"""Utilities to scrape external storm-report catalogs."""
from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import geopandas as gpd
import pandas as pd
import requests
from shapely.geometry import Polygon, shape

from .aoi_utils import aoi_bounds_wgs84, load_aoi

LOGGER = logging.getLogger(__name__)

IEM_LSR_URL = "https://mesonet.agron.iastate.edu/geojson/lsr.php"
BOM_WARNINGS_URL = "https://api.weather.bom.gov.au/v1/warnings"

KNOWN_HAZARD_SYNONYMS = {
    "tornado": {"tornado", "tor", "to"},
    "hail": {"hail", "hi"},
    "wind": {"wind", "tstmwind", "thunderstormwind", "ws", "ds", "gw"},
    "flood": {"flood", "ff", "fl"},
}


@dataclass
class ScrapeResult:
    """Container describing scraped storm reports."""

    events: gpd.GeoDataFrame
    provider: str
    parameters: dict


def fetch_bom_warnings(
    *,
    start: datetime | None = None,
    end: datetime | None = None,
    bounds: Sequence[float] | None = None,
    hazards: Iterable[str] | None = None,
    states: Iterable[str] | None = None,
    phenomena: Iterable[str] | None = None,
    session: requests.Session | None = None,
    timeout: float = 30.0,
) -> gpd.GeoDataFrame:
    """Fetch active Bureau of Meteorology warnings as a GeoDataFrame.

    Parameters
    ----------
    start, end:
        Optional UTC datetimes used to filter warnings locally by their
        ``issuedAt``/``expiresAt`` timestamps. ``start`` must be earlier than
        ``end`` when both are provided.
    bounds:
        Optional ``(min_lon, min_lat, max_lon, max_lat)`` bounding box used to
        retain warnings that intersect the area of interest.
    hazards:
        Optional iterable of hazard keywords used to filter warnings by
        phenomenon/headline text. Comparison is case-insensitive and tokenised
        similarly to :func:`fetch_iem_local_storm_reports`.
    states:
        Optional iterable of Australian state or territory abbreviations (e.g.,
        ``"QLD"``) that the warning must reference.
    phenomena:
        Optional iterable of exact phenomenon names published by the API (for
        example ``"Severe Thunderstorm Warning"``).
    session:
        Optional :class:`requests.Session` used for HTTP requests.
    timeout:
        Timeout (in seconds) passed to :func:`requests.get`.
    """

    if start and end and start >= end:
        raise ValueError("start must be earlier than end")

    http = session or requests
    headers = {
        "Accept": "application/geo+json, application/json",
        "User-Agent": "tornado-ground-scoring/1.0",
    }

    LOGGER.info("Fetching active BOM warnings")
    response = http.get(BOM_WARNINGS_URL, headers=headers, timeout=timeout)
    response.raise_for_status()
    payload = response.json()

    features = _extract_bom_features(payload)
    if not features:
        return _empty_events_gdf()

    hazard_tokens = _normalise_hazard_filters(hazards)
    state_filter = {_normalise_state(code) for code in states or [] if code}
    phenomenon_filter = {
        value.strip().lower() for value in phenomena or [] if isinstance(value, str)
    }
    bbox_geometry = _bounds_to_geometry(bounds) if bounds else None

    records: list[dict] = []
    geometries: list = []

    for feature in features:
        properties = feature.get("properties") or {}

        issued = _parse_bom_datetime(
            properties.get("issuedAt")
            or properties.get("issueTime")
            or properties.get("sent")
        )
        expires = _parse_bom_datetime(properties.get("expiresAt"))

        if start and issued and issued < start:
            continue
        if end and issued and issued > end:
            continue
        if end and expires and expires < start:
            continue

        feature_state = _extract_bom_state(properties)
        if state_filter and feature_state and feature_state not in state_filter:
            continue

        phenomenon_name = _extract_bom_phenomenon(properties)
        if phenomenon_filter and phenomenon_name.lower() not in phenomenon_filter:
            continue

        hazard_label = _extract_bom_hazard_label(properties, phenomenon_name)
        if hazard_tokens:
            hazard_text = hazard_label.lower()
            hazard_compact = hazard_text.replace(" ", "").replace("-", "")
            if not any(
                token in hazard_text or token in hazard_compact
                for token in hazard_tokens
            ):
                continue

        geometry = _parse_bom_geometry(feature)
        if geometry is None:
            continue

        if bbox_geometry is not None and not geometry.intersects(bbox_geometry):
            continue

        centroid = geometry.centroid
        latitude = centroid.y
        longitude = centroid.x

        record = {
            "event_time_utc": issued if issued else expires,
            "expires_time_utc": expires,
            "hazard": hazard_label,
            "latitude": latitude,
            "longitude": longitude,
            "event_date": (
                issued.date() if isinstance(issued, datetime) else None
            ),
            "warning_id": feature.get("id") or properties.get("id"),
            "headline": properties.get("headline") or properties.get("title"),
            "description": properties.get("description")
            or properties.get("instruction"),
            "phenomenon": phenomenon_name,
            "state": feature_state,
            "source": "bom_warnings",
            "area_name": properties.get("area") or properties.get("areaName"),
            "area_sq_km": _calculate_area_sqkm(geometry),
        }
        record["details"] = record.get("description")
        records.append(record)
        geometries.append(geometry)

    if not records:
        return _empty_events_gdf()

    df = pd.DataFrame.from_records(records)
    df["event_time_utc"] = pd.to_datetime(
        df["event_time_utc"], utc=True, errors="coerce"
    )
    df["event_date"] = df["event_time_utc"].dt.date
    geometry_series = gpd.GeoSeries(geometries, crs="EPSG:4326")
    gdf = gpd.GeoDataFrame(df, geometry=geometry_series, crs="EPSG:4326")
    return gdf.reset_index(drop=True)


def fetch_iem_local_storm_reports(
    *,
    start: datetime,
    end: datetime,
    bounds: Sequence[float],
    hazards: Iterable[str] | None = None,
    session: requests.Session | None = None,
    timeout: float = 30.0,
) -> gpd.GeoDataFrame:
    """Fetch Local Storm Reports from the IEM GeoJSON feed.

    Parameters
    ----------
    start, end:
        UTC datetimes bounding the requested reporting window. ``end`` must be
        greater than ``start``.
    bounds:
        Four element sequence ``(min_lon, min_lat, max_lon, max_lat)`` limiting
        the geographic extent of the query in WGS84.
    hazards:
        Optional iterable of hazard names used to filter reports locally. The
        comparison is case-insensitive and matches against the ``typetext``
        attribute provided by the API (for example, ``"Tornado"``).
    session:
        Optional :class:`requests.Session` to reuse connections.
    timeout:
        Timeout (in seconds) passed to :func:`requests.get`.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame containing storm reports. It includes at least the
        ``event_time_utc``, ``hazard``, ``latitude``, ``longitude``,
        ``event_date``, and ``details`` columns used by :mod:`storm_filter`.
    """

    if start >= end:
        raise ValueError("start must be earlier than end")
    if len(bounds) != 4:
        raise ValueError("bounds must be a four element sequence")

    http = session or requests
    params = {
        "sts": _format_utc_iso(start),
        "ets": _format_utc_iso(end),
        "bbox": ",".join(f"{value:.6f}" for value in bounds),
        "fmt": "geojson",
    }

    LOGGER.info(
        "Fetching storm reports from IEM between %s and %s within %s",
        params["sts"],
        params["ets"],
        params["bbox"],
    )

    response = http.get(IEM_LSR_URL, params=params, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    features = payload.get("features", [])

    hazard_filter = _normalise_hazard_filters(hazards)
    if not hazard_filter:
        hazard_filter = None

    records: list[dict] = []
    for feature in features:
        geometry = feature.get("geometry") or {}
        coordinates = geometry.get("coordinates") or []
        if len(coordinates) < 2:
            continue
        lon, lat = float(coordinates[0]), float(coordinates[1])

        properties = feature.get("properties") or {}
        hazard_name = _extract_hazard(properties)
        if hazard_filter:
            hazard_tokens = _tokenize_hazard(hazard_name)
            phenomenon = str(properties.get("phenomena", ""))
            hazard_tokens.update(_tokenize_hazard(phenomenon))
            if hazard_tokens.isdisjoint(hazard_filter):
                continue

        timestamp = _extract_timestamp(properties)
        if timestamp is None:
            continue

        record = {
            "event_time_utc": timestamp,
            "hazard": hazard_name,
            "latitude": lat,
            "longitude": lon,
            "magnitude": properties.get("magnitude"),
            "source": properties.get("source"),
            "city": properties.get("city"),
            "county": properties.get("county"),
            "state": properties.get("state"),
            "details": properties.get("remark")
            or properties.get("comments")
            or properties.get("event_narrative"),
        }
        records.append(record)

    if not records:
        return _empty_events_gdf()

    df = pd.DataFrame.from_records(records)
    df.sort_values("event_time_utc", inplace=True)
    geometry = gpd.points_from_xy(df["longitude"], df["latitude"], crs="EPSG:4326")
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    gdf["event_date"] = gdf["event_time_utc"].dt.date
    return gdf.reset_index(drop=True)


def scrape_iem_catalog(
    *,
    aoi_path: Path,
    start: datetime,
    end: datetime,
    hazards: Iterable[str] | None,
    bbox_buffer_km: float = 0.0,
    session: requests.Session | None = None,
) -> ScrapeResult:
    """Scrape an AOI-constrained catalog of IEM Local Storm Reports."""

    aoi = load_aoi(str(aoi_path))
    bounds = buffered_bounds(aoi, bbox_buffer_km)
    events = fetch_iem_local_storm_reports(
        start=start, end=end, bounds=bounds, hazards=hazards, session=session
    )
    return ScrapeResult(
        events=events,
        provider="iem_lsr",
        parameters={
            "start": _format_utc_iso(start),
            "end": _format_utc_iso(end),
            "bounds": bounds,
            "hazards": list(hazards or []),
            "bbox_buffer_km": bbox_buffer_km,
        },
    )


def buffered_bounds(aoi: gpd.GeoDataFrame, buffer_km: float = 0.0) -> tuple[float, float, float, float]:
    """Return AOI bounds expanded by ``buffer_km`` on all sides."""

    minx, miny, maxx, maxy = aoi_bounds_wgs84(aoi)
    if buffer_km <= 0:
        return (minx, miny, maxx, maxy)

    lat_center = (miny + maxy) / 2.0
    lat_buffer = buffer_km / 111.0
    lon_scale = math.cos(math.radians(lat_center))
    lon_scale = max(lon_scale, 0.1)
    lon_buffer = buffer_km / (111.320 * lon_scale)

    return (
        minx - lon_buffer,
        miny - lat_buffer,
        maxx + lon_buffer,
        maxy + lat_buffer,
    )


def write_catalog_csv(events: gpd.GeoDataFrame, path: Path) -> None:
    """Persist scraped events to a CSV compatible with :mod:`storm_filter`."""

    path.parent.mkdir(parents=True, exist_ok=True)
    df = events.copy()
    if "geometry" in df.columns:
        df = df.drop(columns=["geometry"])
    if "event_time_utc" in df.columns:
        times = pd.to_datetime(df["event_time_utc"], utc=True, errors="coerce")
        formatted = times.dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        df["event_time_utc"] = formatted.where(times.notna(), None)
    df.to_csv(path, index=False)
    LOGGER.info("Wrote %s scraped events to %s", len(df), path)


def save_metadata(metadata: dict, path: Path) -> None:
    """Write scrape metadata for reproducibility."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def _format_utc_iso(value: datetime) -> str:
    ts = pd.to_datetime(value, utc=True)
    return ts.strftime("%Y-%m-%dT%H:%M:%SZ")


def _extract_bom_features(payload: dict) -> list[dict]:
    """Return the feature list from a BOM warnings payload."""

    if not isinstance(payload, dict):
        return []
    for key in ["data", "warnings", "features"]:
        features = payload.get(key)
        if isinstance(features, list):
            return features
    return []


def _normalise_state(value: str | None) -> str | None:
    if not isinstance(value, str):
        return None
    text = value.strip().upper()
    if not text:
        return None
    # Accept both short codes (NSW) and longer state names.
    mapping = {
        "NEW SOUTH WALES": "NSW",
        "VICTORIA": "VIC",
        "QUEENSLAND": "QLD",
        "WESTERN AUSTRALIA": "WA",
        "SOUTH AUSTRALIA": "SA",
        "TASMANIA": "TAS",
        "NORTHERN TERRITORY": "NT",
        "AUSTRALIAN CAPITAL TERRITORY": "ACT",
    }
    return mapping.get(text, text)


def _bounds_to_geometry(bounds: Sequence[float]):
    minx, miny, maxx, maxy = bounds
    return Polygon(
        [
            (minx, miny),
            (minx, maxy),
            (maxx, maxy),
            (maxx, miny),
            (minx, miny),
        ]
    )


def _parse_bom_datetime(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    else:
        parsed = parsed.astimezone(timezone.utc)
    return parsed


def _extract_bom_state(properties: dict) -> str | None:
    for key in ["state", "stateCode", "warningState", "region", "areaState"]:
        value = properties.get(key)
        normalized = _normalise_state(value)
        if normalized:
            return normalized
    return None


def _extract_bom_phenomenon(properties: dict) -> str:
    for key in ["phenomenon", "event", "productType", "warningType"]:
        value = properties.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return "Unknown"


def _extract_bom_hazard_label(properties: dict, default: str) -> str:
    for key in ["headline", "event", "phenomenon", "title", "warningType"]:
        value = properties.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return default or "Unknown"


def _parse_bom_geometry(feature: dict):
    geometry = feature.get("geometry")
    if not geometry:
        areas = feature.get("properties", {}).get("areas")
        if isinstance(areas, list):
            for area in areas:
                geom = area.get("geometry") if isinstance(area, dict) else None
                if geom:
                    geometry = geom
                    break
    if not geometry:
        return None
    try:
        geom = shape(geometry)
    except Exception:  # pragma: no cover - defensive
        return None
    if geom.is_empty:
        return None
    if geom.geom_type == "GeometryCollection":
        parts = [part for part in geom.geoms if not part.is_empty]
        if not parts:
            return None
        geom = parts[0]
    return geom


def _calculate_area_sqkm(geometry) -> float | None:
    if geometry is None or geometry.is_empty:
        return None
    try:
        series = gpd.GeoSeries([geometry], crs="EPSG:4326")
        projected = series.to_crs(3577)
        area_sqkm = float(projected.area.iloc[0] / 1_000_000.0)
    except Exception:  # pragma: no cover - reprojection guard
        return None
    if math.isfinite(area_sqkm):
        return area_sqkm
    return None


def _extract_timestamp(properties: dict) -> pd.Timestamp | None:
    """Parse the best available timestamp from LSR properties."""

    for key in ["valid", "valid_utc", "issue", "event_time_utc", "time"]:
        candidate = properties.get(key)
        if not candidate:
            continue
        ts = pd.to_datetime(candidate, utc=True, errors="coerce")
        if pd.notna(ts):
            return ts
    return None


def _extract_hazard(properties: dict) -> str:
    """Return the best hazard label from the GeoJSON properties."""

    for key in ["typetext", "event_type", "phenomena", "hazard"]:
        value = properties.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return "Unknown"


def _empty_events_gdf() -> gpd.GeoDataFrame:
    df = pd.DataFrame(
        columns=[
            "event_time_utc",
            "hazard",
            "latitude",
            "longitude",
            "event_date",
            "details",
        ]
    )
    geometry = gpd.GeoSeries([], dtype="geometry", crs="EPSG:4326")
    return gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")


def _normalise_hazard_filters(hazards: Iterable[str] | None) -> set[str]:
    values: set[str] = set()
    for raw in hazards or []:
        if not isinstance(raw, str):
            continue
        text = raw.strip().lower()
        if not text:
            continue
        values.add(text)
        values.add(text.replace(" ", ""))
        values.add(text.replace("-", ""))
        values.update(KNOWN_HAZARD_SYNONYMS.get(text, set()))
    values = {token.lower() for token in values if token}
    return values


def _tokenize_hazard(value: str) -> set[str]:
    if not isinstance(value, str):
        return set()
    text = value.strip().lower()
    if not text:
        return set()
    tokens = {
        text,
        text.replace(" ", ""),
        text.replace("-", ""),
    }
    tokens.update(KNOWN_HAZARD_SYNONYMS.get(text, set()))
    return {token.lower() for token in tokens if token}


__all__ = [
    "ScrapeResult",
    "buffered_bounds",
    "fetch_bom_warnings",
    "fetch_iem_local_storm_reports",
    "save_metadata",
    "scrape_iem_catalog",
    "write_catalog_csv",
]

