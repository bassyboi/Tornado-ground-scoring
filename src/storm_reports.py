"""Utilities to scrape external storm-report catalogs."""
from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

import geopandas as gpd
import pandas as pd
import requests

from .aoi_utils import aoi_bounds_wgs84, load_aoi

LOGGER = logging.getLogger(__name__)

IEM_LSR_URL = "https://mesonet.agron.iastate.edu/geojson/lsr.php"

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
    "fetch_iem_local_storm_reports",
    "save_metadata",
    "scrape_iem_catalog",
    "write_catalog_csv",
]

