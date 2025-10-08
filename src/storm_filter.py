"""Utilities to filter storm-day events before running the change pipeline."""
from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from typing import Iterable, Sequence

import geopandas as gpd
import pandas as pd

LOGGER = logging.getLogger(__name__)


@dataclass
class CatalogSchema:
    """Column mapping for a tabular storm-event catalog."""

    datetime_column: str = "event_time_utc"
    hazard_column: str = "hazard"
    latitude_column: str = "latitude"
    longitude_column: str = "longitude"


@dataclass
class EventWindow:
    """Container describing the selected storm-event window."""

    events: gpd.GeoDataFrame
    window_start: date
    window_end: date
    shift_days: int = 0


def load_catalog(path: Path, schema: CatalogSchema) -> gpd.GeoDataFrame:
    """Load a CSV catalog of storm events into a GeoDataFrame."""

    if not path.exists():
        raise FileNotFoundError(f"Storm catalog not found: {path}")

    df = _read_catalog(path)
    for column in [
        schema.datetime_column,
        schema.latitude_column,
        schema.longitude_column,
    ]:
        if column not in df.columns:
            raise ValueError(f"Storm catalog missing required column '{column}'")

    dt = pd.to_datetime(df[schema.datetime_column], errors="coerce", utc=True)
    if dt.isna().all():
        raise ValueError(
            "Failed to parse storm catalog datetimes; check the datetime column format"
        )

    df = df.copy()
    df["event_time_utc"] = dt
    df["event_date"] = df["event_time_utc"].dt.date

    hazard_values = df.get(schema.hazard_column, None)
    if hazard_values is not None:
        df["hazard"] = hazard_values.astype(str)
    else:
        df["hazard"] = "Unknown"

    geometry = gpd.points_from_xy(
        pd.to_numeric(df[schema.longitude_column], errors="coerce"),
        pd.to_numeric(df[schema.latitude_column], errors="coerce"),
        crs="EPSG:4326",
    )
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    return gdf


def _read_catalog(path: Path) -> pd.DataFrame:
    """Return a DataFrame from ``path`` with resilient CSV parsing."""

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration as exc:  # pragma: no cover - empty file guard
            raise ValueError("Storm catalog is empty") from exc

        rows = []
        header_len = len(header)
        for row in reader:
            if not row:
                continue
            if len(row) > header_len:
                trimmed = row[: header_len - 1]
                trimmed.append(
                    ",".join(cell for cell in row[header_len - 1 :] if cell)
                )
                row = trimmed
            elif len(row) < header_len:
                row = row + [""] * (header_len - len(row))
            rows.append(row)

    return pd.DataFrame(rows, columns=header)


def relevant_events(
    *,
    catalog: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    window_start: date,
    window_end: date,
    hazards: Sequence[str] | None,
    days_before: int,
    days_after: int,
    distance_km: float,
) -> gpd.GeoDataFrame:
    """Return events inside the AOI (optionally buffered) and date window."""

    if catalog.empty:
        return catalog.iloc[0:0].copy()

    filtered = _filter_by_hazard_and_distance(
        catalog=catalog,
        aoi=aoi,
        hazards=hazards,
        distance_km=distance_km,
    )
    if filtered.empty:
        return filtered

    start_date = window_start - timedelta(days=max(days_before, 0))
    end_date = window_end + timedelta(days=max(days_after, 0))

    mask_date = filtered["event_date"].between(start_date, end_date)
    return filtered.loc[mask_date].reset_index(drop=True)


def resolve_event_window(
    *,
    catalog: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    window_start: date,
    window_end: date,
    hazards: Sequence[str] | None,
    days_before: int,
    days_after: int,
    distance_km: float,
    auto_backfill_max_days: int | str | None = 0,
    auto_backfill_step_days: int = 1,
    auto_backfill_directions: Iterable[str] | str | None = None,
) -> EventWindow:
    """Return the first window with qualifying events, expanding if needed.

    The search begins with ``window_start``/``window_end``. When the filtered
    catalog is empty and ``auto_backfill_max_days`` is greater than zero, the
    function shifts the window by ``auto_backfill_step_days`` (in days) in the
    requested directions until it finds qualifying events or exhausts the
    search distance. The returned :class:`EventWindow` captures the discovered
    events, the resolved window bounds, and the signed day shift applied.
    """

    resolved_directions = _normalize_backfill_directions(auto_backfill_directions)
    max_days = _resolve_max_backfill_days(
        catalog=catalog,
        aoi=aoi,
        window_start=window_start,
        window_end=window_end,
        hazards=hazards,
        days_before=days_before,
        days_after=days_after,
        distance_km=distance_km,
        step_days=auto_backfill_step_days,
        configured_max_days=auto_backfill_max_days,
    )
    step_days = max(int(auto_backfill_step_days or 1), 1)

    # Always test the base window first.
    candidate_offsets = [0]
    if max_days > 0:
        for delta in range(step_days, max_days + 1, step_days):
            if "backward" in resolved_directions:
                candidate_offsets.append(-delta)
            if "forward" in resolved_directions:
                candidate_offsets.append(delta)

    for offset in candidate_offsets:
        attempt_start = window_start + timedelta(days=offset)
        attempt_end = window_end + timedelta(days=offset)
        events = relevant_events(
            catalog=catalog,
            aoi=aoi,
            window_start=attempt_start,
            window_end=attempt_end,
            hazards=hazards,
            days_before=days_before,
            days_after=days_after,
            distance_km=distance_km,
        )
        if not events.empty:
            return EventWindow(
                events=events,
                window_start=attempt_start,
                window_end=attempt_end,
                shift_days=offset,
            )

    empty = catalog.iloc[0:0].copy()
    return EventWindow(
        events=empty,
        window_start=window_start,
        window_end=window_end,
        shift_days=0,
    )


def _resolve_max_backfill_days(
    *,
    catalog: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    window_start: date,
    window_end: date,
    hazards: Sequence[str] | None,
    days_before: int,
    days_after: int,
    distance_km: float,
    step_days: int,
    configured_max_days: int | str | None,
) -> int:
    """Resolve the auto-backfill span in days.

    ``configured_max_days`` accepts an integer number of days, the string
    ``"auto"`` (or ``"catalog"``/``"events"``) to expand to the furthest
    qualifying catalog event, or ``None`` to apply the same behaviour. Any
    other string is parsed as an integer. Negative values are clamped to zero.
    """

    if step_days <= 0:
        step_days = 1

    option = configured_max_days
    if isinstance(option, str):
        text = option.strip().lower()
        if not text:
            option = 0
        elif text in {"auto", "catalog", "events", "full"}:
            option = None
        else:
            try:
                option = int(text)
            except ValueError as exc:  # pragma: no cover - defensive
                raise ValueError(
                    "auto_backfill_max_days must be an integer or 'auto'"
                ) from exc

    if option is None:
        return _auto_backfill_span_from_catalog(
            catalog=catalog,
            aoi=aoi,
            window_start=window_start,
            window_end=window_end,
            hazards=hazards,
            days_before=days_before,
            days_after=days_after,
            distance_km=distance_km,
            step_days=step_days,
        )

    return max(int(option or 0), 0)


def _auto_backfill_span_from_catalog(
    *,
    catalog: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    window_start: date,
    window_end: date,
    hazards: Sequence[str] | None,
    days_before: int,
    days_after: int,
    distance_km: float,
    step_days: int,
) -> int:
    """Return the minimum span (rounded to ``step_days``) covering catalog events."""

    filtered = _filter_by_hazard_and_distance(
        catalog=catalog,
        aoi=aoi,
        hazards=hazards,
        distance_km=distance_km,
    )
    if filtered.empty or "event_date" not in filtered:
        return 0

    event_dates = pd.to_datetime(filtered["event_date"], errors="coerce")
    event_dates = event_dates.dropna().dt.normalize()
    if event_dates.empty:
        return 0

    base_start = pd.Timestamp(window_start) - pd.Timedelta(days=max(days_before, 0))
    base_end = pd.Timestamp(window_end) + pd.Timedelta(days=max(days_after, 0))

    earliest = event_dates.min()
    latest = event_dates.max()

    backward_span = abs(min((earliest - base_end).days, 0))
    forward_span = max((latest - base_start).days, 0)
    max_span = max(backward_span, forward_span)
    if max_span <= 0:
        return 0

    return int(math.ceil(max_span / step_days) * step_days)


def _filter_by_hazard_and_distance(
    *,
    catalog: gpd.GeoDataFrame,
    aoi: gpd.GeoDataFrame,
    hazards: Sequence[str] | None,
    distance_km: float,
) -> gpd.GeoDataFrame:
    """Return catalog rows matching the requested hazards and AOI buffer."""

    if catalog.empty:
        return catalog.iloc[0:0].copy()

    filtered = catalog.copy()
    if hazards:
        normalized_hazards = {h.upper() for h in hazards if isinstance(h, str)}
        hazard_series = filtered["hazard"].astype(str).str.upper()
        filtered = filtered.loc[hazard_series.isin(normalized_hazards)]
    if filtered.empty:
        return filtered.iloc[0:0].copy()

    if distance_km > 0:
        aoi_geom = _buffer_aoi(aoi, distance_km)
        mask_geom = filtered.geometry.within(aoi_geom)
        filtered = filtered.loc[mask_geom]

    return filtered.reset_index(drop=True)


def _normalize_backfill_directions(
    directions: Iterable[str] | str | None,
) -> set[str]:
    """Return a validated set of backfill directions."""

    if directions is None:
        return {"backward"}
    if isinstance(directions, str):
        candidates = [directions]
    else:
        candidates = list(directions)

    normalized = {value.strip().lower() for value in candidates if value}
    if not normalized:
        return {"backward"}
    if "both" in normalized:
        return {"backward", "forward"}

    valid = {"backward", "forward"}
    invalid = normalized - valid
    if invalid:
        raise ValueError(f"Invalid auto_backfill_directions: {sorted(invalid)}")
    return normalized


def _buffer_aoi(aoi: gpd.GeoDataFrame, distance_km: float):
    """Return the AOI geometry buffered by ``distance_km`` (in kilometers)."""

    if aoi.crs is None:
        aoi = aoi.set_crs(4326)
    aoi_wgs84 = aoi.to_crs(4326)
    geom = aoi_wgs84.unary_union
    if distance_km <= 0:
        return geom

    buffer_series = gpd.GeoSeries([geom], crs=4326).to_crs(3857)
    buffered = buffer_series.buffer(distance_km * 1000.0)
    return buffered.to_crs(4326).iloc[0]


def export_events(events: gpd.GeoDataFrame, path: Path) -> None:
    """Write relevant storm events to GeoJSON for downstream inspection."""

    if events.empty:
        collection = {"type": "FeatureCollection", "features": []}
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(collection))
        LOGGER.info("Wrote empty storm-event log to %s", path)
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    if events.crs is None:
        events = events.set_crs(4326)
    events = _stringify_event_fields(events)
    events_wgs84 = events.to_crs(4326)
    geojson = events_wgs84.to_json()
    path.write_text(geojson, encoding="utf-8")
    LOGGER.info("Wrote %s storm events to %s", len(events), path)


def _stringify_event_fields(events: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Convert datetime/date columns to ISO strings for GeoJSON export."""

    events = events.copy()
    if "event_time_utc" in events.columns:
        times = pd.to_datetime(events["event_time_utc"], utc=True, errors="coerce")
        events["event_time_utc"] = times.dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    if "event_date" in events.columns:
        events["event_date"] = events["event_date"].apply(_date_to_iso)
    return events


def _date_to_iso(value: object) -> str | None:
    """Return an ISO 8601 date string or ``None`` when value is missing."""

    if pd.isna(value):
        return None
    if isinstance(value, date):
        return value.isoformat()
    if isinstance(value, pd.Timestamp):
        if value.tzinfo:
            return value.tz_convert("UTC").date().isoformat()
        return value.date().isoformat()
    if isinstance(value, str):
        try:
            parsed = pd.to_datetime(value, errors="coerce")
        except Exception:  # pragma: no cover - defensive
            return value
        if pd.isna(parsed):
            return value
        if parsed.tzinfo:
            parsed = parsed.tz_convert("UTC")
        return parsed.date().isoformat()
    return str(value)
