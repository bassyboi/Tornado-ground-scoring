"""Utilities to filter storm-day events before running the change pipeline."""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from typing import Sequence

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


def load_catalog(path: Path, schema: CatalogSchema) -> gpd.GeoDataFrame:
    """Load a CSV catalog of storm events into a GeoDataFrame."""

    if not path.exists():
        raise FileNotFoundError(f"Storm catalog not found: {path}")

    df = pd.read_csv(path)
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
        df[schema.longitude_column], df[schema.latitude_column], crs="EPSG:4326"
    )
    gdf = gpd.GeoDataFrame(df, geometry=geometry)
    return gdf


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

    if hazards:
        normalized_hazards = [h.upper() for h in hazards]
        mask_hazard = catalog["hazard"].str.upper().isin(normalized_hazards)
    else:
        mask_hazard = pd.Series(True, index=catalog.index)

    start_date = window_start - timedelta(days=max(days_before, 0))
    end_date = window_end + timedelta(days=max(days_after, 0))

    mask_date = catalog["event_date"].between(start_date, end_date)
    filtered = catalog.loc[mask_hazard & mask_date].copy()
    if filtered.empty:
        return filtered

    aoi_geom = _buffer_aoi(aoi, distance_km)
    mask_geom = filtered.geometry.within(aoi_geom)
    return filtered.loc[mask_geom].reset_index(drop=True)


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
    events.to_crs(4326).to_file(path, driver="GeoJSON")
    LOGGER.info("Wrote %s storm events to %s", len(events), path)
