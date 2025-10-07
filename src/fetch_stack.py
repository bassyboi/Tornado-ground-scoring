"""STAC search and mosaic helpers for Sentinel data."""
from __future__ import annotations

import logging
from datetime import datetime, timedelta
from typing import Callable, List, Sequence, Tuple

import geopandas as gpd
import numpy as np
import xarray as xr
from pystac_client import Client
import rioxarray  # noqa: F401  # required for rio accessor
import stackstac
from pyproj import CRS

from .aoi_utils import aoi_bounds_wgs84, raster_grid_from_aoi

LOGGER = logging.getLogger(__name__)

DEFAULT_STACK_EPSG = 3857


def _open_catalog(url: str) -> Client | None:
    try:
        return Client.open(url)
    except Exception as exc:  # pragma: no cover - defensive
        LOGGER.warning("Failed to open STAC catalog %s: %s", url, exc)
        return None


def search_sentinel2(
    aoi_gdf: gpd.GeoDataFrame,
    date_from: str,
    date_to: str,
    catalogs: Sequence[str],
) -> List:
    """Search Sentinel-2 L2A scenes across catalogs."""

    bbox = aoi_bounds_wgs84(aoi_gdf)
    items = []
    for url in catalogs:
        client = _open_catalog(url)
        if client is None:
            continue
        search = client.search(
            collections=["sentinel-2-l2a"],
            bbox=bbox,
            datetime=f"{date_from}/{date_to}",
            query={"eo:cloud_cover": {"lt": 90}},
        )
        try:
            items.extend(list(search.items()))
        except Exception as exc:  # pragma: no cover - defensive
            LOGGER.warning("Sentinel-2 search failed for %s: %s", url, exc)
    LOGGER.info("Found %d Sentinel-2 items", len(items))
    return items


def search_sentinel1(
    aoi_gdf: gpd.GeoDataFrame,
    date_from: str,
    date_to: str,
    catalogs: Sequence[str],
) -> List:
    """Search Sentinel-1 RTC scenes."""

    bbox = aoi_bounds_wgs84(aoi_gdf)
    items = []
    for url in catalogs:
        client = _open_catalog(url)
        if client is None:
            continue
        search = client.search(
            collections=["sentinel-1-rtc"],
            bbox=bbox,
            datetime=f"{date_from}/{date_to}",
        )
        try:
            items.extend(list(search.items()))
        except Exception as exc:  # pragma: no cover - defensive
            LOGGER.warning("Sentinel-1 search failed for %s: %s", url, exc)
    LOGGER.info("Found %d Sentinel-1 items", len(items))
    return items


def mosaic_s2(
    items: Sequence,
    bands: Sequence[str],
    aoi_gdf: gpd.GeoDataFrame,
    resolution: float = 10.0,
    target_epsg: int | None = None,
) -> xr.DataArray:
    """Create a median Sentinel-2 mosaic for the requested bands."""

    if not items:
        LOGGER.warning("No Sentinel-2 items found. Using placeholder mosaic.")
        return _placeholder_mosaic(bands, aoi_gdf, resolution, target_epsg)

    try:
        epsg = target_epsg or DEFAULT_STACK_EPSG
        stack = stackstac.stack(items, assets=bands, resolution=resolution, epsg=epsg)
        data = stack.median(dim="time", skipna=True)
        if "band" not in data.dims:
            data = data.expand_dims({"band": list(bands)})
        data = data.transpose("band", "y", "x")
        data = data.assign_coords({"band": list(bands)})
        data = data.rio.write_crs(stack.rio.crs)
        data = data.rio.write_transform(stack.rio.transform())
        return data
    except Exception as exc:  # pragma: no cover - defensive
        LOGGER.warning("Failed to mosaic Sentinel-2 stack: %s", exc)
        return _placeholder_mosaic(bands, aoi_gdf, resolution, target_epsg)


def mosaic_s1(
    items: Sequence,
    bands: Sequence[str],
    aoi_gdf: gpd.GeoDataFrame,
    resolution: float = 10.0,
    target_epsg: int | None = None,
) -> xr.DataArray:
    """Create a median Sentinel-1 RTC mosaic."""

    if not items:
        LOGGER.warning("No Sentinel-1 items found. Using placeholder mosaic.")
        return _placeholder_mosaic(bands, aoi_gdf, resolution, target_epsg)

    try:
        epsg = target_epsg or DEFAULT_STACK_EPSG
        stack = stackstac.stack(items, assets=bands, resolution=resolution, epsg=epsg)
        data = stack.median(dim="time", skipna=True)
        if "band" not in data.dims:
            data = data.expand_dims({"band": list(bands)})
        data = data.transpose("band", "y", "x")
        data = data.assign_coords({"band": list(bands)})
        data = data.rio.write_crs(stack.rio.crs)
        data = data.rio.write_transform(stack.rio.transform())
        return data
    except Exception as exc:  # pragma: no cover - defensive
        LOGGER.warning("Failed to mosaic Sentinel-1 stack: %s", exc)
        return _placeholder_mosaic(bands, aoi_gdf, resolution, target_epsg)


def search_sentinel2_with_fallback(
    aoi_gdf: gpd.GeoDataFrame,
    date_from: str,
    date_to: str,
    catalogs: Sequence[str],
    *,
    max_expansion_days: int = 0,
    step_days: int = 1,
    window_label: str = "",
) -> Tuple[List, str, str]:
    """Search Sentinel-2 scenes and optionally expand the date window."""

    items = search_sentinel2(aoi_gdf, date_from, date_to, catalogs)
    if items or max_expansion_days <= 0:
        return items, date_from, date_to

    return _expand_date_window(
        search_fn=lambda start, end: search_sentinel2(aoi_gdf, start, end, catalogs),
        sensor_label="Sentinel-2",
        window_label=window_label,
        date_from=date_from,
        date_to=date_to,
        max_expansion_days=max_expansion_days,
        step_days=step_days,
    )


def search_sentinel1_with_fallback(
    aoi_gdf: gpd.GeoDataFrame,
    date_from: str,
    date_to: str,
    catalogs: Sequence[str],
    *,
    max_expansion_days: int = 0,
    step_days: int = 1,
    window_label: str = "",
) -> Tuple[List, str, str]:
    """Search Sentinel-1 RTC scenes and optionally expand the date window."""

    items = search_sentinel1(aoi_gdf, date_from, date_to, catalogs)
    if items or max_expansion_days <= 0:
        return items, date_from, date_to

    return _expand_date_window(
        search_fn=lambda start, end: search_sentinel1(aoi_gdf, start, end, catalogs),
        sensor_label="Sentinel-1",
        window_label=window_label,
        date_from=date_from,
        date_to=date_to,
        max_expansion_days=max_expansion_days,
        step_days=step_days,
    )


def _expand_date_window(
    *,
    search_fn: Callable[[str, str], Sequence],
    sensor_label: str,
    window_label: str,
    date_from: str,
    date_to: str,
    max_expansion_days: int,
    step_days: int,
) -> Tuple[List, str, str]:
    """Expand the search window until imagery is found or limits are hit."""

    try:
        start_date = datetime.fromisoformat(date_from).date()
        end_date = datetime.fromisoformat(date_to).date()
    except ValueError as exc:  # pragma: no cover - defensive
        LOGGER.warning(
            "Unable to parse %s window dates %s–%s: %s", sensor_label, date_from, date_to, exc
        )
        return [], date_from, date_to

    step_days = max(1, step_days)
    max_expansion_days = max(0, max_expansion_days)

    for pad_days in range(step_days, max_expansion_days + step_days, step_days):
        expanded_from = (start_date - timedelta(days=pad_days)).isoformat()
        expanded_to = (end_date + timedelta(days=pad_days)).isoformat()
        try:
            items = list(search_fn(expanded_from, expanded_to))
        except Exception as exc:  # pragma: no cover - defensive
            LOGGER.warning(
                "Expanded %s search to %s–%s but it failed: %s",
                sensor_label,
                expanded_from,
                expanded_to,
                exc,
            )
            continue
        if items:
            label = f"{sensor_label} {window_label}".strip()
            LOGGER.info(
                "Expanded %s window to %s–%s (±%s days) to find %d items",
                label or sensor_label,
                expanded_from,
                expanded_to,
                pad_days,
                len(items),
            )
            return list(items), expanded_from, expanded_to

    LOGGER.warning(
        "No %s imagery found within ±%s days of requested window %s–%s",
        sensor_label,
        max_expansion_days,
        date_from,
        date_to,
    )
    return [], date_from, date_to


def _placeholder_mosaic(
    bands: Sequence[str],
    aoi_gdf: gpd.GeoDataFrame,
    resolution: float,
    target_epsg: int | None,
) -> xr.DataArray:
    """Build a zero-filled mosaic so downstream steps can proceed."""

    if target_epsg:
        projected_crs = CRS.from_epsg(target_epsg)
    else:
        projected_crs = CRS.from_epsg(DEFAULT_STACK_EPSG)
    grid = raster_grid_from_aoi(aoi_gdf, resolution=resolution, projected_crs=projected_crs)
    data = np.zeros((len(bands), grid.height, grid.width), dtype=np.float32)
    x_coords = grid.transform.c + (np.arange(grid.width) + 0.5) * grid.transform.a
    y_coords = grid.transform.f + (np.arange(grid.height) + 0.5) * grid.transform.e
    arr = xr.DataArray(
        data,
        dims=("band", "y", "x"),
        coords={"band": list(bands), "y": y_coords, "x": x_coords},
        name="placeholder",
    )
    arr = arr.rio.write_crs(grid.crs)
    arr = arr.rio.write_transform(grid.transform)
    arr = arr.rio.write_nodata(np.nan)
    return arr
