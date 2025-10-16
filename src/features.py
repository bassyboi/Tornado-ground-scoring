"""Polygon feature extraction utilities."""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Iterable, Sequence

import geopandas as gpd
import numpy as np
import rioxarray
import xarray as xr
from rasterio import features
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon, mapping
from shapely.geometry.base import BaseGeometry

from . import change_core

LOGGER = logging.getLogger(__name__)


def extract(
    polygons: gpd.GeoDataFrame,
    *,
    target_epsg: int,
    pre_path: Path,
    post_path: Path,
    change_path: Path,
    min_hole_area: float = 0.0,
) -> gpd.GeoDataFrame:
    """Return polygons with geometric and raster-derived attributes."""

    if polygons.empty:
        return polygons.copy()

    metric = polygons.to_crs(target_epsg)
    metric["geometry"] = metric.geometry.apply(lambda geom: _clean_geometry(geom, min_hole_area))
    metric = metric.dropna(subset=["geometry"]).reset_index(drop=True)
    if metric.empty:
        LOGGER.warning("All polygons dropped after geometry fixes")
        return gpd.GeoDataFrame(geometry=[], crs=polygons.crs or 4326)

    oriented = [_oriented_metrics(geom) for geom in metric.geometry]
    metric["length_m"] = [vals[0] for vals in oriented]
    metric["width_m"] = [vals[1] for vals in oriented]
    metric["elongation"] = [_safe_divide(vals[0], vals[1]) for vals in oriented]
    metric["bearing_deg"] = [vals[2] for vals in oriented]
    metric["area_m2"] = metric.geometry.area
    metric["perimeter_m"] = metric.geometry.length
    metric["compactness"] = [_compactness(geom) for geom in metric.geometry]
    metric["solidity"] = [_solidity(geom) for geom in metric.geometry]

    centroids = metric.to_crs(4326).geometry.centroid
    metric["centroid_lat"] = centroids.y
    metric["centroid_lon"] = centroids.x

    pre = _open(pre_path)
    post = _open(post_path)
    change = _open(change_path)
    if pre.rio.crs != post.rio.crs:
        raise ValueError("Pre and post rasters have mismatched CRS")

    polygons_raster = metric.to_crs(pre.rio.crs)
    pre_ndvi = change_core.ndvi(_band(pre, "nir"), _band(pre, "red"))
    post_ndvi = change_core.ndvi(_band(post, "nir"), _band(post, "red"))
    delta_ndvi = post_ndvi - pre_ndvi
    pre_bright = change_core.brightness(_band(pre, "red"), _band(pre, "green"), _band(pre, "blue"))
    post_bright = change_core.brightness(_band(post, "red"), _band(post, "green"), _band(post, "blue"))
    delta_bright = post_bright - pre_bright

    stats = {
        "pre_ndvi": _zonal_stats(pre_ndvi, polygons_raster.geometry),
        "post_ndvi": _zonal_stats(post_ndvi, polygons_raster.geometry),
        "delta_ndvi": _zonal_stats(delta_ndvi, polygons_raster.geometry),
        "delta_brightness": _zonal_stats(delta_bright, polygons_raster.geometry),
        "change_score": _zonal_stats(_ensure_2d(change), polygons_raster.geometry),
    }
    for prefix, (means, medians) in stats.items():
        metric[f"{prefix}_mean"] = means
        metric[f"{prefix}_median"] = medians

    return metric.to_crs(polygons.crs or 4326)


def _open(path: Path) -> xr.DataArray:
    if not path.exists():
        raise FileNotFoundError(path)
    data = rioxarray.open_rasterio(path)
    return data if "band" in data.dims else data.expand_dims({"band": [1]})


def _ensure_2d(arr: xr.DataArray) -> xr.DataArray:
    return arr.isel(band=0, drop=True) if "band" in arr.dims else arr


def _band(data: xr.DataArray, name: str) -> xr.DataArray:
    if "band" not in data.coords:
        raise KeyError("Raster missing 'band' coordinate")
    values = [str(v) for v in data.coords["band"].values]
    lookup = {v.lower(): original for v, original in zip(values, data.coords["band"].values)}
    key = lookup.get(name.lower())
    if key is None:
        raise KeyError(f"Band '{name}' not found")
    return data.sel(band=key)


def _zonal_stats(arr: xr.DataArray, geoms: Sequence[BaseGeometry]) -> tuple[list[float], list[float]]:
    raster = _ensure_2d(arr)
    data = raster.values.astype(float)
    transform = raster.rio.transform()
    shape = data.shape
    means: list[float] = []
    medians: list[float] = []
    for geom in geoms:
        if not geom or geom.is_empty:
            means.append(float("nan"))
            medians.append(float("nan"))
            continue
        mask = features.geometry_mask([mapping(geom)], out_shape=shape, transform=transform, invert=True)
        masked = np.where(mask, data, np.nan)
        valid = masked[np.isfinite(masked)]
        means.append(float(np.nanmean(valid)) if valid.size else float("nan"))
        medians.append(float(np.nanmedian(valid)) if valid.size else float("nan"))
    return means, medians


def _clean_geometry(geom: BaseGeometry, min_hole_area: float) -> BaseGeometry | None:
    if not geom or geom.is_empty:
        return None
    fixed = geom.buffer(0)
    if fixed.is_empty:
        return None
    if isinstance(fixed, Polygon):
        return _remove_small_holes(fixed, min_hole_area)
    if isinstance(fixed, MultiPolygon):
        parts = [_remove_small_holes(part, min_hole_area) for part in fixed.geoms]
        valid = [part for part in parts if part and not part.is_empty]
        if not valid:
            return None
        merged = GeometryCollection(valid).unary_union
        return merged if isinstance(merged, (Polygon, MultiPolygon)) else None
    return None


def _remove_small_holes(poly: Polygon, min_hole_area: float) -> Polygon:
    if min_hole_area <= 0:
        return poly
    holes = [ring for ring in poly.interiors if Polygon(ring).area >= min_hole_area]
    cleaned = Polygon(poly.exterior, holes).buffer(0)
    return cleaned if cleaned and not cleaned.is_empty else poly


def _oriented_metrics(geom: BaseGeometry) -> tuple[float, float, float]:
    if not geom or geom.is_empty:
        return 0.0, 0.0, 0.0
    rect = geom.minimum_rotated_rectangle
    coords = list(rect.exterior.coords)
    if len(coords) < 5:
        return 0.0, 0.0, 0.0
    edges = [((p2[0] - p1[0]), (p2[1] - p1[1])) for p1, p2 in zip(coords[:-1], coords[1:])]
    lengths = [math.hypot(dx, dy) for dx, dy in edges if dx or dy]
    if not lengths:
        return 0.0, 0.0, 0.0
    length = max(lengths)
    width = min(lengths)
    dx, dy = max(edges, key=lambda item: math.hypot(item[0], item[1]))
    angle = math.degrees(math.atan2(dx, dy))
    return length, width, (angle + 360.0) % 180.0


def _compactness(geom: BaseGeometry) -> float:
    if not geom or geom.is_empty:
        return float("nan")
    perimeter = geom.length
    return 4.0 * math.pi * geom.area / (perimeter ** 2) if perimeter > 0 else float("nan")


def _solidity(geom: BaseGeometry) -> float:
    if not geom or geom.is_empty:
        return float("nan")
    hull_area = geom.convex_hull.area
    return geom.area / hull_area if hull_area > 0 else float("nan")


def _safe_divide(a: float, b: float) -> float:
    return float("inf") if b == 0 and a > 0 else (a / b if b else float("nan"))
