"""Core change-detection operations."""
from __future__ import annotations

import logging
from typing import Dict

import geopandas as gpd
import numpy as np
import xarray as xr
from rasterio import features
from shapely.geometry import shape
from skimage.filters import threshold_otsu

LOGGER = logging.getLogger(__name__)


def ndvi(nir: xr.DataArray, red: xr.DataArray) -> xr.DataArray:
    """Compute NDVI."""

    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi_arr = (nir - red) / (nir + red)
    return ndvi_arr.fillna(0)


def brightness(red: xr.DataArray, green: xr.DataArray, blue: xr.DataArray) -> xr.DataArray:
    """Simple brightness metric (mean of RGB)."""

    return (red + green + blue) / 3.0


def scl_cloud_mask(scl: xr.DataArray) -> xr.DataArray:
    """Return mask of valid (non-cloud) pixels from the Sentinel-2 SCL band."""

    good_classes = [2, 4, 5, 6, 7, 11]
    mask = xr.zeros_like(scl, dtype=bool)
    for cls in good_classes:
        mask = mask | (scl == cls)
    return mask


def normalize_local(arr: xr.DataArray, win: int = 201) -> xr.DataArray:
    """Local z-score normalization using a square rolling window."""

    if win % 2 == 0:
        win += 1
    arr = arr.where(np.isfinite(arr))
    arr_filled = arr.fillna(0)

    rolling_kwargs = dict(y=win, x=win, center=True, min_periods=1)

    valid = xr.where(arr.notnull(), 1, 0)
    count = valid.rolling(**rolling_kwargs).sum()
    count = count.where(count > 0, 1)

    total = arr_filled.rolling(**rolling_kwargs).sum()
    mean = total / count

    squared = (arr_filled ** 2).rolling(**rolling_kwargs).sum()
    variance = (squared / count) - mean ** 2
    variance = variance.where(variance > 0, 0)
    std = np.sqrt(variance)
    std = std.where(std > 1e-6, 1e-6)

    norm = (arr_filled - mean) / std
    norm = norm.where(arr.notnull(), 0)
    return norm.fillna(0)


def change_score_s2(
    pre: xr.DataArray,
    post: xr.DataArray,
    weights: Dict[str, float],
) -> xr.DataArray:
    """Compute the Sentinel-2-based change score."""

    def band(data: xr.DataArray, name: str) -> xr.DataArray:
        return data.sel(band=name)

    mask = scl_cloud_mask(band(pre, "scl")) & scl_cloud_mask(band(post, "scl"))

    ndvi_pre = ndvi(band(pre, "nir"), band(pre, "red"))
    ndvi_post = ndvi(band(post, "nir"), band(post, "red"))
    d_ndvi = ndvi_pre - ndvi_post

    bright_pre = brightness(band(pre, "red"), band(pre, "green"), band(pre, "blue"))
    bright_post = brightness(band(post, "red"), band(post, "green"), band(post, "blue"))
    d_bright = bright_post - bright_pre

    ndvi_norm = normalize_local(d_ndvi)
    bright_norm = normalize_local(d_bright)

    score = weights.get("d_ndvi", 0.6) * ndvi_norm + weights.get("d_brightness", 0.3) * bright_norm
    score = score.where(mask)
    score = _normalize_01(score)
    return score


def add_s1_logratio(
    score: xr.DataArray,
    s1_pre: xr.DataArray,
    s1_post: xr.DataArray,
    weight: float = 0.1,
) -> xr.DataArray:
    """Blend Sentinel-1 log-ratio magnitude into the change score."""

    eps = 1e-6

    def band(data: xr.DataArray, name: str) -> xr.DataArray:
        return data.sel(band=name)

    ratios = []
    for pol in ("vv", "vh"):
        pre = band(s1_pre, pol).clip(min=eps)
        post = band(s1_post, pol).clip(min=eps)
        ratios.append(np.abs(np.log10(post / pre)))
    log_mag = sum(ratios) / len(ratios)
    log_norm = normalize_local(log_mag)
    combined = score + weight * log_norm
    combined = _normalize_01(combined)
    return combined


def resolve_threshold_value(
    score: xr.DataArray,
    method: str = "otsu",
    numeric: float | None = None,
) -> float:
    """Resolve the numeric threshold for the provided score array."""

    if method == "numeric" and numeric is None:
        raise ValueError("Numeric threshold requested but no value provided")

    valid = score.values[np.isfinite(score.values)]
    if valid.size == 0:
        LOGGER.warning("No valid pixels found for thresholding")
        return float("nan")

    if method == "otsu":
        try:
            thresh = float(threshold_otsu(valid))
        except ValueError:
            thresh = float(np.nanmean(valid)) if valid.size else 0.5
    elif method == "numeric":
        thresh = float(numeric)
    else:
        raise ValueError(f"Unknown threshold method: {method}")
    return thresh


def threshold_score(
    score: xr.DataArray,
    method: str = "otsu",
    numeric: float | None = None,
) -> xr.DataArray:
    """Threshold the change score to produce a binary mask."""

    thresh = resolve_threshold_value(score, method=method, numeric=numeric)
    if not np.isfinite(thresh):
        return xr.zeros_like(score, dtype=bool)
    mask = score > thresh
    return mask.fillna(False)


def polygons_from_mask(mask: xr.DataArray, min_area_m2: float) -> gpd.GeoDataFrame:
    """Vectorise the binary mask into polygons filtered by area."""

    if not hasattr(mask, "rio"):
        raise ValueError("Mask DataArray must have rio spatial metadata")

    transform = mask.rio.transform()
    crs = mask.rio.crs
    arr = mask.fillna(False).astype(np.uint8).values

    shapes = list(features.shapes(arr, mask=arr.astype(bool), transform=transform))
    geoms = []
    for geom, value in shapes:
        if int(value) != 1:
            continue
        poly = shape(geom)
        if poly.is_empty:
            continue
        geoms.append(poly)

    if not geoms:
        return gpd.GeoDataFrame(geometry=[], crs=crs).to_crs(4326)

    gdf = gpd.GeoDataFrame(geometry=geoms, crs=crs)
    gdf["area_m2"] = gdf.to_crs(3857).area
    gdf = gdf[gdf["area_m2"] >= min_area_m2]
    if gdf.empty:
        return gpd.GeoDataFrame(geometry=[], crs=crs).to_crs(4326)
    return gdf.to_crs(4326)


def _normalize_01(arr: xr.DataArray) -> xr.DataArray:
    minimum = arr.min(skipna=True)
    maximum = arr.max(skipna=True)
    if np.isnan(minimum) or np.isnan(maximum) or float(maximum - minimum) < 1e-6:
        return xr.zeros_like(arr)
    normed = (arr - minimum) / (maximum - minimum)
    return normed.clip(0, 1).fillna(0)
