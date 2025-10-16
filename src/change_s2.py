"""Sentinel-2 change metrics for the thin-slice pipeline."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
import xarray as xr

from .config import WeightConfig


@dataclass
class ChangeProducts:
    """Container for raster outputs generated during change detection."""

    pre_mosaic_path: Path
    post_mosaic_path: Path
    score_path: Path


def compute_change_score(
    pre: xr.Dataset,
    post: xr.Dataset,
    weights: WeightConfig,
    artifacts_dir: Path,
) -> Tuple[xr.DataArray, ChangeProducts]:
    """Compute the composite change score and persist intermediate rasters."""

    artifacts_dir.mkdir(parents=True, exist_ok=True)

    pre_path = artifacts_dir / "pre_mosaic.tif"
    post_path = artifacts_dir / "post_mosaic.tif"
    score_path = artifacts_dir / "change_score.tif"

    _write_multiband(pre, pre_path)
    _write_multiband(post, post_path)

    pre_ndvi = _ndvi(pre["nir"], pre["red"])
    post_ndvi = _ndvi(post["nir"], post["red"])
    delta_ndvi = (pre_ndvi - post_ndvi).clip(min=0)

    pre_brightness = _brightness(pre)
    post_brightness = _brightness(post)
    delta_brightness = (post_brightness - pre_brightness).clip(min=0)

    norm_ndvi = _normalize(delta_ndvi)
    norm_brightness = _normalize(delta_brightness)

    total_weight = weights.total
    score = (weights.ndvi * norm_ndvi + weights.brightness * norm_brightness) / total_weight
    score = score.clip(min=0.0, max=1.0)
    score = score.where(np.isfinite(score), 0.0).astype("float32")
    score.rio.write_crs(pre.rio.crs, inplace=True)

    score.rio.write_nodata(0.0, inplace=True)
    score.rio.to_raster(score_path)

    return score, ChangeProducts(pre_mosaic_path=pre_path, post_mosaic_path=post_path, score_path=score_path)


def _ndvi(nir: xr.DataArray, red: xr.DataArray) -> xr.DataArray:
    numerator = nir - red
    denominator = nir + red
    ndvi = numerator / xr.where(denominator == 0, np.nan, denominator)
    return ndvi


def _brightness(mosaic: xr.Dataset) -> xr.DataArray:
    return (mosaic["red"] + mosaic["green"] + mosaic["blue"]) / 3.0


def _normalize(arr: xr.DataArray) -> xr.DataArray:
    finite = arr.where(np.isfinite(arr))
    min_val = finite.min()
    max_val = finite.max()
    spread = max_val - min_val
    normalized = (finite - min_val) / xr.where(spread == 0, 1, spread)
    return normalized.fillna(0.0).astype("float32")


def _write_multiband(mosaic: xr.Dataset, path: Path) -> None:
    array = mosaic.to_array(dim="band")
    array.rio.write_crs(mosaic.rio.crs, inplace=True)
    array.rio.write_nodata(np.nan, inplace=True)
    array.rio.to_raster(path)
