"""Thresholding helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
import xarray as xr
from skimage.filters import threshold_otsu


def apply_threshold(score: xr.DataArray, threshold: Union[str, float], output_path: Path) -> xr.DataArray:
    """Apply a threshold to the change score and write a mask GeoTIFF."""

    if isinstance(threshold, str):
        if threshold.lower() != "otsu":
            raise ValueError("Only 'otsu' string threshold is supported")
        flat = score.values[np.isfinite(score.values)]
        if flat.size == 0:
            raise ValueError("Cannot compute Otsu threshold on an empty array")
        cutoff = float(threshold_otsu(flat))
    else:
        cutoff = float(threshold)

    mask = (score >= cutoff).astype("uint8")
    mask.rio.write_crs(score.rio.crs, inplace=True)
    mask.rio.write_nodata(0, inplace=True)
    mask.rio.to_raster(output_path)
    return mask
