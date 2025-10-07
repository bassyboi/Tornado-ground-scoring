"""Ground-Scour Score utilities."""
from __future__ import annotations

import math
from typing import Sequence

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon


def major_axis_orientation(poly: Polygon) -> float:
    """Return the major-axis orientation (0-180 degrees from North)."""

    rect = poly.minimum_rotated_rectangle
    coords = list(rect.exterior.coords)[:-1]
    edges = []
    for idx in range(len(coords)):
        p1 = coords[idx]
        p2 = coords[(idx + 1) % len(coords)]
        length = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
        edges.append((length, p1, p2))
    if not edges:
        return 0.0
    major = max(edges, key=lambda e: e[0])
    dx = major[2][0] - major[1][0]
    dy = major[2][1] - major[1][1]
    angle = math.degrees(math.atan2(dx, dy))
    bearing = abs(angle) % 180
    return bearing


def elongation_ratio(poly: Polygon) -> float:
    """Compute the elongation ratio from the minimum rotated rectangle."""

    rect = poly.minimum_rotated_rectangle
    coords = list(rect.exterior.coords)[:-1]
    if len(coords) < 4:
        return 0.0
    lengths = []
    for idx in range(len(coords)):
        p1 = coords[idx]
        p2 = coords[(idx + 1) % len(coords)]
        lengths.append(math.hypot(p2[0] - p1[0], p2[1] - p1[1]))
    lengths = sorted(lengths, reverse=True)
    if len(lengths) < 2 or lengths[1] == 0:
        return 0.0
    return lengths[0] / lengths[1]


def apply_elongation_filter(
    gdf: gpd.GeoDataFrame,
    preferred_bearing: float,
    ang_tol: float = 30.0,
    min_ratio: float = 2.0,
) -> gpd.GeoDataFrame:
    """Filter polygons by elongation and orientation."""

    if gdf.empty:
        return gdf

    gdf = gdf.copy()
    gdf["bearing_deg"] = gdf.geometry.apply(major_axis_orientation)
    gdf["elongation_ratio"] = gdf.geometry.apply(elongation_ratio)

    def ang_diff(angle: float) -> float:
        diff = abs((angle % 180) - (preferred_bearing % 180))
        return min(diff, 180 - diff)

    mask = (gdf["elongation_ratio"] >= min_ratio) & (gdf["bearing_deg"].apply(ang_diff) <= ang_tol)
    return gdf.loc[mask]


def _as_float_array(values: Sequence[float] | np.ndarray) -> np.ndarray:
    """Return the input values as a 1-D float array."""

    if isinstance(values, np.ndarray):
        return values.astype(float)
    return np.asarray(list(values), dtype=float)


def resolve_gss_thresholds(
    values: Sequence[float], breaks: Sequence[float] | str = "quantile"
) -> np.ndarray:
    """Resolve the five Ground-Scour Score thresholds for the provided scores."""

    arr = _as_float_array(values)

    if isinstance(breaks, str):
        if breaks != "quantile":
            raise ValueError(f"Unknown breaks mode: {breaks}")
        valid = arr[arr > 0]
        if valid.size >= 5:
            probs = np.linspace(0, 1, 7)[1:-1]
            thresholds = np.quantile(valid, probs)
        else:
            thresholds = np.linspace(0.1, 0.9, 5)
    else:
        thresholds = np.asarray(list(breaks), dtype=float)
        if thresholds.size != 5:
            raise ValueError("Break list must contain exactly five thresholds")
    return np.sort(thresholds)


def score_to_gss(values: Sequence[float], breaks: Sequence[float] | str = "quantile") -> np.ndarray:
    """Map change scores to Ground-Scour Scores (0-5)."""

    arr = _as_float_array(values)
    thresholds = resolve_gss_thresholds(arr, breaks=breaks)
    gss_values = np.zeros(arr.shape, dtype=int)

    for idx, thresh in enumerate(thresholds):
        gss_values[arr > thresh] = idx + 1
    gss_values[arr > thresholds[-1]] = 5
    return gss_values
