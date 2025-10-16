"""AOI helpers for reading and projecting geometries."""
from __future__ import annotations

from pathlib import Path
from typing import Tuple

import geopandas as gpd
from pyproj import CRS
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union


class AOILoadError(RuntimeError):
    """Raised when the AOI file cannot be read or contains no features."""


def load_aoi(
    path: Path, target_epsg: int
) -> Tuple[BaseGeometry, BaseGeometry, Tuple[float, float, float, float], CRS]:
    """Load an AOI file and return projected and geographic geometries.

    Parameters
    ----------
    path:
        Path to a GeoJSON or GeoPackage describing the AOI features.
    target_epsg:
        EPSG code for the projected CRS that downstream rasters should use.

    Returns
    -------
    projected_geom, geographic_geom, wgs84_bounds, target_crs
        ``projected_geom`` is dissolved into a single polygon in the target CRS.
        ``geographic_geom`` is the unioned polygon in WGS84, suitable for STAC searches.
        ``wgs84_bounds`` provides the bounding box (minx, miny, maxx, maxy) in WGS84.
        ``target_crs`` is the CRS object associated with ``target_epsg``.
    """

    if not Path(path).exists():  # pragma: no cover - defensive
        raise FileNotFoundError(f"AOI file not found: {path}")

    try:
        gdf = gpd.read_file(path)
    except Exception as exc:  # pragma: no cover - defensive
        raise AOILoadError(f"Failed to read AOI file {path}: {exc}") from exc

    if gdf.empty:
        raise AOILoadError(f"AOI file {path} contains no features")

    target_crs = CRS.from_epsg(int(target_epsg))
    wgs84 = CRS.from_epsg(4326)

    gdf = gdf.to_crs(wgs84)
    geographic_geom = unary_union(gdf.geometry)
    if geographic_geom.is_empty:
        raise AOILoadError("AOI geometry is empty after dissolve")
    geographic_geom = geographic_geom.buffer(0)

    projected = gdf.to_crs(target_crs)
    projected_geom = unary_union(projected.geometry).buffer(0)
    if projected_geom.is_empty:
        raise AOILoadError("AOI geometry is empty after projection")

    bounds = geographic_geom.bounds
    return projected_geom, geographic_geom, bounds, target_crs
