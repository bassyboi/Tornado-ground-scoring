"""Utilities for loading and manipulating Areas of Interest (AOIs)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import geopandas as gpd
import numpy as np
from pyproj import CRS
from rasterio.transform import from_origin

DEFAULT_PROJECTED_CRS = CRS.from_epsg(3857)


@dataclass
class RasterGrid:
    """Definition of a raster grid derived from an AOI."""

    height: int
    width: int
    transform: any
    crs: CRS

    @property
    def shape(self) -> Tuple[int, int]:
        return self.height, self.width


def load_aoi(path: str) -> gpd.GeoDataFrame:
    """Load an AOI GeoJSON file as a GeoDataFrame."""

    gdf = gpd.read_file(path)
    if gdf.empty:
        raise ValueError("AOI GeoJSON does not contain any features")
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    return gdf


def aoi_bounds_wgs84(gdf: gpd.GeoDataFrame) -> Tuple[float, float, float, float]:
    """Return the AOI bounds in WGS84 (minx, miny, maxx, maxy)."""

    if gdf.crs is None or gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)
    return tuple(gdf.total_bounds)


def raster_grid_from_aoi(
    gdf: gpd.GeoDataFrame,
    resolution: float = 10.0,
    projected_crs: CRS = DEFAULT_PROJECTED_CRS,
) -> RasterGrid:
    """Create a simple raster grid covering the AOI at the desired resolution.

    The grid is computed in a projected CRS (default Web Mercator) to provide
    approximate metric units for placeholder rasters when STAC data are not
    available locally.
    """

    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    gdf_proj = gdf.to_crs(projected_crs)
    minx, miny, maxx, maxy = gdf_proj.total_bounds
    span_x = max(maxx - minx, resolution)
    span_y = max(maxy - miny, resolution)
    width = max(1, int(np.ceil(span_x / resolution)))
    height = max(1, int(np.ceil(span_y / resolution)))
    transform = from_origin(minx, maxy, resolution, resolution)
    return RasterGrid(height=height, width=width, transform=transform, crs=projected_crs)
