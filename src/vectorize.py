"""Raster to vector helpers for the change mask."""
from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
from rasterio import features
from shapely.geometry import Polygon, shape

DEFAULT_MIN_AREA_M2 = 2000.0


def mask_to_polygons(
    mask: xr.DataArray,
    score: xr.DataArray,
    out_geojson: Path,
    out_kml: Path,
    min_area_m2: float = DEFAULT_MIN_AREA_M2,
) -> gpd.GeoDataFrame:
    """Vectorise a binary mask and write GeoJSON + KML outputs."""

    out_geojson.parent.mkdir(parents=True, exist_ok=True)
    mask_np = mask.values
    transform = mask.rio.transform()
    shapes = features.shapes(mask_np.astype(np.uint8), mask=mask_np.astype(bool), transform=transform)

    geoms: list[Polygon] = []
    for geom, value in shapes:
        if int(value) != 1:
            continue
        poly = shape(geom)
        if poly.is_empty:
            continue
        geoms.append(poly)

    if not geoms:
        gdf = gpd.GeoDataFrame(columns=["geometry", "area_m2", "perimeter_m", "mean_score"], geometry="geometry", crs=mask.rio.crs)
        gdf.to_file(out_geojson, driver="GeoJSON")
        _write_kml(gdf.to_crs(4326), out_kml)
        return gdf

    gdf = gpd.GeoDataFrame(geometry=geoms, crs=mask.rio.crs)
    gdf["area_m2"] = gdf.geometry.area
    gdf["perimeter_m"] = gdf.geometry.length

    if min_area_m2 > 0:
        gdf = gdf[gdf["area_m2"] >= min_area_m2]

    if gdf.empty:
        empty = gpd.GeoDataFrame(columns=gdf.columns, geometry="geometry", crs=mask.rio.crs)
        empty.to_file(out_geojson, driver="GeoJSON")
        _write_kml(empty.to_crs(4326), out_kml)
        return empty

    # Compute mean change score per polygon using raster sampling.
    score_np = score.values
    mean_scores: list[float] = []
    for poly in gdf.geometry:
        mask_raster = features.rasterize(
            [(poly, 1)],
            out_shape=mask_np.shape,
            transform=transform,
            fill=0,
            dtype="uint8",
        )
        values = score_np[mask_raster == 1]
        mean_scores.append(float(np.nanmean(values)) if values.size else float("nan"))
    gdf["mean_score"] = mean_scores

    gdf = gdf.explode(index_parts=False).reset_index(drop=True)
    gdf = gdf.to_crs(4326)

    gdf.to_file(out_geojson, driver="GeoJSON")
    _write_kml(gdf, out_kml)
    return gdf


def _write_kml(gdf: gpd.GeoDataFrame, path: Path) -> None:
    from simplekml import Kml

    path.parent.mkdir(parents=True, exist_ok=True)
    kml = Kml()
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            polys = [geom]
        elif geom.geom_type == "MultiPolygon":
            polys = list(geom.geoms)
        else:
            continue
        for poly in polys:
            coords = [(lon, lat, 0.0) for lon, lat in poly.exterior.coords]
            pol = kml.newpolygon(name="Change Polygon", outerboundaryis=coords)
            pol.style.polystyle.color = "7dff7800"  # semi-transparent orange
            pol.style.linestyle.color = "ff783d00"
            pol.style.linestyle.width = 2
            description = json.dumps(
                {
                    "area_m2": row.get("area_m2"),
                    "perimeter_m": row.get("perimeter_m"),
                    "mean_score": row.get("mean_score"),
                },
                indent=2,
            )
            pol.description = description
    kml.save(path)
