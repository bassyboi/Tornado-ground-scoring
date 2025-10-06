"""End-to-end pipeline CLI."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import click
import geopandas as gpd
import numpy as np
import xarray as xr
import yaml
from rasterio import features
from shapely.geometry import mapping

from . import change_core, fetch_stack, gss
from .aoi_utils import load_aoi

LOGGER = logging.getLogger(__name__)

S2_BANDS = ["nir", "red", "green", "blue", "scl"]
S1_BANDS = ["vv", "vh"]


@click.command()
@click.option("--config", "config_path", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--export-csv", "export_csv", type=click.Path(dir_okay=False), default=None, help="Optional CSV output of polygon stats")
def main(config_path: str, export_csv: Optional[str]) -> None:
    """Run the ground-scour mapping pipeline."""

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    cfg = _load_config(Path(config_path))

    aoi = load_aoi(cfg["aoi_geojson"])
    catalogs = cfg.get("stac_catalogs", [])

    pre_items = fetch_stack.search_sentinel2(aoi, cfg["pre_from"], cfg["pre_to"], catalogs)
    post_items = fetch_stack.search_sentinel2(aoi, cfg["post_from"], cfg["post_to"], catalogs)

    pre_mosaic = fetch_stack.mosaic_s2(pre_items, S2_BANDS, aoi)
    post_mosaic = fetch_stack.mosaic_s2(post_items, S2_BANDS, aoi)

    score = change_core.change_score_s2(pre_mosaic, post_mosaic, cfg.get("weights", {}))
    score = score.rio.write_crs(pre_mosaic.rio.crs)
    score = score.rio.write_transform(pre_mosaic.rio.transform())

    if cfg.get("use_sentinel1_grd"):
        s1_pre_items = fetch_stack.search_sentinel1(aoi, cfg["pre_from"], cfg["pre_to"], catalogs)
        s1_post_items = fetch_stack.search_sentinel1(aoi, cfg["post_from"], cfg["post_to"], catalogs)
        s1_pre = fetch_stack.mosaic_s1(s1_pre_items, S1_BANDS, aoi)
        s1_post = fetch_stack.mosaic_s1(s1_post_items, S1_BANDS, aoi)
        score = change_core.add_s1_logratio(score, s1_pre, s1_post, cfg.get("weights", {}).get("s1_logratio", 0.1))
        score = score.rio.write_crs(pre_mosaic.rio.crs)
        score = score.rio.write_transform(pre_mosaic.rio.transform())

    score_path = Path("artifacts/change_score.tif")
    _ensure_directory(score_path.parent)
    _save_dataarray(score.astype(np.float32), score_path)

    threshold_cfg = cfg.get("threshold", "otsu")
    if isinstance(threshold_cfg, str):
        mask = change_core.threshold_score(score, method=threshold_cfg)
    else:
        mask = change_core.threshold_score(score, method="numeric", numeric=float(threshold_cfg))
    mask = mask.rio.write_crs(score.rio.crs)
    mask = mask.rio.write_transform(score.rio.transform())

    mask_path = Path("artifacts/change_mask.tif")
    _save_dataarray(mask.astype(np.uint8), mask_path)

    polygons = change_core.polygons_from_mask(mask, cfg.get("min_blob_area_m2", 500.0))
    polygons = polygons.reset_index(drop=True)

    if not polygons.empty:
        polygons = _attach_polygon_stats(polygons, score)
        polygons["bearing_deg"] = polygons.geometry.apply(gss.major_axis_orientation)
        polygons["elongation_ratio"] = polygons.geometry.apply(gss.elongation_ratio)
        if cfg.get("elongation_filter", False):
            tol = cfg.get("elongation_tolerance_deg", 30.0)
            min_ratio = cfg.get("elongation_min_ratio", 2.0)
            polygons = gss.apply_elongation_filter(
                polygons,
                preferred_bearing=cfg.get("preferred_track_bearing_deg", 0.0),
                ang_tol=tol,
                min_ratio=min_ratio,
            )
            polygons = polygons.reset_index(drop=True)
        gss_breaks = cfg.get("gss", {}).get("breaks", "quantile")
        polygons["gss"] = gss.score_to_gss(polygons["mean_score"].values, breaks=gss_breaks)
    else:
        polygons["mean_score"] = []
        polygons["max_score"] = []
        polygons["gss"] = []
        polygons["bearing_deg"] = []
        polygons["elongation_ratio"] = []

    geojson_path = Path("docs/data/changes.geojson")
    _ensure_directory(geojson_path.parent)
    _write_geojson(polygons, geojson_path)

    if export_csv and not polygons.empty:
        _ensure_directory(Path(export_csv).parent or Path("."))
        polygons.drop(columns="geometry").to_csv(export_csv, index=False)

    _write_leaflet_page(Path("docs/index.html"), cfg.get("web", {}))

    click.echo("Done. Outputs written to artifacts/ and docs/.")


def _load_config(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config must be a mapping")
    return cfg


def _ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _save_dataarray(arr: xr.DataArray, path: Path) -> None:
    if np.issubdtype(arr.dtype, np.integer):
        arr = arr.rio.write_nodata(0)
    else:
        arr = arr.rio.write_nodata(np.nan)
    arr.rio.to_raster(path)
    LOGGER.info("Wrote %s", path)


def _attach_polygon_stats(polygons: gpd.GeoDataFrame, score: xr.DataArray) -> gpd.GeoDataFrame:
    data = score.values
    transform = score.rio.transform()
    crs = score.rio.crs
    means = []
    maxs = []
    for geom in polygons.geometry.to_crs(crs):
        mask = features.geometry_mask([mapping(geom)], out_shape=data.shape, transform=transform, invert=True)
        masked = np.where(mask, data, np.nan)
        mean_raw = np.nanmean(masked)
        max_raw = np.nanmax(masked)
        mean_val = float(mean_raw) if np.isfinite(mean_raw) else 0.0
        max_val = float(max_raw) if np.isfinite(max_raw) else 0.0
        means.append(mean_val)
        maxs.append(max_val)
    polygons = polygons.copy()
    polygons["mean_score"] = means
    polygons["max_score"] = maxs
    polygons["area_m2"] = polygons.to_crs(3857).area
    return polygons


def _write_geojson(polygons: gpd.GeoDataFrame, path: Path) -> None:
    if polygons.empty:
        collection = {"type": "FeatureCollection", "features": []}
    else:
        collection = json.loads(polygons.to_json())
    with path.open("w", encoding="utf-8") as f:
        json.dump(collection, f)
    LOGGER.info("Wrote %s", path)


def _write_leaflet_page(path: Path, web_cfg: Dict[str, Any]) -> None:
    title = web_cfg.get("title", "Ground Scour Map")
    description = web_cfg.get("description", "Auto-built from Sentinel pre/post windows.")
    center = web_cfg.get("center", [0, 0])
    zoom = web_cfg.get("zoom", 6)

    html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>{title}</title>
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <link rel=\"stylesheet\" href=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.css\" integrity=\"sha512-sA+q5ms5F2yZefRg2fzGEylh4G6dkprdFMn/hTyBC0bY4Z1cdq9VHtV6nRvWmvMRGsiE232zfoFM5x3Dm7G8eA==\" crossorigin=\"\" />
  <link rel=\"stylesheet\" href=\"styles.css\" />
</head>
<body>
  <header>
    <h1>{title}</h1>
    <p>{description}</p>
  </header>
  <div id=\"map\"></div>
  <div id=\"legend\"></div>
  <script src=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.js\" integrity=\"sha512-o9N1j7kG6r0s+P0LlwM651op01qmPvvrLpzjAU6Rz6U02zszbmWzxubUANkqG0x74pJ1gzhS+M4LMEnM08JdKw==\" crossorigin=\"\"></script>
  <script>
    const map = L.map('map').setView([{center[1]}, {center[0]}], {zoom});
    L.tileLayer('https://tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
      maxZoom: 19,
      attribution: '&copy; OpenStreetMap contributors'
    }}).addTo(map);

    const colors = {{
      0: '#f7fbff',
      1: '#c6dbef',
      2: '#9ecae1',
      3: '#6baed6',
      4: '#3182bd',
      5: '#08519c'
    }};

    function style(feature) {{
      const gss = feature.properties?.gss ?? 0;
      return {{
        color: '#333',
        weight: 1,
        fillColor: colors[gss] || '#cccccc',
        fillOpacity: 0.6
      }};
    }}

    fetch('data/changes.geojson')
      .then(resp => resp.json())
      .then(data => {{
        const layer = L.geoJSON(data, {{ style }}).addTo(map);
        if (layer.getLayers().length) {{
          map.fitBounds(layer.getBounds());
        }}
        buildLegend();
      }})
      .catch(err => {{
        console.error('Failed to load GeoJSON', err);
        buildLegend();
      }});

    function buildLegend() {{
      const legend = document.getElementById('legend');
      legend.innerHTML = '<h3>Ground-Scour Score</h3>';
      for (let i = 0; i <= 5; i++) {{
        const row = document.createElement('div');
        row.className = 'legend-row';
        row.innerHTML = `<span class="swatch" style="background:${{colors[i]}}"></span> GSS ${'{'}i{'}'}`;
        legend.appendChild(row);
      }}
    }}
  </script>
</body>
</html>
"""

    with path.open("w", encoding="utf-8") as f:
        f.write(html)
    LOGGER.info("Wrote %s", path)


if __name__ == "__main__":
    main()
