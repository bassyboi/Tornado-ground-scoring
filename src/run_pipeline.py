"""Entry-point for the thin-slice Sentinel-2 change detection pipeline."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict

from . import aoi, change_s2, config, threshold, vectorize, webmap
from .stac_fetch import fetch_s2_mosaic

LOGGER = logging.getLogger(__name__)


def run_pipeline(config_path: Path) -> Dict[str, Any]:
    """Execute the Sentinel-2 change workflow described in the README."""

    pipeline_config = config.PipelineConfig.load(config_path)
    LOGGER.info("Loaded configuration from %s", config_path)

    projected_geom, geographic_geom, _, _ = aoi.load_aoi(
        pipeline_config.aoi_path, pipeline_config.stack_epsg
    )
    bounds_projected = projected_geom.bounds

    LOGGER.info("Fetching pre-event Sentinel-2 mosaic")
    pre_result = fetch_s2_mosaic(
        pipeline_config.stac.catalogs,
        pipeline_config.pre_from,
        pipeline_config.pre_to,
        geographic_geom,
        bounds_projected,
        pipeline_config.stack_epsg,
    )

    LOGGER.info("Fetching post-event Sentinel-2 mosaic")
    post_result = fetch_s2_mosaic(
        pipeline_config.stac.catalogs,
        pipeline_config.post_from,
        pipeline_config.post_to,
        geographic_geom,
        bounds_projected,
        pipeline_config.stack_epsg,
    )

    artifacts_dir = Path("artifacts")
    score, change_products = change_s2.compute_change_score(
        pre_result.data, post_result.data, pipeline_config.weights, artifacts_dir
    )

    mask_path = artifacts_dir / "change_mask.tif"
    mask = threshold.apply_threshold(score, pipeline_config.threshold, mask_path)

    docs_data_dir = Path("docs") / "data"
    geojson_path = docs_data_dir / "changes.geojson"
    kml_path = docs_data_dir / "changes.kml"

    gdf = vectorize.mask_to_polygons(mask, score, geojson_path, kml_path)

    meta_path = docs_data_dir / "meta.json"
    mosaic_metadata = {
        "pre": {"catalog": pre_result.catalog_url, "items": pre_result.item_count},
        "post": {"catalog": post_result.catalog_url, "items": post_result.item_count},
    }
    webmap.write_metadata(pipeline_config, mosaic_metadata, len(gdf), meta_path)

    index_path = Path("docs") / "index.html"
    webmap.ensure_leaflet_page(index_path, pipeline_config)

    outputs = {
        "pre_mosaic": change_products.pre_mosaic_path,
        "post_mosaic": change_products.post_mosaic_path,
        "score_raster": change_products.score_path,
        "mask_raster": mask_path,
        "geojson": geojson_path,
        "kml": kml_path,
        "meta": meta_path,
        "index": index_path,
    }
    LOGGER.info("Pipeline complete")
    for key, value in outputs.items():
        LOGGER.info("%s -> %s", key, value)
    return outputs


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        dest="config_path",
        type=Path,
        required=True,
        help="Path to the YAML configuration file",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = _parse_args(argv or sys.argv[1:])
    try:
        run_pipeline(args.config_path)
    except Exception as exc:  # pragma: no cover - runtime errors
        LOGGER.error("Pipeline failed: %s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
