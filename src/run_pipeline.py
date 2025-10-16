"""End-to-end pipeline CLI."""
from __future__ import annotations

import hashlib
import json
import logging
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Sequence
import xml.etree.ElementTree as ET
from xml.dom import minidom

import click
import geopandas as gpd
import numpy as np
import xarray as xr
import yaml

from . import change_core, fetch_stack, features, scoring, storm_filter, storm_reports
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
    cfg_path = Path(config_path)
    cfg = _load_config(cfg_path)
    run_started = datetime.utcnow().replace(microsecond=0)

    aoi = load_aoi(cfg["aoi_geojson"])
    catalogs = cfg.get("stac_catalogs", [])

    search_cfg = cfg.get("search", {})
    expand_days = int(search_cfg.get("auto_expand_max_days", 0) or 0)
    expand_step = int(search_cfg.get("auto_expand_step_days", 2) or 1)

    stack_epsg = cfg.get("stack_epsg")
    stack_epsg_int = int(stack_epsg) if stack_epsg is not None else None

    post_from = _parse_date(cfg["post_from"])
    post_to = _parse_date(cfg["post_to"])

    storm_cfg = cfg.get("storm_filter", {})
    if storm_cfg.get("enabled"):
        schema = storm_filter.CatalogSchema(
            datetime_column=storm_cfg.get("datetime_column", "event_time_utc"),
            hazard_column=storm_cfg.get("hazard_column", "hazard"),
            latitude_column=storm_cfg.get("latitude_column", "latitude"),
            longitude_column=storm_cfg.get("longitude_column", "longitude"),
        )
        catalog = None
        scrape_cfg = storm_cfg.get("scrape")
        if scrape_cfg:
            provider = (scrape_cfg.get("provider") or "iem_lsr").strip().lower()
            hazard_overrides = scrape_cfg.get("hazards")
            hazard_list = (
                hazard_overrides
                if hazard_overrides is not None
                else storm_cfg.get("hazards")
            )
            if hazard_list:
                hazard_list = [str(value) for value in hazard_list]
            lookback_days = int(scrape_cfg.get("lookback_days", 0) or 0)
            lookahead_days = int(scrape_cfg.get("lookahead_days", 0) or 0)
            bbox_buffer_km = float(scrape_cfg.get("bbox_buffer_km", 0.0) or 0.0)

            scrape_start_date = post_from - timedelta(days=lookback_days)
            scrape_end_date = post_to + timedelta(days=lookahead_days)
            start_dt = datetime.combine(
                scrape_start_date, datetime.min.time()
            ).replace(tzinfo=timezone.utc)
            end_dt = datetime.combine(
                scrape_end_date + timedelta(days=1), datetime.min.time()
            ).replace(tzinfo=timezone.utc)
            bounds = storm_reports.buffered_bounds(aoi, buffer_km=bbox_buffer_km)

            def _coerce_iterable(raw):
                if raw is None:
                    return None
                if isinstance(raw, (list, tuple, set)):
                    return [str(item) for item in raw]
                return [str(raw)]

            states = (
                _coerce_iterable(scrape_cfg.get("states"))
                if provider == "bom_warnings"
                else None
            )
            phenomena = (
                _coerce_iterable(scrape_cfg.get("phenomena"))
                if provider == "bom_warnings"
                else None
            )

            if provider == "iem_lsr":
                catalog = storm_reports.fetch_iem_local_storm_reports(
                    start=start_dt,
                    end=end_dt,
                    bounds=bounds,
                    hazards=hazard_list,
                )
            elif provider == "bom_warnings":
                catalog = storm_reports.fetch_bom_warnings(
                    start=start_dt,
                    end=end_dt,
                    bounds=bounds,
                    hazards=hazard_list,
                    states=states,
                    phenomena=phenomena,
                )
            else:
                raise ValueError(
                    "Unsupported storm_filter.scrape.provider "
                    f"'{provider}' (expected 'iem_lsr' or 'bom_warnings')"
                )

            LOGGER.info(
                "Scraped %s storm reports between %s and %s using %s",
                len(catalog),
                scrape_start_date,
                scrape_end_date,
                provider,
            )
            export_csv = scrape_cfg.get("export_csv")
            if export_csv:
                storm_reports.write_catalog_csv(catalog, Path(export_csv))
            metadata_path = scrape_cfg.get("metadata_path")
            if metadata_path:
                metadata = {
                    "provider": provider,
                    "start": start_dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "end": end_dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    "bounds": list(bounds),
                    "hazards": list(hazard_list or []),
                    "count": int(len(catalog)),
                }
                if provider == "bom_warnings":
                    metadata.update(
                        {
                            "states": list(states) if states else [],
                            "phenomena": list(phenomena) if phenomena else [],
                        }
                    )
                storm_reports.save_metadata(metadata, Path(metadata_path))
            if catalog.empty:
                LOGGER.info(
                    "Storm report scrape returned no events for the requested window"
                )
        if catalog is None:
            catalog_path = storm_cfg.get("catalog")
            if not catalog_path:
                raise KeyError(
                    "storm_filter.scrape is enabled but no catalog was scraped. "
                    "Provide a live scrape configuration or an explicit catalog fallback."
                )
            catalog = storm_filter.load_catalog(Path(catalog_path), schema)
        backfill_dirs = storm_cfg.get("auto_backfill_directions")
        selection = storm_filter.resolve_event_window(
            catalog=catalog,
            aoi=aoi,
            window_start=post_from,
            window_end=post_to,
            hazards=storm_cfg.get("hazards"),
            days_before=int(storm_cfg.get("days_before", 0)),
            days_after=int(storm_cfg.get("days_after", 0)),
            distance_km=float(storm_cfg.get("distance_km", 0.0)),
            auto_backfill_max_days=storm_cfg.get("auto_backfill_max_days", 0),
            auto_backfill_step_days=int(storm_cfg.get("auto_backfill_step_days", 1) or 1),
            auto_backfill_directions=backfill_dirs,
        )
        events = selection.events
        post_from = selection.window_start
        post_to = selection.window_end
        if selection.shift_days:
            LOGGER.info(
                "Storm filter auto-adjusted post window by %+d days to %s–%s",
                selection.shift_days,
                post_from,
                post_to,
            )
        export_path = storm_cfg.get("export_geojson")
        if export_path:
            storm_filter.export_events(events, Path(export_path))
        if events.empty:
            LOGGER.info(
                "Storm filter found no qualifying events between %s and %s; skipping",
                selection.window_start,
                selection.window_end,
            )
            palette = _resolve_palette(cfg.get("webmap", {}).get("palette"))
            _write_empty_outputs(
                Path("docs/data/changes.geojson"),
                Path("docs/data/changes.kml"),
                palette,
            )
            empty_summary = {
                "generated_at": datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
                "config_hash": _config_hash(cfg_path),
                "title": cfg.get("web", {}).get("title", "Ground Scour Map"),
                "description": cfg.get("web", {}).get(
                    "description", "Auto-built from Sentinel pre/post windows."
                ),
                "timestamps": {
                    "started": run_started.isoformat() + "Z",
                    "finished": datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
                },
                "total": 0,
                "total_input": 0,
                "per_gss": {str(level): 0 for level in range(6)},
                "bearing_cluster_deg": None,
                "mean_elongation": None,
                "filters": {},
                "score_rules": {},
                "webmap_palette": palette,
            }
            _write_summary(Path("docs/data/summary.json"), empty_summary)
            click.echo("No storm events in window; outputs cleared.")
            return
        LOGGER.info("Storm filter retained %s events", len(events))

    LOGGER.info(
        "Searching Sentinel-2 %s imagery between %s and %s",
        "pre-event",
        cfg["pre_from"],
        cfg["pre_to"],
    )
    pre_items, pre_from_str, pre_to_str = fetch_stack.search_sentinel2_with_fallback(
        aoi,
        cfg["pre_from"],
        cfg["pre_to"],
        catalogs,
        max_expansion_days=expand_days,
        step_days=expand_step,
        window_label="pre-event",
    )
    LOGGER.info(
        "Searching Sentinel-2 %s imagery between %s and %s",
        "post-event",
        post_from.isoformat(),
        post_to.isoformat(),
    )
    post_items, post_from_str, post_to_str = fetch_stack.search_sentinel2_with_fallback(
        aoi,
        post_from.isoformat(),
        post_to.isoformat(),
        catalogs,
        max_expansion_days=expand_days,
        step_days=expand_step,
        window_label="post-event",
    )
    LOGGER.info("Building Sentinel-2 mosaics")
    pre_mosaic = fetch_stack.mosaic_s2(pre_items, S2_BANDS, aoi, target_epsg=stack_epsg_int)
    post_mosaic = fetch_stack.mosaic_s2(post_items, S2_BANDS, aoi, target_epsg=stack_epsg_int)

    pre_path = Path("artifacts/pre_mosaic.tif")
    post_path = Path("artifacts/post_mosaic.tif")
    _ensure_directory(pre_path.parent)
    LOGGER.info("Writing pre-event mosaic to %s", pre_path)
    _save_dataarray(pre_mosaic.astype(np.float32), pre_path)
    LOGGER.info("Writing post-event mosaic to %s", post_path)
    _save_dataarray(post_mosaic.astype(np.float32), post_path)

    LOGGER.info("Computing change score")
    score = change_core.change_score_s2(pre_mosaic, post_mosaic, cfg.get("weights", {}))
    score = score.rio.write_crs(pre_mosaic.rio.crs)
    score = score.rio.write_transform(pre_mosaic.rio.transform())

    if cfg.get("use_sentinel1_grd"):
        LOGGER.info("Searching Sentinel-1 GRD collections")
        s1_pre_items, _, _ = fetch_stack.search_sentinel1_with_fallback(
            aoi,
            pre_from_str,
            pre_to_str,
            catalogs,
            max_expansion_days=expand_days,
            step_days=expand_step,
            window_label="pre-event",
        )
        s1_post_items, _, _ = fetch_stack.search_sentinel1_with_fallback(
            aoi,
            post_from_str,
            post_to_str,
            catalogs,
            max_expansion_days=expand_days,
            step_days=expand_step,
            window_label="post-event",
        )
        s1_pre = fetch_stack.mosaic_s1(s1_pre_items, S1_BANDS, aoi, target_epsg=stack_epsg_int)
        s1_post = fetch_stack.mosaic_s1(s1_post_items, S1_BANDS, aoi, target_epsg=stack_epsg_int)
        score = change_core.add_s1_logratio(score, s1_pre, s1_post, cfg.get("weights", {}).get("s1_logratio", 0.1))
        score = score.rio.write_crs(pre_mosaic.rio.crs)
        score = score.rio.write_transform(pre_mosaic.rio.transform())

    score_path = Path("artifacts/change_score.tif")
    _ensure_directory(score_path.parent)
    LOGGER.info("Writing change score raster to %s", score_path)
    _save_dataarray(score.astype(np.float32), score_path)

    threshold_cfg = cfg.get("threshold", "otsu")
    if isinstance(threshold_cfg, str):
        threshold_method = threshold_cfg
        threshold_value = change_core.resolve_threshold_value(score, method=threshold_method)
        mask = change_core.threshold_score(score, method=threshold_method)
    else:
        threshold_method = "numeric"
        numeric_value = float(threshold_cfg)
        threshold_value = change_core.resolve_threshold_value(
            score, method=threshold_method, numeric=numeric_value
        )
        mask = change_core.threshold_score(score, method=threshold_method, numeric=numeric_value)
    mask = mask.rio.write_crs(score.rio.crs)
    mask = mask.rio.write_transform(score.rio.transform())

    mask_path = Path("artifacts/change_mask.tif")
    LOGGER.info("Writing change mask raster to %s", mask_path)
    _save_dataarray(mask.astype(np.uint8), mask_path)

    LOGGER.info("Extracting change polygons")
    polygons = change_core.polygons_from_mask(mask, cfg.get("min_blob_area_m2", 500.0))
    polygons = polygons.reset_index(drop=True)

    filters_cfg = cfg.get("filters", {})
    webmap_cfg = cfg.get("webmap", {})
    palette = _resolve_palette(webmap_cfg.get("palette"))

    target_epsg = stack_epsg_int or pre_mosaic.rio.crs.to_epsg()
    if target_epsg is None:
        raise ValueError("Unable to resolve target EPSG for feature extraction")

    min_hole_area = float(filters_cfg.get("min_hole_area", 0.0) or 0.0)
    LOGGER.info("Computing polygon features")
    features_gdf = features.extract(
        polygons,
        target_epsg=int(target_epsg),
        pre_path=pre_path,
        post_path=post_path,
        change_path=score_path,
        min_hole_area=min_hole_area,
    )

    LOGGER.info("Applying scoring rules")
    scoring_result = scoring.apply(
        features_gdf,
        filters_cfg=filters_cfg,
        rules_cfg=cfg.get("score_rules", {}),
    )
    scored_polygons = scoring_result.polygons
    summary = scoring_result.summary
    summary.setdefault("filters", {})["min_hole_area"] = float(min_hole_area)
    summary.setdefault("webmap_palette", list(palette))
    summary["threshold_method"] = threshold_method
    summary["threshold_value"] = float(threshold_value) if np.isfinite(threshold_value) else None

    features_path = Path("docs/data/changes_features.geojson")
    scored_path = Path("docs/data/changes_scored.geojson")
    legacy_path = Path("docs/data/changes.geojson")
    _ensure_directory(features_path.parent)
    LOGGER.info("Writing polygon features to %s", features_path)
    _write_geojson(features_gdf, features_path)
    LOGGER.info("Writing scored polygons to %s", scored_path)
    _write_geojson(scored_polygons, scored_path)
    # Maintain legacy output name for compatibility.
    _write_geojson(scored_polygons, legacy_path)
    _write_kml(scored_polygons, Path("docs/data/changes.kml"), palette)

    if export_csv and not scored_polygons.empty:
        _ensure_directory(Path(export_csv).parent or Path("."))
        scored_polygons.drop(columns="geometry").to_csv(export_csv, index=False)

    run_finished = datetime.utcnow().replace(microsecond=0)
    summary["generated_at"] = run_finished.isoformat() + "Z"
    summary["config_hash"] = _config_hash(cfg_path)
    summary["title"] = cfg.get("web", {}).get("title", "Ground Scour Map")
    summary["description"] = cfg.get("web", {}).get(
        "description", "Auto-built from Sentinel pre/post windows."
    )
    summary["timestamps"] = {
        "started": run_started.isoformat() + "Z",
        "finished": run_finished.isoformat() + "Z",
    }
    summary_path = Path("docs/data/summary.json")
    LOGGER.info("Writing summary to %s", summary_path)
    _write_summary(summary_path, summary)

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


def _write_summary(path: Path, summary: Dict[str, Any]) -> None:
    _ensure_directory(path.parent)
    with path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)


def _config_hash(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hasher.update(chunk)
    return hasher.hexdigest()[:16]


def _resolve_palette(raw: Any) -> list[str]:
    default = ["#d0d0d0", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594"]
    if isinstance(raw, (list, tuple)):
        palette = [str(color) for color in raw if color]
    else:
        palette = []
    if not palette:
        palette = default.copy()
    while len(palette) < 6:
        palette.append(default[len(palette)])
    return palette[:6]


def _save_dataarray(arr: xr.DataArray, path: Path) -> None:
    if np.issubdtype(arr.dtype, np.integer):
        arr = arr.rio.write_nodata(0)
    else:
        arr = arr.rio.write_nodata(np.nan)
    arr.rio.to_raster(path)
    LOGGER.info("Wrote %s", path)


def _write_geojson(polygons: gpd.GeoDataFrame, path: Path) -> None:
    if polygons.empty:
        collection = {"type": "FeatureCollection", "features": []}
    else:
        collection = json.loads(polygons.to_json())
    with path.open("w", encoding="utf-8") as f:
        json.dump(collection, f)
    LOGGER.info("Wrote %s", path)


def _write_kml(polygons: gpd.GeoDataFrame, path: Path, palette: Sequence[str]) -> None:
    colors = {idx: palette[idx] if idx < len(palette) else palette[-1] for idx in range(6)}

    kml = ET.Element("kml", xmlns="http://www.opengis.net/kml/2.2")
    document = ET.SubElement(kml, "Document")
    ET.SubElement(document, "name").text = "Ground Scour Change"

    for gss_value, hex_color in colors.items():
        style = ET.SubElement(document, "Style", id=f"gss-{gss_value}")
        line_style = ET.SubElement(style, "LineStyle")
        ET.SubElement(line_style, "color").text = "ff333333"
        ET.SubElement(line_style, "width").text = "1"
        poly_style = ET.SubElement(style, "PolyStyle")
        ET.SubElement(poly_style, "color").text = _hex_to_kml_color(hex_color, 0.6)
        ET.SubElement(poly_style, "outline").text = "1"

    if not polygons.empty:
        polygons_wgs84 = polygons.to_crs(4326)
        for _, row in polygons_wgs84.iterrows():
            placemark = ET.SubElement(document, "Placemark")
            gss_value = int(row.get("gss", 0))
            ET.SubElement(placemark, "name").text = f"GSS {gss_value}"
            description_parts = []
            for label, key, fmt in [
                ("Change score", "change_score_mean", "{:.3f}"),
                ("ΔNDVI", "delta_ndvi_mean", "{:.3f}"),
                ("Δbrightness", "delta_brightness_mean", "{:.3f}"),
                ("Elongation", "elongation", "{:.2f}×"),
                ("Length", "length_m", "{:.0f} m"),
                ("Width", "width_m", "{:.0f} m"),
                ("Area", "area_m2", "{:.2f} m²"),
            ]:
                value = row.get(key)
                if value is None:
                    continue
                try:
                    numeric = float(value)
                except (TypeError, ValueError):
                    continue
                if np.isfinite(numeric):
                    description_parts.append(f"{label}: {fmt.format(numeric)}")
            explain = row.get("gss_explain")
            if explain:
                description_parts.append(str(explain))
            if description_parts:
                ET.SubElement(placemark, "description").text = "\n".join(description_parts)
            ET.SubElement(placemark, "styleUrl").text = f"#gss-{gss_value}"
            _append_geometry(placemark, row.geometry)

    path.parent.mkdir(parents=True, exist_ok=True)
    xml_bytes = ET.tostring(kml, encoding="utf-8")
    pretty = minidom.parseString(xml_bytes).toprettyxml(indent="  ")
    path.write_text(pretty, encoding="utf-8")
    LOGGER.info("Wrote %s", path)


def _append_geometry(parent: ET.Element, geometry) -> None:
    if geometry.is_empty:
        return
    geom_type = geometry.geom_type
    if geom_type == "Polygon":
        _append_polygon(parent, geometry)
    elif geom_type == "MultiPolygon":
        multi = ET.SubElement(parent, "MultiGeometry")
        for geom in geometry.geoms:
            _append_polygon(multi, geom)
    else:
        # Fallback to centroid representation for unexpected geometries
        point = geometry.centroid
        point_elem = ET.SubElement(parent, "Point")
        ET.SubElement(point_elem, "coordinates").text = _coords_to_string(point.coords)


def _append_polygon(parent: ET.Element, polygon) -> None:
    poly_elem = ET.SubElement(parent, "Polygon")
    outer = ET.SubElement(poly_elem, "outerBoundaryIs")
    outer_ring = ET.SubElement(outer, "LinearRing")
    ET.SubElement(outer_ring, "coordinates").text = _coords_to_string(polygon.exterior.coords)
    for interior in polygon.interiors:
        inner = ET.SubElement(poly_elem, "innerBoundaryIs")
        inner_ring = ET.SubElement(inner, "LinearRing")
        ET.SubElement(inner_ring, "coordinates").text = _coords_to_string(interior.coords)


def _coords_to_string(coords) -> str:
    return " ".join(f"{x:.6f},{y:.6f},0" for x, y, *_ in coords)


def _hex_to_kml_color(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) != 6:
        raise ValueError("Hex colors must be 6 digits")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    a = max(0, min(255, int(round(alpha * 255))))
    return f"{a:02x}{b:02x}{g:02x}{r:02x}"


def _write_leaflet_page(path: Path, web_cfg: Dict[str, Any]) -> None:
    if path.exists():
        LOGGER.info("Web map already present at %s; skipping overwrite", path)
        return

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
  <link rel=\"stylesheet\" href=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.css\" />
  <style>
    body {{ margin: 0; font-family: sans-serif; }}
    #map {{ height: 90vh; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <p>{description}</p>
  <div id=\"map\"></div>
  <script src=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.js\"></script>
  <script>
    const map = L.map('map').setView([{center[1]}, {center[0]}], {zoom});
    L.tileLayer('https://tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
      maxZoom: 19,
      attribution: '&copy; OpenStreetMap contributors'
    }}).addTo(map);
    fetch('data/changes_scored.geojson')
      .then(resp => resp.json())
      .then(data => {{
        const layer = L.geoJSON(data).addTo(map);
        if (layer.getLayers().length) {{
          map.fitBounds(layer.getBounds());
        }}
      }});
  </script>
</body>
</html>
"""

    with path.open("w", encoding="utf-8") as f:
        f.write(html)
    LOGGER.info("Created default web map at %s", path)


def _write_empty_outputs(
    geojson_path: Path, kml_path: Path, palette: Sequence[str] | None = None
) -> None:
    empty = gpd.GeoDataFrame(geometry=[], crs=4326)
    _ensure_directory(geojson_path.parent)
    palette_list = _resolve_palette(palette) if palette is not None else _resolve_palette(None)
    _write_geojson(empty, geojson_path)
    _write_geojson(empty, geojson_path.parent / "changes_scored.geojson")
    _write_geojson(empty, geojson_path.parent / "changes_features.geojson")
    _write_kml(empty, kml_path, palette_list)


def _parse_date(value: str) -> date:
    try:
        return datetime.strptime(value, "%Y-%m-%d").date()
    except ValueError as exc:
        raise ValueError(f"Dates must be YYYY-MM-DD strings: {value}") from exc


if __name__ == "__main__":
    main()
