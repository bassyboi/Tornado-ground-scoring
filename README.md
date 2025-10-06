# Tornado Ground-Scour Mapper

Detect tornado-driven ground change from Sentinel imagery, score the severity, and publish an interactive web map — all from a single Python pipeline.

## Features
- Fetch Sentinel-2 (and optional Sentinel-1 RTC) imagery through public STAC APIs.
- Build median pre/post mosaics and compute a composite change score using vegetation loss, brightness increase, and SAR backscatter change.
- Threshold the change map, filter polygons by elongation, and assign a Ground-Scour Score (GSS 0–5).
- Export GeoTIFF rasters, GeoJSON polygons, and a ready-to-host Leaflet map under `docs/` for GitHub Pages.
- Emit a Google Earth-ready KML file alongside the GeoJSON export for quick sharing.
- Continuous integration workflow that runs the pipeline and publishes run artifacts on every push.

## Quick start
1. **Install dependencies** (Python 3.11):
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```
2. **Configure the analysis** by editing `data/config.yaml`. Update the AOI path, date windows, weights, threshold, and web map metadata as needed.
3. **Run the pipeline**:
   ```bash
   python -m src.run_pipeline --config data/config.yaml
   ```
   The run will create GeoTIFFs in `artifacts/` and a `docs/data/changes.geojson` file consumed by the Leaflet app at `docs/index.html`.
   Add `--export-csv outputs.csv` to capture polygon statistics as a table.
4. **Serve the map locally** (optional):
   ```bash
   python -m http.server --directory docs 8000
   ```
   Navigate to `http://localhost:8000` to view the map.

## Configuration reference
See `data/config.yaml` for the default options:
- **STAC search**: Provide one or more public catalog URLs. Sentinel-2 L2A is required; Sentinel-1 RTC is optional.
- **Weights**: Tune the component weights to emphasize vegetation loss (ΔNDVI), brightness increase, or SAR log-ratio response.
- **Threshold**: Use `"otsu"` for automatic selection or set a numeric cutoff between 0 and 1.
- **Elongation filter**: Enable to keep polygons aligned with an expected tornado-track bearing and an elongation ratio ≥ 2.0. Optional keys `elongation_tolerance_deg` and `elongation_min_ratio` further tune the filter.
- **GSS breaks**: Either `"quantile"` (default quintiles) or an explicit list of five monotonically increasing thresholds.
- **Web map**: Customize the title, description, and initial map viewpoint.
- **Storm filter**: Optionally gate the entire run on a catalog of historical storm reports. Provide a CSV file with at least
  the event timestamp, hazard type, and coordinates. When enabled, the pipeline loads the catalog, keeps events that intersect
  the AOI (with an optional distance buffer) within the specified post-event window, exports them to GeoJSON, and aborts early
  when no qualifying storm days are found.

### Storm-filter workflow

The `storm_filter` block in `data/config.yaml` demonstrates the expected schema:

```yaml
storm_filter:
  enabled: true
  catalog: data/storm_reports_example.csv
  datetime_column: event_time_utc
  hazard_column: hazard
  latitude_column: latitude
  longitude_column: longitude
  hazards: [Tornado]
  distance_km: 50
  days_before: 0
  days_after: 1
  export_geojson: docs/data/storm_events.geojson
```

Populate the catalog with historical or forecast storm-day observations. When the post-analysis window (plus optional lead/lag
days) contains at least one entry that falls within the AOI buffer, the full Sentinel processing pipeline runs and writes both
GeoJSON (`docs/data/changes.geojson`) and KML (`docs/data/changes.kml`) outputs. If no storm events are found, the command
skips the expensive remote sensing steps, clears previous change layers by writing empty GeoJSON/KML stubs, and exits with a
message indicating that no storm day was detected.

## GitHub Pages workflow
The Leaflet app is designed to live under `docs/`. Enable Pages for the repository and target the `docs/` folder to publish the latest change polygons.

## Running on AWS
For an opinionated example that packages the pipeline into a container image and runs it on AWS Batch with S3-backed inputs/outputs, see [`infra/aws-batch-example`](infra/aws-batch-example/README.md). The walkthrough covers building and pushing the image to ECR, wiring Batch compute environments and job definitions, and triggering jobs automatically on storm days with EventBridge.

## Tuning tips
- Tighten the pre/post acquisition windows around the event to avoid seasonal noise.
- Increase the SAR weight and enable `use_sentinel1_grd` for debris detection over bare ground or low vegetation.
- Adjust `min_blob_area_m2` and the elongation filter to suppress false positives from small fields or urban change.

## Future work
- Incorporate Sentinel-1 SLC coherence for finer debris detection.
- Fuse drone or aerial roughness metrics where available.
- Explore translating GSS to EF-scale intensity estimates.

## License
MIT License. See `LICENSE` for details.
