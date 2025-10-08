# Tornado Ground-Scour Mapper

Detect tornado-driven ground change from Sentinel imagery, score the severity, and publish an interactive web map — all from a single Python pipeline.

## Features
- Fetch Sentinel-2 (and optional Sentinel-1 RTC) imagery through public STAC APIs.
- Build median pre/post mosaics and compute a composite change score using vegetation loss, brightness increase, and SAR backscatter change.
- Threshold the change map, filter polygons by elongation, and assign a Ground-Scour Score (GSS 0–5).
- Export GeoTIFF rasters, GeoJSON polygons, and a ready-to-host Leaflet map under `docs/` for GitHub Pages.
- Emit a Google Earth-ready KML file alongside the GeoJSON export for quick sharing.
- Continuous integration workflow that runs the pipeline and publishes a GitHub Pages site on every push.

## Quick start
1. **Install dependencies** (Python 3.11):
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```
2. **Configure the analysis** by editing `data/config.yaml`. Update the AOI path, date windows, target projection, weights, threshold, and web map metadata as needed.
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
- **Search fallback**: The optional `search` block lets you expand the pre/post date windows when no Sentinel imagery is found. Set `auto_expand_max_days` to the maximum ± days of padding to try, and `auto_expand_step_days` to the increment between attempts.
- **Web map**: Customize the title, description, and initial map viewpoint.
- **Stack projection**: Set `stack_epsg` to the EPSG code you want mosaics resampled into. The sample configuration uses EPSG:3577 (Australian Albers) so Sentinel imagery is reprojected to metric units across the continent.
- **Storm filter**: Optionally gate the entire run on a catalog of historical storm reports. Provide a CSV file with at least
  the event timestamp, hazard type, and coordinates. When enabled, the pipeline loads the catalog, keeps events that intersect
  the AOI (with an optional distance buffer) within the specified post-event window, exports them to GeoJSON, and aborts early
  when no qualifying storm days are found. Add the optional `auto_backfill_*` keys to slide the post-event window forward or
  backward (e.g., for backfilling newly reported storms) until qualifying events are discovered. Set
  `auto_backfill_max_days: auto` to let the filter expand as far as the catalog contains qualifying events. The storm filter can
  also auto-scrape NOAA/IEM Local Storm Reports so the catalog stays fresh without manual downloads.

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
  auto_backfill_max_days: auto
  auto_backfill_step_days: 1
  auto_backfill_directions: [backward, forward]
  export_geojson: docs/data/storm_events.geojson
  scrape:
    provider: iem_lsr
    lookback_days: 7
    lookahead_days: 1
    bbox_buffer_km: 25
    hazards: [Tornado]
    export_csv: data/storm_reports_scraped.csv
    metadata_path: data/storm_reports_scrape.json
```

Populate the catalog with historical or forecast storm-day observations. When the post-analysis window (plus optional lead/lag
days) contains at least one entry that falls within the AOI buffer, the full Sentinel processing pipeline runs and writes both
GeoJSON (`docs/data/changes.geojson`) and KML (`docs/data/changes.kml`) outputs. If no storm events are found, the command
skips the expensive remote sensing steps, clears previous change layers by writing empty GeoJSON/KML stubs, and exits with a
message indicating that no storm day was detected. In `auto` mode the backfill search keeps extending the post window until it
reaches the earliest/latest catalog event that also satisfies the AOI and hazard filters, so new reports are picked up as soon
as they appear in the scraped feed.

When the optional `scrape` block is configured, the pipeline automatically queries the Iowa Environmental Mesonet (IEM) Local
Storm Reports feed (`provider: iem_lsr`) for the AOI bounds before each run. The `lookback_days` and `lookahead_days` keys
extend the scrape window around the configured `post_from`/`post_to` dates, while `bbox_buffer_km` grows the geographic search
radius beyond the AOI outline. Successful scrapes write the updated catalog to `export_csv` and capture provenance details in
`metadata_path`. If the scrape fails (e.g., due to a network outage) the pipeline falls back to the static `catalog` so
previously downloaded events can still trigger the run.

To refresh the catalog outside of the pipeline, run the standalone scraper CLI:

```bash
python -m src.fetch_storm_reports \
  --aoi data/aoi_example.geojson \
  --start 2023-04-15 \
  --end 2023-04-24 \
  --hazard Tornado \
  --bbox-buffer-km 25 \
  --output data/storm_reports_scraped.csv \
  --metadata data/storm_reports_scrape.json
```

The command writes a CSV compatible with the storm filter and an optional JSON metadata log documenting the scrape window,
bounds, hazards, and record count.

## GitHub Pages workflow

The Leaflet app is designed to live under `docs/`. A GitHub Actions workflow (`.github/workflows/build.yaml`) runs the pipeline on every push to `main`, packages the `docs/` folder, and deploys it to GitHub Pages. A daily scheduled run (06:00 UTC) also executes the pipeline so the storm filter can automatically publish new tracks whenever qualifying Australian storm reports appear in the catalog. To get automated map updates online:

1. Open the repository settings → **Pages** and choose the "GitHub Actions" build option.
2. Push your configuration changes to `main`. The workflow will fetch imagery, run the change-detection pipeline, and upload the rendered `docs/` directory.
3. Once the `Deploy to GitHub Pages` job completes, GitHub provides a public URL (surfaced in the workflow summary) that serves the latest map along with the generated storm-day report.

If you prefer to publish the map manually, run the pipeline locally, commit the updated `docs/` folder, and serve it with GitHub Pages or any static hosting provider.

### Scoring metadata export

Every pipeline execution now writes a machine-readable scoring summary to `docs/data/scoring.json`. The file captures the component weights, the resolved change-score threshold, the five Ground-Scour Score breakpoints, and polygon statistics (counts, area, mean-score range). The GitHub Actions workflow publishes this file alongside the Leaflet app, so the repository always documents the scoring system used for the latest run. If you customise the configuration, rerun the pipeline to regenerate the scoring metadata before committing.

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
