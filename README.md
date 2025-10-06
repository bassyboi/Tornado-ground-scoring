# Tornado Ground-Scour Mapper

Detect tornado-driven ground change from Sentinel imagery, score the severity, and publish an interactive web map — all from a single Python pipeline.

## Features
- Fetch Sentinel-2 (and optional Sentinel-1 RTC) imagery through public STAC APIs.
- Build median pre/post mosaics and compute a composite change score using vegetation loss, brightness increase, and SAR backscatter change.
- Threshold the change map, filter polygons by elongation, and assign a Ground-Scour Score (GSS 0–5).
- Export GeoTIFF rasters, GeoJSON polygons, and a ready-to-host Leaflet map under `docs/` for GitHub Pages.
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

## GitHub Pages workflow
The Leaflet app is designed to live under `docs/`. Enable Pages for the repository and target the `docs/` folder to publish the latest change polygons.

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
