# Tornado Ground-Scour Mapper

Detect tornado-driven ground change from Sentinel-2 imagery and publish an interactive Leaflet map — all from a single Python pipeline.

## Features
- Fetch median pre- and post-event Sentinel-2 Level-2A mosaics from public STAC APIs (Microsoft Planetary Computer with an Element84 fallback).
- Compute vegetation loss (ΔNDVI) and brightness increase composites, then blend them into a single change score.
- Derive an automatic Otsu threshold (or apply a manual cutoff), rasterise the change mask, and vectorise polygons.
- Export GeoTIFF rasters to `artifacts/` plus GeoJSON, KML, and metadata to `docs/data/`.
- Ship a static Leaflet map (`docs/index.html`) that loads the generated GeoJSON, ready for `python -m http.server` or GitHub Pages hosting.

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
   The run will create GeoTIFFs in `artifacts/` and populate `docs/data/` with GeoJSON, KML, and metadata consumed by the Leaflet app at `docs/index.html`.
4. **Serve the map locally** (optional):
   ```bash
   python -m http.server --directory docs 8000
   ```
   Navigate to `http://localhost:8000` to view the map.

## Configuration reference
See `data/config.yaml` for the default options:
- **AOI**: Path to a GeoJSON/GeoPackage polygon describing the area of interest.
- **STAC search**: Provide one or more public catalog URLs. The pipeline tries each in order until scenes are found.
- **Acquisition windows**: `pre_from`/`pre_to` and `post_from`/`post_to` define the Sentinel-2 date ranges.
- **Weights**: Tune the component weights to emphasise vegetation loss (ΔNDVI) or brightness increase.
- **Threshold**: Use `"otsu"` for automatic selection or set a numeric cutoff between 0 and 1.
- **Stack projection**: Set `stack_epsg` to the projected CRS used for mosaics and polygon areas. The demo configuration uses EPSG:3577 (Australian Albers).
- **Web map**: Customise the title, description, and initial map viewpoint in the generated Leaflet page.

## Outputs
Each successful run produces:

- `artifacts/pre_mosaic.tif` and `artifacts/post_mosaic.tif`: four-band (blue, green, red, NIR) composites in the target projection.
- `artifacts/change_score.tif`: blended ΔNDVI/brightness change score (0–1).
- `artifacts/change_mask.tif`: binary mask after thresholding.
- `docs/data/changes.geojson`: vectorised polygons in WGS84 with area, perimeter, and mean score attributes.
- `docs/data/changes.kml`: Google Earth-friendly version of the polygons.
- `docs/data/meta.json`: run metadata (date windows, weights, catalog provenance, polygon count).

Serve the map locally with:

```bash
python -m http.server --directory docs 8000
```

Then open `http://localhost:8000` in a browser to inspect the generated polygons on top of OpenStreetMap tiles.

## Roadmap
- Add Sentinel-1 SAR coherence and elongation/bearing filters.
- Layer in Ground-Scour Score (0–5) classification buckets.
- Automate Bureau of Meteorology / IEM storm-report scraping and GitHub Pages deployment.

## License
MIT License. See `LICENSE` for details.
