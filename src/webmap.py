"""Helpers for preparing static web map assets."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict

from .config import PipelineConfig

LEAFLET_TEMPLATE = """<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>{title}</title>
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <link rel=\"stylesheet\" href=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.css\" />
  <style>
    body {{ margin: 0; font-family: system-ui, -apple-system, BlinkMacSystemFont, sans-serif; color: #1f2933; }}
    header {{ padding: 1.5rem; background: #0f172a; color: white; }}
    header h1 {{ margin: 0 0 0.5rem 0; font-size: 1.75rem; }}
    header p {{ margin: 0; max-width: 720px; line-height: 1.4; }}
    #map {{ height: calc(100vh - 140px); width: 100%; }}
    .legend {{ position: absolute; bottom: 20px; right: 20px; background: rgba(255, 255, 255, 0.92); padding: 0.75rem 1rem; border-radius: 8px; box-shadow: 0 4px 10px rgba(15, 23, 42, 0.2); }}
    .legend h3 {{ margin-top: 0; font-size: 1rem; margin-bottom: 0.5rem; }}
    .legend p {{ margin: 0; font-size: 0.85rem; }}
  </style>
</head>
<body>
  <header>
    <h1>{title}</h1>
    <p>{description}</p>
  </header>
  <div id=\"map\"></div>
  <div class=\"legend\">
    <h3>Layers</h3>
    <p><span style=\"display:inline-block;width:14px;height:14px;background:#ff7800;opacity:0.45;border:2px solid #b45309;margin-right:6px;vertical-align:middle;\"></span> Change polygons</p>
  </div>
  <script src=\"https://unpkg.com/leaflet@1.9.4/dist/leaflet.js\"></script>
  <script>
    const metaUrl = 'data/meta.json';
    const dataUrl = 'data/changes.geojson';

    const map = L.map('map');
    const osm = L.tileLayer('https://tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
      attribution: '&copy; OpenStreetMap contributors',
      maxZoom: 19
    }});
    osm.addTo(map);

    fetch(metaUrl)
      .then(response => response.json())
      .then(meta => {{
        const center = meta.webmap?.center || [0, 0];
        const zoom = meta.webmap?.zoom || 6;
        map.setView(center, zoom);
        return fetch(dataUrl).then(resp => resp.json()).then(data => {{
          const layer = L.geoJSON(data, {{
            style: () => ({
              color: '#b45309',
              weight: 2,
              fillColor: '#ff7800',
              fillOpacity: 0.45
            }}),
            onEachFeature: (feature, layer) => {{
              const props = feature.properties || {{}};
              const area = props.area_m2 ? `${{(props.area_m2 / 10000).toFixed(2)}} ha` : 'n/a';
              const score = props.mean_score != null ? props.mean_score.toFixed(2) : 'n/a';
              layer.bindPopup(`<strong>Change polygon</strong><br/>Area: ${{area}}<br/>Mean score: ${{score}}`);
            }}
          }});
          layer.addTo(map);
          if (layer.getBounds().isValid()) {{
            map.fitBounds(layer.getBounds().pad(0.1));
          }}
        }});
      }})
      .catch(error => {{
        console.error('Failed to load change data', error);
      }});
  </script>
</body>
</html>
"""


def write_metadata(
    config: PipelineConfig,
    mosaic_info: dict[str, Any],
    polygon_count: int,
    meta_path: Path,
) -> None:
    """Persist run metadata for consumption by the static map."""

    meta_path.parent.mkdir(parents=True, exist_ok=True)
    metadata: Dict[str, Any] = {
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "pre_window": [config.pre_from.isoformat(), config.pre_to.isoformat()],
        "post_window": [config.post_from.isoformat(), config.post_to.isoformat()],
        "weights": {"ndvi": config.weights.ndvi, "brightness": config.weights.brightness},
        "threshold": config.threshold,
        "stac": mosaic_info,
        "webmap": {
            "title": config.webmap.title,
            "description": config.webmap.description,
            "center": config.webmap.center,
            "zoom": config.webmap.zoom,
        },
        "polygons": {"count": int(polygon_count)},
    }
    meta_path.write_text(json.dumps(metadata, indent=2))


def ensure_leaflet_page(html_path: Path, config: PipelineConfig) -> None:
    """Create or overwrite the Leaflet index page to match the configuration."""

    html_path.parent.mkdir(parents=True, exist_ok=True)
    html = LEAFLET_TEMPLATE.replace("{title}", config.webmap.title).replace(
        "{description}", config.webmap.description
    )
    html_path.write_text(html)
