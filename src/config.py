"""Configuration loading utilities for the thin-slice pipeline."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import List, Sequence, Tuple, Union

import yaml


@dataclass(frozen=True)
class StacConfig:
    """Configuration for the upstream STAC catalogs."""

    catalogs: Tuple[str, ...]

    @classmethod
    def from_raw(cls, raw: Sequence[str] | None) -> "StacConfig":
        if not raw:
            raise ValueError("stac.catalogs must contain at least one catalog URL")
        catalogs = tuple(str(url).strip() for url in raw if str(url).strip())
        if not catalogs:
            raise ValueError("stac.catalogs must contain at least one valid URL")
        return cls(catalogs=catalogs)


@dataclass(frozen=True)
class WeightConfig:
    """Weighting applied to the change composite."""

    ndvi: float
    brightness: float

    @property
    def total(self) -> float:
        return self.ndvi + self.brightness


@dataclass(frozen=True)
class WebMapConfig:
    """Metadata for the static Leaflet web map."""

    title: str
    description: str
    center: Tuple[float, float]
    zoom: int


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level configuration consumed by the pipeline."""

    aoi_path: Path
    stac: StacConfig
    pre_from: date
    pre_to: date
    post_from: date
    post_to: date
    stack_epsg: int
    weights: WeightConfig
    threshold: Union[str, float]
    webmap: WebMapConfig

    @classmethod
    def load(cls, path: Path) -> "PipelineConfig":
        raw = yaml.safe_load(Path(path).read_text())
        if not isinstance(raw, dict):
            raise TypeError("Configuration file must contain a mapping at the top level")

        aoi_path = Path(raw.get("aoi_path"))
        if not aoi_path:
            raise KeyError("aoi_path must be provided in the configuration")
        if not aoi_path.exists():
            raise FileNotFoundError(f"AOI file not found: {aoi_path}")

        stac_cfg = StacConfig.from_raw(raw.get("stac", {}).get("catalogs"))

        pre_from = _parse_date(raw.get("pre_from"), "pre_from")
        pre_to = _parse_date(raw.get("pre_to"), "pre_to")
        post_from = _parse_date(raw.get("post_from"), "post_from")
        post_to = _parse_date(raw.get("post_to"), "post_to")

        if pre_from > pre_to:
            raise ValueError("pre_from must be earlier than or equal to pre_to")
        if post_from > post_to:
            raise ValueError("post_from must be earlier than or equal to post_to")
        if pre_to >= post_from:
            raise ValueError("pre_to must be earlier than post_from to avoid overlap")

        stack_epsg = int(raw.get("stack_epsg") or 3577)

        weights_raw = raw.get("weights") or {}
        try:
            weight_cfg = WeightConfig(
                ndvi=float(weights_raw.get("ndvi")),
                brightness=float(weights_raw.get("brightness")),
            )
        except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
            raise ValueError("weights.ndvi and weights.brightness must be numeric") from exc
        if weight_cfg.total <= 0:
            raise ValueError("The sum of weights.ndvi and weights.brightness must be > 0")

        threshold = raw.get("threshold", "otsu")
        if isinstance(threshold, str):
            if threshold.lower() != "otsu":
                raise ValueError("threshold must be 'otsu' or a numeric value in [0, 1]")
            threshold_value: Union[str, float] = "otsu"
        else:
            threshold_value = float(threshold)
            if not 0.0 <= threshold_value <= 1.0:
                raise ValueError("threshold numeric value must be between 0 and 1")

        webmap_raw = raw.get("webmap") or {}
        try:
            center_raw = webmap_raw.get("center")
            if not isinstance(center_raw, (list, tuple)) or len(center_raw) != 2:
                raise ValueError
            center = (float(center_raw[0]), float(center_raw[1]))
        except (TypeError, ValueError):  # pragma: no cover - defensive
            raise ValueError("webmap.center must be a [lat, lon] pair")
        zoom = int(webmap_raw.get("zoom", 9))
        webmap_cfg = WebMapConfig(
            title=str(webmap_raw.get("title") or "Tornado Ground-Scour"),
            description=str(
                webmap_raw.get("description")
                or "ΔNDVI + Δbrightness from Sentinel-2; Otsu-threshold polygons."
            ),
            center=center,
            zoom=zoom,
        )

        return cls(
            aoi_path=aoi_path,
            stac=stac_cfg,
            pre_from=pre_from,
            pre_to=pre_to,
            post_from=post_from,
            post_to=post_to,
            stack_epsg=stack_epsg,
            weights=weight_cfg,
            threshold=threshold_value,
            webmap=webmap_cfg,
        )


def _parse_date(value: str | None, field_name: str) -> date:
    if not value:
        raise ValueError(f"{field_name} must be provided in the configuration")
    try:
        return datetime.strptime(str(value), "%Y-%m-%d").date()
    except ValueError as exc:  # pragma: no cover - defensive
        raise ValueError(f"{field_name} must be formatted as YYYY-MM-DD") from exc
