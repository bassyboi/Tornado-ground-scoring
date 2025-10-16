"""Ground Scour Score filtering and scoring."""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any, Mapping

import geopandas as gpd
import numpy as np

LOGGER = logging.getLogger(__name__)


@dataclass
class ScoringResult:
    """Container for scored polygons and summary statistics."""

    polygons: gpd.GeoDataFrame
    summary: dict[str, Any]


def apply(
    features: gpd.GeoDataFrame,
    *,
    filters_cfg: Mapping[str, Any],
    rules_cfg: Mapping[str, Any],
) -> ScoringResult:
    """Filter polygons and compute Ground Scour Scores."""

    total_input = int(len(features))
    filters = {
        "min_area_m2": _as_float(filters_cfg, "min_area_m2", 0.0),
        "min_elongation": _as_float(filters_cfg, "min_elongation", 0.0),
        "max_compactness": _as_float(filters_cfg, "max_compactness", float("inf")),
        "min_score_mean": _as_float(filters_cfg, "min_score_mean", -float("inf")),
    }

    rules = {
        "e2_elongation": _as_float(rules_cfg, "e2_elongation", 0.0),
        "e2_min_length_m": _as_float(rules_cfg, "e2_min_length_m", 0.0),
        "ndvi_loss": _as_float(rules_cfg, "ndvi_loss", 0.0),
        "bright_increase": _as_float(rules_cfg, "bright_increase", 0.0),
        "max_compactness_for_bonus": _as_float(rules_cfg, "max_compactness_for_bonus", float("inf")),
        "align_tol_deg": _as_float(rules_cfg, "align_tol_deg", 180.0),
    }

    if features.empty:
        summary = _build_summary(gpd.GeoDataFrame(geometry=[], crs=features.crs), total_input, None)
        summary["filters"] = {k: float(v) for k, v in filters.items()}
        summary["score_rules"] = {k: float(v) for k, v in rules.items()}
        return ScoringResult(features.copy(), summary)

    keep_indices: list[int] = []
    for idx, row in features.iterrows():
        reasons: list[str] = []
        area = float(row.get("area_m2", float("nan")))
        if not np.isfinite(area) or area < filters["min_area_m2"]:
            reasons.append(f"area {area:.1f} < {filters['min_area_m2']:.1f}")
        elongation = float(row.get("elongation", float("nan")))
        if not np.isfinite(elongation) or elongation < filters["min_elongation"]:
            reasons.append(f"elongation {elongation:.2f} < {filters['min_elongation']:.2f}")
        compactness = float(row.get("compactness", float("nan")))
        if np.isfinite(compactness) and compactness > filters["max_compactness"]:
            reasons.append(f"compactness {compactness:.3f} > {filters['max_compactness']:.3f}")
        elif not np.isfinite(compactness):
            reasons.append("compactness missing")
        score_mean = float(row.get("change_score_mean", float("nan")))
        if not np.isfinite(score_mean) or score_mean < filters["min_score_mean"]:
            reasons.append(f"mean score {score_mean:.3f} < {filters['min_score_mean']:.3f}")

        if reasons:
            LOGGER.info("Dropping polygon %s: %s", idx, "; ".join(reasons))
        else:
            keep_indices.append(idx)

    scored = features.loc[keep_indices].copy() if keep_indices else gpd.GeoDataFrame(geometry=[], crs=features.crs)
    scored = scored.reset_index(drop=True)
    LOGGER.info("Kept %s of %s polygons after filters", len(scored), total_input)

    if scored.empty:
        summary = _build_summary(scored, total_input, None)
        summary["filters"] = {k: float(v) for k, v in filters.items()}
        summary["score_rules"] = {k: float(v) for k, v in rules.items()}
        return ScoringResult(scored, summary)

    cluster = _bearing_cluster(scored["bearing_deg"].to_numpy(dtype=float))
    if cluster is not None:
        LOGGER.info("Bearing cluster: %.2f°", cluster)
    else:
        LOGGER.info("Bearing cluster undefined (insufficient polygons)")

    gss_values: list[int] = []
    explanations: list[str] = []
    alignments: list[float] = []

    for _, row in scored.iterrows():
        score = 0
        parts: list[str] = []

        elongation = float(row.get("elongation", float("nan")))
        length = float(row.get("length_m", float("nan")))
        if np.isfinite(elongation) and np.isfinite(length):
            if elongation >= rules["e2_elongation"] and length >= rules["e2_min_length_m"]:
                score += 2
                parts.append(f"+2 elongation ({elongation:.2f}×, {length:.0f} m)")
        ndvi_delta = float(row.get("delta_ndvi_mean", float("nan")))
        if np.isfinite(ndvi_delta) and ndvi_delta <= rules["ndvi_loss"]:
            score += 1
            parts.append(f"+1 ΔNDVI {ndvi_delta:.3f}")
        bright_delta = float(row.get("delta_brightness_mean", float("nan")))
        if np.isfinite(bright_delta) and bright_delta >= rules["bright_increase"]:
            score += 1
            parts.append(f"+1 Δbrightness {bright_delta:.3f}")
        compactness = float(row.get("compactness", float("nan")))
        if np.isfinite(compactness) and compactness <= rules["max_compactness_for_bonus"]:
            score += 1
            parts.append(f"+1 compactness {compactness:.3f}")

        if cluster is not None and np.isfinite(row.get("bearing_deg", float("nan"))):
            align = _angular_distance(float(row["bearing_deg"]), cluster)
            alignments.append(align)
            if align <= rules["align_tol_deg"]:
                score += 1
                parts.append(f"+1 alignment Δ{align:.1f}°")
        else:
            alignments.append(float("nan"))

        gss = int(max(0, min(score, 5)))
        gss_values.append(gss)
        explanations.append("; ".join(parts) if parts else "Baseline score")

    scored["gss"] = gss_values
    scored["gss_explain"] = explanations
    scored["bearing_align_deg"] = alignments

    summary = _build_summary(scored, total_input, cluster)
    summary["filters"] = {k: float(v) for k, v in filters.items()}
    summary["score_rules"] = {k: float(v) for k, v in rules.items()}
    return ScoringResult(scored, summary)


def _as_float(cfg: Mapping[str, Any], key: str, default: float) -> float:
    value = cfg.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _bearing_cluster(bearings: np.ndarray) -> float | None:
    finite = bearings[np.isfinite(bearings)]
    if finite.size == 0:
        return None
    angles = np.deg2rad(finite * 2.0)
    sin_mean = np.sin(angles).mean()
    cos_mean = np.cos(angles).mean()
    if sin_mean == 0 and cos_mean == 0:
        return None
    mean = 0.5 * math.degrees(math.atan2(sin_mean, cos_mean))
    return (mean + 180.0) % 180.0


def _angular_distance(angle: float, reference: float) -> float:
    diff = abs((angle % 180.0) - (reference % 180.0))
    return min(diff, 180.0 - diff)


def _build_summary(polygons: gpd.GeoDataFrame, total_input: int, cluster: float | None) -> dict[str, Any]:
    counts = {str(level): 0 for level in range(6)}
    if not polygons.empty and "gss" in polygons:
        for level, count in polygons["gss"].value_counts().items():
            counts[str(int(level))] = int(count)
    mean_elongation = float(polygons["elongation"].mean()) if "elongation" in polygons and len(polygons) else None
    summary: dict[str, Any] = {
        "total": int(len(polygons)),
        "total_input": total_input,
        "per_gss": counts,
        "bearing_cluster_deg": float(cluster) if cluster is not None else None,
        "mean_elongation": mean_elongation if mean_elongation is None or np.isfinite(mean_elongation) else None,
    }
    if "bearing_align_deg" in polygons and len(polygons):
        mean_align = polygons["bearing_align_deg"].dropna()
        summary["mean_alignment_deg"] = float(mean_align.mean()) if not mean_align.empty else None
    return summary
