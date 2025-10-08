"""CLI to scrape storm reports into a local CSV catalog."""
from __future__ import annotations

import logging
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable

import click

from .aoi_utils import load_aoi
from .storm_reports import (
    buffered_bounds,
    fetch_iem_local_storm_reports,
    save_metadata,
    write_catalog_csv,
)


@click.command()
@click.option("--aoi", "aoi_path", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--start", "start_str", required=True, help="UTC start datetime (YYYY-MM-DD or ISO 8601).")
@click.option("--end", "end_str", required=True, help="UTC end datetime (inclusive, YYYY-MM-DD or ISO 8601).")
@click.option(
    "--hazard",
    "hazards",
    multiple=True,
    help="Optional hazard names to retain (e.g., Tornado). Repeat for multiple filters.",
)
@click.option(
    "--bbox-buffer-km",
    default=0.0,
    type=float,
    show_default=True,
    help="Expand the AOI bounds before scraping (kilometres).",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False),
    required=True,
    help="Destination CSV path for the scraped catalog.",
)
@click.option(
    "--metadata",
    type=click.Path(dir_okay=False),
    default=None,
    help="Optional JSON metadata path describing the scrape parameters.",
)
def main(
    aoi_path: str,
    start_str: str,
    end_str: str,
    hazards: Iterable[str],
    bbox_buffer_km: float,
    output: str,
    metadata: str | None,
) -> None:
    """Scrape IEM Local Storm Reports into a CSV compatible with storm_filter."""

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    start_dt = _parse_datetime(start_str)
    end_dt = _parse_datetime(end_str)
    if end_dt < start_dt:
        raise click.BadParameter("--end must be on or after --start")

    # Make the end timestamp inclusive by extending to the start of the next day
    end_inclusive = end_dt + timedelta(days=1)

    aoi = load_aoi(aoi_path)
    bounds = buffered_bounds(aoi, buffer_km=bbox_buffer_km)
    hazard_list = [str(h) for h in hazards] if hazards else None

    reports = fetch_iem_local_storm_reports(
        start=start_dt,
        end=end_inclusive,
        bounds=bounds,
        hazards=hazard_list,
    )
    write_catalog_csv(reports, Path(output))

    if metadata:
        metadata_path = Path(metadata)
        payload = {
            "provider": "iem_lsr",
            "start": start_dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "end": end_inclusive.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "bounds": list(bounds),
            "hazards": list(hazard_list or []),
            "count": int(len(reports)),
        }
        save_metadata(payload, metadata_path)


def _parse_datetime(value: str) -> datetime:
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise click.BadParameter(f"Invalid datetime: {value}") from exc
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    else:
        parsed = parsed.astimezone(timezone.utc)
    return parsed


if __name__ == "__main__":
    main()

