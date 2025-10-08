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
    fetch_bom_warnings,
    fetch_iem_local_storm_reports,
    save_metadata,
    write_catalog_csv,
)


@click.command()
@click.option("--aoi", "aoi_path", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--start", "start_str", required=True, help="UTC start datetime (YYYY-MM-DD or ISO 8601).")
@click.option("--end", "end_str", required=True, help="UTC end datetime (inclusive, YYYY-MM-DD or ISO 8601).")
@click.option(
    "--provider",
    type=click.Choice(["iem_lsr", "bom_warnings"], case_sensitive=False),
    default="iem_lsr",
    show_default=True,
    help="Storm report source to query.",
)
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
    "--state",
    "states",
    multiple=True,
    help="Filter BOM warnings by Australian state/territory (e.g., QLD).",
)
@click.option(
    "--phenomenon",
    "phenomena",
    multiple=True,
    help="Filter BOM warnings by phenomenon name (e.g., Severe Thunderstorm Warning).",
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
    provider: str,
    hazards: Iterable[str],
    bbox_buffer_km: float,
    states: Iterable[str],
    phenomena: Iterable[str],
    output: str,
    metadata: str | None,
) -> None:
    """Scrape storm reports into a CSV compatible with :mod:`storm_filter`."""

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

    provider_key = provider.strip().lower()
    if provider_key not in {"iem_lsr", "bom_warnings"}:
        raise click.BadParameter(
            f"Unsupported provider '{provider}'. Expected 'iem_lsr' or 'bom_warnings'."
        )

    state_list = [str(code) for code in states] if states else None
    phenomenon_list = [str(name) for name in phenomena] if phenomena else None

    if provider_key == "iem_lsr":
        reports = fetch_iem_local_storm_reports(
            start=start_dt,
            end=end_inclusive,
            bounds=bounds,
            hazards=hazard_list,
        )
    else:
        reports = fetch_bom_warnings(
            start=start_dt,
            end=end_inclusive,
            bounds=bounds,
            hazards=hazard_list,
            states=state_list,
            phenomena=phenomenon_list,
        )

    write_catalog_csv(reports, Path(output))

    if metadata:
        metadata_path = Path(metadata)
        payload = {
            "provider": provider_key,
            "start": start_dt.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "end": end_inclusive.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "bounds": list(bounds),
            "hazards": list(hazard_list or []),
            "count": int(len(reports)),
        }
        if provider_key == "bom_warnings":
            payload.update(
                {
                    "states": list(state_list or []),
                    "phenomena": list(phenomenon_list or []),
                }
            )
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

