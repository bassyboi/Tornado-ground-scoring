"""Utilities for fetching Sentinel-2 mosaics from STAC APIs."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import date
from typing import Iterable, List, Sequence

import numpy as np
import stackstac
import xarray as xr
from pystac import Item
from pystac_client import Client
from shapely.geometry import mapping

try:  # pragma: no cover - optional dependency
    import planetary_computer
except Exception:  # noqa: S110 - broad except to handle import/runtime errors
    planetary_computer = None

LOGGER = logging.getLogger(__name__)

S2_COLLECTION = "sentinel-2-l2a"
S2_ASSETS = {"B02": "blue", "B03": "green", "B04": "red", "B08": "nir"}
DEFAULT_RESOLUTION = 10


@dataclass
class MosaicResult:
    """Container for a projected Sentinel-2 mosaic."""

    data: xr.Dataset
    catalog_url: str
    item_count: int


def fetch_s2_mosaic(
    catalogs: Sequence[str],
    start: date,
    end: date,
    geometry_wgs84,
    bounds_projected,
    target_epsg: int,
) -> MosaicResult:
    """Fetch a median Sentinel-2 mosaic for the supplied window.

    Parameters
    ----------
    catalogs:
        Ordered list of STAC catalog URLs to attempt.
    start, end:
        Inclusive start/end dates for the acquisition window.
    geometry_wgs84:
        AOI geometry in WGS84 coordinates.
    bounds_projected:
        Bounding box (minx, miny, maxx, maxy) in the projected CRS; used to clip
        the mosaic before writing to disk.
    target_epsg:
        EPSG code of the target projected CRS.
    """

    timerange = f"{start.isoformat()}/{end.isoformat()}"
    last_error: Exception | None = None
    for catalog_url in catalogs:
        if "planetarycomputer" in catalog_url and planetary_computer is None:
            LOGGER.warning(
                "planetary-computer client unavailable; skipping catalog %s",
                catalog_url,
            )
            continue
        try:
            client = Client.open(catalog_url)
        except Exception as exc:  # pragma: no cover - network issues
            LOGGER.warning("Failed to open STAC catalog %s: %s", catalog_url, exc)
            last_error = exc
            continue

        search = client.search(
            collections=[S2_COLLECTION],
            datetime=timerange,
            intersects=mapping(geometry_wgs84),
        )
        items = list(search.get_items())
        if not items:
            LOGGER.info("No Sentinel-2 items found in %s for %s", timerange, catalog_url)
            continue

        LOGGER.info(
            "Found %s Sentinel-2 items between %s and %s using %s",
            len(items),
            start,
            end,
            catalog_url,
        )
        if "planetarycomputer" in catalog_url:
            if planetary_computer is None:
                LOGGER.warning(
                    "planetary-computer client unavailable; skipping catalog %s",
                    catalog_url,
                )
                continue
            items = [planetary_computer.sign(item) for item in items]

        mosaic = _stack_items(items, bounds_projected, target_epsg)
        return MosaicResult(data=mosaic, catalog_url=catalog_url, item_count=len(items))

    if last_error:
        raise RuntimeError(
            "Failed to fetch Sentinel-2 data from all catalogs; last error was"
            f" {last_error}"
        )
    raise RuntimeError(
        "No Sentinel-2 scenes were found for the requested time range in any catalog"
    )


def _stack_items(
    items: Iterable[Item], bounds_projected, target_epsg: int
) -> xr.Dataset:
    """Stack and median-composite the provided STAC items."""

    stack = stackstac.stack(
        items,
        assets=list(S2_ASSETS.keys()),
        resolution=DEFAULT_RESOLUTION,
        epsg=target_epsg,
        bounds=bounds_projected,
        bounds_crs=f"EPSG:{target_epsg}",
        chunks={},
        fill_value=np.nan,
    )
    stack = stack.where(np.isfinite(stack))
    composite = stack.median(dim="time", skipna=True)
    dataset = composite.to_dataset(dim="band").rename(S2_ASSETS)
    dataset = dataset.astype("float32") / 10000.0
    dataset = dataset.compute()
    dataset.rio.write_crs(f"EPSG:{target_epsg}", inplace=True)
    return dataset
