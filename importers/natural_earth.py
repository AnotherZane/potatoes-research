"""Read Natural Earth 1:110m Admin 0 countries as corpus rows.

A MultiPolygon country becomes one row per part, largest first: "New Zealand_1",
"New Zealand_2", ... with source_id "NZL_1", "NZL_2", ...
"""

import os

import geopandas as gpd
from shapely import wkb as shapely_wkb


METADATA_COLUMNS = [
    "ISO_A3", "ISO_A2", "ISO_N3",
    "SOVEREIGNT", "ADMIN", "NAME_LONG", "FORMAL_EN",
    "CONTINENT", "REGION_UN", "SUBREGION", "REGION_WB",
    "POP_EST", "POP_YEAR", "GDP_MD", "GDP_YEAR",
    "ECONOMY", "INCOME_GRP",
    "NE_ID", "WIKIDATAID",
    "featurecla", "scalerank", "LABELRANK", "TYPE",
]


def build_metadata(row: dict) -> dict:
    """Extract a curated metadata dict from a GeoDataFrame row."""
    meta = {}
    for col in METADATA_COLUMNS:
        val = row.get(col)
        if val is not None and not (isinstance(val, float) and val != val):  # skip NaN
            try:
                meta[col] = val.item()
            except AttributeError:
                meta[col] = val
    return meta


def read_countries(shapefile: str) -> list[dict]:
    """One row per polygon part, in shapefile order.

    Each row: name, source, source_id, geometry_wkb, crs, coordinate_space, metadata
    (a dict).  MultiPolygon countries are split into one row per part, largest first.
    """
    os.environ.setdefault("SHAPE_RESTORE_SHX", "YES")
    gdf = gpd.read_file(shapefile)
    gdf = gdf.set_crs("EPSG:4326", allow_override=True) if gdf.crs is None else gdf

    rows = []
    for _, feature in gdf.iterrows():
        geom_raw = feature.geometry
        if geom_raw is None or geom_raw.is_empty:
            continue

        base_name: str = str(feature.get("NAME") or feature.get("ADMIN") or "Unknown")

        if geom_raw.geom_type == "MultiPolygon":
            geoms = list(geom_raw.geoms)
            if len(geoms) > 1:
                geoms.sort(key=lambda p: p.area, reverse=True)  # largest first
                is_multi = True
            else:
                geoms = [geoms[0]]
                is_multi = False
        elif geom_raw.geom_type == "Polygon":
            geoms = [geom_raw]
            is_multi = False
        else:
            raise ValueError(f"{base_name!r}: unexpected geometry type {geom_raw.geom_type!r}")

        total_parts = len(geoms)

        # Base source_id (ISO_A3 preferred, then NE_ID)
        iso_a3 = feature.get("ISO_A3")
        if iso_a3 and str(iso_a3).strip() not in ("", "-99", "-1"):
            base_source_id = str(iso_a3)
        else:
            ne_id = feature.get("NE_ID")
            base_source_id = str(int(ne_id)) if ne_id else ""

        for part_idx, geom in enumerate(geoms, start=1):
            if is_multi:
                part_name = f"{base_name}_{part_idx}"
                part_num = part_idx
            else:
                part_name = base_name
                part_num = 1

            row_data = feature.to_dict()
            row_data.pop("geometry", None)
            metadata = build_metadata(row_data)

            metadata["original_geometry_type"] = geom_raw.geom_type
            if geom_raw.geom_type == "MultiPolygon":
                metadata["original_polygon_count"] = len(list(geom_raw.geoms))

            if is_multi:
                metadata["part_number"] = part_num
                metadata["total_parts"] = total_parts
                metadata["base_name"] = base_name

            # source_id is unique per part
            if is_multi:
                source_id = f"{base_source_id}_{part_num}" if base_source_id else f"{base_name.replace(' ', '_')}_{part_num}"
            else:
                source_id = base_source_id
            if not source_id:
                raise ValueError(f"{part_name!r}: no ISO_A3 or NE_ID to identify it by")

            rows.append({
                "name": part_name,
                "source": "natural_earth",
                "source_id": source_id,
                "geometry_wkb": shapely_wkb.dumps(geom, include_srid=False),
                "crs": "EPSG:4326",
                "coordinate_space": "geographic",
                "metadata": metadata,
            })
    return rows
