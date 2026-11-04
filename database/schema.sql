-- The polygon corpus, built by reproduce/build_corpus.py.  Results are keyed by id,
-- assigned in (source, source_id) order.

CREATE TABLE IF NOT EXISTS polygons (
    id               INTEGER PRIMARY KEY,
    name             TEXT NOT NULL,     -- not unique ("1" exists in ten fonts)
    source           TEXT NOT NULL,     -- 'natural_earth' or 'font'
    source_id        TEXT NOT NULL,     -- unique within source: "NZL_2", "Roboto_Bold_0041_1"
    geometry_wkb     BLOB NOT NULL,     -- Polygon, WKB
    crs              TEXT,              -- 'EPSG:4326' for natural_earth; NULL for font (font units)
    coordinate_space TEXT NOT NULL,     -- 'geographic' or 'cartesian'
    metadata_json    TEXT NOT NULL,     -- source fields and import parameters

    -- Descriptors of the exterior ring.  Exact values are rationals as "num/den" text.
    n_vertices       INTEGER NOT NULL,  -- n
    n_reflex         INTEGER NOT NULL,  -- reflex corners
    area_exact       TEXT NOT NULL,     -- |shoelace| / 2, exact
    area             REAL NOT NULL,     -- area_exact as a float
    hull_area_exact  TEXT NOT NULL,     -- area of the convex hull of the vertices, exact
    hull_area        REAL NOT NULL,     -- hull_area_exact as a float
    orientation      TEXT NOT NULL,     -- 'ccw' or 'cw' (sign of the exact shoelace sum)
    min_x            REAL NOT NULL,     -- bounding box, in the polygon's own coordinates
    min_y            REAL NOT NULL,
    max_x            REAL NOT NULL,
    max_y            REAL NOT NULL,
    poly_sha         TEXT NOT NULL,     -- SHA-1 of the vertex list, as in result rows

    CHECK (source IN ('natural_earth', 'font')),
    CHECK (orientation IN ('ccw', 'cw')),
    CHECK (n_reflex BETWEEN 1 AND n_vertices - 3),
    UNIQUE (source, source_id)
);
