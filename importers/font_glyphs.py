"""Read font glyph outlines (TTF/OTF) as corpus rows.

Each character's outline is flattened (flatten_steps segments per Bezier curve),
scaled to an em of em_size units and simplified (Douglas-Peucker, simplify_tol em
units; 0 = off).  Only outer boundaries are kept, so holes are dropped.  A glyph with
several components ("i", "j") becomes one row per component, largest first: "i_1",
"i_2", ...  Coordinates are font units, baseline at y = 0, +y up.
"""

import os
import re

from fontTools.ttLib import TTFont
from fontTools.pens.basePen import BasePen

from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
from shapely import wkb as shapely_wkb


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def slugify(text: str) -> str:
    """Filesystem/id-safe slug: alnum runs joined by single underscores."""
    slug = re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_")
    return slug or "font"


# ---------------------------------------------------------------------------
# Glyph outline extraction
# ---------------------------------------------------------------------------

def _quad_point(p0, c, p1, t):
    """Point on a quadratic Bezier at parameter t (TrueType curves)."""
    mt = 1.0 - t
    return (mt * mt * p0[0] + 2 * mt * t * c[0] + t * t * p1[0],
            mt * mt * p0[1] + 2 * mt * t * c[1] + t * t * p1[1])


def _cubic_point(p0, c1, c2, p1, t):
    """Point on a cubic Bezier at parameter t (CFF/OpenType curves)."""
    mt = 1.0 - t
    a, b, c, d = mt ** 3, 3 * mt * mt * t, 3 * mt * t * t, t ** 3
    return (a * p0[0] + b * c1[0] + c * c2[0] + d * p1[0],
            a * p0[1] + b * c1[1] + c * c2[1] + d * p1[1])


class FlattenPen(BasePen):
    """
    A fontTools pen that flattens a glyph's outline into a list of contours,
    each a list of (x, y) points. Handles both quadratic (TrueType) and cubic
    (CFF) segments; BasePen takes care of TrueType implied on-curve points.
    """

    def __init__(self, glyph_set, steps: int = 12):
        super().__init__(glyph_set)
        self.contours = []
        self._current = None
        self._steps = steps

    def _moveTo(self, pt):
        self._current = [pt]

    def _lineTo(self, pt):
        self._current.append(pt)

    def _curveToOne(self, bcp1, bcp2, pt):
        p0 = self._current[-1]
        for i in range(1, self._steps + 1):
            self._current.append(_cubic_point(p0, bcp1, bcp2, pt, i / self._steps))

    def _qCurveToOne(self, bcp, pt):
        p0 = self._current[-1]
        for i in range(1, self._steps + 1):
            self._current.append(_quad_point(p0, bcp, pt, i / self._steps))

    def _closePath(self):
        if self._current:
            self.contours.append(self._current)
            self._current = None

    def _endPath(self):
        self._closePath()


def _dedupe(points, eps: float = 1e-9):
    """Drop consecutive duplicate points and any trailing closing duplicate."""
    out = []
    for p in points:
        if not out or abs(p[0] - out[-1][0]) > eps or abs(p[1] - out[-1][1]) > eps:
            out.append((float(p[0]), float(p[1])))
    if len(out) > 1 and abs(out[0][0] - out[-1][0]) <= eps and abs(out[0][1] - out[-1][1]) <= eps:
        out.pop()
    return out


def extract_outer_polygons(glyph_set, glyph_name, steps, scale, simplify_tol):
    """
    Return (parts, stats) for one glyph.

    parts : list of shapely Polygons (exterior-only, oriented CCW), one per
            disconnected outer component, sorted by descending area.
    stats : dict with original_contour_count and holes_dropped.
    """
    pen = FlattenPen(glyph_set, steps)
    glyph_set[glyph_name].draw(pen)

    contours = []
    for raw in pen.contours:
        pts = _dedupe([(x * scale, y * scale) for (x, y) in raw])
        if len(pts) >= 3:
            contours.append(pts)

    if not contours:
        return [], {"original_contour_count": 0, "holes_dropped": 0}

    # Build helper polygons purely to decide nesting (outer vs hole/counter).
    helpers = []
    for pts in contours:
        poly = Polygon(pts)
        if not poly.is_valid:
            poly = poly.buffer(0)
        helpers.append(poly)

    def is_inner(i):
        pi = helpers[i]
        if pi.is_empty:
            return True
        rp = pi.representative_point()
        for j, pj in enumerate(helpers):
            if j == i or pj.is_empty:
                continue
            if pj.area > pi.area and pj.contains(rp):
                return True
        return False

    outer_idx = [i for i in range(len(contours)) if not is_inner(i)]

    parts = []
    for i in outer_idx:
        # Store the raw outer ring (exterior only, no holes) to keep it simple.
        poly = orient(Polygon(contours[i]), sign=1.0)  # CCW
        if simplify_tol and simplify_tol > 0:
            poly = orient(poly.simplify(simplify_tol, preserve_topology=True), sign=1.0)
        if poly.is_valid and not poly.is_empty and poly.exterior is not None:
            parts.append(poly)

    parts.sort(key=lambda p: p.area, reverse=True)  # largest first

    stats = {
        "original_contour_count": len(contours),
        "holes_dropped": len(contours) - len(outer_idx),
    }
    return parts, stats


# ---------------------------------------------------------------------------
# Font metadata
# ---------------------------------------------------------------------------

def font_identity(font: TTFont, font_path: str) -> dict:
    """Pull family/subfamily/postscript names and a slug for the font."""
    def name(nid):
        try:
            val = font["name"].getDebugName(nid)
        except Exception:
            val = None
        return val

    family = name(16) or name(1) or os.path.splitext(os.path.basename(font_path))[0]
    subfamily = name(17) or name(2) or "Regular"
    ps_name = name(6) or f"{family}-{subfamily}"
    slug = slugify(ps_name)
    return {
        "font_family": family,
        "font_subfamily": subfamily,
        "font_postscript_name": ps_name,
        "font_slug": slug,
        "units_per_em": int(font["head"].unitsPerEm),
    }


def char_class(ch: str) -> str:
    if "A" <= ch <= "Z":
        return "uppercase"
    if "a" <= ch <= "z":
        return "lowercase"
    if "0" <= ch <= "9":
        return "digit"
    return "other"


def default_charset() -> str:
    upper = "".join(chr(c) for c in range(ord("A"), ord("Z") + 1))
    lower = "".join(chr(c) for c in range(ord("a"), ord("z") + 1))
    digit = "".join(chr(c) for c in range(ord("0"), ord("9") + 1))
    return upper + lower + digit


# ---------------------------------------------------------------------------
# Core import logic
# ---------------------------------------------------------------------------

def read_font(
        font_path: str,
        chars: str = None,
        em_size: int = 1000,
        flatten_steps: int = 12,
        simplify_tol: float = 2.0,
) -> list[dict]:
    """One row per character (per disconnected component), in charset order.

    A character the font lacks, or whose outline does not survive simplification,
    yields no row.

    Each row: name, source, source_id, geometry_wkb, crs, coordinate_space, metadata
    (a dict).
    """
    chars = chars if chars else default_charset()

    font = TTFont(font_path)
    ident = font_identity(font, font_path)
    glyph_set = font.getGlyphSet()
    cmap = font.getBestCmap()

    upm = ident["units_per_em"]
    scale = (em_size / upm) if em_size and em_size > 0 else 1.0
    rows = []
    for ch in chars:
        cp = ord(ch)
        if cp not in cmap:
            continue            # not in the font

        glyph_name = cmap[cp]
        parts, stats = extract_outer_polygons(
            glyph_set, glyph_name, flatten_steps, scale, simplify_tol
        )

        if not parts:
            continue            # nothing survives simplification

        is_multi = len(parts) > 1
        total_parts = len(parts)
        base_source_id = f"{ident['font_slug']}_{cp:04X}"

        for part_idx, poly in enumerate(parts, start=1):
            if is_multi:
                part_name = f"{ch}_{part_idx}"
                source_id = f"{base_source_id}_{part_idx}"
            else:
                part_name = ch
                source_id = base_source_id

            metadata = {
                "character": ch,
                "unicode": f"U+{cp:04X}",
                "codepoint": cp,
                "glyph_name": glyph_name,
                "char_class": char_class(ch),
                "font_family": ident["font_family"],
                "font_subfamily": ident["font_subfamily"],
                "font_postscript_name": ident["font_postscript_name"],
                "font_file": os.path.basename(font_path),
                "units_per_em": upm,
                "em_size": em_size,
                "scale_factor": scale,
                "flatten_steps": flatten_steps,
                "simplify_tolerance": simplify_tol,
                "original_contour_count": stats["original_contour_count"],
                "holes_dropped": stats["holes_dropped"],
                "vertex_count": len(poly.exterior.coords) - 1,
                "bbox": list(poly.bounds),  # (minx, miny, maxx, maxy)
                "is_simple": bool(poly.is_valid and poly.exterior.is_simple),
            }
            if is_multi:
                metadata["part_number"] = part_idx
                metadata["total_parts"] = total_parts
                metadata["base_name"] = ch

            rows.append({
                "name": part_name,
                "source": "font",
                "source_id": source_id,
                "geometry_wkb": shapely_wkb.dumps(poly, include_srid=False),
                "crs": None,                       # font units: no CRS
                "coordinate_space": "cartesian",
                "metadata": metadata,
            })
    return rows
