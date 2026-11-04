"""Shared helpers for the reproduce/ scripts."""
from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                           # reproduce -> repo
DATA = os.path.join(ROOT, "data")


def environment(argv: list[str] | None = None) -> dict:
    """Machine, versions, git commit and command line, recorded beside each result."""
    def run(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=20).stdout.strip() or None
        except Exception:
            return None

    cpu = None
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip(); break
    except Exception:
        pass
    mem = None
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem = line.split()[1] + " kB"; break
    except Exception:
        pass
    virt = None
    try:
        virt = "hypervisor" in open("/proc/cpuinfo").read()
    except Exception:
        pass

    return {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "cpu": cpu,
        "cpu_count": os.cpu_count(),
        "virtualised_guest": virt,
        "memory_total": mem,
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "compiler": (run(["g++", "--version"]) or "").split("\n")[0] or None,
        "cmake": (run(["cmake", "--version"]) or "").split("\n")[0] or None,
        "git_commit": run(["git", "-C", ROOT, "rev-parse", "HEAD"]),
        "git_dirty": bool(run(["git", "-C", ROOT, "status", "--porcelain"])),
        "argv": argv if argv is not None else sys.argv,
    }


def write_sidecar(result_path: str, env: dict, summary: dict | None = None) -> str:
    side = result_path.rsplit(".", 1)[0] + ".json"
    with open(side, "w") as fh:
        json.dump({"environment": env, "summary": summary or {}}, fh, indent=1)
    return side


def polygon_stats(coords):
    """(n, n_reflex, area, sha1) of a polygon, exactly.  sha1 identifies the vertex list."""
    import hashlib
    from fractions import Fraction as Fr
    pts = [(Fr(str(x)), Fr(str(y))) for x, y in coords]
    n = len(pts)
    two_a = sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1]
                for i in range(n))
    ccw = two_a > 0
    nref = 0
    for i in range(n):
        ax, ay = pts[(i - 1) % n]; bx, by = pts[i]; cx, cy = pts[(i + 1) % n]
        z = (bx - ax) * (cy - by) - (by - ay) * (cx - bx)
        if (z < 0) if ccw else (z > 0):
            nref += 1
    h = hashlib.sha1(repr([(str(x), str(y)) for x, y in pts]).encode()).hexdigest()
    return n, nref, abs(two_a) / 2, h


def convex_hull_area(coords):
    """Exact area of the convex hull of the vertices (monotone chain)."""
    from fractions import Fraction as Fr
    pts = sorted(set((Fr(str(x)), Fr(str(y))) for x, y in coords))
    if len(pts) < 3:
        return Fr(0)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    n = len(hull)
    two_a = sum(hull[i][0] * hull[(i + 1) % n][1] - hull[(i + 1) % n][0] * hull[i][1]
                for i in range(n))
    return abs(two_a) / 2
