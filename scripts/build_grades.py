#!/usr/bin/env python3
"""Grade along every trail, measured the same way for each.

    .venv/bin/python scripts/build_grades.py            # every trail in trails.yaml
    .venv/bin/python scripts/build_grades.py tmb jmt    # only these

Reads data/osm/<trail>/ (or the trail's track_file), data/fabdem/ and the settings in trails.yaml. Steps:
  1. Track      the route relation's main-route ways, joined end to end into as few pieces as they allow, and
                pointed the way hikers usually go: pieces are chained from the trail's start point, and a loop is
                turned to its usual direction
  2. Points     one every spacing_m meters along each piece, measured on the WGS84 ellipsoid (the track's own
                vertices are never used directly). Pieces shorter than 3 spacings are dropped.
  3. Elevation  FABDEM at each point, bilinear
  4. Smooth     centered moving average of elevation over smooth_points points
  5. Grade      rise over run to the next point, in percent, positive uphill in the hiking direction. Points are
                evenly spaced, so every grade stands for the same distance: a share of points is a share of distance.

Writes data/grades/:
  <trail>.csv.gz   one row per point at the main settings: piece, km along the trail (gaps between pieces not
                   counted), lon, lat, elev_m (smoothed), grade_pct
  summary.csv      one row per trail and setting (the main one and every check): lengths and grade shares
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import shapely
import yaml
from scipy.ndimage import map_coordinates

from get_data import read_track_file, tile_name

ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "trails.yaml").read_text())
SETTINGS = CONFIG["settings"]
OSM_DIR = ROOT / "data" / "osm"
FABDEM_DIR = ROOT / "data" / "fabdem"
OUT = ROOT / "data" / "grades"
A, E2 = 6378137.0, 6.69437999014e-3          # WGS84 semi-major axis (m) and eccentricity squared


# ---- 1. Track -------------------------------------------------------------------------------------------------------

def load_ways(key):
    """Track lines for a trail, each an (n, 2) lon/lat array. OSM: the relation's ways, following section relations
    and members whose role is in keep_roles (alternates and spurs excluded); water and ferry ways dropped.
    track_file: the file's lines as recorded."""
    if "track_file" in CONFIG["trails"][key]:
        return read_track_file(CONFIG["trails"][key]["track_file"])
    rels, ways = {}, {}
    for f in sorted((OSM_DIR / key).glob("rel_*.json")):
        for e in json.loads(f.read_bytes())["elements"]:
            (rels if e["type"] == "relation" else ways)[e["id"]] = e
    tops = CONFIG["trails"][key]["osm_relation"]                  # one relation, or a list of them
    tops = tops if isinstance(tops, list) else [tops]
    if any(top not in rels for top in tops):
        sys.exit(f"No track for {key} in {OSM_DIR / key}: run scripts/get_data.py {key} first")
    keep, todo, seen = set(), list(tops), set(CONFIG["trails"][key].get("skip_relations", []))
    while todo:
        rel = todo.pop()
        if rel in seen:
            continue
        seen.add(rel)
        for m in rels[rel]["members"] if rel in rels else []:
            if m["role"] not in SETTINGS["keep_roles"]:
                continue
            if m["type"] == "way":
                keep.add(m["ref"])
            elif m["type"] == "relation" and m["ref"] not in seen:
                todo.append(m["ref"])
    walked = lambda tags: "waterway" not in tags and tags.get("route") != "ferry"   # river and ferry legs aren't walked
    return [np.array([(p["lon"], p["lat"]) for p in ways[w]["geometry"]])
            for w in sorted(keep) if w in ways and len(ways[w]["geometry"]) >= 2 and walked(ways[w].get("tags", {}))]


def join_pieces(ways):
    """Join ways that share end points into continuous pieces (on lon/lat, which is fine for topology)."""
    merged = shapely.line_merge(shapely.MultiLineString(ways))
    return [np.asarray(p.coords) for p in getattr(merged, "geoms", [merged])]


def orient_pieces(pieces, start, loop=None, end=None):
    """Put the pieces in trail order and point each the way hikers go: the order (and direction of every piece)
    that makes the shortest total jump from the start point, through each piece, to the end point. This works for
    trails that wind, where no single compass direction would. A loop is then turned to its usual direction."""
    out = chain_pieces(pieces, start, end)
    if loop:
        xy = np.concatenate(out)
        x = xy[:, 0] * np.cos(np.radians(xy[:, 1].mean()))
        counterclockwise = np.sum(x * np.roll(xy[:, 1], -1) - np.roll(x, -1) * xy[:, 1]) > 0   # shoelace area
        if counterclockwise != (loop == "counterclockwise"):
            out = [p[::-1] for p in out[::-1]]
        if len(out) == 1 and np.allclose(out[0][0], out[0][-1]):   # one closed ring: start it nearest `start`
            ring = out[0][:-1]
            i = int(np.argmin(np.hypot((ring[:, 0] - start[0]) * np.cos(np.radians(start[1])), ring[:, 1] - start[1])))
            out = [np.vstack([ring[i:], ring[:i], ring[i:i + 1]])]
    return out


def chain_pieces(pieces, start, end=None):
    """Order and direction of the pieces that make the shortest total jump from `start`, end to start through every
    piece, to `end` (if given): a nearest-end chain from the start, then improved by moving single pieces (either
    way round) and reversing runs of pieces until nothing shortens it. Distances are on a flat lon/lat grid scaled
    by latitude, plenty for telling near from far."""
    scale = np.cos(np.radians(start[1]))
    start = np.asarray(start, dtype=float)
    end = None if end is None else np.asarray(end, dtype=float)

    def gap(a, b):
        return 0.0 if a is None or b is None else float(np.hypot((a[0] - b[0]) * scale, a[1] - b[1]))

    # Nearest-end chain from the start. seq holds (piece, first point, last point) in trail order.
    left, seq, here = list(pieces), [], start
    while left:
        i = int(np.argmin([min(gap(here, p[0]), gap(here, p[-1])) for p in left]))
        p = left.pop(i)
        p = p if gap(here, p[0]) <= gap(here, p[-1]) else p[::-1]
        seq.append((p, p[0], p[-1]))
        here = p[-1]

    def before(s, i):                      # the point the trail arrives from at position i
        return start if i == 0 else s[i - 1][2]

    def after(s, i):                       # the point it goes on to after position i
        return end if i == len(s) - 1 else s[i + 1][1]

    improved = True
    while improved and len(seq) > 1:
        improved = False
        # Move one piece to another place in the order, either way round
        for i in range(len(seq)):
            p, a0, a1 = seq[i]
            saved = gap(before(seq, i), a0) + gap(a1, after(seq, i)) - gap(before(seq, i), after(seq, i))
            rest = seq[:i] + seq[i + 1:]
            for j in range(len(rest) + 1):
                x = start if j == 0 else rest[j - 1][2]
                y = end if j == len(rest) else rest[j][1]
                for first, last, piece in ((a0, a1, p), (a1, a0, p[::-1])):
                    if gap(x, first) + gap(last, y) - gap(x, y) < saved - 1e-9:
                        seq, improved = rest[:j] + [(piece, first, last)] + rest[j:], True
                        break
                if improved:
                    break
            if improved:
                break
        if improved:
            continue
        # Reverse a run of pieces (order and direction)
        for i in range(len(seq)):
            for j in range(i + 1, len(seq)):
                a, b = before(seq, i), after(seq, j)
                if gap(a, seq[j][2]) + gap(seq[i][1], b) < gap(a, seq[i][1]) + gap(seq[j][2], b) - 1e-9:
                    seq = seq[:i] + [(p[::-1], l, f) for p, f, l in seq[i:j + 1][::-1]] + seq[j + 1:]
                    improved = True
                    break
            if improved:
                break
    return [p for p, _, _ in seq]


# ---- 2. Points every spacing_m, on the ellipsoid --------------------------------------------------------------------

def step_m(lon, lat):
    """Length of each step between consecutive points, from the WGS84 radii of curvature at that latitude (accurate
    to well under 0.1% for steps of a few meters to a few hundred meters)."""
    phi = np.radians((lat[:-1] + lat[1:]) / 2)
    w = np.sqrt(1 - E2 * np.sin(phi) ** 2)
    dy = np.radians(np.diff(lat)) * A * (1 - E2) / w ** 3       # meridian radius
    dx = np.radians(np.diff(lon)) * A / w * np.cos(phi)         # prime vertical radius
    return np.hypot(dx, dy)


def resample(xy, spacing):
    """Points every `spacing` m along one piece (straight between its vertices): distance along it, lon, lat."""
    d = np.concatenate([[0.0], np.cumsum(step_m(xy[:, 0], xy[:, 1]))])
    keep = np.concatenate([[True], np.diff(d) > 0])
    at = np.arange(0, d[-1], spacing)
    return at, np.interp(at, d[keep], xy[keep, 0]), np.interp(at, d[keep], xy[keep, 1])


# ---- 3-5. Elevation, smoothing, grade -------------------------------------------------------------------------------

def fabdem(lon, lat):
    """FABDEM at each point, bilinear between the four nearest cell centers. Points within half a cell (15 m) of a
    tile's edge use the edge cells rather than the next tile's."""
    z = np.full(len(lon), np.nan)
    tile_lat, tile_lon = np.floor(lat).astype(int), np.floor(lon).astype(int)
    for la, lo in sorted(set(zip(tile_lat, tile_lon))):
        here = (tile_lat == la) & (tile_lon == lo)
        path = FABDEM_DIR / tile_name(la, lo)
        if not path.exists():
            print(f"  missing {path.name}: {here.sum()} points get no elevation")
            continue
        with rasterio.open(path) as src:
            band = src.read(1).astype(float)
            if src.nodata is not None:
                band[band == src.nodata] = np.nan
            col, row = ~src.transform * (lon[here], lat[here])
        z[here] = map_coordinates(band, [row - 0.5, col - 0.5], order=1, mode="nearest")
    return z


def smooth(z, k):
    """Centered moving average over k points (k odd). The first and last k // 2 points have no average."""
    out = np.full(len(z), np.nan)
    if k <= 1:
        return z.copy()
    if len(z) >= k:
        out[k // 2:len(z) - k // 2] = np.convolve(z, np.ones(k) / k, mode="valid")
    return out


def summarize(grade, spacing, k, track_km, official_km):
    a = np.abs(grade)
    return {
        "spacing_m": spacing, "smooth_points": k,
        "track_km": round(track_km, 1), "official_km": official_km,
        "track_vs_official": round(track_km / official_km, 3) if official_km else None,
        "km_graded": round(len(a) * spacing / 1000, 1),
        "median_abs_grade": round(float(np.median(a)), 2), "mean_abs_grade": round(float(a.mean()), 2),
        "pct_within_5": round(100 * float(np.mean(a <= 5)), 1),     # between -5% and 5% grade
        "pct_within_10": round(100 * float(np.mean(a <= 10)), 1),   # between -10% and 10% grade
        "pct_over_10": round(100 * float(np.mean(a > 10)), 1),
        "pct_over_15": round(100 * float(np.mean(a > 15)), 1),
        "pct_over_20": round(100 * float(np.mean(a > 20)), 1),
        "p95_abs_grade": round(float(np.percentile(a, 95)), 1),
    }


def build_trail(key):
    t = CONFIG["trails"][key]
    ways = load_ways(key)
    start = t.get("start", ways[0][0].tolist())                    # own track without a start: its first point
    pieces = orient_pieces(join_pieces(ways), start, t.get("loop"), t.get("end"))
    lengths = [step_m(p[:, 0], p[:, 1]).sum() for p in pieces]
    offsets = np.concatenate([[0.0], np.cumsum(lengths)[:-1]])     # where each piece starts along the trail
    track_km = sum(lengths) / 1000
    print(f"{t['name']}: {len(pieces)} pieces, {track_km:,.0f} km")
    main = (SETTINGS["spacing_m"], SETTINGS["smooth_points"])
    spacings = sorted(set(SETTINGS["checks"]["spacing_m"]) | {main[0]})
    smooths = sorted(set(SETTINGS["checks"]["smooth_points"]) | {main[1]})
    rows, points = [], None
    for spacing in spacings:
        parts = [(i, *resample(p, spacing)) for i, p in enumerate(pieces, 1) if lengths[i - 1] >= 3 * spacing]
        piece = np.concatenate([np.full(len(d), i) for i, d, *_ in parts])
        dist = np.concatenate([d + offsets[i - 1] for i, d, *_ in parts])
        lon = np.concatenate([x for *_, x, _ in parts])
        lat = np.concatenate([y for *_, y in parts])
        raw = fabdem(lon, lat)
        bounds = np.flatnonzero(np.diff(piece)) + 1                 # smooth and grade within each piece only
        for k in smooths:
            elev = np.concatenate([smooth(z, k) for z in np.split(raw, bounds)])
            grade = np.concatenate([np.append(np.diff(z) / spacing * 100, np.nan) for z in np.split(elev, bounds)])
            ok = ~np.isnan(grade)
            rows.append({"trail": key, "name": t["name"], "pieces": len(parts),
                         **summarize(grade[ok], spacing, k, track_km, t.get("official_km"))})
            if (spacing, k) == main:
                points = pd.DataFrame({"piece": piece, "km": (dist / 1000).round(3), "lon": lon.round(5),
                                       "lat": lat.round(5), "elev_m": elev.round(1), "grade_pct": grade.round(2)})[ok]
    return points, rows


def main():
    keys = sys.argv[1:] or list(CONFIG["trails"])
    OUT.mkdir(parents=True, exist_ok=True)
    summary_path = OUT / "summary.csv"
    old = pd.read_csv(summary_path) if summary_path.exists() else pd.DataFrame(columns=["trail"])
    rows = []
    for key in keys:
        points, trail_rows = build_trail(key)
        points.to_csv(OUT / f"{key}.csv.gz", index=False)
        rows += trail_rows
    summary = pd.concat([old[~old.trail.isin(keys)], pd.DataFrame(rows)], ignore_index=True)
    summary = summary.astype({"spacing_m": int, "smooth_points": int, "pieces": int})
    summary.to_csv(summary_path, index=False)
    main_rows = summary[(summary.spacing_m == SETTINGS["spacing_m"]) & (summary.smooth_points == SETTINGS["smooth_points"])]
    print(main_rows[["name", "track_km", "official_km", "pct_within_5", "pct_within_10", "pct_over_10",
                     "pct_over_15", "pct_over_20"]].to_string(index=False))


if __name__ == "__main__":
    main()
