#!/usr/bin/env python3
"""Download each trail's track and the elevation tiles it crosses.

    .venv/bin/python scripts/get_data.py              # every trail in trails.yaml
    .venv/bin/python scripts/get_data.py tmb jmt      # only these trails
    .venv/bin/python scripts/get_data.py --tracks     # tracks only, no elevation tiles

Sources:
  tracks     OpenStreetMap route relations, main OSM API (ODbL). A trail with `track_file` in trails.yaml uses that
             GPX or GeoJSON file instead; nothing is downloaded for its track.
  elevation  FABDEM v1-2 (University of Bristol), 30 m bare-earth DEM (CC BY-NC-SA 4.0)

Output:
  data/osm/<trail>/rel_<id>.json   one file per route relation: the relation (member roles) and its ways (geometry)
  data/fabdem/<tile>.tif           FABDEM 1 x 1 degree tiles crossed by the tracks, ~25 MB each (not in git)

Re-running skips anything already downloaded. FABDEM is published as 10 x 10 degree zips of 1-2 GB; the server
accepts byte-range requests, so this reads each zip's directory and pulls out only the tiles it needs.
"""

import io
import json
import math
import sys
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "trails.yaml").read_text())
OSM_DIR = ROOT / "data" / "osm"
FABDEM_DIR = ROOT / "data" / "fabdem"
OSM_API = "https://api.openstreetmap.org/api/0.6/relation/{rel}{full}.json"
FABDEM_URL = "https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/{zip}"
USER_AGENT = "comparing-trail-grade-distributions/0.1"
EDGE_DEG = 0.01            # a point this close to a tile edge also gets the neighboring tile
WORKERS = 4                # FABDEM tiles downloaded at once (the server gives about 0.5 MB/s per connection)


def log(msg):
    print(msg, flush=True)


def fetch(url, headers=None, tries=5):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as e:
            if attempt == tries - 1:
                raise
            log(f"    {type(e).__name__} on {url}; retrying")
            time.sleep(10 * (attempt + 1))


# ---- Tracks --------------------------------------------------------------------------------------------------------

def osm_relation(rel, full=True):
    """One relation from the OSM API. With full=True, also its member ways with their node coordinates. Returned as
    {"elements": [relation, ways...]}, each way with a "geometry" list of {lat, lon}."""
    url = OSM_API.format(rel=rel, full="/full" if full else "")
    els = json.loads(fetch(url))["elements"]
    nodes = {e["id"]: {"lat": e["lat"], "lon": e["lon"]} for e in els if e["type"] == "node"}
    out = []
    for e in els:
        if e["type"] == "relation":
            out.append({"type": "relation", "id": e["id"], "tags": e.get("tags", {}),
                        "members": [{"type": m["type"], "ref": m["ref"], "role": m.get("role", "")}
                                    for m in e["members"]]})
        elif e["type"] == "way":
            out.append({"type": "way", "id": e["id"], "tags": e.get("tags", {}),
                        "geometry": [nodes[n] for n in e["nodes"] if n in nodes]})
    return {"source": url, "elements": out}


def main_route_relations(rel, keep_roles):
    """The relation and, for a superroute, every section relation under it on the main route (sections with another
    role, such as "alternative", are left out), depth first."""
    top = osm_relation(rel, full=False)["elements"][0]
    found = [rel]
    for m in top["members"]:
        if m["type"] == "relation" and m["role"] in keep_roles:
            time.sleep(1)
            found += main_route_relations(m["ref"], keep_roles)
    return found


def read_track_file(path):
    """A GPX (tracks or routes) or GeoJSON (LineString, MultiLineString) file -> list of (n, 2) lon/lat arrays."""
    import numpy as np
    path = ROOT / path
    if path.suffix.lower() == ".gpx":
        import xml.etree.ElementTree as ET
        tag = lambda e: e.tag.rsplit("}", 1)[-1]                   # drop the XML namespace
        lines = []
        for seg in ET.parse(path).iter():
            if tag(seg) in ("trkseg", "rte"):
                pts = [(float(p.get("lon")), float(p.get("lat"))) for p in seg if tag(p) in ("trkpt", "rtept")]
                if len(pts) >= 2:
                    lines.append(np.array(pts))
        return lines
    geo = json.loads(path.read_text())
    feats = geo.get("features", [geo])
    lines = []
    for f in feats:
        g = f.get("geometry", f)
        parts = [g["coordinates"]] if g["type"] == "LineString" else g["coordinates"] if g["type"] == "MultiLineString" else []
        lines += [np.array([c[:2] for c in part]) for part in parts if len(part) >= 2]
    return lines


def get_tracks(keys):
    keep_roles = CONFIG["settings"]["keep_roles"]
    for key in keys:
        if "track_file" in CONFIG["trails"][key]:
            continue                                               # own track: nothing to download
        out_dir = OSM_DIR / key
        out_dir.mkdir(parents=True, exist_ok=True)
        tops = CONFIG["trails"][key]["osm_relation"]              # one relation, or a list of them
        rels = [r for top in (tops if isinstance(tops, list) else [tops]) for r in main_route_relations(top, keep_roles)]
        log(f"{key}: {len(rels)} relation(s)")
        for rel in rels:
            out = out_dir / f"rel_{rel}.json"
            if out.exists():
                continue
            body = json.dumps(osm_relation(rel)).encode()
            out.write_bytes(body)
            log(f"  rel {rel}: {len(body) / 1e6:.1f} MB")
            time.sleep(1)                        # the OSM API asks for gentle use


# ---- Elevation -----------------------------------------------------------------------------------------------------

def track_points(keys):
    for key in keys:
        if "track_file" in CONFIG["trails"][key]:
            for line in read_track_file(CONFIG["trails"][key]["track_file"]):
                yield from ((lat, lon) for lon, lat in line)
            continue
        for f in sorted((OSM_DIR / key).glob("*.json")):
            for el in json.loads(f.read_bytes())["elements"]:
                yield from ((p["lat"], p["lon"]) for p in el.get("geometry", []))


def tile_name(lat, lon):
    return f"{'N' if lat >= 0 else 'S'}{abs(lat):02d}{'E' if lon >= 0 else 'W'}{abs(lon):03d}_FABDEM_V1-2.tif"


def zip_name(lat, lon):
    la, lo = math.floor(lat / 10) * 10, math.floor(lon / 10) * 10
    corner = lambda a, b: f"{'N' if a >= 0 else 'S'}{abs(a):02d}{'E' if b >= 0 else 'W'}{abs(b):03d}"
    east = lo + 10 if lo + 10 < 180 else lo + 10 - 360         # FABDEM names the 180th meridian W180
    return f"{corner(la, lo)}-{corner(la + 10, east)}_FABDEM_V1-2.zip"


class HttpRangeFile(io.RawIOBase):
    """A read-only, seekable file over HTTP byte ranges, so zipfile can read one member of a remote zip."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Length"])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = {io.SEEK_SET: offset, io.SEEK_CUR: self.pos + offset, io.SEEK_END: self.size + offset}[whence]
        return self.pos

    def readinto(self, b):
        if self.pos >= self.size:
            return 0
        end = min(self.pos + len(b), self.size) - 1
        data = fetch(self.url, {"Range": f"bytes={self.pos}-{end}"})
        b[:len(data)] = data
        self.pos += len(data)
        return len(data)


def get_tile(zname, tile):
    remote = io.BufferedReader(HttpRangeFile(FABDEM_URL.format(zip=zname)), buffer_size=8 << 20)
    with zipfile.ZipFile(remote) as z:
        members = {Path(n).name: n for n in z.namelist()}
        if tile not in members:
            return f"{tile}: not in {zname} (no land there)"
        part = FABDEM_DIR / (tile + ".part")
        part.write_bytes(z.read(members[tile]))
        part.rename(FABDEM_DIR / tile)
    return f"{tile}: {(FABDEM_DIR / tile).stat().st_size / 1e6:.0f} MB"


def get_fabdem(keys):
    FABDEM_DIR.mkdir(parents=True, exist_ok=True)
    need = {(math.floor(lat + dlat), math.floor(lon + dlon))
            for lat, lon in track_points(keys)
            for dlat in (-EDGE_DEG, 0, EDGE_DEG) for dlon in (-EDGE_DEG, 0, EDGE_DEG)}
    todo = [(zip_name(*t), tile_name(*t)) for t in sorted(need) if not (FABDEM_DIR / tile_name(*t)).exists()]
    log(f"FABDEM: {len(need)} tiles, {len(todo)} to download")
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for n, msg in enumerate(pool.map(lambda zt: get_tile(*zt), todo), 1):
            log(f"  [{n}/{len(todo)}] {msg}")


if __name__ == "__main__":
    keys = [a for a in sys.argv[1:] if not a.startswith("--")] or list(CONFIG["trails"])
    get_tracks(keys)
    if "--tracks" not in sys.argv:
        get_fabdem(keys)
