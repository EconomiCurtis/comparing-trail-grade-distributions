# Comparing Trail Grade Distributions

Trail grade every 100 m along long-distance hiking trails, computed the same way for every trail: OpenStreetMap
tracks (or your own GPX/GeoJSON), FABDEM bare-earth elevation, a 3-point moving average, rise over run. Twenty-one
trails are included; adding one is a few lines in `trails.yaml`.

Write-up: [How Hard Is the Trail Going to Be?](https://economicurtis.com) (link to come).

## Where things are

| Path | What it is |
| --- | --- |
| `trails.yaml` | Trails (OSM relation or track file, start/end, direction) and settings (spacing, smoothing, checks, grade zones) |
| `scripts/get_data.py` | Downloads OSM route relations and the FABDEM tiles a track crosses |
| `scripts/build_grades.py` | Track -> points every 100 m -> elevation -> smoothing -> grade; writes `data/grades/` |
| `scripts/plots.py` | Density, grade-by-elevation heatmaps, grade-zone bars; writes `figures/` |
| `notebook/walkthrough.qmd`, `.html` | Runs each step on one trail, then all trails and the figures |
| `data/osm/<trail>/` | OSM route relations as downloaded (relation member roles, way geometry) |
| `data/grades/<trail>.csv.gz` | One row per 100 m point: piece, km, lon, lat, elev_m (smoothed), grade_pct |
| `data/grades/summary.csv` | Per trail and setting (spacing 50/100/200 m x smoothing 1/3/5 points): lengths, grade shares |
| `data/fabdem/` | FABDEM tiles, downloaded by `make download` (not in the repo: ~25 MB each, 4 GB in all) |
| `figures/` | Output of `scripts/plots.py` |

The trails: PCT, AT, CDT, Great Divide Trail, JMT, Colorado Trail, Arizona Trail, Pacific Northwest Trail, Wonderland
Trail, Tahoe Rim Trail, Tour du Mont Blanc, Via Alpina Red and Green, Fra li Monti (GR20), Haute Randonnée Pyrénéenne
(HRP), Pyrenean Way (GR10), Senda Pirenaica (GR11), Annapurna Circuit, Te Araroa, Australian Alps Walking Track,
Camino de Santiago (Camino Francés).

## Run

Python 3.10+. The notebook also needs [Quarto](https://quarto.org).

```bash
make setup                 # .venv from requirements.txt
make download              # FABDEM tiles for every trail (~4 GB); OSM tracks are already in data/osm/
make grades figures        # data/grades/, figures/
make notebook              # notebook/walkthrough.html
```

The grades in `data/grades/` are already built, so `make figures` works without downloading anything. Any step takes
`TRAILS="tmb jmt"` to limit it to some trails.

## Add a trail

From OpenStreetMap: find the route relation (e.g. on [Waymarked Trails](https://hiking.waymarkedtrails.org)) and add
an entry to `trails.yaml` with `osm_relation`, `start` (sets the sign of grade), and `end` or `loop`. From your own
track: put a GPX or GeoJSON file in `data/tracks/` and give `track_file` instead of `osm_relation`. Then:

```bash
make download grades TRAILS=my_trail
make figures
```

Check the printed track length against the published length. Extra length or many pieces usually means alternates
mapped as main route; drop those section relations with `skip_relations`. The notebook's last section has examples.

## Method

1. **Track.** OSM route relation, main-route member roles only (`keep_roles`); water and ferry ways dropped. Ways
   joined end to end; pieces ordered and pointed from `start` to `end` (shortest total gap); loops turned to `loop`.
2. **Points.** Every `spacing_m` (100 m) along each piece, distances on the WGS84 ellipsoid. Track vertices unused.
3. **Elevation.** FABDEM v1-2 at each point, bilinear.
4. **Smooth.** Centered moving average over `smooth_points` (3) points.
5. **Grade.** Rise over run to the next point, %. Equal spacing: a share of points is a share of distance.

Every run repeats steps 2-5 at each spacing and smoothing in `settings.checks` (`summary.csv`), to show how much the
results depend on those choices.

## Licenses

Code (scripts, notebook): MIT, see `LICENSE`. The data keep their sources' licenses:

- Tracks: © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors, ODbL 1.0. `data/osm/` and the
  coordinates in `data/grades/` are OSM data and carry the ODbL.
- Elevation: [FABDEM v1-2](https://data.bris.ac.uk/data/dataset/s5hqmjcdj8yo2ibzi9b4ew3sn) (Hawker et al. 2022,
  University of Bristol), CC BY-NC-SA 4.0. The elevations and grades in `data/grades/` derive from it: non-commercial
  use, same license.
