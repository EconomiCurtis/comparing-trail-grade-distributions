#!/usr/bin/env python3
"""Grade-distribution figures from data/grades/ (scripts/build_grades.py).

    .venv/bin/python scripts/plots.py                      # density: pct at cdt annapurna; heatmaps, bars: all
    .venv/bin/python scripts/plots.py --density pct tmb    # pick the density plot's trails

Writes figures/:
  density.png    grade density (Gaussian KDE) per trail over the grade zones
  heatmaps.png   grade x elevation, 1% x 250 ft cells, hue = zone, lightness = log count; 10th/90th percentile lines
  zone-bars.png  share of distance per grade zone, one bar per trail, sorted by share within ±10%
"""

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.stats import gaussian_kde

ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "trails.yaml").read_text())
GRADES = ROOT / "data" / "grades"
FIGURES = ROOT / "figures"
FT = 3.28084
KDE_BW = 0.18                      # KDE bandwidth, x each trail's grade sd
CELL_GRADE, CELL_FT = 1, 250       # heatmap cells
BAND_FT, MIN_IN_BAND = 500, 30     # percentile lines: per 500 ft band with >= 30 points
LINE_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]   # Okabe-Ito; reference trail (first) is a gray fill

# Grade zones (edges from trails.yaml: 3, 7, 10, 16, 25 -> 11 zones). Bar / band color, and (hue, chroma) at fixed
# CIE lightness for the heatmap, so cell lightness encodes count only.
ZONE_STYLE = [
    ("<-25%", "#7F0000", (15, 45)), ("-16 to -25%", "#D7301F", (40, 70)), ("-10 to -16%", "#3182BD", (290, 70)),
    ("-7 to -10%", "#B9E5C2", (155, 22)), ("-3 to -7%", "#439543", (140, 75)), ("±3%", "#77C180", (148, 46)),
    ("3 to 7%", "#FDB813", (95, 70)), ("7 to 10%", "#FD8D3C", (75, 70)), ("10 to 16%", "#E6550D", (58, 70)),
    ("16 to 25%", "#D7301F", (40, 70)), (">25%", "#7F0000", (15, 45)),
]
EDGES = np.array([-e for e in reversed(CONFIG["zones"])] + CONFIG["zones"], dtype=float)


def lch_rgb(L, C, h):
    """CIE LCh (D65) -> sRGB 0-1; chroma reduced until in gamut."""
    for c in np.arange(C, -0.1, -1.0):
        a, b = c * np.cos(np.radians(h)), c * np.sin(np.radians(h))
        fy = (L + 16) / 116
        f = np.array([fy + a / 500, fy, fy - b / 200])
        xyz = np.where(f ** 3 > 0.008856, f ** 3, (f - 16 / 116) / 7.787) * np.array([0.95047, 1.0, 1.08883])
        rgb = np.array([[3.2404542, -1.5371385, -0.4985314], [-0.9692660, 1.8760108, 0.0415560],
                        [0.0556434, -0.2040259, 1.0572252]]) @ xyz
        if rgb.min() >= -1e-4 and rgb.max() <= 1 + 1e-4:
            rgb = np.clip(rgb, 0, 1)
            return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
    return np.full(3, L / 100)


def load(keys=None):
    """{trail key: 100 m points} for trails with grades built."""
    keys = keys or list(CONFIG["trails"])
    return {k: pd.read_csv(GRADES / f"{k}.csv.gz") for k in keys if (GRADES / f"{k}.csv.gz").exists()}


def name(key):
    return CONFIG["trails"][key]["name"]


def zone_shares(grade):
    """% of points in each grade zone."""
    return np.bincount(np.digitize(grade, EDGES), minlength=len(ZONE_STYLE)) / len(grade) * 100


def density(ax, data):
    """KDE of grade per trail; first trail filled gray, others lines. Zone bands behind."""
    x = np.arange(-40, 40.01, 0.25)
    bounds = np.concatenate([[-40], EDGES, [40]])
    for (label, color, _), lo, hi in zip(ZONE_STYLE, bounds[:-1], bounds[1:]):
        ax.axvspan(lo, hi, color=color, alpha=0.18, lw=0)
    for i, (key, d) in enumerate(data.items()):
        y = gaussian_kde(d.grade_pct, bw_method=KDE_BW)(x) * 100        # % of distance per 1% grade
        if i == 0:
            ax.fill_between(x, y, color="#A6A6A6", alpha=0.9, label=name(key), zorder=2)
        else:
            ax.plot(x, y, color=LINE_COLORS[(i - 1) % 4], lw=2, label=name(key), zorder=3)
    ax.axvline(0, color="black", lw=0.6, alpha=0.5)
    ax.set(xlim=(-40, 40), ylim=(0, None), xlabel="Grade (%): downhill < 0 < uphill", ylabel="% of distance per 1% grade")
    ax.legend(frameon=False, loc="upper left")


def heatmap(ax, d, y_max_ft, steps=16):
    """Grade x elevation cells: hue by zone (fixed lightness), lightness by log count (trail's own max)."""
    xe = np.arange(-40, 40 + CELL_GRADE, CELL_GRADE)
    ye = np.arange(0, y_max_ft + CELL_FT, CELL_FT)
    n, *_ = np.histogram2d(d.elev_m * FT, d.grade_pct, bins=[ye, xe])
    zone = np.digitize((xe[:-1] + xe[1:]) / 2, EDGES)
    level = np.where(n > 0, np.log(np.maximum(n, 1)) / np.log(max(n.max(), 2)), 0)
    img = np.zeros(n.shape + (4,))
    for i, j in zip(*np.nonzero(n)):
        hue, chroma = ZONE_STYLE[zone[j]][2]
        img[i, j] = (*lch_rgb(95 - 45 * level[i, j], chroma, hue), 1)
    ax.imshow(img, origin="lower", extent=(xe[0], xe[-1], ye[0], ye[-1]), aspect="auto", interpolation="nearest")
    band = np.floor(d.elev_m * FT / BAND_FT)
    env = [(b * BAND_FT + BAND_FT / 2, *np.percentile(d.grade_pct[band == b], [10, 90]))
           for b in np.unique(band) if (band == b).sum() >= MIN_IN_BAND]
    if env:
        mid, p10, p90 = map(np.array, zip(*env))
        ax.plot(p10, mid, color="black", lw=1.3)
        ax.plot(p90, mid, color="black", lw=1.3)
    for e in EDGES:
        ax.axvline(e, color="black", lw=0.5, ls=":", alpha=0.4)
    ax.set(xlim=(-40, 40), ylim=(0, y_max_ft))


def zone_bars(ax, data):
    """Stacked share of distance per zone, one bar per trail, easiest (most within ±10%) at the top."""
    shares = {k: zone_shares(d.grade_pct) for k, d in data.items()}
    order = sorted(shares, key=lambda k: -np.mean(np.abs(data[k].grade_pct) <= 10))
    for row, key in enumerate(order):
        left = 0
        for (label, color, _), s in zip(ZONE_STYLE, shares[key]):
            ax.barh(row, s, left=left, color=color, edgecolor="white", lw=0.5, height=0.75,
                    label=label if row == 0 else None)
            left += s
    ax.set(yticks=range(len(order)), yticklabels=[name(k) for k in order], xlim=(0, 100), xlabel="% of trail distance")
    ax.invert_yaxis()
    ax.legend(ncol=6, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=8)


def main():
    args = sys.argv[1:]
    pick = args[args.index("--density") + 1:] if "--density" in args else ["pct", "at", "cdt", "annapurna"]
    FIGURES.mkdir(exist_ok=True)
    data = load()

    fig, ax = plt.subplots(figsize=(10, 5))
    density(ax, {k: data[k] for k in pick if k in data})
    fig.tight_layout()
    fig.savefig(FIGURES / "density.png", dpi=150)

    cols = 4
    rows = int(np.ceil(len(data) / cols))
    y_max = np.ceil(max(d.elev_m.max() for d in data.values()) * FT / 1000) * 1000
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.4 * rows), sharex=True, sharey=True, squeeze=False)
    for ax, (key, d) in zip(axes.flat, data.items()):
        heatmap(ax, d, y_max)
        ax.set_title(name(key), fontsize=10)
    for ax in axes.flat[len(data):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Grade (%)")
    for ax in axes[:, 0]:
        ax.set_ylabel("Elevation (ft)")
    fig.tight_layout()
    fig.savefig(FIGURES / "heatmaps.png", dpi=130)

    fig, ax = plt.subplots(figsize=(10, 0.45 * len(data) + 1.5))
    zone_bars(ax, data)
    fig.tight_layout()
    fig.savefig(FIGURES / "zone-bars.png", dpi=150)
    print("Saved", ", ".join(p.name for p in sorted(FIGURES.glob("*.png"))))


if __name__ == "__main__":
    matplotlib.use("Agg")                 # files only; imports (the notebook) keep their backend
    main()
