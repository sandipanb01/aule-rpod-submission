"""Shared figure style: validated categorical palette, recessive chrome, thin marks."""
from __future__ import annotations

import matplotlib as mpl

# Categorical slots (validated reference palette; fixed order, never cycled)
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
# Status palette (reserved for state; always paired with icon/label)
GOOD, WARNING, SERIOUS, CRITICAL = "#0ca30c", "#fab219", "#ec835a", "#d03b3b"
# Chrome & ink
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SHADOW = "#383835"

STATION_COLORS = {"Svalbard": BLUE, "Bengaluru": ORANGE, "Wellington": AQUA, "Los Angeles": VIOLET}
SUN_COLOR = YELLOW
BAND_COLOR = MAGENTA


def apply():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK2,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "legend.fontsize": 8, "lines.linewidth": 1.4, "savefig.dpi": 200, "savefig.bbox": "tight",
        "text.color": INK,
    })
