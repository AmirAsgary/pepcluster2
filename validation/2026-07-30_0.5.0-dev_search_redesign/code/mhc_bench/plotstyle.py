#!/usr/bin/env python3
"""One visual language for every figure in this study.

Colours are assigned from a fixed order and never cycled, so a given series is
the same colour in every plot: blue first, red second, then green, orange,
yellow. Blue and red lead because they are the two most reliably distinguishable
hues; from the third series onward colour alone is not enough, so callers should
also vary line style or marker.
"""

from __future__ import annotations

SERIES = ("#2a78d6", "#d1342f", "#1baf7a", "#eb6834", "#eda100",
          "#7e3ace", "#e87ba4", "#0b8f8f")
INK = "#52514e"
HEADING = "#0b0b0b"
GRID = "#c9c8c3"
FONT = "DejaVu Sans"

# Fixed assignment, so a mode keeps its colour across every figure.
MODE_COLOR = {"separate_kmer_anchor": SERIES[0], "separate_aln_anchor": SERIES[1]}
MODE_LABEL = {"separate_kmer_anchor": "k-mer + anchor",
              "separate_aln_anchor": "alignment + anchor"}
SPLIT_LABEL = {"inner": "tuning folds", "outer": "held-out alleles", "test": "benchmark"}


def apply() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": FONT,
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "axes.titlecolor": HEADING,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.5,
        "legend.frameon": False,
        "figure.titlesize": 11.5,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "axes.edgecolor": GRID,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.0,
        "lines.markersize": 4.5,
    })


def finish(axis, xlabel: str = "", ylabel: str = "", title: str = "") -> None:
    if title:
        axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
