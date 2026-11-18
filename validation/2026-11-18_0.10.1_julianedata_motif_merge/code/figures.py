#!/usr/bin/env python3
"""Figures for the independent Juliane-data benchmark.

The palette, the mode names and the dash contract are imported from the
2026-08-06 study so the two reports can be read side by side: colour is the
pipeline, a dash means the arm was handed the true allele count K, and
MixMHCp is red. GibbsCluster was not run here, so blue is absent.

  fig1_ami_alleles     AMI and BCubed F1 against the number of alleles
  fig2_metrics         five metrics against the number of alleles
  fig3_granularity     recovered motif count against the true allele count
  fig4_cost            median wall-clock against pool size
  fig5_imbalance       AMI against the per-allele imbalance regime
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / "2026-08-06_0.6.0_motif_merge" / "code"))
from figures import MODES, STYLE, METRICS  # noqa: E402

ORDER = [label for _, label in MODES] + ["MixMHCp", "MixMHCp (given K)"]
BALANCE_ORDER = ["eq", "a10", "a3", "a1", "a03"]
BALANCE_LABEL = {"eq": "equal", "a10": "mild\n(alpha 10)",
                 "a3": "moderate\n(alpha 3)", "a1": "strong\n(alpha 1)",
                 "a03": "extreme\n(alpha 0.3)"}

plt.rcParams.update({"font.size": 13, "axes.titlesize": 15,
                     "figure.titlesize": 16, "savefig.dpi": 160,
                     "savefig.bbox": "tight"})


def load(results: Path, manifest: Path) -> pd.DataFrame:
    frames = []
    for key, label in MODES:
        d = pd.read_csv(results / f"variant_{key}.csv")
        d["tool"] = label
        frames.append(d)
    m = pd.read_csv(results / "mixmhcp.csv")
    # "MixMHCp" alone is the default arm, matching the other study's legend.
    m["tool"] = m["setting"].map({"default": "MixMHCp",
                                  "oracle_k": "MixMHCp (given K)"})
    frames.append(m)
    keep = (["tool", "pool", "allele_count", "peptides", "clusters",
             "elapsed_seconds"] + [k for k, _ in METRICS])
    df = pd.concat([f[keep] for f in frames], ignore_index=True)
    man = pd.read_csv(manifest)[["pool", "balance"]]
    return df.merge(man, on="pool", how="left")


def _allele_axis(ax, df):
    """A log x-axis labelled only at the five allele counts.

    A log scale otherwise draws its own minor ticks (2, 4, 20, 40 ...), whose
    labels collide with the five we actually want.
    """
    ax.set_xscale("log")
    ax.set_xticks(sorted(df.allele_count.unique()))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_xaxis().set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("alleles in pool")


def _legend(fig, ax):
    h, l = ax.get_legend_handles_labels()
    fig.legend(h, l, loc="outside right center", frameon=False)


def _series(ax, df, x, y, agg="mean"):
    for tool in ORDER:
        s = df[df.tool == tool]
        if s.empty:
            continue
        g = s.groupby(x)[y].agg(agg)
        ax.plot(g.index, g.values, label=tool, **STYLE[tool])


def fig1(df, out):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for ax, (metric, name) in zip(axes, [("ami", "AMI"),
                                         ("bcubed_f1_macro", "BCubed F1")]):
        _series(ax, df, "allele_count", metric)
        _allele_axis(ax, df)
        ax.set_ylabel(name)
        ax.set_title(name)
        ax.grid(alpha=0.25)
    fig.suptitle("Independent benchmark: quality against pool complexity")
    _legend(fig, axes[0])
    fig.savefig(out / "fig1_ami_alleles.png")
    plt.close(fig)


def fig2(df, out):
    fig, axes = plt.subplots(5, 1, figsize=(10, 19), layout="constrained",
                             sharex=True)
    for ax, (metric, name) in zip(axes, METRICS):
        _series(ax, df, "allele_count", metric)
        ax.set_ylabel(name)
        ax.set_title(name)
        ax.grid(alpha=0.25)
    _allele_axis(axes[-1], df)
    fig.suptitle("Independent benchmark: every metric against pool complexity")
    _legend(fig, axes[0])
    fig.savefig(out / "fig2_metrics.png")
    plt.close(fig)


def fig3(df, out):
    fig, ax = plt.subplots(figsize=(8.5, 6), layout="constrained")
    _series(ax, df, "allele_count", "clusters")
    lo, hi = 3, 50
    ax.plot([lo, hi], [lo, hi], color="#999999", ls=":", lw=1.5, zorder=0,
            label="one motif per allele")
    _allele_axis(ax, df)
    ax.set_yscale("log")
    ax.set_ylabel("motifs reported (mean)")
    ax.set_title("Granularity: motifs reported against alleles present")
    ax.grid(alpha=0.25, which="both")
    _legend(fig, ax)
    fig.savefig(out / "fig3_granularity.png")
    plt.close(fig)


def fig4(df, out):
    """Cost against pool size.

    Pool sizes are continuous here, not a fixed grid, so the pools are binned
    by decade-ish edges and the median taken within each bin; the band is the
    interquartile range, which is what varies between pools of a similar size.
    """
    edges = [900, 2000, 4000, 8000, 16000, 32000, 64000, 100000]
    d = df.copy()
    d["bin"] = pd.cut(d.peptides, edges)
    centre = {iv: np.sqrt(iv.left * iv.right) for iv in d.bin.cat.categories}
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    for tool in ORDER:
        s = d[d.tool == tool]
        if s.empty:
            continue
        g = s.groupby("bin", observed=True)["elapsed_seconds"]
        x = [centre[iv] for iv in g.median().index]
        ax.plot(x, g.median().values, label=tool, **STYLE[tool])
        ax.fill_between(x, g.quantile(0.25).values, g.quantile(0.75).values,
                        color=STYLE[tool]["color"], alpha=0.12, lw=0, zorder=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # Label the bin edges themselves; the log scale's own minor labels
    # (2x10^3, 4x10^3 ...) overlap at this width.
    ax.set_xticks(edges[1:-1])
    ax.set_xticklabels([f"{e//1000}k" for e in edges[1:-1]])
    ax.get_xaxis().set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel("peptides in pool")
    ax.set_ylabel("wall-clock seconds (median, IQR band)")
    ax.set_title("Cost against pool size")
    ax.grid(alpha=0.25, which="both")
    _legend(fig, ax)
    fig.savefig(out / "fig4_cost.png")
    plt.close(fig)


def fig5(df, out):
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    xs = np.arange(len(BALANCE_ORDER))
    for tool in ORDER:
        s = df[df.tool == tool]
        if s.empty:
            continue
        g = s.groupby("balance")["ami"].mean().reindex(BALANCE_ORDER)
        ax.plot(xs, g.values, label=tool, **STYLE[tool])
    ax.set_xticks(xs)
    ax.set_xticklabels([BALANCE_LABEL[b] for b in BALANCE_ORDER])
    ax.set_xlabel("per-allele abundance imbalance")
    ax.set_ylabel("AMI")
    ax.set_title("Quality against abundance imbalance")
    ax.grid(alpha=0.25, axis="y")
    _legend(fig, ax)
    fig.savefig(out / "fig5_imbalance.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--manifest", type=Path, default=Path("pool_manifest.csv"))
    ap.add_argument("--out", type=Path, default=Path("figures"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    df = load(args.results, args.manifest)
    for fn in (fig1, fig2, fig3, fig4, fig5):
        fn(df, args.out)
        print(f"  {fn.__name__} ok")


main()
