#!/usr/bin/env python3
"""Cost against pool size, and cost against accuracy.

Two panels each. CPU seconds is the primary axis because it is the resource
actually consumed and is comparable across tools regardless of how each one
happens to be threaded; wall seconds is shown beside it because it is what a user
waits for. PepCluster2 and MixMHCp are single-threaded here, so for them the two
are the same number. GibbsCluster is run at -k 16, so its wall time is 16x
cheaper-looking than its cost.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

STYLE = {
    "pc2_cluster":  dict(color="#000000", ls=":",  marker="s",
                         label="PepCluster2, clustering only"),
    "pc2_motif":    dict(color="#000000", ls="-",  marker="o",
                         label="PepCluster2 + merge + EM"),
    "pc2_motif_k":  dict(color="#000000", ls=(0, (4, 1.5)), marker="D",
                         label="PepCluster2 + merge + EM, given k"),
    "mixmhcp":      dict(color="#d62728", ls="-",  marker="o", label="MixMHCp"),
    "gibbscluster": dict(color="#1f77b4", ls="-",  marker="o", label="GibbsCluster"),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--speed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})

    frame = pd.read_csv(args.speed)

    def decade(caxis):
        caxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0),
                                                  numticks=30))
        caxis.set_minor_locator(ticker.LogLocator(base=10.0,
                                                  subs=tuple(np.arange(1, 10) * 0.1),
                                                  numticks=100))
        caxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
        caxis.set_minor_formatter(ticker.NullFormatter())

    # ---- fig6: cost against pool size ---------------------------------
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for axis, column, title in ((axes[0], "cpu_seconds", "CPU seconds"),
                                (axes[1], "wall_seconds", "Wall seconds")):
        for arm, style in STYLE.items():
            sub = frame[frame.arm == arm].sort_values("peptides")
            if sub.empty:
                continue
            axis.plot(sub.peptides, sub[column], markersize=4, lw=1.6, **style)
        axis.set_xscale("log")
        axis.set_yscale("log")
        decade(axis.xaxis)
        decade(axis.yaxis)
        axis.set_xlabel("Peptides in pool")
        axis.set_ylabel(title)
        axis.set_title(title)
        axis.grid(alpha=0.25, lw=0.5, which="both")
    axes[0].legend(fontsize=7, loc="upper left")
    figure.suptitle("Cost against pool size (serial runs; PepCluster2 and MixMHCp "
                    "single-threaded, GibbsCluster at -k 16)")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig6_speed_vs_pool_size.{suffix}",
                       bbox_inches="tight")
    plt.close(figure)
    frame.to_csv(args.out / "fig6_speed_vs_pool_size.csv", index=False)

    # ---- fig7: cost against accuracy ----------------------------------
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    summary = []
    for axis, metric, title in ((axes[0], "bcubed_f1_macro", "BCubed F1"),
                                (axes[1], "ami", "AMI")):
        for arm, style in STYLE.items():
            sub = frame[frame.arm == arm]
            if sub.empty:
                continue
            x = sub.cpu_seconds.median()
            y = sub[metric].mean()
            err = sub[metric].std(ddof=1) / max(np.sqrt(len(sub)), 1)
            marker = dict(style)
            marker.pop("ls")
            axis.errorbar([x], [y], yerr=[err], markersize=11, capsize=3,
                          lw=0, elinewidth=1, **marker)
            if metric == "bcubed_f1_macro":
                summary.append(dict(arm=arm, label=style["label"],
                                    median_cpu_seconds=x,
                                    mean_f1=y, mean_ami=sub.ami.mean(),
                                    pools=len(sub)))
        axis.set_xscale("log")
        decade(axis.xaxis)
        axis.set_xlabel("Median CPU seconds per pool")
        axis.set_ylabel(title)
        axis.set_ylim(0, 1)
        axis.set_title(f"{title} against cost")
        axis.grid(alpha=0.25, lw=0.5, which="both")
    axes[0].legend(fontsize=7, loc="upper right")
    figure.suptitle("Accuracy against cost: up and to the left is better")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig7_speed_vs_performance.{suffix}",
                       bbox_inches="tight")
    plt.close(figure)
    table = pd.DataFrame(summary)
    table.to_csv(args.out / "fig7_speed_vs_performance.csv", index=False)

    print(table.round(4).to_string(index=False))
    print("\ncost relative to PepCluster2 + merge + EM:")
    base = table[table.arm == "pc2_motif"].median_cpu_seconds.iloc[0]
    for row in table.itertuples():
        print(f"  {row.label:42s} {row.median_cpu_seconds / base:7.2f}x")
    # Does cost scale differently with pool size across tools?
    print("\nSpearman rho of CPU seconds against pool size:")
    for arm in STYLE:
        sub = frame[frame.arm == arm]
        if len(sub) >= 5:
            rho, p = stats.spearmanr(sub.peptides, sub.cpu_seconds)
            print(f"  {STYLE[arm]['label']:42s} rho {rho:+.3f}  p {p:.4f}")


if __name__ == "__main__":
    main()
