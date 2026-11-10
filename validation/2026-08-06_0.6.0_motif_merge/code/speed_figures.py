#!/usr/bin/env python3
"""Cost against pool size, and cost against accuracy.

Two panels each. CPU seconds is the primary axis because it is the resource
actually consumed and is comparable across tools regardless of how each one
happens to be threaded; wall seconds is shown beside it because it is what a user
waits for. PepCluster2 and MixMHCp are single-threaded here, so for them the two
are the same number. GibbsCluster is run at -k 16, so its wall time is 16x
cheaper-looking than its cost.

Every arm is timed several times per pool and the repeats are reduced to their
median before anything is plotted. At one repeat the cheap PepCluster2 arms are
separated by less than the run-to-run spread, which is how an earlier version of
this figure came to show a given-K arm as faster than clustering alone.

Colour and line style follow `figures.py` exactly, so an arm looks the same in
every figure of this study.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figures import STYLE, ORDER, MODES  # noqa: E402

# speed_benchmark.py writes the variant stem for our arms and a plain name for
# the comparators. Map both onto the display names the rest of the study uses.
ARM_LABEL = {stem: label for stem, label in MODES}
ARM_LABEL.update({"mixmhcp": "MixMHCp", "gibbscluster": "GibbsCluster"})


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
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15,
                         "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12, "figure.titlesize": 16,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 150})

    raw = pd.read_csv(args.speed)
    raw = raw[raw.get("status", "ok").astype(str).eq("ok") if "status" in raw
              else np.ones(len(raw), bool)]
    unknown = sorted(set(raw.arm) - set(ARM_LABEL))
    if unknown:
        print(f"warning: arms with no display name, dropped: {unknown}")
    raw = raw[raw.arm.isin(ARM_LABEL)].copy()
    raw["tool"] = raw.arm.map(ARM_LABEL)
    repeats = raw.groupby(["tool", "pool"]).size()
    print(f"{len(raw)} rows, {raw.tool.nunique()} arms, "
          f"{raw.pool.nunique()} pools, {repeats.min()}-{repeats.max()} "
          f"repeats per arm and pool")

    # One row per (arm, pool): the median over repeats, which is what a reader
    # should compare. The spread across repeats is carried so the figure can say
    # when two arms are not actually separated.
    metrics = [c for c in ("bcubed_f1_macro", "ami") if c in raw.columns]
    frame = raw.groupby(["tool", "pool"], as_index=False).agg(
        peptides=("peptides", "first"), alleles=("alleles", "first"),
        cpu_seconds=("cpu_seconds", "median"),
        cpu_seconds_min=("cpu_seconds", "min"),
        cpu_seconds_max=("cpu_seconds", "max"),
        wall_seconds=("wall_seconds", "median"),
        repeats=("cpu_seconds", "size"),
        **{m: (m, "mean") for m in metrics})
    order = [t for t in ORDER if t in set(frame.tool)]

    def decade(caxis):
        caxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0),
                                                  numticks=30))
        caxis.set_minor_locator(ticker.LogLocator(base=10.0,
                                                  subs=tuple(np.arange(1, 10) * 0.1),
                                                  numticks=100))
        caxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
        caxis.set_minor_formatter(ticker.NullFormatter())

    def outside_legend(figure, axis):
        handles, labels = axis.get_legend_handles_labels()
        index = {label: i for i, label in enumerate(labels)}
        picked = [(handles[index[t]], t) for t in order if t in index]
        figure.legend([h for h, _ in picked], [t for _, t in picked],
                      loc="outside right center", frameon=False, fontsize=12)

    # ---- fig6: cost against pool size ---------------------------------
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.4), constrained_layout=True)
    for axis, column, title in ((axes[0], "cpu_seconds", "CPU seconds"),
                                (axes[1], "wall_seconds", "Wall seconds")):
        for tool in order:
            sub = frame[frame.tool == tool].sort_values("peptides")
            if sub.empty:
                continue
            style = dict(STYLE[tool])
            marker = style.pop("marker")
            if column == "cpu_seconds":
                axis.fill_between(sub.peptides, sub.cpu_seconds_min,
                                  sub.cpu_seconds_max, color=style["color"],
                                  alpha=0.12, lw=0, zorder=style["zorder"] - 0.5)
            axis.plot(sub.peptides, sub[column], marker=marker, markersize=6,
                      label=tool, **style)
        axis.set_xscale("log")
        axis.set_yscale("log")
        decade(axis.xaxis)
        decade(axis.yaxis)
        axis.set_xlabel("Peptides in pool")
        axis.set_ylabel(title)
        axis.set_title(title)
        axis.grid(alpha=0.25, lw=0.5, which="both")
    outside_legend(figure, axes[0])
    figure.suptitle("Cost against pool size (serial runs, median of repeats; "
                    "band = min to max across repeats; PepCluster2 and MixMHCp "
                    "single-threaded, GibbsCluster at -k 16)")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig6_speed_vs_pool_size.{suffix}",
                       bbox_inches="tight")
    plt.close(figure)
    frame.to_csv(args.out / "fig6_speed_vs_pool_size.csv", index=False)

    # ---- fig7: cost against accuracy ----------------------------------
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.4), constrained_layout=True)
    summary = []
    titles = {"bcubed_f1_macro": "BCubed F1", "ami": "AMI"}
    for axis, metric in zip(axes, metrics):
        for tool in order:
            sub = frame[frame.tool == tool]
            if sub.empty:
                continue
            x = sub.cpu_seconds.median()
            y = sub[metric].mean()
            err = sub[metric].std(ddof=1) / max(np.sqrt(len(sub)), 1)
            style = dict(STYLE[tool])
            # A point, not a line: drop the line properties and keep the mark.
            # A scatter cannot show a dash, so the given-K arms are drawn hollow
            # instead; colour still means the pipeline.
            given_k = style.pop("ls") != "-"
            style.pop("lw")
            if given_k:
                style["markerfacecolor"] = "white"
                style["markeredgecolor"] = style["color"]
                style["markeredgewidth"] = 2.0
            axis.errorbar([x], [y], yerr=[err], markersize=13, capsize=3,
                          lw=0, elinewidth=1.4, label=tool, **style)
            if metric == metrics[0]:
                summary.append(dict(tool=tool, median_cpu_seconds=x,
                                    **{f"mean_{m}": sub[m].mean() for m in metrics},
                                    pools=len(sub),
                                    repeats=int(sub.repeats.median())))
        axis.set_xscale("log")
        decade(axis.xaxis)
        axis.set_xlabel("Median CPU seconds per pool")
        axis.set_ylabel(titles.get(metric, metric))
        axis.set_ylim(0, 1)
        axis.set_title(f"{titles.get(metric, metric)} against cost")
        axis.grid(alpha=0.25, lw=0.5, which="both")
    outside_legend(figure, axes[0])
    figure.suptitle("Accuracy against cost: up and to the left is better; hollow = given the true motif count")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig7_speed_vs_performance.{suffix}",
                       bbox_inches="tight")
    plt.close(figure)
    table = pd.DataFrame(summary)
    table.to_csv(args.out / "fig7_speed_vs_performance.csv", index=False)

    print(table.round(4).to_string(index=False))
    base_rows = table[table.tool == "PepCluster2 + merge + EM"]
    if not base_rows.empty:
        base = base_rows.median_cpu_seconds.iloc[0]
        print("\ncost relative to PepCluster2 + merge + EM:")
        for row in table.itertuples():
            print(f"  {row.tool:42s} {row.median_cpu_seconds / base:8.2f}x")

    # Does cost scale differently with pool size across tools?
    print("\nSpearman rho of CPU seconds against pool size:")
    for tool in order:
        sub = frame[frame.tool == tool]
        if len(sub) >= 5:
            rho, p = stats.spearmanr(sub.peptides, sub.cpu_seconds)
            print(f"  {tool:42s} rho {rho:+.3f}  p {p:.4f}")

    # Are the cheap arms actually separated? Paired over pools, on the medians.
    print("\npaired Wilcoxon on CPU seconds, each arm against PepCluster2:")
    base = frame[frame.tool == "PepCluster2"].set_index("pool").cpu_seconds
    for tool in order:
        if tool == "PepCluster2":
            continue
        sub = frame[frame.tool == tool].set_index("pool").cpu_seconds
        shared = base.index.intersection(sub.index)
        if len(shared) >= 6:
            stat, p = stats.wilcoxon(sub[shared], base[shared])
            ratio = (sub[shared] / base[shared]).median()
            print(f"  {tool:42s} {ratio:6.2f}x  p {p:.4f}  n {len(shared)}")


if __name__ == "__main__":
    main()
