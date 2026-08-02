#!/usr/bin/env python3
"""Three figures, one question each, for the two PepCluster2 scoring modes.

1. hyperparameter_effect - how much do the two thresholds move the metrics?
2. overall_performance   - how do the modes compare across tuning, held-out
                           alleles and the benchmark?
3. benchmark_<mode>       - at the selected setting, how does one mode behave
                           across the benchmark pools?

Only these two modes appear anywhere in the output.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plotstyle as PS  # noqa: E402

MODES = ("separate_kmer_anchor", "separate_aln_anchor")
ROOTS = {"separate_kmer_anchor": "mhc_bench_sep_kmer_anchor",
         "separate_aln_anchor": "mhc_bench_sep_aln_anchor"}


def load(root: Path, tag: str) -> pd.DataFrame:
    frames = [pd.read_csv(f) for f in sorted((root / "grid").glob(f"{tag}_*.csv"))]
    frames = [f for f in frames if "status" in f.columns and len(f)]
    frame = pd.concat(frames, ignore_index=True).drop_duplicates(
        ["pool", "method", "representative_order", "primary_threshold", "anchor_threshold"])
    return frame[frame["status"] == "ok"].copy()


def figure_hyperparameter_effect(data: dict, output: Path) -> pd.DataFrame:
    """How the two thresholds move each metric. Graph path, coverage order."""
    import matplotlib.pyplot as plt

    rows = []
    figure, axes = plt.subplots(2, 3, figsize=(12.5, 6.4), constrained_layout=True)
    panels = [("ami", "AMI"), ("adjusted_purity_macro", "Per-allele purity"),
              ("clusters", "Clusters")]
    for column, (metric, title) in enumerate(panels):
        for row, (swept, xlabel) in enumerate(
                [("primary_threshold", "Similarity threshold"),
                 ("anchor_threshold", "Anchor threshold")]):
            axis = axes[row][column]
            for mode in MODES:
                frame = data[mode]
                frame = frame[(frame.method == "graph") &
                              (frame.representative_order == "coverage")]
                if swept == "anchor_threshold":
                    frame = frame[frame.anchor_threshold > 0]
                per_fold = frame.groupby([swept, "outer_fold"])[metric].mean().reset_index()
                stats = per_fold.groupby(swept)[metric].agg(["mean", "std"])
                axis.errorbar(stats.index, stats["mean"], yerr=stats["std"],
                              marker="o", capsize=2.5, elinewidth=0.9,
                              color=PS.MODE_COLOR[mode], label=PS.MODE_LABEL[mode])
                for threshold, record in stats.iterrows():
                    rows.append({"mode": mode, "swept": swept, "threshold": threshold,
                                 "metric": metric, "mean": record["mean"],
                                 "std": record["std"]})
            if metric == "clusters":
                axis.set_yscale("log")
            # Title on the top row, y-label on the bottom row: each panel is
            # named exactly once instead of carrying the metric twice.
            PS.finish(axis, xlabel, title if row == 1 else "",
                      title if row == 0 else "")
    axes[0][0].legend(loc="best")
    figure.suptitle("Effect of the two thresholds on clustering quality")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"hyperparameter_effect.{suffix}", dpi=200,
                       bbox_inches="tight")
    plt.close(figure)
    return pd.DataFrame(rows)


def figure_overall(selected: dict, output: Path) -> pd.DataFrame:
    """Both modes at their selected setting, across all three splits."""
    import matplotlib.pyplot as plt

    rows = []
    for mode, frames in selected.items():
        for split, frame in frames.items():
            for metric in ("ami", "adjusted_purity_macro", "nmi",
                           "singleton_fraction_of_clusters"):
                rows.append({"mode": mode, "split": split, "metric": metric,
                             "mean": frame[metric].mean(), "std": frame[metric].std(),
                             "pools": len(frame)})
    table = pd.DataFrame(rows)

    panels = [("ami", "AMI"), ("adjusted_purity_macro", "Per-allele purity"),
              ("singleton_fraction_of_clusters", "Singleton fraction")]
    splits = ["inner", "outer", "test"]
    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True)
    width = 0.36
    for axis, (metric, title) in zip(axes, panels):
        for offset, mode in zip((-width / 2, width / 2), MODES):
            values = [table[(table["mode"] == mode) & (table.split == s) &
                            (table.metric == metric)]["mean"].iloc[0] for s in splits]
            errors = [table[(table["mode"] == mode) & (table.split == s) &
                            (table.metric == metric)]["std"].iloc[0] for s in splits]
            axis.bar(np.arange(len(splits)) + offset, values, width, yerr=errors,
                     capsize=3, color=PS.MODE_COLOR[mode], label=PS.MODE_LABEL[mode],
                     edgecolor="white", linewidth=0.8,
                     error_kw=dict(elinewidth=0.9, ecolor=PS.INK))
        axis.set_xticks(np.arange(len(splits)))
        axis.set_xticklabels([PS.SPLIT_LABEL[s] for s in splits])
        axis.grid(axis="x", alpha=0)
        PS.finish(axis, "", "", title)
    axes[0].legend(loc="best")
    figure.suptitle("Performance at the selected setting")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"overall_performance.{suffix}", dpi=200,
                       bbox_inches="tight")
    plt.close(figure)
    return table


def figure_benchmark(mode: str, frame: pd.DataFrame, output: Path) -> None:
    """One mode on the benchmark pools: how quality varies with pool composition."""
    import matplotlib.pyplot as plt

    colour = PS.MODE_COLOR[mode]
    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True)

    grouped = frame.groupby("allele_count")[["ami", "adjusted_purity_macro"]].mean()
    axes[0].plot(grouped.index, grouped["ami"], marker="o", color=colour, label="AMI")
    axes[0].plot(grouped.index, grouped["adjusted_purity_macro"], marker="s",
                 linestyle="--", color=PS.SERIES[2], label="Per-allele purity")
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc="best")
    PS.finish(axes[0], "Alleles in pool", "", "Quality against pool complexity")

    axes[1].scatter(frame["peptides"], frame["ami"], s=26, color=colour,
                    alpha=0.8, edgecolor="white", linewidth=0.5)
    axes[1].set_xscale("log")
    axes[1].set_ylim(0, max(0.5, frame["ami"].max() * 1.15))
    PS.finish(axes[1], "Peptides in pool", "AMI", "Quality against pool size")

    axes[2].hist(frame["adjusted_purity_macro"], bins=14, color=colour, alpha=0.85,
                 edgecolor="white", linewidth=0.7)
    # Neutral, not red: red is the alignment mode's own colour on every other panel.
    axes[2].axvline(frame["adjusted_purity_macro"].mean(), color=PS.INK,
                    linestyle="--", linewidth=1.4)
    PS.finish(axes[2], "Per-allele purity", "Pools", "Spread across pools")

    figure.suptitle(f"Benchmark evaluation: {PS.MODE_LABEL[mode]} "
                    f"({len(frame)} pools)")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"benchmark_{stem_of(mode)}.{suffix}", dpi=200,
                       bbox_inches="tight")
    plt.close(figure)


def build_all(runs: Path, output: Path) -> dict:
    """Write the three figures and their backing CSVs; return what the report needs."""
    output.mkdir(parents=True, exist_ok=True)
    PS.apply()

    inner = {mode: load(runs / ROOTS[mode], "inner") for mode in MODES}
    curve = figure_hyperparameter_effect(inner, output)
    curve.to_csv(output / "hyperparameter_effect.csv", index=False)

    choices, selected = {}, {}
    for mode in MODES:
        root = runs / ROOTS[mode]
        choice = pd.read_csv(root / "tables" / "selected_overall.csv")
        choice = choice[choice.method == "graph"].iloc[0]
        choices[mode] = choice
        selected[mode] = {
            split: (lambda f: f[(f.method == "graph") &
                                (f.representative_order == choice.representative_order) &
                                (f.primary_threshold == choice.primary_threshold) &
                                (f.anchor_threshold == choice.anchor_threshold)])(
                load(root, split))
            for split in ("inner", "outer", "test")}

    table = figure_overall(selected, output)
    table.to_csv(output / "overall_performance.csv", index=False)
    for mode in MODES:
        figure_benchmark(mode, selected[mode]["test"], output)
        selected[mode]["test"].to_csv(
            output / f"benchmark_{stem_of(mode)}_pools.csv", index=False)
    return {"curve": curve, "summary": table, "choices": choices, "selected": selected,
            "inner": inner}


def stem_of(mode: str) -> str:
    return "kmer" if "kmer" in mode else "aln"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    built = build_all(args.runs.resolve(), args.output.resolve())
    print(built["summary"].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
