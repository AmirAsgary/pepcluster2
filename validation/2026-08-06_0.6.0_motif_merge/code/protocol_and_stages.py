#!/usr/bin/env python3
"""Two things the benchmark figures do not cover.

fig1_protocol   What the nested selection actually did: the selection surface on
                the tuning folds, the configuration it chose, and the score on
                each split. The point is that the test pools contributed nothing
                to the choice.

stage contribution
                What the merge stage and the EM stage each contribute, measured
                on EVERY metric rather than on AMI alone. An earlier reading
                took only AMI and concluded the merge was worth ~0.02; precision
                and recall move in opposite directions when granularity changes,
                so a single chance-corrected summary can hide a real effect.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = [("ami", "AMI"), ("adjusted_purity_macro", "Adjusted purity"),
           ("bcubed_precision_macro", "BCubed precision"),
           ("bcubed_recall_macro", "BCubed recall"),
           ("bcubed_f1_macro", "BCubed F1")]
CONFIG = ["merge_concentration", "merge_threshold", "em", "em_concentration"]


def load_grid(grid: Path) -> pd.DataFrame:
    frame = pd.concat([pd.read_csv(f) for f in sorted(grid.glob("motif_shard_*.csv"))],
                      ignore_index=True)
    frame = frame[frame.status == "ok"].copy()
    frame["em_concentration"] = frame.em_concentration.fillna(-1.0)
    return frame.drop_duplicates(["pool"] + CONFIG)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})

    frame = load_grid(args.grid)
    choice = pd.read_csv(args.selected).iloc[0]
    inner = frame[frame.split == "inner"]

    # ---------- stage contribution, every metric ------------------------
    # Three arms, each at its own best configuration chosen on the INNER folds,
    # then read off the test pools. Choosing each arm's own optimum is what makes
    # the comparison fair: penalising an arm for a configuration tuned to a
    # different arm would understate it.
    rows = []
    arms = {
        "similarity only": None,
        "merge only (no EM)": inner[~inner.em.astype(bool)],
        "merge + EM": inner[inner.em.astype(bool)],
        "EM, minimal merging": inner[inner.em.astype(bool) &
                                     (inner.merge_concentration ==
                                      inner.merge_concentration.min()) &
                                     (inner.merge_threshold ==
                                      inner.merge_threshold.max())],
    }
    test = frame[frame.split == "test"]
    for name, pool in arms.items():
        if name == "similarity only":
            base = test.drop_duplicates("pool")
            rows.append(dict(arm=name, ami=base.similarity_ami.mean(),
                             bcubed_f1_macro=base.similarity_f1.mean(),
                             clusters=base.similarity_clusters.mean()))
            continue
        best = pool.groupby(CONFIG, dropna=False).ami.mean().idxmax()
        mask = pd.Series(True, index=test.index)
        for column, value in zip(CONFIG, best):
            mask &= test[column] == value
        chosen = test[mask]
        rows.append(dict(arm=name, **{m: chosen[m].mean() for m, _ in METRICS},
                         clusters=chosen.clusters.mean(),
                         config=", ".join(f"{c}={v}" for c, v in zip(CONFIG, best))))
    stages = pd.DataFrame(rows)
    stages.to_csv(args.out / "fig1c_stage_contribution_all_metrics.csv", index=False)
    print("=== stage contribution on the test pools, every metric ===")
    print(stages.round(4).to_string(index=False))

    # deltas between consecutive stages, per metric
    def get(arm, metric):
        row = stages[stages.arm == arm]
        return float(row[metric].iloc[0]) if metric in row and len(row) else np.nan

    deltas = []
    for metric, label in METRICS:
        deltas.append(dict(
            metric=label,
            similarity=get("similarity only", metric),
            merge_only=get("merge only (no EM)", metric),
            merge_em=get("merge + EM", metric),
            em_minimal_merge=get("EM, minimal merging", metric),
            gain_from_merge=get("merge only (no EM)", metric) - get("similarity only", metric),
            gain_from_em=get("merge + EM", metric) - get("merge only (no EM)", metric),
            merge_worth_given_em=get("merge + EM", metric) - get("EM, minimal merging", metric)))
    deltas = pd.DataFrame(deltas)
    deltas.to_csv(args.out / "fig1c_stage_deltas_all_metrics.csv", index=False)
    print("\n=== per-metric deltas ===")
    print(deltas.round(4).to_string(index=False))

    # ---------- fig 1: protocol -----------------------------------------
    figure = plt.figure(figsize=(14, 4.4), constrained_layout=True)
    spec = figure.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0])

    # (a) selection surface on the tuning folds
    axis = figure.add_subplot(spec[0, 0])
    surface = inner[inner.em.astype(bool)].groupby(
        ["em_concentration", "merge_concentration"]).ami.mean().reset_index()
    for conc, group in surface.groupby("merge_concentration"):
        group = group.sort_values("em_concentration")
        axis.plot(group.em_concentration, group.ami, marker="o", markersize=3.5,
                  lw=1.3, label=f"merge conc. {conc:g}")
    axis.axvline(choice.em_concentration, color="black", ls=":", lw=1.2)
    axis.set_xscale("log")
    axis.set_xlabel("EM prior concentration")
    axis.set_ylabel("Mean AMI on tuning folds")
    axis.set_title("(a) Selection surface, tuning folds only")
    axis.legend(fontsize=6.5, ncol=2)
    axis.grid(alpha=0.25, lw=0.5)
    surface.to_csv(args.out / "fig1a_selection_surface_tuning_folds.csv", index=False)

    # (b) per-fold choice reproducibility
    axis = figure.add_subplot(spec[0, 1])
    per_fold = pd.read_csv(args.selected.parent / "motif_selected_per_fold.csv")
    axis.bar(per_fold.held_out_fold.astype(str), per_fold.ami,
             color="black", alpha=0.8, edgecolor="white")
    axis.axhline(choice.ami, color="#d62728", ls="--", lw=1.4,
                 label=f"independent test = {choice.ami:.3f}")
    axis.set_xlabel("Held-out allele fold")
    axis.set_ylabel("AMI on that fold")
    axis.set_ylim(0, 1)
    axis.set_title("(b) Each fold scored by a choice\nmade without it")
    axis.legend(fontsize=7)
    axis.grid(axis="y", alpha=0.25, lw=0.5)
    per_fold.to_csv(args.out / "fig1b_per_fold_evaluation.csv", index=False)

    # (c) the pipeline, stage by stage
    axis = figure.add_subplot(spec[0, 2])
    order = ["similarity only", "merge only (no EM)", "merge + EM"]
    x = np.arange(len(order))
    for offset, (metric, label) in zip((-0.2, 0.2), (METRICS[0], METRICS[4])):
        values = [get(a, metric) for a in order]
        axis.bar(x + offset, values, 0.38, label=label,
                 color="black" if offset < 0 else "#777777",
                 edgecolor="white", linewidth=0.5)
    axis.set_xticks(x)
    axis.set_xticklabels(["similarity", "+ merge", "+ EM"], fontsize=8)
    axis.set_ylabel("Score on independent test")
    axis.set_ylim(0, 1)
    axis.set_title("(c) What each stage adds")
    axis.legend(fontsize=7)
    axis.grid(axis="y", alpha=0.25, lw=0.5)

    figure.suptitle("Nested selection protocol: hyperparameters chosen on tuning "
                    "folds, evaluated once on held-out alleles and independent test")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig1_protocol.{suffix}", bbox_inches="tight")
    plt.close(figure)
    print("\nwrote fig1_protocol and stage tables to", args.out)


if __name__ == "__main__":
    main()
