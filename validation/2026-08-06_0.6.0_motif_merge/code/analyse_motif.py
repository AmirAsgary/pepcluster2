#!/usr/bin/env python3
"""Nested selection over the motif grid.

Protocol, mirroring `analyse.py` of the 0.5.0 study:

  per fold   For each outer fold f, the configuration with the best mean score
             over the inner pools of the OTHER folds is evaluated on fold f's
             held-out-allele pools. Nothing from fold f informs its own choice.
  overall    One configuration selected over all 120 inner pools is evaluated
             once on the 48 independent test pools.

Selection metric. `metrics.objective` averages AMI, NMI and adjusted per-allele
purity. That composite is not appropriate here: adjusted purity is BCubed
precision corrected only against the allele prior, so it still rises with
fragmentation, and the motif layer's whole effect is to change granularity by an
order of magnitude. Selecting on it would systematically prefer configurations
that merge too little. AMI is used instead - chance-corrected over the whole
partition, and the metric the external-tool comparison already leads with.
BCubed F1 is reported alongside, and the selection is repeated under F1 as a
robustness check: if the two disagree, that is reported rather than hidden.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CONFIG = ["merge_concentration", "merge_threshold", "em", "em_concentration"]
REPORT = ["ami", "bcubed_f1_macro", "adjusted_purity_macro", "bcubed_recall_macro",
          "clusters", "merged_groups", "singleton_fraction_of_clusters"]


def load(grid: Path) -> pd.DataFrame:
    files = sorted(grid.glob("motif_shard_*.csv"))
    if not files:
        raise SystemExit(f"no shard files under {grid}")
    frame = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    bad = (frame.status != "ok").sum()
    if bad:
        print(f"warning: {bad} failed runs excluded")
        print(frame[frame.status != "ok"].error.value_counts().head().to_string())
    frame = frame[frame.status == "ok"].copy()
    # em=False rows carry no EM concentration. groupby drops NaN keys, which
    # would silently delete the merge-only arm from every selection.
    frame["em_concentration"] = frame["em_concentration"].fillna(-1.0)
    frame = frame.drop_duplicates(["pool"] + CONFIG)
    return frame


def choose(inner: pd.DataFrame, fold, metric: str) -> pd.Series:
    data = inner if fold is None else inner[inner.outer_fold != fold]
    grouped = data.groupby(CONFIG, dropna=False).agg(
        score=(metric, "mean"), pools=(metric, "size")).reset_index()
    # Deterministic tie-break, so a rerun cannot silently pick a different
    # configuration when two are numerically equal.
    grouped = grouped.sort_values(
        ["score"] + CONFIG, ascending=[False] + [True] * len(CONFIG))
    return grouped.iloc[0]


def evaluate(frame: pd.DataFrame, choice: pd.Series) -> pd.DataFrame:
    mask = pd.Series(True, index=frame.index)
    for column in CONFIG:
        mask &= frame[column] == choice[column]
    return frame[mask]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--metric", default="ami")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    frame = load(args.grid)
    inner = frame[frame.split == "inner"]
    outer = frame[frame.split == "outer"]
    test = frame[frame.split == "test"]
    print(f"{len(frame)} runs | inner {inner.pool.nunique()} pools, "
          f"outer {outer.pool.nunique()}, test {test.pool.nunique()} | "
          f"{inner.groupby(CONFIG, dropna=False).ngroups} configurations")

    # ---- per-fold nested evaluation on held-out alleles -------------------
    rows = []
    for fold in sorted(inner.outer_fold.unique()):
        choice = choose(inner, fold, args.metric)
        held = evaluate(outer[outer.outer_fold == fold], choice)
        if held.empty:
            continue
        rows.append({"held_out_fold": fold,
                     **{c: choice[c] for c in CONFIG},
                     "inner_score": choice["score"], "pools": len(held),
                     **{m: held[m].mean() for m in REPORT}})
    per_fold = pd.DataFrame(rows)
    per_fold.to_csv(args.out / "motif_selected_per_fold.csv", index=False)

    # ---- one overall selection, evaluated once on test --------------------
    overall = choose(inner, None, args.metric)
    on_test = evaluate(test, overall)
    summary = pd.DataFrame([{
        **{c: overall[c] for c in CONFIG},
        "inner_score": overall["score"], "test_pools": len(on_test),
        **{m: on_test[m].mean() for m in REPORT},
        **{f"{m}_std": on_test[m].std() for m in ("ami", "bcubed_f1_macro")},
    }])
    summary.to_csv(args.out / "motif_selected_overall.csv", index=False)

    # robustness: does selecting on F1 pick the same configuration?
    alternative = choose(inner, None, "bcubed_f1_macro")
    agrees = all(alternative[c] == overall[c] for c in CONFIG)

    # baseline: the similarity partition these motifs were built from, and the
    # merge-only arm, both on the same test pools.
    baseline = test.drop_duplicates("pool")
    merge_only = test[~test.em.astype(bool)]
    merge_choice = choose(inner[~inner.em.astype(bool)], None, args.metric)
    merge_only = evaluate(merge_only, merge_choice)

    print("\n=== selected on inner folds, by", args.metric, "===")
    print(summary.round(4).to_string(index=False))
    print("\n=== per held-out fold ===")
    if not per_fold.empty:
        print(per_fold.round(4).to_string(index=False))
    print("\n=== reference, same 48 test pools ===")
    print(f"similarity clustering  AMI {baseline.similarity_ami.mean():.4f}  "
          f"F1 {baseline.similarity_f1.mean():.4f}  "
          f"clusters {baseline.similarity_clusters.mean():.1f}")
    if not merge_only.empty:
        print(f"merge only, no EM      AMI {merge_only.ami.mean():.4f}  "
              f"F1 {merge_only.bcubed_f1_macro.mean():.4f}  "
              f"clusters {merge_only.clusters.mean():.1f}  "
              f"(merge concentration {merge_choice['merge_concentration']}, "
              f"threshold {merge_choice['merge_threshold']})")
    print(f"merge + EM             AMI {on_test.ami.mean():.4f}  "
          f"F1 {on_test.bcubed_f1_macro.mean():.4f}  "
          f"clusters {on_test.clusters.mean():.1f}")
    print("\nMixMHCp default AMI 0.3918 F1 0.3918 | forced k AMI 0.4915 F1 0.4728")
    print(f"\nselection under F1 agrees with selection under {args.metric}: {agrees}")
    if not agrees:
        print("  F1 would choose:",
              {c: alternative[c] for c in CONFIG})

    # per-band breakdown of the selected configuration
    banded = on_test.copy()
    banded["band"] = pd.cut(banded.allele_count, bins=[1, 6, 12, 30],
                            labels=["2-6", "7-12", "13-30"])
    table = banded.groupby("band", observed=True).agg(
        pools=("pool", "size"), alleles=("allele_count", "mean"),
        ami=("ami", "mean"), f1=("bcubed_f1_macro", "mean"),
        purity=("adjusted_purity_macro", "mean"),
        recall=("bcubed_recall_macro", "mean"), clusters=("clusters", "mean"))
    table.to_csv(args.out / "motif_test_by_band.csv")
    print("\n=== selected configuration, test pools by complexity ===")
    print(table.round(4).to_string())


if __name__ == "__main__":
    main()
