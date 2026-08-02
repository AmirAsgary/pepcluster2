#!/usr/bin/env python3
"""Cross-analysis comparison of the two PepCluster2 scoring modes.

`separate_kmer_anchor` and `separate_aln_anchor` share everything except how the
primary similarity is computed: positionwise BLOSUM62 over the two terminal
3-mers, or a constrained alignment. Both are tuned and evaluated under one
protocol on one set of pools - select on the inner folds by AMI subject to
singletons being at most half of all clusters, then evaluate that selection on
the held-out alleles and on the test pools.

Figures come from `figures.py`; this module owns the tables and the report.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures as F  # noqa: E402
import plotstyle as PS  # noqa: E402

REPORTED = ["ami", "nmi", "adjusted_purity_macro", "adjusted_purity_micro",
            "bcubed_precision_macro", "singleton_fraction_of_clusters", "clusters"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    built = F.build_all(args.runs.resolve(), output)

    rows = []
    for mode, frames in built["selected"].items():
        for split, frame in frames.items():
            for record in frame.itertuples():
                rows.append({"mode": mode, "split": split, "pool": record.pool,
                             **{m: getattr(record, m) for m in REPORTED}})
    detail = pd.DataFrame(rows)
    detail.to_csv(output / "comparison_runs.csv", index=False)
    detail.groupby(["mode", "split"])[REPORTED].agg(["mean", "std"]).to_csv(
        output / "comparison_summary.csv")
    pd.DataFrame([
        {"mode": mode, **{k: choice[k] for k in
                          ("method", "representative_order", "primary_threshold",
                           "anchor_threshold", "ami", "singleton_fraction", "clusters")}}
        for mode, choice in built["choices"].items()
    ]).to_csv(output / "comparison_selected.csv", index=False)

    write_report(detail, built, output)
    print(detail.groupby(["mode", "split"])[["ami", "adjusted_purity_macro"]]
          .mean().round(4).to_string())


def write_report(detail: pd.DataFrame, built: dict, output: Path) -> None:
    lines: list[str] = []
    w = lines.append

    w("# PepCluster2: k-mer versus alignment scoring on the peptide-MHC benchmark")
    w("")
    w("Two scoring modes, one protocol, one set of pools. They differ only in the")
    w("primary similarity: `separate_kmer_anchor` scores the two terminal 3-mers")
    w("positionwise with BLOSUM62 and ignores the middle of the peptide;")
    w("`separate_aln_anchor` runs a constrained alignment over the whole peptide.")
    w("Both then apply the same anchor-combination condition as a second, separate")
    w("threshold.")
    w("")
    w("## Selected configurations")
    w("")
    w("| Mode | Similarity | Anchor | Order | Inner AMI | Singletons |")
    w("|---|---:|---:|---|---:|---:|")
    for mode, choice in built["choices"].items():
        w(f"| {PS.MODE_LABEL[mode]} | {choice['primary_threshold']:.2f} | "
          f"{choice['anchor_threshold']:.2f} | {choice['representative_order']} | "
          f"{choice['ami']:.4f} | {choice['singleton_fraction']:.3f} |")
    w("")
    for split, title in (("outer", "Held-out alleles"),
                         ("test", "Test pools")):
        w(f"## {title}")
        w("")
        w("| Mode | Pools | AMI | NMI | Purity (macro) | BCubed-P | Singletons | Clusters |")
        w("|---|---:|---:|---:|---:|---:|---:|---:|")
        subset = detail[detail.split == split]
        for mode in F.MODES:
            s = subset[subset["mode"] == mode]
            if s.empty:
                continue
            w(f"| {PS.MODE_LABEL[mode]} | {len(s)} | "
              f"{s.ami.mean():.4f} ± {s.ami.std():.4f} | {s.nmi.mean():.4f} | "
              f"{s.adjusted_purity_macro.mean():.4f} ± {s.adjusted_purity_macro.std():.4f} | "
              f"{s.bcubed_precision_macro.mean():.4f} | "
              f"{s.singleton_fraction_of_clusters.mean():.3f} | {s.clusters.mean():.0f} |")
        w("")

    # A selection sitting on the lowest threshold swept is not an optimum, it is
    # the edge of the grid; say so rather than quoting it as one.
    lowest = min(built["inner"][F.MODES[0]]["primary_threshold"].unique())
    at_edge = [mode for mode, choice in built["choices"].items()
               if choice["primary_threshold"] <= lowest]
    w("## Grid coverage")
    w("")
    w(f"Similarity thresholds swept: {lowest:g} to "
      f"{max(built['inner'][F.MODES[0]]['primary_threshold'].unique()):g}.")
    w("")
    if at_edge:
        for mode in at_edge:
            w(f"- **{PS.MODE_LABEL[mode]}** selects {lowest:g}, the lowest value swept.")
        w("  That is a boundary solution, not a located optimum: AMI may still be")
        w("  rising where the grid stops, so the setting should be read as \"at most")
        w(f"  {lowest:g}\" until lower values are run.")
    else:
        w("Every mode selects an interior value, so the grid brackets the optimum.")
    w("")
    w("## Figures")
    w("")
    w("- `hyperparameter_effect` - how far each of the two thresholds moves AMI,")
    w("  purity and cluster count. Error bars are the spread across the 5 outer folds.")
    w("- `overall_performance` - both modes at their selected setting on all three")
    w("  splits, so tuning-set optimism is visible as the gap between them.")
    w("- `benchmark_kmer`, `benchmark_aln` - the selected setting on the test pools,")
    w("  broken down by pool complexity and size.")
    w("")
    w("## Caveats")
    w("")
    w("- The test pools are drawn from a separate source file, but they are not")
    w("  disjoint in content: 66% of their peptides and 6 of 20 alleles also appear")
    w("  in the tuning data. They test transfer to new *pool compositions*, not to")
    w("  unseen peptides. The outer folds, which hold out whole alleles, are the")
    w("  stricter estimate.")
    w("- Alleles above the 97% similarity ceiling never share a pool, so these")
    w("  numbers are an upper bound relative to a realistic mixture.")
    w("- BCubed recall is near zero by construction - motif clustering splits each")
    w("  allele across many clusters - so only the precision half is informative.")
    w("")
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
