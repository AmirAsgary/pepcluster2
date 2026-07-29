#!/usr/bin/env python3
"""Compare best-assignment, all-four-dimer, and aligned-3-mer formulas."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path


def collect(root: Path) -> dict[str, float]:
    run_dirs = sorted((root / "runs" / "no_prefilter").glob("sample_*"))
    cluster_counts, edge_counts, singleton_peptide_fractions = [], [], []
    for directory in run_dirs:
        stats = json.loads((directory / "run_stats.json").read_text())
        cluster_counts.append(stats["final_clusters"])
        edge_counts.append(stats["graph_edge_count"])
        with (directory / "cluster_summary.tsv").open(newline="") as handle:
            sizes = [int(row["size"]) for row in csv.DictReader(handle, delimiter="\t")]
        singleton_peptide_fractions.append(sum(size == 1 for size in sizes) / sum(sizes))
    return {
        "datasets": len(run_dirs),
        "mean_eligible_edges": statistics.mean(edge_counts),
        "mean_clusters": statistics.mean(cluster_counts),
        "mean_singleton_peptide_fraction": statistics.mean(singleton_peptide_fractions),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--best-assignment", type=Path, required=True)
    parser.add_argument("--all-four", type=Path, required=True)
    parser.add_argument("--aligned-three-mer", type=Path, required=True)
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--figure-output", type=Path, required=True)
    args = parser.parse_args()
    result = {
        "best_assignment_kmer_formula": collect(args.best_assignment),
        "all_four_dimer_kmer_formula": collect(args.all_four),
        "aligned_three_mer_kmer_formula": collect(args.aligned_three_mer),
    }
    args.json_output.write_text(json.dumps(result, indent=2) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = [
        ("Best assignment", result["best_assignment_kmer_formula"]),
        ("All four dimers", result["all_four_dimer_kmer_formula"]),
        ("Aligned 3-mers", result["aligned_three_mer_kmer_formula"]),
    ]
    panels = [
        ("mean_eligible_edges", "Mean eligible edges", "count"),
        ("mean_clusters", "Mean clusters", "count"),
        ("mean_singleton_peptide_fraction", "Peptides in singleton clusters", "fraction"),
    ]
    figure, axes = plt.subplots(1, 3, figsize=(11, 4))
    for axis, (key, title, kind) in zip(axes, panels):
        values = [group[key] for _, group in groups]
        bars = axis.bar([label for label, _ in groups], values, color=["#718096", "#c05621", "#2b6cb0"])
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=15)
        if kind == "fraction":
            axis.set_ylim(0, 1)
            axis.set_ylabel("Fraction of input peptides")
            labels = [f"{100*x:.1f}%" for x in values]
        else:
            labels = [f"{x:,.0f}" for x in values]
        axis.bar_label(bars, labels=labels, padding=3)
    figure.suptitle("Comparison of PepCluster2 k-mer similarity formulas")
    figure.tight_layout()
    figure.savefig(args.figure_output, dpi=200, bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
