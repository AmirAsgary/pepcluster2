#!/usr/bin/env python3
"""Plot the exact PepCluster2 400 x 400 dimer-similarity distribution."""

from __future__ import annotations

import argparse
import csv
import json
import struct
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MAGIC = b"PC2K2S01"
N_DIMERS = 400
N_SCORES = N_DIMERS * N_DIMERS


def read_scores(path: Path) -> np.ndarray:
    payload = path.read_bytes()
    expected_size = 16 + N_SCORES * 2
    if len(payload) != expected_size or payload[:8] != MAGIC:
        raise ValueError(f"not a PepCluster2 dimer table: {path}")
    if struct.unpack_from("<I", payload, 8)[0] != N_DIMERS:
        raise ValueError(f"unexpected dimer count in {path}")
    return np.frombuffer(payload, dtype="<i2", offset=16).astype(np.float64) / 1000.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    scores = read_scores(args.table)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    edges = np.arange(-0.85, 1.0001, 0.05)
    counts, edges = np.histogram(scores, bins=edges)
    fractions = counts / scores.size
    with (args.output_dir / "dimer_similarity_histogram.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin_lower", "bin_upper", "count", "fraction_of_all_pairs"])
        for lower, upper, count, fraction in zip(edges[:-1], edges[1:], counts, fractions):
            writer.writerow([f"{lower:.2f}", f"{upper:.2f}", int(count), f"{fraction:.8f}"])

    unique, unique_counts = np.unique(scores, return_counts=True)
    with (args.output_dir / "dimer_similarity_exact_frequencies.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["similarity", "count", "fraction_of_all_pairs"])
        for value, count in zip(unique, unique_counts):
            writer.writerow([f"{value:.3f}", int(count), f"{count / scores.size:.8f}"])

    summary = {
        "ordered_dimer_pairs": int(scores.size),
        "unique_dimer_similarity_values": int(unique.size),
        "minimum": float(scores.min()),
        "maximum": float(scores.max()),
        "mean": float(scores.mean()),
        "median": float(np.median(scores)),
        "negative_count": int(np.count_nonzero(scores < 0)),
        "zero_count": int(np.count_nonzero(scores == 0)),
        "positive_count": int(np.count_nonzero(scores > 0)),
        "at_least_0_50_count": int(np.count_nonzero(scores >= 0.50)),
        "at_least_0_50_fraction": float(np.mean(scores >= 0.50)),
        "at_least_0_60_count": int(np.count_nonzero(scores >= 0.60)),
        "at_least_0_60_fraction": float(np.mean(scores >= 0.60)),
        "at_least_0_75_count": int(np.count_nonzero(scores >= 0.75)),
        "at_least_0_75_fraction": float(np.mean(scores >= 0.75)),
        "identical_score_count": int(np.count_nonzero(scores == 1.0)),
    }
    (args.output_dir / "dimer_similarity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    weights = np.full(scores.size, 100.0 / scores.size)
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)

    axes[0].hist(scores, bins=edges, weights=weights, color="#31688e", edgecolor="white", linewidth=0.4)
    axes[0].axvline(0.50, color="#d1495b", linewidth=2, linestyle="--", label="Seed-neighbour cutoff = 0.50")
    axes[0].set_title("A. Complete dimer-similarity distribution")
    axes[0].set_xlabel("Normalized BLOSUM62 dimer similarity")
    axes[0].set_ylabel("Ordered dimer pairs per bin (%)")
    axes[0].legend(frameon=False)

    nonnegative_edges = np.arange(0.0, 1.0001, 0.025)
    axes[1].hist(scores[scores >= 0], bins=nonnegative_edges, weights=np.full(np.count_nonzero(scores >= 0), 100.0 / scores.size), color="#35b779", edgecolor="white", linewidth=0.4)
    axes[1].axvline(0.50, color="#d1495b", linewidth=2, linestyle="--")
    axes[1].set_title("B. Nonnegative tail (percentage of all pairs)")
    axes[1].set_xlabel("Normalized BLOSUM62 dimer similarity")
    axes[1].set_ylabel("Ordered dimer pairs per bin (%)")
    axes[1].text(
        0.98,
        0.95,
        f"≥ 0.50: {summary['at_least_0_50_fraction']:.2%}\n"
        f"≥ 0.60: {summary['at_least_0_60_fraction']:.2%}\n"
        f"≥ 0.75: {summary['at_least_0_75_fraction']:.2%}",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )

    figure.suptitle(
        "PepCluster2 dimer similarities: all 160,000 ordered 2-mer comparisons",
        fontsize=14,
    )
    figure.savefig(args.output_dir / "dimer_similarity_histogram.png", dpi=220)
    figure.savefig(args.output_dir / "dimer_similarity_histogram.pdf")


if __name__ == "__main__":
    main()
