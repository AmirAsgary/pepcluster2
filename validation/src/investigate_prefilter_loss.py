#!/usr/bin/env python3
"""Classify why a scoped prefilter run differs from its non-prefilter baseline."""

from __future__ import annotations

import argparse
import csv
import json
import struct
from collections import Counter, defaultdict
from pathlib import Path

from run_prefilter_comparison import adjusted_rand_index, canonical_partition, read_partition


def read_edges(path: Path) -> dict[tuple[int, int], int]:
    data = path.read_bytes()
    if len(data) % 10:
        raise ValueError(f"truncated edge file: {path}")
    return {(u, v): weight for u, v, weight in struct.iter_unpack("<IIH", data)}


def labels_in_node_order(path: Path) -> tuple[list[str], list[tuple[str, str]], list[tuple[str, str]]]:
    labels, nodes, reps = [], [], []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels.append(row["cluster_id"])
            nodes.append((row["anchor"], row["geometry_mask"]))
            reps.append((row["representative_anchor"], row["representative_geometry_mask"]))
    return labels, nodes, reps


def partition_comparison(first: Path, second: Path) -> dict:
    a, _ = read_partition(first)
    b, _ = read_partition(second)
    keys = sorted(a)
    return {
        "adjusted_rand_index": adjusted_rand_index([a[k] for k in keys], [b[k] for k in keys]),
        "exact_partition": canonical_partition(a) == canonical_partition(b),
        "clusters_first": len(set(a.values())),
        "clusters_second": len(set(b.values())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-output", type=Path, required=True)
    parser.add_argument("--scoped-output", type=Path, required=True)
    parser.add_argument("--full-output", type=Path, required=True)
    parser.add_argument("--baseline-edges", type=Path, required=True)
    parser.add_argument("--scoped-edges", type=Path, required=True)
    parser.add_argument("--full-edges", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    baseline_edges = read_edges(args.baseline_edges)
    scoped_edges = read_edges(args.scoped_edges)
    full_edges = read_edges(args.full_edges)
    baseline_labels, nodes, baseline_reps = labels_in_node_order(args.baseline_output / "node_clusters.tsv")
    scoped_labels, scoped_nodes, _ = labels_in_node_order(args.scoped_output / "node_clusters.tsv")
    if nodes != scoped_nodes:
        raise ValueError("node orders differ")

    provisional = []
    with (args.scoped_output / "prefilter_provisional_clusters.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            provisional.append({
                "cluster": int(row["provisional_cluster"]),
                "representative": int(row["provisional_representative_node"]),
                "assigned": row["assigned_to_nonsingleton"] == "true",
            })
    is_rep = [i == row["representative"] and row["assigned"] for i, row in enumerate(provisional)]
    is_unassigned = [not row["assigned"] for row in provisional]

    missing = set(baseline_edges) - set(scoped_edges)
    categories = Counter()
    missing_within_baseline_cluster = 0
    missing_baseline_rep_member = 0
    for u, v in missing:
        if is_rep[u] or is_rep[v]:
            category = "provisional_representative_to_any"
        elif is_unassigned[u] and is_unassigned[v]:
            category = "unassigned_to_unassigned"
        elif is_unassigned[u] != is_unassigned[v]:
            category = "assigned_nonrepresentative_to_unassigned"
        elif provisional[u]["cluster"] == provisional[v]["cluster"]:
            category = "assigned_nonrepresentatives_same_provisional_cluster"
        else:
            category = "assigned_nonrepresentatives_different_provisional_clusters"
        categories[category] += 1
        if baseline_labels[u] == baseline_labels[v]:
            missing_within_baseline_cluster += 1
            if baseline_reps[u] == nodes[u] or baseline_reps[v] == nodes[v]:
                missing_baseline_rep_member += 1

    spread = defaultdict(set)
    for baseline, scoped in zip(baseline_labels, scoped_labels):
        spread[baseline].add(scoped)
    split_counts = [len(groups) for groups in spread.values()]

    result = {
        "baseline_vs_scoped": partition_comparison(args.baseline_output / "node_clusters.tsv", args.scoped_output / "node_clusters.tsv"),
        "baseline_vs_full_sensitive_after_prefilter": partition_comparison(args.baseline_output / "node_clusters.tsv", args.full_output / "node_clusters.tsv"),
        "scoped_vs_full_sensitive_after_prefilter": partition_comparison(args.scoped_output / "node_clusters.tsv", args.full_output / "node_clusters.tsv"),
        "edge_counts": {
            "baseline": len(baseline_edges),
            "scoped": len(scoped_edges),
            "full_sensitive_after_prefilter": len(full_edges),
            "baseline_and_full_edge_sets_identical": baseline_edges == full_edges,
            "missing_from_scoped": len(missing),
        },
        "missing_edge_categories": dict(sorted(categories.items())),
        "missing_edges_inside_a_baseline_cluster": missing_within_baseline_cluster,
        "missing_baseline_representative_member_edges": missing_baseline_rep_member,
        "baseline_cluster_fragmentation_in_scoped_result": {
            "clusters_split_across_multiple_scoped_clusters": sum(x > 1 for x in split_counts),
            "mean_scoped_clusters_per_baseline_cluster": sum(split_counts) / len(split_counts),
            "maximum_scoped_clusters_for_one_baseline_cluster": max(split_counts),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
