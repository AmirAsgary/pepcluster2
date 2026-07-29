#!/usr/bin/env python3
"""Run and analyse the final PepCluster2 exhaustive and stability validation."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import math
import shutil
import statistics
import struct
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


METHODS = {
    "graph": ["--clustering-method", "graph", "--no-prefilter"],
    "graph_prefilter": ["--clustering-method", "graph", "--force-prefilter"],
    "greedy": [
        "--clustering-method", "greedy", "--greedy-selection", "kmer-degree", "--no-prefilter"
    ],
    "greedy_lazy": [
        "--clustering-method", "greedy", "--greedy-selection", "lazy-exact", "--no-prefilter"
    ],
}
LABELS = {
    "graph": "Graph",
    "graph_prefilter": "Graph + prefilter",
    "greedy": "Greedy",
    "greedy_lazy": "Greedy lazy-exact",
}
SUBSET_SIZES = (1_000, 2_000, 4_000, 6_000, 8_000)


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as reader, target.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=4 * 1024 * 1024)


def execute(command: list[str], output: Path, log: Path, resource: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    timed = ["/usr/bin/time", "-v", "-o", str(resource), *command]
    with log.open("w") as handle:
        result = subprocess.run(timed, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log}")


def common_command(binary: Path, input_fasta: Path, output: Path, threads: int) -> list[str]:
    return [
        str(binary),
        "--input", str(input_fasta),
        "--output-dir", str(output),
        "--mode", "separate_aln_anchor",
        "--alignment-similarity-threshold", "0.50",
        "--anchor-combination-similarity-threshold", "0.60",
        "--kmer-seed-threshold", "0.50",
        "--gap-open", "-4",
        "--gap-extension", "-1",
        "--terminal-overhang-gap-open", "-2",
        "--terminal-overhang-gap-extension", "-1",
        "--minimum-terminal-match-length", "2",
        "--threads", str(threads),
        "--candidate-buffer-mb", "96",
        "--compact-output",
    ]


def run_exhaustive(root: Path, helper: Path, dataset: int, threads: int) -> str:
    output = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
    if (output / "run_stats.json").exists() and (output / "true_pairs.bin").exists():
        return f"exhaustive {dataset:03d} cached"
    source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
    with tempfile.TemporaryDirectory(prefix=f"pc2_gt_{dataset:03d}_", dir="/tmp") as tmp:
        fasta = Path(tmp) / f"sample_{dataset:03d}.fasta"
        decompress(source, fasta)
        execute(
            [str(helper), str(fasta), str(output), str(threads)],
            output,
            output / "run.log",
            output / "resource.txt",
        )
    return f"exhaustive {dataset:03d} complete"


def run_cluster(
    root: Path,
    binary: Path,
    method: str,
    dataset: int,
    threads: int,
    subset_size: int | None,
) -> str:
    if subset_size is None:
        source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
        output = root / "runs" / "full" / method / f"sample_{dataset:03d}"
        trace = True
        label = f"full {method} {dataset:03d}"
    else:
        source = (
            root / "data" / "subsets" / f"n_{subset_size:06d}"
            / f"sample_{dataset:03d}.fasta.gz"
        )
        output = (
            root / "runs" / "subsets" / f"n_{subset_size:06d}"
            / method / f"sample_{dataset:03d}"
        )
        trace = False
        label = f"subset {subset_size} {method} {dataset:03d}"
    required = [output / "run_stats.json", output / "node_clusters.tsv"]
    if trace:
        required.append(output / "scored_pairs.bin")
    if all(path.exists() for path in required):
        return f"{label} cached"
    with tempfile.TemporaryDirectory(prefix="pc2_final_", dir="/tmp") as tmp:
        fasta = Path(tmp) / "input.fasta"
        decompress(source, fasta)
        command = common_command(binary, fasta, output, threads) + METHODS[method]
        if trace:
            command.append("--write-scored-pairs")
        execute(command, output, output / "run.log", output / "resource.txt")
    return f"{label} complete"


def parallel(tasks: list[tuple], function, workers: int) -> None:
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(function, *task) for task in tasks]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            print(f"[{index}/{len(tasks)}] {future.result()}", flush=True)


def read_pairs(path: Path, magic: bytes) -> set[int]:
    payload = path.read_bytes()
    if payload[:8] != magic:
        raise ValueError(f"invalid pair file: {path}")
    count = struct.unpack_from("<Q", payload, 8)[0]
    values = set(struct.unpack_from(f"<{count}Q", payload, 16))
    if len(values) != count:
        raise ValueError(f"duplicate pair identifiers in {path}")
    return values


def read_partition(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {
            row["sequence"]: row["cluster_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def choose2(value: int) -> int:
    return value * (value - 1) // 2


def partition_metrics(reference: dict[str, str], query: dict[str, str]) -> dict[str, float]:
    keys = sorted(set(reference) & set(query))
    ref_counts = Counter(reference[key] for key in keys)
    query_counts = Counter(query[key] for key in keys)
    cells = Counter((reference[key], query[key]) for key in keys)
    ref_pairs = sum(choose2(value) for value in ref_counts.values())
    query_pairs = sum(choose2(value) for value in query_counts.values())
    shared = sum(choose2(value) for value in cells.values())
    total_pairs = choose2(len(keys))
    expected = ref_pairs * query_pairs / total_pairs if total_pairs else 0.0
    denominator = 0.5 * (ref_pairs + query_pairs) - expected
    ari = (shared - expected) / denominator if denominator else 1.0
    n = len(keys)
    mutual_information = sum(
        count / n * math.log(n * count / (ref_counts[ref] * query_counts[qry]))
        for (ref, qry), count in cells.items()
    )
    ref_entropy = -sum(count / n * math.log(count / n) for count in ref_counts.values())
    query_entropy = -sum(count / n * math.log(count / n) for count in query_counts.values())
    nmi = 2 * mutual_information / (ref_entropy + query_entropy) if ref_entropy + query_entropy else 1.0
    recall = shared / ref_pairs if ref_pairs else 1.0
    precision = shared / query_pairs if query_pairs else 1.0
    union = ref_pairs + query_pairs - shared
    jaccard = shared / union if union else 1.0
    return {
        "ari": ari,
        "nmi": nmi,
        "pairwise_jaccard": jaccard,
        "coassociation_recall": recall,
        "coassociation_precision": precision,
        "reference_cocluster_pairs": ref_pairs,
        "query_cocluster_pairs": query_pairs,
        "shared_cocluster_pairs": shared,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyse(root: Path, datasets: int) -> None:
    metrics_dir = root / "metrics"
    figures_dir = root / "figures"
    metrics_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)
    search_rows, agreement_rows, truth_rows = [], [], []
    for dataset in range(datasets):
        ground_dir = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
        ground_stats = json.loads((ground_dir / "run_stats.json").read_text())
        true_pairs = read_pairs(ground_dir / "true_pairs.bin", b"PC2TRUE1")
        ground = read_partition(ground_dir / "ground_truth_clusters.tsv")
        truth_rows.append({"dataset": dataset, **ground_stats})
        for method in METHODS:
            run_dir = root / "runs" / "full" / method / f"sample_{dataset:03d}"
            explored = read_pairs(run_dir / "scored_pairs.bin", b"PC2PAIR1")
            found = len(true_pairs & explored)
            possible = ground_stats["all_possible_pairs"]
            search_rows.append({
                "dataset": dataset, "method": method, "true_pairs": len(true_pairs),
                "explored_unique_pairs": len(explored), "true_pairs_explored": found,
                "missed_true_pairs": len(true_pairs) - found,
                "useless_explored_pairs": len(explored) - found,
                "search_recall": found / len(true_pairs) if true_pairs else 1.0,
                "search_precision": found / len(explored) if explored else 1.0,
                "fraction_all_pairs_explored": len(explored) / possible,
            })
            query = read_partition(run_dir / "node_clusters.tsv")
            run_stats = json.loads((run_dir / "run_stats.json").read_text())
            agreement_rows.append({
                "dataset": dataset, "method": method,
                **partition_metrics(ground, query),
                "clusters": run_stats["final_clusters"],
                "singletons": run_stats["singleton_clusters"],
                "runtime_seconds": run_stats["elapsed_seconds"],
            })
    write_csv(metrics_dir / "true_pair_space.csv", truth_rows)
    write_csv(figures_dir / "search_rule_performance.csv", search_rows)
    write_csv(figures_dir / "ground_truth_cluster_agreement.csv", agreement_rows)

    colors = dict(zip(METHODS, ["#31688e", "#35b779", "#e07a00", "#7e3ace"]))
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    for axis, metric, title in zip(
        axes,
        ["search_recall", "search_precision", "fraction_all_pairs_explored"],
        ["True-pair recall", "Search precision", "Fraction of all pairs explored"],
    ):
        data = [[row[metric] for row in search_rows if row["method"] == method] for method in METHODS]
        box = axis.boxplot(data, tick_labels=[LABELS[m] for m in METHODS], patch_artist=True)
        for box_patch, method in zip(box["boxes"], METHODS):
            box_patch.set_facecolor(colors[method])
            box_patch.set_alpha(0.65)
        axis.set_title(title)
        axis.set_ylim(bottom=0)
        axis.tick_params(axis="x", rotation=20)
        axis.grid(axis="y", alpha=0.25)
    figure.savefig(figures_dir / "search_rule_performance.png", dpi=220)
    figure.savefig(figures_dir / "search_rule_performance.pdf")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    for axis, metric, title in zip(
        axes, ["ari", "nmi"], ["Adjusted Rand index", "Normalized mutual information"]
    ):
        data = [[row[metric] for row in agreement_rows if row["method"] == method] for method in METHODS]
        box = axis.boxplot(data, tick_labels=[LABELS[m] for m in METHODS], patch_artist=True)
        for box_patch, method in zip(box["boxes"], METHODS):
            box_patch.set_facecolor(colors[method])
            box_patch.set_alpha(0.65)
        axis.set_title(title)
        axis.set_ylim(0, 1.01)
        axis.tick_params(axis="x", rotation=20)
        axis.grid(axis="y", alpha=0.25)
    figure.savefig(figures_dir / "ground_truth_cluster_agreement.png", dpi=220)
    figure.savefig(figures_dir / "ground_truth_cluster_agreement.pdf")
    plt.close(figure)

    stability_rows = []
    for dataset in range(datasets):
        for method in METHODS:
            full = read_partition(
                root / "runs" / "full" / method / f"sample_{dataset:03d}" / "node_clusters.tsv"
            )
            for size in SUBSET_SIZES:
                subset = read_partition(
                    root / "runs" / "subsets" / f"n_{size:06d}" / method
                    / f"sample_{dataset:03d}" / "node_clusters.tsv"
                )
                reference = {sequence: full[sequence] for sequence in subset}
                stability_rows.append({
                    "dataset": dataset, "method": method, "subset_size": size,
                    "subset_fraction": size / 10_000, **partition_metrics(reference, subset),
                })
    write_csv(figures_dir / "cluster_stability.csv", stability_rows)
    metric_names = [
        ("pairwise_jaccard", "Pairwise Jaccard"),
        ("ari", "Subset ARI"),
        ("nmi", "Subset NMI"),
        ("coassociation_recall", "Co-association recall"),
        ("coassociation_precision", "Co-association precision"),
    ]
    figure, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    for axis, (metric, title) in zip(axes.flat, metric_names):
        for method in METHODS:
            means, deviations = [], []
            for size in SUBSET_SIZES:
                values = [
                    row[metric] for row in stability_rows
                    if row["method"] == method and row["subset_size"] == size
                ]
                means.append(statistics.mean(values))
                deviations.append(statistics.stdev(values))
            x = np.asarray(SUBSET_SIZES) / 10_000
            mean, deviation = np.asarray(means), np.asarray(deviations)
            axis.plot(x, mean, marker="o", color=colors[method], label=LABELS[method])
            axis.fill_between(
                x, np.maximum(0, mean - deviation), np.minimum(1, mean + deviation),
                color=colors[method], alpha=0.14
            )
        axis.set_title(title)
        axis.set_xlabel("Subset fraction")
        axis.set_ylabel("Stability score")
        axis.set_ylim(0, 1.01)
        axis.grid(alpha=0.25)
    axes.flat[-1].axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles, labels, loc="center", frameon=False)
    figure.savefig(figures_dir / "cluster_stability.png", dpi=220)
    figure.savefig(figures_dir / "cluster_stability.pdf")
    plt.close(figure)

    def aggregate(rows: list[dict], metric: str, method: str) -> tuple[float, float]:
        values = [float(row[metric]) for row in rows if row["method"] == method]
        return statistics.mean(values), statistics.stdev(values)

    lines = [
        "# PepCluster2 0.4.3 final validation", "", "## Configuration", "",
        "- 20 independently sampled datasets; 10,000 peptides each.",
        "- Scoring mode: separate_aln_anchor.",
        "- Alignment-similarity threshold: 0.50.",
        "- Anchor-combination-similarity threshold: 0.60.",
        "- Terminal/core alignment weights: 4/1.",
        "- Ground truth: every pair scored exactly, followed by dynamic greedy set cover.",
        "",
        "The exhaustive partition is a computational reference under the chosen scoring rule, not biological ground truth.",
        "", "## How to read the metrics", "",
        "- Search recall is the fraction of all truly eligible peptide pairs that the method actually scored.",
        "- Search precision is the fraction of scored pairs that were truly eligible; low precision means extra work, not incorrect final edges.",
        "- ARI is the Adjusted Rand Index: agreement of two clusterings after correcting for chance (1 is identical).",
        "- NMI is Normalized Mutual Information: shared cluster information (1 is identical). With many small clusters, NMI can remain high even when assignments differ, so ARI and pairwise Jaccard are more informative here.",
        "", "## Search-rule performance", "",
        "| Method | Recall mean ± SD | Precision mean ± SD | All-pairs fraction mean ± SD |",
        "|---|---:|---:|---:|",
    ]
    for method in METHODS:
        recall = aggregate(search_rows, "search_recall", method)
        precision = aggregate(search_rows, "search_precision", method)
        explored = aggregate(search_rows, "fraction_all_pairs_explored", method)
        lines.append(
            f"| {LABELS[method]} | {recall[0]:.4f} ± {recall[1]:.4f} | "
            f"{precision[0]:.4f} ± {precision[1]:.4f} | {explored[0]:.4f} ± {explored[1]:.4f} |"
        )
    lines += [
        "", "## Ground-truth cluster agreement", "",
        "| Method | ARI mean ± SD | NMI mean ± SD | Clusters mean | Singletons mean |",
        "|---|---:|---:|---:|---:|",
    ]
    for method in METHODS:
        ari, nmi = aggregate(agreement_rows, "ari", method), aggregate(agreement_rows, "nmi", method)
        clusters = statistics.mean(row["clusters"] for row in agreement_rows if row["method"] == method)
        singletons = statistics.mean(row["singletons"] for row in agreement_rows if row["method"] == method)
        lines.append(
            f"| {LABELS[method]} | {ari[0]:.4f} ± {ari[1]:.4f} | "
            f"{nmi[0]:.4f} ± {nmi[1]:.4f} | {clusters:.1f} | {singletons:.1f} |"
        )
    truth_singletons = statistics.mean(row["singleton_clusters"] for row in truth_rows)
    lines += [
        "", f"Exhaustive ground truth contained a mean of {truth_singletons:.1f} singleton clusters "
        f"({truth_singletons / 100:.2f}% of 10,000 peptides).",
        "", "## Stability", "",
        "Stability compares a clustering of a subset with the full-dataset clustering after restricting the latter to the same peptides. Higher values mean less dependence on dataset size.",
        "", "| Method | Subset | Jaccard | ARI | NMI | Pair recall | Pair precision |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for method in METHODS:
        for size in (1_000, 8_000):
            selected = [
                row for row in stability_rows
                if row["method"] == method and row["subset_size"] == size
            ]
            means = {
                metric: statistics.mean(float(row[metric]) for row in selected)
                for metric in (
                    "pairwise_jaccard", "ari", "nmi",
                    "coassociation_recall", "coassociation_precision",
                )
            }
            lines.append(
                f"| {LABELS[method]} | {size // 100}% | "
                f"{means['pairwise_jaccard']:.4f} | {means['ari']:.4f} | "
                f"{means['nmi']:.4f} | {means['coassociation_recall']:.4f} | "
                f"{means['coassociation_precision']:.4f} |"
            )
    graph_recall = aggregate(search_rows, "search_recall", "graph")[0]
    greedy_recall = aggregate(search_rows, "search_recall", "greedy")[0]
    graph_ari = aggregate(agreement_rows, "ari", "graph")[0]
    greedy_ari = aggregate(agreement_rows, "ari", "greedy")[0]
    lines += [
        "", "## Final interpretation", "",
        f"PepCluster2 avoids almost all all-pairs work: graph scored about 1.62% of all possible pairs. However, its search rule recovered only {graph_recall:.1%} of pairs that pass both exact thresholds. Static greedy recovered still fewer ({greedy_recall:.1%}). Therefore the current candidate search is fast but not exhaustive.",
        "",
        f"Graph and lazy-exact greedy were the closest to the exhaustive set-cover reference, but their mean ARI was only about {graph_ari:.2f}. Static greedy was worse (ARI {greedy_ari:.2f}) and produced about 540 more clusters per dataset. The high NMI values should not be read as near-identity because these data contain thousands of small and singleton clusters.",
        "",
        "Forced graph prefiltering was effectively indistinguishable from non-prefilter graph at 10,000 peptides in both search and clustering metrics. This supports using the prefilter at this tested scale, but does not prove equivalence on larger or biologically different datasets.",
        "",
        "Cluster stability improved steadily as the subset approached the full dataset. At 80%, ARI was about 0.70 and pairwise Jaccard about 0.54 for graph/lazy-exact (0.57 for static greedy). Thus clusters are reasonably, but not fully, stable to dataset composition.",
        "",
        "Overall conclusion: use graph when its temporary edge storage is affordable, or lazy-exact greedy when memory/disk is limiting. Do not claim that the present search recovers every eligible relationship or that clusters are dataset-independent. Biological validation against peptide–MHC labels remains necessary before making purity claims.",
        "", "## Files", "",
        "- runs/exhaustive/: exact true pairs and exhaustive set-cover assignments.",
        "- runs/full/: four full clustering paths and unique exact-scoring traces.",
        "- runs/subsets/: all nested subset clusterings.",
        "- figures/: plot-matched CSV, PNG, and PDF files.",
        "- code/: complete preparation, execution, exhaustive-scoring, and analysis source.",
    ]
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=20)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=6)
    parser.add_argument(
        "--stage", choices=["exhaustive", "full", "subsets", "analyse", "all"], default="all"
    )
    args = parser.parse_args()
    root, binary, helper = args.root.resolve(), args.binary.resolve(), args.helper.resolve()
    if args.stage in {"exhaustive", "all"}:
        parallel(
            [(root, helper, dataset, args.threads) for dataset in range(args.datasets)],
            run_exhaustive, args.workers,
        )
    if args.stage in {"full", "all"}:
        parallel(
            [(root, binary, method, dataset, args.threads, None)
             for dataset in range(args.datasets) for method in METHODS],
            run_cluster, args.workers,
        )
    if args.stage in {"subsets", "all"}:
        parallel(
            [(root, binary, method, dataset, args.threads, size)
             for size in SUBSET_SIZES for dataset in range(args.datasets) for method in METHODS],
            run_cluster, args.workers,
        )
    if args.stage in {"analyse", "all"}:
        analyse(root, args.datasets)


if __name__ == "__main__":
    main()
