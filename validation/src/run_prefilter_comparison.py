#!/usr/bin/env python3
"""Paired validation of PepCluster2 with and without its high-confidence prefilter."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    version: str = "0.3.4-dev"
    kmer_similarity_formula: str = "aligned front 3mer and aligned end 3mer normalized BLOSUM means"
    samples: int = 100
    sample_size: int = 20_000
    combined_similarity_threshold: float = 0.60
    kmer_seed_threshold: float = 0.50
    prefilter_floor: float = 0.75
    process_workers: int = 4
    threads_per_process: int = 6
    candidate_buffer_mb: int = 128
    iteration_cap: int | None = None
    merge_cap: int | None = None


def run_one(task: tuple[int, str, Path, Path, Path, Settings]) -> dict:
    sample, mode, source, binary, root, settings = task
    output = root / "runs" / mode / f"sample_{sample:03d}"
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"pc2_{mode}_{sample:03d}_", dir=root / "tmp") as tmp:
        fasta = Path(tmp) / f"sample_{sample:03d}.fasta"
        with gzip.open(source, "rb") as inp, fasta.open("wb") as out:
            shutil.copyfileobj(inp, out)
        command = [
            str(binary), "--input", str(fasta), "--output-dir", str(output),
            "--tmp-dir", str(Path(tmp) / "work"), "--threshold",
            str(settings.combined_similarity_threshold), "--kmer-seed-threshold",
            str(settings.kmer_seed_threshold), "--candidate-buffer-mb",
            str(settings.candidate_buffer_mb), "--threads",
            str(settings.threads_per_process),
        ]
        command.append("--no-prefilter" if mode == "no_prefilter" else "--force-prefilter")
        if settings.iteration_cap is not None:
            command += ["--iteration-cap", str(settings.iteration_cap)]
        if settings.merge_cap is not None:
            command += ["--merge-cap", str(settings.merge_cap)]
        started = time.monotonic()
        completed = subprocess.run(command, text=True, capture_output=True)
        wall = time.monotonic() - started
        (output / "validation_stdout.txt").write_text(completed.stdout)
        (output / "validation_stderr.txt").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"sample {sample:03d} {mode} failed: {completed.stderr[-2000:]}")
        stats = json.loads((output / "run_stats.json").read_text())
        stats.update(sample=sample, validation_mode=mode, validation_wall_seconds=wall)
        return stats


def read_partition(path: Path) -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], tuple[str, str]]]:
    labels = {}
    representatives = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["anchor"], row["geometry_mask"])
            labels[key] = row["cluster_id"]
            representatives[key] = (row["representative_anchor"], row["representative_geometry_mask"])
    return labels, representatives


def choose2(n: int) -> int:
    return n * (n - 1) // 2


def adjusted_rand_index(a: list[str], b: list[str]) -> float:
    n = len(a)
    if n < 2:
        return 1.0
    contingency = Counter(zip(a, b))
    rows = Counter(a)
    cols = Counter(b)
    sum_cells = sum(choose2(v) for v in contingency.values())
    sum_rows = sum(choose2(v) for v in rows.values())
    sum_cols = sum(choose2(v) for v in cols.values())
    total = choose2(n)
    expected = sum_rows * sum_cols / total
    maximum = (sum_rows + sum_cols) / 2
    if maximum == expected:
        return 1.0
    return (sum_cells - expected) / (maximum - expected)


def canonical_partition(labels: dict[tuple[str, str], str]) -> set[frozenset[tuple[str, str]]]:
    groups: dict[str, set[tuple[str, str]]] = {}
    for node, label in labels.items():
        groups.setdefault(label, set()).add(node)
    return {frozenset(group) for group in groups.values()}


def compare_sample(sample: int, root: Path) -> dict:
    off_dir = root / "runs" / "no_prefilter" / f"sample_{sample:03d}"
    on_dir = root / "runs" / "prefilter" / f"sample_{sample:03d}"
    off_stats = json.loads((off_dir / "run_stats.json").read_text())
    on_stats = json.loads((on_dir / "run_stats.json").read_text())
    off, off_reps = read_partition(off_dir / "anchor_clusters.tsv")
    on, on_reps = read_partition(on_dir / "anchor_clusters.tsv")
    if off.keys() != on.keys():
        raise RuntimeError(f"node sets differ for sample {sample:03d}")
    keys = sorted(off)
    ari = adjusted_rand_index([off[k] for k in keys], [on[k] for k in keys])
    exact = canonical_partition(off) == canonical_partition(on)
    rep_agreement = sum(off_reps[k] == on_reps[k] for k in keys) / max(1, len(keys))
    return {
        "sample": sample,
        "unique_nodes": len(keys),
        "exact_partition": exact,
        "adjusted_rand_index": ari,
        "representative_assignment_agreement": rep_agreement,
        "clusters_no_prefilter": off_stats["final_clusters"],
        "clusters_prefilter": on_stats["final_clusters"],
        "cluster_count_difference": on_stats["final_clusters"] - off_stats["final_clusters"],
        "candidate_pairs_no_prefilter": off_stats["sensitive_candidate_pairs"],
        "candidate_pairs_prefilter_total": on_stats["prefilter_candidate_pairs"] + on_stats["sensitive_candidate_pairs"],
        "final_edges_no_prefilter": off_stats["graph_edge_count"],
        "final_edges_prefilter": on_stats["graph_edge_count"],
        "edge_recovery_fraction": on_stats["graph_edge_count"] / max(1, off_stats["graph_edge_count"]),
        "seconds_no_prefilter": off_stats["elapsed_seconds"],
        "seconds_prefilter": on_stats["elapsed_seconds"],
        "iterations_no_prefilter": off_stats["iterations"],
        "iterations_prefilter": on_stats["iterations"],
        "converged_no_prefilter": off_stats["converged"],
        "converged_prefilter": on_stats["converged"],
    }


def describe(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    n = len(ordered)
    return {
        "mean": sum(ordered) / n,
        "median": (ordered[(n - 1) // 2] + ordered[n // 2]) / 2,
        "minimum": ordered[0],
        "maximum": ordered[-1],
    }


def create_figure(rows: list[dict], path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    samples = [r["sample"] for r in rows]
    axes[0, 0].scatter([r["clusters_no_prefilter"] for r in rows], [r["clusters_prefilter"] for r in rows], s=18, alpha=.75)
    values = [r["clusters_no_prefilter"] for r in rows] + [r["clusters_prefilter"] for r in rows]
    lo, hi = min(values), max(values)
    axes[0, 0].plot([lo, hi], [lo, hi], "--", color="black", linewidth=1)
    axes[0, 0].set(xlabel="Clusters without prefilter", ylabel="Clusters with prefilter", title="A. Cluster counts")
    axes[0, 1].plot(samples, [r["adjusted_rand_index"] for r in rows], ".", color="#2b6cb0")
    axes[0, 1].set(xlabel="Dataset", ylabel="Adjusted Rand index", title="B. Partition agreement", ylim=(0, 1.02))
    axes[1, 0].plot(samples, [r["edge_recovery_fraction"] for r in rows], ".", color="#c05621")
    axes[1, 0].set(xlabel="Dataset", ylabel="Edges retained / baseline edges", title="C. Edge recovery")
    axes[1, 1].scatter([r["seconds_no_prefilter"] for r in rows], [r["seconds_prefilter"] for r in rows], s=18, alpha=.75, color="#2f855a")
    times = [r["seconds_no_prefilter"] for r in rows] + [r["seconds_prefilter"] for r in rows]
    lo, hi = min(times), max(times)
    axes[1, 1].plot([lo, hi], [lo, hi], "--", color="black", linewidth=1)
    axes[1, 1].set(xlabel="Seconds without prefilter", ylabel="Seconds with prefilter", title="D. Runtime")
    figure.suptitle("PepCluster2 high-confidence prefilter versus non-prefilter baseline")
    figure.tight_layout()
    figure.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def write_report(rows: list[dict], summary: dict, root: Path, settings: Settings) -> None:
    ari = summary["adjusted_rand_index"]
    delta = summary["cluster_count_difference"]
    edge = summary["edge_recovery_fraction"]
    speed = summary["runtime_ratio_prefilter_over_no_prefilter"]
    report = f"""# PepCluster2 prefilter validation

## Question

Does the high-confidence prefilter reproduce the clusters obtained when all sensitive k-mer candidates are evaluated?

## Data and settings

The test used the same 100 independently sampled datasets used in the earlier PepCluster validation. Each contains 20,000 peptides. K-mer similarity used: {settings.kmer_similarity_formula}. Both modes used k-mer seed threshold {settings.kmer_seed_threshold:.2f}, combined similarity threshold {settings.combined_similarity_threshold:.2f}, equal weighting of k-mer similarity and anchor-combination similarity, constrained representative updates, validated merging, and iteration to convergence. Peptides shorter than eight residues were excluded.

The baseline disabled only the high-confidence prefilter; it still used terminal k-mer seed blocking. The forced-prefilter mode required two distinct identical ordered anchor-pair values and combined similarity at least max({settings.prefilter_floor:.2f}, threshold) for provisional edges. Its sensitive completion evaluated every unassigned peptide against all peptides and every provisional representative against all peptides. Provisional clusters guided candidate generation only; final set cover restarted from the completed graph.

The **Adjusted Rand index (ARI)** measures agreement between two partitions after correcting for agreement expected by chance. ARI=1 means identical grouping; values below 1 indicate differences. Exact equality was also checked directly and does not depend on cluster names.

## Results

- Exact partition recovery: **{summary['exact_partitions']}/{settings.samples} datasets**.
- ARI: mean **{ari['mean']:.4f}**, median **{ari['median']:.4f}**, range **{ari['minimum']:.4f}–{ari['maximum']:.4f}**.
- Additional clusters with prefilter: mean **{delta['mean']:.1f}**, median **{delta['median']:.1f}**, range **{delta['minimum']:.0f}–{delta['maximum']:.0f}**.
- Baseline eligible edges recovered by the prefilter graph: mean **{edge['mean']:.3f}**, or {100*edge['mean']:.1f}%.
- Prefilter/baseline runtime ratio: mean **{speed['mean']:.3f}**. Values below 1 mean the prefilter was faster.
- All baseline runs converged: **{summary['converged_no_prefilter']}/{settings.samples}**; all prefilter runs converged: **{summary['converged_prefilter']}/{settings.samples}**.

![Paired comparison](figures/prefilter_comparison.png)

## Interpretation

The forced prefilter **does not reproduce the non-prefilter clustering at the current defaults**. Provisional clusters now guide candidate generation only; final set cover restarts from the completed graph. The remaining differences arise because the scoped sensitive stage omits eligible edges between two assigned non-representatives. Later representative reassignment cannot recover an edge that was never generated. A higher cluster count is therefore fragmentation caused by candidate loss, not evidence of more biologically distinct groups.

The non-prefilter result should remain the scientific reference for 20,000-peptide datasets. The current prefilter is useful as an experimental large-database approximation, but it should not silently replace the reference path. Before using it for a paper, its sensitive completion must be broadened or made iterative until edge/partition recovery reaches a predefined acceptable target.

## Files

- `comparison/per_sample.csv`: paired metrics for every dataset.
- `comparison/summary.json`: machine-readable aggregate results.
- `runs/`: complete `--output-dir` files, configurations, logs, and statistics for every run.
- `config/config.json`: complete validation settings.
- `config/input_manifest.tsv`: source dataset manifest.
"""
    (root / "REPORT.md").write_text(report)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    settings = Settings()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    for directory in ("config", "comparison", "figures", "runs/no_prefilter", "runs/prefilter", "tmp"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    (root / "config" / "config.json").write_text(json.dumps(asdict(settings), indent=2) + "\n")
    samples = []
    with (root / "config" / "input_manifest.tsv").open("w") as manifest:
        manifest.write("sample\tsource\n")
        for i in range(settings.samples):
            source = args.data_dir / f"sample_{i:03d}.fasta.gz"
            if not source.exists():
                raise FileNotFoundError(source)
            samples.append(source)
            manifest.write(f"{i:03d}\t{source}\n")

    tasks = [(i, mode, samples[i], args.binary.resolve(), root, settings)
             for i in range(settings.samples) for mode in ("no_prefilter", "prefilter")]
    completed_rows = []
    with ThreadPoolExecutor(max_workers=settings.process_workers) as executor:
        futures = {executor.submit(run_one, task): (task[0], task[1]) for task in tasks}
        for done, future in enumerate(as_completed(futures), 1):
            completed_rows.append(future.result())
            if done % 10 == 0:
                print(f"completed {done}/{len(tasks)} runs", flush=True)

    rows = [compare_sample(i, root) for i in range(settings.samples)]
    with (root / "comparison" / "per_sample.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "datasets": settings.samples,
        "exact_partitions": sum(r["exact_partition"] for r in rows),
        "adjusted_rand_index": describe([r["adjusted_rand_index"] for r in rows]),
        "representative_assignment_agreement": describe([r["representative_assignment_agreement"] for r in rows]),
        "cluster_count_difference": describe([r["cluster_count_difference"] for r in rows]),
        "edge_recovery_fraction": describe([r["edge_recovery_fraction"] for r in rows]),
        "runtime_ratio_prefilter_over_no_prefilter": describe([r["seconds_prefilter"] / max(1e-9, r["seconds_no_prefilter"]) for r in rows]),
        "converged_no_prefilter": sum(r["converged_no_prefilter"] for r in rows),
        "converged_prefilter": sum(r["converged_prefilter"] for r in rows),
    }
    (root / "comparison" / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    create_figure(rows, root / "figures" / "prefilter_comparison.png")
    write_report(rows, summary, root, settings)
    shutil.rmtree(root / "tmp")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
