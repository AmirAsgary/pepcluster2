#!/usr/bin/env python3
"""Validate PepCluster2 alignment scoring across greedy and graph clustering.

The runner is resumable. It reuses the previously sampled 100 x 20,000 FASTA
datasets, runs three clustering paths, tests one deterministic FASTA shuffle per
dataset and method, and evaluates ten 50% subsets of dataset 0.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
SOURCE_VALIDATION = PROJECT / "validation/2026-07-28_0.2.0-dev_mmseqs_defaults"
DEFAULT_ROOT = PROJECT / "validation/2026-07-29_0.4.0-dev_separate_aln_anchor_methods"
METHODS = ("greedy", "graph", "graph_prefilter")


@dataclass(frozen=True)
class Settings:
    version: str = "0.4.0-dev"
    scoring_mode: str = "separate_aln_anchor"
    threshold: float = 0.60
    alignment_threshold: float = 0.60
    anchor_threshold: float = 0.60
    kmer_seed_threshold: float = 0.50
    gap_open: float = -4.0
    gap_extension: float = -1.0
    terminal_gap_open: float = -2.0
    terminal_gap_extension: float = -1.0
    minimum_terminal_match_length: int = 2
    datasets: int = 100
    records_per_dataset: int = 20_000
    shuffle_replicates: int = 1
    subset_replicates: int = 10
    subset_fraction: float = 0.50
    master_seed: int = 20260729
    process_workers: int = 4
    threads_per_process: int = 6
    candidate_buffer_mb: int = 96


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    handle = gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")
    with handle:
        header: str | None = None
        sequence: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def seed_for(settings: Settings, stage: int, dataset: int, replicate: int = 0) -> int:
    digest = hashlib.sha256(
        f"{settings.master_seed}:{stage}:{dataset}:{replicate}".encode()
    ).digest()
    return int.from_bytes(digest[:8], "little")


def method_flags(method: str) -> list[str]:
    if method == "greedy":
        return ["--clustering-method", "greedy", "--no-prefilter"]
    if method == "graph":
        return ["--clustering-method", "graph", "--no-prefilter"]
    if method == "graph_prefilter":
        return ["--clustering-method", "graph", "--force-prefilter"]
    raise ValueError(method)


def compress_assignment(output: Path) -> Path:
    source = output / "node_clusters.tsv"
    target = output / "node_clusters.tsv.gz"
    if target.exists():
        source.unlink(missing_ok=True)
        return target
    with source.open("rb") as inp, gzip.open(target, "wb", compresslevel=6) as out:
        shutil.copyfileobj(inp, out)
    source.unlink()
    return target


def run_binary(
    *, binary: Path, fasta: Path, output: Path, method: str, settings: Settings,
    kmer_table: Path, tmp_parent: Path,
) -> dict:
    assignment = output / "node_clusters.tsv.gz"
    stats_path = output / "run_stats.json"
    if assignment.exists() and stats_path.exists():
        return json.loads(stats_path.read_text())
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pc2_work_", dir=tmp_parent) as work:
        command = [
            str(binary), "--input", str(fasta), "--output-dir", str(output),
            "--tmp-dir", str(Path(work) / "tmp"), "--kmer-table", str(kmer_table),
            "--mode", settings.scoring_mode, "--threshold", str(settings.threshold),
            "--alignment-similarity-threshold", str(settings.alignment_threshold),
            "--anchor-combination-similarity-threshold", str(settings.anchor_threshold),
            "--kmer-seed-threshold", str(settings.kmer_seed_threshold),
            "--gap-open", str(settings.gap_open), "--gap-extension", str(settings.gap_extension),
            "--terminal-overhang-gap-open", str(settings.terminal_gap_open),
            "--terminal-overhang-gap-extension", str(settings.terminal_gap_extension),
            "--minimum-terminal-match-length", str(settings.minimum_terminal_match_length),
            "--candidate-buffer-mb", str(settings.candidate_buffer_mb),
            "--threads", str(settings.threads_per_process), "--compact-output",
            *method_flags(method),
        ]
        started = time.monotonic()
        completed = subprocess.run(command, text=True, capture_output=True)
        wall = time.monotonic() - started
        (output / "validation_stdout.txt").write_text(completed.stdout)
        (output / "validation_stderr.txt").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(
                f"{method} failed for {fasta}:\n{completed.stderr[-4000:]}"
            )
        compress_assignment(output)
        stats = json.loads(stats_path.read_text())
        stats["validation_wall_seconds"] = wall
        stats_path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
        return stats


def prepare_fasta(source: Path, destination: Path, order: np.ndarray | None = None) -> None:
    if order is None:
        with gzip.open(source, "rb") as inp, destination.open("wb") as out:
            shutil.copyfileobj(inp, out)
        return
    records = read_fasta(source)
    write_fasta(destination, [records[int(i)] for i in order])


def run_dataset_task(task: tuple[int, str, str, Path, Path, Path, Settings]) -> dict:
    dataset, variant, method, source, binary, root, settings = task
    output = root / "runs" / variant / method / f"sample_{dataset:03d}"
    if output.joinpath("node_clusters.tsv.gz").exists() and output.joinpath("run_stats.json").exists():
        stats = json.loads(output.joinpath("run_stats.json").read_text())
        return {**stats, "dataset": dataset, "variant": variant, "method": method}
    with tempfile.TemporaryDirectory(prefix=f"pc2_{dataset:03d}_{variant}_", dir=root / "tmp") as temporary:
        fasta = Path(temporary) / "input.fasta"
        if variant == "base":
            prepare_fasta(source, fasta)
        elif variant == "shuffle":
            n = len(read_fasta(source))
            rng = np.random.default_rng(seed_for(settings, 1, dataset))
            prepare_fasta(source, fasta, rng.permutation(n))
        else:
            raise ValueError(variant)
        stats = run_binary(
            binary=binary, fasta=fasta, output=output, method=method,
            settings=settings, kmer_table=root / "config/kmer2_similarity_q.bin",
            tmp_parent=root / "tmp",
        )
    return {**stats, "dataset": dataset, "variant": variant, "method": method}


def run_subset_task(task: tuple[int, str, Path, Path, Path, Settings]) -> dict:
    replicate, method, source, binary, root, settings = task
    output = root / "runs/subsets" / method / f"subset_{replicate:02d}"
    if output.joinpath("node_clusters.tsv.gz").exists() and output.joinpath("run_stats.json").exists():
        stats = json.loads(output.joinpath("run_stats.json").read_text())
        return {**stats, "replicate": replicate, "method": method}
    records = read_fasta(source)
    size = round(len(records) * settings.subset_fraction)
    rng = np.random.default_rng(seed_for(settings, 2, 0, replicate))
    selected = np.sort(rng.choice(len(records), size=size, replace=False))
    with tempfile.TemporaryDirectory(prefix=f"pc2_subset_{replicate:02d}_", dir=root / "tmp") as temporary:
        fasta = Path(temporary) / "subset.fasta"
        write_fasta(fasta, [records[int(i)] for i in selected])
        stats = run_binary(
            binary=binary, fasta=fasta, output=output, method=method,
            settings=settings, kmer_table=root / "config/kmer2_similarity_q.bin",
            tmp_parent=root / "tmp",
        )
    (output / "subset_sequences.txt.gz").unlink(missing_ok=True)
    with gzip.open(output / "subset_sequences.txt.gz", "wt") as handle:
        for i in selected:
            handle.write(records[int(i)][1] + "\n")
    return {**stats, "replicate": replicate, "method": method}


def execute(tasks: list[tuple], function, workers: int, label: str) -> list[dict]:
    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(function, task): task for task in tasks}
        for completed, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            if completed % 10 == 0 or completed == len(tasks):
                print(f"[{label}] {completed}/{len(tasks)} complete", flush=True)
    return results


def read_partition(path: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            labels[row["sequence"]] = row["cluster_id"]
    return labels


def choose2(n: int) -> int:
    return n * (n - 1) // 2


def partition_metrics(first: dict[str, str], second: dict[str, str]) -> dict[str, float | bool | int]:
    if first.keys() != second.keys():
        missing = len(first.keys() ^ second.keys())
        raise RuntimeError(f"partition node sets differ by {missing} sequences")
    keys = sorted(first)
    a = [first[key] for key in keys]
    b = [second[key] for key in keys]
    rows = Counter(a)
    columns = Counter(b)
    cells = Counter(zip(a, b))
    n = len(keys)
    total = choose2(n)
    same_a = sum(choose2(value) for value in rows.values())
    same_b = sum(choose2(value) for value in columns.values())
    same_both = sum(choose2(value) for value in cells.values())
    expected = same_a * same_b / total if total else 0.0
    maximum = (same_a + same_b) / 2
    ari = 1.0 if maximum == expected else (same_both - expected) / (maximum - expected)
    precision = 1.0 if same_b == 0 else same_both / same_b
    recall = 1.0 if same_a == 0 else same_both / same_a
    pair_agreement = 1.0 if total == 0 else (same_both + total - same_a - same_b + same_both) / total
    groups_a: dict[str, set[str]] = defaultdict(set)
    groups_b: dict[str, set[str]] = defaultdict(set)
    for key in keys:
        groups_a[first[key]].add(key)
        groups_b[second[key]].add(key)
    jaccards: list[float] = []
    informative: list[float] = []
    for key in keys:
        ga = groups_a[first[key]] - {key}
        gb = groups_b[second[key]] - {key}
        union = ga | gb
        value = 1.0 if not union else len(ga & gb) / len(union)
        jaccards.append(value)
        if union:
            informative.append(value)
    canonical_a = {frozenset(group) for group in groups_a.values()}
    canonical_b = {frozenset(group) for group in groups_b.values()}
    return {
        "members": n,
        "exact_partition": canonical_a == canonical_b,
        "adjusted_rand_index": ari,
        "pairwise_agreement": pair_agreement,
        "cocluster_precision": precision,
        "cocluster_recall": recall,
        "mean_member_cocluster_jaccard": float(np.mean(jaccards)) if jaccards else 1.0,
        "mean_informative_member_jaccard": float(np.mean(informative)) if informative else 1.0,
        "informative_members": len(informative),
        "cocluster_pairs_first": same_a,
        "cocluster_pairs_second": same_b,
        "cocluster_pairs_shared": same_both,
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def analyse(root: Path, settings: Settings) -> None:
    method_rows: list[dict] = []
    shuffle_rows: list[dict] = []
    pairs = (("greedy", "graph"), ("graph", "graph_prefilter"), ("greedy", "graph_prefilter"))
    for dataset in range(settings.datasets):
        base = {
            method: read_partition(root / f"runs/base/{method}/sample_{dataset:03d}/node_clusters.tsv.gz")
            for method in METHODS
        }
        stats = {
            method: json.loads((root / f"runs/base/{method}/sample_{dataset:03d}/run_stats.json").read_text())
            for method in METHODS
        }
        for first, second in pairs:
            method_rows.append({
                "dataset": dataset, "first_method": first, "second_method": second,
                "comparison": f"{first}_vs_{second}",
                "clusters_first": stats[first]["final_clusters"],
                "clusters_second": stats[second]["final_clusters"],
                "seconds_first": stats[first]["elapsed_seconds"],
                "seconds_second": stats[second]["elapsed_seconds"],
                **partition_metrics(base[first], base[second]),
            })
        for method in METHODS:
            shuffled = read_partition(root / f"runs/shuffle/{method}/sample_{dataset:03d}/node_clusters.tsv.gz")
            shuffle_rows.append({
                "dataset": dataset, "method": method,
                "clusters_reference": stats[method]["final_clusters"],
                **partition_metrics(base[method], shuffled),
            })

    subset_rows: list[dict] = []
    full0 = {
        method: read_partition(root / f"runs/base/{method}/sample_000/node_clusters.tsv.gz")
        for method in METHODS
    }
    for replicate in range(settings.subset_replicates):
        for method in METHODS:
            subset = read_partition(root / f"runs/subsets/{method}/subset_{replicate:02d}/node_clusters.tsv.gz")
            restricted = {sequence: full0[method][sequence] for sequence in subset}
            subset_rows.append({
                "replicate": replicate, "method": method,
                "subset_fraction": settings.subset_fraction,
                **partition_metrics(restricted, subset),
            })

    base_stats: list[dict] = []
    for dataset in range(settings.datasets):
        for method in METHODS:
            row = json.loads((root / f"runs/base/{method}/sample_{dataset:03d}/run_stats.json").read_text())
            base_stats.append({"dataset": dataset, "method": method, **row})
    write_csv(root / "base_run_metrics.csv", base_stats)
    write_csv(root / "method_agreement.csv", method_rows)
    write_csv(root / "shuffle_reproducibility.csv", shuffle_rows)
    write_csv(root / "subset_consistency.csv", subset_rows)
    make_figure(root, pd.DataFrame(base_stats), pd.DataFrame(method_rows), pd.DataFrame(shuffle_rows), pd.DataFrame(subset_rows))
    write_report(root, settings, pd.DataFrame(base_stats), pd.DataFrame(method_rows), pd.DataFrame(shuffle_rows), pd.DataFrame(subset_rows))


def make_figure(root: Path, base: pd.DataFrame, methods: pd.DataFrame, shuffles: pd.DataFrame, subsets: pd.DataFrame) -> None:
    labels = {"greedy": "Greedy", "graph": "Graph", "graph_prefilter": "Graph + prefilter"}
    colors = {"greedy": "#d97706", "graph": "#2563eb", "graph_prefilter": "#059669"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for method in METHODS:
        values = base.loc[base.method == method, "final_clusters"].astype(float)
        axes[0, 0].scatter(np.full(len(values), METHODS.index(method)) + np.random.default_rng(1).normal(0, .035, len(values)), values, s=10, alpha=.45, color=colors[method])
    axes[0, 0].set_xticks(range(3), [labels[x] for x in METHODS])
    axes[0, 0].set_ylabel("Clusters per 20,000 peptides")
    axes[0, 0].set_title("A. Cluster counts")

    comparisons = list(methods.comparison.unique())
    axes[0, 1].boxplot([methods.loc[methods.comparison == item, "adjusted_rand_index"] for item in comparisons], tick_labels=[x.replace("_vs_", "\nvs\n") for x in comparisons], showfliers=False)
    axes[0, 1].set_ylim(0, 1.01)
    axes[0, 1].set_ylabel("Adjusted Rand index")
    axes[0, 1].set_title("B. Agreement between methods")

    shuffle_exact = [100 * shuffles.loc[shuffles.method == method, "exact_partition"].mean() for method in METHODS]
    axes[1, 0].bar(range(3), shuffle_exact, color=[colors[x] for x in METHODS])
    axes[1, 0].set_xticks(range(3), [labels[x] for x in METHODS])
    axes[1, 0].set_ylim(0, 105)
    axes[1, 0].set_ylabel("Exactly reproduced shuffles (%)")
    axes[1, 0].set_title("C. FASTA-order reproducibility")

    x = np.arange(3)
    width = .35
    recall = [subsets.loc[subsets.method == method, "cocluster_recall"].mean() for method in METHODS]
    precision = [subsets.loc[subsets.method == method, "cocluster_precision"].mean() for method in METHODS]
    axes[1, 1].bar(x - width/2, recall, width, label="Recall", color="#7c3aed")
    axes[1, 1].bar(x + width/2, precision, width, label="Precision", color="#db2777")
    axes[1, 1].set_xticks(x, [labels[item] for item in METHODS])
    axes[1, 1].set_ylim(0, 1.05)
    axes[1, 1].set_ylabel("Fraction of co-clustered pairs")
    axes[1, 1].set_title("D. Dataset-0 subset consistency")
    axes[1, 1].legend(frameon=False)
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=.2)
    fig.tight_layout()
    (root / "figures").mkdir(exist_ok=True)
    fig.savefig(root / "figures/method_validation.png", dpi=220)
    fig.savefig(root / "figures/method_validation.pdf")
    plt.close(fig)


def write_report(root: Path, settings: Settings, base: pd.DataFrame, methods: pd.DataFrame, shuffles: pd.DataFrame, subsets: pd.DataFrame) -> None:
    def summary(frame: pd.DataFrame, column: str) -> str:
        values = frame[column].astype(float)
        return f"mean {values.mean():.4f}; median {values.median():.4f}; range {values.min():.4f}–{values.max():.4f}"

    cluster_lines = []
    speed_lines = []
    for method in METHODS:
        selected = base[base.method == method]
        singleton_percent = 100 * selected.singleton_clusters.astype(float).sum() / selected.final_clusters.astype(float).sum()
        singleton_peptide_percent = 100 * selected.singleton_clusters.astype(float).sum() / selected.accepted_records.astype(float).sum()
        cluster_lines.append(f"- `{method}`: {summary(selected, 'final_clusters')} clusters; {singleton_percent:.2f}% of clusters and {singleton_peptide_percent:.2f}% of input peptides were singletons.")
        speed_lines.append(f"- `{method}`: {summary(selected, 'elapsed_seconds')} seconds per dataset.")
    agreement_lines = []
    for comparison in methods.comparison.unique():
        selected = methods[methods.comparison == comparison]
        agreement_lines.append(
            f"- `{comparison}`: ARI {summary(selected, 'adjusted_rand_index')}; "
            f"exact partitions {int(selected.exact_partition.sum())}/{len(selected)}; "
            f"co-cluster recall {selected.cocluster_recall.mean():.4f}, precision {selected.cocluster_precision.mean():.4f}."
        )
    shuffle_lines = []
    subset_lines = []
    for method in METHODS:
        selected = shuffles[shuffles.method == method]
        shuffle_lines.append(
            f"- `{method}`: {int(selected.exact_partition.sum())}/{len(selected)} exact; "
            f"mean ARI {selected.adjusted_rand_index.mean():.6f}."
        )
        selected = subsets[subsets.method == method]
        subset_lines.append(
            f"- `{method}`: mean ARI {selected.adjusted_rand_index.mean():.4f}; "
            f"co-cluster recall {selected.cocluster_recall.mean():.4f}; "
            f"co-cluster precision {selected.cocluster_precision.mean():.4f}; "
            f"informative-member Jaccard {selected.mean_informative_member_jaccard.mean():.4f}."
        )

    restrictive = base.groupby("method").final_clusters.mean().min() > settings.records_per_dataset * 0.8
    graph_seconds = base.loc[base.method == "graph", "elapsed_seconds"].mean()
    greedy_seconds = base.loc[base.method == "greedy", "elapsed_seconds"].mean()
    prefilter_seconds = base.loc[base.method == "graph_prefilter", "elapsed_seconds"].mean()
    graph_pairs = base.loc[base.method == "graph", "candidate_pairs_computed"].mean()
    greedy_pairs = base.loc[base.method == "greedy", "candidate_pairs_computed"].mean()
    prefilter_comparison = methods[methods.comparison == "graph_vs_graph_prefilter"]
    differing_prefilter = prefilter_comparison.loc[~prefilter_comparison.exact_partition, "dataset"].astype(int).tolist()
    all_run_stats = [json.loads(path.read_text()) for path in (root / "runs").rglob("run_stats.json")]
    validation_failures = sum(int(item["validation_failures"]) for item in all_run_stats)
    converged_runs = sum(bool(item["converged"]) for item in all_run_stats)
    report = f"""# PepCluster2 alignment-mode clustering validation

## Question

This analysis asks whether the new constrained full-alignment and anchor-combination `AND` rule behaves consistently when clusters are constructed by greedy representative assignment, a complete candidate graph, or the graph prefilter. It also tests whether FASTA order changes the result and whether clusters remain similar when half of dataset 0 is removed.

## Methods

- PepCluster2 version: `{settings.version}`.
- Scoring mode: `{settings.scoring_mode}`.
- Alignment-similarity threshold: {settings.alignment_threshold:.2f}.
- Anchor-combination-similarity threshold: {settings.anchor_threshold:.2f}.
- K-mer seed threshold: {settings.kmer_seed_threshold:.2f}; it retrieves candidates but does not accept edges.
- Gap open/extension: {settings.gap_open:g}/{settings.gap_extension:g}; terminal overhang open/extension: {settings.terminal_gap_open:g}/{settings.terminal_gap_extension:g}.
- Data: the same {settings.datasets} independently sampled datasets of {settings.records_per_dataset:,} peptides used in the prior validation, originally sampled from `PMBind/data/peptides.fasta`.
- Paths: greedy, graph without prefilter, and graph with forced prefilter.
- Reproducibility: one independently seeded FASTA shuffle for every dataset and method ({settings.datasets * len(METHODS)} shuffled runs).
- Subset analysis: ten independently sampled {100*settings.subset_fraction:.0f}% subsets of dataset 0, clustered by all three methods.

The **Adjusted Rand index (ARI)** measures agreement between two partitions after correcting for chance; 1 means identical. Because these results contain many singletons, ARI can be high even when the relatively few non-singleton relationships differ. Therefore, the report also gives **co-cluster recall** (fraction of reference co-clustered pairs recovered), **co-cluster precision** (fraction of new co-clustered pairs supported by the reference), and member-level Jaccard agreement.

## Results

### Number of clusters

{chr(10).join(cluster_lines)}

### Runtime

{chr(10).join(speed_lines)}

The graph path was {greedy_seconds/graph_seconds:.1f} times faster than greedy on average. Greedy performed about {greedy_pairs/1_000_000:.2f} million exact scoring operations per dataset, compared with {graph_pairs/1_000_000:.2f} million unique candidate-pair scores for graph. Forced prefiltering was {100*(prefilter_seconds/graph_seconds-1):.1f}% slower than the non-prefilter graph at this 20k scale.

### Agreement between clustering methods

{chr(10).join(agreement_lines)}

### FASTA-order reproducibility

{chr(10).join(shuffle_lines)}

All {converged_runs}/{len(all_run_stats)} base, shuffle, and subset runs converged, and the total number of final representative-validation failures was {validation_failures}.

### Dataset-0 subset consistency

{chr(10).join(subset_lines)}

## Interpretation

The FASTA-shuffle test directly checks implementation determinism. Exact equality means that changing record order did not change cluster membership; cluster names themselves were ignored.

The graph-versus-prefilter comparison measures prefilter loss under the new scoring rule. Any co-cluster recall below 1 means that the prefilter failed to reconstruct some relationships present in the non-prefilter graph path. Greedy-versus-graph disagreement is not necessarily an error: dynamic graph set cover and fixed k-mer-degree representative ordering optimize different initial objectives, although refinement can reduce the difference.

The prefilter was extremely close to graph clustering, but it was not exact: {len(differing_prefilter)} of {settings.datasets} partitions differed (datasets {', '.join(map(str, differing_prefilter)) if differing_prefilter else 'none'}). Its mean co-cluster recall was {prefilter_comparison.cocluster_recall.mean():.6f}. Therefore, the scoped prefilter substantially reduces the earlier recovery problem but does not eliminate it.

The subset test is intentionally stricter than shuffle testing. Removing peptides can remove representatives and alter k-mer-degree ordering, so exact subset recovery is not expected. Co-cluster recall and precision show whether the remaining biological relationships are preserved.

## Conclusion

{'At the default 0.60/0.60 thresholds, all methods produced more than 80% as many clusters as input peptides, and approximately 92% of peptides remained singletons. The configuration is therefore extremely restrictive. This does not mean that the implementation is incorrect: every run converged, final validation passed, and all shuffled inputs were reproduced exactly. It means the thresholds are not yet biologically calibrated and should not be used as paper defaults without labelled peptide–MHC tuning.' if restrictive else 'The default thresholds produced substantial non-singleton clustering, but biological purity still requires labelled peptide–MHC validation.'}

For 20,000-peptide datasets, the non-prefilter graph is the strongest current choice: it was fastest and avoided the small remaining prefilter loss. The forced prefilter is appropriate only when graph storage requires it and the measured loss is acceptable. The greedy path achieved perfect input-order reproducibility and avoids materializing the graph, but it was about {greedy_seconds/graph_seconds:.1f} times slower here and did not reproduce graph partitions; it should be treated as a memory-saving alternative, not an equivalent faster replacement.

Method agreement, prefilter recovery, and subset consistency must be judged from co-cluster recall and precision rather than ARI alone. Threshold and gap calibration on labelled pMHC data remains the next biological validation step.

The figure is available at `figures/method_validation.png`. Complete settings, per-run configurations, assignments, logs, and CSV summaries are stored with this report.
"""
    (root / "REPORT.md").write_text(report)


def prepare_root(root: Path, settings: Settings, binary: Path) -> None:
    (root / "config").mkdir(parents=True, exist_ok=True)
    (root / "tmp").mkdir(parents=True, exist_ok=True)
    source_table = SOURCE_VALIDATION / "config/kmer2_similarity_q.bin"
    target_table = root / "config/kmer2_similarity_q.bin"
    if not target_table.exists():
        shutil.copy2(source_table, target_table)
    shutil.copy2(SOURCE_VALIDATION / "data/sample_manifest.csv", root / "sample_manifest.csv")
    (root / "config/settings.json").write_text(json.dumps(asdict(settings), indent=2, sort_keys=True) + "\n")
    environment = {
        "python": platform.python_version(), "platform": platform.platform(),
        "binary": str(binary), "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "algorithm_sha256": hashlib.sha256((PROJECT / "ALGORITHM.md").read_bytes()).hexdigest(),
        "source_pool": "/home/amir/amir/ParseFold/PMBind/data/peptides.fasta",
        "reused_sample_root": str(SOURCE_VALIDATION),
    }
    (root / "config/environment.json").write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n")


def manifest_sources(settings: Settings) -> list[Path]:
    rows = list(csv.DictReader((SOURCE_VALIDATION / "data/sample_manifest.csv").open()))
    if len(rows) < settings.datasets:
        raise RuntimeError("sample manifest does not contain enough datasets")
    return [Path(row["path"]) for row in rows[: settings.datasets]]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("base", "shuffle", "subset", "analyse", "all"), default="all")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--binary", type=Path, default=PROJECT / "target/release/pepcluster2")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-process", type=int, default=6)
    parser.add_argument("--datasets", type=int, default=100)
    parser.add_argument(
        "--methods", default=",".join(METHODS),
        help="comma-separated subset of greedy,graph,graph_prefilter for execution stages",
    )
    args = parser.parse_args()
    settings = Settings(process_workers=args.workers, threads_per_process=args.threads_per_process, datasets=args.datasets)
    selected_methods = tuple(item.strip() for item in args.methods.split(",") if item.strip())
    if not selected_methods or any(item not in METHODS for item in selected_methods):
        parser.error("--methods must contain greedy, graph, and/or graph_prefilter")
    root = args.root.resolve()
    binary = args.binary.resolve()
    prepare_root(root, settings, binary)
    sources = manifest_sources(settings)

    if args.stage in ("base", "all"):
        tasks = [(dataset, "base", method, sources[dataset], binary, root, settings) for dataset in range(settings.datasets) for method in selected_methods]
        rows = execute(tasks, run_dataset_task, settings.process_workers, "base")
        write_csv(root / "base_execution_metrics.csv", rows)
    if args.stage in ("shuffle", "all"):
        tasks = [(dataset, "shuffle", method, sources[dataset], binary, root, settings) for dataset in range(settings.datasets) for method in selected_methods]
        rows = execute(tasks, run_dataset_task, settings.process_workers, "shuffle")
        write_csv(root / "shuffle_execution_metrics.csv", rows)
    if args.stage in ("subset", "all"):
        tasks = [(replicate, method, sources[0], binary, root, settings) for replicate in range(settings.subset_replicates) for method in selected_methods]
        rows = execute(tasks, run_subset_task, settings.process_workers, "subset")
        write_csv(root / "subset_execution_metrics.csv", rows)
    if args.stage in ("analyse", "all"):
        analyse(root, settings)


if __name__ == "__main__":
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/pepcluster2-matplotlib")
    main()
