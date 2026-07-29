#!/usr/bin/env python3
"""Reproducible PepCluster2 index, clustering, and shuffle validation."""

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
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    rand_score,
)


ROOT = Path("/home/amir/amir/ParseFold/Pepcluster2")
BINARY = ROOT / "target/release/pepcluster2"
POOL = Path("/home/amir/amir/ParseFold/PMBind/data/peptides.fasta")
SOURCE_DATA = Path(
    "/home/amir/amir/ParseFold/PepCluster/analysis/results/"
    "00_exhaustive_search/2026-07-27_0.1.7/data"
)
DEFAULT_RESULT = ROOT / "validation/2026-07-28_0.2.0-dev_mmseqs_defaults"
CANONICAL = frozenset("ARNDCQEGHILKMFPSTWYV")


@dataclass(frozen=True)
class Settings:
    master_seed: int = 20260727
    n_datasets: int = 100
    sample_size: int = 20_000
    n_shuffles: int = 10
    kmer_seed_threshold: float = 0.5
    edge_threshold: float = 0.6
    iterations: int = 3
    minimum_improvement: float = 0.01
    strict_merge: bool = True
    workers: int = 4
    threads_per_run: int = 6


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(value)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def read_fasta_gz(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
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
                sequence.append(line.upper())
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]], order: Iterable[int] | None = None) -> None:
    indices = range(len(records)) if order is None else order
    with path.open("w", encoding="utf-8") as handle:
        for index in indices:
            header, sequence = records[int(index)]
            handle.write(f">{header}\n{sequence}\n")


def geometry_mask(length: int) -> int:
    combinations = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    result = 0
    for bit, (front, back_local) in enumerate(combinations):
        if length - 3 + back_local >= front + 6:
            result |= 1 << bit
    return result


def node_key(sequence: str) -> tuple[str, int]:
    if len(sequence) < 7 or not set(sequence).issubset(CANONICAL):
        raise ValueError(f"invalid validation peptide: {sequence!r}")
    return sequence[:3] + sequence[-3:], geometry_mask(len(sequence))


def seed_for(master_seed: int, stage: int, dataset_id: int, replicate: int = 0) -> int:
    state = np.random.SeedSequence(
        [master_seed, stage, dataset_id, replicate]
    ).generate_state(1, dtype=np.uint64)
    return int(state[0])


def run_binary(
    input_fasta: Path,
    output_dir: Path,
    temporary_dir: Path,
    table: Path,
    settings: Settings,
    index_only: bool = False,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(BINARY),
        "--input",
        str(input_fasta),
        "--output-dir",
        str(output_dir),
        "--tmp-dir",
        str(temporary_dir),
        "--kmer-table",
        str(table),
        "--kmer-seed-threshold",
        str(settings.kmer_seed_threshold),
        "--threshold",
        str(settings.edge_threshold),
        "--iterations",
        str(settings.iterations),
        "--min-improvement",
        str(settings.minimum_improvement),
        "--threads",
        str(settings.threads_per_run),
        "--compact-output",
    ]
    if index_only:
        command.append("--index-only")
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(
            f"PepCluster2 failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )


def read_node_mapping(path: Path) -> dict[tuple[str, int], str]:
    result: dict[tuple[str, int], str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["anchor"], int(row["geometry_mask"]))] = row["cluster_id"]
    return result


def save_mapping(path: Path, mapping_file: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(fd)
    try:
        with mapping_file.open("rb") as source, gzip.open(temporary, "wb", compresslevel=6) as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def load_saved_mapping(path: Path) -> dict[tuple[str, int], str]:
    result: dict[tuple[str, int], str] = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["anchor"], int(row["geometry_mask"]))] = row["cluster_id"]
    return result


def cluster_raw_path(result_dir: Path, dataset_id: int) -> Path:
    return result_dir / "raw/clusters" / f"sample_{dataset_id:03d}.json"


def mapping_path(result_dir: Path, dataset_id: int) -> Path:
    return result_dir / "raw/reference_assignments" / f"sample_{dataset_id:03d}.tsv.gz"


def stability_raw_path(result_dir: Path, dataset_id: int) -> Path:
    return result_dir / "raw/stability" / f"sample_{dataset_id:03d}.csv"


def cluster_dataset(dataset_id: int, result_dir_text: str, settings_dict: dict) -> str:
    settings = Settings(**settings_dict)
    result_dir = Path(result_dir_text)
    output = cluster_raw_path(result_dir, dataset_id)
    saved_mapping = mapping_path(result_dir, dataset_id)
    if output.exists() and saved_mapping.exists():
        return str(output)
    sample = result_dir / "data" / f"sample_{dataset_id:03d}.fasta.gz"
    records = read_fasta_gz(sample)
    with tempfile.TemporaryDirectory(prefix=f"pc2-cluster-{dataset_id:03d}-") as tmp_text:
        tmp = Path(tmp_text)
        fasta = tmp / "sample.fasta"
        run_dir = tmp / "run"
        write_fasta(fasta, records)
        started = time.perf_counter()
        run_binary(
            fasta,
            run_dir,
            tmp / "work",
            result_dir / "config/kmer2_similarity_q.bin",
            settings,
        )
        elapsed = time.perf_counter() - started
        with (run_dir / "run_stats.json").open(encoding="utf-8") as handle:
            row = json.load(handle)
        row.update(
            {
                "dataset_id": dataset_id,
                "sampling_seed": seed_for(settings.master_seed, 0, dataset_id),
                "wall_seconds_from_python": elapsed,
            }
        )
        save_mapping(saved_mapping, run_dir / "node_clusters.tsv")
    atomic_json(output, row)
    return str(output)


def canonical_labels(labels: list[str]) -> np.ndarray:
    seen: dict[str, int] = {}
    result = np.empty(len(labels), dtype=np.int32)
    next_label = 0
    for index, label in enumerate(labels):
        if label not in seen:
            seen[label] = next_label
            next_label += 1
        result[index] = seen[label]
    return result


def variation_of_information(a: np.ndarray, b: np.ndarray) -> float:
    _, a = np.unique(a, return_inverse=True)
    _, b = np.unique(b, return_inverse=True)
    contingency = np.zeros((int(a.max()) + 1, int(b.max()) + 1), dtype=np.int64)
    np.add.at(contingency, (a, b), 1)
    n = contingency.sum()
    pa = contingency.sum(axis=1) / n
    pb = contingency.sum(axis=0) / n
    rows, cols = np.nonzero(contingency)
    joint = contingency[rows, cols] / n
    mutual_information = np.sum(joint * np.log(joint / (pa[rows] * pb[cols])))
    h_a = -np.sum(pa[pa > 0] * np.log(pa[pa > 0]))
    h_b = -np.sum(pb[pb > 0] * np.log(pb[pb > 0]))
    return float(h_a + h_b - 2 * mutual_information)


def changed_fraction(a: np.ndarray, b: np.ndarray) -> float:
    _, a = np.unique(a, return_inverse=True)
    _, b = np.unique(b, return_inverse=True)
    contingency = np.zeros((int(a.max()) + 1, int(b.max()) + 1), dtype=np.int64)
    np.add.at(contingency, (a, b), 1)
    rows, cols = linear_sum_assignment(contingency, maximize=True)
    return 1.0 - int(contingency[rows, cols].sum()) / len(a)


def stability_dataset(dataset_id: int, result_dir_text: str, settings_dict: dict) -> str:
    settings = Settings(**settings_dict)
    result_dir = Path(result_dir_text)
    output = stability_raw_path(result_dir, dataset_id)
    if output.exists():
        return str(output)
    sample = result_dir / "data" / f"sample_{dataset_id:03d}.fasta.gz"
    records = read_fasta_gz(sample)
    keys = [node_key(sequence) for _, sequence in records]
    reference_mapping = load_saved_mapping(mapping_path(result_dir, dataset_id))
    reference_raw = [reference_mapping[key] for key in keys]
    reference_labels = canonical_labels(reference_raw)
    reference_clusters = len(set(reference_raw))
    rows: list[dict] = []
    with tempfile.TemporaryDirectory(prefix=f"pc2-stability-{dataset_id:03d}-") as tmp_text:
        tmp = Path(tmp_text)
        fasta = tmp / "shuffle.fasta"
        for shuffle_id in range(settings.n_shuffles):
            shuffle_seed = seed_for(settings.master_seed, 1, dataset_id, shuffle_id)
            order = np.random.default_rng(shuffle_seed).permutation(len(records))
            write_fasta(fasta, records, order)
            run_dir = tmp / f"run_{shuffle_id:02d}"
            started = time.perf_counter()
            run_binary(
                fasta,
                run_dir,
                tmp / f"work_{shuffle_id:02d}",
                result_dir / "config/kmer2_similarity_q.bin",
                settings,
            )
            elapsed = time.perf_counter() - started
            observed_mapping = read_node_mapping(run_dir / "node_clusters.tsv")
            observed_raw = [observed_mapping[key] for key in keys]
            observed_labels = canonical_labels(observed_raw)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "shuffle_id": shuffle_id,
                    "shuffle_seed": shuffle_seed,
                    "n_peptides": len(records),
                    "reference_n_clusters": reference_clusters,
                    "shuffled_n_clusters": len(set(observed_raw)),
                    "cluster_count_delta": len(set(observed_raw)) - reference_clusters,
                    "partition_identical": bool(np.array_equal(reference_labels, observed_labels)),
                    "adjusted_rand_index": adjusted_rand_score(reference_labels, observed_labels),
                    "normalized_mutual_information": normalized_mutual_info_score(
                        reference_labels, observed_labels
                    ),
                    "rand_index": rand_score(reference_labels, observed_labels),
                    "variation_of_information_nats": variation_of_information(
                        reference_labels, observed_labels
                    ),
                    "optimal_label_changed_fraction": changed_fraction(
                        reference_labels, observed_labels
                    ),
                    "representative_label_changed_fraction": float(
                        np.mean(np.asarray(reference_raw, dtype=object) != np.asarray(observed_raw, dtype=object))
                    ),
                    "shuffle_elapsed_seconds": elapsed,
                }
            )
            shutil.rmtree(run_dir)
    atomic_csv(output, rows)
    return str(output)


def run_parallel(name: str, worker, result_dir: Path, settings: Settings) -> None:
    settings_dict = asdict(settings)
    pending = []
    for dataset_id in range(settings.n_datasets):
        expected = (
            cluster_raw_path(result_dir, dataset_id)
            if name == "clusters"
            else stability_raw_path(result_dir, dataset_id)
        )
        if not expected.exists() or (
            name == "clusters" and not mapping_path(result_dir, dataset_id).exists()
        ):
            pending.append(dataset_id)
    if not pending:
        print(f"[{name}] already complete", flush=True)
        return
    print(f"[{name}] {len(pending)} datasets with {settings.workers} workers", flush=True)
    with ProcessPoolExecutor(max_workers=settings.workers) as executor:
        futures = {
            executor.submit(worker, dataset_id, str(result_dir), settings_dict): dataset_id
            for dataset_id in pending
        }
        completed = settings.n_datasets - len(pending)
        for future in as_completed(futures):
            dataset_id = futures[future]
            future.result()
            completed += 1
            print(f"[{name}] {completed}/{settings.n_datasets} (sample {dataset_id:03d})", flush=True)


def initialize(result_dir: Path, settings: Settings) -> None:
    for relative in (
        "config",
        "data",
        "index",
        "raw/clusters",
        "raw/reference_assignments",
        "raw/stability",
        "figures",
    ):
        (result_dir / relative).mkdir(parents=True, exist_ok=True)
    manifest = []
    for dataset_id in range(settings.n_datasets):
        source = SOURCE_DATA / f"sample_{dataset_id:03d}.fasta.gz"
        destination = result_dir / "data" / source.name
        if not destination.exists():
            try:
                os.link(source, destination)
            except OSError:
                shutil.copy2(source, destination)
        manifest.append(
            {
                "dataset_id": dataset_id,
                "sampling_seed": seed_for(settings.master_seed, 0, dataset_id),
                "records": settings.sample_size,
                "sha256": sha256(destination),
                "path": str(destination),
            }
        )
    atomic_csv(result_dir / "data/sample_manifest.csv", manifest)
    config = {
        "created_utc": utc_now(),
        "analysis": "PepCluster2 MMseqs-style candidate and shuffle validation",
        "pepcluster2_version": subprocess.check_output([str(BINARY), "--version"], text=True).strip(),
        "binary": str(BINARY),
        "pool": str(POOL),
        "pool_sha256": sha256(POOL),
        "source_pepcluster_analysis": str(SOURCE_DATA.parent),
        "settings": asdict(settings),
        "definitions": {
            "comparison_fraction": "unique candidate node pairs scored / all unordered distinct node pairs",
            "stability_reference": "saved sample order",
            "stability_replicates": "10 deterministic shuffles per each of 100 datasets",
            "candidate_rule": "at least one front and one back 2-mer neighbour at normalized BLOSUM score >= 0.5",
            "edge_rule": "optimal one-to-one valid anchor-combination score >= 0.6",
        },
    }
    existing = result_dir / "config/config.json"
    if existing.exists():
        with existing.open(encoding="utf-8") as handle:
            old = json.load(handle)
        if old["settings"] != config["settings"] or old["pool_sha256"] != config["pool_sha256"]:
            raise RuntimeError(f"existing result configuration differs: {existing}")
    else:
        atomic_json(existing, config)
    environment = (
        f"created_utc={utc_now()}\n"
        f"python={platform.python_version()}\n"
        f"platform={platform.platform()}\n"
        f"numpy={np.__version__}\n"
        f"scipy={__import__('scipy').__version__}\n"
        f"scikit_learn={__import__('sklearn').__version__}\n"
        f"pepcluster2={config['pepcluster2_version']}\n"
    )
    atomic_text(result_dir / "config/environment.txt", environment)
    tracked_sources = [
        ROOT / "Cargo.toml",
        ROOT / "Cargo.lock",
        ROOT / "src/cli.rs",
        ROOT / "src/edge_store.rs",
        ROOT / "src/fasta.rs",
        ROOT / "src/graph.rs",
        ROOT / "src/index.rs",
        ROOT / "src/kmer.rs",
        ROOT / "src/main.rs",
        ROOT / "src/model.rs",
        ROOT / "src/output.rs",
        ROOT / "src/scoring.rs",
        ROOT / "validation/src/run_validation.py",
    ]
    atomic_json(
        result_dir / "config/implementation_checksums.json",
        {
            "release_binary": {"path": str(BINARY), "sha256": sha256(BINARY)},
            "sources": {str(path.relative_to(ROOT)): sha256(path) for path in tracked_sources},
        },
    )
    atomic_text(
        result_dir / "config/last_command.txt",
        " ".join([sys.executable, *sys.argv]) + "\n",
    )


def ensure_kmer_table(result_dir: Path, settings: Settings) -> None:
    table = result_dir / "config/kmer2_similarity_q.bin"
    if not table.exists():
        with tempfile.TemporaryDirectory(prefix="pc2-table-") as tmp_text:
            tmp = Path(tmp_text)
            fasta = tmp / "tiny.fasta"
            write_fasta(fasta, [("a", "FLNVIVHKA"), ("b", "FLNVIVHSA")])
            run_binary(fasta, tmp / "out", tmp / "work", table, settings, index_only=True)
    atomic_json(
        result_dir / "config/kmer_table_metadata.json",
        {
            "path": str(table),
            "format": "PC2K2S01: 16-byte header followed by 400x400 little-endian i16 scores",
            "bytes": table.stat().st_size,
            "sha256": sha256(table),
            "loading": "read-only memory map on Unix; one-time buffered read fallback elsewhere",
        },
    )


def run_index(result_dir: Path, settings: Settings) -> None:
    output = result_dir / "index/full_pool"
    if (output / "run_stats.json").exists():
        print("[index] already complete", flush=True)
        return
    print("[index] full peptide pool", flush=True)
    run_binary(
        POOL,
        output,
        result_dir / "index/tmp",
        result_dir / "config/kmer2_similarity_q.bin",
        settings,
        index_only=True,
    )


def aggregate(result_dir: Path, settings: Settings) -> None:
    cluster_rows = []
    for dataset_id in range(settings.n_datasets):
        with cluster_raw_path(result_dir, dataset_id).open(encoding="utf-8") as handle:
            cluster_rows.append(json.load(handle))
    atomic_csv(result_dir / "cluster_metrics.csv", cluster_rows)
    stability_rows: list[dict] = []
    for dataset_id in range(settings.n_datasets):
        with stability_raw_path(result_dir, dataset_id).open(newline="", encoding="utf-8") as handle:
            stability_rows.extend(csv.DictReader(handle))
    atomic_csv(result_dir / "stability_metrics.csv", stability_rows)


def make_report(result_dir: Path) -> None:
    import pandas as pd

    os.environ.setdefault("MPLCONFIGDIR", "/tmp/pepcluster2-matplotlib")
    import matplotlib.pyplot as plt

    clusters = pd.read_csv(result_dir / "cluster_metrics.csv")
    stability = pd.read_csv(result_dir / "stability_metrics.csv")
    with (result_dir / "index/full_pool/run_stats.json").open(encoding="utf-8") as handle:
        index = json.load(handle)

    fraction = clusters["fraction_all_pairs_computed"]
    cluster_count = clusters["final_clusters"]
    exact = stability["partition_identical"].astype(str).str.lower().eq("true")
    ari = stability["adjusted_rand_index"].astype(float)
    changed = stability["optimal_label_changed_fraction"].astype(float)

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4))
    axes[0].hist(100 * fraction, bins=20, color="#2874a6", edgecolor="white")
    axes[0].axvline(100 * fraction.mean(), color="#922b21", linestyle="--")
    axes[0].set_xlabel("All unordered node pairs scored (%)")
    axes[0].set_ylabel("Datasets")
    axes[0].set_title("Candidate-search workload")
    axes[1].hist(cluster_count, bins=20, color="#239b56", edgecolor="white")
    axes[1].axvline(cluster_count.mean(), color="#922b21", linestyle="--")
    axes[1].set_xlabel("Final clusters in 20,000 peptides")
    axes[1].set_ylabel("Datasets")
    axes[1].set_title("Cluster count")
    if np.allclose(ari, 1.0):
        axes[2].bar([0], [100 * exact.mean()], color="#7d3c98", width=0.55)
        axes[2].set_xticks([0], ["1,000 shuffles"])
        axes[2].set_ylim(0, 110)
        axes[2].set_ylabel("Exactly identical partitions (%)")
    else:
        axes[2].hist(ari, bins=20, color="#7d3c98", edgecolor="white")
        axes[2].set_xlabel("Adjusted Rand Index (1 = identical)")
        axes[2].set_ylabel("Shuffled runs")
    axes[2].set_title("Input-order stability")
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle("PepCluster2 default-setting validation", fontsize=14)
    figure.tight_layout()
    figure.savefig(result_dir / "figures/validation_summary.png", dpi=190, bbox_inches="tight")
    figure.savefig(result_dir / "figures/validation_summary.pdf", bbox_inches="tight")
    plt.close(figure)

    possible_full = index["unique_anchors"] * (index["unique_anchors"] - 1) // 2
    exact_count = int(exact.sum())
    text = f"""# PepCluster2 default-setting validation

## Purpose

This analysis tested the new similarity-neighbour candidate search and six-anchor edge score. It used the same 100 random datasets of 20,000 peptide records used for PepCluster 0.1.7. Each dataset was also clustered after 10 independent input shuffles, producing 1,000 shuffled runs.

## Default settings

- Candidate seed: normalized BLOSUM similarity of at least **0.50** between terminal 2-mers.
- Final graph edge: optimal one-to-one anchor-combination similarity of at least **0.60**.
- Six geometrically valid anchor hypotheses were used for peptides of length nine or more; shorter peptides retained only hypotheses whose residues were at least six positions apart.
- Candidate pairs were written in bounded chunks, externally sorted and deduplicated, and only then scored.
- Clustering used deterministic greedy set cover, three reassignment/refinement iterations and strict merging.

## Results

### Candidate-search workload

Across the 100 datasets, PepCluster2 scored a mean of **{100*fraction.mean():.3f}%** of all unordered pairs of distinct terminal-anchor nodes. The range was **{100*fraction.min():.3f}%–{100*fraction.max():.3f}%**. Thus, the similarity-neighbour search avoided approximately **{100*(1-fraction.mean()):.3f}%** of exhaustive comparisons while recovering conservative k-mer substitutions.

The mean number of unique candidate pairs actually scored was **{clusters['unique_candidate_pairs_computed'].mean():,.0f}** per dataset. External deduplication reduced a mean of **{clusters['candidate_occurrences_before_dedup'].mean():,.0f}** candidate occurrences to those unique pairs before anchor scoring.

The complete 20,000-peptide run took a mean of **{clusters['elapsed_seconds'].mean():.3f} seconds** in the recorded four-process validation workload. Pair deduplication removed **{100*(1-clusters['unique_candidate_pairs_computed'].sum()/clusters['candidate_occurrences_before_dedup'].sum()):.2f}%** of candidate occurrences before the more expensive anchor score.

### Number and size of clusters

The mean final cluster count was **{cluster_count.mean():,.1f}** per 20,000-peptide dataset, with a range of **{cluster_count.min():,}–{cluster_count.max():,}**. The median was **{cluster_count.median():,.0f}**. The mean largest cluster contained **{clusters['largest_cluster_peptides'].mean():.1f}** peptide records.

These cluster counts describe algorithmic behaviour at the provisional thresholds. They do not by themselves demonstrate that the clusters correspond to HLA binding specificity; labelled biological validation remains necessary.

### Stability under FASTA shuffling

**{exact_count} of {len(stability)} shuffled runs ({100*exact.mean():.1f}%) produced exactly the same peptide partition as the unshuffled reference.** The mean Adjusted Rand Index was **{ari.mean():.6f}** and the mean fraction of peptide records changing optimally matched clusters was **{100*changed.mean():.6f}%**.

The Adjusted Rand Index (ARI) compares every pair of peptide records and asks whether the two clusterings agree that they belong together or apart, while correcting for agreement expected by chance. ARI = 1 means identical partitions.

### Full-pool index sizing

The full pool contained **{index['input_records']:,}** records and **{index['unique_anchors']:,}** distinct terminal-anchor/geometry nodes. Exhaustive comparison would require approximately **{possible_full:,}** unordered node pairs. The neighbour index contained **{index['similar_composite_key_relations']:,}** occupied similar-key relations. Its pre-deduplication candidate-occurrence upper bound was **{index['candidate_occurrence_upper_bound']:,}**.

This upper bound includes duplicates and possible self-occurrences across related buckets; it is a sizing statistic, not the number that would actually be scored. A complete full-pool clustering should be attempted only with adequate temporary disk and after partition-level sizing.

If all upper-bound occurrences were materialized as eight-byte node-ID pairs, they would require approximately **{8*index['candidate_occurrence_upper_bound']/1e12:.2f} TB** before chunk-merge overhead. Therefore, disk spilling controls memory but does not yet make the full-pool expansion practical.

## Interpretation and conclusion

The implementation achieves its computational objective on 20,000-peptide datasets: conservative terminal substitutions are included, the expensive anchor score is applied only after external deduplication, and only about **{100*fraction.mean():.2f}%** of possible node pairs are scored. The observed runtime and temporary workload are practical at this scale.

The input-order objective is **{'fully achieved in this experiment' if exact_count == len(stability) else 'not fully achieved'}**. {'All 1,000 shuffles returned exactly the same clusters.' if exact_count == len(stability) else 'Some shuffled inputs changed the resulting partition and require investigation.'}

The much denser graph produced relatively large clusters: the mean largest cluster contained about **{clusters['largest_cluster_peptides'].mean():.0f} peptides**. This may be useful motif recovery, but it may also indicate over-clustering at the provisional 0.60 edge threshold.

The present analysis therefore does **not** establish biological accuracy. The next decisive test should measure whether known binder/binder pairs receive higher scores than binder/non-binder pairs across multiple HLA alleles, and use those labels to calibrate the seed and edge thresholds. Until that test, version 0.2.0-dev should remain a validation build rather than the paper's final clustering method.

![Validation summary](figures/validation_summary.png)
"""
    atomic_text(result_dir / "REPORT.md", text)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT)
    parser.add_argument(
        "--stage",
        choices=("all", "index", "clusters", "stability", "report"),
        default="all",
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads-per-run", type=int, default=6)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    settings = Settings(workers=args.workers, threads_per_run=args.threads_per_run)
    if not BINARY.exists():
        raise FileNotFoundError(f"build the release binary first: {BINARY}")
    initialize(args.result_dir, settings)
    ensure_kmer_table(args.result_dir, settings)
    if args.stage in ("all", "index"):
        run_index(args.result_dir, settings)
    if args.stage in ("all", "clusters"):
        run_parallel("clusters", cluster_dataset, args.result_dir, settings)
    if args.stage in ("all", "stability"):
        missing = [
            dataset_id
            for dataset_id in range(settings.n_datasets)
            if not mapping_path(args.result_dir, dataset_id).exists()
        ]
        if missing:
            raise RuntimeError("reference assignments are incomplete; run --stage clusters first")
        run_parallel("stability", stability_dataset, args.result_dir, settings)
    cluster_complete = all(
        cluster_raw_path(args.result_dir, dataset_id).exists()
        for dataset_id in range(settings.n_datasets)
    )
    stability_complete = all(
        stability_raw_path(args.result_dir, dataset_id).exists()
        for dataset_id in range(settings.n_datasets)
    )
    index_complete = (args.result_dir / "index/full_pool/run_stats.json").exists()
    if cluster_complete and stability_complete:
        aggregate(args.result_dir, settings)
    if args.stage in ("all", "report"):
        if not (cluster_complete and stability_complete and index_complete):
            raise RuntimeError("index, cluster and stability stages must complete before report")
        make_report(args.result_dir)
    atomic_json(
        args.result_dir / "run_status.json",
        {
            "updated_utc": utc_now(),
            "index_complete": index_complete,
            "cluster_datasets_complete": sum(
                cluster_raw_path(args.result_dir, dataset_id).exists()
                for dataset_id in range(settings.n_datasets)
            ),
            "stability_datasets_complete": sum(
                stability_raw_path(args.result_dir, dataset_id).exists()
                for dataset_id in range(settings.n_datasets)
            ),
            "shuffled_runs_complete": sum(
                settings.n_shuffles
                for dataset_id in range(settings.n_datasets)
                if stability_raw_path(args.result_dir, dataset_id).exists()
            ),
        },
    )
    print(f"results: {args.result_dir}", flush=True)


if __name__ == "__main__":
    main()
