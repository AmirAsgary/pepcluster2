#!/usr/bin/env python3
"""Sweep the reassignment hysteresis margin.

Reassignment, not merging, is where subset stability is lost: on the exhaustive
edge set the 80% pairwise Jaccard falls from 0.5145 to 0.4484 (coverage order)
and from 0.7636 to 0.5709 (intrinsic order) when reassignment runs, while
merging is worth +0.003 and -0.035 respectively.

The cause is the argmax: a peptide takes whichever representative scores best,
so a near-tie flips when the dataset composition changes. A margin makes the
choice sticky. This measures what that buys and what it costs.

Stability here is self-referential — a subset run against the same
configuration's full run — so no reference partition is needed and each margin
can be evaluated independently.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A  # noqa: E402
import run_validation as V  # noqa: E402

MARGINS = ("0.00", "0.005", "0.01", "0.02", "0.05")
SUBSET_SIZES = (1_000, 4_000, 8_000)


def run_one(root: Path, binary: Path, margin: str, order: str, dataset: int,
            threads: int, size: int | None) -> str:
    tag = f"margin_{margin}"
    if size is None:
        source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
        output = root / "runs" / "margin" / tag / order / "full" / f"sample_{dataset:03d}"
    else:
        source = (root / "data" / "subsets" / f"n_{size:06d}"
                  / f"sample_{dataset:03d}.fasta.gz")
        output = (root / "runs" / "margin" / tag / order / f"n_{size:06d}"
                  / f"sample_{dataset:03d}")
    label = f"{tag} {order} {size or 'full'} {dataset:03d}"
    if (output / "node_clusters.tsv").exists() and (output / "run_stats.json").exists():
        return f"{label} cached"
    import tempfile
    with tempfile.TemporaryDirectory(prefix="pc2margin_", dir="/tmp") as tmp:
        fasta = Path(tmp) / "in.fasta"
        V.decompress(source, fasta)
        command = V.cluster_command(binary, fasta, output, threads, V.DEFAULT_GEOMETRY,
                                    V.DEFAULT_SEED_THRESHOLD, order)
        command += ["--reassignment-margin", margin,
                    "--tmp-dir", str(Path(tmp) / "tmp")]
        command += V.METHODS["graph"]
        V.execute(command, output, output / "run.log", output / "resource.txt")
    return f"{label} complete"


def analyse(root: Path, datasets: int) -> None:
    rows = []
    for margin in MARGINS:
        tag = f"margin_{margin}"
        for order in V.ORDERS:
            for dataset in range(datasets):
                full_path = (root / "runs" / "margin" / tag / order / "full"
                             / f"sample_{dataset:03d}")
                if not (full_path / "node_clusters.tsv").exists():
                    continue
                full = A.read_partition(full_path / "node_clusters.tsv")
                stats = A.read_stats(full_path / "run_stats.json")
                for size in SUBSET_SIZES:
                    subset_path = (root / "runs" / "margin" / tag / order
                                   / f"n_{size:06d}" / f"sample_{dataset:03d}")
                    if not (subset_path / "node_clusters.tsv").exists():
                        continue
                    subset = A.read_partition(subset_path / "node_clusters.tsv")
                    restricted = {s: full[s] for s in subset}
                    rows.append({
                        "margin": margin, "order": order, "dataset": dataset,
                        "subset_size": size,
                        **A.partition_metrics(restricted, subset),
                        "full_clusters": stats["final_clusters"],
                        "full_singletons": stats["singleton_clusters"],
                        "reassignment_moves": stats["reassignment_moves"],
                        "strict_merges": stats["strict_merges"],
                        "iterations": stats["iterations"],
                        "elapsed_seconds": stats["elapsed_seconds"],
                    })
    A.write_csv(root / "figures" / "margin_sweep.csv", rows)

    def mean(selected, key):
        values = [float(r[key]) for r in selected]
        return statistics.mean(values) if values else float("nan")

    print()
    print("Reassignment hysteresis sweep (graph, seed 0.40, %d datasets)" % datasets)
    for order in V.ORDERS:
        print()
        print(f"=== {order} order ===")
        print("  %-8s %10s %10s %10s %9s %9s %8s" % (
            "margin", "Jacc 10%", "Jacc 40%", "Jacc 80%", "clusters", "moves", "sec"))
        for margin in MARGINS:
            sel = [r for r in rows if r["margin"] == margin and r["order"] == order]
            if not sel:
                continue
            by_size = {s: [r for r in sel if r["subset_size"] == s] for s in SUBSET_SIZES}
            print("  %-8s %10.4f %10.4f %10.4f %9.0f %9.0f %8.2f" % (
                margin,
                mean(by_size[1_000], "pairwise_jaccard"),
                mean(by_size[4_000], "pairwise_jaccard"),
                mean(by_size[8_000], "pairwise_jaccard"),
                mean(sel, "full_clusters"),
                mean(sel, "reassignment_moves"),
                mean(sel, "elapsed_seconds")))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=8)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--analyse-only", action="store_true")
    args = parser.parse_args()
    root, binary = args.root.resolve(), args.binary.resolve()

    if not args.analyse_only:
        tasks = [
            (root, binary, margin, order, dataset, args.threads, size)
            for margin in MARGINS
            for order in V.ORDERS
            for dataset in range(args.datasets)
            for size in (None, *SUBSET_SIZES)
        ]
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, *task) for task in tasks]
            for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
                if index % 25 == 0 or index == len(tasks):
                    print(f"[{index}/{len(tasks)}] {future.result()}", flush=True)
                else:
                    future.result()
    analyse(root, args.datasets)


if __name__ == "__main__":
    main()
