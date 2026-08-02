#!/usr/bin/env python3
"""Run one shard of the peptide-MHC hyperparameter grid.

The job list is deterministic, so a SLURM array task only needs its own index to
know what to do, and a completed shard is skipped on resubmission. Jobs are
assigned round-robin so every shard gets the same mix of pool sizes; blocking
them contiguously would leave one shard with all the 100k pools.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M  # noqa: E402

# 0.15 and 0.25 were added after the first sweep selected 0.35, the lowest value
# then available: AMI was still rising at the boundary, so the optimum lay
# outside the grid.
THRESHOLDS = (0.15, 0.25, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)
# Zero disables the anchor condition entirely, isolating the primary component.
ANCHOR_THRESHOLDS = (0.0,) + THRESHOLDS
# lazy-exact is itself a dynamic set-cover rule, so it has no intrinsic variant.
METHODS = {
    "graph": (["--clustering-method", "graph", "--no-prefilter"], ("coverage", "intrinsic")),
    "greedy_lazy": (
        ["--clustering-method", "greedy", "--greedy-selection", "lazy-exact", "--no-prefilter"],
        ("coverage",),
    ),
}
# Each separate mode thresholds a different primary component against the anchor.
PRIMARY_FLAG = {
    "separate_aln_anchor": "--alignment-similarity-threshold",
    "separate_kmer_anchor": "--kmer-similarity-threshold",
}


def job_list(manifest: Path, splits: tuple[str, ...],
             configs: pd.DataFrame | None = None, mode: str = "separate_aln_anchor") -> list[dict]:
    pools = pd.read_csv(manifest)
    pools = pools[pools["split"].isin(splits)]
    # Largest pools first: with round-robin sharding this keeps the long tail
    # spread out instead of landing at the end of one shard.
    pools = pools.sort_values("peptides", ascending=False)
    if configs is not None:
        # Several folds usually select the same configuration; run each once.
        unique = configs.drop_duplicates(
            ["method", "representative_order", "primary_threshold", "anchor_threshold"])
        grid = [(r.method, r.representative_order,
                 float(r.primary_threshold), float(r.anchor_threshold))
                for r in unique.itertuples()]
    else:
        grid = [(method, order, primary, anchor)
                for method, (_, orders) in METHODS.items()
                for order in orders
                for primary in THRESHOLDS
                for anchor in ANCHOR_THRESHOLDS]
    jobs = []
    for pool in pools.itertuples():
        for method, order, primary, anchor in grid:
            jobs.append({
                "pool": pool.pool, "split": pool.split,
                "outer_fold": pool.outer_fold,
                "allele_count": pool.allele_count, "peptides": pool.peptides,
                "scoring_mode": mode, "method": method, "representative_order": order,
                "primary_threshold": primary, "anchor_threshold": anchor,
            })
    return jobs


def run_job(job: dict, binary: Path, pools: Path, threads: int, tmp_root: Path) -> dict:
    fasta = pools / f"{job['pool']}.fasta"
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="pc2mhc_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [
            str(binary), "--input", str(fasta), "--output-dir", str(output),
            "--mode", job["scoring_mode"],
            PRIMARY_FLAG[job["scoring_mode"]], f"{job['primary_threshold']:.2f}",
            "--anchor-combination-similarity-threshold", f"{job['anchor_threshold']:.2f}",
            "--representative-order", job["representative_order"],
            "--threads", str(threads), "--candidate-buffer-mb", "512",
            "--compact-output", "--tmp-dir", str(Path(tmp) / "tmp"),
            *METHODS[job["method"]][0],
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            return {**job, "status": "failed",
                    "error": (result.stderr or result.stdout)[-300:].replace("\n", " ")}
        assignments = pd.read_csv(output / "node_clusters.tsv", sep="\t")
        stats = json.loads((output / "run_stats.json").read_text())

    merged = labels.merge(assignments[["sequence", "cluster_id"]],
                          left_on="peptide", right_on="sequence", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"clustered {len(merged)} of {len(labels)} peptides"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["cluster_id"].to_numpy())
    return {
        **job, "status": "ok", "error": "", **scores,
        "objective": M.objective(scores),
        "elapsed_seconds": time.time() - started,
        "tool_seconds": stats["elapsed_seconds"],
        "graph_edges": stats.get("graph_edge_count", 0),
    }


LARGE_POOL_PEPTIDES = 20_000
LARGE_POOL_THREADS = 8


def completed_keys(grid: Path) -> set:
    """Every run already recorded, from any shard layout, so the task count can
    change between submissions without recomputing or double-counting."""
    done = set()
    for path in list(grid.glob("*.csv")) + list(grid.glob("*.partial")):
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if frame.empty or "status" not in frame.columns:
            continue
        if "primary_threshold" not in frame.columns:
            continue
        for row in frame[["pool", "method", "representative_order",
                          "primary_threshold", "anchor_threshold"]].itertuples(index=False):
            try:
                done.add((row[0], row[1], row[2],
                          round(float(row[3]), 2), round(float(row[4]), 2)))
            except (TypeError, ValueError):
                continue
    return done


def key_of(job: dict) -> tuple:
    return (job["pool"], job["method"], job["representative_order"],
            round(float(job["primary_threshold"]), 2),
            round(float(job["anchor_threshold"]), 2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--splits", nargs="+", default=["inner"])
    parser.add_argument("--configs", type=Path)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--mode", default="separate_aln_anchor", choices=list(PRIMARY_FLAG))
    parser.add_argument("--todo", type=Path,
                        help="fixed job list to shard; keeps every task consistent")
    parser.add_argument("--write-todo", type=Path,
                        help="compute the outstanding runs, write them here and exit")
    parser.add_argument("--attempt", default="0",
                        help="unique tag per submission, so result files never collide")
    args = parser.parse_args()

    root = args.root.resolve()
    binary = args.binary.resolve()
    pools = root / "pools"
    args.tmp_root.mkdir(parents=True, exist_ok=True)

    configs = pd.read_csv(args.configs) if args.configs else None
    tag = "_".join(args.splits)
    grid = root / "grid"
    grid.mkdir(parents=True, exist_ok=True)

    # The work list must be identical for every task. Computing "what is left" at
    # task start looks like sensible load balancing but is a race: tasks start
    # minutes to hours apart, each sees a different amount of finished work, and
    # their index%shards slices then cover different job sets. The union misses
    # jobs, silently, and every task still exits successfully. So the remaining
    # work is snapshotted to a file once, before the array is submitted, and every
    # task shards that same file.
    if args.write_todo:
        jobs = job_list(root / "pool_manifest.csv", tuple(args.splits), configs, args.mode)
        done = completed_keys(grid)
        remaining = [job for job in jobs if key_of(job) not in done]
        pd.DataFrame(remaining).to_csv(args.write_todo, index=False)
        print(f"wrote {len(remaining)} outstanding runs to {args.write_todo}", flush=True)
        return

    if args.todo:
        jobs = pd.read_csv(args.todo).to_dict("records")
    else:
        jobs = job_list(root / "pool_manifest.csv", tuple(args.splits), configs, args.mode)
    mine = [job for index, job in enumerate(jobs) if index % args.shards == args.shard]
    done = completed_keys(grid)
    mine = [job for job in mine if key_of(job) not in done]
    if not mine:
        print(f"shard {args.shard}: nothing left to do", flush=True)
        return

    destination = grid / f"{tag}_shard_{args.shard:04d}_{args.attempt}.partial"
    handle = destination.open("w", newline="")
    writer = None
    lock = threading.Lock()
    counter = {"done": 0}
    started = time.time()

    def execute(job, threads):
        nonlocal writer
        row = run_job(job, binary, pools, threads, args.tmp_root.resolve())
        with lock:
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
            writer.writerow({k: row.get(k, "") for k in writer.fieldnames})
            handle.flush()
            counter["done"] += 1
            if counter["done"] % 100 == 0:
                rate = (time.time() - started) / counter["done"]
                print(f"[{counter['done']}/{len(mine)}] {rate:.2f}s/job, "
                      f"eta {rate * (len(mine) - counter['done']) / 3600:.1f}h", flush=True)

    # Large pools first, few at a time with several threads each: one 100k run is
    # slow enough to become the critical path if left to the end. Small pools then
    # run single-threaded and wide, where fixed per-run costs dominate.
    large = [job for job in mine if job["peptides"] >= LARGE_POOL_PEPTIDES]
    small = [job for job in mine if job["peptides"] < LARGE_POOL_PEPTIDES]
    print(f"shard {args.shard}: {len(large)} large + {len(small)} small on "
          f"{args.cores} cores", flush=True)

    for batch, threads in ((large, LARGE_POOL_THREADS), (small, 1)):
        if not batch:
            continue
        workers = max(1, args.cores // threads)
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda job: execute(job, threads), batch))

    handle.close()
    destination.rename(grid / f"{tag}_shard_{args.shard:04d}_{args.attempt}.csv")
    print(f"shard {args.shard} complete: {len(mine)} runs in "
          f"{(time.time() - started) / 3600:.2f}h", flush=True)


if __name__ == "__main__":
    main()
