#!/usr/bin/env python3
"""Hyperparameter grid for the motif layer, under the study's nested protocol.

The base clustering configuration is NOT re-tuned here. It was selected by the
earlier nested cross-validation (`runs/mhc_bench_sep_kmer_anchor/tables/
selected_overall.csv`: graph, coverage, primary 0.25, anchor 0.35, mode
separate_kmer_anchor) and is held fixed, so what follows selects the motif
parameters conditional on it. Re-tuning both at once on the same pools would
spend the inner folds twice.

Every pool of every split is run at every configuration. Selection then uses the
inner folds only, exactly as `analyse.py` does: for each outer fold, the
configuration with the best mean objective over the *other* inner folds is
evaluated on that fold's held-out-allele pools; a single overall selection over
all 120 inner pools is evaluated once on the 48 independent test pools.

The three swept parameters are separated deliberately. An earlier exploratory
sweep passed one concentration to both the merge likelihood and the EM
smoothing, which made the two stages' contributions unattributable.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import itertools
import json
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve()
STUDY = HERE.parents[2] / "2026-07-30_0.5.0-dev_search_redesign"
sys.path.insert(0, str(STUDY / "code" / "mhc_bench"))
import metrics as M  # noqa: E402

# Base clustering, fixed. Mirrors selected_overall.csv for the k-mer mode.
BASE = ["--mode", "separate_kmer_anchor",
        "--kmer-similarity-threshold", "0.25",
        "--anchor-combination-similarity-threshold", "0.35",
        "--representative-order", "coverage",
        "--clustering-method", "graph", "--no-prefilter",
        "--compact-output"]

MERGE_CONCENTRATION = (0.3, 1.0, 3.0, 10.0, 30.0)
EM_CONCENTRATION = (0.3, 1.0, 3.0, 10.0, 30.0)
MERGE_THRESHOLD = (0.0, 25.0)


def grid() -> list[dict]:
    out = []
    for a0, t in itertools.product(MERGE_CONCENTRATION, MERGE_THRESHOLD):
        # Merge-only, so the EM contribution can be read off directly. The EM
        # concentration does not apply, recorded as NaN rather than a value that
        # would imply it was used.
        out.append({"merge_concentration": a0, "merge_threshold": t,
                    "em": False, "em_concentration": float("nan")})
        for a0_em in EM_CONCENTRATION:
            out.append({"merge_concentration": a0, "merge_threshold": t,
                        "em": True, "em_concentration": a0_em})
    return out


def run_one(job: dict, binary: Path, pools: Path, tmp_root: Path) -> dict:
    fasta = pools / f"{job['pool']}.fasta"
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="pc2motif_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [str(binary), "--input", str(fasta), "--output-dir", str(output),
                   *BASE, "--threads", "1", "--tmp-dir", str(Path(tmp) / "tmp"),
                   "--merge-motifs",
                   "--motif-prior-concentration", f"{job['merge_concentration']}",
                   "--motif-merge-threshold", f"{job['merge_threshold']}"]
        if job["em"]:
            command += ["--motif-em-prior-concentration", f"{job['em_concentration']}"]
        else:
            command += ["--no-motif-em"]
        proc = subprocess.run(command, capture_output=True, text=True)
        if proc.returncode:
            tail = (proc.stderr or proc.stdout)[-300:].replace("\n", " ")
            return {**job, "status": "failed", "error": tail}
        motifs = pd.read_csv(output / "motif_clusters.tsv", sep="\t")
        stats = json.loads((output / "run_stats.json").read_text())

    merged = labels.merge(motifs[["sequence", "motif_id", "similarity_cluster_id"]],
                          left_on="peptide", right_on="sequence", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"assigned {len(merged)} of {len(labels)} peptides"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["motif_id"].to_numpy())
    similarity = M.evaluate(merged["allele"].to_numpy(),
                            merged["similarity_cluster_id"].to_numpy())
    return {
        **job, "status": "ok", "error": "", **scores,
        "objective": M.objective(scores),
        # The similarity partition is identical across motif configurations, but
        # carrying it per row makes every comparison paired without a join.
        "similarity_clusters": similarity["clusters"],
        "similarity_ami": similarity["ami"],
        "similarity_f1": similarity["bcubed_f1_macro"],
        "merged_groups": stats.get("motif_merged_groups", 0),
        "em_iterations": stats.get("motif_em_iterations", 0),
        "em_converged": stats.get("motif_em_converged", False),
        "elapsed_seconds": round(time.time() - started, 3),
    }


def schema(job_keys) -> list[str]:
    """Fixed output columns.

    Derived from metrics.evaluate rather than hand-listed, so a metric added
    there appears here automatically instead of silently vanishing. A failed run
    fills the metric columns with blanks; taking the columns from whichever row
    happened to be written first is what broke the first attempt.
    """
    import numpy as np
    probe = M.evaluate(np.array(["a", "b"]), np.array([0, 1]))
    extra = ["objective", "similarity_clusters", "similarity_ami", "similarity_f1",
             "merged_groups", "em_iterations", "em_converged", "elapsed_seconds"]
    return list(job_keys) + ["status", "error"] + list(probe) + extra


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--splits", nargs="+", default=["inner", "outer", "test"])
    parser.add_argument("--workers", type=int, default=40)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    args = parser.parse_args()
    if not 0 <= args.shard < args.shards:
        parser.error("--shard must be in [0, --shards)")
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest["split"].isin(args.splits)]
    # Largest pools first so the long tail is not left to the end.
    manifest = manifest.sort_values("peptides", ascending=False)
    configs = grid()
    jobs = [{"pool": p.pool, "split": p.split, "outer_fold": p.outer_fold,
             "allele_count": p.allele_count, "peptides": p.peptides, **c}
            for p in manifest.itertuples() for c in configs]
    # Round-robin, so every shard gets the same mix of pool sizes. Blocking them
    # contiguously would leave one shard with all the 25k-peptide pools.
    jobs = jobs[args.shard::args.shards]

    done = set()
    if args.out.exists():
        previous = pd.read_csv(args.out)
        done = {(r.pool, r.merge_concentration, r.merge_threshold, r.em,
                 r.em_concentration if r.em else None)
                for r in previous.itertuples()}
    todo = [j for j in jobs
            if (j["pool"], j["merge_concentration"], j["merge_threshold"], j["em"],
                j["em_concentration"] if j["em"] else None) not in done]

    print(f"{len(configs)} configurations x {len(manifest)} pools = {len(jobs)} runs; "
          f"{len(todo)} outstanding; {args.workers} workers", flush=True)
    lock = threading.Lock()
    written = [0]
    new = not args.out.exists()
    fields = schema(jobs[0].keys()) if jobs else []
    with open(args.out, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        if new:
            writer.writeheader()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, j, args.binary.resolve(),
                                   args.pools.resolve(), args.tmp_root.resolve())
                       for j in todo]
            for future in concurrent.futures.as_completed(futures):
                row = future.result()
                with lock:
                    writer.writerow({k: row.get(k, "") for k in fields})
                    written[0] += 1
                    # Flush often: a shard that dies late should not lose the
                    # work it already did, since the runner is resumable.
                    if written[0] % 50 == 0:
                        handle.flush()
                        print(f"  {written[0]}/{len(todo)}", flush=True)
    print(f"wrote {written[0]} rows to {args.out}", flush=True)


if __name__ == "__main__":
    main()
