#!/usr/bin/env python3
"""One motif-layer variant over every benchmark pool.

Three supported variants, all sharing the base clustering and the EM
concentration that nested selection chose on the inner folds:

  merge_em   agglomerative merge, then EM seeded from the merged groups. The
             headline configuration.
  em_only    no merge; EM gets one component per similarity cluster and finds its
             own count. Statistically indistinguishable from merge_em on AMI and
             F1, but finer: higher precision, lower recall.
  forced_k   EM seeded from the K similarity clusters farthest apart in profile
             space, then exactly K motifs are returned.

`forced_k` is reported as a separate arm rather than folded into the headline,
because the benchmark's premise is that the allele count is unknown - the same
reason MixMHCp and GibbsCluster carry separate `forced k` arms. In practice a
sample's alleles are usually known from typing, so it is the arm a real user
would run; both readings are kept side by side.
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

import pandas as pd

HERE = Path(__file__).resolve()
STUDY = HERE.parents[2] / "2026-07-30_0.5.0-dev_search_redesign"
sys.path.insert(0, str(STUDY / "code" / "mhc_bench"))
import metrics as M  # noqa: E402

BASE = ["--mode", "separate_kmer_anchor",
        "--kmer-similarity-threshold", "0.25",
        "--anchor-combination-similarity-threshold", "0.35",
        "--representative-order", "coverage",
        "--clustering-method", "graph", "--no-prefilter", "--compact-output"]


VARIANTS = ("merge_em", "em_only", "forced_k")


def variant_flags(variant: str, selected: dict, allele_count: int) -> list[str]:
    flags = ["--merge-motifs",
             "--motif-em-prior-concentration", str(selected["em_concentration"])]
    if variant == "forced_k":
        return flags + ["--motif-count", str(int(allele_count))]
    if variant == "em_only":
        return flags + ["--no-motif-merge"]
    return flags + ["--motif-prior-concentration",
                    str(selected["merge_concentration"]),
                    "--motif-merge-threshold", str(selected["merge_threshold"])]


def run_one(job, binary, pools, tmp_root, selected):
    fasta = pools / f"{job['pool']}.fasta"
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="pc2fk_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [str(binary), "--input", str(fasta), "--output-dir", str(output),
                   *BASE, "--threads", "1", "--tmp-dir", str(Path(tmp) / "tmp"),
                   *variant_flags(job["variant"], selected, job["allele_count"])]
        proc = subprocess.run(command, capture_output=True, text=True)
        if proc.returncode:
            return {**job, "status": "failed",
                    "error": (proc.stderr or proc.stdout)[-300:].replace("\n", " ")}
        motifs = pd.read_csv(output / "motif_clusters.tsv", sep="\t")
        stats = json.loads((output / "run_stats.json").read_text())
    merged = labels.merge(motifs[["sequence", "motif_id"]],
                          left_on="peptide", right_on="sequence", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"assigned {len(merged)} of {len(labels)}"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["motif_id"].to_numpy())
    return {**job, "status": "ok", "error": "", **scores,
            "merged_groups": stats.get("motif_merged_groups", 0),
            "em_iterations": stats.get("motif_em_iterations", 0),
            "elapsed_seconds": round(time.time() - started, 3)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--variant", choices=VARIANTS, required=True)
    parser.add_argument("--workers", type=int, default=40)
    args = parser.parse_args()
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(args.selected).iloc[0].to_dict()
    manifest = pd.read_csv(args.manifest).sort_values("peptides", ascending=False)
    jobs = [{"pool": r.pool, "split": r.split, "outer_fold": r.outer_fold,
             "allele_count": r.allele_count, "peptides": r.peptides,
             "variant": args.variant}
            for r in manifest.itertuples()]

    import numpy as np
    probe = M.evaluate(np.array(["a", "b"]), np.array([0, 1]))
    fields = list(jobs[0]) + ["status", "error"] + list(probe) + \
             ["merged_groups", "em_iterations", "elapsed_seconds"]
    lock = threading.Lock()
    done = [0]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, j, args.binary.resolve(),
                                   args.pools.resolve(), args.tmp_root.resolve(),
                                   selected) for j in jobs]
            for future in concurrent.futures.as_completed(futures):
                row = future.result()
                with lock:
                    writer.writerow({k: row.get(k, "") for k in fields})
                    done[0] += 1
                    if done[0] % 40 == 0:
                        handle.flush()
                        print(f"  {done[0]}/{len(jobs)}", flush=True)
    print(f"wrote {done[0]} rows to {args.out}")


if __name__ == "__main__":
    main()
