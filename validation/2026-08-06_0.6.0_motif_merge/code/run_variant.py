#!/usr/bin/env python3
"""One motif-layer mode over every benchmark pool.

Eight arms, all sharing the base clustering and the merge and EM parameters that
nested selection chose on the inner folds, so the only thing that differs between
them is the stage under test:

  similarity   clustering only; the motif layer is off.
  merge        agglomerative merge of the clusters, no EM.
  merge_k      the merge continued until exactly K groups remain, no EM. With
               `merge_em_k` this is the clean isolation of EM: identical merge,
               identical count, EM off then on.
  merge_em     merge, then EM seeded from the merged groups. The headline
               configuration, and the only one that chooses its own count.
  merge_em_k   merge to exactly K, then EM from those K groups.
  em_only      no merge; EM gets one component per similarity cluster and finds
               its own count.
  em_only_k    no merge; EM seeded from the K similarity clusters farthest apart
               in profile space, then exactly K motifs are returned. The best
               arm in the study.
  merge_em_k_merged
               a second merge-seeded given-K route: merge at the evidence
               threshold, then take K seeds farthest apart among the merged
               groups. Reported in the tables, kept out of the figures, because
               merge_k -> merge_em_k is the nested pair that isolates EM.

The `given K` arms are reported separately from the automatic ones, for every
tool including ours, because the benchmark's premise is that the allele count is
unknown. In practice a sample's alleles are usually known from typing, so they
are the arms a real user would run; both readings are kept side by side.
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


# The seven reported modes, plus merge_em_k_merged, which seeds a requested
# count from groups merged at the evidence threshold rather than from groups
# merged down to exactly K. Both are merge-seeded given-K; the nested pair
# merge_k -> merge_em_k is the one carried in the figures, because it isolates
# the EM stage against an otherwise identical merge.
VARIANTS = ("similarity", "merge", "merge_k", "merge_em", "merge_em_k",
            "em_only", "em_only_k", "merge_em_k_merged")

# Display names, used by every figure, CSV and table downstream.
LABELS = {
    "similarity":        "PepCluster2",
    "merge":             "PepCluster2 + merge",
    "merge_k":           "PepCluster2 + merge (given K)",
    "merge_em":          "PepCluster2 + merge + EM",
    "merge_em_k":        "PepCluster2 + merge + EM (given K)",
    "em_only":           "PepCluster2 + EM",
    "em_only_k":         "PepCluster2 + EM (given K)",
    "merge_em_k_merged": "PepCluster2 + merge + EM (given K, threshold-merged)",
}


def variant_flags(variant: str, selected: dict, allele_count: int) -> list[str]:
    """Flags for one mode. `similarity` returns none: the motif layer is off."""
    if variant == "similarity":
        return []
    flags = ["--merge-motifs",
             "--motif-em-prior-concentration", str(selected["em_concentration"])]
    merge_params = ["--motif-prior-concentration",
                    str(selected["merge_concentration"]),
                    "--motif-merge-threshold", str(selected["merge_threshold"])]
    k = str(int(allele_count))
    if variant == "em_only":
        return flags + ["--no-motif-merge"]
    if variant == "em_only_k":
        return flags + ["--motif-count", k, "--motif-count-from", "clusters"]
    if variant == "merge":
        return flags + merge_params + ["--no-motif-em"]
    if variant == "merge_k":
        # Stopping the agglomeration at K assigns every peptide, so this is the
        # one requested-count mode that needs no EM.
        return flags + merge_params + [
            "--no-motif-em", "--motif-count", k, "--motif-count-from", "merge-to-k"]
    if variant == "merge_em_k":
        return flags + merge_params + [
            "--motif-count", k, "--motif-count-from", "merge-to-k"]
    if variant == "merge_em_k_merged":
        return flags + merge_params + [
            "--motif-count", k, "--motif-count-from", "merged"]
    return flags + merge_params


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
        if job["variant"] == "similarity":
            motifs = pd.read_csv(output / "node_clusters.tsv", sep="\t")
            motifs = motifs.rename(columns={"cluster_id": "motif_id"})
        else:
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
             "variant": args.variant, "label": LABELS[args.variant]}
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
