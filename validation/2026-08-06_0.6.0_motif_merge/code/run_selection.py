#!/usr/bin/env python3
"""Does the representative-selection rule matter against the allele labels?

Greedy clustering picks representatives in some order and the order decides the
partition. Three rules exist and only two have ever been scored here:

  graph        dynamic set cover over the materialised edge graph
  lazy_exact   the same selection without storing the graph, by evaluating exact
               coverage before committing a representative
  kmer_degree  a cheap static ordering by the number of distinct k-mer candidates

`kmer_degree` is the default for `--clustering-method greedy` and is roughly 18x
faster than `lazy_exact`, because `lazy_exact` re-retrieves and rescores a
candidate list every time a node is reinserted into its priority queue. It has
never been benchmarked against the allele labels: every accuracy figure in this
study used `graph` or `lazy_exact`. If it scores comparably, large runs get that
speed for free.

Each run reports both partitions from one invocation: the similarity clustering
itself, and the motif layer built on top of it. The second matters because the
motif layer is where the accuracy is, and a coarser representative order might be
absorbed by it.
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

SCORING = ["--mode", "separate_kmer_anchor",
           "--kmer-similarity-threshold", "0.25",
           "--anchor-combination-similarity-threshold", "0.35"]

SELECTION = {
    "graph": ["--clustering-method", "graph", "--representative-order", "coverage",
              "--no-prefilter"],
    "lazy_exact": ["--clustering-method", "greedy", "--greedy-selection", "lazy-exact"],
    "kmer_degree": ["--clustering-method", "greedy", "--greedy-selection", "kmer-degree"],
}


def run_one(job, binary, pools, tmp_root, selected):
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="pc2sel_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [str(binary), "--input", str(pools / f"{job['pool']}.fasta"),
                   "--output-dir", str(output), *SCORING,
                   *SELECTION[job["selection"]],
                   "--compact-output", "--threads", "1",
                   "--tmp-dir", str(Path(tmp) / "tmp"),
                   "--merge-motifs",
                   "--motif-prior-concentration", str(selected["merge_concentration"]),
                   "--motif-merge-threshold", str(selected["merge_threshold"]),
                   "--motif-em-prior-concentration", str(selected["em_concentration"])]
        proc = subprocess.run(command, capture_output=True, text=True)
        if proc.returncode:
            return {**job, "status": "failed",
                    "error": (proc.stderr or proc.stdout)[-300:].replace("\n", " ")}
        motifs = pd.read_csv(output / "motif_clusters.tsv", sep="\t")
        stats = json.loads((output / "run_stats.json").read_text())
    elapsed = time.time() - started

    merged = labels.merge(motifs[["sequence", "motif_id", "similarity_cluster_id"]],
                          left_on="peptide", right_on="sequence", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"assigned {len(merged)} of {len(labels)}"}
    allele = merged["allele"].to_numpy()
    similarity = M.evaluate(allele, merged["similarity_cluster_id"].to_numpy())
    motif = M.evaluate(allele, merged["motif_id"].to_numpy())
    row = {**job, "status": "ok", "error": "", "elapsed_seconds": round(elapsed, 3),
           "tool_seconds": stats.get("elapsed_seconds", float("nan")),
           "candidate_pairs_computed": stats.get("candidate_pairs_computed", 0),
           "index_candidate_hits": stats.get("index_candidate_hits", 0)}
    for prefix, scores in (("sim", similarity), ("motif", motif)):
        for key in ("ami", "adjusted_purity_macro", "bcubed_precision_macro",
                    "bcubed_recall_macro", "bcubed_f1_macro", "clusters"):
            row[f"{prefix}_{key}"] = scores[key]
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--selections", nargs="+", default=list(SELECTION))
    parser.add_argument("--workers", type=int, default=40)
    args = parser.parse_args()
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(args.selected).iloc[0].to_dict()
    manifest = pd.read_csv(args.manifest).sort_values("peptides", ascending=False)
    jobs = [{"pool": r.pool, "split": r.split, "allele_count": r.allele_count,
             "peptides": r.peptides, "selection": s}
            for r in manifest.itertuples() for s in args.selections]

    fields = list(jobs[0]) + ["status", "error", "elapsed_seconds", "tool_seconds",
                              "candidate_pairs_computed", "index_candidate_hits"]
    for prefix in ("sim", "motif"):
        fields += [f"{prefix}_{k}" for k in
                   ("ami", "adjusted_purity_macro", "bcubed_precision_macro",
                    "bcubed_recall_macro", "bcubed_f1_macro", "clusters")]
    lock = threading.Lock()
    done = [0]
    print(f"{len(jobs)} runs on {args.workers} workers", flush=True)
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
                    if done[0] % 50 == 0:
                        handle.flush()
                        print(f"  {done[0]}/{len(jobs)}", flush=True)
    print(f"wrote {done[0]} rows to {args.out}", flush=True)


if __name__ == "__main__":
    main()
