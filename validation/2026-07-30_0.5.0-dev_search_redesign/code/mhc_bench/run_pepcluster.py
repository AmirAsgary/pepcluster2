#!/usr/bin/env python3
"""Score PepCluster on the identical pools, under the identical protocol.

The published PepCluster figures cannot be compared with ours directly. They were
produced on different pools (143 alleles, natural imbalance, sizes up to the full
mixture), and above all under a different correctness rule: a peptide counted as
correct there when its allele and its cluster label were at least 99% similar, so
near-identical alleles scored as successes. This benchmark instead removes such
alleles from ever sharing a pool. A "soft" accuracy of that kind is necessarily
higher than the strict agreement measured here, and their AMI is computed over a
different label universe.

The only defensible comparison is to run PepCluster on these pools and score it
with these metrics, giving it the same tuning protocol our own modes received:
sweep its threshold on the inner folds, select by AMI under the singleton
constraint, then evaluate the selection on the held-out alleles and the
independent test pools.

The anchor weight is swept as well, because it is the mechanism the PepCluster
authors credit for the configuration working well: weight 1.0 treats the six
terminal positions equally, while 2.0 with `a_22` doubles the second position of
each terminus, which is where the MHC binding anchors sit.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M  # noqa: E402

THRESHOLDS = (0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80)
ANCHOR_WEIGHTS = (1.0, 2.0)
ANCHORS = "2;2"


def run_job(job: dict, pools: Path, tmp_root: Path) -> dict:
    fasta = pools / f"{job['pool']}.fasta"
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="pcold_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [
            "pepcluster", "-i", str(fasta), "-o", str(output),
            "-t", f"{job['threshold']:.2f}",
            "--anchors", job["anchors"],
            "--anchor-weight", f"{job['anchor_weight']:.1f}",
            "--threads", "1", "-q",
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode:
            return {**job, "status": "failed",
                    "error": (result.stderr or result.stdout)[-300:].replace("\n", " ")}
        assignments = pd.read_csv(output / "clusters.tsv", sep="\t")

    merged = labels.merge(assignments[["sequence", "cluster_id"]].drop_duplicates("sequence"),
                          left_on="peptide", right_on="sequence", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"clustered {len(merged)} of {len(labels)} peptides"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["cluster_id"].to_numpy())
    return {**job, "status": "ok", "error": "", **scores,
            "objective": M.objective(scores), "elapsed_seconds": time.time() - started}


def job_list(manifest: Path, splits: tuple[str, ...], configs=None) -> list[dict]:
    pools = pd.read_csv(manifest)
    pools = pools[pools["split"].isin(splits)].sort_values("peptides", ascending=False)
    if configs is not None:
        grid = [(float(r.threshold), float(r.anchor_weight), r.anchors)
                for r in configs.drop_duplicates(
                    ["threshold", "anchor_weight", "anchors"]).itertuples()]
    else:
        grid = [(t, w, ANCHORS) for w in ANCHOR_WEIGHTS for t in THRESHOLDS]
    return [
        {"pool": p.pool, "split": p.split, "outer_fold": p.outer_fold,
         "allele_count": p.allele_count, "peptides": p.peptides,
         "method": "pepcluster", "threshold": t, "anchor_weight": w, "anchors": a}
        for p in pools.itertuples() for t, w, a in grid
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--splits", nargs="+", default=["inner"])
    parser.add_argument("--configs", type=Path)
    parser.add_argument("--cores", type=int, default=64)
    parser.add_argument("--tmp-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    args.tmp_root.mkdir(parents=True, exist_ok=True)

    configs = pd.read_csv(args.configs) if args.configs else None
    jobs = job_list(root / "pool_manifest.csv", tuple(args.splits), configs)
    destination = root / "grid" / f"pepcluster_{'_'.join(args.splits)}.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)

    handle = destination.open("w", newline="")
    writer = None
    lock = threading.Lock()
    counter = {"done": 0}
    started = time.time()

    def execute(job):
        nonlocal writer
        row = run_job(job, root / "pools", args.tmp_root.resolve())
        with lock:
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
            writer.writerow({k: row.get(k, "") for k in writer.fieldnames})
            handle.flush()
            counter["done"] += 1
            if counter["done"] % 200 == 0:
                print(f"[{counter['done']}/{len(jobs)}]", flush=True)

    print(f"{len(jobs)} PepCluster runs on {args.cores} cores", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cores) as pool:
        list(pool.map(execute, jobs))
    handle.close()
    print(f"complete in {(time.time() - started) / 60:.1f} min -> {destination}", flush=True)


if __name__ == "__main__":
    main()
