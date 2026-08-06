#!/usr/bin/env python3
"""Run MixMHCp on the benchmark pools and score it with our metrics.

Two settings, because they answer different questions.

`default` is MixMHCp as documented: `-m 6`, letting its KLD criterion pick the
number of motifs. This is what a user following the README would get.

`oracle_k` forces the number of motifs to the true number of alleles in the pool
(`-m_min = -m = k`). No user could do this without already knowing the answer, so
it is not a fair headline number - it is the tool's best case, included so that a
poor `default` result can be attributed to model selection rather than to the
motif model itself.

Assignment is the argmax over the responsibility columns of the selected model.
MixMHCp's `Trash` column is treated as one further cluster: it is a single
background component of the mixture, not a per-peptide reject.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent / "code" / "mhc_bench"))
import metrics as M  # noqa: E402

MAX_MOTIFS = 50  # MixMHCp's own ceiling


def read_assignments(output: Path) -> pd.DataFrame | None:
    """Peptide -> cluster, from the responsibility table of the selected model."""
    best = output / "KLD" / "best_ncl.txt"
    if not best.exists():
        return None
    first = best.read_text().splitlines()[0]
    ncl = int(first.split("\t")[1])
    table = output / "responsibility" / f"resp_{ncl}.txt"
    if not table.exists():
        return None
    frame = pd.read_csv(table, sep="\t")
    columns = [c for c in frame.columns
               if c == "Trash" or (c.isdigit() and c != "Peptide")]
    if not columns or "Peptide" not in frame.columns:
        return None
    assigned = frame[columns].astype(float).idxmax(axis=1)
    return pd.DataFrame({"peptide": frame["Peptide"], "cluster": assigned.astype(str)})


def run_one(job: dict, pools: Path, binary: Path, tmp_root: Path) -> dict:
    fasta = pools / f"{job['pool']}.fasta"
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="mixp_", dir=tmp_root) as tmp:
        output = Path(tmp) / "out"
        command = [str(binary), "-i", str(fasta), "-o", str(output),
                   "-m", str(job["max_motifs"]), "-m_min", str(job["min_motifs"]),
                   "-l", "0"]
        # The exit code is not trusted: MixMHCp returns non-zero when its optional
        # Rscript logo/length plots fail, long after clustering has succeeded.
        # Completeness of the responsibility table is the real gate.
        proc = subprocess.run(command, capture_output=True, text=True)
        assignments = read_assignments(output)
        if assignments is None:
            tail = (proc.stderr or proc.stdout)[-300:].replace("\n", " ")
            return {**job, "status": "failed", "error": f"no assignments; {tail}"}

    merged = labels.merge(assignments.drop_duplicates("peptide"),
                          on="peptide", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"assigned {len(merged)} of {len(labels)} peptides"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["cluster"].to_numpy())
    return {**job, "status": "ok", "error": "", **scores,
            "elapsed_seconds": round(time.time() - started, 2)}


def build_jobs(manifest: Path, splits: tuple[str, ...], settings: tuple[str, ...]) -> list[dict]:
    pools = pd.read_csv(manifest)
    pools = pools[pools["split"].isin(splits)].sort_values("peptides", ascending=False)
    jobs = []
    for pool in pools.itertuples():
        for setting in settings:
            if setting == "default":
                lo, hi = 1, 6
            else:
                lo = hi = min(int(pool.allele_count), MAX_MOTIFS)
            jobs.append({"tool": "mixmhcp", "setting": setting, "pool": pool.pool,
                         "split": pool.split, "outer_fold": pool.outer_fold,
                         "allele_count": pool.allele_count, "peptides": pool.peptides,
                         "min_motifs": lo, "max_motifs": hi})
    return jobs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--splits", nargs="+", default=["inner", "outer", "test"])
    parser.add_argument("--settings", nargs="+", default=["default", "oracle_k"])
    parser.add_argument("--cores", type=int, default=32)
    parser.add_argument("--tmp-root", type=Path, required=True)
    args = parser.parse_args()

    args.tmp_root.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    jobs = build_jobs(args.manifest, tuple(args.splits), tuple(args.settings))
    print(f"{len(jobs)} MixMHCp runs on {args.cores} cores", flush=True)

    handle = args.out.open("w", newline="")
    writer = None
    lock = threading.Lock()
    done = {"n": 0, "bad": 0}
    started = time.time()

    def execute(job):
        nonlocal writer
        row = run_one(job, args.pools, args.binary, args.tmp_root.resolve())
        with lock:
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row))
                writer.writeheader()
            writer.writerow({k: row.get(k, "") for k in writer.fieldnames})
            handle.flush()
            done["n"] += 1
            done["bad"] += row["status"] != "ok"
            if done["n"] % 25 == 0:
                print(f"  [{done['n']}/{len(jobs)}] failures={done['bad']}", flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.cores) as pool:
        list(pool.map(execute, jobs))
    handle.close()
    print(f"complete in {(time.time()-started)/60:.1f} min, "
          f"{done['bad']} failures -> {args.out}", flush=True)
    if done["bad"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
