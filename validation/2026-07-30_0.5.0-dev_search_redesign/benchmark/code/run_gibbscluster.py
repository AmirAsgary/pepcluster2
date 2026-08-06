#!/usr/bin/env python3
"""Run GibbsCluster on the benchmark pools and score it with our metrics.

Settings follow the GibbsCluster 2.0 preset for MHC class I ligands: core length
9, five seeds, trash cluster enabled at threshold 2, and the sequence-weighting
default. `default` lets GibbsCluster scan 1-6 groups and keep the best by its own
KLD score, which is what the web server does. `oracle_k` fixes the number of
groups to the true allele count, which no user could do in practice and is
reported only to separate the model from its model selection.

NOTE: unlike the MixMHCp runner, this has not been executed end to end, because
GibbsCluster is not freely downloadable - it is behind a DTU academic licence
form. The output parser below targets the documented `*.gibbs.ds.out` core file.
Run `--selftest` once the tool is installed to confirm the parse before trusting
a full sweep.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import re
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent / "code" / "mhc_bench"))
import metrics as M  # noqa: E402


def choose_group_count(output: Path) -> str | None:
    """Group count with the highest total KLD, GibbsCluster's own criterion.

    `images/gibbs.KLDvsClusters.tab` is a matrix, not a list of scores: the
    header row is the group counts tried, and row `n` holds the per-cluster KLD
    of the n-group solution, zero-padded to the widest row. The figure to
    maximise is therefore the row sum, not any single column.

        1   2         3         4
        1   6.380295  0         0         0
        2   3.853824  3.775969  0         0
        ...
    """
    table = output / "images" / "gibbs.KLDvsClusters.tab"
    if not table.exists():
        found = list(output.rglob("gibbs.KLDvsClusters.tab"))
        if not found:
            return None
        table = found[0]
    best, best_total = None, None
    for line in table.read_text().splitlines()[1:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            total = sum(float(x) for x in parts[1:])
        except ValueError:
            continue
        if best_total is None or total > best_total:
            best, best_total = parts[0], total
    return best


def read_assignments(output: Path) -> pd.DataFrame | None:
    """Peptide -> group from the chosen GibbsCluster core file.

    Format notes, all confirmed against gibbscluster 2.0f output rather than
    assumed. The header is the FIRST line and is not comment-prefixed:

        G Gn  Num  Sequence  Core o of ip IP il IL dp DP dl DL Annotation sS ...
        G  0    0  ALNNLLHSL ALNNLLHSL o  0 ip -99 ...

    The group column is `Gn` and is 0-based; `G` is a constant row marker, not
    data. Columns are located by header name so a version change fails loudly
    instead of being mis-parsed silently.

    Peptides can be absent from this file for two distinct reasons: the trash
    cluster removes outliers, and a peptide shorter than the motif length cannot
    host a core. Both leave the peptide unassigned, and the caller decides how to
    score that - this function only reports what the tool assigned.
    """
    root = output if (output / "res").exists() else output
    candidates = sorted(root.rglob("*.ds.out"))
    if not candidates:
        return None
    chosen = choose_group_count(output)
    if chosen:
        for path in candidates:
            if re.search(rf"[^0-9]{chosen}g\.", path.name):
                candidates = [path]
                break
    path = candidates[0]
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split()
    lower = [h.lower() for h in header]
    try:
        pep = next(i for i, h in enumerate(lower) if h == "sequence")
        grp = next(i for i, h in enumerate(lower) if h == "gn")
    except StopIteration:
        return None
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) > max(pep, grp):
            rows.append((parts[pep], parts[grp]))
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["peptide", "cluster"])


def run_one(job: dict, pools: Path, binary: Path, tmp_root: Path) -> dict:
    labels = pd.read_csv(pools / f"{job['pool']}.labels.tsv", sep="\t")
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="gibbs_", dir=tmp_root) as tmp:
        work = Path(tmp)
        peptides = work / "peptides.txt"
        peptides.write_text("\n".join(labels["peptide"].astype(str)) + "\n")
        output = work / "out"
        output.mkdir()
        command = [str(binary), "-f", str(peptides), "-P", "gibbs",
                   "-g", job["groups"], "-l", "9", "-S", "5",
                   "-T", "-j", str(job.get("trash", 2)), "-k", str(job.get("threads", 1)),
                   "-R", str(output)]
        # GibbsCluster is MCMC over a range of group counts with five seeds
        # each; on the largest pools that is hours. A per-run limit keeps one
        # pathological pool from holding a whole shard, and a timeout is recorded
        # as a failure with its own reason so it stays visible in the results.
        try:
            proc = subprocess.run(command, capture_output=True, text=True, cwd=work,
                                  timeout=job.get("timeout", 14400))
        except subprocess.TimeoutExpired:
            return {**job, "status": "failed",
                    "error": f"timed out after {job.get('timeout', 14400)}s"}
        assignments = read_assignments(output)
        if assignments is None:
            tail = (proc.stderr or proc.stdout)[-300:].replace("\n", " ")
            return {**job, "status": "failed", "error": f"no assignments; {tail}"}

    merged = labels.merge(assignments.drop_duplicates("peptide"), on="peptide", how="left")
    assigned = merged["cluster"].notna()
    if assigned.sum() == 0:
        return {**job, "status": "failed", "error": "no peptide was assigned"}

    # Two scorings, because they answer different questions and the difference
    # between them is exactly the tool's coverage limitation.
    #
    # full:     every peptide of the pool, with the peptides GibbsCluster did not
    #           place collected into one "unassigned" cluster. That is the honest
    #           representation of the tool's own output - its trash cluster is a
    #           cluster - and it keeps the denominator identical to every other
    #           tool in the comparison.
    # assigned: only the peptides it placed, which isolates the quality of the
    #           partition from how much of the pool it covers. Not comparable
    #           across tools, since the peptide set differs.
    full_labels = merged["cluster"].fillna("unassigned").to_numpy()
    scores = M.evaluate(merged["allele"].to_numpy(), full_labels)
    subset = merged[assigned]
    on_assigned = M.evaluate(subset["allele"].to_numpy(), subset["cluster"].to_numpy())
    return {**job, "status": "ok", "error": "", **scores,
            "coverage": round(assigned.sum() / len(merged), 6),
            "unassigned_peptides": int((~assigned).sum()),
            "assigned_ami": on_assigned["ami"],
            "assigned_adjusted_purity_macro": on_assigned["adjusted_purity_macro"],
            "assigned_bcubed_recall_macro": on_assigned["bcubed_recall_macro"],
            "assigned_bcubed_f1_macro": on_assigned["bcubed_f1_macro"],
            "assigned_clusters": on_assigned["clusters"],
            "elapsed_seconds": round(time.time() - started, 2)}


def build_jobs(manifest: Path, splits, settings) -> list[dict]:
    pools = pd.read_csv(manifest)
    pools = pools[pools["split"].isin(splits)].sort_values("peptides", ascending=False)
    jobs = []
    for pool in pools.itertuples():
        for setting in settings:
            groups = "1-6" if setting == "default" else str(int(pool.allele_count))
            jobs.append({"tool": "gibbscluster", "setting": setting, "pool": pool.pool,
                         "split": pool.split, "outer_fold": pool.outer_fold,
                         "allele_count": pool.allele_count, "peptides": pool.peptides,
                         "groups": groups})
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
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=14400,
                        help="per-run wall clock limit in seconds")
    parser.add_argument("--trash-threshold", type=int, default=2,
                        help="GibbsCluster -j; the tool's own default is 0")
    parser.add_argument("--gibbs-threads", type=int, default=1,
                        help="passed to GibbsCluster -k; it forks seeds x group "
                             "counts, so values above seeds*(maxg-ming+1) idle")
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--selftest", action="store_true",
                        help="run one small pool and print the parsed assignment")
    args = parser.parse_args()
    args.tmp_root.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(args.manifest, tuple(args.splits), tuple(args.settings))
    for job in jobs:
        job["timeout"] = args.timeout
        job["threads"] = args.gibbs_threads
        job["trash"] = args.trash_threshold
    # Round-robin so each shard gets the same mix of pool sizes.
    if args.shards > 1:
        jobs = jobs[args.shard::args.shards]
    if args.selftest:
        job = min(jobs, key=lambda j: j["peptides"])
        row = run_one(job, args.pools.resolve(), args.binary.resolve(), args.tmp_root.resolve())
        print({k: row[k] for k in ("pool", "peptides", "status", "error")})
        if row["status"] == "ok":
            print(f"  clusters={row['clusters']} ami={row['ami']:.4f} "
                  f"purity={row['adjusted_purity_macro']:.4f}")
            print(f"  coverage={row['coverage']:.4f} "
                  f"({row['unassigned_peptides']} peptides unassigned) "
                  f"assigned_ami={row['assigned_ami']:.4f}")
        raise SystemExit(0 if row["status"] == "ok" else 1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    print(f"{len(jobs)} GibbsCluster runs on {args.cores} cores", flush=True)
    handle = args.out.open("w", newline="")
    writer = None
    lock = threading.Lock()
    done = {"n": 0, "bad": 0}
    started = time.time()

    def execute(job):
        nonlocal writer
        row = run_one(job, args.pools.resolve(), args.binary.resolve(), args.tmp_root.resolve())
        with lock:
            if writer is None:
                # Schema from a successful row, not whichever finished first: a
                # failed run carries only the job keys, and taking the columns
                # from it makes every later success unwritable.
                import numpy as np
                probe = M.evaluate(np.array(["a", "b"]), np.array([0, 1]))
                fields = (list(jobs[0]) + ["status", "error"] + list(probe) +
                          ["coverage", "unassigned_peptides", "assigned_ami",
                           "assigned_adjusted_purity_macro",
                           "assigned_bcubed_recall_macro",
                           "assigned_bcubed_f1_macro", "assigned_clusters",
                           "elapsed_seconds"])
                writer = csv.DictWriter(handle, fieldnames=fields,
                                        extrasaction="ignore")
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
