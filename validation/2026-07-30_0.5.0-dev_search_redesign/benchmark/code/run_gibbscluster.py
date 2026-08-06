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


def read_assignments(output: Path) -> pd.DataFrame | None:
    """Peptide -> group, from the best-scoring GibbsCluster core file.

    GibbsCluster writes `res/<prefix>.<n>g.ds.out` with whitespace-separated
    columns whose header line begins with '#'. The peptide column and the group
    column are located by name so a version change in column order is caught
    rather than silently mis-parsed.
    """
    candidates = sorted((output / "res").glob("*.ds.out")) if (output / "res").exists() \
        else sorted(output.rglob("*.ds.out"))
    if not candidates:
        return None
    # Prefer the file for the number of groups GibbsCluster itself selected.
    best = output / "images" / "gibbs.KLDvsClusters.tab"
    chosen = None
    if best.exists():
        table = [l.split() for l in best.read_text().splitlines() if l.strip()
                 and not l.startswith("#")]
        if table:
            chosen = max(table, key=lambda r: float(r[-1]))[0]
    if chosen:
        for path in candidates:
            if re.search(rf"[^0-9]{chosen}g\.", path.name):
                candidates = [path]
                break

    rows = []
    for path in candidates[:1]:
        header, body = None, []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                header = line.lstrip("#").split()
            elif line.strip():
                body.append(line.split())
        if header is None or not body:
            return None
        lower = [h.lower() for h in header]
        try:
            pep = next(i for i, h in enumerate(lower) if "sequence" in h or "peptide" in h)
            grp = next(i for i, h in enumerate(lower) if h.startswith("cluster")
                       or h.startswith("group"))
        except StopIteration:
            return None
        for parts in body:
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
                   "-T", "-j", "2", "-R", str(output)]
        proc = subprocess.run(command, capture_output=True, text=True, cwd=work)
        assignments = read_assignments(output)
        if assignments is None:
            tail = (proc.stderr or proc.stdout)[-300:].replace("\n", " ")
            return {**job, "status": "failed", "error": f"no assignments; {tail}"}

    merged = labels.merge(assignments.drop_duplicates("peptide"), on="peptide", how="inner")
    if len(merged) != len(labels):
        return {**job, "status": "failed",
                "error": f"assigned {len(merged)} of {len(labels)} peptides"}
    scores = M.evaluate(merged["allele"].to_numpy(), merged["cluster"].to_numpy())
    return {**job, "status": "ok", "error": "", **scores,
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
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--selftest", action="store_true",
                        help="run one small pool and print the parsed assignment")
    args = parser.parse_args()
    args.tmp_root.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(args.manifest, tuple(args.splits), tuple(args.settings))
    if args.selftest:
        job = min(jobs, key=lambda j: j["peptides"])
        row = run_one(job, args.pools, args.binary, args.tmp_root.resolve())
        print({k: row[k] for k in ("pool", "peptides", "status", "error")})
        if row["status"] == "ok":
            print(f"  clusters={row['clusters']} ami={row['ami']:.4f} "
                  f"purity={row['adjusted_purity_macro']:.4f}")
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
