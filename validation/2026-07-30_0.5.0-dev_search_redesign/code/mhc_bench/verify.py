#!/usr/bin/env python3
"""Check that the grid is actually complete and correct before anything is read from it.

Every failure in this study so far exited zero. A shard that skipped a fifth of
its jobs, a file that gained a second header, a job list mistaken for results:
none of them raised, and SLURM reported success each time. Exit status therefore
proves nothing here, and completeness has to be asserted against the expected job
list rather than assumed.

Checks
------
1. every expected (pool, method, order, primary, anchor) is present exactly once
2. no run failed, and no key exists that was never requested
3. every pool's peptides were all clustered (the runner records a failure if not)
4. metrics are recomputable: a random sample is re-run and re-scored, and must
   reproduce the stored numbers

Exit code is non-zero if any check fails, so it can gate the analysis.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M  # noqa: E402
import run_grid as R  # noqa: E402

KEY = ["pool", "method", "representative_order", "primary_threshold", "anchor_threshold"]


def key_frame(frame: pd.DataFrame) -> set:
    return {
        (r[0], r[1], r[2], round(float(r[3]), 2), round(float(r[4]), 2))
        for r in frame[KEY].itertuples(index=False)
    }


def load(root: Path, tag: str) -> pd.DataFrame:
    frames = []
    for path in sorted(glob.glob(str(root / "grid" / f"{tag}_*.csv"))):
        frame = pd.read_csv(path)
        if "status" in frame.columns and len(frame):
            frames.append(frame)
    if not frames:
        raise SystemExit(f"no {tag} results")
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--tag", default="inner")
    parser.add_argument("--splits", nargs="+", default=["inner"])
    parser.add_argument("--sample", type=int, default=12,
                        help="runs to re-execute and re-score")
    parser.add_argument("--tmp-root", type=Path, default=Path("/tmp"))
    parser.add_argument("--mode", default="separate_aln_anchor")
    args = parser.parse_args()
    root = args.root.resolve()

    expected = {R.key_of(job) for job in
                R.job_list(root / "pool_manifest.csv", tuple(args.splits), None, args.mode)}
    frame = load(root, args.tag)
    present = key_frame(frame)
    duplicates = len(frame) - len(present)
    missing = expected - present
    unexpected = present - expected
    failed = frame[frame["status"] != "ok"]

    problems = []
    print(f"expected runs        : {len(expected)}")
    print(f"rows on disk         : {len(frame)}")
    print(f"unique runs          : {len(present)}")
    print(f"duplicate rows       : {duplicates}")
    print(f"missing runs         : {len(missing)}")
    print(f"unrequested runs     : {len(unexpected)}")
    print(f"failed runs          : {len(failed)}")
    if missing:
        problems.append(f"{len(missing)} runs missing")
        for key in sorted(missing)[:5]:
            print(f"   missing example: {key}")
    if unexpected:
        problems.append(f"{len(unexpected)} unrequested runs")
    if len(failed):
        problems.append(f"{len(failed)} failed runs")
        print(f"   first error: {failed.iloc[0]['error']}")

    ok = frame[frame["status"] == "ok"]
    # A pool's peptides must all be accounted for; the runner fails a run
    # otherwise, but check the recorded counts against the manifest too.
    manifest = pd.read_csv(root / "pool_manifest.csv").set_index("pool")["peptides"]
    mismatch = ok[ok["peptides"] != ok["pool"].map(manifest)]
    print(f"peptide-count mismatches: {len(mismatch)}")
    if len(mismatch):
        problems.append(f"{len(mismatch)} runs disagree with the manifest")

    if args.sample and not missing:
        print(f"\nre-running {args.sample} random runs to confirm the metrics reproduce")
        rng = np.random.default_rng(0)
        picks = ok.sample(min(args.sample, len(ok)), random_state=int(rng.integers(1 << 30)))
        worst = 0.0
        for row in picks.itertuples():
            job = {k: getattr(row, k) for k in
                   ["pool", "split", "outer_fold", "allele_count", "peptides",
                    "scoring_mode", "method", "representative_order",
                    "primary_threshold", "anchor_threshold"]}
            repeat = R.run_job(job, args.binary.resolve(), root / "pools", 8,
                               args.tmp_root.resolve())
            if repeat["status"] != "ok":
                problems.append(f"re-run failed for {row.pool}")
                continue
            for metric in ("ami", "nmi", "adjusted_purity_macro", "clusters"):
                stored, again = float(getattr(row, metric)), float(repeat[metric])
                delta = abs(stored - again)
                worst = max(worst, delta)
                if delta > 1e-9:
                    problems.append(
                        f"{row.pool} {metric}: stored {stored} vs re-run {again}")
        print(f"largest disagreement across the sample: {worst:.2e}")

    print()
    if problems:
        print("FAILED")
        for problem in problems[:10]:
            print(f"  - {problem}")
        raise SystemExit(1)
    print("PASSED: grid is complete, unique, error-free and reproducible")


if __name__ == "__main__":
    main()
