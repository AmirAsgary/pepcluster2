#!/usr/bin/env python3
"""Combine the per-shard result files into one table per arm.

The runs are split over 20 array tasks, so each arm arrives as 20 fragments.
This checks that every pool is present exactly once before writing, because a
silently truncated arm would otherwise look like a real difference between
modes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

MODES = ["similarity", "merge", "merge_k", "merge_em", "merge_em_k",
         "em_only", "em_only_k"]


def combine(shards: list[Path], expected: set[str], name: str) -> pd.DataFrame | None:
    if not shards:
        print(f"  {name:12s} MISSING: no shard files")
        return None
    frame = pd.concat([pd.read_csv(f) for f in shards], ignore_index=True)
    got = set(frame.pool)
    missing, extra = expected - got, got - expected
    dupes = frame.pool.duplicated().sum() if "setting" not in frame else 0
    bad = (frame.status != "ok").sum()
    flag = "ok " if not (missing or dupes) else "!! "
    print(f"  {flag}{name:12s} {len(frame):5d} rows  {len(got):4d} pools  "
          f"{bad:3d} failed" + (f"  MISSING {len(missing)}" if missing else "")
          + (f"  DUPLICATED {dupes}" if dupes else "")
          + (f"  UNEXPECTED {len(extra)}" if extra else ""))
    if bad:
        for row in frame[frame.status != "ok"].head(3).itertuples():
            print(f"       {row.pool}: {str(row.error)[:140]}")
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    expected = set(pd.read_csv(args.manifest).pool)
    shard_dir = args.results / "shards"
    print(f"{len(expected)} pools expected\n")

    incomplete = []
    for mode in MODES:
        frame = combine(sorted(shard_dir.glob(f"variant_{mode}_s*.csv")),
                        expected, mode)
        if frame is None or set(frame.pool) != expected:
            incomplete.append(mode)
            continue
        frame.to_csv(args.results / f"variant_{mode}.csv", index=False)

    mix = combine(sorted(shard_dir.glob("mixmhcp_s*.csv")), expected, "mixmhcp")
    if mix is not None:
        # Two settings per pool, so every pool must appear exactly twice.
        per = mix.groupby("setting").pool.nunique().to_dict()
        print(f"       settings: {per}")
        if all(v == len(expected) for v in per.values()) and len(per) == 2:
            mix.to_csv(args.results / "mixmhcp.csv", index=False)
        else:
            incomplete.append("mixmhcp")
    else:
        incomplete.append("mixmhcp")

    if incomplete:
        print(f"\nincomplete, not written: {', '.join(incomplete)}")
        sys.exit(1)
    print(f"\nall arms complete -> {args.results}")


if __name__ == "__main__":
    main()
