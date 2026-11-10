#!/usr/bin/env python3
"""Runtime against pool size, and runtime against accuracy.

Every run here is executed **serially**, one at a time on an otherwise idle node.
The elapsed times recorded during the grid and benchmark sweeps are not usable
for this: those ran 40-64 jobs concurrently, so their wall clocks include
contention and would flatter whichever tool happened to be scheduled loosely.

Threading is the other trap. PepCluster2 and MixMHCp are run single-threaded, so
their wall time is also their CPU cost. GibbsCluster is impractical
single-threaded - a 972-peptide pool did not finish in 35 minutes - so it is run
at `-k 16` and its CPU cost is reported as wall x 16. Both numbers are recorded;
`cpu_seconds` is the fair cross-tool comparison because it is the resource
actually consumed, and `wall_seconds` is what a user waits.

Nine arms. The seven reported PepCluster2 modes, named and flagged exactly as in
`run_variant.py`, so an arm here and an arm in the accuracy figures are the same
thing:

  similarity     clustering only, no motif layer
  merge          + merge, no EM
  merge_k        + merge continued to exactly K
  merge_em       + merge + EM, the automatic count
  merge_em_k     + merge to K + EM
  em_only        + EM seeded from the raw clusters
  em_only_k      + EM seeded from K of the raw clusters

plus the two comparators, each as documented:

  mixmhcp        scanning 1-6 motifs
  gibbscluster   scanning 1-6 groups, five seeds

Each (arm, pool) is timed `--repeats` times. One repeat is not enough: the cheap
PepCluster2 arms differ by less than the run-to-run spread, and an earlier
single-repeat version of this table put a given-K arm below clustering alone on
the median pool. `speed_figures.py` reduces the repeats to their median.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
STUDY = HERE.parents[2] / "2026-07-30_0.5.0-dev_search_redesign"
sys.path.insert(0, str(STUDY / "code" / "mhc_bench"))
import metrics as M  # noqa: E402
sys.path.insert(0, str(STUDY / "benchmark" / "code"))

BASE = ["--mode", "separate_kmer_anchor",
        "--kmer-similarity-threshold", "0.25",
        "--anchor-combination-similarity-threshold", "0.35",
        "--representative-order", "coverage",
        "--clustering-method", "graph", "--no-prefilter", "--compact-output"]
GIBBS_THREADS = 16


def score(labels, assignment, key):
    merged = labels.merge(assignment, left_on="peptide", right_on="sequence"
                          if "sequence" in assignment else "peptide", how="inner")
    if len(merged) < len(labels):          # unassigned peptides form one cluster
        merged = labels.merge(assignment, left_on="peptide",
                              right_on="sequence" if "sequence" in assignment
                              else "peptide", how="left")
        merged[key] = merged[key].fillna("unassigned")
    return M.evaluate(merged["allele"].to_numpy(), merged[key].to_numpy())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--mixmhcp", type=Path, required=True)
    parser.add_argument("--gibbscluster", type=Path, required=True)
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n-pools", type=int, default=8)
    # The previous table was one run per pool. On sub-second PepCluster2 timings
    # that is dominated by run-to-run variation: it made a requested count look
    # cheaper than clustering alone at the median pool, which is impossible,
    # because given-K does all the clustering work and then more.
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--arms", default="all",
                        choices=("all", "fast", "gibbscluster"),
                        help="fast = PepCluster2 and MixMHCp; gibbscluster "
                             "dominates the cost and can be run separately")
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.mixmhcp = args.mixmhcp.resolve()
    args.gibbscluster = args.gibbscluster.resolve()
    args.pools = args.pools.resolve()
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(args.selected).iloc[0].to_dict()

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest.split == "test"].sort_values("peptides")
    # log-spaced in pool size so the cost curve is sampled evenly on the axis it
    # will be plotted on
    targets = np.geomspace(manifest.peptides.min(), manifest.peptides.max(),
                           args.n_pools)
    picks = manifest.iloc[[int(np.abs(manifest.peptides - t).argmin())
                           for t in targets]].drop_duplicates("pool")

    from run_mixmhcp import read_assignments as mix_read
    from run_gibbscluster import read_assignments as gibbs_read
    from run_variant import VARIANTS, LABELS, variant_flags

    do_fast = args.arms in ("all", "fast")
    do_gibbs = args.arms in ("all", "gibbscluster")
    # The same flag builder the accuracy benchmark uses, so a mode cannot be
    # timed under different settings than it was scored under.
    pc2_arms = [v for v in VARIANTS if v != "merge_em_k_merged"]

    rows = []
    for record in picks.itertuples():
        labels = pd.read_csv(args.pools / f"{record.pool}.labels.tsv", sep="\t")
        fasta = args.pools / f"{record.pool}.fasta"

        for repeat in range(args.repeats if do_fast else 0):
          for arm in pc2_arms:
            extra = variant_flags(arm, selected, record.allele_count)
            with tempfile.TemporaryDirectory(dir=args.tmp_root) as tmp:
                out = Path(tmp) / "out"
                started = time.time()
                subprocess.run([str(args.binary), "--input", str(fasta),
                                "--output-dir", str(out), *BASE, "--threads", "1",
                                "--tmp-dir", str(Path(tmp) / "tmp"), *extra],
                               capture_output=True, text=True, check=True)
                wall = time.time() - started
                if arm == "similarity":
                    a = pd.read_csv(out / "node_clusters.tsv", sep="\t")
                    s_ = score(labels, a[["sequence", "cluster_id"]], "cluster_id")
                else:
                    a = pd.read_csv(out / "motif_clusters.tsv", sep="\t")
                    s_ = score(labels, a[["sequence", "motif_id"]], "motif_id")
            rows.append(dict(pool=record.pool, peptides=record.peptides,
                             alleles=record.allele_count, arm=arm,
                             label=LABELS[arm], repeat=repeat, threads=1,
                             wall_seconds=wall, cpu_seconds=wall,
                             ami=s_["ami"], bcubed_f1_macro=s_["bcubed_f1_macro"],
                             clusters=s_["clusters"]))
            print(f"{record.pool:10s} {record.peptides:6d} r{repeat} {arm:18s} "
                  f"{wall:8.2f}s  F1 {s_['bcubed_f1_macro']:.3f}", flush=True)

        for repeat in range(args.repeats if do_fast else 0):
          with tempfile.TemporaryDirectory(dir=args.tmp_root) as tmp:
              out = Path(tmp) / "out"
              started = time.time()
              subprocess.run([str(args.mixmhcp), "-i", str(fasta), "-o", str(out),
                              "-m", "6", "-m_min", "1", "-l", "0"],
                             capture_output=True, text=True)
              wall = time.time() - started
              a = mix_read(out)
          if a is not None:
              s = score(labels, a.drop_duplicates("peptide"), "cluster")
              rows.append(dict(pool=record.pool, peptides=record.peptides,
                               alleles=record.allele_count, arm="mixmhcp",
                               label="MixMHCp", repeat=repeat, threads=1,
                               wall_seconds=wall, cpu_seconds=wall, ami=s["ami"],
                               bcubed_f1_macro=s["bcubed_f1_macro"],
                               clusters=s["clusters"]))
              print(f"{record.pool:10s} {record.peptides:6d} {'mixmhcp':14s} "
                    f"{wall:8.2f}s  F1 {s['bcubed_f1_macro']:.3f}", flush=True)

        for repeat in range(args.repeats if do_gibbs else 0):
          with tempfile.TemporaryDirectory(dir=args.tmp_root) as tmp:
            work = Path(tmp)
            peptides = work / "peptides.txt"
            peptides.write_text("\n".join(labels.peptide.astype(str)) + "\n")
            out = work / "out"
            out.mkdir()
            started = time.time()
            subprocess.run([str(args.gibbscluster), "-f", str(peptides), "-P", "gibbs",
                            "-g", "1-6", "-l", "9", "-S", "5", "-T", "-j", "2",
                            "-k", str(GIBBS_THREADS), "-R", str(out)],
                           capture_output=True, text=True, cwd=work)
            wall = time.time() - started
            a = gibbs_read(out)
        if a is not None:
            merged = labels.merge(a.drop_duplicates("peptide"), on="peptide",
                                  how="left")
            merged["cluster"] = merged["cluster"].fillna("unassigned")
            s = M.evaluate(merged.allele.to_numpy(), merged.cluster.to_numpy())
            rows.append(dict(pool=record.pool, peptides=record.peptides,
                             alleles=record.allele_count, arm="gibbscluster",
                             label="GibbsCluster", repeat=repeat,
                             threads=GIBBS_THREADS, wall_seconds=wall,
                             cpu_seconds=wall * GIBBS_THREADS, ami=s["ami"],
                             bcubed_f1_macro=s["bcubed_f1_macro"],
                             clusters=s["clusters"]))
            print(f"{record.pool:10s} {record.peptides:6d} {'gibbscluster':14s} "
                  f"{wall:8.2f}s wall ({wall * GIBBS_THREADS:.0f} cpu-s)  "
                  f"F1 {s['bcubed_f1_macro']:.3f}", flush=True)

        pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"\nwrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
