#!/usr/bin/env python3
"""Where does our mixture's advantage over MixMHCp at matched k come from?

The initialisation control (`controls.py`) showed the PepCluster2 seed is worth
about +0.017 AMI. What remains is a +0.134 gap between our EM fitted at the true
allele count from a random start (0.625) and MixMHCp forced to the same k
(0.491). Both are mixtures of position weight matrices fitted by EM from random
initialisations, so that gap is the model, not the pipeline, and it was
previously unexplained.

Two candidate explanations, separated here:

  LENGTH   The pools are 61% 9-mers, 8% 8-mers and 31% 10-mers or longer. Our
           frame projects every length onto nine columns and gives an 8-mer a
           gap; MixMHCp handles other lengths its own way. If the gap is length
           handling, it should vanish on pools restricted to 9-mers, where the
           frame is the identity for both.

  SMOOTHING Our Dirichlet pseudocounts may simply be better tuned than
           MixMHCp's defaults. Swept independently on the full pools.

Design. Every arm fits the same number of components as the pool has alleles,
from random initialisations, choosing the maximum-likelihood restart - so no arm
is given information another lacks, and the comparison is at matched k.
"""

from __future__ import annotations

import argparse
import concurrent.futures
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

sys.path.insert(0, str(HERE.parent / "code"))
sys.path.insert(0, str(HERE.parent))
from motif_lib import (ALPHABET, background, count_matrix, encode,  # noqa: E402
                       random_seeded_em)

RESTARTS = 10
SEED_PEPS = 20


def read_mixmhcp_assignments(output: Path):
    """Mirror of run_mixmhcp.read_assignments, imported here to avoid a cycle."""
    sys.path.insert(0, str(STUDY / "benchmark" / "code"))
    from run_mixmhcp import read_assignments  # noqa: E402
    return read_assignments(output)


def run_mixmhcp(fasta: Path, k: int, binary: Path, tmp_root: Path):
    with tempfile.TemporaryDirectory(prefix="dec_mix_", dir=tmp_root) as tmp:
        out = Path(tmp) / "out"
        subprocess.run([str(binary), "-i", str(fasta), "-o", str(out),
                        "-m", str(k), "-m_min", str(k), "-l", "0"],
                       capture_output=True, text=True)
        return read_mixmhcp_assignments(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mixmhcp", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--concentrations", type=float, nargs="+",
                        default=[0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
    args = parser.parse_args()
    args.tmp_root.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest.split == "test"]

    def one_pool(row):
        labels = pd.read_csv(args.pools / f"{row.pool}.labels.tsv", sep="\t")
        k = int(row.allele_count)
        seed = abs(hash(row.pool)) % (2 ** 31)
        out = []

        for restriction in ("full", "9mer"):
            sub = labels if restriction == "full" else \
                labels[labels.peptide.str.len() == 9].reset_index(drop=True)
            if len(sub) < 50 or sub.allele.nunique() < 2:
                continue
            allele = sub.allele.to_numpy()
            X = encode(sub.peptide.values)
            bg = background(X)

            # ours, at the concentration the exploratory sweep favoured
            for a0 in (args.concentrations if restriction == "full" else [3.0]):
                lab, _ = random_seeded_em(X, k, a0 * bg * ALPHABET, seed,
                                          restarts=RESTARTS)
                s = M.evaluate(allele, lab)
                out.append(dict(pool=row.pool, alleles=k, restriction=restriction,
                                arm="ours", concentration=a0, peptides=len(sub),
                                ami=s["ami"], f1=s["bcubed_f1_macro"],
                                clusters=s["clusters"]))

            # MixMHCp at the same k on the same peptides
            with tempfile.TemporaryDirectory(prefix="dec_fa_", dir=args.tmp_root) as tmp:
                fasta = Path(tmp) / "p.fasta"
                fasta.write_text("".join(f">{i}\n{p}\n"
                                         for i, p in enumerate(sub.peptide)))
                assign = run_mixmhcp(fasta, k, args.mixmhcp.resolve(), args.tmp_root)
            if assign is not None:
                merged = sub.merge(assign.drop_duplicates("peptide"),
                                   on="peptide", how="inner")
                if len(merged) == len(sub):
                    s = M.evaluate(merged.allele.to_numpy(), merged.cluster.to_numpy())
                    out.append(dict(pool=row.pool, alleles=k, restriction=restriction,
                                    arm="mixmhcp", concentration=float("nan"),
                                    peptides=len(sub), ami=s["ami"],
                                    f1=s["bcubed_f1_macro"], clusters=s["clusters"]))
        return out

    rows = []
    started = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(one_pool, list(manifest.itertuples())):
            rows.extend(result)
            print(".", end="", flush=True)
    print(f"\n{time.time() - started:.0f}s")
    frame = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.out, index=False)
    print(f"wrote {len(frame)} rows to {args.out}")


if __name__ == "__main__":
    main()
