#!/usr/bin/env python3
"""Build benchmark pools from the mono-allelic MS dataset.

Each `cell_context` is `<study>-<instrument>-<allele>`. The label is the allele,
so the 113 contexts collapse to 91 molecules; the 20 alleles that appear in two
or three contexts are the same protein sampled on different instruments, and
their peptide sets overlap by a median Jaccard of only 0.19, so merging them
adds independent observations rather than duplicates.

Studies are mixed freely. That leaves a batch effect in the pools - mean peptide
length runs from 9.21 to 9.99 across the six studies - which is recorded per
pool (`study_count`, `length_mean`) rather than designed away.

A peptide observed on several alleles is kept once and assigned to one of the
alleles present in the pool, drawn uniformly at random. The metrics need one
true label per peptide, and dropping such peptides instead would cost 30% of the
data at 91 alleles while costing almost nothing at 3, which would confound
ambiguity with allele count. The fraction affected is recorded per pool as
`reassigned_fraction`, so the label noise it injects stays visible.

Pool size is a target, not a guarantee: per-allele draws are capped at what each
allele actually has, the shortfall is redistributed to alleles with headroom,
and the realised size is what the manifest reports.
"""

from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

# Equal shares, then four Dirichlet concentrations from mild to extreme skew.
# Lower alpha puts more of the pool into fewer alleles.
BALANCE = {"eq": None, "a10": 10.0, "a3": 3.0, "a1": 1.0, "a03": 0.3}

# Target sizes per allele count, chosen to stay inside what random draws of that
# many alleles can actually supply (see REPORT.md). The largest cell is 50
# alleles rather than all 91: MixMHCp will not fit more than 50 motifs, so a
# 91-allele pool could not be given its true count and the comparison would
# stop being like-for-like at exactly the point it gets interesting.
GRID = {
    3:  [1000, 2000, 4000, 8000],
    6:  [1000, 2000, 4000, 8000],
    12: [2000, 4000, 8000, 16000],
    24: [4000, 8000, 16000, 32000],
    50: [16000, 32000, 64000, 100000],
}
REPLICATES = {3: 2, 6: 2, 12: 2, 24: 1, 50: 1}

# Every allele in a pool gets at least this many peptides. Without it the
# skewed regimes hand an allele one or two peptides, which cannot form a motif
# and only adds an unclusterable singleton, and they drop alleles entirely:
# a nominal 91-allele pool came out with 81. The floor is well under the
# smallest allele in the dataset (1,101 peptides), so it never binds on supply.
MIN_PER_ALLELE = 20


def pretty_allele(token: str) -> str:
    """A0201 -> HLA-A02:01, the spelling the earlier benchmark used."""
    return f"HLA-{token[0]}{token[1:3]}:{token[3:]}"


def allocate(counts: np.ndarray, target: int, alpha: float | None,
             rng: np.random.Generator) -> np.ndarray:
    """How many peptides to draw from each allele.

    Caps every allele at what it has and pushes the shortfall onto alleles with
    headroom, so a skewed draw that asks a small allele for more than it owns
    still fills the pool instead of silently shrinking it.
    """
    n = len(counts)
    share = np.full(n, 1.0 / n) if alpha is None else rng.dirichlet(np.full(n, alpha))
    # Seat every allele first, then let the draw distribute what is left. The
    # floor is what keeps a nominal N-allele pool an N-allele pool.
    floor = np.minimum(counts, MIN_PER_ALLELE).astype(np.int64)
    if floor.sum() >= target:                  # target too small to seat everyone
        floor = np.minimum(floor, max(1, target // n))
    remaining = max(0, target - int(floor.sum()))
    want = floor + np.floor(share * remaining).astype(np.int64)
    for _ in range(64):
        want = np.minimum(want, counts)
        deficit = target - want.sum()
        if deficit <= 0:
            break
        headroom = counts - want
        if headroom.sum() == 0:
            break
        # Spread the shortfall in proportion to the room each allele has left.
        add = np.floor(deficit * headroom / headroom.sum()).astype(np.int64)
        if add.sum() == 0:                     # rounding stalled: place one by one
            order = np.argsort(-headroom)
            for i in order[:deficit]:
                if headroom[i] > 0:
                    want[i] += 1
            continue
        want = want + add
    return np.minimum(want, counts).astype(int)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261118)
    args = parser.parse_args()
    (args.out / "pools").mkdir(parents=True, exist_ok=True)

    frame = pd.read_parquet(args.data)
    frame["study"] = frame.cell_context.str.rsplit("-", n=1).str[0]
    frame["allele"] = frame.cell_context.str.rsplit("-", n=1).str[1]
    # One row per (peptide, allele): the same peptide seen on one allele in two
    # contexts is one observation of that molecule, not two.
    frame = frame.drop_duplicates(["peptide", "allele"])

    by_allele = {a: g for a, g in frame.groupby("allele")}
    alleles = np.array(sorted(by_allele))
    counts_all = pd.Series({a: len(g) for a, g in by_allele.items()})
    print(f"{len(frame):,} (peptide, allele) pairs over {len(alleles)} alleles")

    rows = []
    for n_alleles in sorted(GRID):
        for target in GRID[n_alleles]:
            for balance, alpha in BALANCE.items():
                for rep in range(REPLICATES[n_alleles]):
                    name = f"a{n_alleles:02d}_s{target:06d}_{balance}_r{rep}"
                    # One seed per pool, derived from its name, so any single
                    # pool can be rebuilt without rebuilding the rest. crc32,
                    # not hash(): str.__hash__ is salted per process, so hash()
                    # would give a different dataset on every run.
                    seed = (args.seed + zlib.crc32(name.encode())) % 2**32
                    rng = np.random.default_rng(seed)

                    picked = (alleles if n_alleles >= len(alleles)
                              else rng.choice(alleles, n_alleles, replace=False))
                    picked = np.array(sorted(picked))
                    want = allocate(counts_all[picked].to_numpy(), target, alpha, rng)

                    parts = []
                    for allele, k in zip(picked, want):
                        if k <= 0:
                            continue
                        pool_of = by_allele[allele]
                        take = rng.choice(len(pool_of), int(k), replace=False)
                        parts.append(pool_of.iloc[take][["peptide", "study"]]
                                     .assign(allele=allele))
                    drawn = pd.concat(parts, ignore_index=True)

                    # A peptide drawn for two alleles keeps one label, chosen at
                    # random, so the pool stays a partition.
                    duplicated = drawn.peptide.duplicated(keep=False)
                    n_ambiguous = drawn.loc[duplicated, "peptide"].nunique()
                    drawn = (drawn.sample(frac=1.0, random_state=int(seed % 2**31))
                                  .drop_duplicates("peptide")
                                  .sort_values("peptide")
                                  .reset_index(drop=True))

                    fasta = args.out / "pools" / f"{name}.fasta"
                    fasta.write_text("".join(
                        f">p{i}\n{p}\n" for i, p in enumerate(drawn.peptide)))
                    labels = drawn[["peptide", "allele"]].copy()
                    labels["allele"] = labels.allele.map(pretty_allele)
                    labels.to_csv(args.out / "pools" / f"{name}.labels.tsv",
                                  sep="\t", index=False)

                    per = drawn.groupby("allele").size()
                    lengths = drawn.peptide.str.len()
                    rows.append(dict(
                        pool=name, split="test", outer_fold=-1,
                        allele_count=int(per.size), peptides=int(len(drawn)),
                        target_peptides=target, balance=balance,
                        alleles=";".join(sorted(pretty_allele(a) for a in per.index)),
                        locus_A=int(sum(a.startswith("A") for a in per.index)),
                        locus_B=int(sum(a.startswith("B") for a in per.index)),
                        locus_C=int(sum(a.startswith("C") for a in per.index)),
                        smallest_allele_peptides=int(per.min()),
                        largest_allele_peptides=int(per.max()),
                        imbalance_ratio=round(per.max() / per.min(), 4),
                        max_allele_share=round(per.max() / len(drawn), 4),
                        reassigned_fraction=round(n_ambiguous / len(drawn), 4),
                        study_count=int(drawn.study.nunique()),
                        length_mean=round(float(lengths.mean()), 3),
                        pct_9mer=round(float(100 * (lengths == 9).mean()), 2),
                        seed=int(seed)))
                    print(f"  {name:28s} {len(drawn):7,d} peptides  "
                          f"{per.size:3d} alleles  imbalance {rows[-1]['imbalance_ratio']:8.2f}  "
                          f"reassigned {100*rows[-1]['reassigned_fraction']:5.2f}%")

    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.out / "pool_manifest.csv", index=False)
    print(f"\n{len(manifest)} pools, {manifest.peptides.sum():,} peptides in total")
    print(f"manifest: {args.out / 'pool_manifest.csv'}")


if __name__ == "__main__":
    main()
