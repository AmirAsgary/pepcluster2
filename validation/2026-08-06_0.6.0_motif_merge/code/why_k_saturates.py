#!/usr/bin/env python3
"""Why does the motif count stop rising at about ten?

Our component count tracks the allele count up to roughly ten alleles and then
flattens and slightly falls, while pools go to twenty. Two explanations, and they
call for different fixes:

  DEGENERACY   Some alleles genuinely share a binding motif - supertype members
               differ in the groove but not in what the groove selects for. Those
               alleles are not separable from peptide sequence by any method, so
               a partition with fewer motifs than alleles is correct and the
               allele label is simply the wrong target.

  UNDER-SPLIT  Our EM smoothing collapses components that the data would support.
               The concentration is a single global constant, so the number of
               components it sustains is roughly fixed regardless of how many
               motifs are present. That is our problem, and it is fixable.

Three measurements separate them:

  1. How distinguishable the TRUE per-allele profiles are. Built from the labels,
     so it is an oracle, and it bounds what any sequence-based method could
     recover.
  2. Whether forcing the true component count improves agreement. If it does, our
     model selection is leaving something on the table.
  3. What the oracle-merged ceiling looks like at each allele count, which
     bounds the similarity clustering's own resolution.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
STUDY = HERE.parents[2] / "2026-07-30_0.5.0-dev_search_redesign"
sys.path.insert(0, str(STUDY / "code" / "mhc_bench"))
import metrics as M  # noqa: E402
sys.path.insert(0, str(HERE.parent))
from motif_lib import (ALPHABET, NP, background, count_matrix, em_from_counts,  # noqa: E402
                       encode)


def profiles_by_allele(X, alleles, pseudo):
    names, index = np.unique(alleles, return_inverse=True)
    counts = count_matrix(X, index, len(names))
    probs = counts + pseudo
    probs /= probs.sum(-1, keepdims=True)
    return names, probs, np.bincount(index, minlength=len(names))


def js_divergence(p, q):
    """Jensen-Shannon divergence in bits, averaged over the nine columns.

    Symmetric, bounded in [0, 1] for base-2 logs, and finite even when a residue
    is absent from one profile - which a KL divergence would not be.
    """
    m = 0.5 * (p + q)
    def kl(a, b):
        mask = a > 0
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))
    return np.mean([0.5 * kl(p[j], m[j]) + 0.5 * kl(q[j], m[j]) for j in range(NP)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pools", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--assign", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--concentration", type=float, default=3.0)
    parser.add_argument("--distinct-threshold", type=float, default=0.15,
                        help="JS divergence in bits below which two allele "
                             "profiles are treated as the same motif")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest.split == "test"]
    rows = []
    for record in manifest.itertuples():
        path = args.assign / f"{record.pool}.tsv"
        if not path.exists():
            continue
        assign = pd.read_csv(path, sep="\t")
        alleles = assign.allele.to_numpy()
        X = encode(assign.peptide.values)
        bg = background(X)
        pseudo = args.concentration * bg * ALPHABET

        names, probs, sizes = profiles_by_allele(X, alleles, pseudo)
        k_true = len(names)

        # 1. how many of the true allele motifs are actually distinguishable
        pairs = [(i, j) for i in range(k_true) for j in range(i + 1, k_true)]
        divergences = np.array([js_divergence(probs[i], probs[j]) for i, j in pairs])
        # Complete linkage, not single. Single linkage chains: with A close to B
        # and B close to C, everything collapses into one group even when A and C
        # are far apart, which on this data collapsed 11 alleles to 1.6 and was
        # plainly wrong. Complete linkage merges only when EVERY pair in the
        # merged group is below threshold.
        #
        # There is no natural threshold to pick: the divergences form a smooth
        # continuum from 0.06 to 0.24 bits with no gap, so any cut is arbitrary.
        # The pair-level fractions below are the honest summary; k_distinct is
        # reported only as a coarse indication at one stated cut.
        import scipy.cluster.hierarchy as sch
        from scipy.spatial.distance import squareform
        if k_true > 1:
            condensed = np.array(divergences)
            linkage = sch.linkage(condensed, method="complete")
            k_distinct = int(sch.fcluster(linkage, args.distinct_threshold,
                                          criterion="distance").max())
        else:
            k_distinct = 1
        frac_close_015 = float(np.mean(divergences < 0.15)) if len(divergences) else np.nan
        frac_close_010 = float(np.mean(divergences < 0.10)) if len(divergences) else np.nan
        min_js = float(np.min(divergences)) if len(divergences) else np.nan

        # 2. does forcing the true component count help?
        #    EM seeded from the true allele partition would be cheating on the
        #    assignment, so components are seeded from a k-means-free split of the
        #    similarity clusters: take the k_true largest clusters as seeds.
        ids, inv = np.unique(assign.cluster_id.values, return_inverse=True)
        counts = count_matrix(X, inv, len(ids))
        order = np.argsort(-np.bincount(inv, minlength=len(ids)))
        for k in (k_true,):
            seed_counts = counts[order[:k]].copy()
            weights = np.full(k, 1.0 / k)
            labels, _ = em_from_counts(X, seed_counts, weights, pseudo)
            forced = M.evaluate(alleles, labels)
        rows.append(dict(
            pool=record.pool, peptides=record.peptides, alleles=k_true,
            distinct_motifs=k_distinct,
            degenerate_fraction=1 - k_distinct / k_true,
            median_pairwise_js=float(np.median(divergences)) if len(divergences) else np.nan,
            min_pairwise_js=min_js, frac_pairs_below_0_15=frac_close_015,
            frac_pairs_below_0_10=frac_close_010,
            forced_k=int(forced["clusters"]), forced_k_ami=forced["ami"],
            forced_k_f1=forced["bcubed_f1_macro"]))
        print(".", end="", flush=True)
    print()
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out / "k_saturation_analysis.csv", index=False)

    frame["band"] = pd.cut(frame.alleles, bins=[3, 6, 9, 12, 16, 20])
    table = frame.groupby("band", observed=True).agg(
        pools=("pool", "size"), true_alleles=("alleles", "mean"),
        distinct_motifs=("distinct_motifs", "mean"),
        degenerate_fraction=("degenerate_fraction", "mean"),
        forced_k_ami=("forced_k_ami", "mean"),
        forced_k_f1=("forced_k_f1", "mean"))
    print("\n=== how many motifs are actually distinguishable? ===")
    print(f"(alleles collapsed when their true profiles differ by less than "
          f"{args.distinct_threshold} bits Jensen-Shannon)")
    print(table.round(3).to_string())
    print(f"\noverall: {frame.alleles.mean():.2f} alleles -> "
          f"{frame.distinct_motifs.mean():.2f} distinguishable motifs "
          f"({100 * frame.degenerate_fraction.mean():.1f}% degenerate)")
    print(f"forcing the true component count: AMI {frame.forced_k_ami.mean():.4f}, "
          f"F1 {frame.forced_k_f1.mean():.4f}")


if __name__ == "__main__":
    main()
