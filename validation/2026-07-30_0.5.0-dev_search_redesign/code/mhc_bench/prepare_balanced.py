#!/usr/bin/env python3
"""Balanced, similarity-filtered peptide-MHC pools (second benchmark design).

Two changes from the first design, both aimed at removing a confounder rather
than making the task easier for its own sake:

Allele balance
    Every allele contributes the same share of a pool, 100/k, jittered by at most
    5% relative. In the first design one allele could supply an order of
    magnitude more peptides than another, which made per-allele averages depend
    on which alleles happened to be abundant.

Motif redundancy
    Alleles more than `--similarity-threshold` percent identical are never placed
    in the same benchmark, selected as a maximum independent set over the
    similarity graph. Near-identical alleles present near-identical motifs, so no
    sequence-based method can separate them; scoring against them measures the
    label, not the method. This makes the benchmark easier than a realistic pool
    and the resulting numbers must be reported as such.

Pool size is capped by balance: a pool of k alleles cannot exceed k times the
smallest contributing allele, so the reachable maximum is well below the nominal
ceiling and small-k pools are necessarily the small ones.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare import (LOCI, MIN_POOL_PEPTIDES, OUTER_FOLDS, POOLS_PER_COUNT,  # noqa: E402
                     load_table, locus_of, stratified_folds, write_pool)

INNER_ALLELE_COUNTS = (2, 5, 10, 20, 30)
MAX_POOL_PEPTIDES = 50_000
BALANCE_TOLERANCE = 0.05


def compatible_alleles(similarity: Path, alleles: list[str], threshold: float,
                       seed: int) -> list[str]:
    """Largest set of alleles that are pairwise below the similarity threshold.

    Maximum independent set is NP-hard, so this is randomised greedy with many
    restarts, taking lowest-degree vertices first.
    """
    table = pd.read_csv(similarity)
    table = table[table.hla_1.isin(alleles) & table.hla_2.isin(alleles)]
    table = table[table.hla_1 != table.hla_2]
    conflicts: dict[str, set] = defaultdict(set)
    for row in table[table.similarity_pct >= threshold].itertuples():
        conflicts[row.hla_1].add(row.hla_2)
        conflicts[row.hla_2].add(row.hla_1)
    rng = np.random.default_rng(seed)
    best: list[str] = []
    for _ in range(4000):
        order = sorted(alleles, key=lambda a: (len(conflicts.get(a, ())), rng.random()))
        chosen, blocked = [], set()
        for allele in order:
            if allele not in blocked:
                chosen.append(allele)
                blocked |= conflicts.get(allele, set())
                blocked.add(allele)
        if len(chosen) > len(best):
            best = chosen
    return sorted(best)


def locus_quota(count: int, index: int) -> dict[str, int]:
    if count == 2:
        pairs = [("A", "B"), ("A", "C"), ("B", "C")]
        first, second = pairs[index % len(pairs)]
        return {l: (1 if l in (first, second) else 0) for l in LOCI}
    base, remainder = divmod(count, len(LOCI))
    quota = {locus: base for locus in LOCI}
    for offset in range(remainder):
        quota[LOCI[(index + offset) % len(LOCI)]] += 1
    return quota


def choose_alleles(available: list[str], count: int, index: int, rng,
                   supply: dict[str, int] | None = None,
                   size_level: int = 0, levels: int = 1) -> list[str]:
    """Locus-balanced sample, biased toward abundant alleles at large size levels.

    Balance caps a pool at k times its smallest allele, so a uniformly random
    choice makes large pools unreachable: one small allele drags the whole pool
    down. Higher size levels therefore draw from a narrowing window of the most
    abundant alleles within each locus, while level 0 still samples the full
    range so small pools stay representative.
    """
    by_locus = {l: sorted(a for a in available if locus_of(a) == l) for l in LOCI}
    if supply is not None and levels > 1:
        fraction = 1.0 - 0.75 * (size_level / (levels - 1))
        for locus, members in by_locus.items():
            ranked = sorted(members, key=lambda a: -supply.get(a, 0))
            by_locus[locus] = ranked[:max(1, int(round(len(ranked) * fraction)))]
    quota = locus_quota(count, index)
    chosen: list[str] = []
    for locus in LOCI:
        pool = by_locus[locus]
        take = min(quota[locus], len(pool))
        if take:
            chosen.extend(rng.choice(pool, take, replace=False).tolist())
    if len(chosen) < count:
        rest = sorted(set(available) - set(chosen))
        take = min(count - len(chosen), len(rest))
        if take:
            chosen.extend(rng.choice(rest, take, replace=False).tolist())
    return sorted(chosen)


def exclusive_peptides(alleles, peptide_alleles, allele_peptides) -> dict[str, list[str]]:
    """Peptides binding exactly one allele of this pool, grouped by that allele."""
    selected = set(alleles)
    result = {}
    for allele in alleles:
        result[allele] = [p for p in allele_peptides[allele]
                          if len(peptide_alleles[p] & selected) == 1]
    return result


def build_balanced(alleles, peptide_alleles, allele_peptides, target, rng):
    exclusive = exclusive_peptides(alleles, peptide_alleles, allele_peptides)
    smallest = min(len(v) for v in exclusive.values())
    if smallest == 0:
        return None
    # Balance caps the pool: no allele may be asked for more than it has.
    share = min(target // len(alleles), smallest)
    if share * len(alleles) < MIN_POOL_PEPTIDES:
        return None
    peptides, labels = [], []
    for allele in alleles:
        jitter = rng.uniform(1.0 - BALANCE_TOLERANCE, 1.0 + BALANCE_TOLERANCE)
        want = min(int(round(share * jitter)), len(exclusive[allele]))
        take = rng.choice(len(exclusive[allele]), want, replace=False)
        for index in sorted(take):
            peptides.append(exclusive[allele][index])
            labels.append(allele)
    order = np.argsort(peptides, kind="stable")
    return [peptides[i] for i in order], [labels[i] for i in order]


def make_pools(spec, peptide_alleles, allele_peptides, directory) -> list[dict]:
    rows = []
    for entry in spec:
        rng = np.random.default_rng(entry["seed"])
        count = min(entry["allele_count"], len(entry["available"]))
        alleles = choose_alleles(entry["available"], count, entry["index"], rng,
                                 supply={a: len(v) for a, v in allele_peptides.items()},
                                 size_level=entry["size_level"], levels=entry["levels"])
        exclusive = exclusive_peptides(alleles, peptide_alleles, allele_peptides)
        ceiling = min(MAX_POOL_PEPTIDES, len(alleles) * min(len(v) for v in exclusive.values()))
        if ceiling < MIN_POOL_PEPTIDES:
            continue
        grid = np.exp(np.linspace(np.log(MIN_POOL_PEPTIDES), np.log(ceiling), entry["levels"]))
        target = int(round(grid[entry["size_level"] % entry["levels"]]))
        built = build_balanced(alleles, peptide_alleles, allele_peptides, target, rng)
        if built is None:
            continue
        peptides, labels = built
        write_pool(directory, entry["name"], peptides, labels)
        sizes = pd.Series(labels).value_counts()
        rows.append({
            "pool": entry["name"], "split": entry["split"], "outer_fold": entry["outer_fold"],
            "allele_count": len(alleles), "peptides": len(peptides),
            "alleles": ";".join(alleles),
            "locus_A": sum(1 for a in alleles if locus_of(a) == "A"),
            "locus_B": sum(1 for a in alleles if locus_of(a) == "B"),
            "locus_C": sum(1 for a in alleles if locus_of(a) == "C"),
            "smallest_allele_peptides": int(sizes.min()),
            "largest_allele_peptides": int(sizes.max()),
            "imbalance_ratio": float(sizes.max() / sizes.min()),
            "max_allele_share": float(sizes.max() / len(peptides)),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--similarity", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--similarity-threshold", type=float, default=97.0)
    parser.add_argument("--discovery-min-peptides", type=int, default=500)
    parser.add_argument("--test-min-peptides", type=int, default=100)
    parser.add_argument("--test-pools", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260731)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    discovery = load_table(args.discovery, args.discovery_min_peptides)
    every_allele = sorted(discovery["allele"].unique())
    alleles = compatible_alleles(args.similarity, every_allele,
                                 args.similarity_threshold, args.seed)
    discovery = discovery[discovery["allele"].isin(alleles)]
    peptide_alleles = {p: frozenset(g) for p, g in discovery.groupby("peptide")["allele"]}
    allele_peptides = {a: sorted(g) for a, g in discovery.groupby("allele")["peptide"]}
    folds = stratified_folds(alleles, OUTER_FOLDS, args.seed)

    manifest: list[dict] = []
    for fold_index, held_out in enumerate(folds):
        training = sorted(set(alleles) - set(held_out))
        spec = []
        for count in INNER_ALLELE_COUNTS:
            for repeat in range(POOLS_PER_COUNT):
                spec.append({
                    "name": f"inner_f{fold_index}_a{count:02d}_{repeat}",
                    "split": "inner", "outer_fold": fold_index,
                    "allele_count": count, "available": training, "index": repeat,
                    "seed": args.seed + 1000 * fold_index + 17 * count + repeat,
                    "size_level": repeat, "levels": POOLS_PER_COUNT,
                })
        outer_counts = [2, 5, 10, min(20, len(held_out)), len(held_out)]
        for repeat, count in enumerate(outer_counts):
            spec.append({
                "name": f"outer_f{fold_index}_a{count:02d}_{repeat}",
                "split": "outer", "outer_fold": fold_index,
                "allele_count": count, "available": held_out, "index": repeat,
                "seed": args.seed + 500_000 + 1000 * fold_index + repeat,
                "size_level": repeat, "levels": len(outer_counts),
            })
        manifest += make_pools(spec, peptide_alleles, allele_peptides, output / "pools")

    test = load_table(args.test, args.test_min_peptides)
    test_alleles = compatible_alleles(args.similarity, sorted(test["allele"].unique()),
                                      args.similarity_threshold, args.seed)
    test = test[test["allele"].isin(test_alleles)]
    test_peptide_alleles = {p: frozenset(g) for p, g in test.groupby("peptide")["allele"]}
    test_allele_peptides = {a: sorted(g) for a, g in test.groupby("allele")["peptide"]}
    rng = np.random.default_rng(args.seed + 99)
    spec = []
    for index in range(args.test_pools):
        count = int(rng.integers(2, min(len(test_alleles), 20) + 1))
        spec.append({
            "name": f"test_{index:02d}", "split": "test", "outer_fold": -1,
            "allele_count": count, "available": test_alleles, "index": index,
            "seed": args.seed + 900_000 + index, "size_level": index, "levels": 10,
        })
    manifest += make_pools(spec, test_peptide_alleles, test_allele_peptides, output / "pools")

    with (output / "pool_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    summary = {
        "similarity_threshold": args.similarity_threshold,
        "alleles_before_similarity_filter": len(every_allele),
        "alleles_after": len(alleles),
        "locus_counts": {l: sum(1 for a in alleles if locus_of(a) == l) for l in LOCI},
        "test_alleles_after": len(test_alleles),
        "pools": len(manifest),
        "max_allele_share_observed": max(r["max_allele_share"] for r in manifest),
        "worst_imbalance_ratio": max(r["imbalance_ratio"] for r in manifest),
    }
    (output / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
