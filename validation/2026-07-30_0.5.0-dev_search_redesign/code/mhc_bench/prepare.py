#!/usr/bin/env python3
"""Build the peptide-MHC benchmark pools for hyperparameter optimisation.

Ground truth
------------
A peptide may bind several alleles, but clustering assigns it to exactly one
cluster, so a pool keeps a peptide only when it binds exactly one of the alleles
in that pool. A peptide promiscuous across the whole dataset is therefore still
usable wherever the pool contains only one of its alleles, and is never given an
arbitrary label.

Design
------
Alleles, not peptides, are split into the outer folds: the hyperparameters are
similarity thresholds, so the question is whether they transfer to alleles the
tuning never saw.

* 5 outer folds over alleles, stratified by locus (A/B/C).
* 25 inner pools per outer fold, drawn from that fold's training alleles:
  five each at 2, 5, 10, 20 and 40 alleles, locus-balanced.
* 5 outer pools per fold from the held-out alleles.
* Pool size is drawn log-uniformly so that dataset size is crossed with allele
  count rather than determined by it.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

VALID_PEPTIDE = re.compile(r"[ARNDCQEGHILKMFPSTWYV]{8,}")
LOCI = ("A", "B", "C")
INNER_ALLELE_COUNTS = (2, 5, 10, 20, 40)
POOLS_PER_COUNT = 5
OUTER_FOLDS = 5
MIN_POOL_PEPTIDES = 1_000
MAX_POOL_PEPTIDES = 100_000


def load_table(path: Path, minimum_peptides: int) -> pd.DataFrame:
    table = pd.read_excel(path)
    table = table[table["MHC class"].astype(str).str.strip() == "I"].copy()
    table["allele"] = table["HLA allele"].astype(str).str.strip()
    table["peptide"] = table["Peptide"].astype(str).str.strip().str.upper()
    table = table[table["peptide"].str.fullmatch(VALID_PEPTIDE)]
    table = table.drop_duplicates(["allele", "peptide"])
    counts = table.groupby("allele")["peptide"].nunique()
    keep = counts[counts >= minimum_peptides].index
    return table[table["allele"].isin(keep)][["allele", "peptide"]].reset_index(drop=True)


def locus_of(allele: str) -> str:
    match = re.search(r"HLA-([ABC])", allele)
    return match.group(1) if match else "?"


def stratified_folds(alleles: list[str], folds: int, seed: int) -> list[list[str]]:
    """Round-robin alleles of each locus across folds, so every fold keeps the
    A/B/C ratio of the source data."""
    rng = np.random.default_rng(seed)
    assignment: list[list[str]] = [[] for _ in range(folds)]
    for locus in LOCI:
        members = sorted(a for a in alleles if locus_of(a) == locus)
        rng.shuffle(members)
        for index, allele in enumerate(members):
            assignment[index % folds].append(allele)
    return [sorted(fold) for fold in assignment]


def locus_quota(count: int, index: int) -> dict[str, int]:
    """Split `count` alleles across A/B/C as evenly as possible. Two-allele pools
    cannot cover three loci, so they rotate through the locus pairs."""
    if count == 2:
        pairs = [("A", "B"), ("A", "C"), ("B", "C")]
        first, second = pairs[index % len(pairs)]
        return {first: 1, second: 1, **{l: 0 for l in LOCI if l not in (first, second)}}
    base, remainder = divmod(count, len(LOCI))
    quota = {locus: base for locus in LOCI}
    for offset in range(remainder):
        quota[LOCI[(index + offset) % len(LOCI)]] += 1
    return quota


def choose_alleles(available: list[str], count: int, index: int, rng) -> list[str]:
    """Locus-balanced sample. Any shortfall in one locus is refilled from the
    others so the requested pool size is still met."""
    by_locus = {locus: sorted(a for a in available if locus_of(a) == locus) for locus in LOCI}
    quota = locus_quota(count, index)
    chosen: list[str] = []
    for locus in LOCI:
        pool = by_locus[locus]
        take = min(quota[locus], len(pool))
        if take:
            chosen.extend(rng.choice(pool, take, replace=False).tolist())
    if len(chosen) < count:
        rest = sorted(set(available) - set(chosen))
        extra = min(count - len(chosen), len(rest))
        chosen.extend(rng.choice(rest, extra, replace=False).tolist())
    return sorted(chosen)


def build_pool(
    alleles: list[str],
    peptide_alleles: dict[str, frozenset],
    allele_peptides: dict[str, list[str]],
    rng,
) -> tuple[list[str], list[str]]:
    """Peptides binding exactly one allele of this pool, with their label."""
    selected = set(alleles)
    peptides, labels = [], []
    for allele in alleles:
        for peptide in allele_peptides[allele]:
            if len(peptide_alleles[peptide] & selected) == 1:
                peptides.append(peptide)
                labels.append(allele)
    order = np.argsort(peptides, kind="stable")
    return [peptides[i] for i in order], [labels[i] for i in order]


def subsample(peptides, labels, target, rng):
    if len(peptides) <= target:
        return peptides, labels
    keep = np.sort(rng.choice(len(peptides), target, replace=False))
    return [peptides[i] for i in keep], [labels[i] for i in keep]


def write_pool(directory: Path, name: str, peptides: list[str], labels: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / f"{name}.fasta").open("w") as handle:
        for index, peptide in enumerate(peptides):
            handle.write(f">p{index}\n{peptide}\n")
    with (directory / f"{name}.labels.tsv").open("w") as handle:
        handle.write("peptide\tallele\n")
        for peptide, allele in zip(peptides, labels):
            handle.write(f"{peptide}\t{allele}\n")


def make_pools(spec, peptide_alleles, allele_peptides, directory, rng_seed) -> list[dict]:
    rows = []
    for entry in spec:
        rng = np.random.default_rng(entry["seed"] + rng_seed)
        available = entry["available"]
        count = min(entry["allele_count"], len(available))
        alleles = choose_alleles(available, count, entry["index"], rng)
        peptides, labels = build_pool(alleles, peptide_alleles, allele_peptides, rng)
        if len(peptides) < MIN_POOL_PEPTIDES:
            continue
        # Size is a designed factor, not a random draw: the repeats at a given
        # allele count are spread evenly on the log scale from the floor to what
        # that allele count can actually supply. Drawing at random left the large
        # sizes almost unsampled, because a draw rarely lands near the ceiling.
        ceiling = min(MAX_POOL_PEPTIDES, len(peptides))
        grid = np.exp(np.linspace(np.log(MIN_POOL_PEPTIDES), np.log(ceiling), entry["levels"]))
        target = int(round(grid[entry["size_level"] % entry["levels"]]))
        peptides, labels = subsample(peptides, labels, target, rng)
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
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--discovery", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--discovery-min-peptides", type=int, default=500)
    parser.add_argument("--test-min-peptides", type=int, default=100)
    parser.add_argument("--test-pools", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260730)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    discovery = load_table(args.discovery, args.discovery_min_peptides)
    peptide_alleles = {
        peptide: frozenset(group)
        for peptide, group in discovery.groupby("peptide")["allele"]
    }
    allele_peptides = {
        allele: sorted(group)
        for allele, group in discovery.groupby("allele")["peptide"]
    }
    alleles = sorted(allele_peptides)
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
                    "allele_count": count, "available": training,
                    "index": repeat, "seed": 1000 * fold_index + 17 * count + repeat,
                    "size_level": repeat, "levels": POOLS_PER_COUNT,
                })
        # Held-out folds hold about a fifth of the alleles, so the largest inner
        # size is not reproducible there; the last pool uses every held-out allele.
        outer_counts = [2, 5, 10, 20, len(held_out)]
        for repeat, count in enumerate(outer_counts):
            spec.append({
                "name": f"outer_f{fold_index}_a{count:02d}_{repeat}",
                "split": "outer", "outer_fold": fold_index,
                "allele_count": count, "available": held_out,
                "index": repeat, "seed": 500_000 + 1000 * fold_index + repeat,
                "size_level": repeat, "levels": len(outer_counts),
            })
        manifest += make_pools(spec, peptide_alleles, allele_peptides,
                               output / "pools", args.seed)

    test = load_table(args.test, args.test_min_peptides)
    test_peptide_alleles = {
        peptide: frozenset(group) for peptide, group in test.groupby("peptide")["allele"]
    }
    test_allele_peptides = {
        allele: sorted(group) for allele, group in test.groupby("allele")["peptide"]
    }
    test_alleles = sorted(test_allele_peptides)
    rng = np.random.default_rng(args.seed + 99)
    spec = []
    for index in range(args.test_pools):
        count = int(rng.integers(2, min(len(test_alleles), 20) + 1))
        spec.append({
            "name": f"test_{index:02d}", "split": "test", "outer_fold": -1,
            "allele_count": count, "available": test_alleles,
            "index": index, "seed": 900_000 + index,
            "size_level": index, "levels": 10,
        })
    manifest += make_pools(spec, test_peptide_alleles, test_allele_peptides,
                           output / "pools", args.seed)

    with (output / "pool_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)

    summary = {
        "discovery_alleles": len(alleles),
        "discovery_peptides": int(discovery["peptide"].nunique()),
        "discovery_locus_counts": {l: sum(1 for a in alleles if locus_of(a) == l) for l in LOCI},
        "test_alleles": len(test_alleles),
        "test_peptides": int(test["peptide"].nunique()),
        "outer_folds": [
            {"fold": i, "held_out_alleles": len(f),
             "locus": {l: sum(1 for a in f if locus_of(a) == l) for l in LOCI}}
            for i, f in enumerate(folds)
        ],
        "pools": len(manifest),
    }
    (output / "prepare_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
