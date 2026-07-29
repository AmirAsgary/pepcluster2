#!/usr/bin/env python3
"""Create one deterministic nested random benchmark sample up to one million peptides."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import random
from pathlib import Path


SIZES = (1_000, 10_000, 20_000, 50_000, 100_000, 500_000, 1_000_000)
CANONICAL = set("ARNDCQEGHILKMFPSTWYV")


def records(path: Path):
    header, sequence = "", []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            else:
                sequence.append(line.upper())
    if header:
        yield header, "".join(sequence)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pool", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260729)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    reservoir: list[tuple[str, str]] = []
    eligible = 0
    for header, sequence in records(args.pool):
        if len(sequence) < 8 or not set(sequence) <= CANONICAL:
            continue
        eligible += 1
        if len(reservoir) < SIZES[-1]:
            reservoir.append((header, sequence))
        else:
            position = rng.randrange(eligible)
            if position < SIZES[-1]:
                reservoir[position] = (header, sequence)
    if len(reservoir) != SIZES[-1]:
        raise RuntimeError(f"only {len(reservoir)} eligible records found")
    rng.shuffle(reservoir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for size in SIZES:
        path = args.output_dir / f"benchmark_{size:07d}.fasta.gz"
        digest = hashlib.sha256()
        with gzip.open(path, "wt", compresslevel=6) as handle:
            for index, (header, sequence) in enumerate(reservoir[:size]):
                text = f">benchmark_{index:07d} {header}\n{sequence}\n"
                handle.write(text)
                digest.update(text.encode())
        manifest.append([size, args.seed, digest.hexdigest(), path])
    with (args.output_dir / "manifest.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["records", "seed", "sha256_uncompressed", "path"])
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
