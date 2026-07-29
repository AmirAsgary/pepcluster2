#!/usr/bin/env python3
"""Create the final 20 x 10k datasets and deterministic nested subsets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import random
from pathlib import Path


SIZES = (1_000, 2_000, 4_000, 6_000, 8_000)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(sequence)))
                header, sequence = line[1:], []
            else:
                sequence.append(line)
    if header:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with gzip.open(path, "wt", compresslevel=6) as handle:
        for header, sequence in records:
            text = f">{header}\n{sequence}\n"
            handle.write(text)
            digest.update(text.encode())
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=20)
    parser.add_argument("--records", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260729)
    args = parser.parse_args()

    manifest: list[list[object]] = []
    for dataset in range(args.datasets):
        source = args.source_dir / f"sample_{dataset:03d}.fasta.gz"
        records = read_fasta(source)[: args.records]
        if len(records) != args.records:
            raise RuntimeError(f"{source} contains only {len(records)} records")
        full_path = args.output_dir / "full" / f"sample_{dataset:03d}.fasta.gz"
        digest = write_fasta(full_path, records)
        manifest.append(["full", dataset, args.records, args.seed + dataset, digest, full_path])

        order = list(range(args.records))
        random.Random(args.seed + dataset).shuffle(order)
        for size in SIZES:
            subset_records = [records[index] for index in order[:size]]
            subset_path = (
                args.output_dir / "subsets" / f"n_{size:06d}" / f"sample_{dataset:03d}.fasta.gz"
            )
            subset_digest = write_fasta(subset_path, subset_records)
            manifest.append(
                ["subset", dataset, size, args.seed + dataset, subset_digest, subset_path]
            )

    with (args.output_dir / "manifest.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["kind", "dataset", "records", "sampling_seed", "sha256_uncompressed", "path"])
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
