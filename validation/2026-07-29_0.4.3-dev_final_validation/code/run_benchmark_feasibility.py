#!/usr/bin/env python3
"""Run index-only resource/disk feasibility checks for benchmark sizes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


SIZES = (1_000, 10_000, 20_000, 50_000, 100_000, 500_000, 1_000_000)


def peak_memory_kb(path: Path) -> int:
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", path.read_text())
    return int(match.group(1)) if match else 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args()
    rows = []
    for size in SIZES:
        output = args.output_dir / f"n_{size:07d}"
        output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="pc2_index_", dir="/tmp") as tmp:
            fasta = Path(tmp) / "input.fasta"
            with gzip.open(args.data_dir / f"benchmark_{size:07d}.fasta.gz", "rb") as source:
                with fasta.open("wb") as target:
                    shutil.copyfileobj(source, target, 4 * 1024 * 1024)
            resource = output / "resource.txt"
            command = [
                "/usr/bin/time", "-v", "-o", str(resource),
                str(args.binary), "--input", str(fasta), "--output-dir", str(output),
                "--index-only", "--threads", str(args.threads),
                "--kmer-seed-threshold", "0.50",
            ]
            with (output / "run.log").open("w") as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"index feasibility failed at {size}; see {output / 'run.log'}")
        stats = json.loads((output / "run_stats.json").read_text())
        rows.append({
            "records": size,
            "candidate_occurrence_upper_bound": stats["candidate_occurrence_upper_bound"],
            "estimated_nonprefilter_disk_bytes": stats["estimated_non_prefilter_disk_bytes"],
            "estimated_nonprefilter_disk_gib": stats["estimated_non_prefilter_disk_bytes"] / 1024**3,
            "automatic_prefilter": stats["prefilter_active"],
            "peak_memory_mb": peak_memory_kb(resource) / 1024,
        })
    with (args.output_dir / "benchmark_feasibility.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
