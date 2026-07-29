#!/usr/bin/env python3
"""Run the PepCluster2 speed/peak-memory benchmark on a compute cluster.

The runner is resumable: a successful status.json is treated as cached. Failed
runs remain in the CSV and may be retried with --retry-failed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import matplotlib.pyplot as plt


SIZES = (1_000, 10_000, 20_000, 50_000, 100_000, 500_000, 1_000_000)
METHODS = {
    "graph": ["--clustering-method", "graph", "--no-prefilter"],
    "graph_prefilter": ["--clustering-method", "graph", "--force-prefilter"],
    "greedy": [
        "--clustering-method", "greedy", "--greedy-selection", "kmer-degree",
        "--no-prefilter",
    ],
    "greedy_lazy": [
        "--clustering-method", "greedy", "--greedy-selection", "lazy-exact",
        "--no-prefilter",
    ],
}
LABELS = {
    "graph": "Graph",
    "graph_prefilter": "Graph + prefilter",
    "greedy": "Greedy",
    "greedy_lazy": "Greedy lazy-exact",
}


def parse_peak_kb(path: Path) -> int | None:
    if not path.exists():
        return None
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", path.read_text())
    return int(match.group(1)) if match else None


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as reader, target.open("wb") as writer:
        shutil.copyfileobj(reader, writer, 8 * 1024 * 1024)


def run_one(args: argparse.Namespace, method: str, size: int) -> dict:
    output = args.output_dir / "runs" / method / f"n_{size:07d}"
    output.mkdir(parents=True, exist_ok=True)
    status_path = output / "status.json"
    if status_path.exists() and not args.retry_failed:
        return json.loads(status_path.read_text())
    if status_path.exists() and json.loads(status_path.read_text()).get("status") == "ok":
        return json.loads(status_path.read_text())

    source = args.data_dir / f"benchmark_{size:07d}.fasta.gz"
    if not source.exists():
        raise FileNotFoundError(source)
    run_tmp_parent = args.tmp_root / method / f"n_{size:07d}"
    run_tmp_parent.mkdir(parents=True, exist_ok=True)
    resource = output / "resource.txt"
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="input_", dir=run_tmp_parent) as temporary:
        fasta = Path(temporary) / "input.fasta"
        decompress(source, fasta)
        command = [
            str(args.binary),
            "--input", str(fasta),
            "--output-dir", str(output),
            "--tmp-dir", str(run_tmp_parent / "pepcluster_tmp"),
            "--mode", "separate_aln_anchor",
            "--alignment-similarity-threshold", "0.50",
            "--anchor-combination-similarity-threshold", "0.60",
            "--kmer-seed-threshold", "0.50",
            "--gap-open", "-4",
            "--gap-extension", "-1",
            "--terminal-overhang-gap-open", "-2",
            "--terminal-overhang-gap-extension", "-1",
            "--minimum-terminal-match-length", "2",
            "--threads", str(args.threads),
            "--candidate-buffer-mb", str(args.candidate_buffer_mb),
            "--compact-output",
            *METHODS[method],
        ]
        if args.max_memory_gb is not None:
            command += ["--max-memory-gb", str(args.max_memory_gb)]
        (output / "benchmark_command.txt").write_text(" ".join(command) + "\n")
        timed = ["/usr/bin/time", "-v", "-o", str(resource), *command]
        with (output / "benchmark.log").open("w") as log:
            result = subprocess.run(timed, stdout=log, stderr=subprocess.STDOUT)

    elapsed = time.time() - started
    peak_kb = parse_peak_kb(resource)
    status = {
        "method": method,
        "records": size,
        "status": "ok" if result.returncode == 0 else "failed",
        "returncode": result.returncode,
        "wall_seconds": elapsed,
        "peak_memory_mb": peak_kb / 1024 if peak_kb is not None else None,
        "output_dir": str(output),
    }
    if (output / "run_stats.json").exists():
        stats = json.loads((output / "run_stats.json").read_text())
        status["pepcluster_elapsed_seconds"] = stats.get("elapsed_seconds")
        status["final_clusters"] = stats.get("final_clusters")
        status["singleton_clusters"] = stats.get("singleton_clusters")
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    return status


def collect(args: argparse.Namespace, methods: list[str], sizes: list[int]) -> list[dict]:
    rows = []
    for method in methods:
        for size in sizes:
            status = args.output_dir / "runs" / method / f"n_{size:07d}" / "status.json"
            if status.exists():
                rows.append(json.loads(status.read_text()))
    figures = args.output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    csv_path = figures / "resource_benchmark.csv"
    fields = [
        "method", "records", "status", "returncode", "wall_seconds",
        "pepcluster_elapsed_seconds", "peak_memory_mb", "final_clusters",
        "singleton_clusters", "output_dir",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot(args: argparse.Namespace, rows: list[dict], methods: list[str]) -> None:
    successful = [row for row in rows if row["status"] == "ok"]
    if not successful:
        return
    colors = ["#31688e", "#35b779", "#e07a00", "#7e3ace"]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for method, color in zip(methods, colors):
        selected = sorted(
            (row for row in successful if row["method"] == method),
            key=lambda row: row["records"],
        )
        if not selected:
            continue
        x = [row["records"] for row in selected]
        axes[0].plot(x, [row["wall_seconds"] for row in selected], "o-", label=LABELS[method], color=color)
        axes[1].plot(x, [row["peak_memory_mb"] for row in selected], "o-", label=LABELS[method], color=color)
    for axis, ylabel in zip(axes, ["Wall time (seconds)", "Peak memory (MiB)"]):
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Number of peptides")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False)
    figure.savefig(args.output_dir / "figures" / "resource_benchmark.png", dpi=220)
    figure.savefig(args.output_dir / "figures" / "resource_benchmark.pdf")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--sizes", nargs="+", type=int, choices=SIZES, default=list(SIZES))
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--candidate-buffer-mb", type=int, default=512)
    parser.add_argument("--max-memory-gb", type=float)
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    args.binary = args.binary.resolve()
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.tmp_root = args.tmp_root.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    if not args.plot_only:
        for method in args.methods:
            for size in args.sizes:
                row = run_one(args, method, size)
                print(json.dumps(row), flush=True)
    rows = collect(args, args.methods, args.sizes)
    plot(args, rows, args.methods)


if __name__ == "__main__":
    main()

