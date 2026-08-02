#!/usr/bin/env python3
"""Speed, peak-memory and temporary-disk scaling of the redesigned candidate search.

The 10,000-peptide validation cannot show whether the extra index expansion is
affordable at scale, because candidate volume grows with the square of the
dataset while the sound anchor bound removes a roughly constant fraction. This
runs the historical and the redesigned configuration over the existing 1k-1M
benchmark datasets and reports the cost decomposition for each.

Resumable: an existing status.json is treated as cached unless --retry-failed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

SIZES = (1_000, 10_000, 20_000, 50_000, 100_000, 500_000, 1_000_000)
METHODS = {
    "graph": ["--clustering-method", "graph", "--no-prefilter"],
    "graph_prefilter": ["--clustering-method", "graph", "--force-prefilter"],
    "greedy_lazy": [
        "--clustering-method", "greedy", "--greedy-selection", "lazy-exact", "--no-prefilter",
    ],
}
# (label, terminal seed, k-mer seed threshold)
CONFIGURATIONS = (
    ("legacy_0.50", "contiguous", "0.50"),
    ("redesign_0.40", "all-column-pairs", "0.40"),
)
LABELS = {
    "graph": "Graph",
    "graph_prefilter": "Graph + prefilter",
    "greedy_lazy": "Greedy lazy-exact",
}


def peak_kbytes(path: Path) -> int | None:
    if not path.exists():
        return None
    match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", path.read_text())
    return int(match.group(1)) if match else None


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as reader, target.open("wb") as writer:
        shutil.copyfileobj(reader, writer, 8 * 1024 * 1024)


def directory_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def run_one(args, configuration: tuple[str, str, str], method: str, size: int) -> dict:
    label, geometry, seed_threshold = configuration
    output = args.output_dir / "runs" / label / method / f"n_{size:07d}"
    output.mkdir(parents=True, exist_ok=True)
    status_path = output / "status.json"
    if status_path.exists():
        previous = json.loads(status_path.read_text())
        if previous.get("status") == "ok" or not args.retry_failed:
            return previous

    source = args.data_dir / f"benchmark_{size:07d}.fasta.gz"
    if not source.exists():
        raise FileNotFoundError(source)
    run_tmp = args.tmp_root / label / method / f"n_{size:07d}"
    run_tmp.mkdir(parents=True, exist_ok=True)
    pepcluster_tmp = run_tmp / "pepcluster_tmp"
    resource = output / "resource.txt"
    # Node load before and after. These nodes are shared, so a timing taken
    # under heavy background load should be visible as such rather than silently
    # compared against a quiet one.
    load_before = os.getloadavg()[0]
    started = time.time()
    with tempfile.TemporaryDirectory(prefix="input_", dir=run_tmp) as temporary:
        fasta = Path(temporary) / "input.fasta"
        decompress(source, fasta)
        command = [
            str(args.binary),
            "--input", str(fasta),
            "--output-dir", str(output),
            "--tmp-dir", str(pepcluster_tmp),
            "--mode", "separate_aln_anchor",
            "--alignment-similarity-threshold", "0.50",
            "--anchor-combination-similarity-threshold", "0.60",
            "--kmer-seed-threshold", seed_threshold,
            "--terminal-seed", geometry,
            "--gap-open", "-4",
            "--gap-extension", "-1",
            "--terminal-overhang-gap-open", "-2",
            "--terminal-overhang-gap-extension", "-1",
            "--minimum-terminal-match-length", "2",
            "--threads", str(args.threads),
            "--candidate-buffer-mb", str(args.candidate_buffer_mb),
            "--compact-output",
            "--keep-tmp",
            *METHODS[method],
        ]
        if args.max_memory_gb is not None:
            command += ["--max-memory-gb", str(args.max_memory_gb)]
        (output / "benchmark_command.txt").write_text(" ".join(command) + "\n")
        timed = ["/usr/bin/time", "-v", "-o", str(resource), *command]
        with (output / "benchmark.log").open("w") as log:
            result = subprocess.run(timed, stdout=log, stderr=subprocess.STDOUT)
        temporary_bytes = directory_bytes(pepcluster_tmp)

    elapsed = time.time() - started
    kb = peak_kbytes(resource)
    status = {
        "configuration": label,
        "terminal_seed": geometry,
        "kmer_seed_threshold": seed_threshold,
        "method": method,
        "records": size,
        "status": "ok" if result.returncode == 0 else "failed",
        "returncode": result.returncode,
        "wall_seconds": elapsed,
        "peak_memory_mb": kb / 1024 if kb is not None else None,
        "temporary_disk_mb": temporary_bytes / 1024 / 1024,
        "node_load_before": load_before,
        "node_load_after": os.getloadavg()[0],
    }
    stats_path = output / "run_stats.json"
    if stats_path.exists():
        stats = json.loads(stats_path.read_text())
        for key in ("elapsed_seconds", "final_clusters", "singleton_clusters",
                    "index_candidate_hits", "anchor_bound_rejected",
                    "candidate_pairs_computed", "alignment_evaluations",
                    "graph_edge_count", "fraction_all_pairs_computed"):
            status[key] = stats.get(key)
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    # The candidate/edge spill can be very large; keep the measurement, not the data.
    if not args.keep_tmp:
        shutil.rmtree(pepcluster_tmp, ignore_errors=True)
    return status


def collect(args, configurations, methods, sizes) -> list[dict]:
    rows = []
    for label, _, _ in configurations:
        for method in methods:
            for size in sizes:
                path = (args.output_dir / "runs" / label / method / f"n_{size:07d}"
                        / "status.json")
                if path.exists():
                    rows.append(json.loads(path.read_text()))
    figures = args.output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    fields = [
        "configuration", "terminal_seed", "kmer_seed_threshold", "method", "records",
        "status", "returncode", "wall_seconds", "elapsed_seconds", "peak_memory_mb",
        "temporary_disk_mb", "node_load_before", "node_load_after",
        "index_candidate_hits", "anchor_bound_rejected",
        "candidate_pairs_computed", "alignment_evaluations", "graph_edge_count",
        "fraction_all_pairs_computed", "final_clusters", "singleton_clusters",
    ]
    with (figures / "scaling_benchmark.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot(args, rows, configurations, methods) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    successful = [row for row in rows if row["status"] == "ok"]
    if not successful:
        return
    colors = {"graph": "#31688e", "graph_prefilter": "#35b779", "greedy_lazy": "#7e3ace"}
    panels = [
        ("wall_seconds", "Wall time (s)"),
        ("peak_memory_mb", "Peak memory (MiB)"),
        ("temporary_disk_mb", "Temporary disk (MiB)"),
        ("candidate_pairs_computed", "Pairs exactly scored"),
    ]
    figure, axes = plt.subplots(1, 4, figsize=(19, 4.2), constrained_layout=True)
    for axis, (key, ylabel) in zip(axes, panels):
        for label, _, _ in configurations:
            for method in methods:
                selected = sorted(
                    (r for r in successful
                     if r["configuration"] == label and r["method"] == method
                     and r.get(key) is not None),
                    key=lambda r: r["records"],
                )
                if not selected:
                    continue
                axis.plot([r["records"] for r in selected], [r[key] for r in selected],
                          "o-" if label.startswith("redesign") else "o--",
                          color=colors[method],
                          label=f"{LABELS[method]} ({label})")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("Number of peptides")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=0.25)
    axes[0].legend(frameon=False, fontsize=7)
    figure.savefig(args.output_dir / "figures" / "scaling_benchmark.png", dpi=200)
    figure.savefig(args.output_dir / "figures" / "scaling_benchmark.pdf")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--tmp-root", type=Path, required=True)
    parser.add_argument("--methods", nargs="+", choices=list(METHODS), default=list(METHODS))
    parser.add_argument("--sizes", nargs="+", type=int, choices=SIZES, default=list(SIZES))
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--candidate-buffer-mb", type=int, default=4096)
    parser.add_argument("--max-memory-gb", type=float)
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--keep-tmp", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    for name in ("binary", "data_dir", "output_dir", "tmp_root"):
        setattr(args, name, getattr(args, name).resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.tmp_root.mkdir(parents=True, exist_ok=True)
    if not args.plot_only:
        # Configuration is the innermost loop so the historical and redesigned
        # settings for a given size run back to back. These nodes are shared, so
        # background load drifts over hours; running the two configurations
        # minutes apart keeps the comparison between them meaningful even when
        # the absolute seconds are not perfectly reproducible.
        for size in sorted(args.sizes):
            for method in args.methods:
                for configuration in CONFIGURATIONS:
                    print(json.dumps(run_one(args, configuration, method, size)), flush=True)
    rows = collect(args, CONFIGURATIONS, args.methods, args.sizes)
    plot(args, rows, CONFIGURATIONS, args.methods)


if __name__ == "__main__":
    main()
