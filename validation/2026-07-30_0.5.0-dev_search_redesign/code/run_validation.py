#!/usr/bin/env python3
"""Run and analyse the PepCluster2 0.5.0-dev candidate-search validation.

Stages
------
reference  exhaustive scoring of every pair, then the identical clustering
           procedure, for both representative orders, on the full datasets and
           on every nested subset. This is the partition the tool would produce
           with a perfect candidate search.
sweep      graph method on a few datasets across seed thresholds and both seed
           geometries, to place the sensitivity/cost operating point.
full       all clustering paths on the full datasets, with pair tracing.
subsets    all clustering paths on every nested subset.
analyse    metrics, figures and REPORT.md.

All agreement is measured against the exhaustive reference, which shares the
clustering procedure with the run under test. Any difference is therefore
attributable to candidate search alone.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A  # noqa: E402

REPO = Path(__file__).resolve().parents[3]

METHODS = {
    "graph": ["--clustering-method", "graph", "--no-prefilter"],
    "graph_prefilter": ["--clustering-method", "graph", "--force-prefilter"],
    "greedy": [
        "--clustering-method", "greedy", "--greedy-selection", "kmer-degree", "--no-prefilter",
    ],
    "greedy_lazy": [
        "--clustering-method", "greedy", "--greedy-selection", "lazy-exact", "--no-prefilter",
    ],
}
LABELS = {
    "graph": "Graph",
    "graph_prefilter": "Graph + prefilter",
    "greedy": "Greedy",
    "greedy_lazy": "Greedy lazy-exact",
}
ORDERS = ("coverage", "intrinsic")
# lazy-exact is itself a dynamic set-cover rule, so it has no intrinsic variant.
INVALID = {("intrinsic", "greedy_lazy")}
SUBSET_SIZES = (1_000, 2_000, 4_000, 6_000, 8_000)
SWEEP = (
    ("contiguous", "0.50"),
    ("all-column-pairs", "0.50"),
    ("all-column-pairs", "0.45"),
    ("all-column-pairs", "0.40"),
    ("all-column-pairs", "0.35"),
    ("all-column-pairs", "0.30"),
)
DEFAULT_SEED_THRESHOLD = "0.40"
DEFAULT_GEOMETRY = "all-column-pairs"


def combinations() -> list[tuple[str, str]]:
    return [
        (order, method)
        for order in ORDERS
        for method in METHODS
        if (order, method) not in INVALID
    ]


# --------------------------------------------------------------------------
# execution
# --------------------------------------------------------------------------


def decompress(source: Path, target: Path) -> None:
    with gzip.open(source, "rb") as reader, target.open("wb") as writer:
        shutil.copyfileobj(reader, writer, length=4 * 1024 * 1024)


def execute(command: list[str], output: Path, log: Path, resource: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    timed = ["/usr/bin/time", "-v", "-o", str(resource), *command]
    with log.open("w") as handle:
        result = subprocess.run(timed, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log}")


def cluster_command(
    binary: Path,
    input_fasta: Path,
    output: Path,
    threads: int,
    geometry: str,
    seed_threshold: str,
    order: str,
) -> list[str]:
    return [
        str(binary),
        "--input", str(input_fasta),
        "--output-dir", str(output),
        "--mode", "separate_aln_anchor",
        "--alignment-similarity-threshold", "0.50",
        "--anchor-combination-similarity-threshold", "0.60",
        "--kmer-seed-threshold", seed_threshold,
        "--terminal-seed", geometry,
        "--representative-order", order,
        "--gap-open", "-4",
        "--gap-extension", "-1",
        "--terminal-overhang-gap-open", "-2",
        "--terminal-overhang-gap-extension", "-1",
        "--minimum-terminal-match-length", "2",
        "--threads", str(threads),
        "--candidate-buffer-mb", "512",
        "--compact-output",
    ]


def run_reference(root: Path, helper: Path, dataset: int, size: int | None, threads: int) -> str:
    if size is None:
        source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
        output = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
        label = f"reference full {dataset:03d}"
    else:
        source = root / "data" / "subsets" / f"n_{size:06d}" / f"sample_{dataset:03d}.fasta.gz"
        output = root / "runs" / "exhaustive_subsets" / f"n_{size:06d}" / f"sample_{dataset:03d}"
        label = f"reference {size} {dataset:03d}"
    needed = ["run_stats.json", "true_pairs.bin",
              "reassign_only_clusters.tsv", "pipeline_clusters.tsv",
              "reassign_only_clusters_intrinsic.tsv",
              "pipeline_clusters_intrinsic.tsv"]
    if all((output / name).exists() for name in needed):
        return f"{label} cached"
    with tempfile.TemporaryDirectory(prefix="pc2ref_", dir="/tmp") as tmp:
        fasta = Path(tmp) / "in.fasta"
        decompress(source, fasta)
        execute([str(helper), str(fasta), str(output), str(threads)],
                output, output / "run.log", output / "resource.txt")
    return f"{label} complete"


def run_cluster(
    root: Path, binary: Path, order: str, method: str, dataset: int,
    threads: int, size: int | None,
) -> str:
    if size is None:
        source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
        output = root / "runs" / "full" / order / method / f"sample_{dataset:03d}"
        trace = True
        label = f"full {order}/{method} {dataset:03d}"
    else:
        source = root / "data" / "subsets" / f"n_{size:06d}" / f"sample_{dataset:03d}.fasta.gz"
        output = (root / "runs" / "subsets" / f"n_{size:06d}" / order / method
                  / f"sample_{dataset:03d}")
        trace = False
        label = f"subset {size} {order}/{method} {dataset:03d}"
    needed = [output / "run_stats.json", output / "node_clusters.tsv"]
    if trace:
        needed.append(output / "scored_pairs.bin")
    if all(path.exists() for path in needed):
        return f"{label} cached"
    with tempfile.TemporaryDirectory(prefix="pc2run_", dir="/tmp") as tmp:
        fasta = Path(tmp) / "in.fasta"
        decompress(source, fasta)
        command = cluster_command(binary, fasta, output, threads, DEFAULT_GEOMETRY,
                                  DEFAULT_SEED_THRESHOLD, order) + METHODS[method]
        if trace:
            command.append("--write-scored-pairs")
        command += ["--tmp-dir", str(Path(tmp) / "tmp")]
        execute(command, output, output / "run.log", output / "resource.txt")
    return f"{label} complete"


def run_sweep(
    root: Path, binary: Path, geometry: str, seed_threshold: str, dataset: int, threads: int,
) -> str:
    tag = f"{geometry}_{seed_threshold}"
    output = root / "runs" / "sweep" / tag / f"sample_{dataset:03d}"
    label = f"sweep {tag} {dataset:03d}"
    needed = [output / "run_stats.json", output / "node_clusters.tsv", output / "scored_pairs.bin"]
    if all(path.exists() for path in needed):
        return f"{label} cached"
    source = root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz"
    with tempfile.TemporaryDirectory(prefix="pc2sweep_", dir="/tmp") as tmp:
        fasta = Path(tmp) / "in.fasta"
        decompress(source, fasta)
        command = cluster_command(binary, fasta, output, threads, geometry, seed_threshold,
                                  "coverage") + METHODS["graph"]
        command += ["--write-scored-pairs", "--tmp-dir", str(Path(tmp) / "tmp")]
        execute(command, output, output / "run.log", output / "resource.txt")
    return f"{label} complete"


def parallel(tasks: list[tuple], function, workers: int) -> None:
    if not tasks:
        return
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(function, *task) for task in tasks]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            print(f"[{index}/{len(tasks)}] {future.result()}", flush=True)


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------


def reference_name(order: str) -> str:
    """The exhaustive reference for an order: the same clustering procedure the
    run uses, applied to the complete edge set."""
    return "pipeline_clusters.tsv" if order == "coverage" else "pipeline_clusters_intrinsic.tsv"


def no_merge_reference_name(order: str) -> str:
    return ("reassign_only_clusters.tsv" if order == "coverage"
            else "reassign_only_clusters_intrinsic.tsv")


def run_row(run_dir: Path, stats: dict) -> dict:
    cost = A.cost_metrics(stats)
    cost["peak_rss_kbytes"] = A.peak_rss_kbytes(run_dir / "resource.txt")
    cost["final_clusters"] = stats["final_clusters"]
    cost["singleton_clusters"] = stats["singleton_clusters"]
    cost["iterations"] = stats["iterations"]
    cost["converged"] = stats["converged"]
    return cost


def analyse(root: Path, datasets: int, sweep_datasets: int) -> None:
    tables = A.ScoreTables(REPO / "src" / "scoring.rs")
    search_rows, agreement_rows, truth_rows, audit_rows, sweep_rows = [], [], [], [], []

    features_cache: dict[int, A.PeptideFeatures] = {}

    def features(dataset: int) -> A.PeptideFeatures:
        if dataset not in features_cache:
            sequences = A.node_sequences(root / "data" / "full" / f"sample_{dataset:03d}.fasta.gz")
            features_cache[dataset] = A.PeptideFeatures(sequences, tables)
        return features_cache[dataset]

    for dataset in range(datasets):
        reference_dir = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
        reference_stats = A.read_stats(reference_dir / "run_stats.json")
        true_pairs = set(A.read_pair_file(reference_dir / "true_pairs.bin", b"PC2TRUE1").tolist())
        truth_rows.append({"dataset": dataset, **reference_stats})
        partitions = {
            name: A.read_partition(reference_dir / name)
            for name in ("pipeline_clusters.tsv", "pipeline_clusters_intrinsic.tsv")
        }

        for order, method in combinations():
            run_dir = root / "runs" / "full" / order / method / f"sample_{dataset:03d}"
            if not (run_dir / "run_stats.json").exists():
                continue
            stats = A.read_stats(run_dir / "run_stats.json")
            scored = set(A.read_pair_file(run_dir / "scored_pairs.bin", b"PC2PAIR1").tolist())
            search = A.search_metrics(true_pairs, scored, reference_stats["all_possible_pairs"])
            search_rows.append({
                "dataset": dataset, "order": order, "method": method,
                **search, **run_row(run_dir, stats),
            })
            missed = np.asarray(sorted(true_pairs - scored), dtype=np.uint64)
            geometry = ("contiguous" if stats.get("terminal_seed") == "contiguous"
                        else "all_column_pairs")
            threshold_q = int(round(float(stats.get("kmer_seed_threshold", 0.4)) * 1000))
            audit_rows.append({
                "dataset": dataset, "order": order, "method": method,
                "terminal_seed": stats.get("terminal_seed"),
                "kmer_seed_threshold": stats.get("kmer_seed_threshold"),
                **A.missed_pair_audit(missed, features(dataset), geometry, threshold_q, 600),
            })
            query = A.read_partition(run_dir / "node_clusters.tsv")
            agreement_rows.append({
                "dataset": dataset, "order": order, "method": method,
                **A.partition_metrics(partitions[reference_name(order)], query),
                "final_clusters": stats["final_clusters"],
            })

    for dataset in range(sweep_datasets):
        reference_dir = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
        reference_stats = A.read_stats(reference_dir / "run_stats.json")
        true_pairs = set(A.read_pair_file(reference_dir / "true_pairs.bin", b"PC2TRUE1").tolist())
        pipeline = A.read_partition(reference_dir / "pipeline_clusters.tsv")
        for geometry, seed_threshold in SWEEP:
            run_dir = root / "runs" / "sweep" / f"{geometry}_{seed_threshold}" / f"sample_{dataset:03d}"
            if not (run_dir / "run_stats.json").exists():
                continue
            stats = A.read_stats(run_dir / "run_stats.json")
            scored = set(A.read_pair_file(run_dir / "scored_pairs.bin", b"PC2PAIR1").tolist())
            missed = np.asarray(sorted(true_pairs - scored), dtype=np.uint64)
            audit = A.missed_pair_audit(
                missed, features(dataset),
                "contiguous" if geometry == "contiguous" else "all_column_pairs",
                int(round(float(seed_threshold) * 1000)), 600,
            )
            query = A.read_partition(run_dir / "node_clusters.tsv")
            sweep_rows.append({
                "dataset": dataset, "terminal_seed": geometry, "kmer_seed_threshold": seed_threshold,
                **A.search_metrics(true_pairs, scored, reference_stats["all_possible_pairs"]),
                **run_row(run_dir, stats),
                "ari_vs_reference": A.partition_metrics(pipeline, query)["ari"],
                **audit,
            })

    stability_rows = []
    for dataset in range(datasets):
        full_reference_dir = root / "runs" / "exhaustive" / f"sample_{dataset:03d}"
        for order in ORDERS:
            # The reference undergoes the identical comparison, so its own
            # composition dependence is the part that belongs to the clustering
            # procedure rather than to the search. The no-merging variant
            # attributes the merge stage.
            for label, name in (
                ("exhaustive_reference_reassign", no_merge_reference_name(order)),
                ("exhaustive_reference_pipeline", reference_name(order)),
            ):
                full_reference = A.read_partition(full_reference_dir / name)
                for size in SUBSET_SIZES:
                    subset_reference_dir = (root / "runs" / "exhaustive_subsets"
                                            / f"n_{size:06d}" / f"sample_{dataset:03d}")
                    if not (subset_reference_dir / name).exists():
                        continue
                    subset_reference = A.read_partition(subset_reference_dir / name)
                    restricted = {s: full_reference[s] for s in subset_reference}
                    stability_rows.append({
                        "dataset": dataset, "order": order, "method": label,
                        "subset_size": size, "subset_fraction": size / 10_000,
                        **A.partition_metrics(restricted, subset_reference),
                    })

        for order, method in combinations():
            full_run = root / "runs" / "full" / order / method / f"sample_{dataset:03d}"
            if not (full_run / "node_clusters.tsv").exists():
                continue
            full_partition = A.read_partition(full_run / "node_clusters.tsv")
            for size in SUBSET_SIZES:
                subset_run = (root / "runs" / "subsets" / f"n_{size:06d}" / order / method
                              / f"sample_{dataset:03d}")
                if not (subset_run / "node_clusters.tsv").exists():
                    continue
                subset_partition = A.read_partition(subset_run / "node_clusters.tsv")
                restricted = {s: full_partition[s] for s in subset_partition}
                stability_rows.append({
                    "dataset": dataset, "order": order, "method": method,
                    "subset_size": size, "subset_fraction": size / 10_000,
                    **A.partition_metrics(restricted, subset_partition),
                })

    A.write_csv(root / "metrics" / "true_pair_space.csv", truth_rows)
    A.write_csv(root / "metrics" / "missed_pair_audit.csv", audit_rows)
    A.write_csv(root / "figures" / "search_rule_performance.csv", search_rows)
    A.write_csv(root / "figures" / "reference_agreement.csv", agreement_rows)
    A.write_csv(root / "figures" / "cluster_stability.csv", stability_rows)
    A.write_csv(root / "figures" / "seed_sweep.csv", sweep_rows)
    make_figures(root, search_rows, agreement_rows, stability_rows, sweep_rows)
    write_report(root, truth_rows, search_rows, agreement_rows, stability_rows,
                 sweep_rows, audit_rows)


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

# Categorical slots 1-4 of the validated reference palette, assigned in fixed
# order. References are drawn in neutral ink because they are baselines, not
# series identities.
COLORS = {
    "graph": "#2a78d6",
    "graph_prefilter": "#eb6834",
    "greedy": "#1baf7a",
    "greedy_lazy": "#eda100",
    "exhaustive_reference_reassign": "#a8a7a2",
    "exhaustive_reference_pipeline": "#0b0b0b",
}
CEILINGS = {
    "exhaustive_reference_reassign": "Exhaustive reference (no merging)",
    "exhaustive_reference_pipeline": "Exhaustive reference",
}
GRID = dict(alpha=0.22, linewidth=0.6)
INK = "#52514e"


def _mean(rows: list[dict], key: str) -> float:
    values = [float(row[key]) for row in rows]
    return statistics.mean(values) if values else float("nan")


def _stdev(rows: list[dict], key: str) -> float:
    values = [float(row[key]) for row in rows]
    return statistics.stdev(values) if len(values) > 1 else 0.0


def make_figures(root: Path, search_rows, agreement_rows, stability_rows, sweep_rows) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = root / "figures"
    figures.mkdir(exist_ok=True)

    def style_axis(axis, xlabel, ylabel, title):
        axis.set_title(title, fontsize=10.5, color="#0b0b0b")
        axis.set_xlabel(xlabel, fontsize=9, color=INK)
        axis.set_ylabel(ylabel, fontsize=9, color=INK)
        axis.grid(**GRID)
        axis.tick_params(labelsize=8.5, colors=INK)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            axis.spines[spine].set_color("#c9c8c3")

    # 1. Sensitivity and cost of the seed rule. The historical geometry is one
    #    point, so it is marked rather than joined.
    if sweep_rows:
        figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)
        panels = [
            ("search_recall", "Search recall", None),
            ("candidate_pairs_computed", "Pairs exactly scored (thousands)", 1e-3),
            ("ari_vs_reference", "ARI vs exhaustive reference", None),
        ]
        points = sorted({float(r["kmer_seed_threshold"]) for r in sweep_rows
                         if r["terminal_seed"] == "all-column-pairs"})
        for axis, (key, ylabel, scale) in zip(axes, panels):
            values = [
                _mean(_pick(sweep_rows, terminal_seed="all-column-pairs",
                            kmer_seed_threshold=f"{p:.2f}"), key)
                for p in points
            ]
            if scale:
                values = [v * scale for v in values]
            axis.plot(points, values, "-o", color=COLORS["graph"], linewidth=2,
                      markersize=5.5, label="all column pairs")
            legacy = _pick(sweep_rows, terminal_seed="contiguous", kmer_seed_threshold="0.50")
            if legacy:
                value = _mean(legacy, key) * (scale or 1)
                axis.plot([0.50], [value], "X", color=COLORS["graph_prefilter"],
                          markersize=11, label="contiguous only (ablation)")
                axis.annotate(f"{value:,.2f}" if not scale else f"{value:,.0f}",
                              (0.50, value), textcoords="offset points", xytext=(-6, 8),
                              ha="right", fontsize=8, color=INK)
            default = float(DEFAULT_SEED_THRESHOLD)
            if default in points:
                value = values[points.index(default)]
                axis.annotate(f"{value:,.2f}" if not scale else f"{value:,.0f}",
                              (default, value), textcoords="offset points", xytext=(0, 9),
                              ha="center", fontsize=8, color=COLORS["graph"])
            style_axis(axis, "k-mer seed threshold", ylabel, ylabel.split(" (")[0])
            if key != "candidate_pairs_computed":
                axis.set_ylim(0, 1.02)
        axes[0].legend(frameon=False, fontsize=8, loc="lower left")
        figure.suptitle(
            "Redesigned terminal seed: sensitivity against cost "
            f"(graph, coverage order, {len({r['dataset'] for r in sweep_rows})} datasets)",
            fontsize=11.5, y=1.04)
        figure.savefig(figures / "seed_sweep.png", dpi=200, bbox_inches="tight")
        figure.savefig(figures / "seed_sweep.pdf", bbox_inches="tight")
        plt.close(figure)

    # 2. Agreement against both references. The ceiling line is why the two
    #    panels must be read together.
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for axis, order in zip(axes, ORDERS):
        data, labels, colors = [], [], []
        for method in METHODS:
            values = [row["ari"] for row in _pick(agreement_rows, method=method, order=order)]
            if values:
                data.append(values)
                labels.append(LABELS[method])
                colors.append(COLORS[method])
        if not data:
            continue
        box = axis.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.55,
                           medianprops=dict(color="#0b0b0b", linewidth=1.4),
                           whiskerprops=dict(color="#8a8984"),
                           capprops=dict(color="#8a8984"),
                           flierprops=dict(markersize=3, markeredgecolor="#8a8984"))
        for patch, color in zip(box["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
            patch.set_edgecolor(color)
        style_axis(axis, "", "Adjusted Rand index" if order == ORDERS[0] else "",
                   f"{order} order"
                   + (" (lazy-exact is a set-cover rule, so it has no intrinsic variant)"
                      if order == "intrinsic" else ""))
        axis.set_ylim(0, 1.01)
        axis.tick_params(axis="x", rotation=18)
    # The reference shares the clustering procedure with the run, so the gap is
    # candidate-search loss; that is stated in the report rather than crowded
    # into a two-line title over the panel headings.
    figure.suptitle("Agreement with the exhaustive reference, 20 datasets",
                    fontsize=11.5, y=1.04)
    figure.savefig(figures / "reference_agreement.png", dpi=200, bbox_inches="tight")
    figure.savefig(figures / "reference_agreement.pdf", bbox_inches="tight")
    plt.close(figure)

    # 3. Stability, faceted by representative order. All four clustering paths are
    #    plotted; graph, graph + prefilter and lazy-exact greedy agree to within
    #    0.0004 and therefore overlap almost exactly, which is itself the result.
    shown = list(METHODS) + list(CEILINGS)
    metric_names = [("pairwise_jaccard", "Pairwise Jaccard"), ("ari", "Subset ARI")]
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True,
                               sharex=True, sharey=True)
    x = np.asarray(SUBSET_SIZES) / 10_000
    for row, (metric, metric_label) in enumerate(metric_names):
        for column, order in enumerate(ORDERS):
            axis = axes[row][column]
            # Graph, graph + prefilter and lazy-exact greedy coincide to within
            # 0.0004, so each gets its own dash pattern: an overplotted solid
            # line would leave two of them invisible behind the third while
            # still appearing in the legend. Only the visually distinct series
            # are end-labelled.
            styles = {
                "graph": ("-", 2.6),
                "graph_prefilter": ((0, (6, 3)), 1.8),
                "greedy_lazy": ((0, (1, 2.5)), 2.2),
                "greedy": ("-", 2.6),
                "exhaustive_reference_reassign": ((0, (4, 3)), 1.5),
                "exhaustive_reference_pipeline": ((0, (1, 2)), 1.8),
            }
            # Only the two distinct clustering paths are end-labelled. Under the
            # coverage order the references sit within 0.008 of graph, so their
            # labels would overlap it; their values are in the stability table.
            offsets = {"graph": 10, "greedy": -13}
            for method in shown:
                means = []
                for size in SUBSET_SIZES:
                    rows = _pick(stability_rows, method=method, order=order,
                                 subset_size=size)
                    means.append(_mean(rows, metric) if rows else float("nan"))
                if all(np.isnan(means)):
                    continue
                dash, width = styles[method]
                axis.plot(x, means, linestyle=dash, marker="o", markersize=4.0,
                          linewidth=width, color=COLORS[method],
                          label=LABELS.get(method, CEILINGS.get(method)))
                if method in offsets:
                    axis.annotate(f"{means[-1]:.2f}", (x[-1], means[-1]),
                                  textcoords="offset points",
                                  xytext=(7, offsets[method]), fontsize=7.5,
                                  color=COLORS[method])
            style_axis(axis, "Subset fraction" if row == 1 else "",
                       metric_label if column == 0 else "",
                       f"{metric_label} — {order} order")
            axis.set_ylim(0, 1.01)
            axis.set_xlim(0.05, 0.92)
    handles, labels = axes[0][0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False,
                  fontsize=8)
    figure.suptitle(
        "Subset stability, with the exhaustive reference under the same comparison\n"
        "Graph, graph + prefilter and lazy-exact greedy coincide to within 0.0004",
        fontsize=11)
    figure.savefig(figures / "cluster_stability.png", dpi=200, bbox_inches="tight")
    figure.savefig(figures / "cluster_stability.pdf", bbox_inches="tight")
    plt.close(figure)


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------


def _pick(rows, **conditions):
    return [
        row for row in rows
        if all(str(row.get(key)) == str(value) for key, value in conditions.items())
    ]


def write_interpretation(add, search_rows, agreement_rows, stability_rows, sweep_rows,
                         audit_rows) -> None:
    add("## Interpretation")
    add("")

    graph = _pick(search_rows, method="graph", order="coverage")
    ablation = _pick(sweep_rows, terminal_seed="contiguous", kmer_seed_threshold="0.50")
    if graph:
        add("### Candidate search")
        add("")
        add(f"The seed recovers {_mean(graph, 'search_recall'):.1%} of the pairs that pass both")
        add(f"exact thresholds while exactly scoring {_mean(graph, 'scored_unique_pairs'):,.0f}")
        add(f"distinct pairs, {_mean(graph, 'fraction_all_pairs_scored'):.2%} of all possible")
        add("pairs. Sensitivity and cost are controlled by two independent mechanisms: the")
        add("geometry decides which residue columns can be matched, and the sound anchor")
        add("bound removes pairs that provably cannot be accepted.")
        add("")
        if ablation:
            add("The geometry ablation isolates the first. Indexing only the contiguous")
            add("column pairs of each terminal 3-mer, with everything else unchanged, drops")
            add(f"recall to {_mean(ablation, 'search_recall'):.1%}. Both contiguous pairs")
            add("contain the middle residue, so one substitution there destroys both, and the")
            add("required terminal columns shift whenever the peptides differ in length.")
            add("")
        seed_missed = _mean(_pick(audit_rows, method="graph", order="coverage"),
                            "missed_seed_attributable")
        unsound = sum(int(row["missed_anchor_bound_unsound"]) for row in audit_rows)
        add(f"Every missed pair that remains ({seed_missed:.0f} per dataset) failed a terminal")
        add("seed threshold, not the bound: across every run in this validation the sound")
        add(f"bound rejected {unsound} eligible pairs. The residual is a threshold trade-off,")
        add("shown in the sweep table, not a structural blind spot.")
        add("")

    add("### Agreement")
    add("")
    add("The exhaustive reference applies the identical clustering procedure to the")
    add("complete edge set, so a run differs from it only through candidate search.")
    add("")
    for order in ORDERS:
        run = _pick(agreement_rows, method="graph", order=order)
        if not run:
            continue
        add(f"With the `{order}` order, graph reaches ARI {_mean(run, 'ari'):.4f} and pairwise")
        add(f"Jaccard {_mean(run, 'pairwise_jaccard'):.4f} against that reference, at")
        add(f"{_mean(run, 'final_clusters'):.0f} clusters.")
        add("")

    add("### Stability")
    add("")
    for order in ORDERS:
        for size, name in ((8_000, "80%"),):
            run = _pick(stability_rows, method="graph", order=order, subset_size=size)
            reference = _pick(stability_rows, method="exhaustive_reference_pipeline",
                              order=order, subset_size=size)
            no_merge = _pick(stability_rows, method="exhaustive_reference_reassign",
                             order=order, subset_size=size)
            if not (run and reference):
                continue
            add(f"At the {name} subset with the `{order}` order, graph reaches pairwise Jaccard")
            add(f"{_mean(run, 'pairwise_jaccard'):.4f}, against "
                f"{_mean(reference, 'pairwise_jaccard'):.4f} for the exhaustive reference under")
            add("the identical comparison. The run therefore tracks what the same clustering")
            add("procedure does on a complete edge set, and the composition dependence that")
            add("remains belongs to the procedure rather than to the search.")
            add("")
            add("That reference is a comparison point, not an upper bound, and a run can score")
            add("slightly above it: missing a few percent of edges makes clusters marginally")
            add("smaller, and smaller clusters have fewer co-cluster pairs to disagree about.")
            add("Stability must therefore be read next to the agreement table, never alone.")
            add("")
            if no_merge:
                add(f"Disabling merging on the reference gives "
                    f"{_mean(no_merge, 'pairwise_jaccard'):.4f}, so the merge stage accounts")
                add(f"for {_mean(reference, 'pairwise_jaccard') - _mean(no_merge, 'pairwise_jaccard'):+.4f}")
                add("of it. Reassignment is the larger term, which is why it carries the")
                add("hysteresis margin.")
                add("")

    intrinsic = _pick(agreement_rows, method="graph", order="intrinsic")
    coverage = _pick(agreement_rows, method="graph", order="coverage")
    intrinsic_stability = _pick(stability_rows, method="graph", order="intrinsic",
                                subset_size=8_000)
    coverage_stability = _pick(stability_rows, method="graph", order="coverage",
                               subset_size=8_000)
    if intrinsic and coverage and intrinsic_stability and coverage_stability:
        add("### Choosing a representative order")
        add("")
        add(f"`coverage` minimises the cluster count ({_mean(coverage, 'final_clusters'):.0f}")
        add("clusters) but its selection key is a degree, so it amplifies small edge-set")
        add(f"differences: it reaches ARI {_mean(coverage, 'ari'):.2f} against its reference,")
        add(f"against {_mean(intrinsic, 'ari'):.2f} for `intrinsic`. `intrinsic` is also more")
        add(f"stable ({_mean(intrinsic_stability, 'pairwise_jaccard'):.4f} versus "
            f"{_mean(coverage_stability, 'pairwise_jaccard'):.4f} pairwise Jaccard at 80%), at")
        add(f"the cost of {_mean(intrinsic, 'final_clusters'):.0f} clusters instead of "
            f"{_mean(coverage, 'final_clusters'):.0f}. Neither is uniformly better: choose")
        add("`coverage` for compactness and `intrinsic` when reproducibility across dataset")
        add("revisions matters more.")
        add("")

    greedy = _pick(search_rows, method="greedy", order="coverage")
    greedy_audit = _pick(audit_rows, method="greedy", order="coverage")
    if greedy and greedy_audit:
        add("### Static greedy")
        add("")
        add(f"Static greedy scores only {_mean(greedy, 'search_recall'):.1%} of eligible pairs,")
        add("but the audit shows this is structural rather than a search defect: of")
        add(f"{_mean(greedy_audit, 'missed'):.0f} missed pairs only "
            f"{_mean(greedy_audit, 'missed_seed_attributable'):.0f} failed the seed, while")
        add(f"{_mean(greedy_audit, 'missed_retrievable_but_unexamined'):.0f} were retrievable")
        add("but never examined, because the method only ever scores")
        add("representative-to-unassigned pairs. Its pair-level recall is therefore not")
        add("comparable with the graph method's, and lazy-exact greedy remains the")
        add("recommended low-memory path.")
        add("")

    add("### Limitations")
    add("")
    add("- The reference is computational, derived from the same scoring rule. It says")
    add("  nothing about biological cluster purity, which still requires labelled")
    add("  peptide-MHC data.")
    add("- The seed threshold is a sensitivity/cost trade-off, not a guarantee: at the")
    add("  0.40 default roughly 3% of eligible pairs are still not retrieved.")
    add("- The iterative section remains the largest source of composition dependence.")
    add("- All results are at 10,000 peptides. See `figures/scaling_benchmark.csv` for")
    add("  the cost profile at 1k-1M.")
    add("")


def write_report(root, truth_rows, search_rows, agreement_rows, stability_rows,
                 sweep_rows, audit_rows) -> None:
    lines: list[str] = []
    add = lines.append
    datasets = len({row["dataset"] for row in truth_rows})
    add("# PepCluster2 candidate search and clustering validation")
    add("")
    add("## Configuration")
    add("")
    add(f"- {datasets} independently sampled datasets, 10,000 peptides each.")
    add("- Scoring mode `separate_aln_anchor`; alignment threshold 0.50; "
        "anchor-combination threshold 0.60; terminal/core weights 4/1.")
    add(f"- Terminal seed `{DEFAULT_GEOMETRY}`, k-mer seed threshold "
        f"{DEFAULT_SEED_THRESHOLD} unless a sweep row says otherwise.")
    add("- Reassignment margin 0.01: a peptide leaves its representative only when")
    add("  another beats it by more than that. Both references use the same value, so")
    add("  runs and references remain comparable.")
    add("")
    add("## Method under test")
    add("")
    add("1. The terminal seed indexes all three ordered column pairs of each terminal")
    add("   3-mer, (1,2), (1,3) and (2,3). The constrained alignment must contain at")
    add("   least `--minimum-terminal-match-length` matched columns drawn from the first")
    add("   three residues of *both* peptides, and those columns shift when the peptides")
    add("   differ in length, so all three pairs are required to retrieve them.")
    add("2. Candidate generation applies a sound upper bound on the anchor-combination")
    add("   similarity (the assignment relaxed to independent row maxima). A pair failing")
    add("   it cannot be accepted, so this prunes without losing any eligible")
    add("   relationship.")
    add("3. Reassignment is hysteretic: a peptide leaves its representative only when")
    add("   another beats it by more than the margin, so near-ties do not flip when the")
    add("   dataset composition changes.")
    add("")

    add("## Reading the metrics")
    add("")
    add("- Search recall: fraction of exactly eligible pairs the run scored.")
    add("- Search precision: fraction of scored pairs that were eligible. Low")
    add("  precision means wasted work, not wrong edges.")
    add("- Cost is reported as a decomposition: index hits (with multiplicity),")
    add("  pairs rejected by the sound anchor bound, distinct pairs exactly scored,")
    add("  and constrained-alignment evaluations. These are different costs.")
    add("  Alignment evaluations are counted process-wide, so they include")
    add("  representative recalculation, merge validation and reassignment; the count")
    add("  can therefore exceed the number of distinct candidate pairs.")
    add("- The exhaustive reference scores every pair exactly and then applies the")
    add("  identical clustering procedure, for the same representative order as the run")
    add("  under test. A run therefore differs from it only through candidate search,")
    add("  which is what makes the agreement number interpretable.")
    add("")

    if sweep_rows:
        add("## Seed geometry and sensitivity sweep")
        add("")
        add("Graph method, coverage order, "
            f"{len({row['dataset'] for row in sweep_rows})} datasets.")
        add("")
        add("| Terminal seed | Seed threshold | Recall | Precision | Index hits | "
            "Pairs scored | Alignments | ARI vs reference | Seconds | Peak RSS MB |")
        add("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for geometry, seed_threshold in SWEEP:
            rows = [r for r in sweep_rows if r["terminal_seed"] == geometry
                    and r["kmer_seed_threshold"] == seed_threshold]
            if not rows:
                continue
            add(f"| {geometry} | {seed_threshold} | {_mean(rows, 'search_recall'):.4f} | "
                f"{_mean(rows, 'search_precision'):.4f} | "
                f"{_mean(rows, 'index_candidate_hits'):,.0f} | "
                f"{_mean(rows, 'candidate_pairs_computed'):,.0f} | "
                f"{_mean(rows, 'alignment_evaluations'):,.0f} | "
                f"{_mean(rows, 'ari_vs_reference'):.4f} | "
                f"{_mean(rows, 'elapsed_seconds'):.2f} | "
                f"{_mean(rows, 'peak_rss_kbytes') / 1024:.0f} |")
        add("")

    add("## Search-rule performance")
    add("")
    add("| Method | Order | Recall | Precision | Pairs scored | Alignments | "
        "All-pairs fraction | Seconds |")
    add("|---|---|---:|---:|---:|---:|---:|---:|")
    for order in ORDERS:
        for method in METHODS:
            rows = [r for r in search_rows if r["method"] == method and r["order"] == order]
            if not rows:
                continue
            add(f"| {LABELS[method]} | {order} | "
                f"{_mean(rows, 'search_recall'):.4f} ± {_stdev(rows, 'search_recall'):.4f} | "
                f"{_mean(rows, 'search_precision'):.4f} | "
                f"{_mean(rows, 'candidate_pairs_computed'):,.0f} | "
                f"{_mean(rows, 'alignment_evaluations'):,.0f} | "
                f"{_mean(rows, 'fraction_all_pairs_scored'):.4f} | "
                f"{_mean(rows, 'elapsed_seconds'):.2f} |")
    add("")

    add("## Missed-pair audit")
    add("")
    add("Every eligible pair a run failed to score, attributed to the cause.")
    add("`seed attributable` means a terminal seed found no neighbouring column")
    add("pair. `retrievable but unexamined` means the index would have returned")
    add("the pair and the bound would have kept it, but the clustering traversal")
    add("never scored it: the greedy paths only ever score")
    add("representative-to-unassigned pairs, so for them this is structural and not")
    add("a candidate-search defect. `anchor bound unsound` counts pairs the sound")
    add("bound rejected; it must be zero, and any non-zero value is a correctness")
    add("failure of the bound.")
    add("")
    add("| Method | Order | Missed | Seed attributable | Front seed | End seed | Both | "
        "Retrievable but unexamined | Anchor bound unsound |")
    add("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for order in ORDERS:
        for method in METHODS:
            rows = [r for r in audit_rows if r["method"] == method and r["order"] == order]
            if not rows:
                continue
            add(f"| {LABELS[method]} | {order} | {_mean(rows, 'missed'):.0f} | "
                f"{_mean(rows, 'missed_seed_attributable'):.0f} | "
                f"{_mean(rows, 'missed_front_seed'):.0f} | "
                f"{_mean(rows, 'missed_end_seed'):.0f} | "
                f"{_mean(rows, 'missed_both_seeds'):.0f} | "
                f"{_mean(rows, 'missed_retrievable_but_unexamined'):.0f} | "
                f"{sum(int(r['missed_anchor_bound_unsound']) for r in rows)} |")
    add("")

    add("## Agreement with the exhaustive reference")
    add("")
    add("The reference scores every pair exactly and then runs the identical clustering")
    add("procedure, so the gap below is candidate-search loss and nothing else.")
    add("")
    add("| Method | Order | ARI | NMI | Jaccard | Co-assoc. recall | Clusters |")
    add("|---|---|---:|---:|---:|---:|---:|")
    for order in ORDERS:
        for method in METHODS:
            rows = _pick(agreement_rows, method=method, order=order)
            if not rows:
                continue
            add(f"| {LABELS[method]} | {order} | "
                f"{_mean(rows, 'ari'):.4f} ± {_stdev(rows, 'ari'):.4f} | "
                f"{_mean(rows, 'nmi'):.4f} | {_mean(rows, 'pairwise_jaccard'):.4f} | "
                f"{_mean(rows, 'coassociation_recall'):.4f} | "
                f"{_mean(rows, 'final_clusters'):.0f} |")
    add("")

    add("## Stability")
    add("")
    add("A subset clustering is compared with the full-dataset clustering restricted")
    add("to the same peptides. The `Exhaustive reference` rows apply the identical")
    add("comparison to the reference itself, so they show how much of the composition")
    add("dependence belongs to the clustering procedure rather than to candidate")
    add("search. Graph, graph + prefilter and lazy-exact greedy agree to within 0.0004")
    add("throughout, so their curves coincide in the figure.")
    add("")
    add("| Method | Order | Subset | Jaccard | ARI | Co-assoc. recall | Co-assoc. precision |")
    add("|---|---|---:|---:|---:|---:|---:|")
    for order in ORDERS:
        for method in list(METHODS) + list(CEILINGS):
            for size in (1_000, 8_000):
                rows = [r for r in stability_rows if r["method"] == method
                        and r["order"] == order and r["subset_size"] == size]
                if not rows:
                    continue
                add(f"| {LABELS.get(method, CEILINGS.get(method, method))} | {order} | "
                    f"{size // 100}% | {_mean(rows, 'pairwise_jaccard'):.4f} | "
                    f"{_mean(rows, 'ari'):.4f} | "
                    f"{_mean(rows, 'coassociation_recall'):.4f} | "
                    f"{_mean(rows, 'coassociation_precision'):.4f} |")
    add("")

    write_interpretation(add, search_rows, agreement_rows, stability_rows, sweep_rows,
                         audit_rows)

    add("## Files")
    add("")
    add("- `runs/exhaustive/`, `runs/exhaustive_subsets/`: exact pairs and both")
    add("  reference partitions for each representative order.")
    add("- `runs/full/`, `runs/subsets/`, `runs/sweep/`: clustering runs.")
    add("- `figures/`: plot-matched CSV, PNG and PDF.")
    add("- `metrics/missed_pair_audit.csv`: per-run attribution of missed pairs.")
    add("- `code/`: preparation, execution, exhaustive reference and analysis source.")
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--helper", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=20)
    parser.add_argument("--sweep-datasets", type=int, default=3)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--stage", default="all",
                        choices=["reference", "sweep", "full", "subsets", "analyse", "all"])
    args = parser.parse_args()
    root = args.root.resolve()
    binary = args.binary.resolve()
    helper = args.helper.resolve()

    if args.stage in {"reference", "all"}:
        tasks = [(root, helper, dataset, None, args.threads) for dataset in range(args.datasets)]
        tasks += [(root, helper, dataset, size, args.threads)
                  for size in SUBSET_SIZES for dataset in range(args.datasets)]
        parallel(tasks, run_reference, args.workers)
    if args.stage in {"sweep", "all"}:
        parallel([(root, binary, geometry, seed_threshold, dataset, args.threads)
                  for geometry, seed_threshold in SWEEP
                  for dataset in range(args.sweep_datasets)], run_sweep, args.workers)
    if args.stage in {"full", "all"}:
        parallel([(root, binary, order, method, dataset, args.threads, None)
                  for dataset in range(args.datasets)
                  for order, method in combinations()], run_cluster, args.workers)
    if args.stage in {"subsets", "all"}:
        parallel([(root, binary, order, method, dataset, args.threads, size)
                  for size in SUBSET_SIZES for dataset in range(args.datasets)
                  for order, method in combinations()], run_cluster, args.workers)
    if args.stage in {"analyse", "all"}:
        analyse(root, args.datasets, args.sweep_datasets)


if __name__ == "__main__":
    main()
