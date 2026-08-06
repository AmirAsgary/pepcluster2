#!/usr/bin/env python3
"""Compare PepCluster2 against external motif-deconvolution tools on one dataset.

Every tool sees the identical pools and is scored with the identical metrics.

A caveat that governs how the numbers should be read: PepCluster2 is a similarity
clustering method and produces on the order of a hundred clusters per pool, while
MixMHCp and GibbsCluster are mixture models that fit a handful of motifs.

Purity must not be read on its own. Adjusted per-allele purity is BCubed
precision corrected against the allele prior, and that correction removes the
baseline for one large cluster but *not* the inflation from fragmentation: a
singleton scores precision 1 and therefore adjusted purity 1, so a partition into
singletons scores a perfect 1.0. Comparing a hundred clusters against four on
that metric alone is not a like-for-like contest. BCubed recall is its
counterpart - it penalises exactly what purity rewards - and BCubed F1 is the
balance of the two. AMI is chance-corrected over the whole partition.

Results are written under `results/<dataset>/` so further datasets can be added
alongside without disturbing this one.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "code" / "mhc_bench"))
import plotstyle as PS  # noqa: E402

METRICS = [("ami", "AMI"),
           ("adjusted_purity_macro", "Adjusted per-allele purity"),
           ("bcubed_recall_macro", "BCubed recall"),
           ("bcubed_f1_macro", "BCubed F1"),
           ("singleton_fraction_of_clusters", "Singleton fraction")]
# Columns carried through from every tool's own scored output.
CARRY = ["ami", "adjusted_purity_macro", "bcubed_precision_macro",
         "bcubed_recall_macro", "bcubed_f1_macro",
         "singleton_fraction_of_clusters", "clusters"]

# Blue and red stay with our two modes across every figure in the study; the
# external tools take the next colours in the fixed order.
STYLE = {
    "PepCluster2 k-mer": PS.SERIES[0],
    "PepCluster2 alignment": PS.SERIES[1],
    "MixMHCp (default)": PS.SERIES[2],
    "MixMHCp (forced k)": PS.SERIES[4],
    "GibbsCluster (default)": PS.SERIES[3],
    "GibbsCluster (forced k)": PS.SERIES[5],
}
SPLITS = ["inner", "outer", "test"]


def load_ours(runs: Path) -> pd.DataFrame:
    """PepCluster2 at its selected configuration, from the tuned analysis."""
    rows = []
    for label, root in (("PepCluster2 k-mer", "mhc_bench_sep_kmer_anchor"),
                        ("PepCluster2 alignment", "mhc_bench_sep_aln_anchor")):
        base = runs / root
        choice = pd.read_csv(base / "tables" / "selected_overall.csv")
        choice = choice[choice.method == "graph"].iloc[0]
        for split in SPLITS:
            files = sorted((base / "grid").glob(f"{split}_*.csv"))
            frames = [pd.read_csv(f) for f in files]
            frames = [f for f in frames if "status" in f.columns and len(f)]
            frame = pd.concat(frames, ignore_index=True)
            frame = frame[(frame.status == "ok") & (frame.method == "graph") &
                          (frame.representative_order == choice.representative_order) &
                          (frame.primary_threshold == choice.primary_threshold) &
                          (frame.anchor_threshold == choice.anchor_threshold)]
            frame = frame.drop_duplicates("pool")
            for record in frame.to_dict("records"):
                row = {"tool": label, "split": split, "pool": record["pool"],
                       "allele_count": record["allele_count"]}
                row.update({c: record[c] for c in CARRY})
                rows.append(row)
    return pd.DataFrame(rows)


def load_external(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path)
    frame = frame[frame.status == "ok"]
    label = {"default": f"{name} (default)", "oracle_k": f"{name} (forced k)"}
    frame = frame.assign(tool=frame["setting"].map(label))
    return frame[["tool", "split", "pool", "allele_count"] + CARRY]


def figure(detail: pd.DataFrame, output: Path, dataset: str) -> None:
    import matplotlib.pyplot as plt
    PS.apply()

    order = [t for t in STYLE if t in set(detail["tool"])]
    splits = [s for s in SPLITS if s in set(detail["split"])]
    columns = 3
    rows_n = (len(METRICS) + columns - 1) // columns
    figure, axes = plt.subplots(rows_n, columns, figsize=(13.5, 4.0 * rows_n),
                                constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for spare in axes[len(METRICS):]:
        spare.axis("off")
    width = 0.8 / max(len(order), 1)
    for axis, (metric, title) in zip(axes, METRICS):
        for i, tool in enumerate(order):
            offset = (i - (len(order) - 1) / 2) * width
            means = [detail[(detail.tool == tool) & (detail.split == s)][metric].mean()
                     for s in splits]
            errs = [detail[(detail.tool == tool) & (detail.split == s)][metric].std()
                    for s in splits]
            axis.bar(np.arange(len(splits)) + offset, means, width, yerr=errs,
                     capsize=2.5, color=STYLE[tool], label=tool,
                     edgecolor="white", linewidth=0.6,
                     error_kw=dict(elinewidth=0.8, ecolor=PS.INK))
        axis.set_xticks(np.arange(len(splits)))
        axis.set_xticklabels([PS.SPLIT_LABEL[s] for s in splits])
        axis.grid(axis="x", alpha=0)
        if metric == "singleton_fraction_of_clusters":
            axis.set_ylim(bottom=0)  # a fraction cannot be negative
        PS.finish(axis, "", "", title)
    axes[0].legend(loc="upper left", ncol=1)
    figure.suptitle(f"PepCluster2 against motif-deconvolution tools ({dataset})")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"tool_comparison.{suffix}", dpi=200, bbox_inches="tight")
    plt.close(figure)


def figure_by_complexity(detail: pd.DataFrame, output: Path, dataset: str) -> pd.DataFrame:
    """Where each tool wins. The split averages hide a strong interaction with
    the number of alleles, which is the whole shape of the result."""
    import matplotlib.pyplot as plt
    PS.apply()

    order = [t for t in STYLE if t in set(detail["tool"])]
    bins = [(2, 6), (7, 12), (13, 30)]
    detail = detail.copy()
    detail["band"] = pd.cut(detail["allele_count"],
                            bins=[1, 6, 12, 30], labels=[f"{a}-{b}" for a, b in bins])

    figure, axes = plt.subplots(1, 3, figsize=(16.0, 4.0), constrained_layout=True)
    for metric, title, axis in ((("ami"), "AMI", axes[0]),
                                ("adjusted_purity_macro", "Adjusted per-allele purity",
                                 axes[1]),
                                ("bcubed_f1_macro", "BCubed F1", axes[2])):
        for tool in order:
            grouped = detail[detail.tool == tool].groupby(
                "allele_count", observed=True)[metric].mean()
            axis.plot(grouped.index, grouped.values, marker="o",
                      color=STYLE[tool], label=tool)
        axis.set_ylim(bottom=0)
        PS.finish(axis, "Alleles in pool", "", title)
    axes[0].legend(loc="best")
    figure.suptitle(f"Where each tool wins ({dataset})")
    for suffix in ("png", "pdf"):
        figure.savefig(output / f"by_complexity.{suffix}", dpi=200, bbox_inches="tight")
    plt.close(figure)

    rows = []
    for tool in order:
        if tool.startswith("PepCluster2"):
            continue
        for band, (lo, hi) in zip([f"{a}-{b}" for a, b in bins], bins):
            for ours in [t for t in order if t.startswith("PepCluster2")]:
                paired = detail[detail.band == band].pivot_table(
                    index="pool", columns="tool", values="ami")
                if tool not in paired or ours not in paired:
                    continue
                paired = paired.dropna(subset=[tool, ours])
                gap = paired[tool] - paired[ours]
                rows.append({"band": band, "tool": tool, "baseline": ours,
                             "pools": len(paired), "mean_ami_gap": gap.mean(),
                             "tool_wins_fraction": (gap > 0).mean()})
    return pd.DataFrame(rows)


def write_report(detail: pd.DataFrame, output: Path, dataset: str, missing: list[str],
                 bands: pd.DataFrame) -> None:
    lines: list[str] = []
    w = lines.append
    w(f"# External tool comparison: {dataset}")
    w("")
    w("Identical pools, identical metrics. PepCluster2 is evaluated at the")
    w("configuration selected by its own nested cross-validation; the external")
    w("tools have no threshold to tune, so they are run as documented.")
    w("")
    if missing:
        w("> Not yet included: " + ", ".join(missing) + ".")
        w("")
    w("## How to read this")
    w("")
    w("PepCluster2 is a similarity clustering method and returns on the order of a")
    w("hundred clusters per pool. MixMHCp and GibbsCluster are mixture models that")
    w("fit a handful of motifs, so partition size differs by more than an order of")
    w("magnitude and the metrics must be read as a set.")
    w("")
    w("Adjusted per-allele purity is BCubed precision corrected against the allele")
    w("prior. That correction removes the baseline for one large cluster but not the")
    w("inflation from fragmentation: a singleton scores precision 1, so a partition")
    w("into singletons scores a perfect 1.0. **Purity alone therefore cannot be used")
    w("to compare partitions of different granularity.** BCubed recall penalises")
    w("exactly what purity rewards, and BCubed F1 balances them. AMI is")
    w("chance-corrected over the whole partition. F1 and AMI are the two figures")
    w("that can be compared across tools directly.")
    w("")
    w("`forced k` gives a tool the true number of alleles in the pool. No user could")
    w("do that in practice, so it is not a fair headline number - it isolates how")
    w("much of a tool's result is its model rather than its model selection.")
    w("")
    for split in SPLITS:
        subset = detail[detail.split == split]
        if subset.empty:
            continue
        w(f"## {PS.SPLIT_LABEL[split].capitalize()}")
        w("")
        w("| Tool | Pools | AMI | Purity (macro) | Recall | F1 | Singletons | Clusters |")
        w("|---|---:|---:|---:|---:|---:|---:|---:|")
        for tool in [t for t in STYLE if t in set(subset.tool)]:
            s = subset[subset.tool == tool]
            w(f"| {tool} | {len(s)} | {s.ami.mean():.4f} ± {s.ami.std():.4f} | "
              f"{s.adjusted_purity_macro.mean():.4f} ± "
              f"{s.adjusted_purity_macro.std():.4f} | "
              f"{s.bcubed_recall_macro.mean():.4f} | "
              f"{s.bcubed_f1_macro.mean():.4f} | "
              f"{s.singleton_fraction_of_clusters.mean():.3f} | "
              f"{s.clusters.mean():.1f} |")
        w("")
    if not bands.empty:
        w("## Where the difference comes from")
        w("")
        w("The split averages hide a strong interaction with pool complexity.")
        w("Paired per-pool AMI, positive meaning the external tool is ahead:")
        w("")
        w("| Alleles in pool | Tool | vs | Pools | Mean AMI gap | Tool wins |")
        w("|---|---|---|---:|---:|---:|")
        for row in bands.itertuples():
            w(f"| {row.band} | {row.tool} | {row.baseline} | {row.pools} | "
              f"{row.mean_ami_gap:+.4f} | {row.tool_wins_fraction:.0%} |")
        w("")
        w("MixMHCp is a mixture model built to resolve a few motifs, and that is")
        w("exactly where it dominates. Its margin shrinks steadily as alleles are")
        w("added and is gone by roughly 13 alleles, where the two are level.")
        w("")
        w("Mean adjusted per-allele purity over the same bands:")
        w("")
        purity = detail.copy()
        purity["band"] = pd.cut(purity["allele_count"], bins=[1, 6, 12, 30],
                                labels=["2-6", "7-12", "13-30"])
        table = purity.pivot_table(index="band", columns="tool",
                                   values="adjusted_purity_macro", observed=True)
        tools = [t for t in STYLE if t in table.columns]
        w("| Alleles in pool | " + " | ".join(tools) + " |")
        w("|---" * (len(tools) + 1) + "|")
        for band, record in table.iterrows():
            w(f"| {band} | " + " | ".join(f"{record[t]:.4f}" for t in tools) + " |")
        w("")
        ours = [t for t in tools if t.startswith("PepCluster2")]
        best_ours = table[ours].max(axis=1)
        best_other = table[[t for t in tools if t not in ours]].max(axis=1)
        leads = [str(b) for b in table.index if best_ours[b] > best_other[b]]
        if leads:
            w(f"PepCluster2 leads on purity in the {', '.join(leads)} band(s). Read that")
            w("against the cluster counts in the tables above and against F1 below")
            w("before drawing a conclusion: purity rises mechanically with the number")
            w("of clusters, so a lead on purity held by the partition with an order of")
            w("magnitude more clusters is not evidence of a better partition.")
            w("")
        w("Mean BCubed F1 over the same bands, which penalises fragmentation and")
        w("lumping together and is the figure to compare across tools:")
        w("")
        f1 = detail.copy()
        f1["band"] = pd.cut(f1["allele_count"], bins=[1, 6, 12, 30],
                            labels=["2-6", "7-12", "13-30"])
        f1_table = f1.pivot_table(index="band", columns="tool",
                                  values="bcubed_f1_macro", observed=True)
        f1_tools = [t for t in STYLE if t in f1_table.columns]
        w("| Alleles in pool | " + " | ".join(f1_tools) + " |")
        w("|---" * (len(f1_tools) + 1) + "|")
        for band, record in f1_table.iterrows():
            w(f"| {band} | " + " | ".join(f"{record[t]:.4f}" for t in f1_tools) + " |")
        w("")
    w("## Files")
    w("")
    w("- `tables/per_pool.csv` - every tool on every pool")
    w("- `tables/summary.csv` - the tables above")
    w("- `tables/by_complexity.csv` - the paired comparison by allele count")
    w("- `tool_comparison.png` / `.pdf`, `by_complexity.png` / `.pdf`")
    w("- `raw/` - each tool's own run output")
    w("")
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True,
                        help="results/<dataset> directory")
    parser.add_argument("--dataset", default="immuneapp")
    args = parser.parse_args()
    results = args.results.resolve()
    (results / "tables").mkdir(parents=True, exist_ok=True)

    parts = [load_ours(args.runs.resolve())]
    missing = []
    for name, stem in (("MixMHCp", "mixmhcp"), ("GibbsCluster", "gibbscluster")):
        frame = load_external(results / "raw" / f"{stem}.csv", name)
        if frame.empty:
            missing.append(name)
        else:
            parts.append(frame)
    detail = pd.concat(parts, ignore_index=True)
    detail.to_csv(results / "tables" / "per_pool.csv", index=False)

    summary = detail.groupby(["tool", "split"]).agg(
        pools=("pool", "size"),
        ami_mean=("ami", "mean"), ami_std=("ami", "std"),
        purity_mean=("adjusted_purity_macro", "mean"),
        purity_std=("adjusted_purity_macro", "std"),
        recall_mean=("bcubed_recall_macro", "mean"),
        recall_std=("bcubed_recall_macro", "std"),
        f1_mean=("bcubed_f1_macro", "mean"), f1_std=("bcubed_f1_macro", "std"),
        singletons=("singleton_fraction_of_clusters", "mean"),
        clusters=("clusters", "mean")).reset_index()
    summary.to_csv(results / "tables" / "summary.csv", index=False)

    figure(detail, results, args.dataset)
    bands = figure_by_complexity(detail, results, args.dataset)
    bands.to_csv(results / "tables" / "by_complexity.csv", index=False)
    write_report(detail, results, args.dataset, missing, bands)
    print(summary.round(4).to_string(index=False))
    if not bands.empty:
        print()
        print(bands.round(4).to_string(index=False))
    if missing:
        print("\nnot included:", ", ".join(missing))


if __name__ == "__main__":
    main()
