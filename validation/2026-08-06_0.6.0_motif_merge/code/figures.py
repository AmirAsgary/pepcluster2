#!/usr/bin/env python3
"""Benchmark figures and tables across all tools, splits and pool properties.

Builds one per-pool table holding every tool at every pool, then draws:

  fig2_granularity         cluster count and mean cluster size against pool size
                           and against allele count
  fig3_performance_alleles five metrics against the number of alleles
  fig4_performance_size    five metrics against the number of peptides
  fig5_all_splits          every metric on the independent test split

Seven PepCluster2 modes are carried through every figure, so the reader can see
what each stage of the pipeline buys:

  PepCluster2                          similarity clustering alone
  PepCluster2 + merge                  agglomerative merge of the clusters
  PepCluster2 + merge (given K)        merge continued down to exactly K groups
  PepCluster2 + merge + EM             EM seeded from the merged groups
  PepCluster2 + merge + EM (given K)   EM seeded from K merged groups
  PepCluster2 + EM                     EM seeded from the raw clusters
  PepCluster2 + EM (given K)           EM seeded from K of the raw clusters

Colour encodes the pipeline, line style encodes whether the allele count was
supplied. Black is clustering alone; grey is the merge stage; purple and violet
are the four modes that run EM. MixMHCp is red and GibbsCluster is blue
throughout. Solid is each tool as a user would actually run it; dashed is the
variant handed the true allele count, which no user could supply in practice.

Correlations use Spearman's rho. The quantities are counts spanning more than an
order of magnitude with no reason to be linear or normal, so a rank statistic is
the honest choice. Every test in the family is corrected together by
Benjamini-Hochberg.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

METRICS = [("ami", "AMI"),
           ("adjusted_purity_macro", "Adjusted purity"),
           ("bcubed_precision_macro", "BCubed precision"),
           ("bcubed_recall_macro", "BCubed recall"),
           ("bcubed_f1_macro", "BCubed F1")]

BLACK, GREY, PURPLE, VIOLET = "#000000", "#808080", "#6A1B9A", "#C77DFF"
RED, BLUE = "#D62728", "#1F77B4"


def dash(phase):
    """A 5-2 dash with a per-series phase offset.

    Every given-K arm returns exactly K clusters, so in the granularity panels
    all five dashed lines sit on the same curve. Offsetting the phase makes them
    interleave instead of hiding one another.
    """
    return (phase, (5, 2))

# variant file stem -> display name. The order here is the order everywhere.
MODES = [
    ("similarity", "PepCluster2"),
    ("merge",      "PepCluster2 + merge"),
    ("merge_k",    "PepCluster2 + merge (given K)"),
    ("merge_em",   "PepCluster2 + merge + EM"),
    ("merge_em_k", "PepCluster2 + merge + EM (given K)"),
    ("em_only",    "PepCluster2 + EM"),
    ("em_only_k",  "PepCluster2 + EM (given K)"),
]
# Reported in the tables, kept out of the figures: a second merge-seeded given-K
# route that seeds from groups merged at the evidence threshold rather than from
# groups merged down to exactly K. The nested pair merge_k -> merge_em_k is the
# one carried in the figures, because it isolates EM against an identical merge.
SUPPLEMENTARY = [("merge_em_k_merged",
                  "PepCluster2 + merge + EM (given K, threshold-merged)")]

# Comparators, renamed from "forced k" so the whole figure uses one phrase.
COMPARATORS = {"MixMHCp (default)": "MixMHCp",
               "MixMHCp (forced k)": "MixMHCp (given K)",
               "GibbsCluster (default)": "GibbsCluster",
               "GibbsCluster (forced k)": "GibbsCluster (given K)"}

STYLE = {
    "PepCluster2":                        dict(color=BLACK,  ls="-",  marker="o", lw=2.2, zorder=11),
    "PepCluster2 + merge":                dict(color=GREY,   ls="-",  marker="s", lw=2.0, zorder=10),
    "PepCluster2 + merge (given K)":      dict(color=GREY,   ls=dash(0.0), marker="s", lw=2.0, zorder=9),
    "PepCluster2 + merge + EM":           dict(color=PURPLE, ls="-",  marker="o", lw=2.0, zorder=8),
    "PepCluster2 + merge + EM (given K)": dict(color=PURPLE, ls=dash(1.4), marker="o", lw=2.0, zorder=7),
    "PepCluster2 + EM":                   dict(color=VIOLET, ls="-",  marker="D", lw=2.0, zorder=6),
    "PepCluster2 + EM (given K)":         dict(color=VIOLET, ls=dash(2.8), marker="D", lw=2.0, zorder=5),
    "MixMHCp":                            dict(color=RED,    ls="-",  marker="^", lw=1.8, zorder=4),
    "MixMHCp (given K)":                  dict(color=RED,    ls=dash(4.2), marker="^", lw=1.8, zorder=3),
    "GibbsCluster":                       dict(color=BLUE,   ls="-",  marker="v", lw=1.8, zorder=2),
    "GibbsCluster (given K)":             dict(color=BLUE,   ls=dash(5.6), marker="v", lw=1.8, zorder=1),
}
ORDER = [label for _, label in MODES] + list(COMPARATORS.values())
TABLE_ORDER = ([label for _, label in MODES] +
               [label for _, label in SUPPLEMENTARY] + list(COMPARATORS.values()))
SPLITS = ["inner", "outer", "test"]
SPLIT_LABEL = {"inner": "tuning folds", "outer": "held-out alleles",
               "test": "independent test"}

KEEP = ["tool", "split", "pool", "allele_count"] + [m for m, _ in METRICS] + \
       ["singleton_fraction_of_clusters", "clusters"]


def build_table(bench: Path, study: Path, results: Path) -> pd.DataFrame:
    """Every tool on every pool, with the pool's own properties attached.

    The PepCluster2 arms all come from one variant sweep, so the only thing that
    differs between them is the motif stage under test. The comparators come from
    the benchmark run. The `PepCluster2 k-mer` row of the benchmark table is the
    same numbers as the `similarity` variant and is dropped to avoid a duplicate.
    """
    parts = []
    for stem, label in MODES + SUPPLEMENTARY:
        path = results / f"variant_{stem}.csv"
        if not path.exists():
            raise SystemExit(f"missing variant file: {path}")
        frame = pd.read_csv(path)
        frame = frame.loc[:, ~frame.columns.duplicated()]
        bad = (frame.status != "ok").sum()
        if bad:
            print(f"warning: {stem} has {bad} non-ok pools, dropped")
        frame = frame[frame.status == "ok"].copy()
        frame["tool"] = label
        parts.append(frame[KEEP])

    per_pool = pd.read_csv(bench / "results/immuneapp/tables/per_pool.csv")
    per_pool = per_pool[per_pool.tool.isin(COMPARATORS)].copy()
    per_pool["tool"] = per_pool.tool.map(COMPARATORS)
    parts.append(per_pool[KEEP])

    frame = pd.concat(parts, ignore_index=True)
    manifest = pd.read_csv(study / "runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv")
    frame = frame.merge(manifest[["pool", "peptides"]], on="pool", how="left")
    frame["mean_cluster_size"] = frame.peptides / frame.clusters
    return frame


def correlations(frame: pd.DataFrame, out: Path) -> pd.DataFrame:
    """Spearman rho for every (tool, split, response, predictor), BH-corrected.

    One family, corrected together: testing many pairs and reporting the small
    p-values without adjustment is how spurious associations get published.
    """
    rows = []
    responses = ["clusters", "mean_cluster_size"] + [m for m, _ in METRICS]
    for split in SPLITS:
        for tool in TABLE_ORDER:
            sub = frame[(frame.split == split) & (frame.tool == tool)]
            if len(sub) < 8:
                continue
            for response in responses:
                for predictor in ("peptides", "allele_count"):
                    x = sub[predictor].to_numpy(float)
                    y = sub[response].to_numpy(float)
                    ok = np.isfinite(x) & np.isfinite(y)
                    if ok.sum() < 8 or len(np.unique(x[ok])) < 3:
                        continue
                    rho, p = stats.spearmanr(x[ok], y[ok])
                    rows.append(dict(split=split, tool=tool, response=response,
                                     predictor=predictor, n=int(ok.sum()),
                                     spearman_rho=rho, p_raw=p))
    table = pd.DataFrame(rows)
    if not table.empty:
        table["p_bh"] = stats.false_discovery_control(table.p_raw.to_numpy())
        table["significant_bh_0.05"] = table.p_bh < 0.05
    table.to_csv(out, index=False)
    return table


def _decade_ticks(caxis):
    """Ticks at 1, 2 and 5 of every decade, labelled as plain numbers.

    A log axis labelled only at 10, 100, 1000 makes it hard to read off a value
    that sits between decades, which is most of them here.
    """
    from matplotlib import ticker
    caxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0),
                                              numticks=30))
    caxis.set_minor_locator(ticker.LogLocator(base=10.0,
                                              subs=tuple(np.arange(1, 10) * 0.1),
                                              numticks=100))
    caxis.set_major_formatter(ticker.FuncFormatter(
        lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
    caxis.set_minor_formatter(ticker.NullFormatter())


def _panel(axis, frame, response, xkey, bins, logx, logy=False):
    """Binned mean with standard error, one line per tool.

    Returns the plotted points so the panel can be written out verbatim: the
    figure and its CSV are then guaranteed to be the same numbers.
    """
    collected = []
    for tool in ORDER:
        sub = frame[frame.tool == tool]
        if sub.empty:
            continue
        cut = pd.cut(sub[xkey], bins=bins)
        grouped = sub.groupby(cut, observed=True).agg(
            x=(xkey, "mean"), y=(response, "mean"),
            err=(response, lambda v: v.std(ddof=1) / max(np.sqrt(len(v)), 1)),
            n=(response, "size")).dropna(subset=["x", "y"])
        grouped = grouped[grouped["n"] >= 3]
        if grouped.empty:
            continue
        style = dict(STYLE[tool])
        marker = style.pop("marker")
        axis.errorbar(grouped["x"], grouped["y"], yerr=grouped["err"].fillna(0),
                      marker=marker, markersize=6, capsize=2.5, label=tool,
                      **style)
        collected.append(grouped.reset_index(drop=True).assign(
            tool=tool, response=response, predictor=xkey))
    if logx:
        axis.set_xscale("log")
        _decade_ticks(axis.xaxis)
    if logy:
        axis.set_yscale("log")
        _decade_ticks(axis.yaxis)
    axis.grid(alpha=0.25, lw=0.5, which="both")
    return (pd.concat(collected, ignore_index=True) if collected
            else pd.DataFrame())


def _outside_legend(figure, axis):
    """One legend for the whole figure, to the right of the axes.

    Placed outside so no series is ever hidden behind the key, and ordered by
    ORDER rather than by draw order so the pipeline reads top to bottom.
    """
    handles, labels = axis.get_legend_handles_labels()
    index = {label: i for i, label in enumerate(labels)}
    picked = [(handles[index[t]], t) for t in ORDER if t in index]
    figure.legend([h for h, _ in picked], [t for _, t in picked],
                  loc="outside right center", frameon=False,
                  fontsize=plt.rcParams["font.size"] - 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", type=Path, required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True,
                        help="directory holding variant_<mode>.csv")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    global plt
    import matplotlib.pyplot as plt
    # Every size is the old value plus four points.
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 15,
                         "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12, "figure.titlesize": 16,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.dpi": 150})

    frame = build_table(args.bench, args.study, args.results)
    frame.to_csv(args.out / "per_pool_all_tools.csv", index=False)
    print(f"per-pool table: {len(frame)} rows, {frame.tool.nunique()} tools, "
          f"{frame.pool.nunique()} pools")

    # Headline tables, generated here so they cannot drift from the figures.
    def block(subset):
        rows = []
        for tool in TABLE_ORDER:
            s = subset[subset.tool == tool]
            if s.empty:
                continue
            row = {"tool": tool, "pools": len(s)}
            for metric, label in METRICS:
                row[label] = round(s[metric].mean(), 4)
                row[f"{label}_sd"] = round(s[metric].std(ddof=1), 4)
            row["clusters"] = round(s.clusters.mean(), 2)
            row["mean_cluster_size"] = round(s.mean_cluster_size.mean(), 1)
            rows.append(row)
        return pd.DataFrame(rows)

    test = frame[frame.split == "test"]
    block(test).to_csv(args.out / "table1_benchmark_test.csv", index=False)
    banded = test.assign(allele_band=pd.cut(test.allele_count, bins=[3, 6, 12, 20],
                                            labels=["4-6", "7-12", "13-20"]))
    pd.concat([block(g).assign(allele_band=b)
               for b, g in banded.groupby("allele_band", observed=True)],
              ignore_index=True).to_csv(
        args.out / "table2_test_by_allele_band.csv", index=False)
    pd.concat([block(frame[frame.split == s]).assign(split=s) for s in SPLITS],
              ignore_index=True).to_csv(args.out / "table3_all_splits.csv",
                                        index=False)

    corr = correlations(frame, args.out / "correlations_spearman.csv")
    print(f"correlations: {len(corr)} tests, "
          f"{int(corr['significant_bh_0.05'].sum())} significant after BH")

    bands = pd.cut(frame.allele_count, bins=[1, 6, 12, 30],
                   labels=["2-6", "7-12", "13-30"])
    by_band = frame.assign(band=bands).groupby(["split", "band", "tool"],
                                               observed=True).agg(
        pools=("pool", "size"),
        **{m: (m, "mean") for m, _ in METRICS},
        clusters=("clusters", "mean")).reset_index()
    by_band.to_csv(args.out / "benchmark_by_allele_band.csv", index=False)

    size_bins = [0, 1500, 2200, 3200, 5000, 12000]
    allele_bins = [3, 5, 7, 9, 12, 16, 20]

    # ---- fig 2: granularity -------------------------------------------
    figure, axes = plt.subplots(2, 2, figsize=(15, 10.5), constrained_layout=True)
    panel_names = {("clusters", "peptides"): "fig2a_clusters_vs_pool_size",
                   ("clusters", "allele_count"): "fig2b_clusters_vs_allele_count",
                   ("mean_cluster_size", "peptides"):
                       "fig2c_mean_cluster_size_vs_pool_size",
                   ("mean_cluster_size", "allele_count"):
                       "fig2d_mean_cluster_size_vs_allele_count"}
    for row, (response, ylabel) in enumerate(
            [("clusters", "Clusters returned"),
             ("mean_cluster_size", "Mean peptides per cluster")]):
        data = _panel(axes[row][0], test, response, "peptides", size_bins, True,
                      logy=True)
        data.to_csv(args.out / f"{panel_names[(response, 'peptides')]}.csv",
                    index=False)
        axes[row][0].set_xlabel("Peptides in pool")
        axes[row][0].set_ylabel(ylabel)
        data = _panel(axes[row][1], test, response, "allele_count", allele_bins,
                      False, logy=True)
        data.to_csv(args.out / f"{panel_names[(response, 'allele_count')]}.csv",
                    index=False)
        axes[row][1].set_xlabel("Alleles in pool")
        axes[row][1].set_ylabel(ylabel)
    _outside_legend(figure, axes[0][0])
    figure.suptitle("Granularity of the returned partition (independent test pools)")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig2_granularity.{suffix}", bbox_inches="tight")
    plt.close(figure)

    # ---- fig 3 and 4: performance -------------------------------------
    for name, xkey, bins, xlabel, logx in (
            ("fig3_performance_alleles", "allele_count", allele_bins,
             "Alleles in pool", False),
            ("fig4_performance_size", "peptides", size_bins,
             "Peptides in pool", True)):
        figure, axes = plt.subplots(1, 5, figsize=(24, 5.2),
                                    constrained_layout=True)
        for letter, (axis, (metric, title)) in zip("abcde", zip(axes, METRICS)):
            data = _panel(axis, test, metric, xkey, bins, logx)
            number = name.split("_")[0]
            data.to_csv(args.out / f"{number}{letter}_{metric}_vs_"
                                   f"{'allele_count' if xkey == 'allele_count' else 'pool_size'}.csv",
                        index=False)
            axis.set_xlabel(xlabel)
            axis.set_title(title)
            axis.set_ylim(-0.02, 1.0)
        axes[0].set_ylabel("Score")
        _outside_legend(figure, axes[0])
        figure.suptitle(f"Performance against {xlabel.lower()} "
                        "(independent test pools, mean ± s.e.m.)")
        for suffix in ("png", "pdf"):
            figure.savefig(args.out / f"{name}.{suffix}", bbox_inches="tight")
        plt.close(figure)

    # ---- fig 5: every metric on the independent test split ------------
    # Independent test split only: the tuning folds informed the choice of
    # configuration, so showing them beside the test set invites reading a number
    # that is not held out. Stacked rows share one set of x labels, which is the
    # only layout in which eleven full tool names stay legible without a legend.
    figure, axes = plt.subplots(5, 1, figsize=(15, 19), sharex=True,
                                constrained_layout=True)
    positions = np.arange(len(ORDER))
    for letter, (axis, (metric, title)) in zip("abcde", zip(axes, METRICS)):
        panel_rows, means, errs, colors, hatches = [], [], [], [], []
        for tool in ORDER:
            sub = frame[(frame.tool == tool) & (frame.split == "test")][metric]
            mean = sub.mean()
            err = sub.std(ddof=1) / max(np.sqrt(len(sub)), 1)
            means.append(mean)
            errs.append(err)
            style = STYLE[tool]
            colors.append(style["color"])
            hatches.append("" if style["ls"] == "-" else "///")
            panel_rows.append(dict(tool=tool, split="test", mean=mean, sem=err,
                                   pools=len(sub)))
        bars = axis.bar(positions, means, 0.72, yerr=errs, capsize=3,
                        color=colors, edgecolor="white", linewidth=0.8)
        for bar, hatch in zip(bars, hatches):
            bar.set_hatch(hatch)
        pd.DataFrame(panel_rows).to_csv(
            args.out / f"fig5{letter}_{metric}_test.csv", index=False)
        axis.set_ylabel(title)
        axis.set_ylim(0, 1.0)
        axis.grid(axis="y", alpha=0.25, lw=0.5)
    axes[-1].set_xticks(positions)
    axes[-1].set_xticklabels(ORDER, rotation=45, ha="right",
                             rotation_mode="anchor")
    figure.suptitle("All metrics, independent test pools (mean ± s.e.m.); "
                    "hatched = variant given the true allele count")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig5_all_splits.{suffix}", bbox_inches="tight")
    plt.close(figure)

    print("wrote figures and tables to", args.out)


if __name__ == "__main__":
    main()
