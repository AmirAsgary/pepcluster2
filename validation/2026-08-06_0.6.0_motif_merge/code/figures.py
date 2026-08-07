#!/usr/bin/env python3
"""Benchmark figures and tables across all tools, splits and pool properties.

Builds one per-pool table holding every tool at every pool, then draws:

  fig1_protocol            how the nested selection was run and what it chose
  fig2_granularity         cluster count and mean cluster size against pool size
                           and against allele count
  fig3_performance_alleles five metrics against the number of alleles
  fig4_performance_size    five metrics against the number of peptides

Colour is fixed by tool family throughout: black for PepCluster2, red for
MixMHCp, blue for GibbsCluster. Solid is each tool as documented; dashed is the
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

# black = ours, red = MixMHCp, blue = GibbsCluster.
STYLE = {
    "PepCluster2 + motif":     dict(color="#000000", ls="-",  marker="o", lw=2.0, zorder=6),
    "PepCluster2 + motif (given k)": dict(color="#000000", ls=(0, (4, 1.5)),
                                          marker="D", lw=1.6, alpha=0.85, zorder=5),
    "PepCluster2 + EM only":   dict(color="#555555", ls="-.", marker="v", lw=1.6,
                                    zorder=4),
    "PepCluster2 (similarity)": dict(color="#000000", ls=":",  marker="s", lw=1.4, zorder=4),
    "MixMHCp (default)":       dict(color="#d62728", ls="-",  marker="o", lw=1.6, zorder=3),
    "MixMHCp (forced k)":      dict(color="#d62728", ls=(0, (5, 2)), marker="^", lw=1.6, alpha=0.85, zorder=3),
    "GibbsCluster (default)":  dict(color="#1f77b4", ls="-",  marker="o", lw=1.6, zorder=2),
    "GibbsCluster (forced k)": dict(color="#1f77b4", ls=(2.5, (5, 2)), marker="v", lw=1.6, alpha=0.85, zorder=2),
}
ORDER = list(STYLE)
SPLITS = ["inner", "outer", "test"]
SPLIT_LABEL = {"inner": "tuning folds", "outer": "held-out alleles",
               "test": "independent test"}


VARIANT_LABEL = {"em_only": "PepCluster2 + EM only",
                 "forced_k": "PepCluster2 + motif (given k)"}


def build_table(bench: Path, study: Path, grid: Path, selected: Path,
                variants=()) -> pd.DataFrame:
    """Every tool on every pool, with the pool's own properties attached."""
    per_pool = pd.read_csv(bench / "results/immuneapp/tables/per_pool.csv")
    per_pool = per_pool[per_pool.tool != "PepCluster2 alignment"].copy()
    per_pool["tool"] = per_pool.tool.replace(
        {"PepCluster2 k-mer": "PepCluster2 (similarity)"})

    # our motif layer at the configuration nested selection chose
    choice = pd.read_csv(selected).iloc[0]
    shards = pd.concat([pd.read_csv(f) for f in sorted(grid.glob("motif_shard_*.csv"))],
                       ignore_index=True)
    shards = shards[shards.status == "ok"].copy()
    shards["em_concentration"] = shards.em_concentration.fillna(-1.0)
    mask = ((shards.merge_concentration == choice.merge_concentration) &
            (shards.merge_threshold == choice.merge_threshold) &
            (shards.em.astype(bool) == bool(choice.em)) &
            (shards.em_concentration == choice.em_concentration))
    ours = shards[mask].drop_duplicates("pool").copy()
    ours["tool"] = "PepCluster2 + motif"
    keep = ["tool", "split", "pool", "allele_count"] + [m for m, _ in METRICS] + \
           ["singleton_fraction_of_clusters", "clusters"]
    parts = [per_pool[keep], ours[keep]]
    for path in variants:
        path = Path(path)
        if not path.exists():
            continue
        name = path.stem.replace("variant_", "")
        extra = pd.read_csv(path)
        extra = extra[extra.status == "ok"].copy()
        extra["tool"] = VARIANT_LABEL.get(name, f"PepCluster2 ({name})")
        parts.append(extra[keep])
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
        for tool in ORDER:
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
                      marker=marker,
                      markersize=4, capsize=2, label=tool, **style)
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bench", type=Path, required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True)
    parser.add_argument("--variant", type=Path, nargs="*", default=(),
                        help="variant_<name>.csv files, each added as an arm")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})

    frame = build_table(args.bench, args.study, args.grid, args.selected,
                        args.variant)
    frame.to_csv(args.out / "per_pool_all_tools.csv", index=False)
    print(f"per-pool table: {len(frame)} rows, {frame.tool.nunique()} tools, "
          f"{frame.pool.nunique()} pools")

    # Headline tables, generated here so they cannot drift from the figures.
    def block(subset):
        rows = []
        for tool in ORDER:
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

    test = frame[frame.split == "test"]
    size_bins = [0, 1500, 2200, 3200, 5000, 12000]
    allele_bins = [3, 5, 7, 9, 12, 16, 20]

    # ---- fig 2: granularity -------------------------------------------
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
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
    axes[0][0].legend(fontsize=7, loc="upper left")
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
        figure, axes = plt.subplots(1, 5, figsize=(19, 3.6), constrained_layout=True)
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
        axes[0].legend(fontsize=7, loc="upper right")
        figure.suptitle(f"Performance against {xlabel.lower()} "
                        "(independent test pools, mean ± s.e.m.)")
        for suffix in ("png", "pdf"):
            figure.savefig(args.out / f"{name}.{suffix}", bbox_inches="tight")
        plt.close(figure)

    # ---- fig 5: every metric on every split ---------------------------
    # Independent test split only: the tuning folds informed the choice of
    # configuration, so showing them beside the test set invites reading a number
    # that is not held out.
    figure, axes = plt.subplots(1, 5, figsize=(17, 3.8), constrained_layout=True)
    width = 0.8 / len(ORDER)
    for letter, (axis, (metric, title)) in zip("abcde", zip(axes, METRICS)):
        panel_rows = []
        for i, tool in enumerate(ORDER):
            offset = (i - (len(ORDER) - 1) / 2) * width
            means, errs = [], []
            for split in ["test"]:
                sub = frame[(frame.tool == tool) & (frame.split == split)][metric]
                means.append(sub.mean())
                errs.append(sub.std(ddof=1) / max(np.sqrt(len(sub)), 1))
            panel_rows += [dict(tool=tool, split=s, mean=m, sem=e)
                           for s, m, e in zip(["test"], means, errs)]
            style = STYLE[tool]
            axis.bar(np.arange(1) + offset, means, width, yerr=errs,
                     capsize=2, color=style["color"], label=tool,
                     alpha=1.0 if style["ls"] == "-" else 0.45,
                     hatch="" if style["ls"] == "-" else "///",
                     edgecolor="white", linewidth=0.5)
        pd.DataFrame(panel_rows).to_csv(
            args.out / f"fig5{letter}_{metric}_test.csv", index=False)
        axis.set_xticks([])
        axis.set_title(title)
        axis.set_ylim(0, 1.0)
        axis.grid(axis="y", alpha=0.25, lw=0.5)
    axes[0].set_ylabel("Score")
    axes[0].legend(fontsize=6.5, loc="upper right")
    figure.suptitle("All metrics, independent test pools (mean ± s.e.m.); "
                    "hatched = variant given the true allele count")
    for suffix in ("png", "pdf"):
        figure.savefig(args.out / f"fig5_all_splits.{suffix}", bbox_inches="tight")
    plt.close(figure)

    print("wrote figures and tables to", args.out)


if __name__ == "__main__":
    main()
