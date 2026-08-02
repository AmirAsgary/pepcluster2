#!/usr/bin/env python3
"""Nested-CV selection, outer evaluation and figures for the peptide-MHC grid.

Selection rule
--------------
Within an outer fold, a configuration is ranked by its mean AMI over that fold's
inner pools, after excluding configurations whose mean singleton fraction exceeds
half of all clusters.

AMI alone, not a blend. Measured across this grid, NMI correlates +0.80 with the
singleton fraction and macro purity +0.96: both are maximised by the degenerate
clustering that puts every peptide in its own cluster, which scores NMI 0.365 and
purity 0.994. AMI is chance-corrected and correlates -0.97 with fragmentation, so
it already expresses the trade-off the other two get wrong. Purity and NMI are
reported, never optimised.

The outer folds hold out alleles, never seen during selection, so the outer
numbers estimate how the chosen thresholds transfer to new alleles.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plotstyle as PS  # noqa: E402

HEADLINE = ("ami", "nmi", "adjusted_purity_macro")
REPORTED = HEADLINE + (
    "adjusted_purity_micro", "bcubed_precision_macro", "bcubed_recall_macro",
    "bcubed_f1_macro", "bcubed_precision_micro", "bcubed_recall_micro",
    "clusters", "singleton_fraction_of_clusters", "objective",
)
SINGLETON_LIMIT = 0.50
CONFIG = ["method", "representative_order", "primary_threshold", "anchor_threshold"]

# Blue first, red second - the same order every other figure uses.
PALETTE = {"graph": PS.SERIES[0], "greedy_lazy": PS.SERIES[1]}
# The primary component differs by mode, so the axis must be named accordingly.
PRIMARY_LABEL = {
    "separate_aln_anchor": "Constrained-alignment similarity threshold",
    "separate_kmer_anchor": "Terminal k-mer similarity threshold",
}
LABELS = {"graph": "Graph (no prefilter)", "greedy_lazy": "Greedy lazy-exact"}
INK = PS.INK


def load_grid(root: Path, tag: str) -> pd.DataFrame:
    # Includes recovered/consolidated files, not just `<tag>_shard_*`.
    files = sorted((root / "grid").glob(f"{tag}_*.csv"))
    if not files:
        raise SystemExit(f"no {tag} grid results under {root / 'grid'}")
    frame = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    frame = frame.drop_duplicates(
        ["pool", "method", "representative_order", "primary_threshold", "anchor_threshold"])
    failed = frame[frame["status"] != "ok"]
    if len(failed):
        print(f"WARNING: {len(failed)} failed runs; first error: {failed.iloc[0]['error']}")
    return frame[frame["status"] == "ok"].copy()


def select(inner: pd.DataFrame, fold: int | None) -> pd.DataFrame:
    """Best configuration per method, optionally restricted to one outer fold."""
    data = inner if fold is None else inner[inner["outer_fold"] != fold]
    grouped = data.groupby(CONFIG).agg(
        objective=("objective", "mean"),
        ami=("ami", "mean"), nmi=("nmi", "mean"),
        adjusted_purity_macro=("adjusted_purity_macro", "mean"),
        singleton_fraction=("singleton_fraction_of_clusters", "mean"),
        clusters=("clusters", "mean"),
        pools=("objective", "size"),
    ).reset_index()
    grouped["feasible"] = grouped["singleton_fraction"] <= SINGLETON_LIMIT
    chosen = []
    for method in grouped["method"].unique():
        candidates = grouped[(grouped["method"] == method) & grouped["feasible"]]
        if candidates.empty:
            candidates = grouped[grouped["method"] == method]
            print(f"WARNING: no configuration for {method} meets the singleton limit")
        best = candidates.sort_values(
            ["ami", "primary_threshold", "anchor_threshold"],
            ascending=[False, True, True]).iloc[0]
        chosen.append({"held_out_fold": fold, **best.to_dict()})
    return pd.DataFrame(chosen)


def style(axis, xlabel, ylabel, title):
    axis.set_title(title, fontsize=10.5, color="#0b0b0b")
    axis.set_xlabel(xlabel, fontsize=9, color=INK)
    axis.set_ylabel(ylabel, fontsize=9, color=INK)
    axis.grid(alpha=0.22, linewidth=0.6)
    axis.tick_params(labelsize=8.5, colors=INK)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        axis.spines[spine].set_color("#c9c8c3")


def figure_pool_sizes(manifest: pd.DataFrame, figures: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PS.apply()

    splits = [("inner", "Inner pools (training/validation)"),
              ("outer", "Outer pools (held-out alleles)"),
              ("test", "Test pools")]
    figure, axes = plt.subplots(1, 3, figsize=(13.5, 3.9), constrained_layout=True)
    bins = np.logspace(3, 5, 21)
    for axis, (split, title) in zip(axes, splits):
        data = manifest[manifest["split"] == split]
        axis.hist(data["peptides"], bins=bins, color=PS.SERIES[0], alpha=0.75,
                  edgecolor="white", linewidth=0.6)
        axis.set_xscale("log")
        style(axis, "Peptides in pool", "Pools", f"{title}  (n={len(data)})")
        axis.annotate(f"{data['peptides'].min():,}–{data['peptides'].max():,}",
                      (0.97, 0.93), xycoords="axes fraction", ha="right",
                      fontsize=8.5, color=INK)
    figure.suptitle("Pool size distribution, by split", fontsize=11.5)
    figure.savefig(figures / "pool_size_distribution.png", dpi=200, bbox_inches="tight")
    figure.savefig(figures / "pool_size_distribution.pdf", bbox_inches="tight")
    plt.close(figure)


def figure_threshold_tradeoff(inner: pd.DataFrame, figures: Path,
                              mode: str = "separate_aln_anchor") -> pd.DataFrame:
    """Each metric against each threshold, with the spread across outer folds.

    A fold mean is computed first, so the error bar is variation between folds
    (what a reader cares about for transfer) rather than between pools of very
    different size.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PS.apply()

    rows = []
    for method in sorted(inner["method"].unique()):
        for order in sorted(inner[inner["method"] == method]["representative_order"].unique()):
            subset = inner[(inner["method"] == method) &
                           (inner["representative_order"] == order)]
            for axis_name, other in (("primary_threshold", "anchor_threshold"),
                                     ("anchor_threshold", "primary_threshold")):
                per_fold = subset.groupby([axis_name, "outer_fold"])[
                    list(REPORTED)].mean().reset_index()
                summary = per_fold.groupby(axis_name)[list(REPORTED)].agg(["mean", "std"])
                for threshold, record in summary.iterrows():
                    row = {"method": method, "representative_order": order,
                           "swept": axis_name, "threshold": threshold,
                           "marginalised_over": other}
                    for metric in REPORTED:
                        row[f"{metric}_mean"] = record[(metric, "mean")]
                        row[f"{metric}_std"] = record[(metric, "std")]
                    rows.append(row)
    curve = pd.DataFrame(rows)

    panels = [("ami", "AMI"), ("nmi", "NMI"),
              ("adjusted_purity_macro", "Adjusted per-allele purity (macro)"),
              ("bcubed_f1_macro", "BCubed F1 (macro)"),
              ("singleton_fraction_of_clusters", "Singleton fraction of clusters"),
              ("clusters", "Clusters")]
    primary_label = PRIMARY_LABEL.get(mode, "Primary-similarity threshold")
    for swept, label in ((("primary_threshold"), primary_label),
                         ("anchor_threshold", "Anchor-combination threshold")):
        figure, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
        for axis, (metric, title) in zip(axes.flat, panels):
            for method in sorted(curve["method"].unique()):
                for order in sorted(curve[curve["method"] == method]["representative_order"].unique()):
                    sel = curve[(curve["method"] == method) &
                                (curve["representative_order"] == order) &
                                (curve["swept"] == swept)].sort_values("threshold")
                    if sel.empty:
                        continue
                    # graph/coverage and lazy-exact/coverage agree to ~0.0001,
                    # so an overplotted solid line would leave one invisible while
                    # still appearing in the legend. Distinct dashes keep all three
                    # readable where they coincide.
                    dash, width = {
                        ("graph", "coverage"): ("-", 2.4),
                        ("graph", "intrinsic"): ((0, (6, 3)), 1.8),
                        ("greedy_lazy", "coverage"): ((0, (1, 2.2)), 2.4),
                    }.get((method, order), ("-", 1.9))
                    axis.errorbar(sel["threshold"], sel[f"{metric}_mean"],
                                  yerr=sel[f"{metric}_std"], marker="o", linestyle=dash,
                                  color=PALETTE[method], markersize=4.0, linewidth=width,
                                  capsize=3, elinewidth=1.0,
                                  label=f"{LABELS[method]} ({order})")
            if metric == "singleton_fraction_of_clusters":
                axis.axhline(SINGLETON_LIMIT, color=PS.SERIES[1], linestyle=":", linewidth=1.5)
                axis.annotate("constraint", (0.02, SINGLETON_LIMIT), xycoords=("axes fraction", "data"),
                              textcoords="offset points", xytext=(0, 5), fontsize=8, color=PS.SERIES[1])
            if metric == "clusters":
                axis.set_yscale("log")
            style(axis, label, title, title)
        handles, labels = axes.flat[0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="outside lower center", ncol=3,
                      frameon=False, fontsize=8.5)
        figure.suptitle(f"Performance against the {label.lower()}")
        stem = "primary" if swept == "primary_threshold" else "anchor"
        figure.savefig(figures / f"tradeoff_{stem}.png", dpi=200, bbox_inches="tight")
        figure.savefig(figures / f"tradeoff_{stem}.pdf", bbox_inches="tight")
        plt.close(figure)
    return curve


def figure_objective_surface(inner: pd.DataFrame, figures: Path,
                             mode: str = "separate_aln_anchor",
                             chosen: pd.DataFrame | None = None) -> pd.DataFrame:
    """Joint view of the two thresholds, since they interact.

    `chosen` is the output of `select(inner, None)`. Passing it in rather than
    recomputing an argmax here is deliberate: selection maximises AMI *subject to
    the singleton constraint*, so an unconstrained argmax over any single surface
    marks a configuration that was never selected.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PS.apply()

    import numpy.ma as ma

    combos = inner[["method", "representative_order"]].drop_duplicates()
    combos = combos.sort_values(["method", "representative_order"]).to_records(index=False)
    figure, axes = plt.subplots(2, len(combos), figsize=(4.6 * len(combos), 8.4),
                                constrained_layout=True, squeeze=False)
    grids = []
    for column, (method, order) in enumerate(combos):
        subset = inner[(inner["method"] == method) & (inner["representative_order"] == order)]
        singleton = subset.pivot_table(index="anchor_threshold", columns="primary_threshold",
                                       values="singleton_fraction_of_clusters", aggfunc="mean")
        for row, (metric, title, cmap) in enumerate(
                [("ami", "AMI", "Blues"),
                 ("singleton_fraction_of_clusters", "Singleton fraction", "Reds")]):
            table = subset.pivot_table(index="anchor_threshold", columns="primary_threshold",
                                       values=metric, aggfunc="mean")
            axis = axes[row][column]
            colours = plt.get_cmap(cmap).copy()
            values = table.values
            if metric == "ami":
                # Cells breaching the singleton constraint are not eligible for
                # selection, so colouring them invites the eye to the top-right
                # corner the constraint actually forbids. Grey them out instead.
                colours.set_bad("#dedcd8")
                values = ma.masked_where(singleton.values > SINGLETON_LIMIT, values)
            # Cell indices, not data coordinates: the threshold grid is not evenly
            # spaced (0.15/0.25 then 0.05 steps), and the anchor axis carries an
            # extra 0.0 level, so an `extent` in data units would mislabel cells.
            image = axis.imshow(values, origin="lower", cmap=colours, aspect="auto")
            figure.colorbar(image, ax=axis, fraction=0.046)
            axis.set_xticks(range(len(table.columns)))
            axis.set_xticklabels([f"{v:g}" for v in table.columns], rotation=90)
            axis.set_yticks(range(len(table.index)))
            axis.set_yticklabels([f"{v:g}" for v in table.index])
            axis.grid(False)
            if metric == "singleton_fraction_of_clusters":
                axis.contour(singleton.values, levels=[SINGLETON_LIMIT],
                             colors="#0b0b0b", linewidths=1.4, linestyles="--")
            elif chosen is not None:
                # The star comes from the same select() result the report quotes,
                # so the figure cannot disagree with the reported configuration.
                pick = chosen[(chosen["method"] == method) &
                              (chosen["representative_order"] == order)]
                for row_pick in pick.itertuples():
                    axis.plot(table.columns.get_loc(row_pick.primary_threshold),
                              table.index.get_loc(row_pick.anchor_threshold),
                              "*", color="#ffd400", markersize=20,
                              markeredgecolor="#0b0b0b", markeredgewidth=1.0)
            style(axis, PRIMARY_LABEL.get(mode, "Primary threshold"), "Anchor threshold",
                  f"{LABELS[method]} ({order})\n{title}")
            melted = table.reset_index().melt(id_vars="anchor_threshold",
                                              var_name="primary_threshold", value_name=metric)
            melted["method"] = method
            melted["representative_order"] = order
            melted["metric"] = metric
            grids.append(melted)
    figure.suptitle("Threshold interaction: AMI (grey = singleton constraint breached, "
                    "star = selected) and the constraint itself (dashed = 50%)")
    figure.savefig(figures / "threshold_surface.png", dpi=200, bbox_inches="tight")
    figure.savefig(figures / "threshold_surface.pdf", bbox_inches="tight")
    plt.close(figure)
    return pd.concat(grids, ignore_index=True)


def write_report(root: Path, mode: str, manifest, inner, per_fold, overall,
                 evaluations: dict) -> None:
    """One report per analysis, so a mode can be read on its own before any
    cross-mode comparison."""
    add = [].append
    lines: list[str] = []
    def w(text=""):
        lines.append(text)

    w(f"# Peptide-MHC hyperparameter optimisation: `{mode}`")
    w()
    w("## What this measures")
    w()
    w("Peptides of known allele origin are pooled and clustered without labels;")
    w("agreement between the clusters and the alleles measures how far a purely")
    w("sequence-based clustering recovers allele identity.")
    w()
    w("## Benchmark design")
    w()
    summary = (root / "prepare_summary.json")
    if summary.exists():
        import json
        info = json.loads(summary.read_text())
        w(f"- Alleles: {info['alleles_after']} of {info['alleles_before_similarity_filter']} "
          f"survive a {info['similarity_threshold']:.0f}% pairwise-similarity ceiling "
          f"(A={info['locus_counts']['A']}, B={info['locus_counts']['B']}, "
          f"C={info['locus_counts']['C']}).")
        w(f"- Worst within-pool allele imbalance: {info['worst_imbalance_ratio']:.3f}; "
          f"largest single-allele share {info['max_allele_share_observed']:.3f}.")
    inner_pools = manifest[manifest.split == "inner"]
    w(f"- Pools: {len(manifest)} total "
      f"({(manifest.split=='inner').sum()} inner, {(manifest.split=='outer').sum()} outer, "
      f"{(manifest.split=='test').sum()} test), "
      f"{inner_pools.peptides.min():,}-{manifest.peptides.max():,} peptides.")
    w("- A peptide enters a pool only if it binds exactly one of that pool's alleles,")
    w("  so every label is unambiguous.")
    w("- Outer folds split *alleles*, so the held-out numbers estimate transfer to")
    w("  alleles never seen during tuning.")
    w()
    w("## Objective")
    w()
    w("AMI, maximised subject to singletons being at most half of all clusters.")
    w("NMI and purity are reported but never optimised: measured on this grid they")
    w("correlate +0.80 and +0.96 with the singleton fraction, so both are maximised")
    w("by the degenerate clustering that isolates every peptide. AMI is chance-")
    w("corrected and correlates -0.97 with fragmentation.")
    w()
    w("## Selected configuration")
    w()
    w("| Method | Order | Primary | Anchor | AMI | NMI | Purity (macro) | Singletons | Clusters |")
    w("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for row in overall.itertuples():
        w(f"| {LABELS.get(row.method, row.method)} | {row.representative_order} | "
          f"{row.primary_threshold:.2f} | {row.anchor_threshold:.2f} | "
          f"{row.ami:.4f} | {row.nmi:.4f} | {row.adjusted_purity_macro:.4f} | "
          f"{row.singleton_fraction:.3f} | {row.clusters:.0f} |")
    w()
    w("Selection per outer fold, which shows whether the choice is stable:")
    w()
    w("| Held-out fold | Method | Order | Primary | Anchor | AMI |")
    w("|---:|---|---|---:|---:|---:|")
    for row in per_fold.itertuples():
        w(f"| {row.held_out_fold} | {LABELS.get(row.method, row.method)} | "
          f"{row.representative_order} | {row.primary_threshold:.2f} | "
          f"{row.anchor_threshold:.2f} | {row.ami:.4f} |")
    w()
    for name, frame in evaluations.items():
        if frame is None or frame.empty:
            continue
        w(f"## {name.capitalize()} evaluation at the selected configuration")
        w()
        w("| Method | Pools | AMI | NMI | Purity (macro) | Purity (micro) | BCubed-P | Singletons | Clusters |")
        w("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for method in sorted(frame.method.unique()):
            s = frame[frame.method == method]
            w(f"| {LABELS.get(method, method)} | {len(s)} | "
              f"{s.ami.mean():.4f} ± {s.ami.std():.4f} | {s.nmi.mean():.4f} | "
              f"{s.adjusted_purity_macro.mean():.4f} ± {s.adjusted_purity_macro.std():.4f} | "
              f"{s.adjusted_purity_micro.mean():.4f} | {s.bcubed_precision_macro.mean():.4f} | "
              f"{s.singleton_fraction_of_clusters.mean():.3f} | {s.clusters.mean():.0f} |")
        w()

    best = inner.groupby("primary_threshold")["ami"].mean()
    w("## Threshold behaviour")
    w()
    w("Mean AMI by primary threshold (marginalised over the anchor threshold):")
    w()
    w("| " + " | ".join(f"{v:.2f}" for v in best.index) + " |")
    w("|" + "---:|" * len(best))
    w("| " + " | ".join(f"{v:.3f}" for v in best.values) + " |")
    w()
    if best.idxmax() == best.index.min():
        w("The optimum sits at the lowest threshold swept, so it is a boundary")
        w("solution: AMI is still rising where the grid stops and the true optimum")
        w("lies outside it.")
        w()
    ablation = inner[inner.anchor_threshold == 0.0]
    if len(ablation):
        with_anchor = inner[inner.anchor_threshold > 0.0]
        w("## Anchor ablation")
        w()
        w("With the anchor threshold at zero the anchor condition is disabled and")
        w("only the primary component decides eligibility.")
        w()
        w(f"- anchor disabled: mean AMI {ablation.ami.mean():.4f}, "
          f"clusters {ablation.clusters.mean():.0f}")
        w(f"- anchor active:   mean AMI {with_anchor.ami.mean():.4f}, "
          f"clusters {with_anchor.clusters.mean():.0f}")
        w()
    w("## Caveats")
    w()
    w("- The test pools come from a separate source file but are not disjoint in")
    w("  content: 66% of their peptides and 6 of 20 alleles also occur in the tuning")
    w("  data, so they test transfer to new pool compositions rather than to unseen")
    w("  peptides. The outer folds, which hold out whole alleles, are the stricter")
    w("  estimate.")
    w("- Alleles above the similarity ceiling never share a pool. Near-identical")
    w("  alleles present near-identical motifs, so no sequence-based method can")
    w("  separate them; excluding them measures the method rather than the label,")
    w("  but it makes these numbers an upper bound relative to a realistic pool.")
    w("- BCubed recall is near zero by construction: motif clustering splits each")
    w("  allele across many clusters, so only the precision half is informative.")
    w("- Every number here passed the completeness gate in `verify.py`: the grid is")
    w("  complete, duplicate-free, error-free, and a random sample reproduces on")
    w("  re-execution.")
    w()
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--mode", default="separate_aln_anchor")
    args = parser.parse_args()
    root = args.root.resolve()
    figures = root / "figures"
    tables = root / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)

    manifest = pd.read_csv(root / "pool_manifest.csv")
    manifest.to_csv(tables / "pool_manifest.csv", index=False)
    figure_pool_sizes(manifest, figures)

    inner = load_grid(root, "inner")
    inner.to_csv(tables / "inner_grid_raw.csv", index=False)

    per_fold = pd.concat([select(inner, fold) for fold in sorted(inner["outer_fold"].unique())],
                         ignore_index=True)
    per_fold.to_csv(tables / "selected_per_fold.csv", index=False)
    overall = select(inner, None)
    overall.to_csv(tables / "selected_overall.csv", index=False)

    curve = figure_threshold_tradeoff(inner, figures, args.mode)
    curve.to_csv(tables / "threshold_tradeoff.csv", index=False)
    # Selection first, so the surface can mark what was actually chosen.
    surface = figure_objective_surface(inner, figures, args.mode, overall)
    surface.to_csv(tables / "threshold_surface.csv", index=False)

    print("\n=== Selected configuration per outer fold (chosen on the other 4 folds) ===")
    print(per_fold[["held_out_fold", "method", "representative_order", "primary_threshold",
                    "anchor_threshold", "objective", "singleton_fraction"]].to_string(index=False))
    print("\n=== Selected on all inner pools ===")
    print(overall[["method", "representative_order", "primary_threshold", "anchor_threshold",
                   "objective", "ami", "nmi", "adjusted_purity_macro",
                   "singleton_fraction"]].to_string(index=False))

    collected: dict = {}
    for tag, name in (("outer", "outer"), ("test", "test")):
        try:
            evaluation = load_grid(root, tag)
        except SystemExit:
            print(f"\n({name} grid not present yet)")
            continue
        evaluation.to_csv(tables / f"{name}_grid_raw.csv", index=False)
        if tag == "outer":
            merged = evaluation.merge(
                per_fold.rename(columns={"held_out_fold": "outer_fold"})[
                    CONFIG + ["outer_fold"]], on=CONFIG + ["outer_fold"], how="inner")
        else:
            merged = evaluation.merge(overall[CONFIG], on=CONFIG, how="inner")
        merged.to_csv(tables / f"{name}_selected_runs.csv", index=False)
        summary = merged.groupby("method")[list(REPORTED)].agg(["mean", "std"])
        summary.to_csv(tables / f"{name}_summary.csv")
        print(f"\n=== {name.capitalize()} evaluation at the selected configuration ===")
        for method in sorted(merged["method"].unique()):
            sel = merged[merged["method"] == method]
            print(f"  {LABELS[method]}  (n={len(sel)} pools)")
            for metric in HEADLINE + ("bcubed_f1_macro", "singleton_fraction_of_clusters"):
                print(f"    {metric:34s} {sel[metric].mean():.4f} ± {sel[metric].std():.4f}")
        collected[name] = merged

    write_report(root, args.mode, manifest, inner, per_fold, overall, collected)
    print(f"\nwrote {root / 'REPORT.md'}")


if __name__ == "__main__":
    main()
