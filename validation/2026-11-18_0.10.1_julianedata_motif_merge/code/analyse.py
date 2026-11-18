"""Summarise the independent Juliane-data benchmark.

No tuning happens here: every arm is the frozen 0.11.0 setting from the
2026-08-06 study, and this script only reduces the 1,440 runs to tables.
"""
import argparse
import pandas as pd
from scipy.stats import wilcoxon

MODES = [
    ("similarity", "PepCluster2"),
    ("merge",      "PepCluster2 + merge"),
    ("merge_k",    "PepCluster2 + merge (given K)"),
    ("merge_em",   "PepCluster2 + merge + EM"),
    ("merge_em_k", "PepCluster2 + merge + EM (given K)"),
    ("em_only",    "PepCluster2 + EM"),
    ("em_only_k",  "PepCluster2 + EM (given K)"),
]
METRICS = ["ami", "bcubed_precision_macro", "bcubed_recall_macro",
           "bcubed_f1_macro", "adjusted_purity_macro",
           "singleton_fraction_of_peptides"]


def load(results):
    frames = []
    for key, label in MODES:
        d = pd.read_csv(f"{results}/variant_{key}.csv")
        d["arm"] = label
        frames.append(d)
    m = pd.read_csv(f"{results}/mixmhcp.csv")
    m["arm"] = m["setting"].map({"default": "MixMHCp (default)",
                                 "oracle_k": "MixMHCp (given K)"})
    frames.append(m)
    cols = ["arm", "pool", "allele_count", "peptides", "status", "clusters",
            "elapsed_seconds"] + METRICS
    return pd.concat([f[cols] for f in frames], ignore_index=True)


ORDER = [l for _, l in MODES] + ["MixMHCp (default)", "MixMHCp (given K)"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--manifest", default="pool_manifest.csv")
    args = ap.parse_args()

    df = load(args.results)
    man = pd.read_csv(args.manifest)[["pool", "balance", "imbalance_ratio"]]
    df = df.merge(man, on="pool", how="left")
    bad = df[df.status != "ok"]
    print(f"{len(df)} runs, {len(bad)} not ok")
    if len(bad):
        print(bad[["arm", "pool", "status"]].to_string(index=False))

    print("\n== headline (mean over all 160 pools) ==")
    h = (df.groupby("arm")[METRICS + ["clusters", "elapsed_seconds"]]
           .mean().reindex(ORDER).round(3))
    print(h.to_string())

    print("\n== AMI by allele count ==")
    p = df.pivot_table(index="arm", columns="allele_count",
                       values="ami", aggfunc="mean").reindex(ORDER).round(3)
    print(p.to_string())

    print("\n== BCubed F1 (macro) by allele count ==")
    p = df.pivot_table(index="arm", columns="allele_count",
                       values="bcubed_f1_macro", aggfunc="mean").reindex(ORDER).round(3)
    print(p.to_string())

    print("\n== AMI by imbalance regime ==")
    p = df.pivot_table(index="arm", columns="balance",
                       values="ami", aggfunc="mean").reindex(ORDER).round(3)
    print(p[["eq", "a10", "a3", "a1", "a03"]].to_string())

    print("\n== clusters recovered vs true alleles ==")
    p = df.pivot_table(index="arm", columns="allele_count",
                       values="clusters", aggfunc="mean").reindex(ORDER).round(1)
    print(p.to_string())

    print("\n== median seconds by allele count ==")
    p = df.pivot_table(index="arm", columns="allele_count",
                       values="elapsed_seconds", aggfunc="median").reindex(ORDER).round(2)
    print(p.to_string())

    print("\n== paired Wilcoxon on AMI, 160 pools ==")
    wide = df.pivot_table(index="pool", columns="arm", values="ami")
    pairs = [("PepCluster2 + merge + EM", "MixMHCp (default)"),
             ("PepCluster2 + merge + EM (given K)", "MixMHCp (given K)"),
             ("PepCluster2 + merge", "PepCluster2"),
             ("PepCluster2 + merge + EM", "PepCluster2 + merge"),
             ("PepCluster2 + merge + EM", "PepCluster2 + EM"),
             ("PepCluster2 + merge + EM (given K)", "PepCluster2 + merge + EM")]
    for a, b in pairs:
        w = wide[[a, b]].dropna()
        st, pv = wilcoxon(w[a], w[b])
        d = (w[a] - w[b]).mean()
        print(f"  {a:38s} vs {b:34s} delta={d:+.4f} p={pv:.2e} n={len(w)}")

    df.to_csv(f"{args.results}/all_arms_long.csv", index=False)
    print(f"\nwrote {args.results}/all_arms_long.csv")


main()
