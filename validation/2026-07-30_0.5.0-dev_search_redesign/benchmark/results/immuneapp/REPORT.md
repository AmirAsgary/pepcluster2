# External tool comparison: immuneapp

Identical pools, identical metrics. PepCluster2 is evaluated at the
configuration selected by its own nested cross-validation; the external
tools have no threshold to tune, so they are run as documented.

## How to read this

PepCluster2 is a similarity clustering method and returns on the order of a
hundred clusters per pool. MixMHCp and GibbsCluster are mixture models that
fit a handful of motifs, so partition size differs by more than an order of
magnitude and the metrics must be read as a set.

Adjusted per-allele purity is BCubed precision corrected against the allele
prior. That correction removes the baseline for one large cluster but not the
inflation from fragmentation: a singleton scores precision 1, so a partition
into singletons scores a perfect 1.0. **Purity alone therefore cannot be used
to compare partitions of different granularity.** BCubed recall penalises
exactly what purity rewards, and BCubed F1 balances them. AMI is
chance-corrected over the whole partition. F1 and AMI are the two figures
that can be compared across tools directly.

`forced k` gives a tool the true number of alleles in the pool. No user could
do that in practice, so it is not a fair headline number - it isolates how
much of a tool's result is its model rather than its model selection.

## Tuning folds

| Tool | Pools | AMI | Purity (macro) | Recall | F1 | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 120 | 0.2753 ± 0.0505 | 0.4884 ± 0.2120 | 0.0521 | 0.0918 | 0.073 | 195.4 |
| PepCluster2 alignment | 120 | 0.2520 ± 0.0497 | 0.4263 ± 0.2092 | 0.0500 | 0.0873 | 0.032 | 154.3 |
| MixMHCp (default) | 120 | 0.5012 ± 0.2372 | 0.3853 ± 0.3451 | 0.7947 | 0.4949 | 0.003 | 4.5 |
| MixMHCp (forced k) | 120 | 0.5474 ± 0.2272 | 0.5013 ± 0.3099 | 0.5524 | 0.5449 | 0.000 | 13.7 |
| GibbsCluster (default) | 120 | 0.2780 ± 0.1978 | 0.2442 ± 0.2908 | 0.5524 | 0.3631 | 0.000 | 3.6 |
| GibbsCluster (forced k) | 120 | 0.3127 ± 0.1716 | 0.2971 ± 0.2713 | 0.3406 | 0.3570 | 0.000 | 13.7 |

## Held-out alleles

| Tool | Pools | AMI | Purity (macro) | Recall | F1 | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 25 | 0.3071 ± 0.0369 | 0.5649 ± 0.1913 | 0.0564 | 0.1024 | 0.074 | 163.0 |
| PepCluster2 alignment | 25 | 0.2835 ± 0.0368 | 0.5074 ± 0.2074 | 0.0520 | 0.0943 | 0.039 | 131.8 |
| MixMHCp (default) | 25 | 0.5965 ± 0.1542 | 0.4431 ± 0.2668 | 0.8205 | 0.5955 | 0.000 | 4.9 |
| MixMHCp (forced k) | 25 | 0.6876 ± 0.1092 | 0.6689 ± 0.1631 | 0.7256 | 0.7089 | 0.000 | 8.7 |
| GibbsCluster (default) | 25 | 0.3604 ± 0.1908 | 0.2918 ± 0.3036 | 0.5916 | 0.4352 | 0.000 | 3.8 |
| GibbsCluster (forced k) | 25 | 0.3738 ± 0.1736 | 0.3517 ± 0.2710 | 0.3996 | 0.4161 | 0.000 | 8.7 |

## Benchmark

| Tool | Pools | AMI | Purity (macro) | Recall | F1 | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 48 | 0.3333 ± 0.0208 | 0.4974 ± 0.1212 | 0.0624 | 0.1098 | 0.078 | 174.9 |
| PepCluster2 alignment | 48 | 0.3185 ± 0.0204 | 0.4457 ± 0.1203 | 0.0605 | 0.1057 | 0.031 | 139.2 |
| MixMHCp (default) | 48 | 0.3918 ± 0.2289 | 0.2215 ± 0.2252 | 0.8282 | 0.3918 | 0.000 | 4.1 |
| MixMHCp (forced k) | 48 | 0.4915 ± 0.2033 | 0.4178 ± 0.2576 | 0.4843 | 0.4728 | 0.000 | 12.4 |
| GibbsCluster (default) | 48 | 0.1799 ± 0.1393 | 0.1000 ± 0.1330 | 0.6046 | 0.2656 | 0.000 | 2.9 |
| GibbsCluster (forced k) | 48 | 0.2522 ± 0.0905 | 0.1861 ± 0.1226 | 0.2576 | 0.2628 | 0.000 | 12.4 |

## Where the difference comes from

The split averages hide a strong interaction with pool complexity.
Paired per-pool AMI, positive meaning the external tool is ahead:

| Alleles in pool | Tool | vs | Pools | Mean AMI gap | Tool wins |
|---|---|---|---:|---:|---:|
| 2-6 | MixMHCp (default) | PepCluster2 k-mer | 72 | +0.4112 | 97% |
| 2-6 | MixMHCp (default) | PepCluster2 alignment | 72 | +0.4304 | 97% |
| 7-12 | MixMHCp (default) | PepCluster2 k-mer | 55 | +0.1553 | 82% |
| 7-12 | MixMHCp (default) | PepCluster2 alignment | 55 | +0.1799 | 84% |
| 13-30 | MixMHCp (default) | PepCluster2 k-mer | 66 | -0.0151 | 53% |
| 13-30 | MixMHCp (default) | PepCluster2 alignment | 66 | +0.0056 | 62% |
| 2-6 | MixMHCp (forced k) | PepCluster2 k-mer | 72 | +0.4682 | 100% |
| 2-6 | MixMHCp (forced k) | PepCluster2 alignment | 72 | +0.4874 | 100% |
| 7-12 | MixMHCp (forced k) | PepCluster2 k-mer | 55 | +0.2416 | 89% |
| 7-12 | MixMHCp (forced k) | PepCluster2 alignment | 55 | +0.2662 | 91% |
| 13-30 | MixMHCp (forced k) | PepCluster2 k-mer | 66 | +0.0417 | 53% |
| 13-30 | MixMHCp (forced k) | PepCluster2 alignment | 66 | +0.0624 | 55% |
| 2-6 | GibbsCluster (default) | PepCluster2 k-mer | 72 | +0.1825 | 82% |
| 2-6 | GibbsCluster (default) | PepCluster2 alignment | 72 | +0.2017 | 83% |
| 7-12 | GibbsCluster (default) | PepCluster2 k-mer | 55 | -0.1053 | 7% |
| 7-12 | GibbsCluster (default) | PepCluster2 alignment | 55 | -0.0807 | 11% |
| 13-30 | GibbsCluster (default) | PepCluster2 k-mer | 66 | -0.1979 | 0% |
| 13-30 | GibbsCluster (default) | PepCluster2 alignment | 66 | -0.1772 | 0% |
| 2-6 | GibbsCluster (forced k) | PepCluster2 k-mer | 72 | +0.1959 | 92% |
| 2-6 | GibbsCluster (forced k) | PepCluster2 alignment | 72 | +0.2151 | 94% |
| 7-12 | GibbsCluster (forced k) | PepCluster2 k-mer | 55 | -0.0705 | 4% |
| 7-12 | GibbsCluster (forced k) | PepCluster2 alignment | 55 | -0.0459 | 16% |
| 13-30 | GibbsCluster (forced k) | PepCluster2 k-mer | 66 | -0.1208 | 0% |
| 13-30 | GibbsCluster (forced k) | PepCluster2 alignment | 66 | -0.1001 | 0% |

MixMHCp is a mixture model built to resolve a few motifs, and that is
exactly where it dominates. Its margin shrinks steadily as alleles are
added and is gone by roughly 13 alleles, where the two are level.

Mean adjusted per-allele purity over the same bands:

| Alleles in pool | PepCluster2 k-mer | PepCluster2 alignment | MixMHCp (default) | MixMHCp (forced k) | GibbsCluster (default) | GibbsCluster (forced k) |
|---|---|---|---|---|---|---|
| 2-6 | 0.6982 | 0.6332 | 0.6794 | 0.7958 | 0.4840 | 0.5381 |
| 7-12 | 0.4691 | 0.4080 | 0.2625 | 0.4918 | 0.0949 | 0.1751 |
| 13-30 | 0.3111 | 0.2605 | 0.0697 | 0.1906 | 0.0201 | 0.0759 |

PepCluster2 leads on purity in the 13-30 band(s). Read that
against the cluster counts in the tables above and against F1 below
before drawing a conclusion: purity rises mechanically with the number
of clusters, so a lead on purity held by the partition with an order of
magnitude more clusters is not evidence of a better partition.

Mean BCubed F1 over the same bands, which penalises fragmentation and
lumping together and is the figure to compare across tools:

| Alleles in pool | PepCluster2 k-mer | PepCluster2 alignment | MixMHCp (default) | MixMHCp (forced k) | GibbsCluster (default) | GibbsCluster (forced k) |
|---|---|---|---|---|---|---|
| 2-6 | 0.0955 | 0.0913 | 0.7769 | 0.8339 | 0.6185 | 0.6096 |
| 7-12 | 0.1006 | 0.0951 | 0.4465 | 0.5471 | 0.2707 | 0.2523 |
| 13-30 | 0.0976 | 0.0926 | 0.1908 | 0.2375 | 0.1178 | 0.1227 |

## Files

- `tables/per_pool.csv` - every tool on every pool
- `tables/summary.csv` - the tables above
- `tables/by_complexity.csv` - the paired comparison by allele count
- `tool_comparison.png` / `.pdf`, `by_complexity.png` / `.pdf`
- `raw/` - each tool's own run output

