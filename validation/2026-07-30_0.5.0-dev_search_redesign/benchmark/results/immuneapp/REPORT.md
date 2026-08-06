# External tool comparison: immuneapp

Identical pools, identical metrics. PepCluster2 is evaluated at the
configuration selected by its own nested cross-validation; the external
tools have no threshold to tune, so they are run as documented.

> Not yet included: GibbsCluster.

## How to read this

PepCluster2 is a similarity clustering method and returns on the order of a
hundred clusters per pool. MixMHCp and GibbsCluster are mixture models that
fit a handful of motifs. The two metrics therefore disagree by construction:
AMI is chance-corrected and penalises the finer partition, while purity
rewards it. Both are reported; neither on its own decides the question.

`forced k` gives a tool the true number of alleles in the pool. No user could
do that in practice, so it is not a fair headline number - it isolates how
much of a tool's result is its model rather than its model selection.

## Tuning folds

| Tool | Pools | AMI | Purity (macro) | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 120 | 0.2753 ± 0.0505 | 0.4884 ± 0.2120 | 0.073 | 195.4 |
| PepCluster2 alignment | 120 | 0.2520 ± 0.0497 | 0.4263 ± 0.2092 | 0.032 | 154.3 |
| MixMHCp (default) | 120 | 0.5012 ± 0.2372 | 0.3853 ± 0.3451 | 0.003 | 4.5 |
| MixMHCp (forced k) | 120 | 0.5474 ± 0.2272 | 0.5013 ± 0.3099 | 0.000 | 13.7 |

## Held-out alleles

| Tool | Pools | AMI | Purity (macro) | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 25 | 0.3071 ± 0.0369 | 0.5649 ± 0.1913 | 0.074 | 163.0 |
| PepCluster2 alignment | 25 | 0.2835 ± 0.0368 | 0.5074 ± 0.2074 | 0.039 | 131.8 |
| MixMHCp (default) | 25 | 0.5965 ± 0.1542 | 0.4431 ± 0.2668 | 0.000 | 4.9 |
| MixMHCp (forced k) | 25 | 0.6876 ± 0.1092 | 0.6689 ± 0.1631 | 0.000 | 8.7 |

## Benchmark

| Tool | Pools | AMI | Purity (macro) | Singletons | Clusters |
|---|---:|---:|---:|---:|---:|
| PepCluster2 k-mer | 48 | 0.3333 ± 0.0208 | 0.4974 ± 0.1212 | 0.078 | 174.9 |
| PepCluster2 alignment | 48 | 0.3185 ± 0.0204 | 0.4457 ± 0.1203 | 0.031 | 139.2 |
| MixMHCp (default) | 48 | 0.3918 ± 0.2289 | 0.2215 ± 0.2252 | 0.000 | 4.1 |
| MixMHCp (forced k) | 48 | 0.4915 ± 0.2033 | 0.4178 ± 0.2576 | 0.000 | 12.4 |

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

MixMHCp is a mixture model built to resolve a few motifs, and that is
exactly where it dominates. Its margin shrinks steadily as alleles are
added and is gone by roughly 13 alleles, where the two are level.

Mean adjusted per-allele purity over the same bands:

| Alleles in pool | PepCluster2 k-mer | PepCluster2 alignment | MixMHCp (default) | MixMHCp (forced k) |
|---|---|---|---|---|
| 2-6 | 0.6982 | 0.6332 | 0.6794 | 0.7958 |
| 7-12 | 0.4691 | 0.4080 | 0.2625 | 0.4918 |
| 13-30 | 0.3111 | 0.2605 | 0.0697 | 0.1906 |

PepCluster2 leads on purity in the 13-30 band(s); it does
not lead everywhere. So neither tool is simply better: MixMHCp wins
decisively on simple pools, PepCluster2 holds up better as pools grow,
and which matters depends on how many alleles a sample actually mixes.

## Files

- `tables/per_pool.csv` - every tool on every pool
- `tables/summary.csv` - the tables above
- `tables/by_complexity.csv` - the paired comparison by allele count
- `tool_comparison.png` / `.pdf`, `by_complexity.png` / `.pdf`
- `raw/` - each tool's own run output

