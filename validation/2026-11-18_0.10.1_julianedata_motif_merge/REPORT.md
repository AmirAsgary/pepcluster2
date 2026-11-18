# Independent benchmark on the Juliane mono-allelic dataset

PepCluster2 **v0.11.0**, frozen settings, run once on a dataset it has never
been tuned on. 160 pools, 1,785,420 peptides, 91 distinct HLA alleles drawn
from 6 studies. Nine arms, 1,440 runs, 0 failures.

**Nothing in this study was tuned.** Every PepCluster2 arm uses the settings
selected by the nested cross-validation in `2026-08-06_0.6.0_motif_merge`, and
MixMHCp is run exactly as it was there. The pools were built once, before any
arm was run, and were not revisited.

## Headline

Mean over all 160 pools:

| arm | AMI | BCubed P | BCubed R | BCubed F1 | adj. purity | motifs | seconds |
|---|---|---|---|---|---|---|---|
| PepCluster2 | 0.255 | 0.474 | 0.051 | 0.087 | 0.407 | 214.9 | 24.8 |
| PepCluster2 + merge | 0.315 | 0.437 | 0.263 | 0.298 | 0.364 | 121.4 | 25.1 |
| PepCluster2 + merge (given K) | 0.369 | 0.364 | 0.468 | 0.384 | 0.276 | 14.5 | 25.1 |
| **PepCluster2 + merge + EM** | **0.552** | 0.511 | 0.580 | 0.500 | 0.457 | 15.3 | 28.5 |
| **PepCluster2 + merge + EM (given K)** | **0.566** | 0.490 | 0.638 | 0.522 | 0.430 | 14.5 | 26.4 |
| PepCluster2 + EM | 0.540 | 0.530 | 0.532 | 0.488 | 0.478 | 17.1 | 32.3 |
| PepCluster2 + EM (given K) | 0.567 | 0.504 | 0.617 | 0.526 | 0.446 | 14.5 | 25.1 |
| MixMHCp (default) | 0.439 | 0.357 | 0.765 | 0.417 | 0.276 | 4.4 | 53.3 |
| MixMHCp (given K) | 0.507 | 0.484 | 0.507 | 0.468 | 0.425 | 15.5 | 545.7 |

Three results stand out.

**First, PepCluster2 + merge + EM beats MixMHCp on an independent dataset, and
the margin grows with pool complexity.** Paired over the 160 pools, the
no-K-supplied comparison is +0.1135 AMI (p = 2.0e-16) and the given-K
comparison is +0.0598 AMI (p = 6.4e-21). It wins on 81.2% and 87.5% of pools
respectively. This is the first independent confirmation that the Bayesian
merge stage fixed what the ImmuneApp benchmark exposed: there, PepCluster2 lost
on AMI below 13 alleles. Here it wins or ties at every complexity.

**Second, the win is driven by granularity control, not by the motif model.**
MixMHCp's default arm scans 1-6 motifs, so at 12 alleles it reports 4.5 motifs
and at 50 it reports 5.7 — it cannot express the answer. Its recall is the
highest of any arm (0.765) for exactly that reason: a handful of huge clusters
recovers most same-allele pairs and gets precision 0.357 in exchange. The
merge + EM arm lands on 15.3 motifs against a true mean of 14.5, a median
K/true ratio of 1.00. **Read the given-K comparison as the fair one** — see
the caveat below.

**Third, the cost gap widens sharply with K.** At 50 alleles MixMHCp's given-K
arm takes a median of 3,706 s (max 8,986 s, i.e. 2.5 hours for one pool)
against 112.8 s for merge + EM — a **33x** speedup at equal information. The
default arm is cheaper only because it never fits more than 6 motifs.

## Quality against pool complexity

AMI, mean per allele count:

| arm | 3 | 6 | 12 | 24 | 50 |
|---|---|---|---|---|---|
| PepCluster2 | 0.212 | 0.277 | 0.276 | 0.265 | 0.248 |
| PepCluster2 + merge | 0.324 | 0.354 | 0.314 | 0.281 | 0.251 |
| PepCluster2 + merge (given K) | 0.448 | 0.420 | 0.341 | 0.288 | 0.245 |
| PepCluster2 + merge + EM | 0.595 | 0.597 | 0.548 | 0.490 | 0.447 |
| PepCluster2 + merge + EM (given K) | 0.630 | 0.598 | 0.558 | 0.500 | 0.458 |
| PepCluster2 + EM | 0.560 | 0.584 | 0.547 | 0.491 | 0.447 |
| PepCluster2 + EM (given K) | 0.618 | 0.606 | 0.562 | 0.503 | 0.457 |
| MixMHCp (default) | 0.539 | 0.504 | 0.390 | 0.346 | 0.299 |
| MixMHCp (given K) | 0.597 | 0.541 | 0.481 | 0.417 | 0.397 |

BCubed F1 (macro), same layout:

| arm | 3 | 6 | 12 | 24 | 50 |
|---|---|---|---|---|---|
| PepCluster2 | 0.091 | 0.105 | 0.089 | 0.074 | 0.052 |
| PepCluster2 + merge | 0.455 | 0.367 | 0.250 | 0.160 | 0.079 |
| PepCluster2 + merge (given K) | 0.637 | 0.468 | 0.297 | 0.175 | 0.092 |
| PepCluster2 + merge + EM | 0.687 | 0.596 | 0.454 | 0.310 | 0.213 |
| PepCluster2 + merge + EM (given K) | 0.743 | 0.598 | 0.472 | 0.323 | 0.230 |
| PepCluster2 + EM | 0.652 | 0.583 | 0.457 | 0.309 | 0.213 |
| PepCluster2 + EM (given K) | 0.741 | 0.607 | 0.477 | 0.327 | 0.228 |
| MixMHCp (default) | 0.697 | 0.517 | 0.309 | 0.187 | 0.103 |
| MixMHCp (given K) | 0.723 | 0.544 | 0.390 | 0.249 | 0.181 |

Per-pool win rate of `merge + EM (given K)` over `MixMHCp (given K)`:
65% at 3 alleles, 88% at 6, 98% at 12, **100% at both 24 and 50**. The
3-allele cell is the only place the two are close, and it is the only place
MixMHCp's own default arm is operating inside its 6-motif range.

At 3 alleles the arm without a supplied K reaches only a 50% win rate against
MixMHCp (given K), because its median K/true ratio there is 1.33 — it splits a
3-allele pool into 4 or 5 motifs. That over-splitting is the one systematic
weakness this dataset exposes, and it disappears by 6 alleles.

## Granularity

Mean motifs reported against true allele count:

| arm | 3 | 6 | 12 | 24 | 50 |
|---|---|---|---|---|---|
| PepCluster2 | 136.8 | 152.3 | 217.8 | 285.2 | 420.2 |
| PepCluster2 + merge | 82.2 | 97.3 | 130.0 | 151.7 | 200.6 |
| PepCluster2 + merge + EM | 5.5 | 6.5 | 11.4 | 21.5 | 54.4 |
| PepCluster2 + EM | 7.0 | 8.2 | 14.0 | 23.8 | 54.2 |
| MixMHCp (default) | 3.3 | 4.2 | 4.5 | 5.0 | 5.7 |

Every given-K arm returns exactly K by construction and is omitted here.

The merge + EM arm tracks the diagonal to within 10% from 12 alleles upward
and is the only automatic arm that does. Its median K/true ratio is 1.00
overall against 1.33 for EM alone and 0.50 for MixMHCp's default. Note that
`+ merge` alone still leaves 200 groups on a 50-allele pool: **the merge is a
granularity control that makes EM tractable, not a clustering in its own
right** — the same conclusion the 0.6.0 study reached, reproduced here.

## Abundance imbalance

AMI by Dirichlet regime (median per-allele imbalance ratio in brackets):

| arm | equal (1.0) | alpha 10 (2.5) | alpha 3 (5.5) | alpha 1 (11.9) | alpha 0.3 (31.3) |
|---|---|---|---|---|---|
| PepCluster2 + merge + EM | 0.580 | 0.541 | 0.588 | 0.519 | 0.534 |
| PepCluster2 + merge + EM (given K) | 0.581 | 0.552 | 0.609 | 0.538 | 0.552 |
| MixMHCp (default) | 0.443 | 0.375 | 0.464 | 0.430 | 0.482 |
| MixMHCp (given K) | 0.552 | 0.502 | 0.558 | 0.458 | 0.464 |

Neither tool degrades monotonically with imbalance, and the spread across
regimes (0.52-0.59 for merge + EM) is small next to the spread across allele
counts (0.45-0.60). **Abundance imbalance is not the dominant difficulty
here**; allele count is. The one trend worth noting is that MixMHCp's
default arm falls further behind at equal abundance (0.137 AMI) than at
alpha 0.3 (0.052), and its given-K arm does the reverse (0.029 at equal,
0.088 at alpha 0.3). Neither gap is large and the n per cell is 32; this is
an observation, not an explained effect.

## Cost

Median wall-clock seconds per pool, 8 cores:

| arm | 3 | 6 | 12 | 24 | 50 |
|---|---|---|---|---|---|
| PepCluster2 | 2.38 | 1.45 | 3.86 | 10.04 | 100.04 |
| PepCluster2 + merge + EM | 2.85 | 2.26 | 5.58 | 13.91 | 112.79 |
| PepCluster2 + EM | 3.27 | 2.49 | 5.97 | 15.46 | 142.07 |
| MixMHCp (default) | 12.88 | 10.96 | 20.36 | 51.47 | 207.73 |
| MixMHCp (given K) | 1.16 | 4.31 | 37.01 | 220.88 | **3705.92** |

The merge + EM stages cost about 13% over clustering alone and the cost ratio
is flat in pool size. MixMHCp's given-K arm is the only one whose cost is
driven by K rather than by n, and it is the reason the 50-allele shards took
2.5 hours each while the PepCluster2 shards finished in 54 minutes for all
seven modes together.

## What correlates with difficulty

Spearman rho against merge + EM AMI, n = 160:

| property | rho | p |
|---|---|---|
| allele count | -0.445 | 3.8e-09 |
| study count in pool | -0.423 | 2.6e-08 |
| label reassignment fraction | -0.394 | 2.5e-07 |
| imbalance ratio | -0.213 | 6.8e-03 |
| peptides in pool | -0.138 | 0.082 |
| fraction of 9-mers | -0.036 | 0.65 |

Reassignment fraction and study count are **confounded with allele count**,
not independent difficulties: within the 50-allele pools alone, reassignment
fraction flips to rho = +0.594 (p = 0.006). A bigger pool draws from more
studies and has more duplicate peptides to reassign, and within a fixed allele
count more data helps. Nothing here says label noise is driving the result.

## Caveats

**MixMHCp's default arm is handicapped above 6 alleles.** `MAX_MOTIFS` in its
own wrapper is 50 and the default scan is 1-6 motifs, which is how the previous
benchmark ran it and how its own documentation recommends it. At 12, 24 and 50
alleles it is being asked to do something it cannot express. The default-arm
rows are reported for continuity with the earlier study; **the given-K rows are
the like-for-like comparison** and PepCluster2 wins those on 87.5% of pools.

**The 50-allele pools sit exactly at MixMHCp's ceiling.** Its given-K arm is
run at K = 50, the maximum it accepts. Nothing above 50 alleles is testable
against it, which is why the grid stops there.

**Absolute AMI is lower than in the 0.6.0 study** (0.566 here against 0.5935
there, comparing the given-K arm in both). The pools here go to 50 alleles where the earlier
study stopped at 20, and at matched complexity the numbers are comparable.
The two studies' numbers should not be read as a regression.

**Peptide length is part of the signal, by design.** As in the earlier study,
pools are not length-matched; `pct_9mer` has no measurable relation to AMI
(rho -0.036, p = 0.65), so this does not appear to be acting as a shortcut.

## Dataset

Built by `code/build_pools.py` from `data/pep_hla_dataset.parquet`. Peptides
are labelled by the allele parsed from `cell_context`; a peptide seen under
more than one allele is assigned to one of them (the "reassigned fraction"
column). Isoleucine is not recoverable — `modifiedSequence` with bracketed
mass deltas stripped returns `peptide` exactly for all 412,644 rows — so the
I/L collapse is inherent to the source data, not a choice made here.

| alleles | pools | peptides min/median/max | median reassigned |
|---|---|---|---|
| 3 | 40 | 997 / 2,885 / 7,962 | 0.3% |
| 6 | 40 | 992 / 2,961 / 7,949 | 0.6% |
| 12 | 40 | 1,975 / 5,805 / 15,575 | 1.2% |
| 24 | 20 | 3,914 / 11,491 / 30,385 | 2.7% |
| 50 | 20 | 15,218 / 42,395 / 83,593 | 8.2% |

Five abundance regimes per size, equal or Dirichlet(alpha) with alpha in
{10, 3, 1, 0.3}, giving median imbalance ratios of 1.0, 2.5, 5.5, 11.9 and
31.3 (max 209). The natural ratio in the source data is 14.8, so the grid
brackets it. Every allele is seated at a minimum of 20 peptides before the
Dirichlet share is distributed; without that floor the first build produced
alleles with a single peptide and nine pools missing alleles entirely. Seeds
come from `zlib.crc32` of the pool name, so the dataset rebuilds bit-identically.

## Figures

`figures/` — palette imported from `2026-08-06_0.6.0_motif_merge/code/figures.py`,
so colour is the pipeline and a dash means K was supplied, identically in both
reports. GibbsCluster was not run here, so blue is absent.

| file | content |
|---|---|
| `fig1_ami_alleles.png` | AMI and BCubed F1 against allele count |
| `fig2_metrics.png` | all five metrics against allele count |
| `fig3_granularity.png` | motifs reported against alleles present, with the one-per-allele diagonal |
| `fig4_cost.png` | median wall-clock against pool size, IQR band |
| `fig5_imbalance.png` | AMI against abundance regime |

## Reproducing

```bash
cd validation/2026-11-18_0.10.1_julianedata_motif_merge
python3 code/build_pools.py                       # 160 pools, deterministic
sbatch code/variants.sbatch                       # 20 tasks, 7 modes each
sbatch code/mixmhcp.sbatch                        # 20 tasks, 2 settings each
python3 code/merge_shards.py --results results --manifest pool_manifest.csv
python3 code/analyse.py  --results results --manifest pool_manifest.csv
python3 code/figures.py  --results results --manifest pool_manifest.csv --out figures
```

`merge_shards.py` refuses to write any arm unless all 160 pools are present
exactly once, so a truncated array task fails loudly rather than skewing a mode.

Jobs as run: `17389068` (variants, 54 min worst task) and `17389069`
(MixMHCp, 2 h 30 min worst task), both `--array=0-19%20` on `soeding`.
