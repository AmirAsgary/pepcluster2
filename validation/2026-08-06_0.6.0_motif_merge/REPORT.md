# Motif merge and EM refinement (0.6.0)

Development and calibration record for the optional motif layer, Section 15 of
`ALGORITHM.md`.

## Summary

The similarity clustering fragments each binding motif into roughly 175 clusters
per pool. Those clusters are enriched for their allele but far too small: BCubed
recall 0.06. Merging their amino-acid profiles and refining the result with EM
raises recall to 0.72 at about 7 motifs, and beats MixMHCp on every metric and in
every complexity band.

Two findings qualify that. First, **EM does nearly all of the work**: the
Bayesian merge stage is worth about +0.02 AMI once EM runs, while the EM
smoothing concentration alone moves AMI by 0.34 across its swept range. Second,
**the PepCluster2 seed contributes about +0.017 AMI**, so the pipeline's
advantage is its mixture model, not its initialisation.

## Headline, nested cross-validated

Configuration selected on the inner folds only, evaluated once on the 48
independent test pools. Selection metric AMI; selecting on BCubed F1 instead
chooses the same configuration.

Selected: `--motif-prior-concentration 10 --motif-merge-threshold 25
--motif-em-prior-concentration 3`.

| | AMI | BCubed F1 | Purity | Recall | Clusters |
|---|---:|---:|---:|---:|---:|
| Similarity clustering alone | 0.3333 | 0.1098 | 0.4974 | 0.062 | 174.9 |
| + merge, no EM | 0.4295 | 0.4275 | — | — | 10.8 |
| **+ merge + EM** | **0.5964** | **0.5680** | 0.4617 | 0.7163 | 7.4 |
| MixMHCp (default) | 0.3918 | 0.3918 | 0.2215 | 0.828 | 4.1 |
| MixMHCp (forced k) | 0.4915 | 0.4728 | 0.4178 | 0.484 | 12.4 |

By pool complexity, selected configuration on test:

| Alleles | Pools | AMI | F1 | Purity | Recall | Clusters |
|---|---:|---:|---:|---:|---:|---:|
| 2-6 | 12 | 0.7348 | 0.7978 | 0.7229 | 0.8383 | 5.2 |
| 7-12 | 15 | 0.6321 | 0.6232 | 0.5201 | 0.7298 | 7.8 |
| 13-30 | 21 | 0.4918 | 0.3972 | 0.2707 | 0.6369 | 8.4 |

Per held-out-allele fold, each selected without seeing that fold: AMI 0.607,
0.660, 0.752, 0.722, 0.758.

The earlier exploratory figure of 0.626 was selected on the test split. The
honest number is **0.596**, so that shortcut was worth about +0.03 — a useful
measure of how much such a shortcut inflates a result.

## What each stage actually contributes

Mean inner-fold AMI with EM enabled, at the best EM concentration for each cell:

| merge concentration | threshold 0 | threshold 25 |
|---|---:|---:|
| 0.3 | 0.5657 | 0.5684 |
| 1 | 0.5699 | 0.5681 |
| 3 | 0.5770 | 0.5675 |
| 10 | 0.5841 | **0.5899** |
| 30 | 0.5660 | 0.5823 |

The whole surface spans 0.024. The cell that merges least (concentration 0.3,
threshold 25) leaves 195 of 195 clusters essentially untouched and still reaches
0.5684 once EM runs, against 0.5899 for the best cell. **The merge stage is worth
about +0.02 AMI.**

Without EM it matters much more - merge-only inner AMI ranges 0.276 to 0.392 -
but merge-only is not competitive with MixMHCp anyway.

The EM concentration is the dominant parameter. Within the least-merging cell,
varying it alone gives:

| EM concentration | 0.3 | 1 | 3 | 10 | 30 |
|---|---:|---:|---:|---:|---:|
| AMI | 0.3176 | 0.4184 | **0.5684** | 0.4911 | 0.2454 |
| Clusters | 193.2 | 54.3 | 12.0 | 4.3 | 2.0 |

A 0.34 range, against 0.02 for every merge setting combined. EM also performs its
own model selection by emptying components, so the merge stage is not even
required to choose the component count.

This is a negative result about the part of the design that took the most effort.
The merge stage earns its place on grounds other than accuracy: it is
deterministic, it reduces the number of EM components to fit, and it produces an
interpretable intermediate. It is not what makes the method work.

## Initialisation control

`code/controls.py`. Identical EM, frame and pseudocounts, four initialisations,
48 test pools:

| Condition | k | AMI |
|---|---:|---:|
| merge -> EM | 10.5 | 0.6263 |
| similarity clusters -> EM, no merge | 10.0 | 0.6045 |
| random init at matched k, max-likelihood of 10 restarts | 24 | 0.6094 |
| random init at true allele count, max-likelihood of 10 | 11.4 | 0.6253 |

Restart standard deviation within a pool: 0.027.

The paired gap between the full pipeline and random initialisation is +0.017
(33/48 pools, Wilcoxon p = 4.2e-4) - significant but small - and the pipeline
loses to random best-of-10. **The hypothesis that a deterministic similarity
clustering is what makes deconvolution work at high allele count is not
supported.** An earlier argument for it cited MixMHCp's ±0.203 as evidence of
restart instability; that is variance across pools, not across restarts, and the
two were conflated.

What the seed does buy is determinism - one run instead of ten restarts plus a
selection step - and a component count chosen without being told the answer.

Note these four numbers share the test-set-tuned concentration and so are jointly
optimistic; the comparison between them is unaffected, since all four used
identical hyperparameters.

## Where the advantage over MixMHCp comes from

`code/decompose.py`. Both arms fit the same number of components as the pool has
alleles, from random initialisations, keeping the maximum-likelihood restart, on
identical peptides.

| Peptide set | Ours | MixMHCp | Gap |
|---|---:|---:|---:|
| Full pools | 0.6245 | 0.4915 | +0.133 |
| 9-mers only | 0.5857 | 0.5025 | +0.083 |

Restricting to 9-mers removes 38% of the gap and 62% survives. So the nine-column
frame with its 8-mer gap is worth roughly 0.05 AMI, and the majority of the
advantage is the model itself, persisting where both tools face an identical
9-position problem. MixMHCp improves slightly when non-9-mers are removed, so
those peptides were hurting it.

Smoothing is a powerful lever within our own model - AMI 0.157 to 0.625 across
concentrations 30 to 3 at matched k - but that demonstrates smoothing matters, not
that it is specifically what differs from MixMHCp. Establishing the latter needs
MixMHCp's own pseudocount scheme, which has not been read.

## Pool composition

Relevant because it determines whether length handling is exercised at all.

| Length | 8 | 9 | 10 | 11 | 12-15 |
|---|---:|---:|---:|---:|---:|
| Share | 8.17% | 60.81% | 18.03% | 10.53% | 2.46% |

935,447 peptides over 193 pools. 8-mers occur in **every** pool, up to 25% of one;
no pool is purely 9-mers. Per-pool table in `results/pool_lengths.csv`.

## A metric correction to the 0.5.0 benchmark

`compare_tools.py` reported AMI and adjusted per-allele purity. That purity is
BCubed precision corrected against the allele prior; the correction removes the
baseline for one large cluster but **not** the inflation from fragmentation,
since a singleton scores precision 1 and so adjusted purity 1. A partition into
singletons scores a perfect 1.0. Comparing 175 clusters against 4 on that metric
is not a like-for-like contest, and the earlier reading that PepCluster2 led at
13-30 alleles rested on it.

`metrics.py` already computed BCubed recall and F1; they were absent from
`METRICS` in `compare_tools.py` and so never reached the report. Now added,
together with a corrected description of what purity can and cannot compare.

## Status of these numbers

1. **The headline is properly nested.** Selection used inner folds only;
   test pools were touched once. The compiled defaults are the selected values.
2. **`--motif-merge-threshold 25` sits at the top of its grid.** The surface is
   flat there (0.5899 against 0.5841 at threshold 0), so the boundary is unlikely
   to matter, but it is untested above 25.
3. **The concentration grids are interior.** Merge concentration peaks at 10 of
   {0.3, 1, 3, 10, 30}; EM concentration peaks at 3 and collapses to 0.245 by 30.
   The earlier worry that the optimum lay outside the grid is resolved.
4. **Single dataset**, one label universe, 48 test pools.
5. **No background or outlier component**, so contaminant peptides are forced
   into a real motif.
6. **Positions are assumed independent** given the motif, the standard PWM
   assumption and the reason the stage recovers motifs rather than homologues.

## Rust implementation

`src/motif.rs`, cross-checked against `code/bhc.py` on pools test_01-03 at
concentration 1.0, threshold 0: cluster counts agree exactly at every stage
(106->23->16, 159->46->33, 151->41->32). Post-EM AMI agrees to 0.0005 on test_03
and differs by up to 0.023 on test_01, because the prototype capped EM at 60
iterations while the Rust runs to convergence; the Rust result is the
better-converged one. Output is bit-identical between 1 and 4 threads.

Largest pool (24,888 peptides): 54.5 s single-threaded, 562 MB temporary spill,
407 MB resident.

## Reproducing

```bash
V=validation/2026-08-06_0.6.0_motif_merge
sbatch $V/code/motif_grid.sbatch          # 11,580 runs, ~6 min on 8x64 cores
python3 $V/code/analyse_motif.py --grid $V/grid --out $V/results
python3 $V/code/decompose.py  ...         # see --help
python3 $V/code/controls.py               # initialisation control
```

`runs/assign/` holds per-peptide similarity assignments for the Python analyses;
it is gitignored and regenerates with `code/oracle.py` in about 7 seconds.

Files: `results/motif_selected_overall.csv`, `motif_selected_per_fold.csv`,
`motif_test_by_band.csv`, `decompose.csv`, `controls.csv`, `ceiling.csv`,
`full.csv` (exploratory, superseded by the nested sweep).
