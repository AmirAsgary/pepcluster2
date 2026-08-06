# Motif merge and EM refinement (0.6.0)

Development record for the optional motif layer, Section 15 of `ALGORITHM.md`.

**Read the caveats in "Status of these numbers" before quoting anything here.**

## Question

The peptide-MHC benchmark in `../2026-07-30_0.5.0-dev_search_redesign/benchmark/`
showed PepCluster2 losing to MixMHCp on AMI at low allele count and drawing level
only above about 13 alleles. Diagnosing why produced this stage.

## What the benchmark comparison actually said

`compare_tools.py` reported AMI and an adjusted per-allele purity. The purity is
BCubed precision corrected against the allele prior only; that correction removes
the baseline for one large cluster but **not** the inflation from fragmentation,
since a singleton scores precision 1 and therefore adjusted purity 1. Comparing
175 clusters against 4 on that metric is not a like-for-like contest.

`metrics.py` already computes the counterpart that penalises fragmentation, but it
was not in `compare_tools.py`'s `METRICS` list and so never reached the report.
Recomputed on the 48 test pools (`code/fair.py`):

| Tool | AMI | Adj. purity | BCubed recall | BCubed F1 | Clusters |
|---|---:|---:|---:|---:|---:|
| MixMHCp (default) | 0.392 | 0.222 | 0.828 | 0.392 | 4.1 |
| MixMHCp (forced k) | 0.491 | 0.418 | 0.484 | 0.473 | 12.4 |
| PepCluster2 k-mer | 0.333 | 0.497 | 0.062 | 0.110 | 174.9 |
| PepCluster2 alignment | 0.319 | 0.446 | 0.061 | 0.106 | 139.2 |

Recall 0.062 is roughly `1 / (clusters per allele)`. The clusters are enriched for
their allele — precision runs about eight times chance at high complexity — but
far too small. On F1 MixMHCp leads in every complexity band, so the earlier
reading that PepCluster2 led at 13–30 alleles does not hold.

**Action taken:** add `bcubed_recall_macro` and `bcubed_f1_macro` to `METRICS` in
`benchmark/code/compare_tools.py`. Not yet done.

## Ceiling on merging

`code/ceiling.py` assigns each similarity cluster to its majority allele and
merges clusters sharing one. This uses the labels, so it is not a method; it
bounds any procedure that only coarsens.

| | AMI | BCubed F1 | Clusters |
|---|---:|---:|---:|
| PepCluster2 as-is | 0.333 | 0.110 | 174.9 |
| Oracle-merged | 0.494 | 0.498 | 11.4 |

Enough headroom to justify building the stage.

## Merge and refinement

`code/full.py` sweeps the Dirichlet concentration and merge threshold, with EM
after each. Best configuration found (`prior_concentration = 3`, threshold 0):

| Alleles | MixMHCp def | MixMHCp forced k | PC2 as-is | + merge | + merge + EM | Oracle |
|---|---:|---:|---:|---:|---:|---:|
| 2–6 | 0.573 | 0.712 | 0.337 | 0.496 | 0.734 | 0.592 |
| 7–12 | 0.441 | 0.535 | 0.348 | 0.446 | 0.658 | 0.510 |
| 13–30 | 0.253 | 0.335 | 0.321 | 0.376 | 0.542 | 0.425 |
| overall | 0.392 | 0.491 | 0.333 | 0.428 | 0.626 | 0.494 |

Recall rises 0.062 → 0.656 with purity roughly held; F1 0.110 → 0.602. EM exceeds
the merge-only ceiling on 45 of 48 pools, as it must, since merging cannot move a
peptide between clusters and EM can.

## Chinese-restaurant partition prior: rejected

`code/bhc.py` implements Bayesian hierarchical clustering (Heller & Ghahramani
2005) with the Dirichlet-process prior. It fails here:

| | merged k | merged AMI | after EM |
|---|---:|---:|---:|
| BHC, alpha = 1e2 | 10.8 | 0.162 | 0.459 |
| BHC, alpha = 1e6 | 21.4 | 0.162 | 0.474 |
| plain Bayes factor | 24.0 | 0.428 | 0.626 |

Cluster count looks reasonable while the partition carries no information, which
is order corruption rather than a bad cut point. The DP prior contributes about
`n log 2` nats toward merging two clusters of size `n`, against a likelihood term
of `n` times the per-peptide divergence between the two motifs; over nine columns
that divergence is often below `log 2` for alleles of one supertype, so the prior
outvotes the data. BHC was validated at n in the tens to hundreds; at n ≈ 10,000
`Gamma(n_k)` pins `pi_k` at 1.

The shipped criterion uses a flat per-cluster penalty instead: requiring
`log BF > t` is a partition prior proportional to `exp(-t k)`, which does not
scale with cluster size. Because `t` is constant across candidate pairs it never
perturbs the merge order.

## Initialisation control: the seed is not where the accuracy comes from

`code/controls.py`. Identical EM, frame and pseudocounts from four inits:

| Condition | k | AMI | F1 |
|---|---:|---:|---:|
| A — merge → EM | 10.5 | 0.626 | 0.602 |
| C — similarity clusters → EM, no merge | 10.0 | 0.604 | 0.577 |
| B — random init at A's k, max-likelihood of 10 restarts | 24 → | 0.609 | 0.567 |
| B — best of 10 | | 0.633 | |
| D — random init at true allele count, max-likelihood of 10 | 11.4 | 0.625 | 0.565 |

Within-pool restart standard deviation: 0.027.

A beats B's max-likelihood restart by +0.017 (33/48 pools, Wilcoxon p = 4.2e-4) —
significant but small — and loses to B's best-of-10. Condition D, with no
PepCluster2 anywhere, reaches 0.625 against A's 0.626.

So the hypothesis that a deterministic similarity clustering is what makes
mixture deconvolution work at high allele count is **not supported**. An earlier
argument for it cited MixMHCp's ±0.203 as evidence of restart instability; that
figure is variance across pools, not across restarts, and the two were conflated.

What the seed does buy is determinism — one run against ten restarts plus a
selection step — and an automatic component count, since B needed k supplied and
D needed the true allele count.

The gain over MixMHCp is therefore in the mixture model, not the pipeline around
it: D against MixMHCp forced-k is +0.134 at matched k and matched random
initialisation. That difference is **not yet decomposed**. Candidates are the
pseudocount strength, the nine-column frame, and the 8-mer gap handling.

## Status of these numbers

1. **The configuration was selected on the test split.** `prior_concentration = 3,
   threshold = 0` won a 9-cell grid scored on the same 48 pools reported above,
   and the control experiment inherited those values. Every absolute figure here
   is optimistic. A nested cross-validated selection on the inner folds is
   required before any of this is quotable, and the compiled defaults are marked
   provisional for that reason.
2. **The optimum sits at the grid boundary.** Performance rose monotonically over
   concentrations 0.3 → 1 → 3, and 3 was the largest tested.
3. **The concentration is confounded.** The sweep passed one value to both the
   merge likelihood and the EM smoothing, so the gain cannot be attributed
   between the two stages. They are separate flags in the shipped code
   (`--motif-prior-concentration`, `--motif-em-prior-concentration`) and must be
   swept independently.
4. **Single dataset**, 48 test pools, one label universe.
5. **8-mer handling differs from the prototype.** `code/merge.py` dropped 8-mers
   from the profiles entirely; `code/bhc.py` and the Rust implementation give them
   a gap at column 5. Numbers from `merge.py` are therefore slightly pessimistic.
6. **GibbsCluster is still missing** from the comparison; see the benchmark
   README for the licence step.

## Rust implementation

`src/motif.rs`. Cross-checked against `code/bhc.py` on pools test_01–03 at
concentration 1.0, threshold 0: cluster counts agree exactly at every stage
(106→23→16, 159→46→33, 151→41→32). Post-EM AMI agrees to 0.0005 on test_03 and
differs by up to 0.023 on test_01, because the prototype capped EM at 60
iterations while the Rust runs to convergence (156 iterations on test_01); the
Rust result is the better-converged one. Output is bit-identical between 1 and 4
threads.

## Reproducing

`runs/assign/` holds the per-peptide similarity assignments the Python analysis
reads. It is gitignored; regenerate it with `code/oracle.py` in about 7 seconds.
Then `code/full.py`, `code/controls.py`, `code/report.py`.
