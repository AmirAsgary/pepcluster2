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

## What "PepCluster2 + motif" is

The full pipeline, three stages, run in one command
(`pepcluster2 --merge-motifs ...`):

1. **Similarity clustering.** Sections 1-14 of `ALGORITHM.md`, unchanged.
   Produces ~175 clusters per pool.
2. **Merge.** Each cluster is summarised as amino-acid counts on nine columns and
   pairs are merged greedily while a Dirichlet-multinomial marginal likelihood
   prefers one shared profile to two. Produces ~112 groups at the selected
   settings.
3. **EM.** A mixture of position weight matrices, one component per surviving
   group, fitted by expectation-maximization and **seeded from stage 2 rather
   than at random**. Components that lose all their weight disappear. Produces
   ~7 motifs.

**No k is supplied at any point.** The number of motifs is a consequence of the
merge threshold and of EM emptying components; nothing in the pipeline is told
how many alleles a pool contains. The only arms in this study that receive k are
the `forced k` variants of MixMHCp and GibbsCluster, which exist precisely to
show how much of those tools' results come from their model rather than their
model selection.

"Motif layer" means stages 2 and 3 together. It is reported as a second,
separate partition: it does not satisfy the representative-to-member invariant
that stage 1 guarantees, so it is written to its own files.

### The two concentration parameters

Both are Dirichlet pseudocount totals per column, spread over the dataset's
background residue frequencies as `alpha_a = concentration * 20 * background_a`.
They differ only in where they are used:

| Flag | Used by | Effect |
|---|---|---|
| `--motif-prior-concentration` | stage 2, the merge Bayes factor | higher smooths the profiles being compared, so more pairs look alike and more merges happen |
| `--motif-em-prior-concentration` | stage 3, smoothing the PWMs each EM iteration | higher pulls components toward background, collapsing more of them; this is the dominant parameter |

They were swept independently precisely because an earlier exploratory sweep tied
them to one value and could not attribute the result between the two stages.

### It is a Dirichlet prior, not a Dirichlet mixture

A single Dirichlet per column, scaled by background composition. A Dirichlet
*mixture* prior - several components encoding amino-acid classes, which would
make the merge aware that D and E are chemically alike - was considered and
deliberately **not** implemented. The criterion is therefore blind to amino-acid
similarity: it treats the 20 residues as unordered categories, so an all-D column
and an all-E column are as different as all-D and all-W. That is a known
limitation, not an oversight.

### The Chinese-restaurant prior was removed

A size-scaled partition prior of Chinese-restaurant form was implemented,
measured, and rejected: it corrupted the merge order and collapsed the
dendrogram (numbers below). **It never entered the Rust implementation**, and its
prototype code has now been deleted. The rationale is retained in `ALGORITHM.md`
Section 15.4 so the same design is not attempted again; the evidence is in
`results/full.csv`.

## Three variants

All share the base clustering, the frame and the EM stage; they differ only in
how EM is seeded. Test split, 48 independent pools:

| Variant | Flags | AMI | Purity | Precision | Recall | F1 | Motifs |
|---|---|---:|---:|---:|---:|---:|---:|
| Merge and refine | `--merge-motifs` | 0.5964 | 0.4617 | 0.5115 | 0.7163 | 0.5680 | 7.4 |
| Refine only | `+ --no-motif-merge` | 0.6056 | 0.5177 | 0.5618 | 0.6519 | 0.5785 | 10.0 |
| Given count | `+ --motif-count K` | 0.6201 | 0.5106 | 0.5564 | 0.6897 | 0.5998 | exactly K |

Merge-and-refine and refine-only are statistically indistinguishable on AMI and
F1 (p = 0.21 and 0.22 paired); they differ in granularity, refine-only trading
recall for precision. Both are supported and neither is declared correct.

Reproducibility was verified for all three: output is bit-identical across 1, 4
and 16 threads and across repeated runs, and a requested count is returned
exactly for every K tried from 2 to 40.

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

### Corrected on the test pools, and on every metric

The `+0.02` above is an inner-fold figure read off AMI alone. Measured properly -
paired on the 48 test pools, selected configuration against near-zero merging,
both with EM - it does not replicate, and the picture on the other metrics is
different in kind:

| Metric | merge + EM | EM, minimal merging | Delta | Wilcoxon p |
|---|---:|---:|---:|---:|
| AMI | 0.5964 | 0.6074 | -0.0110 | 0.21 |
| BCubed F1 | 0.5680 | 0.5801 | -0.0121 | 0.22 |
| Adjusted purity | 0.4617 | 0.5195 | **-0.0578** | <0.0001 |
| BCubed precision | 0.5115 | 0.5632 | **-0.0517** | <0.0001 |
| BCubed recall | 0.7163 | 0.6538 | **+0.0625** | <0.0001 |
| Clusters | 7.42 | 9.92 | -2.50 | <0.0001 |

So the merge stage does **not** improve overall agreement - AMI and F1 are
statistically indistinguishable, and both point slightly the wrong way. What it
does, highly significantly, is trade precision for recall by returning a coarser
partition. **It is a granularity control, not a quality improvement**, and the
earlier "+0.02 AMI" reading was inner-fold noise that reversed sign on the test
pools.

Reporting this on AMI alone would have hidden it in both directions: AMI shows no
effect where precision and recall each move by more than 0.05, because a
chance-corrected summary is insensitive to a shift that moves both in
compensating directions.

The merge stage still earns its place on determinism, on the number of components
EM must fit, and on producing an interpretable intermediate. It is not what makes
the method work, and it does not make it more accurate.

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

What the seed buys, stated precisely, is three things:

1. **Accuracy, slightly.** +0.017 AMI over the max-likelihood of ten random
   restarts, paired, p = 4.2e-4.
2. **Compute.** Measured on three test pools spanning 995 to 11,656 peptides:
   the merge costs 0.23 s and one EM run 0.70 s, so the whole motif stage is
   0.93 s. Matching it from random initialisation needs ten EM restarts plus a
   likelihood comparison, about 7.0 s - roughly **7.5x the motif-stage compute
   for a slightly worse result**. (Timings in
   `figures/fig1c_stage_contribution_all_metrics.csv` context; raw in the
   commit message.)
3. **A component count, unprompted.** The random-init arms had to be *given* k -
   either ours to match, or the true allele count. Nothing supplies it here.

So the earlier phrasing "the seed is not where the accuracy comes from" is right
about accuracy and understated everything else. It is better read as: the seed is
not what makes the method *work*, but it makes it cheaper, reproducible, and
free of a parameter the alternatives require.

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

## GibbsCluster

Now included. It had never been run: the tool sits behind a DTU licence form, and
the parser written against its documentation could not have worked. Three
independent defects, all fixed:

- it looked for a comment-prefixed header, but the header is the uncommented
  first line;
- it searched for a column named `cluster` or `group`, but the group column is
  `Gn` and is 0-based;
- it chose the number of groups from the last column of
  `images/gibbs.KLDvsClusters.tab`. That file is a zero-padded matrix - rows are
  group counts, columns are per-cluster KLD - so the last column is usually
  padding. The criterion is the row sum. On a test case the old rule would have
  chosen 4 groups (row sum 6.42) over the correct 2 (7.63).

Runtime was the other obstacle. A 972-peptide pool did not finish in 35 minutes,
because `-k` defaults to 1. GibbsCluster forks over `seeds x group counts`, which
is 30 independent fits for the default preset, and with `-k 16` the mean run over
all 386 is **106 s** (max 2219 s). All 386 completed; none timed out.

Test split, 48 pools:

| Tool | AMI | Purity | Recall | F1 | Clusters |
|---|---:|---:|---:|---:|---:|
| GibbsCluster (default) | 0.1799 | 0.1000 | 0.6046 | 0.2656 | 2.9 |
| GibbsCluster (forced k) | 0.2522 | 0.1861 | 0.2576 | 0.2628 | 12.4 |
| MixMHCp (default) | 0.3918 | 0.2215 | 0.8282 | 0.3918 | 4.1 |
| MixMHCp (forced k) | 0.4915 | 0.4178 | 0.4843 | 0.4728 | 12.4 |
| PepCluster2, similarity only | 0.3333 | 0.4974 | 0.0624 | 0.1098 | 174.9 |
| **PepCluster2 + motif layer** | **0.5964** | 0.4617 | 0.7163 | **0.5680** | 7.4 |

AMI by complexity:

| Alleles | Gibbs def | Gibbs forced k | MixMHCp def | MixMHCp forced k | PC2 sim | PC2 + motif |
|---|---:|---:|---:|---:|---:|---:|
| 2-6 | 0.348 | 0.378 | 0.573 | 0.712 | 0.337 | 0.735 |
| 7-12 | 0.185 | 0.249 | 0.441 | 0.535 | 0.348 | 0.632 |
| 13-30 | 0.080 | 0.183 | 0.253 | 0.335 | 0.321 | 0.492 |

GibbsCluster trails MixMHCp throughout and degrades faster with allele count.
Against the *similarity* clustering alone it wins below 7 alleles and loses above,
the same crossover MixMHCp shows but at a lower level.

**Coverage caveat.** At motif length 9 GibbsCluster cannot place a core in an
8-mer, and its trash cluster removes outliers, so mean coverage is 0.849 (minimum
0.554). Every run is therefore scored twice: over the full pool with unassigned
peptides collected into one cluster, and over assigned peptides only. The two
barely differ - assigned-only AMI is 0.186 (default) and 0.266 (forced k) against
0.180 and 0.252 - so coverage is *not* what places it behind MixMHCp.

**Trash threshold, settled.** The `-j 2` used here is more aggressive than the
tool's own default of 0. Paired over all 96 test runs:

| Setting | AMI at -j 2 | AMI at -j 0 | Coverage -j 2 | Coverage -j 0 |
|---|---:|---:|---:|---:|
| default | 0.1799 | 0.1568 | 0.771 | 0.855 |
| forced k | 0.2522 | 0.2521 | 0.912 | 0.925 |

Wilcoxon p = 0.068, and the tool's own default is if anything *worse*. The choice
does not affect any conclusion, and `--trash-threshold` is now a flag rather than
a hard-coded constant.

## How granularity responds to the pool

Spearman rho between the number of clusters returned and the pool's properties,
on the test pools; the whole family of 252 tests is Benjamini-Hochberg corrected
together (`figures/correlations_spearman.csv`).

| Tool | vs peptides | vs alleles |
|---|---:|---:|
| PepCluster2 + motif | 0.766 *** | **0.444 **** |
| PepCluster2 (similarity) | 0.906 *** | 0.437 ** |
| MixMHCp (default) | 0.678 *** | -0.298 (ns) |
| GibbsCluster (default) | 0.163 (ns) | **-0.664 **** |
| MixMHCp / GibbsCluster (forced k) | 0.103 (ns) | 1.000 *** |

The forced-k rows returning exactly 1.000 against allele count is the sanity
check that the harness does what it claims.

The informative row is the last unforced one. **GibbsCluster's own model
selection returns *fewer* motifs as a pool gains alleles**, strongly and
significantly so, and MixMHCp's is flat. Only the motif layer's component count
rises with the number of alleles present without being told it. That is a
statement about model selection, not about partition quality, and it is the
mechanism behind the collapse of both external tools in the 13-30 band.

Rank correlation is used because these are counts spanning an order of magnitude
with no reason to be linear or normally distributed. Pools are treated as
independent, which is approximate: they are drawn from a shared peptide universe
and some alleles recur across pools, so the p-values are mildly optimistic.

## Figures and tables

In `figures/`:

- `fig1_protocol` - the selection surface on the tuning folds, each held-out
  fold scored by a choice made without it, and what each stage adds
- `fig2_granularity` - clusters returned and mean cluster size, against pool size
  and against allele count
- `fig3_performance_alleles`, `fig4_performance_size` - AMI, purity, precision,
  recall and F1 against allele count and pool size
- `fig5_all_splits` - every metric on every split
- `per_pool_all_tools.csv` - the per-pool table behind all of it
- `table1_benchmark_test.csv`, `table2_test_by_allele_band.csv`,
  `table3_all_splits.csv`, `correlations_spearman.csv`,
  `stage_contribution_all_metrics.csv`, `stage_deltas_all_metrics.csv`

Colour is fixed throughout: black PepCluster2, red MixMHCp, blue GibbsCluster;
solid as documented, dashed given the true allele count.

## Why the motif count saturates near ten

The component count tracks the allele count to about ten and then flattens and
falls: 5.2 motifs at 5 alleles, 9.3 at 10, 8.8 at 15, 8.0 at 19. Two causes, and
both are real.

**Allele motifs genuinely converge as pools grow.** Jensen-Shannon divergence
between the *true* per-allele profiles, built from the labels, over 3,470 pairs:

| percentile | p01 | p10 | p50 | p90 | p99 |
|---|---:|---:|---:|---:|---:|
| bits | 0.057 | 0.085 | 0.140 | 0.193 | 0.238 |

It is a smooth continuum with no gap, so there is no natural cut between "same
motif" and "different motif" and any count of distinct motifs is threshold-
dependent. What is not threshold-dependent is the trend: median divergence falls
from 0.170 bits at 5 alleles to 0.130 at 19, and the closest pair in a pool falls
from 0.110 to 0.058. Larger pools necessarily contain more closely related
alleles, so part of the saturation is the target genuinely becoming less
separable. At a 0.15-bit complete-linkage cut the distinguishable count is 3.3 at
5 alleles rising only to 4.6 at 19 - the same saturation shape our method shows.

An earlier version of this analysis used *single* linkage and reported 1.6
distinguishable motifs from 11 alleles. That was wrong: single linkage chains, so
one path of near-neighbours collapsed everything. Complete linkage is used
instead, and the pair-level percentiles above are the honest summary.

**Our model selection also under-splits, and that part is fixable.** Seeding EM
with the true allele count instead of letting the merge choose:

| Alleles | Automatic k | Automatic AMI | Forced-k AMI | Automatic F1 | Forced-k F1 |
|---|---:|---:|---:|---:|---:|
| 4-6 | 5.2 | 0.735 | 0.725 | 0.798 | 0.789 |
| 7-9 | 6.5 | 0.619 | 0.627 | 0.625 | 0.633 |
| 10-12 | 9.3 | 0.647 | 0.646 | 0.621 | 0.625 |
| 13-16 | 8.8 | 0.527 | 0.555 | 0.447 | 0.482 |
| 17-20 | 8.0 | 0.453 | 0.493 | 0.342 | 0.393 |

Overall +0.013 AMI (p = 0.13) and +0.018 F1 (p = 0.057) - not significant. Above
twelve alleles, where the automatic count has saturated, +0.033 AMI and +0.043
F1. Below twelve it costs nothing. So the automatic count is adequate for simple
pools and leaves something on the table for complex ones.

`--motif-count K` is therefore provided. It seeds EM with one component per each
of the K largest similarity clusters and bypasses the merge; every other peptide
is placed by the first E-step. EM may still empty a component, so the reported
count can come out below K.

### K is strict

`--motif-count K` returns exactly K motifs, on 48 of 48 test pools. A user who
asks for K wants K; the tool does not get to return fewer because it prefers its
own answer.

Delivering that needs one step beyond seeding, because EM on its own merges
components the data does not separate: a request averaging 11.4 converged to 8.8.
Any component EM empties therefore reclaims the peptide with the highest
likelihood under it, drawn from a component that can spare one. Ties resolve to
the lowest node index, so the result stays reproducible.

The cost of strictness was expected to be real and turned out to be negligible:
AMI 0.6201 strict against 0.6209 non-strict, F1 0.5998 against 0.5991 - well
inside noise. The reclaimed motifs are ones the likelihood did not support alone,
so on data where the requested count genuinely exceeds the number of separable
motifs some returned motifs will resemble each other. That is the honest price of
a strict count, and it is the caller's choice to pay it: the automatic path never
invokes this and still reports whatever the data supports.

Four seedings were implemented and measured. The differences are large, so this
is the part of the feature that needed care:

| Seeding (before strictness was added) | Delivered of 11.4 | AMI | F1 |
|---|---:|---:|---:|
| Merge down to K, then EM from those groups | 7.2 | 0.5946 | 0.5643 |
| Top-K largest, remainder folded into last component | 6.1 | 0.5213 | 0.5031 |
| Top-K largest, remainder unseeded | 8.4 | 0.6054 | 0.5811 |
| **Farthest-apart clusters, remainder unseeded** | **8.8** | **0.6209** | **0.5991** |
| the above, plus strict reclaim | **11.4 = requested** | 0.6201 | 0.5998 |

Merging first is worst because the merge has already blended its groups toward
one another. Folding the remainder into the last component makes that component a
bin of everything unlike the others.

Seeding by *size* is the subtle failure. The largest similarity clusters are not
the most distinct ones: measured on a 20-allele pool, the twenty largest clusters
had median pairwise Jensen-Shannon divergence 0.098 bits against 0.153 for the
true allele profiles, and 54% of seed pairs fell below 0.10 bits. Handing EM
near-duplicate seeds guarantees it will merge them. Choosing seeds farthest apart
in profile space instead - a farthest-first traversal over clusters above the
median size - is worth +0.014 AMI and +0.016 F1 over the size rule.

Against the automatic count, over 48 test pools: AMI 0.5964 -> 0.6201
(p = 0.010), F1 0.5680 -> 0.5998 (p = 0.002), and +0.043 AMI / +0.057 F1 above
twelve alleles.

This is not an oracle setting in practice. A sample's alleles are normally known
from HLA typing, so supplying the count is ordinary use - unlike the `forced k`
arms of MixMHCp and GibbsCluster in the benchmark, which are reported separately
precisely because the benchmark treats the allele count as unknown.

## Cost

Serial runs on an **exclusive** node, nine test pools log-spaced from 995 to
11,656 peptides. The elapsed times recorded during the sweeps are not usable for
this: those ran 40-64 jobs at once, so their wall clocks include contention.

A first attempt at this table was itself contaminated. It ran on a shared login
node while a GibbsCluster sweep was using about 48 cores, and one pool came out
at 13.5 s against 2.2 s on a clean re-run - a spike that looked like a scaling
cliff and was not. The numbers below come from a job holding a whole node. Cost
is now monotone in pool size for every PepCluster2 arm, Spearman rho = 1.000.

PepCluster2 and MixMHCp are single-threaded here, so wall time is also CPU cost.
GibbsCluster is impractical single-threaded and is run at `-k 16`, so its cost is
wall x 16. CPU seconds is the comparable figure; wall seconds is what a user
waits.

| Arm | Median CPU s | Relative | Mean F1 | Mean AMI |
|---|---:|---:|---:|---:|
| PepCluster2, clustering only | 0.88 | 0.57x | 0.099 | 0.316 |
| **PepCluster2 + merge + EM** | **1.55** | **1.00x** | **0.519** | **0.560** |
| PepCluster2 + merge + EM, given k | 0.85 | 0.55x | 0.570 | 0.599 |
| MixMHCp | 7.87 | 5.06x | 0.335 | 0.332 |
| GibbsCluster | 459.52 | 296x | 0.229 | 0.173 |

The motif layer roughly doubles the cost of clustering alone and takes F1 from
0.099 to 0.519. Against the external tools the pipeline is 5x cheaper than
MixMHCp and 296x cheaper than GibbsCluster while scoring above both. Supplying k
is *cheaper* than the automatic count - it skips the agglomeration entirely, and
EM converges faster from diverse seeds than from merged groups.

These F1 and AMI means are over the nine timing pools only, so they differ from
the 48-pool benchmark figures elsewhere in this report; they are here to place
each arm on the cost axis, not to restate accuracy.

Cost rises with pool size for every tool (Spearman rho 0.90 to 1.00, all
significant). Nine pools is too few to separate the scaling exponents, so no
claim is made about asymptotic complexity - only about cost in this range. The
MixMHCp curve is not monotone because it selects its own motif count, and how
many it settles on changes its runtime.

`figures/fig6_speed_vs_pool_size` and `fig7_speed_vs_performance`.

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

`src/motif.rs`, cross-checked against `results/full.csv` on pools test_01-03 at
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
