# Changelog

## 0.7.0

The motif layer settles into three supported variants, and the benchmark is
completed with GibbsCluster and a cost comparison.

### Added

- `--no-motif-merge`: skip the agglomerative merge and give EM one component per
  similarity cluster. Statistically indistinguishable from merging first on AMI
  and BCubed F1 (p = 0.21 and 0.22 paired over 48 independent pools), but finer -
  precision 0.562 against 0.512, recall 0.652 against 0.716. Neither is declared
  correct; the choice is the user's.

- `--motif-count K`: return **exactly** K motifs. Seeds are the K similarity
  clusters farthest apart in profile space, and any component EM empties reclaims
  the peptide that fits it best. Helps above ~12 alleles, where the automatic
  count saturates, and costs nothing below: AMI 0.596 -> 0.620, F1 0.568 ->
  0.600, and +0.043 AMI / +0.057 F1 above twelve alleles.

### Benchmark

GibbsCluster is now included; its parser had never been executed and carried
three independently fatal defects, all fixed. `compare_tools.py` reports BCubed
recall and F1 alongside purity: adjusted purity is maximised by total
fragmentation and cannot compare partitions of different granularity on its own.

Test split, AMI / BCubed F1: merge and refine 0.596 / 0.568, refine only 0.606 /
0.579, given count 0.620 / 0.600, against MixMHCp 0.392 / 0.392 as documented and
0.492 / 0.473 given the true count, and GibbsCluster 0.180 / 0.266 and 0.252 /
0.263.

Cost, serial on an exclusive node: the pipeline is 5x cheaper than MixMHCp and
296x cheaper than GibbsCluster.

### Changed

- Defaults are now calibrated by nested cross-validation rather than provisional.

### Notes

- `--no-motif-em` is retained but marked diagnostic. Merging without refinement
  reaches AMI 0.430 against 0.596 and is not a supported variant; the flag exists
  so the published hyperparameter grid stays reproducible.

- Refinement, not merging, carries the method. The merge is worth about +0.02 AMI
  on the tuning folds and nothing on the independent test set; what it does
  significantly is trade precision for recall by returning a coarser partition.

- The deterministic seed is worth +0.017 AMI over ten random restarts and about
  7.5x less refinement compute, and needs no component count. It is not the
  source of the method's accuracy.

Full methodology in `ALGORITHM.md` Section 15; measurements, figures and
per-panel CSVs in `validation/2026-08-06_0.6.0_motif_merge/`.

## 0.6.0

Adds an optional motif layer above the similarity clustering. It is **off by
default** and does not change any existing output.

### Why

A PepCluster2 cluster is a similarity ball around a representative. A binding
motif is a product of per-position residue preferences: narrow at the anchors,
near-flat elsewhere, and therefore strongly anisotropic in sequence space. A ball
cannot cover such a region, and lowering the threshold widens it along every axis
at once instead of only the tolerant ones. One motif consequently fragments into
many clusters, and no single threshold repairs it.

Measured against MHC allele labels on the 48 held-out benchmark pools, the
similarity clustering alone reaches BCubed recall 0.06 at ~175 clusters: the
clusters are enriched for their allele but far too small. Merging their profiles
and refining with EM raises recall to 0.72 at ~7 motifs, AMI from 0.333 to 0.596
and BCubed F1 from 0.110 to 0.568 - ahead of MixMHCp forced to the true allele
count (0.491 / 0.473) in every complexity band.

### Added

- `--merge-motifs`, building a motif partition from the finished similarity
  clusters. Clusters are summarised as amino-acid profiles on a nine-column frame
  and merged greedily while a Dirichlet-multinomial marginal likelihood prefers
  one shared profile to two separate ones. EM refinement of a mixture of position
  weight matrices then follows, seeded from the merged partition.

- `--motif-prior-concentration`, `--motif-merge-threshold`, `--no-motif-em`,
  `--motif-em-prior-concentration`, `--motif-em-max-iterations`,
  `--motif-em-tolerance`.

- Outputs `motif_clusters.tsv` (motif and similarity cluster per sequence, with
  the framed peptide) and `motif_profiles.tsv` (the fitted profile and mixing
  weight per motif). `run_summary.txt`, `run_config.txt` and `run_stats.json`
  gain the corresponding fields.

### Frame

Nine columns. For `L >= 9` they are peptide positions 1–4 and `L-4..L`, so a
9-mer maps identically and the centre of a longer peptide — which bulges out of
the binding groove and carries little allele-specific signal — is dropped. An
8-mer fills columns 1–4 and 6–9 and leaves column 5 unobserved, rather than
shifting its C-terminal residues into the wrong columns.

### Limits, read these before using it

- **The motif partition does not satisfy the representative-to-member
  invariant.** Two peptides in one motif need not pass the scoring rule against
  any common representative. That is the intent, not a defect, but it is why the
  motif layer is written to separate files and never overwrites `clusters.tsv`.

- **EM does nearly all of the work.** The Bayesian merge is worth about +0.02
  AMI once EM runs; the EM smoothing concentration alone moves AMI by 0.34 across
  its range, and EM performs its own model selection by emptying components. The
  merge stage earns its place on determinism, cost and interpretability, not on
  accuracy. Numbers in `validation/2026-08-06_0.6.0_motif_merge/REPORT.md`.

- **The defaults are calibrated on one dataset.** They were selected by nested
  cross-validation on the inner folds of the peptide-MHC benchmark and evaluated
  once on its 48 held-out test pools. That is a valid protocol but a single
  label universe, so treat them as a starting point on data of a different
  character.

- **Positions are assumed independent** given the motif, so two clusters with
  matching per-position marginals merge even if their joint residue
  distributions are disjoint. This is the standard position weight matrix
  assumption and the reason the stage recovers motifs rather than homologues.

- **Merging can only coarsen.** Contamination already present in a similarity
  cluster propagates; only the EM stage can move a peptide out of it.

- There is no background or "trash" component, so contaminant peptides are forced
  into a real motif.

### Determinism

Preserved. The agglomeration argmax breaks ties by cluster index, EM is seeded
from the merge rather than at random, and parallel accumulation uses fixed chunk
boundaries combined in index order. Output is bit-identical across thread counts.

## 0.5.0

Candidate search redesigned. Recall against an exhaustive all-versus-all
reference rises from 0.64 to about 0.97 while *fewer* candidate pairs are
written to disk than in 0.4.3.

### Fixed

- **Terminal seeding missed roughly a third of eligible pairs.** The index built
  keys from contiguous dimers only (`N1N2`, `N2N3` and `C1C2`, `C2C3`), but a
  pair is accepted on a constrained alignment whose required terminal columns may
  sit at any two of the three positions in each terminal corner — including the
  `(1,3)` combination the seed never indexed. Each terminus was therefore about
  81% sensitive, and a pair had to pass both, giving 0.81² ≈ 0.64 overall. The
  seed now enumerates all three column-pairs per terminus. Peptides of differing
  length were hit hardest, which is 87% of true pairs.

- **Reads of the spilled pair file could return short.** `Read::read` is not
  obliged to fill its buffer; the pair reader assumed it did and aborted with
  "truncated pair record". This affected every run above roughly 500k peptides.
  Records are now read in a fill loop, with a regression test.

### Added

- `separate_kmer_anchor` scoring mode, with `--kmer-similarity-threshold`.
  Similarity is positionwise BLOSUM62 over the first and last three residues,
  normalised per position as `B62(x,y) / sqrt(B62(x,x)·B62(y,y))`; the middle of
  the peptide does not contribute. No alignment is performed, so it is markedly
  faster than `separate_aln_anchor`.

- An anchor upper bound applied before a candidate is spilled to disk. Relaxing
  the one-to-one anchor assignment to a per-row maximum drops injectivity and so
  can only over-estimate the exact score, making rejection provably safe. Costs
  about 36 table lookups and no dynamic programming. Search precision rises from
  1.4% to 3.4%. A test asserts every eligible pair passes the bound.

- `--terminal-seed contiguous|all-column-pairs` (default `all-column-pairs`).
  `contiguous` reproduces 0.4.3 seeding.

- `--representative-order coverage|intrinsic`. `coverage` keeps the dynamic
  greedy set cover. `intrinsic` orders on peptide properties alone, so the
  representative order of a subset is the restriction of the full-data order,
  which makes subset partitions nested-stable at some cost in agreement.

- `--reassignment-margin`, a hysteresis threshold on refinement moves.

### Changed

- `--kmer-seed-threshold` default 0.50 → **0.40**, which is where measured recall
  reaches about 0.97.

Both default changes alter clustering output relative to 0.4.3. To reproduce
0.4.3 exactly, pass `--terminal-seed contiguous --kmer-seed-threshold 0.50`.

### Validation

`validation/2026-07-30_0.5.0-dev_search_redesign/` holds the supporting runs,
including a peptide-MHC benchmark: nested cross-validation over 56,160
configurations per scoring mode, with allele-holdout outer folds. See
`START_HERE.md` there.

## 0.4.3

Previous release.
