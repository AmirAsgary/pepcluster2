# PepCluster2 candidate search and clustering validation

## Configuration

- 20 independently sampled datasets, 10,000 peptides each.
- Scoring mode `separate_aln_anchor`; alignment threshold 0.50; anchor-combination threshold 0.60; terminal/core weights 4/1.
- Terminal seed `all-column-pairs`, k-mer seed threshold 0.40 unless a sweep row says otherwise.
- Reassignment margin 0.01: a peptide leaves its representative only when
  another beats it by more than that. Both references use the same value, so
  runs and references remain comparable.

## Method under test

1. The terminal seed indexes all three ordered column pairs of each terminal
   3-mer, (1,2), (1,3) and (2,3). The constrained alignment must contain at
   least `--minimum-terminal-match-length` matched columns drawn from the first
   three residues of *both* peptides, and those columns shift when the peptides
   differ in length, so all three pairs are required to retrieve them.
2. Candidate generation applies a sound upper bound on the anchor-combination
   similarity (the assignment relaxed to independent row maxima). A pair failing
   it cannot be accepted, so this prunes without losing any eligible
   relationship.
3. Reassignment is hysteretic: a peptide leaves its representative only when
   another beats it by more than the margin, so near-ties do not flip when the
   dataset composition changes.

## Reading the metrics

- Search recall: fraction of exactly eligible pairs the run scored.
- Search precision: fraction of scored pairs that were eligible. Low
  precision means wasted work, not wrong edges.
- Cost is reported as a decomposition: index hits (with multiplicity),
  pairs rejected by the sound anchor bound, distinct pairs exactly scored,
  and constrained-alignment evaluations. These are different costs.
  Alignment evaluations are counted process-wide, so they include
  representative recalculation, merge validation and reassignment; the count
  can therefore exceed the number of distinct candidate pairs.
- The exhaustive reference scores every pair exactly and then applies the
  identical clustering procedure, for the same representative order as the run
  under test. A run therefore differs from it only through candidate search,
  which is what makes the agreement number interpretable.

## Seed geometry and sensitivity sweep

Graph method, coverage order, 5 datasets.

| Terminal seed | Seed threshold | Recall | Precision | Index hits | Pairs scored | Alignments | ARI vs reference | Seconds | Peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| contiguous | 0.50 | 0.6430 | 0.0694 | 1,003,038 | 166,900 | 180,870 | 0.4298 | 0.76 | 36 |
| all-column-pairs | 0.50 | 0.9423 | 0.0474 | 4,176,869 | 358,303 | 283,741 | 0.7489 | 1.10 | 57 |
| all-column-pairs | 0.45 | 0.9464 | 0.0465 | 4,436,888 | 366,282 | 281,476 | 0.7563 | 1.10 | 60 |
| all-column-pairs | 0.40 | 0.9695 | 0.0358 | 9,253,657 | 487,334 | 320,923 | 0.8065 | 1.46 | 84 |
| all-column-pairs | 0.35 | 0.9865 | 0.0304 | 16,350,032 | 584,385 | 346,585 | 0.8728 | 1.77 | 109 |
| all-column-pairs | 0.30 | 0.9917 | 0.0276 | 25,082,015 | 646,097 | 367,473 | 0.9127 | 1.85 | 140 |

## Search-rule performance

| Method | Order | Recall | Precision | Pairs scored | Alignments | All-pairs fraction | Seconds |
|---|---|---:|---:|---:|---:|---:|---:|
| Graph | coverage | 0.9701 ± 0.0014 | 0.0358 | 486,814 | 312,769 | 0.0097 | 1.17 |
| Graph + prefilter | coverage | 0.9701 ± 0.0013 | 0.0295 | 752,835 | 335,737 | 0.0118 | 1.36 |
| Greedy | coverage | 0.7198 ± 0.0061 | 0.0128 | 2,272,308 | 1,236,388 | 0.0202 | 12.41 |
| Greedy lazy-exact | coverage | 0.9705 ± 0.0013 | 0.0163 | 2,716,011 | 1,311,917 | 0.0215 | 16.41 |
| Graph | intrinsic | 0.9701 ± 0.0014 | 0.0358 | 486,814 | 344,267 | 0.0097 | 1.28 |
| Graph + prefilter | intrinsic | 0.9701 ± 0.0014 | 0.0295 | 752,835 | 367,248 | 0.0118 | 1.42 |
| Greedy | intrinsic | 0.7329 ± 0.0065 | 0.0108 | 2,519,583 | 1,129,933 | 0.0244 | 15.08 |

## Missed-pair audit

Every eligible pair a run failed to score, attributed to the cause.
`seed attributable` means a terminal seed found no neighbouring column
pair. `retrievable but unexamined` means the index would have returned
the pair and the bound would have kept it, but the clustering traversal
never scored it: the greedy paths only ever score
representative-to-unassigned pairs, so for them this is structural and not
a candidate-search defect. `anchor bound unsound` counts pairs the sound
bound rejected; it must be zero, and any non-zero value is a correctness
failure of the bound.

| Method | Order | Missed | Seed attributable | Front seed | End seed | Both | Retrievable but unexamined | Anchor bound unsound |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Graph | coverage | 537 | 537 | 271 | 266 | 0 | 0 | 0 |
| Graph + prefilter | coverage | 537 | 537 | 271 | 266 | 0 | 0 | 0 |
| Greedy | coverage | 5037 | 527 | 266 | 261 | 0 | 4509 | 0 |
| Greedy lazy-exact | coverage | 529 | 529 | 267 | 262 | 0 | 0 | 0 |
| Graph | intrinsic | 537 | 537 | 271 | 266 | 0 | 0 | 0 |
| Graph + prefilter | intrinsic | 537 | 537 | 271 | 266 | 0 | 0 | 0 |
| Greedy | intrinsic | 4801 | 526 | 266 | 260 | 0 | 4275 | 0 |

## Agreement with the exhaustive reference

The reference scores every pair exactly and then runs the identical clustering
procedure, so the gap below is candidate-search loss and nothing else.

| Method | Order | ARI | NMI | Jaccard | Co-assoc. recall | Clusters |
|---|---|---:|---:|---:|---:|---:|
| Graph | coverage | 0.8053 ± 0.0206 | 0.9741 | 0.6746 | 0.7948 | 3894 |
| Graph + prefilter | coverage | 0.8063 ± 0.0209 | 0.9742 | 0.6760 | 0.7956 | 3894 |
| Greedy | coverage | 0.4361 ± 0.0136 | 0.9324 | 0.2791 | 0.3699 | 4384 |
| Greedy lazy-exact | coverage | 0.8055 ± 0.0206 | 0.9741 | 0.6750 | 0.7951 | 3893 |
| Graph | intrinsic | 0.9480 ± 0.0081 | 0.9955 | 0.9013 | 0.9402 | 4776 |
| Graph + prefilter | intrinsic | 0.9480 ± 0.0080 | 0.9955 | 0.9012 | 0.9401 | 4776 |
| Greedy | intrinsic | 0.9493 ± 0.0078 | 0.9956 | 0.9036 | 0.9422 | 4775 |

## Stability

A subset clustering is compared with the full-dataset clustering restricted
to the same peptides. The `Exhaustive reference` rows apply the identical
comparison to the reference itself, so they show how much of the composition
dependence belongs to the clustering procedure rather than to candidate
search. Graph, graph + prefilter and lazy-exact greedy agree to within 0.0004
throughout, so their curves coincide in the figure.

| Method | Order | Subset | Jaccard | ARI | Co-assoc. recall | Co-assoc. precision |
|---|---|---:|---:|---:|---:|---:|
| Graph | coverage | 10% | 0.2365 | 0.3808 | 0.4171 | 0.3527 |
| Graph | coverage | 80% | 0.4590 | 0.6288 | 0.6444 | 0.6142 |
| Graph + prefilter | coverage | 10% | 0.2360 | 0.3803 | 0.4166 | 0.3521 |
| Graph + prefilter | coverage | 80% | 0.4591 | 0.6289 | 0.6445 | 0.6143 |
| Greedy | coverage | 10% | 0.2656 | 0.4186 | 0.5125 | 0.3552 |
| Greedy | coverage | 80% | 0.5266 | 0.6894 | 0.7164 | 0.6646 |
| Greedy lazy-exact | coverage | 10% | 0.2364 | 0.3807 | 0.4168 | 0.3527 |
| Greedy lazy-exact | coverage | 80% | 0.4588 | 0.6287 | 0.6443 | 0.6140 |
| Exhaustive reference (no merging) | coverage | 10% | 0.2345 | 0.3787 | 0.4181 | 0.3480 |
| Exhaustive reference (no merging) | coverage | 80% | 0.4514 | 0.6216 | 0.6373 | 0.6070 |
| Exhaustive reference | coverage | 10% | 0.2357 | 0.3803 | 0.4176 | 0.3509 |
| Exhaustive reference | coverage | 80% | 0.4535 | 0.6237 | 0.6394 | 0.6091 |
| Graph | intrinsic | 10% | 0.2787 | 0.4340 | 0.5991 | 0.3414 |
| Graph | intrinsic | 80% | 0.5651 | 0.7219 | 0.7595 | 0.6880 |
| Graph + prefilter | intrinsic | 10% | 0.2785 | 0.4337 | 0.5991 | 0.3411 |
| Graph + prefilter | intrinsic | 80% | 0.5652 | 0.7219 | 0.7595 | 0.6881 |
| Greedy | intrinsic | 10% | 0.2780 | 0.4332 | 0.5981 | 0.3407 |
| Greedy | intrinsic | 80% | 0.5650 | 0.7218 | 0.7594 | 0.6879 |
| Exhaustive reference (no merging) | intrinsic | 10% | 0.2886 | 0.4455 | 0.6430 | 0.3418 |
| Exhaustive reference (no merging) | intrinsic | 80% | 0.5989 | 0.7490 | 0.7798 | 0.7205 |
| Exhaustive reference | intrinsic | 10% | 0.2755 | 0.4301 | 0.5981 | 0.3368 |
| Exhaustive reference | intrinsic | 80% | 0.5600 | 0.7177 | 0.7550 | 0.6841 |

## Interpretation

### Candidate search

The seed recovers 97.0% of the pairs that pass both
exact thresholds while exactly scoring 486,814
distinct pairs, 0.97% of all possible
pairs. Sensitivity and cost are controlled by two independent mechanisms: the
geometry decides which residue columns can be matched, and the sound anchor
bound removes pairs that provably cannot be accepted.

The geometry ablation isolates the first. Indexing only the contiguous
column pairs of each terminal 3-mer, with everything else unchanged, drops
recall to 64.3%. Both contiguous pairs
contain the middle residue, so one substitution there destroys both, and the
required terminal columns shift whenever the peptides differ in length.

Every missed pair that remains (537 per dataset) failed a terminal
seed threshold, not the bound: across every run in this validation the sound
bound rejected 0 eligible pairs. The residual is a threshold trade-off,
shown in the sweep table, not a structural blind spot.

### Agreement

The exhaustive reference applies the identical clustering procedure to the
complete edge set, so a run differs from it only through candidate search.

With the `coverage` order, graph reaches ARI 0.8053 and pairwise
Jaccard 0.6746 against that reference, at
3894 clusters.

With the `intrinsic` order, graph reaches ARI 0.9480 and pairwise
Jaccard 0.9013 against that reference, at
4776 clusters.

### Stability

At the 80% subset with the `coverage` order, graph reaches pairwise Jaccard
0.4590, against 0.4535 for the exhaustive reference under
the identical comparison. The run therefore tracks what the same clustering
procedure does on a complete edge set, and the composition dependence that
remains belongs to the procedure rather than to the search.

That reference is a comparison point, not an upper bound, and a run can score
slightly above it: missing a few percent of edges makes clusters marginally
smaller, and smaller clusters have fewer co-cluster pairs to disagree about.
Stability must therefore be read next to the agreement table, never alone.

Disabling merging on the reference gives 0.4514, so the merge stage accounts
for +0.0022
of it. Reassignment is the larger term, which is why it carries the
hysteresis margin.

At the 80% subset with the `intrinsic` order, graph reaches pairwise Jaccard
0.5651, against 0.5600 for the exhaustive reference under
the identical comparison. The run therefore tracks what the same clustering
procedure does on a complete edge set, and the composition dependence that
remains belongs to the procedure rather than to the search.

That reference is a comparison point, not an upper bound, and a run can score
slightly above it: missing a few percent of edges makes clusters marginally
smaller, and smaller clusters have fewer co-cluster pairs to disagree about.
Stability must therefore be read next to the agreement table, never alone.

Disabling merging on the reference gives 0.5989, so the merge stage accounts
for -0.0389
of it. Reassignment is the larger term, which is why it carries the
hysteresis margin.

### Choosing a representative order

`coverage` minimises the cluster count (3894
clusters) but its selection key is a degree, so it amplifies small edge-set
differences: it reaches ARI 0.81 against its reference,
against 0.95 for `intrinsic`. `intrinsic` is also more
stable (0.5651 versus 0.4590 pairwise Jaccard at 80%), at
the cost of 4776 clusters instead of 3894. Neither is uniformly better: choose
`coverage` for compactness and `intrinsic` when reproducibility across dataset
revisions matters more.

### Static greedy

Static greedy scores only 72.0% of eligible pairs,
but the audit shows this is structural rather than a search defect: of
5037 missed pairs only 527 failed the seed, while
4509 were retrievable
but never examined, because the method only ever scores
representative-to-unassigned pairs. Its pair-level recall is therefore not
comparable with the graph method's, and lazy-exact greedy remains the
recommended low-memory path.

### Limitations

- The reference is computational, derived from the same scoring rule. It says
  nothing about biological cluster purity, which still requires labelled
  peptide-MHC data.
- The seed threshold is a sensitivity/cost trade-off, not a guarantee: at the
  0.40 default roughly 3% of eligible pairs are still not retrieved.
- The iterative section remains the largest source of composition dependence.
- All results are at 10,000 peptides. See `figures/scaling_benchmark.csv` for
  the cost profile at 1k-1M.

## Files

- `runs/exhaustive/`, `runs/exhaustive_subsets/`: exact pairs and both
  reference partitions for each representative order.
- `runs/full/`, `runs/subsets/`, `runs/sweep/`: clustering runs.
- `figures/`: plot-matched CSV, PNG and PDF.
- `metrics/missed_pair_audit.csv`: per-run attribution of missed pairs.
- `code/`: preparation, execution, exhaustive reference and analysis source.
