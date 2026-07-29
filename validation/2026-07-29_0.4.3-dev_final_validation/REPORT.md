# PepCluster2 0.4.3 final validation

## Configuration

- 20 independently sampled datasets; 10,000 peptides each.
- Scoring mode: separate_aln_anchor.
- Alignment-similarity threshold: 0.50.
- Anchor-combination-similarity threshold: 0.60.
- Terminal/core alignment weights: 4/1.
- Ground truth: every pair scored exactly, followed by dynamic greedy set cover.

The exhaustive partition is a computational reference under the chosen scoring rule, not biological ground truth.

## How to read the metrics

- Search recall is the fraction of all truly eligible peptide pairs that the method actually scored.
- Search precision is the fraction of scored pairs that were truly eligible; low precision means extra work, not incorrect final edges.
- ARI is the Adjusted Rand Index: agreement of two clusterings after correcting for chance (1 is identical).
- NMI is Normalized Mutual Information: shared cluster information (1 is identical). With many small clusters, NMI can remain high even when assignments differ, so ARI and pairwise Jaccard are more informative here.

## Search-rule performance

| Method | Recall mean ± SD | Precision mean ± SD | All-pairs fraction mean ± SD |
|---|---:|---:|---:|
| Graph | 0.6428 ± 0.0039 | 0.0142 ± 0.0002 | 0.0162 ± 0.0002 |
| Graph + prefilter | 0.6428 ± 0.0039 | 0.0142 ± 0.0002 | 0.0162 ± 0.0002 |
| Greedy | 0.5398 ± 0.0037 | 0.0136 ± 0.0002 | 0.0143 ± 0.0002 |
| Greedy lazy-exact | 0.6454 ± 0.0039 | 0.0143 ± 0.0002 | 0.0163 ± 0.0002 |

## Ground-truth cluster agreement

| Method | ARI mean ± SD | NMI mean ± SD | Clusters mean | Singletons mean |
|---|---:|---:|---:|---:|
| Graph | 0.2999 ± 0.0069 | 0.9132 ± 0.0013 | 4837.9 | 2630.1 |
| Graph + prefilter | 0.2999 ± 0.0069 | 0.9132 ± 0.0013 | 4837.8 | 2630.1 |
| Greedy | 0.1773 ± 0.0064 | 0.9047 ± 0.0017 | 5377.1 | 2877.6 |
| Greedy lazy-exact | 0.3008 ± 0.0073 | 0.9133 ± 0.0013 | 4834.1 | 2630.2 |

Exhaustive ground truth contained a mean of 2313.0 singleton clusters (23.13% of 10,000 peptides).

## Stability

Stability compares a clustering of a subset with the full-dataset clustering after restricting the latter to the same peptides. Higher values mean less dependence on dataset size.

| Method | Subset | Jaccard | ARI | NMI | Pair recall | Pair precision |
|---|---:|---:|---:|---:|---:|---:|
| Graph | 10% | 0.3243 | 0.4884 | 0.9885 | 0.5267 | 0.4581 |
| Graph | 80% | 0.5391 | 0.7002 | 0.9732 | 0.7148 | 0.6863 |
| Graph + prefilter | 10% | 0.3244 | 0.4885 | 0.9885 | 0.5270 | 0.4581 |
| Graph + prefilter | 80% | 0.5391 | 0.7002 | 0.9732 | 0.7148 | 0.6863 |
| Greedy | 10% | 0.3517 | 0.5184 | 0.9912 | 0.6430 | 0.4365 |
| Greedy | 80% | 0.5713 | 0.7268 | 0.9807 | 0.7569 | 0.6991 |
| Greedy lazy-exact | 10% | 0.3240 | 0.4881 | 0.9885 | 0.5260 | 0.4580 |
| Greedy lazy-exact | 80% | 0.5385 | 0.6996 | 0.9731 | 0.7146 | 0.6855 |

## Final interpretation

PepCluster2 avoids almost all all-pairs work: graph scored about 1.62% of all possible pairs. However, its search rule recovered only 64.3% of pairs that pass both exact thresholds. Static greedy recovered still fewer (54.0%). Therefore the current candidate search is fast but not exhaustive.

Graph and lazy-exact greedy were the closest to the exhaustive set-cover reference, but their mean ARI was only about 0.30. Static greedy was worse (ARI 0.18) and produced about 540 more clusters per dataset. The high NMI values should not be read as near-identity because these data contain thousands of small and singleton clusters.

Forced graph prefiltering was effectively indistinguishable from non-prefilter graph at 10,000 peptides in both search and clustering metrics. This supports using the prefilter at this tested scale, but does not prove equivalence on larger or biologically different datasets.

Cluster stability improved steadily as the subset approached the full dataset. At 80%, ARI was about 0.70 and pairwise Jaccard about 0.54 for graph/lazy-exact (0.57 for static greedy). Thus clusters are reasonably, but not fully, stable to dataset composition.

Overall conclusion: use graph when its temporary edge storage is affordable, or lazy-exact greedy when memory/disk is limiting. Do not claim that the present search recovers every eligible relationship or that clusters are dataset-independent. Biological validation against peptide–MHC labels remains necessary before making purity claims.

## Files

- runs/exhaustive/: exact true pairs and exhaustive set-cover assignments.
- runs/full/: four full clustering paths and unique exact-scoring traces.
- runs/subsets/: all nested subset clusterings.
- figures/: plot-matched CSV, PNG, and PDF files.
- code/: complete preparation, execution, exhaustive-scoring, and analysis source.
