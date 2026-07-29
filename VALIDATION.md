# PepCluster2 0.4.3 validation

The final computational validation is stored in:

```text
validation/2026-07-29_0.4.3-dev_final_validation/
```

Despite the historical directory suffix, all final clustering runs used the
release-identical scoring and clustering implementation now labelled `0.4.3`.
The release change only removed the development version suffix and added
packaging metadata.

The validation used 20 independently sampled datasets of 10,000 peptides,
alignment-similarity threshold 0.50, anchor-combination-similarity threshold
0.60, and terminal/core weights 4:1. It completed:

- exhaustive scoring of every pair in all 20 datasets;
- 80 full-data runs across graph, forced-prefilter graph, static greedy, and
  lazy-exact greedy;
- 400 nested-subset runs at 1k, 2k, 4k, 6k, and 8k peptides.

Main findings:

- graph and lazy-exact greedy scored about 1.63% of all possible pairs and
  recovered about 64.5% of all exactly eligible pairs;
- static greedy recovered about 54.0% of eligible pairs;
- graph and forced-prefilter graph were effectively indistinguishable at this
  tested scale;
- graph/lazy-exact agreement with exhaustive greedy set cover had mean ARI
  about 0.30; static greedy had mean ARI about 0.18;
- at an 80% subset, graph/lazy-exact stability had ARI about 0.70 and pairwise
  Jaccard about 0.54.

Read `REPORT.md` in the validation directory for definitions, full results,
interpretation, and conclusions. Plot-matched CSV, PNG, and PDF files are in
`figures/`; exact source is in `code/`.

The resource benchmark data are already generated under `data/benchmark/`.
Run the cluster-only benchmark using `CLUSTER_RESOURCE_BENCHMARK.md` and
`code/run_resource_benchmark.py`.
