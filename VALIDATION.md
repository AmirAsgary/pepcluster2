# PepCluster2 validation

`validation/2026-07-30_0.5.0-dev_search_redesign/REPORT.md` holds the current
results, definitions, figures and limitations. It covers exhaustive references
for 20 independently sampled datasets of 10,000 peptides and all 100 nested
subsets, a seed geometry and threshold sweep, 140 full-data runs and 700 subset
runs across four clustering paths and two representative orders, at
alignment-similarity threshold 0.50, anchor-combination-similarity threshold
0.60, and terminal/core alignment weights 4:1.

## The reference

Every pair is scored exactly, which removes candidate search from the
comparison, and the complete edge set then feeds the identical clustering
procedure the tool runs. A run therefore differs from the reference only through
candidate search, which is what makes the agreement number interpretable. A
reference that stopped earlier in the procedure would not be reachable by a run
that completes it, and agreement against such a reference would confound search
loss with the effect of the later stages.

The same reference is emitted for both representative orders, so a run is always
compared against a reference using its own selection rule. A variant with merging
disabled is also emitted, which attributes the merge stage.

## Headline results (graph, 20 datasets)

| | coverage order | intrinsic order |
|---|---:|---:|
| Search recall | 0.9701 | 0.9701 |
| Search precision | 0.0358 | 0.0358 |
| Distinct pairs exactly scored | 486,814 | 486,814 |
| Fraction of all pairs scored | 0.97% | 0.97% |
| ARI vs exhaustive reference | 0.8053 | **0.9480** |
| Pairwise Jaccard, 80% subset | 0.4590 | 0.5651 |
| Clusters | 3,894 | 4,776 |

Graph, forced-prefilter graph and lazy-exact greedy agree with each other to
within 0.001 on every metric. Static greedy is the outlier at ARI 0.4361, which
the missed-pair audit attributes to pairs it never examines rather than to the
seed.

## What the numbers do and do not say

- Search recall is a threshold trade-off, not a structural blind spot: every
  missed pair failed a terminal seed threshold, the sound anchor bound rejected
  none, and the sweep shows recall 0.9865 at seed 0.35 and 0.9917 at 0.30.
- Stability tracks what the same procedure achieves on a complete edge set. At
  the 80% subset the graph run reaches pairwise Jaccard 0.4590 against 0.4535 for
  the reference under the identical comparison. That reference is a comparison
  point, not an upper bound: missing a few percent of edges yields slightly
  smaller clusters, which have fewer co-cluster pairs to disagree about, so a run
  can score marginally above it. Stability must be read next to the agreement
  table, never alone.
- Reassignment, not merging, dominates composition dependence, which is why
  reassignment carries the `--reassignment-margin` hysteresis.
- The references are computational, derived from the same scoring rule. They say
  nothing about biological cluster purity, which still requires labelled
  peptide-MHC data.

## Reproducing

```bash
module load rust/1.85.0
cargo build --release
cd validation/2026-07-30_0.5.0-dev_search_redesign
(cd code/exhaustive_reference && cargo build --release)
sbatch code/full_validation.sbatch      # references, sweep, full, subsets, analysis
sbatch code/scaling_graph.sbatch        # 1k-1M speed, memory and temporary disk
sbatch code/scaling_greedy.sbatch
python3 code/margin_sweep.py --root . --binary ../../target/release/pepcluster2
```

An earlier validation round is archived under
`validation/2026-07-29_0.4.3-dev_final_validation/` together with the shared
benchmark datasets in its `data/` directory, which the current round reuses
byte-for-byte.

All thresholds and gap parameters remain subject to biological calibration on
labelled peptide-MHC data.
