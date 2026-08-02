# Where to read the results

Three reports, in the order they are worth reading.

| Read | File | Answers |
|---|---|---|
| 1 | `runs/mhc_bench_comparison/REPORT.md` | Which scoring mode wins, at which settings, and how it scores on the benchmark. Start here. |
| 2 | `runs/mhc_bench_sep_kmer_anchor/REPORT.md` | The k-mer mode on its own: full threshold behaviour, per-fold selection stability, anchor ablation. |
| 3 | `runs/mhc_bench_sep_aln_anchor/REPORT.md` | The same for the alignment mode. |

`REPORT.md` at this level covers the earlier search-redesign work (candidate
generation, seed geometry, stability) and is independent of the MHC benchmark.

## Figures

All in `runs/mhc_bench_comparison/`, as `.png` and `.pdf`:

- `hyperparameter_effect` - how far each threshold moves AMI, purity and cluster count.
- `overall_performance` - both modes at their selected setting, on all three splits.
- `benchmark_kmer`, `benchmark_aln` - the selected setting on the test pools.

Per-mode figures (threshold sweeps, the two-threshold interaction surface, pool
sizes) live in each mode's `figures/` directory.

## Numbers behind the figures

Every figure has a CSV of the exact plotted values next to it. Raw per-run
results are in each mode's `tables/*_grid_raw.csv`, and the untouched grid output
in `grid/`.

## Reproducing

```
code/mhc_bench/run_grid.py    # the sweep (sharded, run under SLURM)
code/mhc_bench/verify.py      # completeness gate - run before trusting any grid
code/mhc_bench/analyse.py     # per-mode selection, figures, REPORT.md
code/mhc_bench/compare.py     # cross-mode figures and REPORT.md
```

`verify.py` is not optional: a shard that dies still exits zero, so completeness
is asserted against the expected job list rather than against exit status.

## Not reported

PepCluster, the earlier method, is excluded entirely: its results and figures
have been deleted. The scripts that produced them (`code/mhc_bench/
run_pepcluster.py` and `pepcluster.sbatch`) are kept, so the run can be
reproduced if it is ever wanted again.
