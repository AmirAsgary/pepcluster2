# Benchmark inputs

`discovery_immuneapp.xlsx` and `test_immuneapp.xlsx` are ImmuneApp supplementary
tables. They are not redistributed here — they are third-party data — so the two
files are gitignored and must be placed in this directory before the peptide-MHC
benchmark can be rebuilt.

Everything downstream of them is reproducible from the tracked scripts:

```
code/mhc_bench/prepare_balanced.py   # pools, splits, similarity filtering
code/mhc_bench/run_grid.py           # the sweep
code/mhc_bench/verify.py             # completeness gate
code/mhc_bench/analyse.py            # per-mode selection and figures
code/mhc_bench/compare.py            # cross-mode figures and report
```

`hla.fasta` and `hla_similarity.csv` are tracked, since they are derived here
rather than taken from another publication.
