# Cluster resource benchmark

The seven nested benchmark FASTA files are already present in `data/benchmark/`
and their sampling provenance is recorded in `data/benchmark/manifest.csv`.
The runner records wall time with Python and peak resident memory with
`/usr/bin/time -v`. It is resumable and writes one status file per run.

Build PepCluster2 on the compute node, then run:

```bash
cargo build --release --locked
python code/run_resource_benchmark.py \
  --binary ../../target/release/pepcluster2 \
  --data-dir data/benchmark \
  --output-dir benchmark \
  --tmp-root /path/to/large/scratch/pepcluster2 \
  --threads 16
```

Use a scratch filesystem with at least 200 GiB free if the non-prefilter graph
is run at one million peptides. The feasibility estimate for that case is
approximately 149 GiB before safety margin. To split the benchmark across jobs,
pass subsets such as:

```bash
python code/run_resource_benchmark.py ... \
  --methods graph_prefilter greedy greedy_lazy \
  --sizes 500000 1000000
```

After all jobs finish, regenerate the combined CSV and figure without running
anything again:

```bash
python code/run_resource_benchmark.py ... --plot-only
```

Outputs are `benchmark/figures/resource_benchmark.csv`, PNG, and PDF, with
complete logs and exact commands retained below `benchmark/runs/`.
