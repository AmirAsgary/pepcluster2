#!/bin/bash
# GibbsCluster on the benchmark pools: 193 pools x 2 settings = 386 runs.
#
# GibbsCluster is far more expensive than the other tools here. The documented
# MHC-I preset scans 1-6 groups with five seeds each, so the `default` setting is
# 21 clusterings x 5 seeds = 105 MCMC fits per pool, and the pools run to 24,888
# peptides. A 972-peptide pool did not finish in 35 minutes on one core.
#
# There are only 386 runs, so throwing more concurrent slots at it does nothing -
# per-run latency is the binding constraint. Each run therefore gets -k 16 (the
# tool forks over seed x group-count, 30 fits for the default preset and 5 for
# forced k) and only 4 run concurrently per task, filling the 64 cores.
# 16 tasks x 4 = 64 concurrent runs over 386 jobs, ~6 rounds.
#
# Per-run timeout 8 h: anything slower is recorded as a failure with that reason
# rather than silently missing, and the coverage columns keep a partial result
# distinguishable from a complete one.
#SBATCH --job-name=pc2_gibbs
#SBATCH --partition=soeding
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gibbs_%A_%a.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
STUDY=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
BENCH=$STUDY/benchmark
V=$ROOT/validation/2026-08-06_0.6.0_motif_merge

# Node-local NVMe, not the shared mount: GibbsCluster writes a directory tree per
# run and 64 concurrent runs on BeeGFS produced I/O errors for the motif grid.
TMP=${TMPDIR:-/local/jobs/${USER}_${SLURM_JOB_ID}}/gibbs
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$BENCH/results/immuneapp/raw/gibbs_shards"

export PATH="$BENCH/tools/env/bin:$PATH"
echo "host $(hostname) task $SLURM_ARRAY_TASK_ID cores $SLURM_CPUS_PER_TASK started $(date -Is)"
python3 -u "$BENCH/code/run_gibbscluster.py" \
  --pools "$STUDY/runs/mhc_bench_sep_kmer_anchor/pools" \
  --manifest "$STUDY/runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv" \
  --binary "$BENCH/tools/gibbscluster-2.0/gibbscluster" \
  --out "$BENCH/results/immuneapp/raw/gibbs_shards/gibbs_$(printf '%04d' "$SLURM_ARRAY_TASK_ID").csv" \
  --tmp-root "$TMP" \
  --shard "$SLURM_ARRAY_TASK_ID" --shards 16 \
  --cores 4 --gibbs-threads 16 \
  --timeout 28800 \
  --splits inner outer test || true
echo "finished $(date -Is)"
