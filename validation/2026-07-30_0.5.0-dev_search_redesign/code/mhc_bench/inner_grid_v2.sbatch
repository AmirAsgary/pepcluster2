#!/bin/bash
# Inner hyperparameter grid: 125 pools x 300 configurations = 37,500 runs.
#
# 16 tasks of 64 cores = 1024 cores. Few, large tasks rather than many small
# ones: a per-user cap on concurrent jobs makes a wide array unusable, and each
# task fills its own cores with a worker pool, which packs better than one run
# per job because pepcluster2 does not saturate 8 threads on a small pool.
#SBATCH --job-name=pc2_mhc2_inner
#SBATCH --partition=soeding
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/inner2_%A_%a.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
# Node-local /tmp is only ~4 GB; a 100k graph run spills a few hundred megabytes
# and 64 concurrent runs would overflow it.
TMP=/cbscratch/amirasgary2/pc2_mhc_tmp/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

echo "host $(hostname) task $SLURM_ARRAY_TASK_ID cores $SLURM_CPUS_PER_TASK started $(date -Is)"
python3 -u "$NEW/code/mhc_bench/run_grid.py" \
  --root "$NEW/runs/mhc_bench_v2" \
  --binary "$ROOT/target/release/pepcluster2" \
  --shard "$SLURM_ARRAY_TASK_ID" --shards 16 \
  --cores "$SLURM_CPUS_PER_TASK" --splits inner --tmp-root "$TMP" \
  --attempt "$SLURM_ARRAY_JOB_ID"
echo "finished $(date -Is)"
