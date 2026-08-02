#!/bin/bash
# separate_aln_anchor: 120 pools x 330 configurations.
#
# The grid is 10 primary thresholds x 11 anchor thresholds (0.35-0.80 plus zero,
# which disables the anchor condition and isolates the primary component) x
# 2 representative orders for graph and 1 for lazy-exact greedy.
#SBATCH --job-name=pc2_aln_inner
#SBATCH --partition=soeding
#SBATCH --array=0-15%16
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=logs/aln_inner_%A_%a.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
BENCH=$NEW/runs/mhc_bench_sep_aln_anchor
TMP=/cbscratch/amirasgary2/pc2_mhc_tmp/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

echo "host $(hostname) task $SLURM_ARRAY_TASK_ID cores $SLURM_CPUS_PER_TASK started $(date -Is)"
python3 -u "$NEW/code/mhc_bench/run_grid.py" \
  --root "$BENCH" --binary "$ROOT/target/release/pepcluster2" \
  --mode separate_aln_anchor \
  --shard "$SLURM_ARRAY_TASK_ID" --shards 16 \
  --cores "$SLURM_CPUS_PER_TASK" --splits inner --tmp-root "$TMP" \
  --todo "$BENCH/todo_inner.csv" --attempt "$SLURM_ARRAY_JOB_ID"
echo "finished $(date -Is)"
