#!/bin/bash
# The seven motif modes on an independent mono-allelic MS benchmark.
# 160 pools are split round-robin over 20 array tasks, at most 20 at a time;
# each task runs all seven modes over its own eighth of the pools.
#SBATCH --job-name=jd_variants
#SBATCH --partition=soeding
#SBATCH --array=0-19%20
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-11-18_0.10.1_julianedata_motif_merge/logs/variants_%A_%a.out

set -euo pipefail
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
HERE=$ROOT/validation/2026-11-18_0.10.1_julianedata_motif_merge
SEL=$ROOT/validation/2026-08-06_0.6.0_motif_merge
# $TMPDIR is node-local NVMe. The scratch filesystem must not be used for
# --tmp-dir: the edge store writes continuously and BeeGFS cannot keep up.
TMP=${TMPDIR:-/tmp}/jdvar/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

S=$SLURM_ARRAY_TASK_ID
echo "host $(hostname) shard $S/20 cores $SLURM_CPUS_PER_TASK started $(date -Is)"

for MODE in similarity merge merge_k merge_em merge_em_k em_only em_only_k; do
  echo "--- $MODE (shard $S) $(date -Is) ---"
  python3 -u "$SEL/code/run_variant.py" \
    --binary "$ROOT/target/release/pepcluster2" \
    --pools "$HERE/pools" \
    --manifest "$HERE/pool_manifest.csv" \
    --selected "$SEL/results/motif_selected_overall.csv" \
    --out "$HERE/results/shards/variant_${MODE}_s${S}.csv" \
    --tmp-root "$TMP" --variant "$MODE" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --max-memory-gb 18 \
    --shard "$S" --shards 20
done
echo "finished $(date -Is)"
