#!/bin/bash
# The seven reported motif modes plus the threshold-merged given-K variant,
# one array task per mode, every one of the 193 benchmark pools.
#SBATCH --job-name=pc2_variants
#SBATCH --partition=soeding
#SBATCH --array=0-7%8
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-08-06_0.6.0_motif_merge/logs/variant_%A_%a.out

set -euo pipefail
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
HERE=$ROOT/validation/2026-08-06_0.6.0_motif_merge
STUDY=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
V=(similarity merge merge_k merge_em merge_em_k em_only em_only_k merge_em_k_merged)
NAME=${V[$SLURM_ARRAY_TASK_ID]}
TMP=${TMPDIR:-/tmp}/pc2var/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

echo "host $(hostname) mode $NAME started $(date -Is)"
python3 -u "$HERE/code/run_variant.py" \
  --binary "$ROOT/target/release/pepcluster2" \
  --pools "$STUDY/runs/mhc_bench_sep_kmer_anchor/pools" \
  --manifest "$STUDY/runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv" \
  --selected "$HERE/results/motif_selected_overall.csv" \
  --out "$HERE/results/variant_${NAME}.csv" \
  --tmp-root "$TMP" --variant "$NAME" --workers "$SLURM_CPUS_PER_TASK"
echo "finished $(date -Is)"
