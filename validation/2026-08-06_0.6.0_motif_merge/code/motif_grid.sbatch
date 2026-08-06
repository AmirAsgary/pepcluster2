#!/bin/bash
# Motif-layer hyperparameter grid: 193 pools x 60 configurations = 11,580 runs.
#
# Same shape as the 0.5.0 inner grid: few large tasks rather than a wide array,
# because a per-user cap on concurrent jobs makes a wide array unusable, and each
# task fills its own cores with a worker pool. Each pepcluster2 run here is
# single-threaded, so the pool packs one run per core.
#SBATCH --job-name=pc2_motif_grid
#SBATCH --partition=soeding
#SBATCH --array=0-7%8
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/motif_grid_%A_%a.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
STUDY=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
V=$ROOT/validation/2026-08-06_0.6.0_motif_merge

# The graph path spills ~560 MB for the largest pool here, so 64 concurrent runs
# need ~36 GB of scratch. The first attempt put that on shared BeeGFS and every
# shard died with "Remote I/O error (os error 121)". The nodes have 184 GB of
# local NVMe and SLURM already points $TMPDIR at it; /tmp is only 4 GB.
TMP=${TMPDIR:-/local/jobs/$USER_$SLURM_JOB_ID}/pc2_motif
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$V/grid"

echo "host $(hostname) task $SLURM_ARRAY_TASK_ID cores $SLURM_CPUS_PER_TASK started $(date -Is)"
python3 -u "$V/code/run_motif_grid.py" \
  --binary "$ROOT/target/release/pepcluster2" \
  --pools "$STUDY/runs/mhc_bench_sep_kmer_anchor/pools" \
  --manifest "$STUDY/runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv" \
  --out "$V/grid/motif_shard_$(printf '%04d' "$SLURM_ARRAY_TASK_ID")_${SLURM_ARRAY_JOB_ID}.csv" \
  --tmp-root "$TMP" \
  --shard "$SLURM_ARRAY_TASK_ID" --shards 8 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --splits inner outer test
echo "finished $(date -Is)"
