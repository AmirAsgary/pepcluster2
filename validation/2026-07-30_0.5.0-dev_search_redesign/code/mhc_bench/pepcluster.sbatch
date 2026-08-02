#!/bin/bash
#SBATCH --job-name=pc2_pepcluster
#SBATCH --partition=soeding
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --output=logs/pepcluster_%j.out
set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
BENCH=$NEW/runs/mhc_bench_pepcluster
TMP=/cbscratch/amirasgary2/pc2_mhc_tmp/pepcluster_$SLURM_JOB_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
for split in inner outer test; do
  echo "=== $split ==="
  python3 -u "$NEW/code/mhc_bench/run_pepcluster.py" --root "$BENCH" \
    --splits $split --cores 64 --tmp-root "$TMP"
done
echo "finished $(date -Is)"
