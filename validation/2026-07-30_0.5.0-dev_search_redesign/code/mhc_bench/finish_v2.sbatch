#!/bin/bash
# Everything after the inner grid, chained so the study completes unattended:
# select on the inner folds, evaluate the selection on the held-out alleles and
# on the independent test pools, then build every figure and CSV.
#SBATCH --job-name=pc2_mhc2_finish
#SBATCH --partition=soeding
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/finish2_%j.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
BENCH=$NEW/runs/mhc_bench_v2
TMP=/cbscratch/amirasgary2/pc2_mhc_tmp/finish_$SLURM_JOB_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

echo "=== verifying the inner grid is complete before using it ==="
python3 "$NEW/code/mhc_bench/verify.py" --root "$BENCH" \
  --binary "$ROOT/target/release/pepcluster2" --tag inner --splits inner \
  --sample 12 --tmp-root "$TMP"

echo "=== selecting hyperparameters on the inner folds ==="
python3 "$NEW/code/mhc_bench/analyse.py" --root "$BENCH"

rm -f "$BENCH"/grid/outer_*.csv "$BENCH"/grid/test_*.csv
echo "=== outer folds: held-out alleles, per-fold selection ==="
python3 "$NEW/code/mhc_bench/run_grid.py" --root "$BENCH" \
  --binary "$ROOT/target/release/pepcluster2" --shard 0 --shards 1 --cores 64 --attempt "$SLURM_JOB_ID" \
  --splits outer --configs "$BENCH/tables/selected_per_fold.csv" --tmp-root "$TMP"

echo "=== independent test pools ==="
python3 "$NEW/code/mhc_bench/run_grid.py" --root "$BENCH" \
  --binary "$ROOT/target/release/pepcluster2" --shard 0 --shards 1 --cores 64 --attempt "$SLURM_JOB_ID" \
  --splits test --configs "$BENCH/tables/selected_overall.csv" --tmp-root "$TMP"

echo "=== verifying the evaluation runs ==="
python3 "$NEW/code/mhc_bench/verify.py" --root "$BENCH" \
  --binary "$ROOT/target/release/pepcluster2" --tag test --splits test --sample 0 \
  --tmp-root "$TMP" || echo "note: test grid holds only the selected configurations"

echo "=== final tables and figures ==="
python3 "$NEW/code/mhc_bench/analyse.py" --root "$BENCH"
echo "finished $(date -Is)"
