#!/bin/bash
# Verify, select, evaluate and report one scoring-mode analysis.
#SBATCH --job-name=pc2_aln_finish
#SBATCH --partition=soeding
#SBATCH --cpus-per-task=64
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=logs/aln_finish_%j.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
BENCH=$NEW/runs/mhc_bench_sep_aln_anchor
BIN=$ROOT/target/release/pepcluster2
TMP=/cbscratch/amirasgary2/pc2_mhc_tmp/finish_$SLURM_JOB_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

echo "=== gate: the inner grid must be complete before anything is read from it ==="
python3 "$NEW/code/mhc_bench/verify.py" --root "$BENCH" --binary "$BIN" \
  --tag inner --splits inner --mode separate_aln_anchor --sample 12 --tmp-root "$TMP"

echo "=== select on the inner folds ==="
python3 "$NEW/code/mhc_bench/analyse.py" --root "$BENCH" --mode separate_aln_anchor

rm -f "$BENCH"/grid/outer_*.csv "$BENCH"/grid/test_*.csv
echo "=== held-out alleles ==="
python3 "$NEW/code/mhc_bench/run_grid.py" --root "$BENCH" --binary "$BIN" --mode separate_aln_anchor \
  --shard 0 --shards 1 --cores 64 --splits outer \
  --configs "$BENCH/tables/selected_per_fold.csv" --tmp-root "$TMP" --attempt "$SLURM_JOB_ID"

echo "=== independent test pools ==="
python3 "$NEW/code/mhc_bench/run_grid.py" --root "$BENCH" --binary "$BIN" --mode separate_aln_anchor \
  --shard 0 --shards 1 --cores 64 --splits test \
  --configs "$BENCH/tables/selected_overall.csv" --tmp-root "$TMP" --attempt "$SLURM_JOB_ID"

echo "=== final tables, figures and REPORT.md ==="
python3 "$NEW/code/mhc_bench/analyse.py" --root "$BENCH" --mode separate_aln_anchor
echo "finished $(date -Is)"
