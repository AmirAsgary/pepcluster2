#!/bin/bash
# MixMHCp, default and given-K, exactly as the previous benchmark ran it.
# The given-K arm is capped by MixMHCp's own 50-motif ceiling, which is why the
# largest pools in this study hold 50 alleles and not all 91.
#SBATCH --job-name=jd_mixmhcp
#SBATCH --partition=soeding
#SBATCH --array=0-19%20
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-11-18_0.10.1_julianedata_motif_merge/logs/mixmhcp_%A_%a.out

set -euo pipefail
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
HERE=$ROOT/validation/2026-11-18_0.10.1_julianedata_motif_merge
BENCH=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign/benchmark
TMP=${TMPDIR:-/tmp}/jdmix/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

S=$SLURM_ARRAY_TASK_ID
echo "host $(hostname) shard $S/20 started $(date -Is)"
python3 -u "$BENCH/code/run_mixmhcp.py" \
  --binary "$BENCH/tools/MixMHCp/MixMHCp" \
  --pools "$HERE/pools" \
  --manifest "$HERE/pool_manifest.csv" \
  --out "$HERE/results/shards/mixmhcp_s${S}.csv" \
  --splits test --settings default oracle_k \
  --cores "$SLURM_CPUS_PER_TASK" --tmp-root "$TMP" \
  --shard "$S" --shards 20
echo "finished $(date -Is)"
