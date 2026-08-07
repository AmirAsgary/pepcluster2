#!/bin/bash
# Runtime comparison. Requires an EXCLUSIVE node: the first attempt at this ran
# on a shared login node while a GibbsCluster sweep was using ~48 cores, and one
# pool came out 6x slower than a re-run of the same pool. Serial timing is only
# meaningful with nothing else competing.
#SBATCH --job-name=pc2_speed
#SBATCH --partition=soeding
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --output=logs/speed_%j.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
STUDY=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
V=$ROOT/validation/2026-08-06_0.6.0_motif_merge
TMP=${TMPDIR:-/local/jobs/${USER}_${SLURM_JOB_ID}}/speed
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
export PATH="$STUDY/benchmark/tools/env/bin:$PATH"

echo "host $(hostname) exclusive; started $(date -Is)"
python3 -u "$V/code/speed_benchmark.py" \
  --binary "$ROOT/target/release/pepcluster2" \
  --mixmhcp "$STUDY/benchmark/tools/MixMHCp/MixMHCp" \
  --gibbscluster "$STUDY/benchmark/tools/gibbscluster-2.0/gibbscluster" \
  --pools "$STUDY/runs/mhc_bench_sep_kmer_anchor/pools" \
  --manifest "$STUDY/runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv" \
  --selected "$V/results/motif_selected_overall.csv" \
  --tmp-root "$TMP" --out "$V/results/speed_benchmark.csv" --n-pools 10
echo "finished $(date -Is)"
