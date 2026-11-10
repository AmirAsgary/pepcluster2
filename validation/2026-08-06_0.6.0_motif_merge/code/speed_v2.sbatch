#!/bin/bash
# Speed benchmark, rebuilt. Two differences from the first version:
#   * all seven PepCluster2 modes, with flags taken from run_variant.py so a
#     mode cannot be timed under settings it was not scored under;
#   * --repeats, because the first table was one run per pool and sub-second
#     timings made a requested count look cheaper than clustering alone.
# Task 0 is the fast arms (PepCluster2 + MixMHCp); task 1 is GibbsCluster,
# which dominates the cost and is reported separately until it lands.
# --exclusive holds the whole node so nothing else contends, but on this
# partition it allocates a single CPU unless --cpus-per-task is also given, and
# GibbsCluster forks at -k 16. Without the core request its sixteen processes
# share one core and every GibbsCluster number in the table is wrong.
#SBATCH --job-name=pc2_speed2
#SBATCH --partition=soeding
#SBATCH --array=0-1
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-08-06_0.6.0_motif_merge/logs/speed2_%A_%a.out

set -euo pipefail
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
HERE=$ROOT/validation/2026-08-06_0.6.0_motif_merge
STUDY=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
ARMS=(fast gibbscluster)
A=${ARMS[$SLURM_ARRAY_TASK_ID]}
TMP=${TMPDIR:-/tmp}/pc2speed/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

echo "host $(hostname) arms $A cores $(nproc) started $(date -Is)"
python3 -u "$HERE/code/speed_benchmark.py" \
  --binary "$ROOT/target/release/pepcluster2" \
  --mixmhcp "$STUDY/benchmark/tools/MixMHCp/MixMHCp" \
  --gibbscluster "$STUDY/benchmark/tools/gibbscluster-2.0/gibbscluster" \
  --pools "$STUDY/runs/mhc_bench_sep_kmer_anchor/pools" \
  --manifest "$STUDY/runs/mhc_bench_sep_kmer_anchor/pool_manifest.csv" \
  --selected "$HERE/results/motif_selected_overall.csv" \
  --tmp-root "$TMP" --n-pools 9 --repeats 3 --arms "$A" \
  --out "$HERE/results/speed_benchmark_${A}.csv"
echo "finished $(date -Is)"
