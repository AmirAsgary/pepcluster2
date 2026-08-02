#!/bin/bash
# Graph paths, 1k-1M, both seed configurations.
#
# These nodes are always shared, so an exclusive allocation would never
# schedule. Slurm pins this job to its own cores, so the threads are not
# stolen, but memory bandwidth is shared and absolute seconds drift with
# background load. The runner therefore executes sizes sequentially, runs the
# two seed configurations back to back, and records node load with each
# measurement, so the legacy-versus-redesign comparison stays valid even
# though the absolute numbers are not perfectly reproducible.
#SBATCH --job-name=pc2_scale_graph
#SBATCH --partition=soeding
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
TMP=/cbscratch/amirasgary2/pc2_scaling_tmp_$SLURM_JOB_ID

echo "host: $(hostname); cpus: $(nproc); started $(date -Is)"
python3 "$NEW/code/run_scaling_benchmark.py" \
  --binary "$ROOT/target/release/pepcluster2" \
  --data-dir "$NEW/data/benchmark" \
  --output-dir "$NEW/scaling" \
  --tmp-root "$TMP" \
  --methods graph graph_prefilter \
  --threads 48 --candidate-buffer-mb 8192 --max-memory-gb 400 \
  --retry-failed
rm -rf "$TMP"
echo "finished $(date -Is)"
