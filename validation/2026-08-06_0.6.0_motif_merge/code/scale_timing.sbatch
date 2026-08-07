#!/bin/bash
# Scale timing with the OPTIMISED binary at full node parallelism, which is how
# a large run would actually be executed. The earlier table
# (selection_timing_v070_serial.csv) used the pre-optimisation binary at one
# thread and understates what is now achievable.
#SBATCH --job-name=pc2_scale
#SBATCH --partition=soeding
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=logs/scale_%j.out
set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
D=$ROOT/validation/2026-07-29_0.4.3-dev_final_validation/data/benchmark
V=$ROOT/validation/2026-08-06_0.6.0_motif_merge
TMP=${TMPDIR:-/local/jobs/${USER}_${SLURM_JOB_ID}}/scale
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
OUT=$V/results/scale_timing_optimised.csv
echo "peptides,selection,threshold,threads,wall_seconds,peak_rss_mb,clusters,pairs_scored" > "$OUT"
for N in 0100000 0500000 1000000; do
  gunzip -c "$D/benchmark_$N.fasta.gz" > "$TMP/$N.fasta"
  for SEL in lazy-exact kmer-degree; do
    for T in 0.40 0.25; do
      rm -rf "$TMP/out"
      R=$( { /usr/bin/time -f "%e %M" timeout 14400 "$ROOT/target/release/pepcluster2" \
          -i "$TMP/$N.fasta" -o "$TMP/out" \
          --mode separate_kmer_anchor --kmer-similarity-threshold "$T" \
          --anchor-combination-similarity-threshold 0.35 \
          --clustering-method greedy --greedy-selection "$SEL" \
          --compact-output --threads 64 --tmp-dir "$TMP/tmp" >/dev/null; } 2>&1 | tail -1 )
      if [ -f "$TMP/out/run_stats.json" ]; then
        python3 - "$TMP/out/run_stats.json" "$SEL" "$T" "$R" "$OUT" <<'PY'
import json,sys
s=json.load(open(sys.argv[1])); w,m=sys.argv[4].split()
open(sys.argv[5],'a').write(f"{s['unique_sequences']},{sys.argv[2]},{sys.argv[3]},64,"
    f"{w},{float(m)/1024:.0f},{s['final_clusters']},{s['candidate_pairs_computed']}\n")
PY
      else
        echo "$N,$SEL,$T,64,TIMEOUT,,," >> "$OUT"
      fi
      echo "done $N $SEL $T"
    done
  done
done
