#!/bin/bash
# Serial timing of the three representative-selection rules on an EXCLUSIVE node.
# Accuracy is measured separately and in parallel; timing cannot be, because
# concurrent jobs contaminate wall clocks.
#SBATCH --job-name=pc2_seltime
#SBATCH --partition=soeding
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/seltime_%j.out
set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
D=$ROOT/validation/2026-07-29_0.4.3-dev_final_validation/data/benchmark
V=$ROOT/validation/2026-08-06_0.6.0_motif_merge
TMP=${TMPDIR:-/local/jobs/${USER}_${SLURM_JOB_ID}}/seltime
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT
OUT=$V/results/selection_timing.csv
echo "peptides,selection,threshold,wall_seconds,clusters,pairs_scored,index_hits" > "$OUT"
for N in 0010000 0050000 0100000 0500000; do
  gunzip -c "$D/benchmark_$N.fasta.gz" > "$TMP/$N.fasta"
  for SEL in kmer-degree lazy-exact; do
    for T in 0.40 0.25; do
      rm -rf "$TMP/out"
      # lazy-exact is superlinear; cap so one cell cannot eat the whole job
      S=$( { /usr/bin/time -f "%e" timeout 7200 "$ROOT/target/release/pepcluster2" \
          -i "$TMP/$N.fasta" -o "$TMP/out" \
          --mode separate_kmer_anchor --kmer-similarity-threshold "$T" \
          --anchor-combination-similarity-threshold 0.35 \
          --clustering-method greedy --greedy-selection "$SEL" \
          --compact-output --threads 1 --tmp-dir "$TMP/tmp" >/dev/null; } 2>&1 | tail -1 )
      if [ -f "$TMP/out/run_stats.json" ]; then
        python3 - "$TMP/out/run_stats.json" "$N" "$SEL" "$T" "$S" "$OUT" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
with open(sys.argv[6],'a') as h:
    h.write(f"{s['unique_sequences']},{sys.argv[3]},{sys.argv[4]},{sys.argv[5]},"
            f"{s['final_clusters']},{s['candidate_pairs_computed']},{s['index_candidate_hits']}\n")
PY
      else
        echo "$N,$SEL,$T,TIMEOUT,,," >> "$OUT"
      fi
      echo "done $N $SEL $T"
    done
  done
done
