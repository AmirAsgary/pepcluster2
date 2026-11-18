#!/bin/bash
# Time and memory for the largest pool, one mode, before committing the grid.
#SBATCH --job-name=jd_probe
#SBATCH --partition=soeding
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2/validation/2026-11-18_0.10.1_julianedata_motif_merge/logs/probe_%j.out

set -euo pipefail
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
HERE=$ROOT/validation/2026-11-18_0.10.1_julianedata_motif_merge
TMP=${TMPDIR:-/tmp}/jdprobe/$SLURM_JOB_ID
mkdir -p "$TMP"; trap 'rm -rf "$TMP"' EXIT

for POOL in a91_s100000_eq_r0 a24_s032000_eq_r0; do
  N=$(head -1 <(grep -c '^>' "$HERE/pools/$POOL.fasta") || true)
  echo "=== $POOL ($(grep -c '^>' $HERE/pools/$POOL.fasta) peptides) ==="
  for FLAGS in "" "--merge-motifs"; do
    /usr/bin/time -f "  mode[${FLAGS:-clustering}] wall %e s  maxRSS %M KB" \
      "$ROOT/target/release/pepcluster2" \
      --input "$HERE/pools/$POOL.fasta" --output-dir "$TMP/out" \
      --mode separate_kmer_anchor --kmer-similarity-threshold 0.25 \
      --anchor-combination-similarity-threshold 0.35 \
      --representative-order coverage --clustering-method graph \
      --no-prefilter --compact-output --threads 1 \
      --tmp-dir "$TMP/t" --max-memory-gb 100 $FLAGS \
      --motif-prior-concentration 10 --motif-merge-threshold 25 \
      --motif-em-prior-concentration 3 2>&1 | tail -2
    rm -rf "$TMP/out" "$TMP/t"
  done
done
