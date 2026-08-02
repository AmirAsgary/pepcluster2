#!/bin/bash
# Single-dataset smoke comparison of the historical and redesigned terminal seed.
set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
BIN=$ROOT/target/release/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
WORK=$(mktemp -d /tmp/pc2smoke.XXXXXX)
zcat "$NEW/data/full/sample_000.fasta.gz" > "$WORK/in.fasta"

common="--input $WORK/in.fasta --mode separate_aln_anchor \
  --alignment-similarity-threshold 0.50 \
  --anchor-combination-similarity-threshold 0.60 \
  --gap-open -4 --gap-extension -1 \
  --terminal-overhang-gap-open -2 --terminal-overhang-gap-extension -1 \
  --minimum-terminal-match-length 2 --threads ${THREADS:-32} \
  --candidate-buffer-mb 512 --compact-output --write-scored-pairs \
  --no-prefilter --clustering-method graph"

for tag in legacy new050 new040 new035 new030; do
  case $tag in
    legacy) extra="--terminal-seed contiguous --kmer-seed-threshold 0.50";;
    new050) extra="--terminal-seed all-column-pairs --kmer-seed-threshold 0.50";;
    new040) extra="--terminal-seed all-column-pairs --kmer-seed-threshold 0.40";;
    new035) extra="--terminal-seed all-column-pairs --kmer-seed-threshold 0.35";;
    new030) extra="--terminal-seed all-column-pairs --kmer-seed-threshold 0.30";;
  esac
  out=$NEW/runs/smoke/$tag
  mkdir -p "$out"
  /usr/bin/time -v -o "$out/resource.txt" $BIN $common $extra \
    --output-dir "$out" --tmp-dir "$WORK/tmp_$tag" > "$out/run.log" 2>&1
  echo "=== $tag ==="
  grep -E "Index candidate|Rejected by|Candidate pairs scored|Constrained-alignment|Eligible graph edges|^Clusters|Singleton clusters|Elapsed seconds" "$out/run_summary.txt"
  grep -E "Maximum resident" "$out/resource.txt"
done
rm -rf "$WORK"
