set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
W=$(mktemp -d /tmp/pc2v.XXXXXX)
zcat $NEW/data/full/sample_000.fasta.gz > $W/in.fasta
$ROOT/target/release/pepcluster2 --input $W/in.fasta --output-dir $W/out \
  --mode separate_aln_anchor --alignment-similarity-threshold 0.50 \
  --anchor-combination-similarity-threshold 0.60 --kmer-seed-threshold 0.40 \
  --terminal-seed all-column-pairs --representative-order coverage \
  --gap-open -4 --gap-extension -1 --terminal-overhang-gap-open -2 \
  --terminal-overhang-gap-extension -1 --minimum-terminal-match-length 2 \
  --threads 16 --candidate-buffer-mb 512 --compact-output \
  --clustering-method graph --no-prefilter --tmp-dir $W/tmp > /dev/null 2>&1
if cmp -s $W/out/node_clusters.tsv $NEW/runs/full/coverage/graph/sample_000/node_clusters.tsv; then
  echo "PASS: post-fix rerun is byte-identical to the benchmarked partition"
else
  echo "DIFFERS"; diff <(head -5 $W/out/node_clusters.tsv) <(head -5 $NEW/runs/full/coverage/graph/sample_000/node_clusters.tsv)
fi
rm -rf $W
