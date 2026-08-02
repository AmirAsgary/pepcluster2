set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
OLD=$ROOT/validation/2026-07-29_0.4.3-dev_final_validation
HELPER=$NEW/code/exhaustive_reference/target/release/pepcluster2-exhaustive-reference
W=$(mktemp -d /tmp/pc2ref.XXXXXX)
zcat $NEW/data/full/sample_000.fasta.gz > $W/in.fasta
mkdir -p $NEW/runs/exhaustive/sample_000
$HELPER $W/in.fasta $NEW/runs/exhaustive/sample_000 32
cat $NEW/runs/exhaustive/sample_000/run_stats.json
echo "--- ground truth identical to 0.4.3? ---"
if cmp -s $NEW/runs/exhaustive/sample_000/ground_truth_clusters.tsv $OLD/runs/exhaustive/sample_000/ground_truth_clusters.tsv; then
  echo "IDENTICAL"
else
  echo "DIFFERS"; diff <(head -20 $NEW/runs/exhaustive/sample_000/ground_truth_clusters.tsv) <(head -20 $OLD/runs/exhaustive/sample_000/ground_truth_clusters.tsv) | head
fi
echo "--- true_pairs identical? ---"
cmp -s $NEW/runs/exhaustive/sample_000/true_pairs.bin $OLD/runs/exhaustive/sample_000/true_pairs.bin && echo "IDENTICAL" || echo "DIFFERS"
rm -rf $W
