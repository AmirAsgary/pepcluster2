#!/bin/bash
#SBATCH --job-name=pc2_0500dev
#SBATCH --partition=soeding
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out

set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign

echo "host: $(hostname); cpus: $(nproc); started $(date -Is)"
python3 "$NEW/code/run_validation.py" \
  --root "$NEW" \
  --binary "$ROOT/target/release/pepcluster2" \
  --helper "$NEW/code/exhaustive_reference/target/release/pepcluster2-exhaustive-reference" \
  --datasets 20 --sweep-datasets 5 \
  --workers 12 --threads 8 \
  --stage all
echo "finished $(date -Is)"
