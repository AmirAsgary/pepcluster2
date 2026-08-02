#!/bin/bash
set -e
ROOT=/user/a.hajialiasgarynaj01/u14286/cbscratch/amirasgary2/pepcluster2
NEW=$ROOT/validation/2026-07-30_0.5.0-dev_search_redesign
python3 $NEW/code/margin_sweep.py --root "$NEW" --binary $ROOT/target/release/pepcluster2 "$@"
