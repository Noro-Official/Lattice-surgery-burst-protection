#!/usr/bin/env bash
# Step 2: recompute every study from its per-event counts and check it against the stored analysis.
# Writes out/results/<study>.json. Takes about 10 s.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

python tools/recompute.py --out_dir out/results
