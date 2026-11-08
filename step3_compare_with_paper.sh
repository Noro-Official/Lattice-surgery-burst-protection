#!/usr/bin/env bash
# Step 3: compare the recomputed results with the numbers printed in the paper.
# Prints a table and writes it to out/summary.md.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

python tools/compare_with_paper.py --results_dir out/results --summary out/summary.md
