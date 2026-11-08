#!/usr/bin/env bash
# Step 4: draw Figs. 2-5 into out/figures/ (needs step 2 and tools/figures/requirements.txt).
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"

python tools/build_figure_inputs.py --results_dir out/results --out out/figure_data.json
python tools/figures/make_manuscript_figures.py data/count_sweep/RESULTS.json   # Figs. 2-5
python tools/figures/figs_v11.py data/count_sweep/RESULTS.json                  # final Figs. 2 and 5
