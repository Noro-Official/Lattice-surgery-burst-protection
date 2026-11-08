# Code and data for "Circuit-aware protection against burst errors in lattice surgery"

[![reproduce](https://github.com/Noro-Official/Lattice-surgery-burst-protection/actions/workflows/reproduce.yml/badge.svg)](https://github.com/Noro-Official/Lattice-surgery-burst-protection/actions/workflows/reproduce.yml)
[![Code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24046121.svg)](https://doi.org/10.5281/zenodo.24046121)

Eric Kim

This repository contains the simulation data and analysis code behind every number and data figure in the paper.
It recomputes each result from per-event failure counts, checks the result against the stored analysis, and compares
it with the value printed in the paper.

### Abstract

Large-scale quantum computing will require logical operations protected against high-energy impacts that disturb many
superconducting qubits at once. Gap engineering suppresses relaxation induced by these impacts, but protecting only a
subset of qubits requires choosing where that protection is placed. During lattice surgery, measurements change as
patches merge and split, exposing different ancillas and temporal boundaries to burst errors. Here we compare three
static protection layouts for a simulated lattice-surgery CNOT, matching the numbers of protected ancillas and the X:Z
check-type splits. At distance 9, protecting the check type measured during each merge reduces burst-induced excess
risk by 38 (+6/−5)% (95% bootstrap interval) relative to a deep geometric band under correlated background matching,
without giving the decoder the impact parameters. We attribute this difference to report errors at temporal
boundaries and ancilla faults propagating to data, which interact to increase the excess risk of a shallow band. We
also find that the typed layout reaches the shallow-band reference risk with an estimated 49 (+2/−14) fewer protected
ancillas. At distance 13, the advantage is not resolved under oracle decoding, and a five-round burst cap reverses
the shallow-band ordering. These results indicate how the measurements of a logical circuit can guide protection
placement within the simulated decoder and burst-duration regimes tested here.

## Requirements

Python 3.13 (3.11 or later should work) with:

- numpy==2.5.3
- pandas==3.0.5
- scipy==1.18.1
- pyarrow==26.0.1

Drawing the figures (step 4) also needs the packages in [`tools/figures/requirements.txt`](tools/figures/requirements.txt).
The workflow was tested on Linux (Debian and Ubuntu).

## Reproducing the results

Step 0: create an environment and install the dependencies.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Step 1: verify the checksums of every data file, including the archived cohorts against their original package.

```bash
./step1_verify_data.sh
```

Step 2: recompute every study from its per-event counts. Each result must match the stored analysis exactly
(relative tolerance 1e-9); the script exits with an error otherwise. Results are written to `out/results/`. This
takes about 10 seconds.

```bash
./step2_recompute_results.sh
```

Step 3: compare the recomputed results with the numbers printed in the paper. The table is printed and written to
`out/summary.md`.

```bash
./step3_compare_with_paper.sh
```

Step 4 (optional): draw Figs. 2–5 into `out/figures/`.

```bash
pip install -r tools/figures/requirements.txt
./step4_plot_figures.sh
```

The [`reproduce`](.github/workflows/reproduce.yml) workflow runs all four steps on every push.

## Results checked against the paper

Step 3 prints the following table.

| Paper item | Paper | Recomputed | Trace level |
|---|---|---|---|
| Protected ancillas at the shallow-reference risk (Table I, Fig. 2b) | 118.8 (117.0–133.1) | 118.8 (117.0–133.1) | per-event counts |
| Protected-ancilla saving ΔN (abstract) | 49 (35–51) | 49 (35–51) | per-event counts |
| Typed excess risk at 117, 123, 132, 150 sites (Table I, 10⁻⁴) | 8.29, 6.19, 5.10, 3.07 | same | per-event counts |
| Deep/typed ratio η, correlated background matching, d = 9 (Fig. 2a) | 1.62 (1.49–1.77) | 1.62 (1.49–1.77) | per-event counts |
| Typed reduction in excess risk, 1 − 1/η (abstract) | 38% (33–44%) | 38% (33–44%) | per-event counts |
| η for FC, D0 and F (Fig. 2a) | 1.43, 1.36, 1.28 | same | per-event counts |
| Oracle deep/typed ratios at d = 13, F and FC (Fig. 5a) | 0.96 (0.62–1.56), 1.12 (0.48–2.65) | same | per-event counts |
| Five-round cap, shallow/typed, background matching (Fig. 5b) | 0.60 (0.28–0.97) | 0.60 (0.28–0.97) | per-event counts |
| Estimated-burst decoder break-even burst probability | about 4.3% per operation | about 4.3% | per-event counts |
| d = 13, three response models: typed lower under both background decoders (SM S5) | yes, all three | yes, all three | per-event counts |
| d = 13, three response models: oracle comparisons (SM S5) | not resolved | resolved for Fc under V1 and V10 | per-event counts |
| Report × ancilla-Pauli interaction (Fig. 4b, 10⁻⁴) | 1.12 (0.71–1.60) | 1.12 (0.71–1.60) | stored analysis |
| Path ranking of local swaps, Spearman ρ (SM S4) | 0.77 | 0.77 | stored analysis |
| Path cost W at depth one (Fig. 3c) | 7.8 / 10.4 / 22.7 | same | stored analysis |

"Per-event counts" means the value is recomputed from the integer failure counts of every simulated impact.
"Stored analysis" means the value is read from an analysis file whose per-event counts are not included.

## Calibration extension (version 1.1.0, not in the paper)

[`data/calibration_extension/`](data/calibration_extension/RUN_RESULTS.md) reruns the Fig. 2a layouts on the count
sweep's 600-impact bank. It varies three model choices the paper held fixed: impact energy (×1 to ×100, where ×100 is
the scale at which the phonon generation was calibrated to Willow's 2 MHz median peak shift), standard versus echoed
decoupling, and faster trapping or a 10 GHz gap on protected sites. The deep/typed ratio η is resolved above 1 only for
small impacts: 1.43 (1.11–1.90) at the paper's scale and 1.20 (1.09–1.31) at ×3. From ×10 it is indistinguishable from
1, including for Willow-calibrated impacts with an echoed circuit: 0.97 (0.95–1.00). Decoupling and the trap and gap
variants do not change this. Check it with `python data/calibration_extension/recompute_results.py`; the workflow
runs this as an extra step.

## Known differences and caveats

- **SM S5 oracle comparisons.** The paper says these comparisons are not resolved. Under the correlated oracle (Fc),
  the deep-minus-typed difference has a 95% interval that excludes zero for response models V1 and V10. The deep/typed ratio
  intervals are unbounded there because the Fc excess risk is close to zero. Step 3 reports this row as differing.
- **Break-even burst probability.** The 4.3% value rests on 1 harmful and 0 rescued burst-free shots out of 131,072,
  and the false-alarm penalty's 95% interval includes zero. The paper presents it as a point estimate.
- **Figures 3 and 4.** Their plotted values were recovered from the manuscript's vector figures (about 0.2% of each
  axis span) and are kept in [`assets/figure_data.json`](assets/figure_data.json). Figs. 2 and 5 are drawn from
  recomputed results. Fig. 1 (layout masks) is not generated here.
- **Simulator.** This repository contains simulation outputs, not the simulator. The circuit-level simulator, impact
  model and decoders are in [Noro-Official/waveguide-please](https://github.com/Noro-Official/waveguide-please)
  (commit `9848047e`; `scripts/pacer_postv2_campaign.py`, `scripts/pacer_postv2_r3/`). The count sweep's own run
  scripts are in [`data/count_sweep/analysis_source/`](data/count_sweep/analysis_source/).

## Repository layout

- [`step1_verify_data.sh`](step1_verify_data.sh) … [`step4_plot_figures.sh`](step4_plot_figures.sh): the
  reproduction steps above.
- [`data/`](data/): simulation outputs, one folder per study; formats in [`data/README.md`](data/README.md).
    - [`data/count_sweep/`](data/count_sweep/): refined protection-count sweep (Table I, Fig. 2b).
    - [`data/p1_archive/`](data/p1_archive/): original cohorts for Fig. 2a, Fig. 5a and Fig. 5b,c.
    - [`data/r9_estimated_burst_decoder/`](data/r9_estimated_burst_decoder/): realistic decoders, break-even value.
    - [`data/r1_d13_response_models/`](data/r1_d13_response_models/): three response models at d = 13 (SM S5).
    - [`data/d9_mechanism_inputs/`](data/d9_mechanism_inputs/): stored analyses for Figs. 3c and 4b and SM S1, S4.
    - [`data/calibration_extension/`](data/calibration_extension/RUN_RESULTS.md): energy-scale, decoupling and trapping
      extension of Fig. 2a (version 1.1.0; not in the paper).
- [`src/lsburst/`](src/lsburst/): Python package; one module per study, each with `analyze(data_dir) -> dict`.
    - [`count_sweep.py`](src/lsburst/count_sweep.py), [`p1_archive.py`](src/lsburst/p1_archive.py),
      [`r9_estimated_decoder.py`](src/lsburst/r9_estimated_decoder.py),
      [`r1_response_models.py`](src/lsburst/r1_response_models.py): analyses.
    - [`paper.py`](src/lsburst/paper.py): the paper's printed values and the checks against them.
- [`tools/`](tools/): command-line scripts called by the step scripts.
    - [`tools/figures/`](tools/figures/): figure scripts and style (`make_manuscript_figures.py`, `figs_v11.py`).
- [`assets/figure_data.json`](assets/figure_data.json): plotted values for Figs. 3 and 4, and the starting point
  that step 4 updates with recomputed values.
- `out/`: ignored by git; the step scripts write their output here.
- [`SHA256SUMS.txt`](SHA256SUMS.txt): checksums of every file under `data/`.

## Citation

If you use these data or this code, please cite the paper and this archive:

```bibtex
@misc{kim2026burst,
  title  = {Circuit-aware protection against burst errors in lattice surgery},
  author = {Kim, Eric},
  year   = {2026},
  note   = {Manuscript in preparation}
}

@dataset{kim2026burstdata,
  title     = {Code and data for "Circuit-aware protection against burst errors in lattice surgery"},
  author    = {Kim, Eric},
  year      = {2026},
  publisher = {Zenodo},
  version   = {1.1.0},
  doi       = {10.5281/zenodo.24046121},
  url       = {https://doi.org/10.5281/zenodo.24046121}
}
```

The archive is on Zenodo: [10.5281/zenodo.24046121](https://doi.org/10.5281/zenodo.24046121) always resolves to the latest
version; [10.5281/zenodo.24046122](https://doi.org/10.5281/zenodo.24046122) is version 1.0.0. Citation metadata are also
in [`CITATION.cff`](CITATION.cff). The paper entry will be updated with its journal reference.

## License

Code is released under the [MIT License](LICENSE). The data in `data/` are released under
[CC BY 4.0](LICENSE-DATA).
