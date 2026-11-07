# Circuit-aware protection against burst errors in lattice surgery: reproducibility package

Data and code for checking the results in *Circuit-aware protection against burst errors in lattice surgery*
(E. Kim). One command recomputes every result that can be traced to per-event failure counts. It then compares each
one with the value printed in the paper.

```bash
pip install -r requirements.txt
python reproduce.py          # about 10 s; writes results/summary.md
```

`reproduce.py` does four things:
1. Checks the SHA-256 of every data file.
2. Recomputes three studies from per-event counts, with the original estimators and bootstrap seeds, and requires an
   exact match with the stored analyses.
3. Reads the stored analyses for the distance-9 mechanism studies.
4. Prints a paper-versus-data table.

It exits with an error only if a checksum or a recomputation fails. The GitHub Actions workflow runs it on every push.

## What can be checked

| Paper result | Paper | Recomputed here | Trace level |
|---|---|---|---|
| Protected ancillas at the shallow-reference risk (Table I, Fig. 2b) | 118.8 (117.0–133.1) | 118.8 (117.0–133.1) | per-event counts |
| Saving ΔN (abstract) | 49 (35–51) | 49 (35–51) | per-event counts |
| Table I typed excess risk, 117–150 sites (10⁻⁴) | 8.29, 6.19, 5.10, 3.07 | same | per-event counts |
| Estimated-burst decoder break-even burst probability | ≈ 4.3% per operation | 4.3% | per-event counts; see note 1 |
| d = 13, three response models: typed lower under both background decoders (SM S5) | yes | yes | per-event counts |
| d = 13, three response models: oracle comparisons (SM S5) | not resolved | Fc resolved for V1 and V10 | per-event counts; see note 2 |
| Report × ancilla-Pauli interaction (Fig. 4b) | 1.12 (0.71–1.60) × 10⁻⁴ | same | stored analysis |
| Path ranking of local swaps, Spearman ρ (SM S4) | 0.77 | 0.77 | stored analysis |
| Path cost W at depth one (Fig. 3c) | 7.8 / 10.4 / 22.7 | same | stored isolated-site scores |
| Deep/typed ratio η at d = 9, 3,200 impacts (Fig. 2a) | 1.62 (1.49–1.77) | same | original archived analysis |
| Typed reduction, 1 − 1/η (abstract) | 38% (33–44%) | same | original archived analysis |
| η for the other decoders (Fig. 2a) | 1.43, 1.36, 1.28 | same | original archived analysis |
| d = 13 oracle ratios, 1,200 impacts (Fig. 5a) | 0.96 (0.62–1.56), 1.12 (0.48–2.65) | same | original archived analysis |
| Five-round burst cap, shallow/typed (Fig. 5b) | 0.60 (0.28–0.97) | same | original archived analysis |

"Original archived analysis" means the analysis file from the 31 October 2026 publication package the figures were
built from, checked byte for byte against that package's own checksums (`data/p1_archive/`). Those cohorts' 15,600
per-event block files are listed in `data/p1_archive/PACKAGE_SHA256SUMS.txt` but not yet included. Until they are,
these rows confirm the archived analysis rather than recompute it.

### Notes

1. The 4.3% break-even rests on 1 harmful and 0 rescued burst-free shots out of 131,072. The penalty's 95% interval
   includes zero ([`data/r9_estimated_burst_decoder/README.md`](data/r9_estimated_burst_decoder/README.md)).
2. Under the correlated oracle, the deep-minus-typed difference is resolved for two response models, while the ratio
   intervals are unbounded ([`data/r1_d13_response_models/README.md`](data/r1_d13_response_models/README.md)).

## Contents

```
reproduce.py                       runs every check below
requirements.txt                   versions this package was verified with
analysis/
  r9_estimated_burst_decoder.py    R9 analysis (standalone copy of the original)
  r1_d13_response_models.py        R1 analysis (standalone copy of the original)
data/
  count_sweep/                     refined protection-count sweep, d = 9 (Table I, Fig. 2b)
    recompute_results.py           standalone recomputation, NumPy only
    RUN_RESULTS.md                 design, results and provenance
    *PROTOCOL*.json, RUNTIME_AMENDMENT.json, CONTROL_REPRODUCTION.json   pre-registration and audit trail
    analysis_source/               original simulation and analysis code for this run (for reference)
  r9_estimated_burst_decoder/      R9: 1,200 impacts, nine decoders
  r1_d13_response_models/          R1: 2,000 impacts at d = 13, three response models
  d9_mechanism_inputs/             stored analyses: factorial, local-swap exchange test, isolated-site path costs,
                                   initial count sweep (baselines/frontier.json), composition
  p1_archive/                      original Fig. 2a, Fig. 5a and Fig. 5b,c cohorts: analyses, event banks, quiet
                                   baselines, and the package's checksums
figures/                           scripts that draw Figs. 2–5 from figure_data.json and the count-sweep results
SHA256SUMS.txt                     checksums of every file under data/
```

## Data format

Each study stores failures per simulated impact (*event*), input mode (Z or X), layout and decoder as integer counts
out of a fixed number of shots. Burst-free blocks for the quiet baseline are stored separately. Excess risk is the
pool-weighted stratified mean failure fraction minus the quiet baseline, taking the worse input mode. Intervals come from 4,000
paired within-stratum event-bootstrap draws, with Dirichlet source-pool weights and Jeffreys quiet-baseline draws. Each
analysis script documents its estimator in its header.

## Figures

```bash
pip install -r figures/requirements.txt
cd figures
python make_manuscript_figures.py ../data/count_sweep/RESULTS.json   # Figs. 2-5
python figs_v11.py               ../data/count_sweep/RESULTS.json   # final Figs. 2 and 5
```

Output goes to `figures/out/`. The figures use URW Nimbus Sans if you place it in `figures/fonts/`
([urw-base35-fonts](https://github.com/ArtifexSoftware/urw-base35-fonts)); otherwise a system sans-serif is used.

## Simulation code

These files are simulation outputs. The circuit-level simulator, impact model and decoders are in
[Noro-Official/waveguide-please](https://github.com/Noro-Official/waveguide-please) at commit `9848047e`
(`scripts/pacer_postv2_campaign.py`, `scripts/pacer_postv2_r3/`, `requirements-pacer-postv2.txt`).
`data/count_sweep/analysis_source/` keeps the count sweep's own run and analysis scripts. They need that full
simulation tree; `recompute_results.py` does not.

## License and citation

Code: MIT ([LICENSE](LICENSE)). Data: CC BY 4.0 ([LICENSE-DATA](LICENSE-DATA)). Please cite the paper; citation
metadata are in [CITATION.cff](CITATION.cff).
