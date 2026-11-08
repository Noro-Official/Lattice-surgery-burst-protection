# Data

Simulation outputs for "Circuit-aware protection against burst errors in lattice surgery". Every file is listed in
[`../SHA256SUMS.txt`](../SHA256SUMS.txt) and checked by `step1_verify_data.sh`.

## Overview

| Folder | Paper use | Impacts | Distance | Decoders | Shots per cell |
|---|---|---|---|---|---|
| [`count_sweep/`](count_sweep/RUN_RESULTS.md) | Table I, Fig. 2b; ΔN = 49 | 600 | 9 | D0C | 1,024 |
| [`p1_archive/precision_d9/`](p1_archive/README.md) | Fig. 2a; η = 1.62, the 38% reduction | 3,200 | 9 | D0, D0C, F, FC | 4,096 |
| [`p1_archive/precision_d13/`](p1_archive/README.md) | Fig. 5a | 1,200 | 13 | D0, D0C, F, FC | 4,096 |
| [`p1_archive/regime/`](p1_archive/README.md) | Fig. 5b,c (burst-duration caps) | 600 | 9 | D0, F | 1,024 |
| [`r9_estimated_burst_decoder/`](r9_estimated_burst_decoder/README.md) | Estimated-burst decoder break-even (≈ 4.3%) | 1,200 | 9 | D0, D0c, F, Fc, Ff, E, L1, L3, K | 1,024 |
| [`r1_d13_response_models/`](r1_d13_response_models/README.md) | SM S5, three response models | 2,000 | 13 | D0, D0c, F, Fc | 1,024 |
| [`d9_mechanism_inputs/`](d9_mechanism_inputs/README.md) | Fig. 3c, Fig. 4b, SM S1, SM S4 | — | 9 | — | stored analyses only |
| [`calibration_extension/`](calibration_extension/RUN_RESULTS.md) | Not in the paper (v1.1.0): η vs impact energy, decoupling, trapping | 600 | 9 | D0C | 1,024 |

Decoder labels: **D0** background matching; **D0C**/**D0c** background matching with correlations; **F** event-informed
oracle; **FC**/**Fc** oracle with correlations. The R9-only decoders (**Ff**, **E**, **L1**, **L3**, **K**) are
described in [`r9_estimated_burst_decoder/README.md`](r9_estimated_burst_decoder/README.md).

Layout labels: **typed** (`typed`, `typed_paper`) protects the check type measured during each merge; **deep** and
**shallow** are geometric bands with matched protected counts and X:Z splits.

## Shared conventions

All studies share the following conventions.

- **Event.** One simulated impact, drawn from a stratified pool, with a position, onset and footprint. Each study has
  an event bank (`.npz`) that holds `stratum` (per event) and `pool_counts` (per stratum). The stratum weights are
  `pool_counts / pool_counts.sum()`.
- **Cell.** One event × input mode (`Z` or `X`) × layout × decoder (and response model, where present). Each cell
  holds integer failure counts out of a fixed number of shots.
- **Burst-free baseline ("quiet").** Failures in shots with no impact, stored separately for each input mode and
  base decoder. F uses D0's baseline and Fc/FC uses D0c/D0C's.
- **Excess risk.** The stratum-weighted mean per-event failure fraction minus the burst-free baseline, taking the
  worse input mode.
- **Intervals.** 95% percentile intervals from 4,000 paired within-stratum event-bootstrap replicates, with Dirichlet(pool counts + 1)
  stratum weights and Jeffreys Beta(k + ½, n − k + ½) baseline draws. Each analysis module in `src/lsburst/` states its
  seeds.

## File formats

### `count_sweep/`

| File | Format | Contents |
|---|---|---|
| `per_event_counts.csv.gz` | gzip CSV, 8,400 rows | `condition` (e.g. `n132\|typed`, `n168\|shallow`), `input_mode`, `event_index`, `stratum`, `shots`, `any_output_failures`, `histogram_00` … `histogram_11` (shots by the two logical outputs' correctness), `onset_round` |
| `events_d9_seed78001_n20000_300-180-120_fp1.npz` | NumPy archive | `xy` (impact position), `energy_keV`, `u_onset`, `dose`, `stratum`, `pool_index`, `pool_counts`, `footprints`, `x0_fp1` (per-site initial excitation, 600 × 521) |
| `quiet_d9.json` | JSON | `"{T_M}\|{mode}\|{decoder}": {fail, shots, p, se}` |
| `study.json`, `config.json` | JSON | study definition: counts swept, seeds, allocation, shots |
| `RESULTS.json` | JSON | stored analysis; `decoders.D0C` holds `typed_excess`, `n_star`, `n_star_ci`, `saving`, `saving_ci` |
| `sampled_risks.csv` | CSV | per-condition risks as reported |
| `FROZEN_PROTOCOL.json`, `ORIGINAL_FROZEN_PROTOCOL.json`, `ANALYSIS_PLAN.json`, `RUNTIME_AMENDMENT.json`, `CONTROL_REPRODUCTION.json`, `EXECUTION_RESUME.json`, `EXPORT_VERIFICATION.json` | JSON | pre-registration, amendment and audit trail |
| `recompute_results.py` | Python | the export's own standalone check (NumPy only) |
| `analysis_source/` | Python | the sweep's run and analysis scripts; they need the full simulation tree |
| `SHA256SUMS.json` | JSON | the export's checksums |

### `calibration_extension/`

Same formats as `count_sweep/` (same event bank and quiet file), with these differences:

| File | Format | Contents |
|---|---|---|
| `per_event_counts.csv.gz` | gzip CSV, 46,800 rows | as in `count_sweep/`; `condition` = `<layout>\|E<energy scale>\|<decoupling>\|T<trap multiplier>\|G<protected gap, GHz>` with layout `shallow`, `deep` or `typed_paper` |
| `RESULTS.json` | JSON | one row per physics setting: per-layout excess risks and intervals, `eta_deep_over_typed`, `eta_ci`, `shallow_over_typed` |
| `ANCHOR.json` | JSON | pre-registered regression of the ×1 standard setting against `count_sweep/` |
| `STEP5_CHECKS.json`, `STEP5_TIMING.json` | JSON | pre-freeze controls: bit-identity of default channels, peak-shift table, runtime, pooled ×100 saturation |
| `FROZEN_PROTOCOL.json`, `study_calext_*.json` | JSON | pre-registration and study definitions |
| `recompute_results.py` | Python | standalone check (NumPy only) |
| `eta_vs_energy.png`, `.pdf` | figure | η against energy scale |
| `analysis_source/` | Python | run, analysis and export scripts; they need the full simulation tree |

### `p1_archive/`

| File | Format | Contents |
|---|---|---|
| `<cohort>/blocks.jsonl.gz` | gzip JSON Lines | one line per original block file: `{"file", "sha256", "content"}`. `content` is the block's JSON text: `rows` of `{cond, e, mode, shots, onset, hist: {decoder: [n00, n01, n10, n11]}}`, where `n00` counts shots with both logical outputs correct. |
| `<cohort>/study.json` | JSON | conditions (`name`, `mask_name`, `tm`, and for `regime` the varied parameter), seeds, allocation, shots, decoders |
| `<cohort>/analysis.json`, `precision_analysis.json` | JSON | archived analyses; the d = 13 analysis is only in `precision_analysis.json` |
| `events/*.npz` | NumPy archive | event banks, with the same fields as above plus `score_disagree` and `score_typed` for the precision cohorts and `x0_fp{s}` per footprint scale for `regime` |
| `quiet_d9.json`, `quiet_d9_p2.json`, `quiet_d13_p2.json` | JSON | burst-free baselines, `"{T_M}\|{mode}\|{decoder}": {fail, shots, p, se}` |
| `PACKAGE_SHA256SUMS.txt` | text | the original package's checksums for every file of these cohorts, including the 15,600 block files |

### `r9_estimated_burst_decoder/` and `r1_d13_response_models/`

| File | Format | Contents |
|---|---|---|
| `rows.parquet` | Parquet | one row per cell. R9: `mode`, `event`, `stratum`, `layout`, `variant`, `shots`, `fail_<decoder>`, `harm_<decoder>`, `rescue_<decoder>` (harms and rescues are relative to D0), and estimator diagnostics (`E_accepted`, `E_distinct`, `E_onset_abs_err_median`, `E_pos_err_mm_median`, `K_alarms`, `F_eq_Ff`). R1: `exp`, `mode`, `event`, `variant`, `mask`, `shots`, `seed`, `fail_D0`, `fail_D0c`, `fail_F`, `fail_Fc`, `stratum`. |
| `quiet.parquet` | Parquet | burst-free blocks: `mode`, `block`, `shots`, `fail_D0`, `fail_D0c`; R9 adds `{layout}\|{variant}\|{fail,harm,rescue}_{E,K}` and `{layout}\|{E,K}_alarms` |
| `events.npz` | NumPy archive | `x0` (per-site initial excitation), `onset`, `stratum`, `pool_counts`, `pool_index`, `xy`, `heavy_threshold` |
| `RESULTS.json` | JSON | stored analysis |

### `d9_mechanism_inputs/`

These folders hold stored analyses only; their per-event counts are not part of this repository.

| File | Paper use | Key |
|---|---|---|
| `factorial/analysis.json` | Fig. 4b, report × ancilla-Pauli interaction | `effects["shallow - typed_paper\|D0"].RA` |
| `exchange_test/analysis.json` | SM S4, path ranking of local swaps | `D0.spearman_calc` |
| `site_scores_d9_tm9.json` | Fig. 3c, isolated-site path costs W | `W["full\|D0"]`, `W["reports\|D0"]`, `quiet` (burst-free W) |
| `baselines/frontier.json` | SM S1, initial count sweep (132, 168, 204 sites) | — |
| `composition/sealed_analysis.json` | composition comparisons | — |
