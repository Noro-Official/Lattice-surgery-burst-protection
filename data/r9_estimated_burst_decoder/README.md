# R9: estimated-burst and other realistic decoders

Used in the paper for the estimated-burst decoder's break-even burst probability (about 4.3% per operation).

| File | Contents |
|---|---|
| `events.npz` | Event bank: 1,200 impacts, strata (100/200/900) and source-pool counts |
| `rows.parquet` | Failure counts per event × input mode (Z, X) × layout (typed, deep) × response variant (V0, V10), 1,024 shots per cell, for nine decoders (D0, D0c, F, Fc, Ff, E, L1, L3, K) plus estimator diagnostics |
| `quiet.parquet` | Burst-free blocks (128 per mode, 1,024 shots each), with background failures and the false-alarm harms and rescues of decoders E and K |
| `RESULTS.json` | Stored analysis. `python analysis/r9_estimated_burst_decoder.py` recomputes it exactly. |

The break-even value is `quiet_penalty["V0_baseline|typed|E|X"].break_even_burst_prob_per_op`
= dq / (R + dq). Here R = 1.69 × 10⁻⁴ is the burst-risk reduction of E relative to D0, and dq = 7.6 × 10⁻⁶ is the
burst-free penalty.

**Caution.** dq comes from 1 harmful and 0 rescued burst-free shots out of 131,072. Its 95% interval
(−1.1 × 10⁻⁵ to 3.4 × 10⁻⁵) includes zero, so the 4.3% figure is a point estimate. The data are also consistent with no
false-alarm penalty at all.

Source: `data/reports/pacer_round3/r9/` in Noro-Official/waveguide-please. That folder also holds the raw part
files, the calibration library and the completion receipt (`PUBLICATION.json`).
