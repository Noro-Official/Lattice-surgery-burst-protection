# R1: response models at distance 13

The paper's SM S5 cohort: 2,000 impacts at d = 13, three stochastic response models, typed and deep layouts,
four decoders.

| File | Contents |
|---|---|
| `events.npz` | Event bank: 2,000 impacts in five strata, source-pool counts, impact positions and onsets |
| `rows.parquet` | Failure counts per event × input mode × layout × response model (V0_baseline, V1_phase_unsuppressed, V10_phase_echo_measured), 1,024 shots per cell, for D0, D0c, F, Fc |
| `quiet.parquet` | Burst-free blocks (4,096 per mode, 1,024 shots each) |
| `RESULTS.json` | Stored analysis (key `R1`). `src/lsburst/r1_response_models.py` recomputes it exactly (step 2). |
| `FILES.sha256` | Checksums from the original release |

The verdict for each response model and decoder uses the deep-minus-typed excess-risk difference ΔE. It is "superior" when the 95% interval of ΔE lies above zero.

| Decoder | V0 | V1 | V10 |
|---|---|---|---|
| D0 (background) | superior | superior | superior |
| D0c (background, correlated) | superior | superior | superior |
| F (oracle) | unresolved | unresolved | unresolved |
| Fc (oracle, correlated) | unresolved | **superior** | **superior** |

**Differs from the paper.** The paper says the oracle comparisons are not resolved. By the ΔE criterion, Fc favours typed
for V1 and V10. Under Fc the excess risks are close to zero, so the deep/typed *ratio* intervals are very wide and
include negative values. That may be what "not resolved" was meant to describe; if so, the text should say the ratio
is unresolved.

This cohort was run as a pilot (see the release notes in Noro-Official/waveguide-please,
`data/reports/pacer_postv2_release_20261105/README.md`). The event bank was extracted from that release's archive and
its SHA-256 checked against the release manifest.
