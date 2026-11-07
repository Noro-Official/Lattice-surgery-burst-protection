# Data

| Folder | Paper use | Level |
|---|---|---|
| [`count_sweep/`](count_sweep/RUN_RESULTS.md) | Table I, Fig. 2b, abstract (118.8 protected ancillas, ΔN = 49) | per-event counts |
| [`r9_estimated_burst_decoder/`](r9_estimated_burst_decoder/README.md) | Estimated-burst decoder break-even (≈ 4.3%) | per-event counts |
| [`r1_d13_response_models/`](r1_d13_response_models/README.md) | SM S5: three response models at d = 13 | per-event counts |
| [`d9_mechanism_inputs/`](d9_mechanism_inputs/README.md) | Fig. 3 path costs, Fig. 4 factorial, SM S4 local swaps, SM S1 initial count sweep | stored analyses |
| [`p1_archive/`](p1_archive/README.md) | Fig. 2a (η = 1.62), Fig. 5a (d = 13), Fig. 5b,c (burst caps) | per-event counts |

In `d9_mechanism_inputs/`:
- `factorial/analysis.json` holds the Fig. 4 channel factorial. The paper's interaction is `effects["shallow - typed_paper|D0"].RA`.
- `exchange_test/analysis.json` holds the SM S4 local swaps, `D0.spearman_calc`.
- `site_scores_d9_tm9.json` holds the Fig. 3 isolated-site path costs W per channel set; `quiet` gives burst-free W.
- `baselines/frontier.json` holds the SM S1 initial count sweep at 132, 168 and 204 sites.
- `composition/sealed_analysis.json` holds the composition comparisons.

The README inside that folder describes how these files were used to build later layouts. It is kept unchanged.

Every file here is listed in `../SHA256SUMS.txt`.
