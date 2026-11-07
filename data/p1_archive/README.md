# Original P1 cohorts (31 October 2026 publication package)

These files are copied byte for byte from `publication_package_20261031_revised/experiments/pacer_p1/results/`, the
archive the manuscript's figures were built from. `PACKAGE_SHA256SUMS.txt` holds that package's own checksum
entries for these cohorts, and `reproduce.py` checks every file here against it.

| Folder / file | Paper use | Design |
|---|---|---|
| `precision_d9/` | Fig. 2a, abstract (η = 1.62, 38% reduction) | d = 9; 3,200 impacts in five strata from a 60,000-event pool (seed 85001); shallow, deep and typed layouts; D0, D0C, F, FC; 4,096 shots per event, mode and layout |
| `precision_d13/`, `precision_analysis.json` | Fig. 5a (oracle ratios 0.96, 1.12) | d = 13; 1,200 impacts in five strata from a 40,000-event pool (seed 86001); same layouts and decoders; 4,096 shots |
| `regime/` | Fig. 5b,c (burst-duration caps; five-round ratio 0.60) | d = 9; 600 impacts (seed 76001); 72 conditions = six layouts × report strength, burst duration (2, 5, 10, 20 rounds) and footprint scale; D0, F; 1,024 shots |
| `events/` | Event banks for the three cohorts | strata, pool counts, impact parameters |
| `quiet_d9.json`, `quiet_d9_p2.json`, `quiet_d13_p2.json` | Burst-free baselines | `regime` uses `quiet_d9.json`; the precision studies use the `_p2` tables |

`precision_d9/analysis.json` and the `precision_d9` entry of `precision_analysis.json` are identical. The d = 13
analysis is stored only in `precision_analysis.json`.

**Per-event counts not yet included.** The package also holds 15,600 per-event block files for these cohorts
(4,800 for `precision_d9`, 3,600 for `precision_d13`, 7,200 for `regime`), listed with their checksums in
`PACKAGE_SHA256SUMS.txt`. Once they are added, the stored analyses here can be recomputed from raw counts, as for the
other studies in this repository.
