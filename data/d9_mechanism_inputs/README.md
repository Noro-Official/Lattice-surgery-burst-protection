# Manuscript inputs for round 3 (copied from the 4 Nov 2026 manuscript package, figure_inputs/)

- `site_scores_d9_tm9.json`: isolated-site scores, d = 9, T_M = 9. `W["full|D0"][site]` holds 4 values; `sites` are geometry site indices.
- **Priority layout rule (manuscript Methods, verbatim):**
  - "For isolated-site scores, we apply burst channels to one weak site at a reference dose: the median epicentral weak-site dose among the heaviest 5% of development impacts. The site's smallest W_D0 across merges defines its danger score. Priority leaves the least dangerous sites weak within each check type."
  - Reconstruction:
    - danger(site) = min of its W["full|D0"] values (smaller W = more dangerous);
    - within each check type (X-anc, Z-anc), mark weak the sites with the LARGEST min-W;
    - per-type counts equal the typed layout's per-type counts (`cw.budget_counts(72)`, or the typed mask's split);
    - total 72.
  - Verify the per-type split equals typed, and record ties and their tie-break (site index) in DEVIATIONS.md.
- `composition/`, `factorial/`, `baselines/`, `exchange_test/`: the manuscript's analysis JSONs, including the line-model predictions in composition/prediction.json and the factorial cohort analysis.
