# PACER protection-count refinement: completed run
6 November 2026

**Main result:** Under correlated background matching at distance nine, log-risk interpolation estimates **118.8 protected ancillas** [117.0, 133.1] to reach the shallow layout's risk at 168 sites. The estimated saving is **49.2 protected ancillas** [34.9, 51.0]. Intervals are nominal 95% paired-bootstrap intervals.

The crossing is bracketed by the **117- and 123-site sampled points**. At 123 protected sites, the typed point estimate is below the reference; that is an observed point-estimate saving of 45 sites. The fractional 118.8-site crossing is interpolated, not a directly implemented layout or a separately simulated 119-site configuration.

## Sampled results
Excess risk is the larger of the two input-mode increases over quiet operation, averaged over the impact distribution. It is conditional on an impact.

| Layout | Protected ancillas | Excess risk | Nominal 95% interval |
|---|---:|---:|---|
| Typed | 117 | 0.000828674316 | [0.00069111551, 0.00116376053] |
| Typed | 123 | 0.0006190389 | [0.000521728527, 0.000714719786] |
| Typed | 132 | 0.000509989421 | [0.000422272964, 0.000598116894] |
| Typed | 150 | 0.000307027181 | [0.000239359408, 0.000370627962] |
| Typed | 168 | 0.000242411296 | [0.000174087439, 0.000308130738] |
| Typed | 204 | 0.000111063639 | [5.44863893e-05, 0.00016357414] |
| Shallow reference | 168 | 0.000757873535 | [0.000500477046, 0.00105542875] |

## Text for the paper
Under correlated background matching at distance nine, a refined protection-count sweep brackets the typed band's crossing of the shallow band's 168-site reference risk between 117 and 123 protected ancillas. Log-risk interpolation gives an estimated requirement of 118.8 [117.0, 133.1] protected ancillas, corresponding to a saving of 49.2 [34.9, 51.0] sites. The nominal 95% intervals use a paired stratified event bootstrap with quiet-baseline uncertainty. The interval is constrained by the rule's minimum feasible count of 117, which is reached in 28.1% of bootstrap draws. This is a refinement on the existing 600-impact bank rather than an independent impact confirmation.

## Design and scope
- Distance nine; 240 eligible ancillas; five-layer typed band; declared type budgets unchanged.
- Typed counts: 117, 123, 132, 150, 168 and 204; shallow reference: 168. All seven conditions were regenerated in the same runtime.
- Requested counts of 96 and 114 are infeasible for this fixed band and type allocation. They were not simulated. Testing them would require a different layout rule.
- 600 saved impacts, 1,024 trajectories per event, input mode and layout; 8,601,600 retained trajectories in total. All 2,100 checkpoint blocks completed.
- Primary decoder: correlated background matching (D0C). No refined oracle or other-decoder result is claimed. An outcome-independent runtime amendment narrowed the initially broader plan to this primary endpoint.
- Schedule: 12 pre-rounds, nine-round ZZ merge, nine-round gap, nine-round XX merge, helper readout, 12 post-rounds. Background SI1000 parameter 0.001; strong/weak gaps 12/5 GHz.
- Event-pool seed 78001; trajectory seed base 78002; stratum allocation 300/180/120; weights from the saved 20,000-event source pool.
- 4,000 paired within-stratum event-bootstrap draws, Dirichlet source-weight draws and shared Jeffreys quiet-baseline draws. The worse input mode is recomputed in each draw. Piecewise log-risk interpolation uses the first crossing, without a forced monotone fit or extrapolation.
- The observed six-point curve decreases throughout. The lower feasible boundary is reached in 28.1% of draws; the upper grid boundary in 0%. The saving interval's 51-site endpoint is a feasible-domain boundary, not an unconstrained optimum or proof of a universal maximum saving.
- The grid was refined after examining the coarse result and reuses its impacts. Intervals are nominal and have no independent coverage calibration. Protected count is not measured fabrication cost.

## Provenance and verification
The original experiment archive supplied an explicit 12/12 schedule and summaries matching the manuscript inputs byte for byte. A four-event seeded control was not bit-identical in this runtime. Stim permits machine-dependent seeded samples; that does not establish the specific cause of the discrepancy. Therefore, the new curve and its shallow reference were regenerated together rather than pooled with old trajectory estimates.

The exported per-event D0C histograms sum to 1,024 in all 8,400 event/mode/layout records. Weighted means and the point crossing were independently recomputed from these counts. `recompute_results.py` recomputes the full bootstrap result using only NumPy and the exported data. Historical other-decoder values and the archived 38% placement result are separate studies.

## Files
- `paper_insert.tex`: insertion-ready paragraph and numeric table.
- `RESULTS.json`, `sampled_risks.csv`: exact machine-readable estimates.
- `per_event_counts.csv.gz`, event-bank NPZ, `quiet_d9.json`: data for independent recomputation.
- `recompute_results.py`: standalone analysis check (Python and NumPy).
- Protocol, amendment, control diagnostic and original analysis sources: audit trail. The original analysis sources refer to the full simulation tree; use the standalone checker to analyze this export without that tree.
