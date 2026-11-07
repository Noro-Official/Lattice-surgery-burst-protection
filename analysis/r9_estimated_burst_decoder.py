"""Standalone copy of scripts/pacer_postv2_r3/r9/analyze_r9.py from Noro-Official/waveguide-please.

Only the input/output paths and the JSON writer differ from the original; the estimator, bootstrap and seed are
unchanged. Usage: python analysis/r9_estimated_burst_decoder.py [output.json]

R9 analysis (pre-registered with the runner, before any outcome is sampled).

Estimand per (layout L, variant V, decoder X), as the manuscript:
  m_X(L, V, mode) = pool-weighted stratified mean over events of the per-event failure fraction
  E_X(L, V)       = max over input modes (Z, X) of [ m_X(L, V, mode) - q_ref(mode) ]
  q_ref           = burst-free failure of D0 for D0, F, Ff, E, L1, L3, K (they all fall back to / share the uncorrelated
                    background matcher); of D0c for D0c and Fc.  Burst-free penalties of the decoders that can
                    false-alarm (E, K) are reported SEPARATELY (below) and not folded into E_X.
Contrasts:
  dE_X(V)   = E_X(deep, V) - E_X(typed, V);   rho_X(V) = E_X(deep, V) / E_X(typed, V)
  R_X(L, V) = E_D0(L, V) - E_X(L, V)                      (burst-risk reduction relative to D0)
  g_X(L, V) = R_X(L, V) / R_F(L, V)                        (fraction of the oracle F gain retained; F = uncorrelated oracle,
                                                           the same matcher family as every realistic decoder)
  quiet penalty dq_X(L, V, mode) = (harms - rescues) / shots on burst-free shots (paired with D0);
  break-even burst probability per operation lambda* = dq / (R + dq): below lambda*, the decoder's burst-free cost
  outweighs its burst benefit.
Uncertainty: 4,000 paired within-stratum event-bootstrap replicates (shared over layouts, variants, decoders, modes),
Dirichlet(pool counts + 1) stratum weights, Jeffreys Beta(k + 1/2, n - k + 1/2) quiet draws; quiet penalties from a
Dirichlet(harm + 1/2, rescue + 1/2, other + 1/2) draw of the paired burst-free counts.  Percentile 95% intervals.

PRIMARY (declared): under V0, for the syndrome-estimated decoder E: g_E(typed), g_E(deep) and R_E with 95% intervals,
and dE_E(V0).  Everything else is secondary.  No pass mark is set for g; a reduction is "resolved" if the R interval
excludes 0, a layout difference if the dE interval excludes 0.  Oracle-only gains are reported as such.
"""
from __future__ import annotations


import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/r9_estimated_burst_decoder"     # inputs: events.npz, rows.parquet, quiet.parquet
RESULT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/r9_RESULTS.json"
B = 4000
BOOT_SEED = 20351199
LAYOUTS = ("typed", "deep")
VARIANTS = ("V0_baseline", "V10_phase_echo_measured")
DECS = ("D0", "D0c", "F", "Fc", "Ff", "E", "L1", "L3", "K")
REALISTIC = ("E", "L1", "L3", "K")
QREF = {d: ("D0c" if d in ("D0c", "Fc") else "D0") for d in DECS}


def ci(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.percentile(x, [2.5, 97.5])] if len(x) else [None, None]


def main() -> int:
    ev = np.load(OUT / "events.npz")
    st, counts = ev["stratum"], ev["pool_counts"].astype(float)
    n = len(st)
    df = pd.read_parquet(OUT / "rows.parquet")
    q = pd.read_parquet(OUT / "quiet.parquet")
    complete = len(df) == 2 * n * len(LAYOUTS) * len(VARIANTS)
    events_present = np.unique(df.event.values)
    # analysis on the events present in BOTH modes for every layout / variant (complete-case; reported)
    k = df.groupby("event").size()
    use = np.array(sorted(k[k == 2 * len(LAYOUTS) * len(VARIANTS)].index))
    idx = [use[st[use] == h] for h in range(3)]
    if any(len(i) == 0 for i in idx):
        print("a stratum has no complete event; nothing to analyse")
        return 1
    rng = np.random.default_rng(BOOT_SEED)
    p_hat = counts / counts.sum()
    Pd = rng.dirichlet(counts + 1.0, size=B)
    bidx = [ii[rng.integers(0, len(ii), size=(B, len(ii)))] for ii in idx]
    qp, qb = {}, {}
    for mode in ("Z", "X"):
        qm = q[q["mode"] == mode]
        ns = float(qm.shots.sum())
        for d in ("D0", "D0c"):
            kq = float(qm[f"fail_{d}"].sum())
            qp[(d, mode)] = kq / ns
            qb[(d, mode)] = rng.beta(kq + 0.5, ns - kq + 0.5, size=B)
    V = {}
    for (lay, var, mode), g in df.groupby(["layout", "variant", "mode"]):
        g = g.set_index("event")
        for d in DECS:
            v = np.full(n, np.nan)
            v[g.index.values] = g[f"fail_{d}"].values / g["shots"].values
            V[(lay, var, d, mode)] = v

    def E(lay, var, d, boot=False):
        vals = []
        for mode in ("Z", "X"):
            v = V[(lay, var, d, mode)]
            if boot:
                vals.append(sum(Pd[:, h] * v[bidx[h]].mean(1) for h in range(3)) - qb[(QREF[d], mode)])
            else:
                vals.append(float(sum(p_hat[h] * v[idx[h]].mean() for h in range(3))) - qp[(QREF[d], mode)])
        return np.max(vals, axis=0) if boot else max(vals)

    res = {"complete": bool(complete), "events_total": int(n), "events_analysed": int(len(use)), "events_present": int(len(events_present)),
           "events_analysed_per_stratum": [int(len(i)) for i in idx], "pool_counts": counts.astype(int).tolist(),
           "quiet": {f"{d}|{m}": {"rate": qp[(d, m)], "shots": int(q[q['mode'] == m].shots.sum())} for (d, m) in qp},
           "excess": {}, "layout_contrast": {}, "reduction_vs_D0": {}, "quiet_penalty": {}, "diagnostics": {}}
    for var in VARIANTS:
        for d in DECS:
            for lay in LAYOUTS:
                res["excess"][f"{var}|{d}|{lay}"] = {"E": E(lay, var, d), "ci95": ci(E(lay, var, d, True))}
            bt, bd = E("typed", var, d, True), E("deep", var, d, True)
            pt, pd_ = E("typed", var, d), E("deep", var, d)
            with np.errstate(divide="ignore", invalid="ignore"):
                res["layout_contrast"][f"{var}|{d}"] = {"dE_deep_minus_typed": pd_ - pt, "dE_ci95": ci(bd - bt),
                                                        "ratio_deep_over_typed": pd_ / pt if pt != 0 else None, "ratio_ci95": ci(bd / bt),
                                                        "resolved": bool(np.all(np.array(ci(bd - bt)) > 0) or np.all(np.array(ci(bd - bt)) < 0))}
        for lay in LAYOUTS:
            b0, p0 = E(lay, var, "D0", True), E(lay, var, "D0")
            bF, pF = E(lay, var, "F", True), E(lay, var, "F")
            bFc, pFc = E(lay, var, "Fc", True), E(lay, var, "Fc")
            res["reduction_vs_D0"][f"{var}|{lay}|F(oracle)"] = {"R": p0 - pF, "ci95": ci(b0 - bF), "relative": (p0 - pF) / p0 if p0 else None}
            res["reduction_vs_D0"][f"{var}|{lay}|Fc(oracle, corr.)"] = {"R": p0 - pFc, "ci95": ci(b0 - bFc), "relative": (p0 - pFc) / p0 if p0 else None}
            for d in REALISTIC:
                bX, pX = E(lay, var, d, True), E(lay, var, d)
                with np.errstate(divide="ignore", invalid="ignore"):
                    g = (b0 - bX) / (b0 - bF)
                res["reduction_vs_D0"][f"{var}|{lay}|{d}"] = {"R": p0 - pX, "ci95": ci(b0 - bX), "relative": (p0 - pX) / p0 if p0 else None,
                                                              "resolved_positive": bool(ci(b0 - bX)[0] is not None and ci(b0 - bX)[0] > 0),
                                                              "retained_fraction_of_oracle_gain": (p0 - pX) / (p0 - pF) if p0 != pF else None,
                                                              "retained_ci95": ci(g),
                                                              "oracle_gain_resolved": bool(ci(b0 - bF)[0] is not None and ci(b0 - bF)[0] > 0)}
    # burst-free penalty of the decoders that can false-alarm
    for lay in LAYOUTS:
        for var in VARIANTS:
            for d in ("E", "K"):
                for mode in ("Z", "X"):
                    qm = q[q["mode"] == mode]
                    h, r, ns = float(qm[f"{lay}|{var}|harm_{d}"].sum()), float(qm[f"{lay}|{var}|rescue_{d}"].sum()), float(qm.shots.sum())
                    dd = rng.dirichlet([h + 0.5, r + 0.5, ns - h - r + 0.5], size=B)
                    dq = (h - r) / ns
                    R = res["reduction_vs_D0"][f"{var}|{lay}|{d}"]["R"]
                    res["quiet_penalty"][f"{var}|{lay}|{d}|{mode}"] = {
                        "alarms": int(qm[f"{lay}|{d}_alarms"].sum()), "shots": int(ns), "harms": int(h), "rescues": int(r),
                        "dq": dq, "dq_ci95": ci(dd[:, 0] - dd[:, 1]),
                        "break_even_burst_prob_per_op": (dq / (R + dq)) if (dq > 0 and R > 0) else None}
    res["diagnostics"]["F_equals_Ff_arm_fraction"] = float(df.F_eq_Ff.mean())
    res["diagnostics"]["fail_F_minus_fail_Ff_total"] = int((df.fail_F - df.fail_Ff).sum())
    acc = df[df.variant == "V0_baseline"].groupby(["layout", "stratum"]).apply(
        lambda g: pd.Series({"accept_frac": g.E_accepted.sum() / g.shots.sum(), "onset_abs_err_median": float(np.nanmedian(g.E_onset_abs_err_median)),
                             "pos_err_mm_median": float(np.nanmedian(g.E_pos_err_mm_median)), "distinct_estimates_mean": float(g.E_distinct.mean())}), include_groups=False)
    res["diagnostics"]["E_estimator_V0"] = {f"{a}|h{b}": r.to_dict() for (a, b), r in acc.iterrows()}
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(json.loads(json.dumps(res, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x))), indent=1))
    # console table
    print(f"R9 | complete {complete} | events analysed {len(use)} of {n} ({res['events_analysed_per_stratum']})")
    for var in VARIANTS:
        print(f"\n{var}: excess risk given a burst (worst mode) and deep - typed")
        for d in DECS:
            t, dp, c = res["excess"][f"{var}|{d}|typed"], res["excess"][f"{var}|{d}|deep"], res["layout_contrast"][f"{var}|{d}"]
            print("  %-3s typed %.3e  deep %.3e  dE %+.3e [%+.2e, %+.2e]  ratio %s [%s, %s]" % (
                d, t["E"], dp["E"], c["dE_deep_minus_typed"], *c["dE_ci95"], f"{c['ratio_deep_over_typed']:.2f}" if c["ratio_deep_over_typed"] else "-",
                *(f"{x:.2f}" if x is not None else "-" for x in c["ratio_ci95"])))
        for lay in LAYOUTS:
            print(f"  {lay}: reduction vs D0 (retained fraction of oracle F gain)")
            for d in ("F(oracle)", "Fc(oracle, corr.)") + REALISTIC:
                r = res["reduction_vs_D0"][f"{var}|{lay}|{d}"]
                extra = "" if d not in REALISTIC else "  g %s [%s, %s]" % (
                    f"{r['retained_fraction_of_oracle_gain']:.2f}" if r["retained_fraction_of_oracle_gain"] is not None else "-",
                    *(f"{x:.2f}" if x is not None else "-" for x in r["retained_ci95"]))
                print("    %-18s R %+.3e [%+.2e, %+.2e]%s" % (d, r["R"], *r["ci95"], extra))
    print("\nburst-free penalty (harms - rescues per shot):")
    for kk, v in res["quiet_penalty"].items():
        print(f"  {kk}: alarms {v['alarms']}/{v['shots']}, dq {v['dq']:+.2e} [{v['dq_ci95'][0]:+.1e}, {v['dq_ci95'][1]:+.1e}], break-even {v['break_even_burst_prob_per_op']}")
    print("\ndiagnostics:", json.dumps(res["diagnostics"], indent=1, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
