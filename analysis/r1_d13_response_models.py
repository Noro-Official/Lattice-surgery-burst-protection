"""Standalone copy of the R1 part of scripts/pacer_postv2_ms/analyze_round3.py from Noro-Official/waveguide-please.

Only the input/output paths and the JSON writer differ (the unrelated R6 analysis is omitted); the estimator,
bootstrap and seed are unchanged. Usage: python analysis/r1_d13_response_models.py [output.json]

Round-3 analysis (written and committed before any round-3 outcome was read).

R6: the E1 V10 rows get the new Fc column from R6 (D0 and F from R6 must equal E1 exactly; any mismatch is reported
    and the E1 values are kept). Fc ratio and its paired change vs E1 V0/Fc use the shared E1 bootstrap replicates.
R1: d = 13 pilot with 5 strata (pool weights). Per (variant, decoder):
    E_L = max_mode( sum_h w_h * mean_h(fail rate) - q_mode ), Delta E = E_deep - E_typed (primary), rho = E_deep/E_typed.
    4,000 paired within-stratum event-bootstrap replicates, Dirichlet(pool counts + 1) weights, Jeffreys quiet draws.
    Equivalence (margin 1.33) if the rho 95% CI lies within [1/1.33, 1.33]; superiority if the Delta E CI excludes 0.
    Pilot sizing: per-stratum SD of the per-event (deep - typed) rate difference in the worst mode gives the Neyman
    allocation and the total N needed for a Delta E 95% half-width of 0.25 * |Delta E| (reported, not acted on here).
"""
from __future__ import annotations


import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/r1_d13_response_models"     # inputs: events.npz, rows.parquet, quiet.parquet
B = 4000
QMAP = {"D0": "D0", "F": "D0", "D0c": "D0c", "Fc": "D0c"}


def core(df, ev, qt, seed, decs):
    counts = ev["pool_counts"].astype(float)
    st = ev["stratum"]
    n, H = len(st), len(counts)
    idx = [np.nonzero(st == h)[0] for h in range(H)]
    qk = {(d, m): (int(qt[qt["mode"] == m][f"fail_{d}"].sum()), int(qt[qt["mode"] == m].shots.sum())) for d in ("D0", "D0c") for m in ("Z", "X")}
    rng = np.random.default_rng(seed)
    reps = [([rng.choice(i, len(i)) for i in idx], rng.dirichlet(counts + 1), rng.random(4)) for _ in range(B)]
    w0 = counts / counts.sum()

    def arr(dec, v):
        s = df[df.variant == v]
        F = np.zeros((2, 2, n))
        for i, mk in enumerate(("typed", "deep")):
            for j, md in enumerate(("Z", "X")):
                t = s[(s["mask"] == mk) & (s["mode"] == md)].set_index("event")
                F[i, j] = (t[f"fail_{dec}"] / t.shots).reindex(range(n)).values
        return F

    def stat(F, ix, w, q):
        per = np.stack([F[:, :, i].mean(axis=2) for i in ix], -1)
        ex = ((per * w).sum(-1) - q[None, :]).max(1)
        return ex

    out = []
    for dec in decs:
        if f"fail_{dec}" not in df.columns:
            continue
        q0 = np.array([qk[(QMAP[dec], m)][0] / qk[(QMAP[dec], m)][1] for m in ("Z", "X")])
        qr = [np.array([beta.ppf(r[2][j], qk[(QMAP[dec], m)][0] + .5, qk[(QMAP[dec], m)][1] - qk[(QMAP[dec], m)][0] + .5) for j, m in enumerate(("Z", "X"))]) for r in reps]
        vs = [v for v in sorted(df.variant.unique()) if df[df.variant == v][f"fail_{dec}"].notna().all()]
        base = None
        for v in (["V0_baseline"] if "V0_baseline" in vs else []) + [v for v in vs if v != "V0_baseline"]:
            F = arr(dec, v)
            ex = stat(F, idx, w0, q0)
            er = np.array([stat(F, r[0], r[1], q) for r, q in zip(reps, qr)])
            rho, rr = ex[1] / ex[0], er[:, 1] / er[:, 0]
            dE, dr = ex[1] - ex[0], er[:, 1] - er[:, 0]
            row = {"decoder": dec, "variant": v, "E_typed": float(ex[0]), "E_deep": float(ex[1]),
                   "E_typed_ci": np.quantile(er[:, 0], [.025, .975]).tolist(), "E_deep_ci": np.quantile(er[:, 1], [.025, .975]).tolist(),
                   "delta_E": float(dE), "delta_E_ci": np.quantile(dr, [.025, .975]).tolist(), "P_delta_E_le_0": float((dr <= 0).mean()),
                   "rho": float(rho), "rho_ci": np.quantile(rr, [.025, .975]).tolist(),
                   "frac_reps_nonpositive_E": float(((er <= 0).any(1)).mean())}
            lo, hi = row["rho_ci"]
            row["verdict"] = ("superior (Delta E > 0)" if row["delta_E_ci"][0] > 0 else "inferior (Delta E < 0)" if row["delta_E_ci"][1] < 0
                              else "equivalent within 1.33" if (lo >= 1 / 1.33 and hi <= 1.33 and row["frac_reps_nonpositive_E"] == 0) else "unresolved")
            if base is None:
                base = (ex, er)
            else:
                row["change_delta_E_vs_V0"] = float(dE - (base[0][1] - base[0][0]))
                row["change_delta_E_vs_V0_ci"] = np.quantile(dr - (base[1][:, 1] - base[1][:, 0]), [.025, .975]).tolist()
                ok = (er > 0).all(1) & (base[1] > 0).all(1)
                lr = np.log(rr[ok]) - np.log(base[1][ok, 1] / base[1][ok, 0])
                row["change_rho_vs_V0"] = float(rho / (base[0][1] / base[0][0]))
                row["change_rho_vs_V0_ci"] = np.exp(np.quantile(lr, [.025, .975])).tolist() if ok.mean() > .99 else None
            # pilot sizing (worst mode of the deep layout)
            m =int(np.argmax(((np.stack([F[:, :, i].mean(axis=2) for i in idx], -1) * w0).sum(-1) - q0[None, :])[1]))
            dpe = F[1, m] - F[0, m]
            sd = np.array([dpe[i].std(ddof=1) if len(i) > 1 else 0.0 for i in idx])
            ws = (w0 * sd).sum()
            if dE != 0 and ws > 0:
                ntot = (1.96 * ws / (0.25 * abs(dE))) ** 2
                row["pilot_sizing"] = {"worst_mode": "ZX"[m], "stratum_sd": sd.tolist(), "neyman_alloc_frac": (w0 * sd / ws).tolist(),
                                       "N_total_for_halfwidth_25pct_of_deltaE": float(ntot)}
            out.append(row)
    return out


def r1():
    return core(pd.read_parquet(DATA / "rows.parquet"), np.load(DATA / "events.npz"), pd.read_parquet(DATA / "quiet.parquet"),
                20301209, ("D0", "D0c", "F", "Fc"))


def main() -> int:
    import sys
    out = {"R1": r1()}
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results/r1_RESULTS.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json.loads(json.dumps(out, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else float(x))), indent=1))
    for r in out["R1"]:
        print(f"{r['variant']:<26} {r['decoder']:<4} rho {r['rho']:.2f} [{r['rho_ci'][0]:.2f}, {r['rho_ci'][1]:.2f}]  {r['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
