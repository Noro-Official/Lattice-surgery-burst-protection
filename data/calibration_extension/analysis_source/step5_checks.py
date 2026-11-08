"""Step 5 pre-campaign checks (no risk estimates are computed here).

1. Regression: channels() with default keys is bit-identical to the unmodified code path (original design() with
   uniform traps + the backed-up original pacer_continuous.trajectory_channels), 4 events x 3 masks, all 7 arrays.
2. Mask hashes of the Fig. 2a layouts (shallow / deep / typed) at d = 9, f = 0.3, T_M = 9.
3. Standard vs echo at energy_scale 1: total ancilla Pauli probability (bridge output, summed over cycles and ancillas).
4. Peak-shift table per energy scale on the 600-event bank (all sites at 12 GHz, i.e. a Willow-like device; and
   under the typed mask), peak over a 0.5 us - 3 ms grid.
Writes results/step5_checks.json.
"""
from __future__ import annotations

import importlib.util
import json
import sys

import numpy as np

import p1core as p
import studies

D, F, TM = 9, 0.3, 9
SCALES = (1, 3, 10, 30, 100)
BANK = p.HERE / "results" / "events" / "events_d9_seed78001_n20000_300-180-120_fp1.npz"
OUT = p.HERE / "results" / "step5_checks.json"


def orig_pc():
    spec = importlib.util.spec_from_file_location("src.analysis._pc_orig", p.HERE / "pacer_continuous.orig.py",
                                                  submodule_search_locations=None)
    m = importlib.util.module_from_spec(spec)
    m.__package__ = "src.analysis"
    sys.modules["src.analysis._pc_orig"] = m
    spec.loader.exec_module(m)
    return m


def masks():
    fam = studies.family(D, F, TM)
    return {k: fam[k] for k in ("shallow", "deep", "typed_paper")}


def main():
    out = {}
    EV = dict(np.load(BANK))
    x0s = EV["x0_fp1"]
    sch = p.schedule(D, TM)
    P = p.prepared("Z", D, sch)
    M = masks()
    out["mask_hashes"] = {k: p.mask_hash(v) for k, v in M.items()}
    out["mask_weak_counts"] = {k: int(v.sum()) for k, v in M.items()}
    print("mask hashes", out["mask_hashes"], flush=True)

    # 1. regression
    import campaign as cg
    pco = orig_pc()
    ok = True
    for name, w in M.items():
        for e in range(4):
            onset = cg.onset_of(float(EV["u_onset"][e]), D, TM)
            new = p.channels(P, x0s[e], w, onset)
            lay = P["layout"]
            nd = len(lay.data_ids)
            gaps = np.where(np.asarray(w, bool), p.WEAK_GAP_HZ, p.STRONG_GAP_HZ).astype(float)
            gaps[:nd] = p.STRONG_GAP_HZ
            des = {"gaps": gaps, "traps": np.full(len(gaps), p.physics()["pstar"].trapping_hz)}
            old = pco.trajectory_channels(np.asarray(x0s[e], float), des, p.physics()["pstar"], lay, p.adapter(), P["R"], int(onset))
            same = all(np.array_equal(a, b) for a, b in zip(new, old)) and len(new) == len(old) == 7
            ok &= same
    out["regression_bit_identical"] = bool(ok)
    print("regression bit-identical:", ok, flush=True)

    # 3. standard vs echo, total ancilla Pauli probability, x1, all 600 events, typed mask and deep mask
    from src.physics.pacer_qp_bridge_v1 import channels_from_site_design
    pstar = p.physics()["pstar"]
    nd = len(P["layout"].data_ids)
    times = (np.arange(P["R"] + 1) + 0.5) * pstar.cycle_time_s
    se = {}
    for name in ("typed_paper", "deep"):
        des = p.design(D, M[name])
        for sc in (1, 100):
            tot = {}
            for dec in ("standard", "echo"):
                s = 0.0
                for e in range(len(x0s)):
                    ch = channels_from_site_design(x0s[e] * sc, times, pstar, gap_difference_hz=des["gaps"],
                                                   trapping_hz=des["traps"], decoupling=dec)
                    s += float(np.asarray(ch.pauli_ancilla)[:, nd:, :].sum())
                tot[dec] = s
            se[f"{name}|x{sc}"] = {"standard": tot["standard"], "echo": tot["echo"], "rel_diff": (tot["echo"] - tot["standard"]) / tot["standard"]}
            print("std vs echo", name, sc, se[f"{name}|x{sc}"], flush=True)
    out["standard_vs_echo_total_ancilla_pauli"] = se

    # 4. peak-shift table
    from src.physics.pacer_qp_bridge_v1 import Platform, frequency_shift_coefficient, qp_density_with_generation
    tgrid = np.unique(np.concatenate([np.geomspace(0.5e-6, 3e-3, 160), (np.arange(60) + 0.5) * pstar.cycle_time_s]))
    table = {}
    for label, gaps in (("all_12GHz", np.full(x0s.shape[1], p.STRONG_GAP_HZ)), ("typed_mask", p.design(D, M["typed_paper"])["gaps"])):
        a = np.array([frequency_shift_coefficient(pstar, Platform("g", float(g), "standard")) for g in gaps])
        traps = np.full(len(gaps), pstar.trapping_hz)
        for sc in SCALES:
            peak_ev, n200, peak_med_site = [], [], []
            for e in range(len(x0s)):
                x = qp_density_with_generation(x0s[e] * sc, tgrid, pstar, trapping_hz=traps)
                pk = (a[None, :] * pstar.qubit_frequency_hz * x).max(0)        # per-site peak |shift| (Hz)
                peak_ev.append(pk.max())
                n200.append(int((pk > 200e3).sum()))
            peak_ev, n200 = np.array(peak_ev), np.array(n200)
            w = EV["pool_counts"] / EV["pool_counts"].sum()
            # population-weighted median: weight each event by its stratum's pool fraction / stratum sample size
            st = EV["stratum"]
            ew = np.array([w[s] / (st == s).sum() for s in st])

            def wq(v, q):
                o = np.argsort(v)
                c = np.cumsum(ew[o])
                return float(v[o][np.searchsorted(c, q * c[-1])])
            table[f"{label}|x{sc}"] = {"median_peak_khz": wq(peak_ev, 0.5) / 1e3, "p90_peak_khz": wq(peak_ev, 0.9) / 1e3,
                                       "median_sites_above_200khz": wq(n200.astype(float), 0.5),
                                       "unweighted_median_peak_khz": float(np.median(peak_ev)) / 1e3}
            print("peak", label, sc, table[f"{label}|x{sc}"], flush=True)
    out["peak_shift_table"] = table
    out["peak_shift_note"] = ("Population-weighted by stratum pool fractions; peak over a 0.5 us-3 ms grid of "
                              "qp_density_with_generation with base trapping; shift = a(gap) f_q x_qp.")
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
