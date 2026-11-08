"""Calibration-extension analysis (frozen with the protocol).

Estimator = count_sweep/recompute_results.py: stratified event means with pool-fraction weights; 4000 paired
bootstrap replicates (within-stratum event resampling + Dirichlet(pool_counts + 1) source weights), Jeffreys
quiet-baseline draws (D0C, quiet_d9.json), worse input mode taken inside each replicate.  The SAME bootstrap event
indices, weights and quiet draws are used for every condition of every wave (all paired on the 600-impact bank).

  python analyze_calext.py anchor          -> anchor regression only (typed + shallow, E1|standard)
  python analyze_calext.py all             -> every complete setting: excess risks, eta = deep/typed, shallow/typed
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
B, SEED = 4000, 20261108
CS = Path(r"C:\Users\admin\PACER_CALEXT\lsbp\data\count_sweep")
QUIET = json.loads((HERE / "quiet_d9.json").read_text())


def load_wave(name):
    S = json.loads((RES / name / "study.json").read_text())
    H, cnt = {}, {}
    for f in sorted((RES / name / "blocks").glob("*.json")):
        for r in json.loads(f.read_text())["rows"]:
            h = r["hist"]["D0C"]
            k = (r["cond"], r["mode"])
            H.setdefault(k, np.full(600, np.nan))[r["e"]] = 1.0 - h[0] / r["shots"]
            cnt[k] = cnt.get(k, 0) + 1
    return S, H, {k: v == 600 for k, v in cnt.items()}


class Boot:
    def __init__(self, EV):
        counts = EV["pool_counts"].astype(float)
        self.pi = counts / counts.sum()
        rng = np.random.default_rng(SEED)
        self.idx = [np.flatnonzero(EV["stratum"] == h) for h in range(3)]
        self.picks = [ii[rng.integers(0, len(ii), size=(B, len(ii)))] for ii in self.idx]
        self.w = rng.dirichlet(counts + 1, size=B)
        rng = np.random.default_rng(SEED + 1)
        self.q = {m: rng.beta(QUIET[f"9|{m}|D0C"]["fail"] + .5, QUIET[f"9|{m}|D0C"]["shots"] - QUIET[f"9|{m}|D0C"]["fail"] + .5, size=B)
                  for m in ("Z", "X")}

    def risk(self, H, c, boot=False):
        v = []
        for m in ("Z", "X"):
            x = H[c, m]
            if boot:
                v.append(sum(self.w[:, h] * x[self.picks[h]].mean(1) for h in range(3)) - self.q[m])
            else:
                v.append(sum(self.pi[h] * x[self.idx[h]].mean() for h in range(3)) - QUIET[f"9|{m}|D0C"]["p"])
        return np.max(v, axis=0) if boot else float(max(v))


def ci(a):
    return np.nanpercentile(a, [2.5, 97.5]).tolist()


def anchor():
    S, H, comp = load_wave("calext_A0")
    EV = dict(np.load(RES / "events" / S["events_file"]))
    bt = Boot(EV)
    old = {}
    with gzip.open(CS / "per_event_counts.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if r["condition"] in ("n168|typed", "n168|shallow"):
                old.setdefault((r["condition"], r["input_mode"]), np.full(600, np.nan))[int(r["event_index"])] = int(r["any_output_failures"]) / int(r["shots"])
    ref = json.loads((CS / "RESULTS.json").read_text())["decoders"]["D0C"]
    target = {"typed_paper": (ref["typed_excess"][4], ref["typed_excess_ci"][4], "n168|typed"),
              "shallow": (ref["target"], ref["target_ci"], "n168|shallow")}
    out = {"complete": all(comp.values())}
    ok = True
    for mk, (pt, iv, oname) in target.items():
        c = f"{mk}|E1|standard|T1|G12"
        new_pt, new_b = bt.risk(H, c), bt.risk(H, c, True)
        old_pt, old_b = bt.risk(old, oname), bt.risk(old, oname, True)
        diff_ci = ci(new_b - old_b)
        within = iv[0] <= new_pt <= iv[1]
        paired_ok = diff_ci[0] <= 0 <= diff_ci[1]
        ok &= within and paired_ok
        out[mk] = {"new": new_pt, "new_ci": ci(new_b), "count_sweep": pt, "count_sweep_recomputed": old_pt, "count_sweep_ci": iv,
                   "new_within_count_sweep_ci": within, "paired_diff": new_pt - old_pt, "paired_diff_ci": diff_ci,
                   "paired_diff_ci_contains_0": paired_ok}
    out["ANCHOR_PASS"] = bool(ok)
    return out


def all_settings():
    import calext_studies as cs
    H, comp, EVF = {}, {}, None
    for w in ("A0", "A", "B", "C"):
        if not (RES / f"calext_{w}" / "study.json").exists():
            continue
        S, h, c = load_wave(f"calext_{w}")
        H |= h
        comp |= c
        EVF = S["events_file"]
    bt = Boot(dict(np.load(RES / "events" / EVF)))
    rows = []
    for w, sets in cs.SETTINGS.items():
        for es, dec, tr, gap in sets:
            tag = cs.setting_tag(es, dec, tr, gap)
            names = {k: f"{k}|{tag}" for k in cs.MASKS}
            if not all(comp.get((n, m), False) for n in names.values() for m in ("Z", "X")):
                continue
            pt = {k: bt.risk(H, n) for k, n in names.items()}
            bb = {k: bt.risk(H, n, True) for k, n in names.items()}
            with np.errstate(divide="ignore", invalid="ignore"):
                eta_b, sh_b = bb["deep"] / bb["typed_paper"], bb["shallow"] / bb["typed_paper"]
            eci, sci = ci(eta_b), ci(sh_b)
            rows.append({"wave": w, "energy_scale": es, "decoupling": dec, "trap_multiplier": tr, "protected_gap_ghz": gap / 1e9,
                         "excess": pt, "excess_ci": {k: ci(v) for k, v in bb.items()},
                         "eta_deep_over_typed": pt["deep"] / pt["typed_paper"], "eta_ci": eci,
                         "eta_resolved_gt1": eci[0] > 1, "eta_resolved_lt1": eci[1] < 1,
                         "shallow_over_typed": pt["shallow"] / pt["typed_paper"], "shallow_over_typed_ci": sci,
                         "eta_boot_frac_nonfinite": float(np.mean(~np.isfinite(eta_b)))})
    return rows


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = anchor() if what == "anchor" else all_settings()
    (RES / f"analysis_{what}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
