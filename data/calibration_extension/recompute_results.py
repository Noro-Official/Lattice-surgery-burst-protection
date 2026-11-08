"""Standalone check (Python + NumPy): recompute RESULTS.json and ANCHOR.json's new values from per_event_counts.csv.gz."""
from pathlib import Path
import csv, gzip, json
import numpy as np

R = Path(__file__).resolve().parent
A = json.loads((R / "RESULTS.json").read_text())
Q = json.loads((R / "quiet_d9.json").read_text())
E = dict(np.load(R / "events_d9_seed78001_n20000_300-180-120_fp1.npz"))
H = {}
with gzip.open(R / "per_event_counts.csv.gz", "rt") as f:
    for r in csv.DictReader(f):
        H.setdefault((r["condition"], r["input_mode"]), np.full(600, np.nan))[int(r["event_index"])] = int(r["any_output_failures"]) / int(r["shots"])
B, seed = 4000, 20261108
counts = E["pool_counts"].astype(float); pi = counts / counts.sum()
rng = np.random.default_rng(seed)
idx = [np.flatnonzero(E["stratum"] == h) for h in range(3)]
picks = [ii[rng.integers(0, len(ii), size=(B, len(ii)))] for ii in idx]
w = rng.dirichlet(counts + 1, size=B)
rng = np.random.default_rng(seed + 1)
q = {m: rng.beta(Q[f"9|{m}|D0C"]["fail"] + .5, Q[f"9|{m}|D0C"]["shots"] - Q[f"9|{m}|D0C"]["fail"] + .5, size=B) for m in ("Z", "X")}


def risk(c, boot=False):
    v = []
    for m in ("Z", "X"):
        x = H[c, m]
        v.append(sum(w[:, h] * x[picks[h]].mean(1) for h in range(3)) - q[m] if boot else
                 sum(pi[h] * x[idx[h]].mean() for h in range(3)) - Q[f"9|{m}|D0C"]["p"])
    return np.max(v, axis=0) if boot else float(max(v))


out = []
for row in A:
    tag = f"E{row['energy_scale']:g}|{row['decoupling']}|T{row['trap_multiplier']:g}|G{row['protected_gap_ghz']:g}"
    pt = {k: risk(f"{k}|{tag}") for k in ("shallow", "deep", "typed_paper")}
    bb = {k: risk(f"{k}|{tag}", True) for k in ("shallow", "deep", "typed_paper")}
    with np.errstate(divide="ignore", invalid="ignore"):
        eci = np.nanpercentile(bb["deep"] / bb["typed_paper"], [2.5, 97.5])
    assert all(abs(pt[k] - row["excess"][k]) < 1e-15 for k in pt), tag
    assert np.allclose(eci, row["eta_ci"], rtol=0, atol=1e-12), tag
    out.append({"setting": tag, "eta": pt["deep"] / pt["typed_paper"], "eta_ci": eci.tolist()})
print(json.dumps({"settings": out, "verified": True}, indent=1))
