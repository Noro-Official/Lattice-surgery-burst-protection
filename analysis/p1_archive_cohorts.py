"""Recompute the archived P1 cohorts behind Fig. 2a (precision_d9), Fig. 5a (precision_d13) and Fig. 5b,c (regime).

Standalone re-implementation of the loader and estimator in experiments/pacer_p1/{campaign.py, analyze.py,
analyze_precision.py} of the 31 October 2026 publication package (publication_package_20261031_revised). The
estimators, bootstrap sizes and seeds are the originals:

  per-event risk        1 - hist[0] / shots          (hist[0] = shots with every logical output correct)
  stratum-weighted mean sum_h pi_h * mean over events in stratum h, pi = pool counts / total
  bootstrap             4,000 within-stratum event resamples + Dirichlet(pool counts + 1) weights, seed 20261029
  excess risk           worse input mode of (mean - quiet baseline)
  precision studies     quiet uncertainty from Jeffreys Beta draws (seed 20261030 + d) on quiet_d{d}_p2.json
  regime study          quiet point values from quiet_d9.json (as analyze.generic)

Usage: python analysis/p1_archive_cohorts.py [output.json]
Each cohort's blocks.jsonl.gz holds the original block files byte for byte; every line's SHA-256 is checked against
the package's own SHA256SUMS entries (data/p1_archive/PACKAGE_SHA256SUMS.txt).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "p1_archive"
B = 4000
DECS = ("D0", "D0C", "F", "FC")
MARGIN = {9: 1.10, 13: 1.33}


def _ci(x):
    return [float(v) for v in np.nanpercentile(x, [2.5, 97.5])]


def package_sums():
    sums = {}
    for line in (DATA / "PACKAGE_SHA256SUMS.txt").read_text().splitlines():
        h, f = line.split(maxsplit=1)
        sums[f.lstrip("*").removeprefix("./")] = h
    return sums


def load(name, sums):
    """{(cond, decoder, mode): per-event risk array}, as campaign.load; verifies every block against the package."""
    folder = DATA / name
    S = json.loads((folder / "study.json").read_text())
    EV = dict(np.load(DATA / "events" / S["events_file"]))
    n = len(EV["stratum"])
    H, cnt, nblocks = {}, {}, 0
    with gzip.open(folder / "blocks.jsonl.gz", "rt", encoding="utf-8") as fh:
        lines = [json.loads(line) for line in fh]
    for rec in sorted(lines, key=lambda r: r["file"]):          # campaign.load reads blocks in sorted order
        raw = rec["content"].encode("utf-8")
        h = hashlib.sha256(raw).hexdigest()
        key = f"experiments/pacer_p1/results/{name}/blocks/{rec['file']}"
        if h != rec["sha256"] or sums.get(key) != h:
            raise ValueError(f"block {name}/{rec['file']} does not match the package checksum")
        nblocks += 1
        for r in json.loads(raw)["rows"]:
            for dec, hist in r["hist"].items():
                k = (r["cond"], dec, r["mode"])
                H.setdefault(k, np.full(n, np.nan))[r["e"]] = 1.0 - hist[0] / r["shots"]
                cnt[k] = cnt.get(k, 0) + 1
    expected = sum(1 for f in sums if f.startswith(f"experiments/pacer_p1/results/{name}/blocks/"))
    if nblocks != expected:
        raise ValueError(f"{name}: {nblocks} blocks present, package lists {expected}")
    return {"S": S, "EV": EV, "H": H, "n": n, "complete": all(v == n for v in cnt.values()), "blocks": nblocks}


class Estimator:
    def __init__(self, D, B=B, seed=20261029):
        st = D["EV"]["stratum"]
        counts = D["EV"]["pool_counts"].astype(float)
        self.pi = counts / counts.sum()
        rng = np.random.default_rng(seed)
        self.K = len(counts)
        self.idx = [np.nonzero(st == h)[0] for h in range(self.K)]
        self.bidx = [ii[rng.integers(0, len(ii), size=(B, len(ii)))] for ii in self.idx]
        self.Pd = rng.dirichlet(counts + 1.0, size=B)
        self.H = D["H"]

    def mean(self, cond, dec, mode, boot=False):
        v = self.H[(cond, dec, mode)]
        if boot:
            return sum(self.Pd[:, h] * v[self.bidx[h]].mean(1) for h in range(self.K))
        return float(sum(self.pi[h] * v[self.idx[h]].mean() for h in range(self.K)))

    def risk(self, cond, dec, boot=False, quiet=None):
        vals = [self.mean(cond, dec, m, boot) - (0.0 if quiet is None else quiet[(dec, m)]) for m in ("Z", "X")]
        return np.max(vals, axis=0) if boot else max(vals)


def quiet_table(name):
    """{(tm, mode, decoder): (fail, shots)}; F shares D0's baseline and FC shares D0C's (quiet.table)."""
    out = {}
    for k, v in json.loads((DATA / name).read_text()).items():
        tm, mode, dec = k.split("|")
        for d_ in (dec, {"D0": "F", "D0C": "FC"}[dec]):
            out[(int(tm), mode, d_)] = (v["fail"], v["shots"])
    return out


def classify(ci_d, ci_r, m):
    if ci_d[0] > 0:
        return "typed better"
    if ci_d[1] < 0:
        return "comparator better"
    if 1 / m <= ci_r[0] and ci_r[1] <= m:
        return f"equivalent within {m:g}"
    return "unresolved"


def precision(name, sums):
    D = load(name, sums)
    d, tm = D["S"]["d"], D["S"]["conditions"][0]["tm"]
    E = Estimator(D)
    rng = np.random.default_rng(20261030 + d)          # analyze_precision.quiet_draws
    C = quiet_table(f"quiet_d{d}_p2.json")
    qp, qb = {}, {}
    for dec in DECS:
        for m in ("Z", "X"):
            f, n = C[(tm, m, dec)]
            qp[(dec, m)] = f / n
            qb[(dec, m)] = rng.beta(f + 0.5, n - f + 0.5, size=B)
    out = {"events": D["n"], "shots": D["S"]["shots"], "blocks": D["blocks"], "complete": D["complete"], "comparisons": {}}
    for dec in DECS:
        rt, rtb = E.risk("typed_paper", dec, quiet=qp), E.risk("typed_paper", dec, True, qb)
        for a in ("shallow", "deep"):
            ra, rab = E.risk(a, dec, quiet=qp), E.risk(a, dec, True, qb)
            ci_d, ci_r = _ci(rab - rtb), _ci(rab / rtb)
            out["comparisons"][f"{a} vs typed|{dec}"] = {"diff": ra - rt, "diff_ci": ci_d, "ratio": ra / rt if rt > 0 else None,
                                                         "ratio_ci": ci_r, "class": classify(ci_d, ci_r, MARGIN[d])}
    return out


def regime(sums, keys):
    """Contrasts a - b (excess risk, worse mode) with point quiet baselines, as analyze.generic."""
    D = load("regime", sums)
    E = Estimator(D)
    C = quiet_table("quiet_d9.json")
    tm_of = {c["name"]: c["tm"] for c in D["S"]["conditions"]}
    out = {"events": D["n"], "shots": D["S"]["shots"], "blocks": D["blocks"], "complete": D["complete"], "contrast": {}}
    for key in keys:
        pair, dec = key.rsplit("|", 1)
        a, b = pair.split(" - ")
        q = {(dec, m): C[(tm_of[a], m, dec)][0] / C[(tm_of[a], m, dec)][1] for m in ("Z", "X")}
        ra, rb = E.risk(a, dec, quiet=q), E.risk(b, dec, quiet=q)
        ba, bb = E.risk(a, dec, True, q), E.risk(b, dec, True, q)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = ba / bb
        out["contrast"][key] = {"a": ra, "b": rb, "diff": ra - rb, "diff_ci": np.percentile(ba - bb, [2.5, 97.5]).tolist(),
                                "ratio": ra / rb if rb > 0 else float("inf"), "ratio_ci": np.nanpercentile(ratio, [2.5, 97.5]).tolist()}
    return out


def main():
    sums = package_sums()
    stored_regime = json.loads((DATA / "regime" / "analysis.json").read_text())
    res = {"precision_d9": precision("precision_d9", sums),
           "precision_d13": precision("precision_d13", sums),
           "regime": regime(sums, list(stored_regime["contrast"]))}
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "p1_archive_RESULTS.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(res, indent=1, default=float))
    for k in ("precision_d9", "precision_d13"):
        r = res[k]
        print(f"{k}: {r['events']} events x {r['shots']} shots, {r['blocks']} blocks")
        for c, v in r["comparisons"].items():
            print(f"  {c:22s} ratio {v['ratio']:.3f} [{v['ratio_ci'][0]:.3f}, {v['ratio_ci'][1]:.3f}]  {v['class']}")
    v = res["regime"]["contrast"]["dur5|shallow - dur5|typed_paper|D0"]
    print(f"regime: {res['regime']['events']} events, {res['regime']['blocks']} blocks; five-round cap shallow/typed (D0) "
          f"{v['ratio']:.3f} [{v['ratio_ci'][0]:.3f}, {v['ratio_ci'][1]:.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
