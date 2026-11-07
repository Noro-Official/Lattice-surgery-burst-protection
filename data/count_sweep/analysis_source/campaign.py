"""Generic P1 campaign runner: event sets, condition grids, parallel trajectory blocks, stratified paired analysis.

A study (see `studies.py`) declares: distance, event-set spec (pool seed, pool size, stratum allocation), conditions
and decoders.  A condition is a dict
    {"name", "mask": index list (weak sites), "tm": merge duration, "footprint": scale, "gap": protected gap (Hz),
     "keep": factor switches, "report_scale": s, "duration": rounds | None}
Conditions sharing (mask, tm, footprint, gap) share the burst channels of each event; every condition of an event and
input mode uses the SAME trajectory seed (common random numbers).

Event population: a pool of impacts (uniform position, log-uniform energy, onset uniform on the schedule's window via
a stored uniform variate).  Strata by pool rank of total delivered dose: H1 top 5 %, H2 next 15 %, H3 the rest; a
fixed number of events is drawn from each.  The population mean is sum_h pi_h * mean_h with pi_h the pool fractions.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

import p1core as p

RES = p.HERE / "results"
STRATA_Q = (0.95, 0.80)                                                     # dose quantiles separating H1 | H2 | H3


# ----------------------------------------------------------------------------------------------- events
def build_events(d: int, pool_seed: int, n_pool: int, alloc: tuple, footprints=(1.0,), chunk: int = 500) -> Path:
    f = RES / "events" / f"events_d{d}_seed{pool_seed}_n{n_pool}_{'-'.join(map(str, alloc))}_fp{'-'.join(f'{x:g}' for x in footprints)}.npz"
    if f.exists():
        return f
    f.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(pool_seed)
    lo, hi = p.impact_box(d)
    pos = rng.uniform(lo, hi, size=(n_pool, 2))
    ds = p.physics()["config"]["design_study"]
    e_lo, e_hi = (0.01 * float(v) for v in ds["energy_keV"])
    energy = np.exp(rng.uniform(np.log(e_lo), np.log(e_hi), size=n_pool))
    u_onset = rng.random(n_pool)
    dose = np.concatenate([_deliver(d, pos[k:k + chunk], energy[k:k + chunk], 1.0).sum(1) for k in range(0, n_pool, chunk)])
    q1, q2 = np.quantile(dose, STRATA_Q[0]), np.quantile(dose, STRATA_Q[1])
    stratum = np.where(dose >= q1, 0, np.where(dose >= q2, 1, 2))
    counts = np.bincount(stratum, minlength=3)
    sel_rng = np.random.default_rng(pool_seed + 1)
    sel = np.concatenate([np.sort(sel_rng.choice(np.nonzero(stratum == h)[0], size=int(alloc[h]), replace=False)) for h in range(3)])
    out = {"xy": pos[sel], "energy_keV": energy[sel], "u_onset": u_onset[sel], "stratum": stratum[sel], "pool_counts": counts, "pool_index": sel,
           "dose": dose[sel], "footprints": np.array(footprints, float)}
    for fp in footprints:
        out[f"x0_fp{fp:g}"] = np.concatenate([_deliver(d, pos[sel[k:k + chunk]], energy[sel[k:k + chunk]], fp) for k in range(0, len(sel), chunk)])
    np.savez_compressed(f, **out)
    return f


def build_events_scored(d: int, pool_seed: int, n_pool: int, alloc: tuple, arms: dict, typed: str, chunk: int = 500) -> Path:
    """Precision-study event set: five pre-trajectory score strata on the pool (fixed in PROTOCOL.md, Sec. P2).
    S0 / S1: top 1 % / next 4 % of the exposure of sites weak in some but not all `arms` (where the layouts disagree);
    S2: of the rest, the top 5 % (pool fraction) of exposure of the `typed` layout's weak sites; S3: of the rest, the
    next 15 % by total dose; S4: the remainder.  Unbiased for any score: the estimator weights strata by pool counts."""
    f = RES / "events" / f"events_d{d}_seed{pool_seed}_n{n_pool}_scored_{'-'.join(map(str, alloc))}.npz"
    if f.exists():
        return f
    f.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(pool_seed)
    lo, hi = p.impact_box(d)
    pos = rng.uniform(lo, hi, size=(n_pool, 2))
    ds = p.physics()["config"]["design_study"]
    e_lo, e_hi = (0.01 * float(v) for v in ds["energy_keV"])
    energy = np.exp(rng.uniform(np.log(e_lo), np.log(e_hi), size=n_pool))
    u_onset = rng.random(n_pool)
    nsite = len(p.cw.geometry(d)["role"])
    W = {}
    for a, idx in arms.items():
        W[a] = np.zeros(nsite, bool)
        W[a][list(idx)] = True
    anyw, allw = np.any(list(W.values()), 0), np.all(list(W.values()), 0)
    dis_set, typ_set = anyw & ~allw, W[typed]
    dose, dis, typ = (np.empty(n_pool) for _ in range(3))
    for k in range(0, n_pool, chunk):
        x = _deliver(d, pos[k:k + chunk], energy[k:k + chunk], 1.0)
        dose[k:k + chunk], dis[k:k + chunk], typ[k:k + chunk] = x.sum(1), x[:, dis_set].sum(1), x[:, typ_set].sum(1)
    rank = lambda v, m: np.argsort(np.argsort(-np.where(m, v, -np.inf)))  # noqa: E731  descending rank within mask
    st = np.full(n_pool, 4)
    r = rank(dis, np.ones(n_pool, bool))
    st[r < round(0.01 * n_pool)] = 0
    st[(r >= round(0.01 * n_pool)) & (r < round(0.05 * n_pool))] = 1
    rest = st == 4
    st[rest & (rank(typ, rest) < round(0.05 * n_pool))] = 2
    rest = st == 4
    st[rest & (rank(dose, rest) < round(0.15 * n_pool))] = 3
    counts = np.bincount(st, minlength=5)
    sel_rng = np.random.default_rng(pool_seed + 1)
    sel = np.concatenate([np.sort(sel_rng.choice(np.nonzero(st == h)[0], size=int(alloc[h]), replace=False)) for h in range(5)])
    out = {"xy": pos[sel], "energy_keV": energy[sel], "u_onset": u_onset[sel], "stratum": st[sel], "pool_counts": counts, "pool_index": sel,
           "dose": dose[sel], "score_disagree": dis[sel], "score_typed": typ[sel], "footprints": np.array([1.0])}
    out["x0_fp1"] = np.concatenate([_deliver(d, pos[sel[k:k + chunk]], energy[sel[k:k + chunk]], 1.0) for k in range(0, len(sel), chunk)])
    np.savez_compressed(f, **out)
    return f


def events_for(S: dict) -> Path:
    if "score_arms" in S:
        return build_events_scored(S["d"], S["pool_seed"], S["n_pool"], tuple(S["alloc"]), S["score_arms"], S["score_typed"])
    return build_events(S["d"], S["pool_seed"], S["n_pool"], tuple(S["alloc"]), tuple(S.get("footprints", (1.0,))))


def _deliver(d, pos, energy, fp):
    """Per-site dose for given impacts (the p1core.sample_events delivery, without resampling)."""
    from dataclasses import replace
    G = p.cw.geometry(d)
    ph = p.physics()
    W = G["chip"][0]
    cell = float(ph["config"]["design_study"]["cell_mm"])
    pref, used = ph["params_ref"], ph["pstar"]
    if fp != 1.0:
        used = replace(used, substrate_thickness_mm=float(used.substrate_thickness_mm) * float(fp) ** 2)
    x_ref = p.ls.rect_phonon_delivery(W, W, cell, pref, pos, G["pos"], base_loss_hz=pref.flipchip_loss_hz) * energy[:, None] * pref.generation_scale
    x0 = p.ls.rect_phonon_delivery(W, W, cell, used, pos, G["pos"], base_loss_hz=used.flipchip_loss_hz) * energy[:, None] * used.generation_scale
    return x0 * (x_ref.sum(1) / np.maximum(x0.sum(1), 1e-300))[:, None]


def onset_of(u: float, d: int, tm: int) -> int:
    lo, hi = p.onset_window(p.schedule(d, tm))
    return int(lo + min(int(u * (hi - lo + 1)), hi - lo))


# ----------------------------------------------------------------------------------------------- worker
_EV: dict = {}


def _events(path):
    if path not in _EV:
        _EV[path] = dict(np.load(path))
    return _EV[path]


def block(args):
    """args = (study, events_file, d, mode, group, conds, items, shots, decoders, seed_base); group = shared-channel key."""
    study, events_file, d, mode, group, conds, items, shots, decoders, seed_base = args
    t0 = time.perf_counter()
    EV = _events(events_file)
    tm, fp, gap = group["tm"], group["footprint"], group["gap"]
    sch = p.schedule(d, tm)
    P, C, FT = p.worker_setup(mode, d, sch)
    weak = np.zeros(len(P["layout"].site_ids), bool)
    weak[list(group["mask"])] = True
    rows = []
    x0s = EV[f"x0_fp{fp:g}"]
    for e in items:
        onset = onset_of(float(EV["u_onset"][e]), d, tm)
        base = p.channels(P, x0s[e], weak, onset, protected_gap_hz=gap)
        seed = p.task_seed(seed_base, e, mode)
        for c in conds:
            ch = base
            if c.get("keep"):
                ch = p.apply_factors(ch, c["keep"])
            if c.get("report_scale", 1.0) != 1.0:
                ch = p.scale_reports(ch, c["report_scale"])
            if c.get("duration") is not None:
                ch = p.cap_duration(ch, onset, c["duration"])
            tt = time.perf_counter()
            res = p.run_trajectories(P, C, FT, ch, shots, seed, decoders=tuple(c.get("decoders", decoders)))
            rows.append({"cond": c["name"], "e": int(e), "mode": mode, "shots": int(shots), "onset": onset, "hist": res, "sec": time.perf_counter() - tt})
    return {"rows": rows, "seconds": time.perf_counter() - t0}


# ----------------------------------------------------------------------------------------------- driver
def group_key(c):
    return json.dumps({"mask": c["mask_name"], "tm": c["tm"], "footprint": c.get("footprint", 1.0), "gap": c.get("gap", p.STRONG_GAP_HZ)}, sort_keys=True)


def make_jobs(S: dict, events_file: Path, n_events: int, block_events: int = 4) -> list:
    groups: dict = {}
    for c in S["conditions"]:
        groups.setdefault(group_key(c), []).append(c)
    jobs = []
    for gk, conds in groups.items():
        c0 = conds[0]
        group = {"mask": c0["mask"], "tm": c0["tm"], "footprint": c0.get("footprint", 1.0), "gap": c0.get("gap", p.STRONG_GAP_HZ)}
        lite = [{k: v for k, v in c.items() if k not in ("mask",)} for c in conds]
        gid = hashlib.sha256(gk.encode()).hexdigest()[:10]
        for mode in ("Z", "X"):
            for k in range(0, n_events, block_events):
                items = list(range(k, min(n_events, k + block_events)))
                name = f"{mode}_{gid}_{k // block_events:04d}"
                jobs.append((name, (S["name"], str(events_file), S["d"], mode, group, lite, items, S["shots"], tuple(S["decoders"]), S["seed"])))
    return jobs


def run(S: dict, workers: int) -> int:
    from concurrent.futures import ProcessPoolExecutor, as_completed
    ef = events_for(S)
    n = int(sum(S["alloc"]))
    folder = RES / S["name"] / "blocks"
    folder.mkdir(parents=True, exist_ok=True)
    (RES / S["name"] / "study.json").write_text(json.dumps({k: v for k, v in S.items() if k not in ("conditions", "score_arms")} | {"conditions": [{k: v for k, v in c.items() if k != "mask"} | {"mask_hash": p.mask_hash(_mk(S["d"], c["mask"]))} for c in S["conditions"]], "events_file": ef.name}, indent=1, default=str))
    jobs = [j for j in make_jobs(S, ef, n, S.get("block_events", 4)) if not (folder / f"{j[0]}.json").exists()]
    # warm caches (fast tables) once per schedule before forking work
    for tm in sorted({c["tm"] for c in S["conditions"]}):
        for mode in ("Z", "X"):
            p.worker_setup(mode, S["d"], p.schedule(S["d"], tm))
    print(f"{S['name']}: {len(jobs)} blocks, {workers} workers", flush=True)
    t0, failed, done = time.perf_counter(), 0, 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(block, a): name for name, a in jobs}
        for f in as_completed(futs):
            name = futs[f]
            try:
                res = f.result()
            except Exception as exc:  # noqa: BLE001
                failed += 1
                print("block failed", name, repr(exc), flush=True)
                continue
            tmp = folder / f"{name}.tmp"
            tmp.write_text(json.dumps(res))
            tmp.replace(folder / f"{name}.json")
            done += 1
            if done % 50 == 0 or done == len(jobs):
                el = time.perf_counter() - t0
                print(f"  {done}/{len(jobs)} blocks, {el / 60:.1f} min, eta {(len(jobs) - done) * el / done / 60:.1f} min", flush=True)
    print(f"{S['name']}: finished, {failed} failed, {(time.perf_counter() - t0) / 60:.1f} min", flush=True)
    return 1 if failed else 0


def _mk(d, idx):
    w = np.zeros(len(p.cw.geometry(d)["role"]), bool)
    w[list(idx)] = True
    return w


# ----------------------------------------------------------------------------------------------- analysis
def load(name: str) -> dict:
    """{(cond, decoder, mode): array over events of P(any output wrong)} plus event metadata."""
    S = json.loads((RES / name / "study.json").read_text())
    EV = dict(np.load(RES / "events" / S["events_file"]))
    n = len(EV["stratum"])
    H: dict = {}
    cnt: dict = {}
    for f in sorted((RES / name / "blocks").glob("*.json")):
        for r in json.loads(f.read_text())["rows"]:
            for dec, h in r["hist"].items():
                k = (r["cond"], dec, r["mode"])
                H.setdefault(k, np.full(n, np.nan))[r["e"]] = 1.0 - h[0] / r["shots"]
                cnt[k] = cnt.get(k, 0) + 1
    return {"S": S, "EV": EV, "H": H, "n": n, "complete": {k: v == n for k, v in cnt.items()}}


class Estimator:
    """Stratified, event-paired bootstrap; worse-basis maximum taken inside every replicate."""

    def __init__(self, D: dict, B: int = 4000, seed: int = 20261029):
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

    def risk(self, cond, dec, boot=False, quiet: dict | None = None, agg: str = "max"):
        """agg "max": worse input mode (R_max); "mean": input-mode average (R_common)."""
        vals = []
        for mode in ("Z", "X"):
            q = 0.0 if quiet is None else quiet[(dec, mode)]
            vals.append(self.mean(cond, dec, mode, boot) - q)
        if agg == "mean":
            return np.mean(vals, axis=0) if boot else float(np.mean(vals))
        return np.max(vals, axis=0) if boot else max(vals)

    def contrast(self, a, b, dec, quiet=None, agg="max"):
        """risk(a) - risk(b) with 95 % percentile interval, and the ratio a / b (excess risks if quiet given)."""
        ra, rb = self.risk(a, dec, quiet=quiet, agg=agg), self.risk(b, dec, quiet=quiet, agg=agg)
        ba, bb = self.risk(a, dec, True, quiet, agg), self.risk(b, dec, True, quiet, agg)
        diff = ba - bb
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = ba / bb
        return {"a": ra, "b": rb, "diff": ra - rb, "diff_ci": np.percentile(diff, [2.5, 97.5]).tolist(),
                "ratio": ra / rb if rb > 0 else float("inf"), "ratio_ci": np.nanpercentile(ratio, [2.5, 97.5]).tolist()}


if __name__ == "__main__":
    import studies
    S = studies.get(sys.argv[1])
    raise SystemExit(run(S, int(sys.argv[2]) if len(sys.argv) > 2 else 26))
