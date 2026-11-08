"""PACER-P1: one self-contained, reproducible burst/lattice-surgery CNOT experiment (2026-10-29).

Everything a rerun needs is declared here or in `config.json`:
  circuit      `simulation/src/analysis/pacer_cnot_surgery.cnot_circuit` (L-shaped A/C/T CNOT, unchanged), with the
               schedule PRE -> ZZ merge (T_M rounds) -> GAP -> XX merge (T_M) -> A measured -> POST.  T_M is a
               parameter (default d) so the merge duration can be varied independently of the distance.
  noise        SI1000 background p = 1e-3 (unchanged injector) + round-resolved burst channels from the frozen bridge
               (`channels_from_site_design`) + state-conditioned readout flips and reset excess from `ReadoutAdapter`
               (computed from `ReadoutParameters()` defaults; no historical cache is used or claimed).
  events       NEW contract (the historical `pacer_vulnerability_surgery.mapping` is unavailable and is NOT
               reconstructed): impacts uniform over the site bounding box + 1.15 mm on the square chip of
               `pacer_cnot_workload.geometry`, energies log-uniform on [0.5, 10] keV (the 0.01 x config range used by
               `pacer_lattice_surgery.prepare_surgery`), delivery by the rectangular phonon solver, flip-chip loss
               kappa = 0.5 with the per-event reference total, onset uniform on the declared window.
  masks        memory boundary rule at a weak fraction, type-blind bands, the typed band, universal protection.
  decoders     D0 (background-only MWPM), F (event-informed MWPM, symmetric report model), FC (F with PyMatching
               correlated matching), D0C (D0 with correlated matching).  F/FC receive the realised event and are
               oracle-informed interfaces, never deployable decoders.

Statistical unit: the impact event.  All masks share each event and the trajectory seed (common random numbers).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from dataclasses import replace
from pathlib import Path

for _n in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_n, "1")

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SIM = ROOT / "simulation"
if str(SIM) not in sys.path:
    sys.path.insert(0, str(SIM))

from scripts import pacer_cnot_workload as cw  # noqa: E402
from src.analysis import pacer_cnot_surgery as cn, pacer_lattice_surgery as ls  # noqa: E402

CONFIG = json.loads((HERE / "config.json").read_text())
CACHE = HERE / "cache"
P_BG = 1e-3
STRONG_GAP_HZ, WEAK_GAP_HZ = 12.0e9, 5.0e9
KAPPA = 0.5
IMPACT_MARGIN_MM = 1.15
_C: dict = {}


# ----------------------------------------------------------------------------------------------- schedule / circuits
def schedule(d: int, t_merge: int | None = None) -> dict:
    s = CONFIG["schedule"]
    tm = d if t_merge is None else int(t_merge)
    return {"pre": int(s["pre"]), "t_merge": tm, "gap": int(s.get("gap", tm) if s.get("gap") is not None else tm), "post": int(s["post"])}


def phase_names(d: int, sch: dict) -> list:
    return cn.phases(d, sch["pre"], sch["t_merge"], sch["gap"], sch["t_merge"], sch["post"])


def merge_windows(sch: dict) -> dict:
    """Syndrome rounds [start, stop) of the two merges."""
    zz0 = sch["pre"]
    xx0 = zz0 + sch["t_merge"] + sch["gap"]
    return {"zz": (zz0, zz0 + sch["t_merge"]), "xx": (xx0, xx0 + sch["t_merge"])}


def prepared(mode: str, d: int, sch: dict) -> dict:
    """`pacer_cnot_workload.prepared` with an explicit merge duration."""
    key = ("prep", mode, d, tuple(sorted(sch.items())))
    if key not in _C:
        from src.analysis import pacer_bridge_v1_qec as qec, pacer_continuous as pc, pacer_surgery_continuous as sx
        nm = phase_names(d, sch)
        cc = cn.cnot_circuit(d, mode, nm)
        tpl = cn.template_cnot(d, mode, nm, P_BG)
        ch, gr, rid = pc.split_rounds(tpl)
        R = len(nm)
        if max(rid) != R:
            raise AssertionError("slot rounds do not end at the final data measurement")
        bg = cn.noisy_cnot(d, mode, nm, P_BG)
        if bg.num_detectors != tpl.num_detectors or bg.num_measurements != tpl.num_measurements:
            raise AssertionError("template and background circuit disagree")
        fixed = ls.matcher_for(bg)
        _C[key] = {"layout": cc.layout, "chunks": ch, "groups": gr, "rid": rid, "R": R, "phases": nm, "records": sx.record_map(tpl, cc.layout),
                   "observables": cc.observables, "detector_phase": cc.detector_phase, "background": bg, "fixed": fixed,
                   "fmask": qec.observable_support_mask(fixed), "template": tpl, "mode": mode, "d": d, "sch": dict(sch)}
    return _C[key]


def fast_tables(P: dict) -> dict:
    """Calibrated-matcher tables (`pacer_fast_calibrated.build`), cached on disk per (mode, d, schedule)."""
    import pickle
    from src.analysis import pacer_fast_calibrated as fc
    sch = P["sch"]
    key = ("fast", P["mode"], P["d"], tuple(sorted(sch.items())))
    if key in _C:
        return _C[key]
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"fast_{P['mode']}_d{P['d']}_pre{sch['pre']}_tm{sch['t_merge']}_gap{sch['gap']}_post{sch['post']}.pkl"
    if f.exists():
        C = pickle.loads(f.read_bytes())
    else:
        clean = cn.template_cnot(P["d"], P["mode"], P["phases"], 0.0)
        C = fc.build(P, clean)
        tmp = f.with_suffix(".tmp")
        tmp.write_bytes(pickle.dumps(C))
        tmp.replace(f)
    _C[key] = C
    return C


# ----------------------------------------------------------------------------------------------- events
def physics() -> dict:
    if "phys" not in _C:
        config, params = ls.load_params()
        used = replace(params, flipchip_loss_hz=float(params.flipchip_loss_hz) * KAPPA)
        _C["phys"] = {"params_ref": params, "pstar": used, "config": config}
    return _C["phys"]


def adapter():
    if "adapter" not in _C:
        from src.analysis import pacer_readout_model as rm
        rp = rm.ReadoutParameters(**CONFIG.get("readout_overrides", {}))
        base = rm.ReadoutAdapter(rp, cache_dir=CACHE)                       # the table every study before 2026-10-29 15:00 used
        ext = rm.ReadoutAdapter(rp, hi_hz=1e9, cache_dir=CACHE)              # extension for rates above 1e7 Hz (narrow footprints)
        n = len(base.grid)
        assert np.allclose(ext.grid[:n], base.grid)
        # hybrid: the original values on the original grid (rebuilt tables differ by Monte Carlo fallback noise), new
        # values only where a rate exceeds the original range; results inside the original range are unchanged
        ext.flip[:n, :n], ext.reset[:n, :n] = base.flip, base.reset
        ext.base_flip, ext.base_reset = base.base_flip, base.base_reset
        _C["adapter"] = ext
    return _C["adapter"]


def mapping(d: int) -> dict:
    """The fields `pacer_surgery_continuous.design` and `pacer_continuous.trajectory_channels` read."""
    G = cw.geometry(d)
    ph = physics()
    return {"layout": G["layout"], "pstar": ph["pstar"], "params": ph["pstar"], "strong": STRONG_GAP_HZ, "weak": WEAK_GAP_HZ}


def impact_box(d: int) -> tuple:
    G = cw.geometry(d)
    lo = G["pos"].min(0) - IMPACT_MARGIN_MM
    hi = G["pos"].max(0) + IMPACT_MARGIN_MM
    W = G["chip"][0]
    return np.maximum(lo, 0.0), np.minimum(hi, W)


def sample_events(d: int, n: int, seed: int, onset_window: tuple, *, footprint_scale: float = 1.0) -> dict:
    """n impacts: position, energy, onset and per-site delivered dose x0 [n, sites] (layout.site_ids order).

    footprint_scale s rescales the phonon diffusivity by s^2 (the delivery length scales as s) and renormalises each
    event to the reference total, i.e. the same deposited quasiparticle number spread over a wider/narrower region."""
    G = cw.geometry(d)
    ph = physics()
    cfg = ph["config"]["design_study"]
    lo, hi = impact_box(d)
    rng = np.random.default_rng(int(seed))
    pos = rng.uniform(lo, hi, size=(int(n), 2))
    e_lo, e_hi = (0.01 * float(v) for v in cfg["energy_keV"])
    energy = np.exp(rng.uniform(math.log(e_lo), math.log(e_hi), size=int(n)))
    onset = rng.integers(int(onset_window[0]), int(onset_window[1]) + 1, size=int(n))
    W = G["chip"][0]
    cell = float(cfg["cell_mm"])
    pref = ph["params_ref"]

    def deliver(p):
        return ls.rect_phonon_delivery(W, W, cell, p, pos, G["pos"], base_loss_hz=p.flipchip_loss_hz) * energy[:, None] * p.generation_scale
    x_ref = deliver(pref)
    used = ph["pstar"]
    if footprint_scale != 1.0:
        used = replace(used, substrate_thickness_mm=float(used.substrate_thickness_mm) * float(footprint_scale) ** 2)
    x0 = deliver(used)
    x0 = x0 * (x_ref.sum(1) / np.maximum(x0.sum(1), 1e-300))[:, None]
    return {"x0": x0, "xy": pos, "energy_keV": energy, "onset": onset, "seed": int(seed), "d": int(d), "footprint_scale": float(footprint_scale)}


def onset_window(sch: dict) -> tuple:
    """Bursts may begin in any round from 0 up to the last round of the XX merge."""
    w = merge_windows(sch)
    return 0, w["xx"][1] - 1


# ----------------------------------------------------------------------------------------------- masks
def eligible(d: int) -> np.ndarray:
    G = cw.geometry(d)
    return (G["role"] != "data") & ~G["seam"]


def memory_counts(d: int, weak_fraction: float) -> dict:
    """Weak ancillas per patch and type: the memory rule's Z : X split (22 : 29 at d = 13) at a weak fraction of a patch's ancillas."""
    per_patch = d * d - 1
    per = int(round(weak_fraction * per_patch))
    z = int(round(22 * per / 51))
    return {"Z": z, "X": per - z}


def mask_memory(d: int, weak_fraction: float) -> np.ndarray:
    return cw.mask_memory(d, counts=memory_counts(d, weak_fraction))


def _relocate(mem: np.ndarray, band: np.ndarray, d: int, what: str) -> np.ndarray:
    G = cw.geometry(d)
    rank = cw.donor_rank(d)
    out = mem.copy()
    for t in ("Z", "X"):
        typ = G["role"] == f"{t}-anc"
        take = np.nonzero(mem & typ & band)[0]
        out[take] = False
        cand = [int(i) for i in np.argsort(rank, kind="stable") if typ[i] and np.isfinite(rank[i]) and not mem[i] and not band[i]]
        if len(cand) < len(take):
            raise ValueError(f"{what} infeasible for {t}: {len(take)} to move, {len(cand)} donors")
        out[cand[:len(take)]] = True
    assert out.sum() == mem.sum()
    return out


def typed_depth(d: int) -> np.ndarray:
    """`pacer_typed_band_rule.typed_depth` for any distance: steps from the new seam checks through checks of the
    merge's measured type that share a data qubit (99 if the patch is not in that merge, -1 for data/seam-only)."""
    from collections import deque
    key = ("td", d)
    if key in _C:
        return _C[key]
    G = cw.geometry(d)
    lay = G["layout"]
    sites = list(lay.site_ids)
    pos = {q: i for i, q in enumerate(sites)}
    st = lay.stabilizers
    out = np.full(len(sites), -1.0)
    for table, typ in (("zz", "Z"), ("xx", "X")):
        chk = {a: set(s) for a, (t, s) in st[table].items() if t == typ}
        byq: dict = {}
        for a, s in chk.items():
            for q in s:
                byq.setdefault(q, []).append(a)
        dist = {a: 0 for a in chk if a not in st["sep"]}
        dq = deque(dist)
        while dq:
            a = dq.popleft()
            for q in chk[a]:
                for b in byq[q]:
                    if b not in dist:
                        dist[b] = dist[a] + 1
                        dq.append(b)
        for q in sites:
            if lay.role[q] == f"{typ}-anc" and q not in lay.seam_only:
                out[pos[q]] = dist.get(q, 99)
    _C[key] = out
    return out


def mask_typed(d: int, weak_fraction: float, w) -> np.ndarray:
    """Typed band B_w of the paper (Eq. main-band): no weak measured-type check with typed depth <= w.
    (`pacer_typed_band_rule` uses depth < K; K = w + 1.)  w may be a pair (w_ZZ, w_XX): one width per merge."""
    td = typed_depth(d)
    wz, wx = (w, w) if np.isscalar(w) else (int(w[0]), int(w[1]))
    role = cw.geometry(d)["role"]
    band = (td >= 0) & (((role == "Z-anc") & (td <= wz)) | ((role == "X-anc") & (td <= wx)))
    return _relocate(mask_memory(d, weak_fraction), band, d, f"typed band w={w}")


def mask_priority(d: int, weak_fraction: float, score: np.ndarray) -> np.ndarray:
    """Calculation-derived allocation: per type, the weak sites are the eligible checks with the LOWEST calculated
    danger score (ties broken by the memory rule's priority, then index); type-specific budgets as the memory rule."""
    G = cw.geometry(d)
    mem = mask_memory(d, weak_fraction)
    rank = cw.donor_rank(d)
    el = eligible(d)
    out = np.zeros(len(mem), bool)
    for t in ("Z", "X"):
        ids = [i for i in np.nonzero(el & (G["role"] == f"{t}-anc"))[0]]
        ids.sort(key=lambda i: (score[i], -rank[i], i))
        out[ids[:int((mem & (G["role"] == f"{t}-anc")).sum())]] = True
    assert out.sum() == mem.sum()
    return out


def mask_blind(d: int, weak_fraction: float, band_lattice_units: int) -> np.ndarray:
    """Type-blind band: no weak ancilla of either type within the band (lattice units) of either seam."""
    sd = cw.seam_distance(d)
    band = (sd["zz"] <= band_lattice_units) | (sd["xx"] <= band_lattice_units)
    return _relocate(mask_memory(d, weak_fraction), band, d, f"type-blind band {band_lattice_units}")


def max_feasible(fn, d: int, weak_fraction: float, values) -> int | None:
    best = None
    for v in values:
        try:
            fn(d, weak_fraction, v)
            best = v
        except ValueError:
            pass
    return best


def universal(d: int) -> np.ndarray:
    return np.zeros(len(cw.geometry(d)["role"]), bool)


def mask_hash(w: np.ndarray) -> str:
    return hashlib.sha256(np.packbits(np.asarray(w, bool)).tobytes()).hexdigest()[:16]


# ----------------------------------------------------------------------------------------------- trajectories
def design(d: int, weak: np.ndarray, protected_gap_hz: float = STRONG_GAP_HZ, trap_multiplier: float = 1.0) -> dict:
    """Per-site gaps and traps.  Weak sites: WEAK_GAP_HZ, base trapping.  Protected ancillas: protected_gap_hz and
    base trapping x trap_multiplier.  Data qubits: STRONG_GAP_HZ always, and (never weak) the protected trap rate."""
    lay = cw.geometry(d)["layout"]
    nd = len(lay.data_ids)
    weak = np.asarray(weak, bool)
    gaps = np.where(weak, WEAK_GAP_HZ, protected_gap_hz).astype(float)
    gaps[:nd] = STRONG_GAP_HZ
    base = physics()["pstar"].trapping_hz
    traps = np.where(weak, base, base * float(trap_multiplier)).astype(float)
    return {"gaps": gaps, "traps": traps}


def channels(P: dict, x0, weak: np.ndarray, onset: int, *, protected_gap_hz: float = STRONG_GAP_HZ,
             energy_scale: float = 1.0, decoupling: str = "standard", trap_multiplier: float = 1.0) -> tuple:
    """(gd, qd, ga, qa, pr, e0, e1), rows 0..R: the burst channels of one event under one mask.
    energy_scale multiplies the delivered dose x0 (delivery is linear in energy; all nonlinearity is downstream),
    relative to the 0.01 x config bank: 1 = the paper's events, 100 = generation-scale (Willow-calibrated) impacts."""
    from src.analysis import pacer_continuous as pc
    des = design(P["d"], weak, protected_gap_hz, trap_multiplier)
    x = np.asarray(x0, float)
    if float(energy_scale) != 1.0:
        x = x * float(energy_scale)
    return pc.trajectory_channels(x, des, physics()["pstar"], P["layout"], adapter(), P["R"], int(onset),
                                  decoupling=decoupling)


FACTORS = ("data", "anc", "reset", "report")


def apply_factors(ch: tuple, keep: dict) -> tuple:
    """Switch off burst-excess components.  keep[f] False removes: data -> data-qubit Pauli excess; anc -> ancilla Pauli
    excess (the propagating channel); reset -> post-reset excess; report -> state-conditioned record flips."""
    gd, qd, ga, qa, pr, e0, e1 = ch
    z = np.zeros_like
    if not keep.get("data", True):
        gd, qd = z(gd), z(qd)
    if not keep.get("anc", True):
        ga, qa = z(ga), z(qa)
    if not keep.get("reset", True):
        pr = z(pr)
    if not keep.get("report", True):
        e0, e1 = z(e0), z(e1)
    return gd, qd, ga, qa, pr, e0, e1


def cap_duration(ch: tuple, onset: int, rounds: int | None) -> tuple:
    """Burst-duration family: every burst excess is removed from round onset + rounds on (None: unchanged)."""
    if rounds is None:
        return ch
    stop = int(onset) + int(rounds)
    out = []
    for a in ch:
        a = a.copy()
        a[stop:] = 0.0
        out.append(a)
    return tuple(out)


def scale_reports(ch: tuple, s: float) -> tuple:
    """Report-strength family: excess record-flip probabilities multiplied by s (capped at 1/2)."""
    gd, qd, ga, qa, pr, e0, e1 = ch
    return gd, qd, ga, qa, pr, np.minimum(e0 * s, 0.5), np.minimum(e1 * s, 0.5)


def run_trajectories(P: dict, C: dict, FT: dict, ch: tuple, shots: int, seed: int, decoders=("D0", "F")) -> dict:
    """Sample one (event, mask) task; returns per decoder the 4-bin histogram of output-error bitmasks."""
    from src.analysis import pacer_fast_calibrated as fc, pacer_fast_runtime as fr, pacer_surgery_continuous as sx
    gd, qd, ga, qa, pr, e0, e1 = ch
    circ = fr.fast_fill(P, C, FT, gd, qd, ga, qa, pr)
    full, obs, _y = sx.sample(circ, shots, seed, P["records"], e0, e1)
    # D0 / F decode the observable-supporting component only (historical convention); the correlated decoders need
    # the conjugate-basis detectors too, since Y-type correlations are their only extra information
    det = full & P["fmask"][None, :] if P["fmask"] is not None else full
    out = {}
    pm = 0.5 * (e0 + e1)
    for name in decoders:
        if name == "D0":
            pred = P["fixed"].decode_batch(det)
        elif name == "D0C":
            pred = background_correlated(P).decode_batch(full, enable_correlations=True)
        elif name == "F":
            pred = fc.matcher(C, gd, qd, ga, qa, pr, pm).decode_batch(det)
        elif name == "FC":
            pred = event_dem_matcher(P, ch, correlated=True).decode_batch(full, enable_correlations=True)
        else:
            raise ValueError(name)
        f = pred != obs
        out[name] = np.bincount(f[:, 0] + 2 * f[:, 1], minlength=4).astype(int).tolist()
    return out


def background_correlated(P: dict):
    key = ("d0c", P["mode"], P["d"], tuple(sorted(P["sch"].items())))
    if key not in _C:
        import pymatching
        dem = P["background"].detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
        _C[key] = pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
    return _C[key]


def event_dem_matcher(P: dict, ch: tuple, correlated: bool = True):
    """Event-informed model built through stim (needed for correlated matching, which reads hyperedge structure)."""
    import pymatching
    from src.analysis import pacer_continuous as pc
    gd, qd, ga, qa, pr, e0, e1 = ch
    _c, model = pc.fill_rounds_pair(P["chunks"], P["groups"], P["rid"], P["layout"], gd, qd, ga, qa, pr, 0.5 * (e0 + e1))
    dem = model.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    return pymatching.Matching.from_detector_error_model(dem, enable_correlations=correlated)


def worker_setup(mode: str, d: int, sch: dict) -> tuple:
    from src.analysis import pacer_fast_runtime as fr
    P = prepared(mode, d, sch)
    C = fast_tables(P)
    key = ("ft", mode, d, tuple(sorted(sch.items())))
    if key not in _C:
        _C[key] = fr.fill_tables(P, C)
    return P, C, _C[key]


def task_seed(base: int, *parts) -> int:
    h = hashlib.sha256(json.dumps([int(base)] + [str(p) for p in parts]).encode()).digest()
    return int.from_bytes(h[:4], "little") & 0x3FFFFFFF
