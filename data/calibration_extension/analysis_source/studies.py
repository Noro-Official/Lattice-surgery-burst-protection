"""Frozen study definitions for PACER-P1 (see PROTOCOL.md).  Seeds are distinct per study; no study reuses another's events.

Mask families at distance d, weak fraction f (default 0.30 of eligible ancillas, memory rule's Z : X split):
  memory      the memory boundary rule
  shallow     type-blind band of round(6 d / 13) lattice units (the paper's two-interface heuristic, scaled)
  deep        the widest feasible type-blind band (even widths)
  typedK      typed band with width K for both merges; typed_paper = round(7 d / 13)
  priority    calculation-ranked rule: weak sites are those with the largest danger score (results/site_scores_*)
  universal   every ancilla protected (different protection count)
"""
from __future__ import annotations

import json

import numpy as np

import p1core as p

F_DEFAULT = 0.30


def shallow_units(d):
    return int(round(6 * d / 13))


def paper_width(d):
    return int(round(7 * d / 13))


def deep_units(d, f):
    return p.max_feasible(p.mask_blind, d, f, range(2, 4 * d, 2))


def site_scores(d: int, tm: int, variant: str = "full", dec: str = "D0") -> np.ndarray:
    S = json.loads((p.HERE / "results" / f"site_scores_d{d}_tm{tm}.json").read_text())
    sc = np.full(len(p.cw.geometry(d)["role"]), np.inf)
    for i, v in S["W"][f"{variant}|{dec}"].items():
        sc[int(i)] = min(v)
    return sc


def danger(d, tm, variant="full", dec="D0"):
    """Higher = more dangerous to leave weak: quiet distance minus the site's mismatch distance."""
    S = json.loads((p.HERE / "results" / f"site_scores_d{d}_tm{tm}.json").read_text())
    q = min(min(S["quiet"][m][dec] for m in ("Z", "X")), 1e9)
    return q - site_scores(d, tm, variant, dec)


def mask_priority_calc(d, f, tm, variant="full", dec="D0"):
    # mask_priority weakens the LOWEST scores: pass the danger itself
    return p.mask_priority(d, f, danger(d, tm, variant, dec))


def family(d: int, f: float, tm: int) -> dict:
    out = {"memory": p.mask_memory(d, f), "shallow": p.mask_blind(d, f, shallow_units(d)), "deep": p.mask_blind(d, f, deep_units(d, f)),
           "typed_paper": p.mask_typed(d, f, paper_width(d)), "universal": p.universal(d)}
    try:
        out["priority"] = mask_priority_calc(d, f, tm)
    except FileNotFoundError:
        pass
    return out


def cond(name, mask_name, mask, tm, **kw):
    return {"name": name, "mask_name": mask_name, "mask": np.nonzero(mask)[0].tolist(), "tm": int(tm), **kw}


def get(name: str) -> dict:
    d = 9
    if name == "width_sweep":
        # Task 2a: typed width (both merges) at three merge durations; prediction frozen in PROTOCOL.md
        conds = []
        for tm in (5, 9, 13):
            fam = family(d, F_DEFAULT, tm)
            for k in ("memory", "shallow", "deep", "universal", "priority"):
                if k in fam:
                    conds.append(cond(f"tm{tm}|{k}", f"tm{tm}|{k}", fam[k], tm))
            wp = paper_width(d)
            # one merge's width varied, the other held at the scaled paper width
            for side, rng_ in (("zz", range(0, 7)), ("xx", range(0, 8))):
                for w in rng_:
                    if side == "xx" and w == wp:
                        continue                                       # identical to zz{wp}
                    pair = (w, wp) if side == "zz" else (wp, w)
                    try:
                        m = p.mask_typed(d, F_DEFAULT, pair)
                    except ValueError:
                        continue
                    conds.append(cond(f"tm{tm}|{side}{w}", f"tm{tm}|typed{pair[0]}-{pair[1]}", m, tm))
        return {"name": name, "d": d, "pool_seed": 75001, "n_pool": 20000, "alloc": (300, 180, 120), "conditions": conds,
                "decoders": ("D0", "F"), "shots": 1024, "seed": 75002}
    if name == "regime":
        # Task 2b-d: report strength, burst duration and footprint scale, each crossed with the main masks (T_M = d)
        tm = d
        conds = []
        fam = family(d, F_DEFAULT, tm)
        main = ("memory", "shallow", "deep", "typed_paper", "priority", "universal")
        for k in main:
            for s in (0.0, 0.25, 0.5, 1.0, 2.0):
                conds.append(cond(f"rs{s:g}|{k}", k, fam[k], tm, report_scale=s))
            for D in (2, 5, 10, 20):
                conds.append(cond(f"dur{D}|{k}", k, fam[k], tm, duration=D))
            for fp in (0.7, 1.4, 2.0):
                conds.append(cond(f"fp{fp:g}|{k}", f"fp{fp:g}|{k}", fam[k], tm, footprint=fp))
        return {"name": name, "d": d, "pool_seed": 76001, "n_pool": 20000, "alloc": (300, 180, 120), "footprints": (1.0, 0.7, 1.4, 2.0),
                "conditions": conds, "decoders": ("D0", "F"), "shots": 1024, "seed": 76002}
    if name == "factorial":
        # Task 4: 2^3 factorial over report flips (R), ancilla Pauli propagation (A) and reset excess (P); data excess on.
        tm = d
        fam = family(d, F_DEFAULT, tm)
        conds = []
        for k in ("shallow", "deep", "typed_paper", "priority", "memory"):
            for R in (0, 1):
                for A in (0, 1):
                    for Pr in (0, 1):
                        conds.append(cond(f"R{R}A{A}P{Pr}|{k}", k, fam[k], tm, keep={"report": bool(R), "anc": bool(A), "reset": bool(Pr), "data": True}))
            conds.append(cond(f"R0A0P0D0|{k}", k, fam[k], tm, keep={"report": False, "anc": False, "reset": False, "data": False}))
        return {"name": name, "d": d, "pool_seed": 77001, "n_pool": 20000, "alloc": (300, 180, 120), "conditions": conds,
                "decoders": ("D0", "F"), "shots": 1024, "seed": 77002}
    if name == "baselines":
        # Tasks 3 + 5: layout x decoder on identical events, three budgets; FC/D0C are correlated matching
        tm = d
        conds = []
        for f in (0.15, 0.30, 0.45):
            fam = family(d, f, tm)
            for k, m in fam.items():
                if k == "universal" and f != 0.30:
                    continue
                conds.append(cond(f"f{f:g}|{k}", f"f{f:g}|{k}", m, tm))
        return {"name": name, "d": d, "pool_seed": 78001, "n_pool": 20000, "alloc": (300, 180, 120), "conditions": conds,
                "decoders": ("D0", "D0C", "F", "FC"), "shots": 1024, "seed": 78002}
    if name.startswith("exchange_test"):
        X = json.loads((p.HERE / "results" / "exchange_predictions.json").read_text())
        conds = [cond(c["name"], c["name"], _from_idx(d, c["weak"]), X["tm"]) for c in X["layouts"]]
        return {"name": name, "d": d, "pool_seed": 79101, "n_pool": 20000, "alloc": (300, 180, 120), "conditions": conds,
                "decoders": ("D0", "F"), "shots": 1024, "seed": 79102}
    if name == "d13_confirm":
        d13 = 13
        tm = 13
        fam = family(d13, F_DEFAULT, tm)
        conds = [cond(k, k, fam[k], tm) for k in ("memory", "shallow", "deep", "typed_paper", "priority", "universal") if k in fam]
        for w in (4, 5, 6, 8):
            conds.append(cond(f"typed{w}", f"typed{w}", p.mask_typed(d13, F_DEFAULT, w), tm))
        return {"name": name, "d": d13, "pool_seed": 80001, "n_pool": 20000, "alloc": (240, 120, 60), "conditions": conds,
                "decoders": ("D0", "F"), "shots": 1024, "seed": 80002, "block_events": 2}
    if name in ("precision_d9", "precision_d13"):
        # Sec. P2 precision replication: the three central layouts at f = 0.30 on fresh events, five pre-trajectory score
        # strata with allocation and shots fixed from the pilot variance (vr_size.py), all four decoders
        dd = 9 if name == "precision_d9" else 13
        fam = family(dd, F_DEFAULT, dd)
        arms = ("shallow", "deep", "typed_paper")
        conds = [cond(k, k, fam[k], dd) for k in arms]
        spec = {"precision_d9": dict(pool_seed=85001, n_pool=60000, alloc=(352, 897, 289, 713, 949), seed=85002, block_events=4),
                "precision_d13": dict(pool_seed=86001, n_pool=40000, alloc=(73, 120, 101, 562, 344), seed=86002, block_events=2)}[name]
        return {"name": name, "d": dd, **spec, "conditions": conds, "decoders": ("D0", "D0C", "F", "FC"), "shots": 4096,
                "score_arms": {k: np.nonzero(fam[k])[0].tolist() for k in arms}, "score_typed": "typed_paper", "quiet_tag": "p2"}
    raise KeyError(name)


def _from_idx(d, idx):
    w = np.zeros(len(p.cw.geometry(d)["role"]), bool)
    w[list(idx)] = True
    return w
