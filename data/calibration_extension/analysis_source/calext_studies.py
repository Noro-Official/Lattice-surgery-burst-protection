"""Calibration-extension study grid (frozen in FROZEN_PROTOCOL.json; do not edit after freezing).

Physics setting = (energy_scale, decoupling, trap_multiplier, protected gap).  Every setting runs the three Fig. 2a
layouts (shallow / deep / typed, d = 9, weak fraction 0.3, T_M = 9) on the 78001 bank, D0C, 1024 shots, trajectory
seed base 78002 (all settings share trajectory seeds per (event, input mode): common random numbers).
"""
from __future__ import annotations

import numpy as np

import p1core as p
import studies

D, F, TM = 9, 0.3, 9
MASKS = ("shallow", "deep", "typed_paper")
EXPECTED_HASH = {"shallow": "d8ba0ac6af1b5bce", "deep": "1603afb01c8bb608", "typed_paper": "feb78d5b7189d641"}
JJONLY_GAP_HZ = 10.0e9

# (wave, energy_scale, decoupling, trap_multiplier, protected_gap_hz)
SETTINGS = {
    "A0": [(1, "standard", 1, p.STRONG_GAP_HZ)],                                         # anchor (regression)
    "A": [(1, "echo", 1, p.STRONG_GAP_HZ), (10, "standard", 1, p.STRONG_GAP_HZ), (10, "echo", 1, p.STRONG_GAP_HZ),
          (100, "standard", 1, p.STRONG_GAP_HZ), (100, "echo", 1, p.STRONG_GAP_HZ)],
    "B": [(3, "standard", 1, p.STRONG_GAP_HZ), (3, "echo", 1, p.STRONG_GAP_HZ),
          (30, "standard", 1, p.STRONG_GAP_HZ), (30, "echo", 1, p.STRONG_GAP_HZ)],
    "C": [(1, "echo", 10, p.STRONG_GAP_HZ), (100, "echo", 10, p.STRONG_GAP_HZ), (1, "echo", 1, JJONLY_GAP_HZ)],
}


def setting_tag(es, dec, tr, gap) -> str:
    return f"E{es:g}|{dec}|T{tr:g}|G{gap / 1e9:g}"


def cond_name(mask, es, dec, tr, gap) -> str:
    return f"{mask}|{setting_tag(es, dec, tr, gap)}"


def get(wave: str) -> dict:
    fam = studies.family(D, F, TM)
    for k in MASKS:
        h = p.mask_hash(fam[k])
        assert h == EXPECTED_HASH[k], (k, h)
    conds = []
    for es, dec, tr, gap in SETTINGS[wave]:
        for k in MASKS:
            c = {"name": cond_name(k, es, dec, tr, gap), "mask_name": k, "mask": np.nonzero(fam[k])[0].tolist(), "tm": TM}
            if gap != p.STRONG_GAP_HZ:
                c["gap"] = float(gap)
            if es != 1:
                c["energy_scale"] = float(es)
            if dec != "standard":
                c["decoupling"] = dec
            if tr != 1:
                c["trap_multiplier"] = float(tr)
            conds.append(c)
    return {"name": f"calext_{wave}", "d": D, "pool_seed": 78001, "n_pool": 20000, "alloc": [300, 180, 120],
            "decoders": ["D0C"], "shots": 1024, "seed": 78002, "block_events": 4, "conditions": conds}


if __name__ == "__main__":
    import sys

    import campaign as cg
    S = get(sys.argv[1])
    raise SystemExit(cg.run(S, int(sys.argv[2]) if len(sys.argv) > 2 else 20))
