"""Step 5 runtime control and x100 saturation check (pre-freeze; no layout contrast is computed or printed).

Timing: a 4-event typed control in Z mode at energy_scale 1 and 100, D0C, 1024 shots (seconds per task, measured
inside block(); run in parallel with the saturation blocks, i.e. under campaign-like load).
Saturation: 20 events (every 30th of the 600-impact bank) at x100, standard decoupling, deep and typed masks, both
input modes.  Reported POOLED over the two layouts only: the fraction of tasks with any-output failure >= 0.5 and
the pooled mean failure.  Writes results/step5_timing.json.
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import campaign as cg
import p1core as p
import studies

D, F, TM = 9, 0.3, 9
EVF = str(p.HERE / "results" / "events" / "events_d9_seed78001_n20000_300-180-120_fp1.npz")
OUT = p.HERE / "results" / "step5_timing.json"


def job(mask_idx, mode, items, scale, tag):
    group = {"mask": mask_idx, "tm": TM, "footprint": 1.0, "gap": p.STRONG_GAP_HZ}
    if scale != 1:
        group["energy_scale"] = float(scale)
    conds = [{"name": tag}]
    return ("probe", EVF, D, mode, group, conds, items, 1024, ("D0C",), 78002)


def main(workers: int):
    fam = studies.family(D, F, TM)
    idx = {k: np.nonzero(fam[k])[0].tolist() for k in ("deep", "typed_paper")}
    for mode in ("Z", "X"):
        p.worker_setup(mode, D, p.schedule(D, TM))
    jobs = [("timing_x1", job(idx["typed_paper"], "Z", [0, 1, 2, 3], 1, "typed")),
            ("timing_x100", job(idx["typed_paper"], "Z", [0, 1, 2, 3], 100, "typed"))]
    sel = list(range(0, 600, 30))
    for k in ("deep", "typed_paper"):
        for mode in ("Z", "X"):
            for b in range(0, 20, 4):
                jobs.append((f"sat|{k}|{mode}", job(idx[k], mode, sel[b:b + 4], 100, k)))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        res = list(pool.map(cg.block, [a for _, a in jobs]))
    wall = time.perf_counter() - t0
    out = {"workers": workers, "wall_s": wall}
    for (name, _), r in zip(jobs[:2], res[:2]):
        out[name] = {"sec_per_task": [row["sec"] for row in r["rows"]], "block_seconds": r["seconds"]}
    fails = [1 - row["hist"]["D0C"][0] / row["shots"] for r in res[2:] for row in r["rows"]]
    secs = [row["sec"] for r in res[2:] for row in r["rows"]]
    out["x100_saturation_pooled"] = {"tasks": len(fails), "frac_tasks_fail_ge_0.5": float(np.mean(np.array(fails) >= 0.5)),
                                     "frac_tasks_fail_ge_0.1": float(np.mean(np.array(fails) >= 0.1)),
                                     "pooled_mean_fail": float(np.mean(fails)), "median_sec_per_task": float(np.median(secs))}
    OUT.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 22)
