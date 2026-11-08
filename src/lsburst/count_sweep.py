"""Refined protection-count sweep at d = 9 (Table I, Fig. 2b).

The typed layout protects the n highest-scoring ancillas for n in {117, 123, 132, 150, 168}; the shallow layout
protects 168. For each count this module estimates the typed layout's burst-attributable excess risk under correlated
background matching (D0C), then finds the count at which typed reaches the shallow layout's risk by log-linear
interpolation between swept counts.

Estimator (identical to `data/count_sweep/recompute_results.py`, the export's own check):
    per-event risk    p_e = (shots with any logical output wrong) / shots
    stratified mean   m(c, mode) = sum_h pi_h * mean_{e in h} p_e,  pi = pool counts / total
    excess risk       R(c) = max over input modes of [m(c, mode) - q(mode)],  q = burst-free D0C failure rate
    bootstrap         4,000 within-stratum event resamples and Dirichlet(pool counts + 1) weights (seed 20261106);
                      Jeffreys Beta(k + 1/2, n - k + 1/2) quiet-baseline draws (seed 20261107)
    interval          2.5th and 97.5th percentiles of the per-replicate crossing count
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import numpy as np

from lsburst._io import DATA_DIR, load_json

STUDY_DIR = DATA_DIR / "count_sweep"
N_BOOT = 4000
SEED = 20261106
SHALLOW_COUNT = 168
STRATA = 3


def _per_event_risk(path: Path, n_events: int) -> dict[tuple[str, str], np.ndarray]:
    """Reads per_event_counts.csv.gz into {(condition, input mode): per-event failure fraction}."""
    risk: dict[tuple[str, str], np.ndarray] = {}
    with gzip.open(path, "rt") as f:
        for row in csv.DictReader(f):
            key = (row["condition"], row["input_mode"])
            risk.setdefault(key, np.empty(n_events))[int(row["event_index"])] = (
                int(row["any_output_failures"]) / int(row["shots"]))
    return risk


def _crossing(typed_risk: np.ndarray, target: float, counts: np.ndarray) -> float:
    """Count at which typed risk first falls to `target`, by interpolation in log risk.

    Returns the lowest swept count if typed is already below the target there, and the highest if it never gets there.
    """
    log_risk = np.log(np.maximum(typed_risk, 1e-7))
    log_target = np.log(max(target, 1e-7))
    if log_risk[0] <= log_target:
        return float(counts[0])
    for k in range(1, len(counts)):
        if log_risk[k] <= log_target:
            frac = (log_risk[k - 1] - log_target) / (log_risk[k - 1] - log_risk[k])
            return float(counts[k - 1] + (counts[k] - counts[k - 1]) * frac)
    return float(counts[-1])


def analyze(data_dir: Path = STUDY_DIR) -> dict:
    """Recomputes the count sweep from its per-event counts.

    Args:
        data_dir: Folder holding the count-sweep export.

    Returns:
        `typed_excess` (one value per swept count), `n_star` (required protected count), `n_star_ci`, `saving`
        (168 - n_star), `saving_ci` and `lower_edge_fraction` (share of replicates at the lowest swept count).
    """
    stored = load_json(data_dir / "RESULTS.json")
    study = load_json(data_dir / "study.json")
    events = dict(np.load(data_dir / study["events_file"]))
    quiet = load_json(data_dir / "quiet_d9.json")
    risk = _per_event_risk(data_dir / "per_event_counts.csv.gz", len(events["stratum"]))

    # Bootstrap draws, in the same order as the original script so the random streams match.
    rng = np.random.default_rng(SEED)
    pool = events["pool_counts"].astype(float)
    pi = pool / pool.sum()
    idx = [np.flatnonzero(events["stratum"] == h) for h in range(STRATA)]
    picks = [ii[rng.integers(0, len(ii), size=(N_BOOT, len(ii)))] for ii in idx]
    weights = rng.dirichlet(pool + 1, size=N_BOOT)
    rng = np.random.default_rng(SEED + 1)
    quiet_draws = {}
    for base in ("D0", "D0C"):
        for mode in ("Z", "X"):
            q = quiet[f"9|{mode}|{base}"]
            quiet_draws[base, mode] = rng.beta(q["fail"] + 0.5, q["shots"] - q["fail"] + 0.5, size=N_BOOT)

    def excess(condition: str, boot: bool = False):
        per_mode = []
        for mode in ("Z", "X"):
            x = risk[condition, mode]
            if boot:
                mean = sum(weights[:, h] * x[picks[h]].mean(1) for h in range(STRATA))
                per_mode.append(mean - quiet_draws["D0C", mode])
            else:
                mean = sum(pi[h] * x[idx[h]].mean() for h in range(STRATA))
                per_mode.append(mean - quiet[f"9|{mode}|D0C"]["p"])
        return np.max(per_mode, axis=0) if boot else max(per_mode)

    counts = np.array(stored["protected_counts"])
    typed = np.array([excess(f"n{n}|typed") for n in counts])
    typed_boot = np.array([excess(f"n{n}|typed", boot=True) for n in counts])
    target = excess(f"n{SHALLOW_COUNT}|shallow")
    target_boot = excess(f"n{SHALLOW_COUNT}|shallow", boot=True)

    n_star = _crossing(typed, target, counts)
    n_star_boot = np.array([_crossing(typed_boot[:, j], target_boot[j], counts) for j in range(N_BOOT)])
    return {
        "typed_excess": typed.tolist(),
        "shallow_excess": float(target),
        "n_star": n_star,
        "n_star_ci": np.percentile(n_star_boot, [2.5, 97.5]).tolist(),
        "saving": SHALLOW_COUNT - n_star,
        "saving_ci": np.percentile(SHALLOW_COUNT - n_star_boot, [2.5, 97.5]).tolist(),
        "lower_edge_fraction": float(np.mean(n_star_boot == counts[0])),
    }


def stored_results(data_dir: Path = STUDY_DIR) -> dict:
    """The export's stored values for the same quantities (correlated background matching, D0C)."""
    d = load_json(data_dir / "RESULTS.json")["decoders"]["D0C"]
    return {"typed_excess": d["typed_excess"], "n_star": d["n_star"], "n_star_ci": d["n_star_ci"],
            "saving_ci": d["saving_ci"]}
