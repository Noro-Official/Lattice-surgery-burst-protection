"""Original cohorts behind Fig. 2a, Fig. 5a and Fig. 5b,c (31 October 2026 publication package).

Re-implements the loader and estimator of `experiments/pacer_p1/{campaign.py, analyze.py, analyze_precision.py}` from
`publication_package_20261031_revised`, with the original bootstrap sizes and seeds, so that the archived
`analysis.json` files are reproduced exactly from the per-event block files.

Cohorts:
    precision_d9   d = 9, 3,200 impacts, 4,096 shots; shallow, deep, typed; D0, D0C, F, FC     (Fig. 2a)
    precision_d13  d = 13, 1,200 impacts, 4,096 shots; same layouts and decoders                (Fig. 5a)
    regime         d = 9, 600 impacts, 1,024 shots; six layouts x 12 regimes; D0, F             (Fig. 5b,c)

Estimator:
    per-event risk      1 - hist[0] / shots, where hist[0] counts shots with every logical output correct
    stratified mean     sum_h pi_h * mean over events in stratum h, pi = pool counts / total
    excess risk         worse input mode of (mean - burst-free baseline)
    bootstrap           4,000 within-stratum event resamples + Dirichlet(pool counts + 1) weights, seed 20261029
    precision studies   baseline uncertainty from Jeffreys Beta draws (seed 20261030 + d) on quiet_d{d}_p2.json
    regime study        baseline point values from quiet_d9.json, as `analyze.generic` in the original
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from lsburst._io import DATA_DIR, load_json, read_checksum_list

ARCHIVE_DIR = DATA_DIR / "p1_archive"
PACKAGE_PREFIX = "experiments/pacer_p1/results/"
N_BOOT = 4000
DECODERS = ("D0", "D0C", "F", "FC")
EQUIVALENCE_MARGIN = {9: 1.10, 13: 1.33}


def _ci(x) -> list[float]:
    return [float(v) for v in np.nanpercentile(x, [2.5, 97.5])]


def load_cohort(name: str, data_dir: Path, package_sums: dict[str, str]) -> dict:
    """Loads one cohort's per-event risks, checking every block against the package checksums.

    Args:
        name: Cohort folder name (`precision_d9`, `precision_d13` or `regime`).
        data_dir: The `data/p1_archive` folder.
        package_sums: The package's checksum list, from `PACKAGE_SHA256SUMS.txt`.

    Returns:
        A dict with the study definition `S`, event bank `EV`, per-event risks `H` keyed by
        (condition, decoder, input mode), the number of events `n`, a `complete` flag and the block count.

    Raises:
        ValueError: If a block's bytes or the number of blocks differ from the package.
    """
    folder = data_dir / name
    study = load_json(folder / "study.json")
    events = dict(np.load(data_dir / "events" / study["events_file"]))
    n = len(events["stratum"])
    with gzip.open(folder / "blocks.jsonl.gz", "rt", encoding="utf-8") as f:
        records = [json.loads(line) for line in f]

    risk: dict[tuple[str, str, str], np.ndarray] = {}
    filled: dict[tuple[str, str, str], int] = {}
    for rec in sorted(records, key=lambda r: r["file"]):     # the original loader reads blocks in sorted order
        raw = rec["content"].encode("utf-8")
        digest = hashlib.sha256(raw).hexdigest()
        if digest != rec["sha256"] or package_sums.get(f"{PACKAGE_PREFIX}{name}/blocks/{rec['file']}") != digest:
            raise ValueError(f"block {name}/{rec['file']} does not match the package checksum")
        for row in json.loads(raw)["rows"]:
            for decoder, hist in row["hist"].items():
                key = (row["cond"], decoder, row["mode"])
                risk.setdefault(key, np.full(n, np.nan))[row["e"]] = 1.0 - hist[0] / row["shots"]
                filled[key] = filled.get(key, 0) + 1

    expected = sum(1 for f in package_sums if f.startswith(f"{PACKAGE_PREFIX}{name}/blocks/"))
    if len(records) != expected:
        raise ValueError(f"{name}: {len(records)} blocks present, the package lists {expected}")
    return {"S": study, "EV": events, "H": risk, "n": n, "blocks": len(records),
            "complete": all(v == n for v in filled.values())}


class Estimator:
    """Stratified, event-paired bootstrap; the worse input mode is taken inside every replicate."""

    def __init__(self, cohort: dict, n_boot: int = N_BOOT, seed: int = 20261029):
        strata = cohort["EV"]["stratum"]
        pool = cohort["EV"]["pool_counts"].astype(float)
        self.pi = pool / pool.sum()
        rng = np.random.default_rng(seed)
        self.n_strata = len(pool)
        self.idx = [np.nonzero(strata == h)[0] for h in range(self.n_strata)]
        self.boot_idx = [ii[rng.integers(0, len(ii), size=(n_boot, len(ii)))] for ii in self.idx]
        self.boot_weights = rng.dirichlet(pool + 1.0, size=n_boot)
        self.risk_by_event = cohort["H"]

    def mean(self, cond: str, decoder: str, mode: str, boot: bool = False):
        v = self.risk_by_event[(cond, decoder, mode)]
        if boot:
            return sum(self.boot_weights[:, h] * v[self.boot_idx[h]].mean(1) for h in range(self.n_strata))
        return float(sum(self.pi[h] * v[self.idx[h]].mean() for h in range(self.n_strata)))

    def risk(self, cond: str, decoder: str, boot: bool = False, quiet: dict | None = None):
        """Excess risk (risk if `quiet` is None) for the worse input mode."""
        per_mode = [self.mean(cond, decoder, m, boot) - (0.0 if quiet is None else quiet[(decoder, m)])
                    for m in ("Z", "X")]
        return np.max(per_mode, axis=0) if boot else max(per_mode)


def _quiet_counts(path: Path) -> dict[tuple[int, str, str], tuple[int, int]]:
    """{(T_M, mode, decoder): (failures, shots)}; F shares D0's baseline and FC shares D0C's."""
    out = {}
    for key, v in load_json(path).items():
        tm, mode, decoder = key.split("|")
        for d in (decoder, {"D0": "F", "D0C": "FC"}[decoder]):
            out[(int(tm), mode, d)] = (v["fail"], v["shots"])
    return out


def _classify(diff_ci: list[float], ratio_ci: list[float], margin: float) -> str:
    if diff_ci[0] > 0:
        return "typed better"
    if diff_ci[1] < 0:
        return "comparator better"
    if 1 / margin <= ratio_ci[0] and ratio_ci[1] <= margin:
        return f"equivalent within {margin:g}"
    return "unresolved"


def analyze_precision(name: str, data_dir: Path, package_sums: dict[str, str]) -> dict:
    """Deep and shallow versus typed, as `analyze_precision.study` (Fig. 2a, Fig. 5a)."""
    cohort = load_cohort(name, data_dir, package_sums)
    d, tm = cohort["S"]["d"], cohort["S"]["conditions"][0]["tm"]
    est = Estimator(cohort)
    rng = np.random.default_rng(20261030 + d)
    counts = _quiet_counts(data_dir / f"quiet_d{d}_p2.json")
    q_point, q_boot = {}, {}
    for decoder in DECODERS:
        for mode in ("Z", "X"):
            fail, shots = counts[(tm, mode, decoder)]
            q_point[(decoder, mode)] = fail / shots
            q_boot[(decoder, mode)] = rng.beta(fail + 0.5, shots - fail + 0.5, size=N_BOOT)

    out = {"events": cohort["n"], "shots": cohort["S"]["shots"], "blocks": cohort["blocks"],
           "complete": cohort["complete"], "comparisons": {}}
    for decoder in DECODERS:
        typed = est.risk("typed_paper", decoder, quiet=q_point)
        typed_boot = est.risk("typed_paper", decoder, True, q_boot)
        for arm in ("shallow", "deep"):
            other = est.risk(arm, decoder, quiet=q_point)
            other_boot = est.risk(arm, decoder, True, q_boot)
            diff_ci, ratio_ci = _ci(other_boot - typed_boot), _ci(other_boot / typed_boot)
            out["comparisons"][f"{arm} vs typed|{decoder}"] = {
                "diff": other - typed, "diff_ci": diff_ci,
                "ratio": other / typed if typed > 0 else None, "ratio_ci": ratio_ci,
                "class": _classify(diff_ci, ratio_ci, EQUIVALENCE_MARGIN[d])}
    return out


def analyze_regime(data_dir: Path, package_sums: dict[str, str], keys: list[str]) -> dict:
    """Contrasts `a - b|decoder` with point baselines, as `analyze.generic` (Fig. 5b,c).

    Args:
        data_dir: The `data/p1_archive` folder.
        package_sums: The package's checksum list.
        keys: Contrast keys of the form "a - b|decoder"; the archived analysis lists all 144.
    """
    cohort = load_cohort("regime", data_dir, package_sums)
    est = Estimator(cohort)
    counts = _quiet_counts(data_dir / "quiet_d9.json")
    tm_of = {c["name"]: c["tm"] for c in cohort["S"]["conditions"]}
    out = {"events": cohort["n"], "shots": cohort["S"]["shots"], "blocks": cohort["blocks"],
           "complete": cohort["complete"], "contrast": {}}
    for key in keys:
        pair, decoder = key.rsplit("|", 1)
        a, b = pair.split(" - ")
        q = {(decoder, m): counts[(tm_of[a], m, decoder)][0] / counts[(tm_of[a], m, decoder)][1] for m in ("Z", "X")}
        ra, rb = est.risk(a, decoder, quiet=q), est.risk(b, decoder, quiet=q)
        ba, bb = est.risk(a, decoder, True, q), est.risk(b, decoder, True, q)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = ba / bb
        out["contrast"][key] = {"a": ra, "b": rb, "diff": ra - rb,
                                "diff_ci": np.percentile(ba - bb, [2.5, 97.5]).tolist(),
                                "ratio": ra / rb if rb > 0 else float("inf"),
                                "ratio_ci": np.nanpercentile(ratio, [2.5, 97.5]).tolist()}
    return out


def analyze(data_dir: Path = ARCHIVE_DIR) -> dict:
    """Recomputes all three cohorts. Keys: `precision_d9`, `precision_d13`, `regime`."""
    sums = read_checksum_list(data_dir / "PACKAGE_SHA256SUMS.txt")
    regime_keys = list(load_json(data_dir / "regime" / "analysis.json")["contrast"])
    return {"precision_d9": analyze_precision("precision_d9", data_dir, sums),
            "precision_d13": analyze_precision("precision_d13", data_dir, sums),
            "regime": analyze_regime(data_dir, sums, regime_keys)}


def stored_results(data_dir: Path = ARCHIVE_DIR) -> dict:
    """The archived analyses, restricted to the fields `analyze` recomputes."""
    precision = load_json(data_dir / "precision_analysis.json")
    keys = ("diff", "diff_ci", "ratio", "ratio_ci", "class")
    out = {s: {"comparisons": {c: {k: v[k] for k in keys} for c, v in precision[s]["comparisons"].items()}}
           for s in ("precision_d9", "precision_d13")}
    out["regime"] = {"contrast": load_json(data_dir / "regime" / "analysis.json")["contrast"]}
    return out
