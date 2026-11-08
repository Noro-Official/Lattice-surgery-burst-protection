"""Write FROZEN_PROTOCOL.json for calibration_extension_20261108 (run once, before any campaign block exists)."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np

import calext_studies as cs
import p1core as p

HERE = Path(__file__).resolve().parent
OUT = HERE / "FROZEN_PROTOCOL.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert not OUT.exists(), "protocol already frozen"
    for w in cs.SETTINGS:
        assert not (HERE / "results" / f"calext_{w}" / "blocks").exists(), f"blocks exist for {w}: cannot freeze after running"
    import pymatching
    import stim
    timing = json.loads((HERE / "results" / "step5_timing.json").read_text())
    checks = json.loads((HERE / "results" / "step5_checks.json").read_text())
    sim = p.SIM
    code = {k: sha(v) for k, v in {
        "p1core.py": HERE / "p1core.py", "campaign.py": HERE / "campaign.py", "studies.py": HERE / "studies.py",
        "calext_studies.py": HERE / "calext_studies.py", "analyze_calext.py": HERE / "analyze_calext.py",
        "config.json": HERE / "config.json", "quiet_d9.json": HERE / "quiet_d9.json",
        "simulation/src/analysis/pacer_continuous.py": sim / "src/analysis/pacer_continuous.py",
        "simulation/src/physics/pacer_qp_bridge_v1.py": sim / "src/physics/pacer_qp_bridge_v1.py",
        "simulation/src/physics/nonmarkov.py": sim / "src/physics/nonmarkov.py",
        "simulation/src/analysis/pacer_lattice_surgery.py": sim / "src/analysis/pacer_lattice_surgery.py",
        "simulation/src/analysis/pacer_readout_model.py": sim / "src/analysis/pacer_readout_model.py",
    }.items()}
    ph = p.physics()
    grid = []
    for w, sets in cs.SETTINGS.items():
        for es, dec, tr, gap in sets:
            grid.append({"wave": w, "tag": cs.setting_tag(es, dec, tr, gap), "energy_scale_rel_bank": es,
                         "energy_scale_rel_generation_calibration": es / 100, "decoupling": dec,
                         "trap_multiplier_protected": tr, "protected_ancilla_gap_ghz": gap / 1e9})
    proto = {
        "study": "calibration_extension_20261108",
        "date": "2026-11-08",
        "status": "frozen before any campaign trajectory; only the Step 5 controls (regression, std-vs-echo channel totals, "
                  "peak-shift table, 4-event timing, pooled x100 saturation) were run, none of which computes a layout contrast",
        "question": "Does the typed band's advantage over the deep (and shallow) type-blind band survive (1) impact energies up to the "
                    "generation-scale calibration (Willow median 2 MHz peak shift), (2) echoed instead of standard decoupling, "
                    "(3) M1-style faster trapping on protected sites, (4) a JJ-Only 10 GHz protected gap?",
        "model_changes": {
            "energy_scale": "x0 (delivered dose) multiplied by energy_scale before the bridge; delivery is linear in energy, all "
                            "nonlinearity (recombination, tunnelling) is downstream. Relative to the paper bank (0.01 x config "
                            "energy range); x100 = config range on which generation_scale was solved.",
            "decoupling": "passed to channels_from_site_design (Platform.decoupling); 'echo' uses the frozen power-law echo residual.",
            "trap_multiplier": "trapping_hz of every non-weak site (protected ancillas AND data qubits) = base x multiplier; weak sites base.",
            "protected_gap": "protected ancilla gap; data qubits stay 12 GHz. Hybrid (Willow f_q, mean gap), not an MIT device model.",
            "defaults_bit_identical": checks["regression_bit_identical"],
        },
        "frozen_draw": {"generation_scale": ph["params_ref"].generation_scale, "trapping_hz": ph["pstar"].trapping_hz,
                        "echo_model": ph["pstar"].echo_model, "kappa": p.KAPPA},
        "grid": grid,
        "layouts": {k: {"hash": cs.EXPECTED_HASH[k], "weak_count": checks["mask_weak_counts"][k]} for k in cs.MASKS},
        "layout_definitions": "studies.family(d=9, f=0.3, T_M=9): shallow, deep (widest feasible type-blind band), typed_paper (= count sweep n168|typed)",
        "bank": {"file": "events_d9_seed78001_n20000_300-180-120_fp1.npz",
                 "sha256": sha(HERE / "results/events/events_d9_seed78001_n20000_300-180-120_fp1.npz"),
                 "pool_seed": 78001, "n_pool": 20000, "alloc": [300, 180, 120]},
        "trajectories": {"shots": 1024, "seed_base": 78002, "seed_rule": "p1core.task_seed(78002, event, mode): shared by every "
                         "condition (common random numbers across layouts AND physics settings)", "decoder": "D0C",
                         "block_events": 4, "tasks_per_setting": 3600, "tasks_total": 3600 * len(grid)},
        "endpoints": {
            "primary": "eta = R(deep) / R(typed) per setting, R = impact-conditioned excess risk (worse input mode, stratified mean "
                       "minus D0C quiet), 95% paired bootstrap percentile interval",
            "secondary": ["shallow/typed ratio per setting", "excess risk per layout per setting", "eta vs energy scale curve (std and echo)"],
            "resolved": "95% interval excludes 1",
            "estimator": "analyze_calext.py: B=4000, seed 20261108, Dirichlet(pool_counts+1), Jeffreys quiet draws (quiet_d9.json D0C), "
                         "identical bootstrap indices for all conditions",
        },
        "order_and_looks": [
            "1. Wave A0 (E1|standard, all three layouts). Look ONLY at the anchor: analyze_calext.py anchor.",
            "2. Anchor PASS iff, for both typed and shallow, (a) the new point estimate lies in the count sweep's nominal 95% interval "
            "(typed [1.74e-4, 3.08e-4], shallow [5.00e-4, 1.06e-3]) AND (b) the 95% paired bootstrap interval of new minus count-sweep "
            "excess risk (same events, same estimator) contains 0. Bit-identity is NOT required: the count sweep ran on macOS arm64 and "
            "Stim seeded sampling is machine-dependent (CONTROL_REPRODUCTION.json).",
            "3. Anchor FAIL -> stop and debug; no other setting is run or read.",
            "4. Anchor PASS -> run waves A, B, C in that order without interim looks; analyse all settings together (analyze_calext.py all).",
        ],
        "decision_rules": {
            "eta_resolved_gt1_at_x100_echo": "Typed placement helps for Willow-scale impacts on a well-decoupled circuit (strongest claim).",
            "eta_gt1_at_x1_to_x10_not_resolved_at_x100": "Claim restricted to sub-threshold (small/medium) impacts; argue from their rate.",
            "eta_gt1_at_x100_standard_not_echo": "Large-impact advantage came from unechoed phase errors; the echo result is the device-relevant one.",
            "trap_x10_changes_eta": "Trapping and placement interact: 'changes sharply' = the trap-x10 eta interval excludes the matching "
                                    "trap-x1 echo point estimate (E1 and E100 separately).",
            "jj_only_gap": "Reported as a sensitivity point; same rule as trap_x10 against E1|echo.",
            "all_outcomes_reported": True,
        },
        "known_limits": [
            "Post hoc extension on the existing 600-impact bank; no independent impact cohort.",
            "The bank's strata were defined on x1 dose; scaling is uniform, so strata and weights remain valid, but at x100 the "
            "stratification is no longer variance-optimal.",
            "Quiet baseline unchanged (decoupling and traps act only on burst channels in this model).",
            "Step 5 x100 saturation (pooled over deep+typed, 20 events): " + json.dumps(timing["x100_saturation_pooled"]),
            "Total ancilla Pauli mass std vs echo differs by ~6% even at x100 (T1/excitation dominate the summed channel).",
            "Recovery-class consistency check of the trap prior failed in calibration (0.86-13 ms); trap x10 departs from it on purpose.",
        ],
        "runtime": {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__,
                    "stim": stim.__version__, "pymatching": pymatching.__version__,
                    "numpy_trapz": "simulation/src/physics/nonmarkov.py already falls back to np.trapezoid; no pin needed",
                    "measured_sec_per_task": {"x1": timing["timing_x1"]["sec_per_task"], "x100": timing["timing_x100"]["sec_per_task"],
                                              "x100_saturation_median": timing["x100_saturation_pooled"]["median_sec_per_task"],
                                              "workers_during_measurement": timing["workers"]}},
        "code_sha256": code,
        "step5_checks": checks,
    }
    OUT.write_text(json.dumps(proto, indent=1))
    print("frozen:", sha(OUT))


if __name__ == "__main__":
    main()
