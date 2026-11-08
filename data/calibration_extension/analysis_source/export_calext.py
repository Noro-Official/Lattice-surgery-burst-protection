"""Export the calibration extension to <data repo>/data/calibration_extension/ (count_sweep layout)."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
DEST = Path(sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\admin\PACER_CALEXT\lsbp\data\calibration_extension")
WAVES = ("A0", "A", "B", "C")


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    EV = dict(np.load(RES / "events" / "events_d9_seed78001_n20000_300-180-120_fp1.npz"))
    rows = []
    for w in WAVES:
        for f in sorted((RES / f"calext_{w}" / "blocks").glob("*.json")):
            for r in json.loads(f.read_text())["rows"]:
                h = r["hist"]["D0C"]
                rows.append([r["cond"], r["mode"], r["e"], int(EV["stratum"][r["e"]]), r["shots"], r["shots"] - h[0], *h, r["onset"]])
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    keys = {(r[0], r[1], r[2]) for r in rows}
    assert len(keys) == len(rows), "duplicate cells"
    assert all(sum(r[6:10]) == r[4] for r in rows)
    buf = io.StringIO()
    wr = csv.writer(buf, lineterminator="\n")
    wr.writerow(["condition", "input_mode", "event_index", "stratum", "shots", "any_output_failures",
                 "histogram_00", "histogram_01", "histogram_10", "histogram_11", "onset_round"])
    wr.writerows(rows)
    with gzip.GzipFile(DEST / "per_event_counts.csv.gz", "wb", mtime=0) as g:
        g.write(buf.getvalue().encode())
    copies = {
        "FROZEN_PROTOCOL.json": HERE / "FROZEN_PROTOCOL.json", "quiet_d9.json": HERE / "quiet_d9.json",
        "config.json": HERE / "config.json", "RESULTS.json": RES / "analysis_all.json", "ANCHOR.json": RES / "analysis_anchor.json",
        "STEP5_CHECKS.json": RES / "step5_checks.json", "STEP5_TIMING.json": RES / "step5_timing.json",
        "eta_vs_energy.png": RES / "eta_vs_energy.png", "eta_vs_energy.pdf": RES / "eta_vs_energy.pdf",
        "events_d9_seed78001_n20000_300-180-120_fp1.npz": RES / "events" / "events_d9_seed78001_n20000_300-180-120_fp1.npz",
        "recompute_results.py": HERE / "recompute_results.py", "RUN_RESULTS.md": HERE / "RUN_RESULTS.md",
    }
    for k, v in copies.items():
        shutil.copyfile(v, DEST / k)
    for w in WAVES:
        shutil.copyfile(RES / f"calext_{w}" / "study.json", DEST / f"study_calext_{w}.json")
    src = DEST / "analysis_source"
    src.mkdir(exist_ok=True)
    for k in ("p1core.py", "campaign.py", "studies.py", "calext_studies.py", "analyze_calext.py", "step5_checks.py",
              "step5_timing.py", "freeze_protocol.py", "figure_eta.py", "export_calext.py", "run_ABC.sh"):
        shutil.copyfile(HERE / k, src / k)
    shutil.copyfile(HERE / "pacer_continuous.modified.py", src / "pacer_continuous.py")
    sums = {str(p.relative_to(DEST)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(DEST.rglob("*")) if p.is_file() and p.name != "SHA256SUMS.json"}
    (DEST / "SHA256SUMS.json").write_text(json.dumps(sums, indent=1))
    print(f"exported {len(rows)} cells, {len(sums)} files -> {DEST}")


if __name__ == "__main__":
    main()
