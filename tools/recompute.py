#!/usr/bin/env python3
"""Recomputes each study from its per-event counts and checks the result against the stored analysis.

Writes one JSON file per study to `--out_dir`. Exits with status 1 if any recomputed value differs from the stored
analysis (floats compared to a relative tolerance of 1e-9).

Studies:
    count_sweep           refined protection-count sweep, d = 9          data/count_sweep/
    p1_archive            original Fig. 2a, 5a and 5b,c cohorts           data/p1_archive/
    r9_estimated_decoder  estimated-burst and realistic decoders          data/r9_estimated_burst_decoder/
    r1_response_models    three response models at d = 13                data/r1_d13_response_models/

Usage:
    python tools/recompute.py --out_dir out/results
    python tools/recompute.py --out_dir out/results --studies count_sweep r9_estimated_decoder
"""

from __future__ import annotations

import argparse
import importlib
import pathlib
import sys
import time

src_path = pathlib.Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(src_path))

from lsburst._io import tree_differences, write_json  # noqa: E402

STUDIES = ("count_sweep", "p1_archive", "r9_estimated_decoder", "r1_response_models")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out_dir", type=pathlib.Path, default=pathlib.Path("out/results"))
    parser.add_argument("--studies", nargs="+", choices=STUDIES, default=list(STUDIES))
    args = parser.parse_args()

    failed = False
    for name in args.studies:
        module = importlib.import_module(f"lsburst.{name}")
        t0 = time.perf_counter()
        result = module.analyze()
        diffs = tree_differences(module.stored_results(), result)
        write_json(args.out_dir / f"{name}.json", result)
        status = "matches the stored analysis" if not diffs else f"DIFFERS from the stored analysis ({len(diffs)})"
        print(f"{name:<22} {status}  [{time.perf_counter() - t0:.1f} s]")
        for d in diffs[:5]:
            print(f"    {d}", file=sys.stderr)
        failed |= bool(diffs)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
