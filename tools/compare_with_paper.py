#!/usr/bin/env python3
"""Compares the recomputed results with the numbers printed in the paper.

Reads the JSON files written by `tools/recompute.py`, prints one line per reported number and writes the same table
as Markdown to `--summary`. Disagreements are listed but do not change the exit status; a missing results file does.

Usage:
    python tools/compare_with_paper.py --results_dir out/results --summary out/summary.md
"""

from __future__ import annotations

import argparse
import pathlib
import sys

src_path = pathlib.Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(src_path))

from lsburst._io import load_json  # noqa: E402
from lsburst.paper import all_checks  # noqa: E402

STUDIES = ("count_sweep", "p1_archive", "r9_estimated_decoder", "r1_response_models")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results_dir", type=pathlib.Path, default=pathlib.Path("out/results"))
    parser.add_argument("--summary", type=pathlib.Path, default=pathlib.Path("out/summary.md"))
    args = parser.parse_args()

    missing = [s for s in STUDIES if not (args.results_dir / f"{s}.json").exists()]
    if missing:
        print(f"missing results for {', '.join(missing)}; run step2_recompute_results.sh first", file=sys.stderr)
        return 1
    checks = all_checks({s: load_json(args.results_dir / f"{s}.json") for s in STUDIES})

    width = max(len(c.item) for c in checks)
    for c in checks:
        print(f"{'agrees' if c.agrees else 'DIFFERS':<8} {c.item:<{width}}  paper: {c.paper:<36} data: {c.data}")
    notes = [c for c in checks if c.note]
    for c in notes:
        print(f"\nnote ({c.item}): {c.note}")
    n_agree = sum(c.agrees for c in checks)
    print(f"\n{n_agree} of {len(checks)} reported values agree with the data.")

    lines = ["| Paper item | Paper | This repository | Trace level | Source | Agrees |", "|---|---|---|---|---|---|"]
    lines += [f"| {c.item} | {c.paper} | {c.data} | {c.level} | `{c.source}` | {'yes' if c.agrees else '**no**'} |"
              for c in checks]
    lines += ["", f"{n_agree} of {len(checks)} reported values agree with the data.", ""]
    lines += [f"- {c.item}: {c.note}" for c in notes]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
