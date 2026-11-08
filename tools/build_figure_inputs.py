#!/usr/bin/env python3
"""Builds the plotted values for Figs. 2-5 from the recomputed results.

Starts from `assets/figure_data.json` and replaces the panels that this repository recomputes from per-event counts:

    fig2a   deep/typed ratio and 95% interval, d = 9, rows D0, D0C, F, FC        p1_archive / precision_d9
    fig5a   the same at d = 13                                                  p1_archive / precision_d13
    fig5b   shallow - typed excess risk (10^-4) for burst caps of 2, 5, 10 and   p1_archive / regime
            20 rounds and the uncapped reference, rows D0 ("background") and
            F ("oracle")
    fig5c   the same for deep - typed                                           p1_archive / regime

Fig. 2b is drawn directly from `data/count_sweep/RESULTS.json`. Figs. 3 and 4 keep the values in
`assets/figure_data.json`, which were recovered from the manuscript's vector figures (about 0.2% of each axis span).

Usage:
    python tools/build_figure_inputs.py --results_dir out/results --out out/figure_data.json
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

src_path = pathlib.Path(__file__).resolve().parent.parent / "src"
sys.path.insert(0, str(src_path))

from lsburst._io import REPO_ROOT, load_json  # noqa: E402

ROW_DECODERS = ("D0", "D0C", "F", "FC")          # row order of fig2a and fig5a
CAP_CONDITIONS = ("dur2", "dur5", "dur10", "dur20", "rs1")   # x positions of fig5b,c; rs1 is the uncapped reference
CAP_DECODERS = {"background": "D0", "oracle": "F"}


def ratio_rows(comparisons: dict) -> list[list[float]]:
    return [[c["ratio"], *c["ratio_ci"]] for c in (comparisons[f"deep vs typed|{d}"] for d in ROW_DECODERS)]


def cap_rows(contrast: dict, arm: str) -> dict[str, list[list[float]]]:
    out = {}
    for label, decoder in CAP_DECODERS.items():
        rows = []
        for cond in CAP_CONDITIONS:
            c = contrast[f"{cond}|{arm} - {cond}|typed_paper|{decoder}"]
            rows.append([1e4 * c["diff"], 1e4 * c["diff_ci"][0], 1e4 * c["diff_ci"][1]])
        out[label] = rows
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results_dir", type=pathlib.Path, default=pathlib.Path("out/results"))
    parser.add_argument("--out", type=pathlib.Path, default=pathlib.Path("out/figure_data.json"))
    args = parser.parse_args()

    data = load_json(REPO_ROOT / "assets" / "figure_data.json")
    p1 = load_json(args.results_dir / "p1_archive.json")
    data["fig2a"] = ratio_rows(p1["precision_d9"]["comparisons"])
    data["fig5a"] = ratio_rows(p1["precision_d13"]["comparisons"])
    data["fig5b"] = cap_rows(p1["regime"]["contrast"], "shallow")
    data["fig5c"] = cap_rows(p1["regime"]["contrast"], "deep")
    data["_provenance"] = {
        "fig2a, fig5a, fig5b, fig5c": "recomputed from per-event counts (tools/recompute.py, p1_archive)",
        "fig2b": "data/count_sweep/RESULTS.json",
        "fig3, fig4": data["_provenance"]["fig3,fig4,fig5"].replace("replace with frozen summaries when available",
                                                                    "assets/figure_data.json"),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(data, indent=1))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
