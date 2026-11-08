"""The numbers printed in the paper, and their comparison with recomputed results.

Each check formats the recomputed value the way the paper prints it and compares the two strings, so a check passes
only if the paper's rounding of the data gives the printed number.
"""

from __future__ import annotations

import dataclasses

from lsburst._io import DATA_DIR, load_json

PER_EVENT = "recomputed from per-event counts"
STORED = "stored analysis (no per-event counts)"


@dataclasses.dataclass(frozen=True)
class Check:
    """One reported number, compared with the data.

    Attributes:
        item: Where the number appears in the paper.
        paper: The value as printed.
        data: The value from this repository, formatted like the paper.
        level: How far back toward raw data the value is recomputed.
        source: The study or file the value comes from.
        note: Optional context shown under the table.
    """
    item: str
    paper: str
    data: str
    level: str
    source: str
    note: str = ""

    @property
    def agrees(self) -> bool:
        return self.paper == self.data


def _ratio(v: dict) -> str:
    return f"{v['ratio']:.2f} ({v['ratio_ci'][0]:.2f}–{v['ratio_ci'][1]:.2f})"


def count_sweep_checks(r: dict) -> list[Check]:
    lo, hi = r["n_star_ci"]
    slo, shi = r["saving_ci"]
    return [
        Check("Protected ancillas at the shallow-reference risk (Table I, Fig. 2b)", "118.8 (117.0–133.1)",
              f"{r['n_star']:.1f} ({lo:.1f}–{hi:.1f})", PER_EVENT, "count_sweep"),
        Check("Protected-ancilla saving ΔN (abstract)", "49 (35–51)",
              f"{r['saving']:.0f} ({slo:.0f}–{shi:.0f})", PER_EVENT, "count_sweep"),
        Check("Typed excess risk at 117, 123, 132, 150 sites (Table I, 10⁻⁴)", "8.29, 6.19, 5.10, 3.07",
              ", ".join(f"{1e4 * x:.2f}" for x in r["typed_excess"][:4]), PER_EVENT, "count_sweep"),
    ]


def p1_archive_checks(r: dict) -> list[Check]:
    d9 = r["precision_d9"]["comparisons"]
    d13 = r["precision_d13"]["comparisons"]
    eta = d9["deep vs typed|D0C"]
    reduction = [round(100 * (1 - 1 / x)) for x in (eta["ratio"], *eta["ratio_ci"])]
    cap = r["regime"]["contrast"]["dur5|shallow - dur5|typed_paper|D0"]
    return [
        Check("Deep/typed ratio η, correlated background matching, d = 9 (Fig. 2a)", "1.62 (1.49–1.77)",
              _ratio(eta), PER_EVENT, "p1_archive/precision_d9"),
        Check("Typed reduction in excess risk, 1 − 1/η (abstract)", "38% (33–44%)",
              f"{reduction[0]}% ({reduction[1]}–{reduction[2]}%)", PER_EVENT, "p1_archive/precision_d9"),
        Check("η for FC, D0 and F (Fig. 2a)", "1.43, 1.36, 1.28",
              ", ".join(f"{d9[f'deep vs typed|{k}']['ratio']:.2f}" for k in ("FC", "D0", "F")),
              PER_EVENT, "p1_archive/precision_d9"),
        Check("Oracle deep/typed ratios at d = 13, F and FC (Fig. 5a)", "0.96 (0.62–1.56), 1.12 (0.48–2.65)",
              f"{_ratio(d13['deep vs typed|F'])}, {_ratio(d13['deep vs typed|FC'])}", PER_EVENT,
              "p1_archive/precision_d13"),
        Check("Five-round cap, shallow/typed, background matching (Fig. 5b)", "0.60 (0.28–0.97)",
              _ratio(cap), PER_EVENT, "p1_archive/regime"),
    ]


def r9_checks(r: dict) -> list[Check]:
    q = r["quiet_penalty"]["V0_baseline|typed|E|X"]
    return [Check(
        "Estimated-burst decoder break-even burst probability", "about 4.3% per operation",
        f"about {100 * q['break_even_burst_prob_per_op']:.1f}% per operation", PER_EVENT, "r9_estimated_decoder",
        note=(f"The break-even rests on {q['harms']} harmful and {q['rescues']} rescued burst-free shots out of "
              f"{q['shots']:,}; the penalty's 95% interval includes zero."))]


def r1_checks(r: dict) -> list[Check]:
    verdict = {(row["decoder"], row["variant"]): row["verdict"] for row in r["R1"]}
    models = sorted({v for _, v in verdict})
    background = all(verdict[(d, m)].startswith("superior") for d in ("D0", "D0c") for m in models)
    resolved = sorted(f"{d} {m.split('_')[0]}" for (d, m), s in verdict.items()
                      if d in ("F", "Fc") and s != "unresolved")
    return [
        Check("d = 13, three response models: typed lower under both background decoders (SM S5)",
              "yes, all three models", "yes, all three models" if background else "no", PER_EVENT,
              "r1_response_models"),
        Check("d = 13, three response models: oracle comparisons (SM S5)", "not resolved",
              ("resolved for " + ", ".join(resolved)) if resolved else "not resolved", PER_EVENT,
              "r1_response_models",
              note=("Resolved here means the deep-minus-typed difference ΔE has a 95% interval excluding 0. The "
                    "deep/typed ratio intervals for Fc are unbounded because Fc's excess risk is close to zero.")),
    ]


def d9_mechanism_checks() -> list[Check]:
    """Values read from the stored distance-9 mechanism analyses (`data/d9_mechanism_inputs`)."""
    d = DATA_DIR / "d9_mechanism_inputs"
    ra = load_json(d / "factorial" / "analysis.json")["effects"]["shallow - typed_paper|D0"]["RA"]
    rho = load_json(d / "exchange_test" / "analysis.json")["D0"]["spearman_calc"]
    scores = load_json(d / "site_scores_d9_tm9.json")
    pairs = sorted({(round(f, 1), round(g, 1)) for site in scores["W"]["full|D0"]
                    for f, g in zip(scores["W"]["full|D0"][site], scores["W"]["reports|D0"][site])
                    if f is not None and g is not None and round(f, 1) == 7.8 and round(g, 1) == 10.4})
    w_quiet = scores["quiet"]["X"]["D0"]
    return [
        Check("Report × ancilla-Pauli interaction (Fig. 4b, 10⁻⁴)", "1.12 (0.71–1.60)",
              f"{1e4 * ra['effect']:.2f} ({1e4 * ra['ci'][0]:.2f}–{1e4 * ra['ci'][1]:.2f})", STORED,
              "d9_mechanism_inputs/factorial"),
        Check("Path ranking of local swaps, Spearman ρ (SM S4)", "0.77", f"{rho:.2f}", STORED,
              "d9_mechanism_inputs/exchange_test"),
        Check("Path cost W at depth one: all channels / reports only / no burst excess (Fig. 3c)",
              "7.8 / 10.4 / 22.7",
              (f"{pairs[0][0]} / {pairs[0][1]} / {w_quiet:.1f}" if pairs else f"not found / not found / {w_quiet:.1f}"),
              STORED, "d9_mechanism_inputs/site_scores_d9_tm9.json"),
    ]


def all_checks(results: dict) -> list[Check]:
    """All checks, given recomputed results keyed by study name (see `tools/recompute.py`)."""
    return (count_sweep_checks(results["count_sweep"]) + p1_archive_checks(results["p1_archive"])
            + r9_checks(results["r9_estimated_decoder"]) + r1_checks(results["r1_response_models"])
            + d9_mechanism_checks())
