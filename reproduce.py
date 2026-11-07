"""Check every number in the paper that this repository can reproduce.

    python reproduce.py

Steps
  1. Verify SHA-256 checksums of every data file (SHA256SUMS.txt).
  2. Recompute the protection-count sweep, the R9 decoder study and the R1 response-model study from per-event
     failure counts, and require the recomputed analyses to match the stored ones exactly.
  3. Read the stored analysis outputs for the distance-9 mechanism studies.
  4. Compare each result with the value printed in the paper.

The exit code is 1 if a checksum fails or a recomputation does not match the stored analysis. Paper statements that
disagree with the data are listed separately; they are reported, not counted as failures.
A Markdown summary is written to results/summary.md.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
DATA = ROOT / "data"

rows: list[tuple[str, str, str, str, str]] = []   # (paper item, paper value, this repository, trace level, status)
problems: list[str] = []                          # reproduction failures (exit code 1)
notes: list[str] = []                             # paper statements that disagree with the data


def record(item, paper, got, level, ok):
    rows.append((item, paper, got, level, "match" if ok else "DIFFERS"))
    if not ok:
        notes.append(f"{item}: paper says {paper}, data give {got}.")


def load(path):
    return json.loads(Path(path).read_text())


def run(script, *args):
    out = subprocess.run([sys.executable, "-I", str(script), *map(str, args)], cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        problems.append(f"{script.relative_to(ROOT)} exited with {out.returncode}:\n{out.stderr[-2000:]}")
    return out.stdout


def same(a, b, path=""):
    """Differences between two JSON trees (floats to 1e-9 relative)."""
    diffs = []
    if isinstance(a, dict):
        for k, v in a.items():
            if k in ("provenance", "authenticity", "source"):
                continue
            if not isinstance(b, dict) or k not in b:
                diffs.append(f"{path}/{k} missing")
            else:
                diffs += same(v, b[k], f"{path}/{k}")
    elif isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b):
            return [f"{path} length"]
        for i, (x, y) in enumerate(zip(a, b)):
            diffs += same(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and isinstance(b, (int, float)):
        if not (math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-15) or (math.isnan(a) and math.isnan(b))):
            diffs.append(f"{path}: {a} vs {b}")
    elif a != b:
        diffs.append(f"{path}: {a!r} vs {b!r}")
    return diffs


def checksums():
    bad, n = [], 0
    for line in (ROOT / "SHA256SUMS.txt").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        n += 1
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            bad.append(name)
    if bad:
        problems.append("checksum mismatch: " + ", ".join(bad))
    print(f"[1/4] checksums: {n - len(bad)}/{n} files verified")


def count_sweep():
    out = run(DATA / "count_sweep" / "recompute_results.py")   # asserts equality with RESULTS.json internally
    if not out:
        return
    r = json.loads(out)
    lo, hi = r["count_interval"]
    slo, shi = r["saving_interval"]
    record("Protected ancillas at shallow-reference risk (Table I, Fig. 2b)", "118.8 (117.0–133.1)",
           f"{r['required_count']:.1f} ({lo:.1f}–{hi:.1f})", "recomputed from per-event counts",
           f"{r['required_count']:.1f} ({lo:.1f}–{hi:.1f})" == "118.8 (117.0–133.1)")
    record("Saving ΔN versus shallow (abstract)", "49 (35–51)", f"{r['saving']:.0f} ({slo:.0f}–{shi:.0f})",
           "recomputed from per-event counts", f"{r['saving']:.0f} ({slo:.0f}–{shi:.0f})" == "49 (35–51)")
    stored = load(DATA / "count_sweep" / "RESULTS.json")["decoders"]["D0C"]["typed_excess"]
    paper = [8.29, 6.19, 5.10, 3.07]          # Table I, typed 117/123/132/150, excess risk in 1e-4
    got = [round(x * 1e4, 2) for x in stored[:4]]
    record("Table I typed excess risk, 117–150 sites (10⁻⁴)", ", ".join(f"{x:.2f}" for x in paper),
           ", ".join(f"{x:.2f}" for x in got), "recomputed from per-event counts", got == paper)
    print("[2/4] count sweep recomputed and matches RESULTS.json")


def r9():
    new = RESULTS / "r9_RESULTS.json"
    run(ROOT / "analysis" / "r9_estimated_burst_decoder.py", new)
    if not new.exists():
        return
    diffs = same(load(DATA / "r9_estimated_burst_decoder" / "RESULTS.json"), load(new))
    if diffs:
        problems.append("R9 recomputation differs from stored RESULTS.json: " + "; ".join(diffs[:5]))
    q = load(new)["quiet_penalty"]["V0_baseline|typed|E|X"]
    be = q["break_even_burst_prob_per_op"]
    record("Estimated-burst decoder break-even burst probability", "about 4.3% per operation", f"{100 * be:.1f}%",
           "recomputed from per-event counts", round(100 * be, 1) == 4.3)
    rows.append(("  … false-alarm evidence behind it", "", f"{q['harms']} harm, {q['rescues']} rescues in "
                 f"{q['shots']:,} burst-free shots; penalty interval includes 0", "", "see README"))
    print(f"[2/4] R9 recomputed: {'matches' if not diffs else 'DIFFERS FROM'} stored RESULTS.json")


def r1():
    new = RESULTS / "r1_RESULTS.json"
    run(ROOT / "analysis" / "r1_d13_response_models.py", new)
    if not new.exists():
        return
    got = load(new)["R1"]
    diffs = same(load(DATA / "r1_d13_response_models" / "RESULTS.json")["R1"], got)
    if diffs:
        problems.append("R1 recomputation differs from stored RESULTS.json: " + "; ".join(diffs[:5]))
    v = {(r["decoder"], r["variant"]): r["verdict"] for r in got}
    bg = all(v[(d, m)].startswith("superior") for d in ("D0", "D0c") for m in {k[1] for k in v})
    record("d = 13, three response models: typed lower under both background decoders (SM S5)",
           "yes, all three models", "yes, all three models" if bg else "no",
           "recomputed from per-event counts", bg)
    resolved = sorted(f"{d} {m.split('_')[0]}" for (d, m), s in v.items() if d in ("F", "Fc") and s != "unresolved")
    record("d = 13, three response models: oracle comparisons (SM S5)", "not resolved",
           "resolved (ΔE interval excludes 0) for " + ", ".join(resolved) if resolved else "not resolved",
           "recomputed from per-event counts", not resolved)
    print(f"[2/4] R1 recomputed: {'matches' if not diffs else 'DIFFERS FROM'} stored results")


def stored_d9():
    d = DATA / "d9_mechanism_inputs"
    ra = load(d / "factorial" / "analysis.json")["effects"]["shallow - typed_paper|D0"]["RA"]
    got = f"{ra['effect'] * 1e4:.2f} ({ra['ci'][0] * 1e4:.2f}–{ra['ci'][1] * 1e4:.2f})"
    record("Report × ancilla-Pauli interaction (Fig. 4b), 10⁻⁴", "1.12 (0.71–1.60)", got,
           "stored analysis (no per-event counts)", got == "1.12 (0.71–1.60)")
    rho = load(d / "exchange_test" / "analysis.json")["D0"]["spearman_calc"]
    record("Path ranking of local swaps, Spearman ρ (SM S4)", "0.77", f"{rho:.2f}",
           "stored analysis (no per-event counts)", f"{rho:.2f}" == "0.77")
    s = load(d / "site_scores_d9_tm9.json")
    w_quiet = s["quiet"]["X"]["D0"]
    pairs = sorted({(round(f, 1), round(r, 1)) for site in s["W"]["full|D0"]
                    for f, r in zip(s["W"]["full|D0"][site], s["W"]["reports|D0"][site])
                    if f is not None and r is not None and round(f, 1) == 7.8 and round(r, 1) == 10.4})
    record("Path cost W at depth one: all channels / reports only / no burst excess (Fig. 3c)", "7.8 / 10.4 / 22.7",
           f"{pairs[0][0]} / {pairs[0][1]} / {w_quiet:.1f}" if pairs else f"no 7.8/10.4 pair / {w_quiet:.1f}",
           "stored isolated-site scores", bool(pairs) and f"{w_quiet:.1f}" == "22.7")
    print("[3/4] distance-9 mechanism outputs read")


def archive_p1():
    """Stored analyses from the 31 Oct 2026 publication package, verified against that package's own SHA256SUMS."""
    d = DATA / "p1_archive"
    sums = {}
    for line in (d / "PACKAGE_SHA256SUMS.txt").read_text().splitlines():
        h, f = line.split(maxsplit=1)
        sums[f] = h
    bad = [str(p.relative_to(d)) for p in d.rglob("*") if p.is_file() and p.name != "PACKAGE_SHA256SUMS.txt"
           and p.suffix != ".md" and p.name != "blocks.jsonl.gz"     # blocks are checked line by line below
           and sums.get("experiments/pacer_p1/results/" + str(p.relative_to(d))) != hashlib.sha256(p.read_bytes()).hexdigest()]
    if bad:
        problems.append("files differ from the original package checksums: " + ", ".join(bad))

    # Recompute from the 15,600 per-event blocks (each checked against the package checksum) and require an exact
    # match with the archived analyses.
    new = RESULTS / "p1_archive_RESULTS.json"
    run(ROOT / "analysis" / "p1_archive_cohorts.py", new)
    if not new.exists():
        return
    R = load(new)
    stored = load(d / "precision_analysis.json")
    keys = ("diff", "diff_ci", "ratio", "ratio_ci", "class")
    diffs = [f"{s}/{c}: {x}" for s in ("precision_d9", "precision_d13") for c, v in R[s]["comparisons"].items()
             for x in same({k: stored[s]["comparisons"][c][k] for k in keys}, v)]
    stored_rg = load(d / "regime" / "analysis.json")["contrast"]
    diffs += [f"regime/{c}: {x}" for c, v in R["regime"]["contrast"].items() for x in same(stored_rg[c], v)]
    if diffs:
        problems.append("P1 recomputation differs from the archived analyses: " + "; ".join(diffs[:5]))
    n_cmp = sum(len(R[s]["comparisons"]) for s in ("precision_d9", "precision_d13")) + len(R["regime"]["contrast"])
    print(f"[2/4] P1 cohorts recomputed from {sum(R[s]['blocks'] for s in R):,} blocks: "
          f"{'all' if not diffs else 'NOT all'} {n_cmp} archived comparisons reproduced")

    level = "recomputed from per-event counts"
    fmt = lambda v: f"{v['ratio']:.2f} ({v['ratio_ci'][0]:.2f}–{v['ratio_ci'][1]:.2f})"
    c9 = R["precision_d9"]["comparisons"]
    got = fmt(c9["deep vs typed|D0C"])
    record("Deep/typed ratio η at d = 9, correlated background matching (Fig. 2a)", "1.62 (1.49–1.77)", got, level,
           got == "1.62 (1.49–1.77)")
    eta = c9["deep vs typed|D0C"]
    red = [round(100 * (1 - 1 / x)) for x in (eta["ratio"], eta["ratio_ci"][0], eta["ratio_ci"][1])]
    got = f"{red[0]}% ({red[1]}–{red[2]}%)"
    record("Typed reduction in excess risk, 1 − 1/η (abstract)", "38% (33–44%)", got, level, got == "38% (33–44%)")
    got = ", ".join(f"{c9[f'deep vs typed|{k}']['ratio']:.2f}" for k in ("FC", "D0", "F"))
    record("η for the other three decoders (Fig. 2a)", "1.43, 1.36, 1.28", got, level, got == "1.43, 1.36, 1.28")
    c13 = R["precision_d13"]["comparisons"]
    got = f"{fmt(c13['deep vs typed|F'])}, {fmt(c13['deep vs typed|FC'])}"
    record("d = 13 oracle ratios, 1,200 impacts (Fig. 5a)", "0.96 (0.62–1.56), 1.12 (0.48–2.65)", got, level,
           got == "0.96 (0.62–1.56), 1.12 (0.48–2.65)")
    rg = R["regime"]["contrast"]
    got = fmt(rg["dur5|shallow - dur5|typed_paper|D0"])
    record("Five-round burst cap, shallow/typed ratio, background matching (Fig. 5b)", "0.60 (0.28–0.97)", got, level,
           got == "0.60 (0.28–0.97)")


NOT_HERE = []


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    checksums()
    count_sweep()
    r9()
    r1()
    stored_d9()
    archive_p1()
    for item, paper in NOT_HERE:
        rows.append((item, paper, "not in this repository", "figures/figure_data.json only", "not reproducible"))

    width = max(len(r[0]) for r in rows)
    print("\n[4/4] comparison with the paper\n")
    for item, paper, got, level, status in rows:
        print(f"  {status:<16} {item:<{width}}  paper: {paper:<36} here: {got}")
    md = ["| Paper item | Paper | This repository | Trace level | Status |", "|---|---|---|---|---|"]
    md += [f"| {' | '.join(r)} |" for r in rows]
    if notes:
        print("\nPaper statements that disagree with the data:")
        md += ["", "**Paper statements that disagree with the data**", ""]
        for n in notes:
            print("  - " + n)
            md.append("- " + n)
    if problems:
        print("\nREPRODUCTION FAILED:")
        md += ["", "**Reproduction failed**", ""]
        for p in problems:
            print("  - " + p)
            md.append("- " + p)
    else:
        print("\nAll checksums verified and all recomputations match the stored analyses.")
    (RESULTS / "summary.md").write_text("\n".join(md) + "\n")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
