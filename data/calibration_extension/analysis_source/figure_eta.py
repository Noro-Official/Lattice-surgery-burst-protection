"""eta = R(deep)/R(typed) vs impact energy scale; reads results/analysis_all.json, writes results/eta_vs_energy.{png,pdf}."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
rows = json.loads((HERE / "results" / "analysis_all.json").read_text())

INK, INK2, GRID = "#1f1f1e", "#5c5b55", "#e4e3de"
SERIES = [  # (label, selector, color, marker, x-dodge factor) - fixed categorical order
    ("Standard decoupling", lambda r: r["decoupling"] == "standard" and r["trap_multiplier"] == 1 and r["protected_gap_ghz"] == 12, "#2a78d6", "o", 0.94),
    ("Echo", lambda r: r["decoupling"] == "echo" and r["trap_multiplier"] == 1 and r["protected_gap_ghz"] == 12, "#eb6834", "s", 1.0),
    ("Echo + traps ×10 on protected sites", lambda r: r["trap_multiplier"] == 10, "#1baf7a", "^", 1.06),
    ("Echo + 10 GHz protected gap (JJ-Only)", lambda r: r["protected_gap_ghz"] == 10, "#eda100", "D", 1.12),
]

plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=200)
ax.axhline(1.0, color=INK2, lw=1, ls=(0, (4, 3)), zorder=1)
ax.text(140, 1.005, "no placement benefit", color=INK2, fontsize=8, va="bottom", ha="right")
for label, sel, col, mk, dx in SERIES:
    pts = sorted((r for r in rows if sel(r)), key=lambda r: r["energy_scale"])
    if not pts:
        continue
    x = [r["energy_scale"] * dx for r in pts]
    y = [r["eta_deep_over_typed"] for r in pts]
    lo = [r["eta_deep_over_typed"] - r["eta_ci"][0] for r in pts]
    hi = [r["eta_ci"][1] - r["eta_deep_over_typed"] for r in pts]
    if len(pts) > 2:                      # only the full energy series are joined; sparse variants stay as points
        ax.plot(x, y, color=col, lw=2, zorder=2)
    ax.errorbar(x, y, yerr=[lo, hi], fmt=mk, ms=7, color=col, ecolor=col, elinewidth=1.5, capsize=0,
                mec="white", mew=1.5, zorder=3, label=label)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xticks([1, 3, 10, 30, 100])
ax.set_xticklabels(["×1\npaper", "×3", "×10", "×30", "×100\nWillow-calibrated"])
ax.set_xlim(0.7, 150)
ax.set_xlabel("Impact energy scale (relative to the paper's event bank)")
ax.set_ylabel("η = excess risk, deep ÷ typed")
ax.grid(True, which="major", axis="y", color=GRID, lw=0.8)
ax.set_yticks([0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_ylim(0.9, 2.0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="best")
ax.set_title("Typed placement vs impact energy, d = 9 lattice-surgery CNOT (D0C)", fontsize=10, color=INK, loc="left")
fig.text(0.01, 0.005, "Whiskers: 95% paired stratified bootstrap, 600 impacts, 1,024 shots per event, layout and input mode.",
         fontsize=7, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(HERE / "results" / f"eta_vs_energy.{ext}")
print("wrote figure")
