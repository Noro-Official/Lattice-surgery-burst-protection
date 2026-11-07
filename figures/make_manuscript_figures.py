"""PACER manuscript Figs. 2-5 in the publication style of pacer_style.py.

Panel structure and numbering follow PACER_updated_current_20261106.tex.
Fig. 1 (distance-13 layout masks) is kept as supplied.

Inputs
  figure_data.json : Figs. 2a, 3, 4, 5 (provenance recorded inside the file)
  RESULTS.json     : Fig. 2b, refined protection-count sweep (pass path as argv[1])

Toolchain (all open source): SciencePlots base style, cmcrameri colour ramps,
matplotlib-label-lines, highlight_text, URW Nimbus Sans (bundled), colorspacious
palette check, scour (SVG) and pikepdf (PDF metadata, font audit).
"""
import json
import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullFormatter, NullLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pacer_style as S                                   # noqa: E402
from pacer_style import (ax_in, letter, errpt, null_line, direction_cue,  # noqa: E402
                         hl_xlabel, INK, MUTED, FAINT, HAIR, DEEP, SHALLOW,
                         SHALLOW_DARK, ORACLE, ZZ, XX, labelLines)
from highlight_text import ax_text                         # noqa: E402

W = 518.4 / 72                         # match the manuscript's figure width (7.2 in)
D = json.load(open(os.path.join(HERE, "figure_data.json")))

DECODERS = ["Background matching", "Background + correlations",
            "Event-informed oracle", "Oracle + correlations"]
ROW_Y = np.array([3.0, 2.0, 0.8, -0.2])


def decoder_axis(ax, labels=True):
    ax.set_ylim(-0.75, 3.6)
    ax.set_yticks(ROW_Y)
    ax.set_yticklabels(DECODERS if labels else [])
    ax.tick_params(axis="y", length=0, pad=4)
    ax.spines["left"].set_visible(False)


def header(ax, text, colors=(), x=0.0, y=1.04, ha="left"):
    props = [dict(color=c, fontweight="bold") for c in colors]
    ax_text(x, y, text, ax=ax, ha=ha, va="bottom", fontsize=6.5, color=INK,
            highlight_textprops=props or None,
            annotationbbox_kw=dict(xycoords=ax.transAxes))


# ---------------------------------------------------------------------------
def fig2(results_path):
    R = json.load(open(results_path))
    V = R["decoders"]["D0C"]
    counts = np.array(R["protected_counts"], float)
    fig = plt.figure(figsize=(W, 2.25))
    B, H = 0.40, 1.42

    a = ax_in(fig, 1.62, B, 2.05, H)
    for y, (p, lo, hi) in zip(ROW_Y, D["fig2a"]):
        errpt(a, p, y, lo, hi, DEEP, filled=lo > 1)
    a.set_xlim(0.92, 1.85)
    a.set_xticks([1.0, 1.2, 1.4, 1.6, 1.8])
    null_line(a, x=1.0)
    decoder_axis(a)
    hl_xlabel(a, "<Deep> / <typed> excess risk", [DEEP, INK])
    direction_cue(a, 1.0, "deep better", "typed better", y=1.07)
    letter(fig, 0.02, "a")

    b = ax_in(fig, 4.48, B, 2.58, H)
    sc = 1e4
    t = V["target"] * sc
    tlo, thi = (v * sc for v in V["target_ci"])
    xmin, xmax = 108, 212
    b.axhspan(tlo, thi, color="#efefef", lw=0, zorder=0)
    b.axhline(t, color=SHALLOW, lw=0.8, ls=(0, (3.5, 1.8)), zorder=1)
    ax_text(xmax - 1.5, thi * 1.04, "<shallow>, 168 sites", ax=b, ha="right",
            va="bottom", fontsize=6, color=MUTED,
            highlight_textprops=[dict(color=SHALLOW_DARK, fontweight="bold")])
    b.add_patch(Rectangle((xmin, 0.2), 117 - xmin, 40, facecolor="white",
                          edgecolor=HAIR, hatch="//////", lw=0, zorder=0.5))
    b.axvline(117, color=FAINT, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
    b.text((xmin + 117) / 2, 2.0, "infeasible", rotation=90, ha="center", va="center",
           fontsize=5.5, color=MUTED,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    ys = np.array(V["typed_excess"]) * sc
    cis = np.array(V["typed_excess_ci"]) * sc
    b.plot(counts, ys, color=INK, lw=0.8, zorder=3)
    for x, y, (lo, hi) in zip(counts, ys, cis):
        errpt(b, x, y, lo, hi, INK, filled=True, vertical=True, ms=3.6)
    n, (nlo, nhi) = V["n_star"], V["n_star_ci"]
    yb = 0.6
    b.plot([nlo, nhi], [yb, yb], color=INK, lw=1.4, solid_capstyle="butt", zorder=4)
    b.plot(n, yb, "D", ms=3.6, mfc="white", mec=INK, mew=0.9, zorder=5)
    b.plot([n, n], [yb * 1.08, t], color=MUTED, lw=0.45, ls=(0, (1, 1.2)), zorder=2)
    b.text(nhi + 3, yb, f"crossing {n:.1f} [{nlo:.0f}, {nhi:.0f}]   "
           f"saving {V['saving']:.0f} [{V['saving_ci'][0]:.0f}, {V['saving_ci'][1]:.0f}]",
           va="center", ha="left", fontsize=5.8, color=INK)
    b.set_yscale("log")
    b.set_ylim(0.45, 15)
    b.set_yticks([1, 2, 5, 10]); b.set_yticklabels(["1", "2", "5", "10"])
    b.yaxis.set_minor_formatter(NullFormatter())
    b.set_xlim(xmin, xmax)
    b.set_xticks([117, 132, 150, 168, 204]); b.set_xticks([123], minor=True)
    hl_xlabel(b, "Protected ancillas, <typed> band", [INK])
    b.set_ylabel("Excess risk (10$^{-4}$)")
    letter(fig, 3.95, "b")
    return fig


# ---------------------------------------------------------------------------
def fig3():
    fig = plt.figure(figsize=(W, 2.05))
    B, H, Wd = 0.40, 1.36, 1.74
    lefts = [0.66, 2.98, 5.30]
    a, b, c = (ax_in(fig, l, B, Wd, H) for l in lefts)

    def quiet(ax, label=False):
        ax.axhline(0, color=FAINT, lw=0.5, ls=(0, (2, 1.5)), zorder=0)
        if label:
            ax.text(ax.get_xlim()[1], 0.6, "quiet", ha="right", va="bottom",
                    fontsize=5.5, color=MUTED)

    for (k, ys), col, mk in zip(D["fig3a"].items(), ZZ, ["^", "o", "s"]):
        a.plot(np.arange(1, len(ys) + 1), ys, color=col, marker=mk, ms=2.6, mew=0,
               lw=0.9, zorder=3)
    a.text(0.75, D["fig3a"]["5"][0], "5", color=ZZ[0], fontsize=6, ha="right", va="center")
    a.text(0.75, D["fig3a"]["9"][0] + 0.7, "9", color=ZZ[1], fontsize=6, ha="right", va="center")
    a.text(0.75, D["fig3a"]["13"][0] - 0.6, "13", color=ZZ[2], fontsize=6, ha="right", va="center")
    a.text(0.15, -3.2, "merge\nrounds", fontsize=5.5, color=MUTED, ha="center",
           va="center", linespacing=0.95)
    a.set_xlim(-0.6, 8.4); a.set_xticks([1, 2, 4, 6, 8])
    a.set_ylim(-18.8, 2); a.set_yticks([-15, -10, -5, 0])
    quiet(a, label=True)
    header(a, "$ZZ$ merge, $d$ = 9")
    a.set_ylabel("Path-cost change, $W - W_{\\rm quiet}$")
    a.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 0.02, "a")

    lines = []
    for (dd, ys), col in zip(D["fig3b"].items(), XX):
        ln, = b.plot(np.arange(1, len(ys) + 1), ys, color=col, marker="o", ms=2.2,
                     mew=0, lw=0.9, zorder=3, label=f"$d$ = {dd}")
        lines.append(ln)
    labelLines(lines, xvals=[4.5, 5.5, 7.5, 8.5], fontsize=5.5, align=False,
               outline_width=0, zorder=6,
               bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    b.annotate("plateau set by\nreadout of A", xy=(5.2, -28.6), xytext=(6.8, -30.4),
               fontsize=5.5, color=MUTED, va="center", ha="left", linespacing=0.95,
               arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED, shrinkA=1,
                               shrinkB=1.5))
    b.set_xlim(0.3, 12.6); b.set_xticks([1, 4, 8, 12])
    b.set_ylim(-33, 2.2); b.set_yticks([-30, -20, -10, 0])
    quiet(b)
    header(b, "$XX$ merge, $d$ = 7–13")
    b.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 2.52, "b")

    sty = {"All channels": (ZZ[2], "-", "o"),
           "Reports only": (ZZ[1], (0, (3.2, 1.4)), "s"),
           "Pauli + reset": (ZZ[0], (0, (1, 1.1)), "^")}
    for k, ys in D["fig3c"].items():
        colr, ls, mk = sty[k]
        c.plot(np.arange(1, len(ys) + 1), ys, color=colr, ls=ls, marker=mk, ms=2.6,
               mew=0, lw=0.95, label=k, zorder=3)
    c.legend(loc="lower right", handlelength=2.6, borderaxespad=0.1,
             labelspacing=0.3, handletextpad=0.5, fontsize=6)
    c.set_xlim(0.4, 8.4); c.set_xticks([1, 2, 4, 6, 8])
    c.set_ylim(-18.8, 2); c.set_yticks([-15, -10, -5, 0])
    quiet(c)
    header(c, "Helper $Z$ checks, reference dose")
    c.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 4.84, "c")
    return fig


# ---------------------------------------------------------------------------
FIG4_ROWS = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
             (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
FIG4B_LABELS = ["Report", "Ancilla Pauli", "Reset",
                "Report × ancilla", "Report × reset", "Ancilla × reset"]


def fig4():
    fig = plt.figure(figsize=(W, 2.1))
    B, H = 0.40, 1.30
    ys = np.arange(8)[::-1].astype(float)
    m = ax_in(fig, 0.42, B, 0.50, H)
    a = ax_in(fig, 1.00, B, 2.95, H, sharey=m)
    b = ax_in(fig, 5.15, B, 1.92, H)
    for ax in (m, a):
        ax.axhline(3.5, color=HAIR, lw=0.6, zorder=0)
    for y, (r, an, p) in zip(ys, FIG4_ROWS):
        on = [i for i, v in enumerate((r, an, p)) if v]
        if len(on) > 1:
            m.plot([on[0], on[-1]], [y, y], color=INK, lw=0.8, zorder=2)
        for i, v in enumerate((r, an, p)):
            m.plot(i, y, "o", ms=3.6, mfc=INK if v else "white",
                   mec=INK if v else FAINT, mew=0.6, zorder=3)
    m.set_xlim(-0.55, 2.55); m.set_ylim(-0.6, 7.6)
    m.set_xticks([0, 1, 2])
    m.set_xticklabels(["Report", "Ancilla", "Reset"], rotation=45, ha="left",
                      rotation_mode="anchor", fontsize=6)
    m.xaxis.tick_top(); m.tick_params(axis="x", length=0, pad=1.5)
    m.tick_params(axis="y", length=0, labelleft=False)
    for sp in m.spines.values():
        sp.set_visible(False)
    m.text(-0.75, 1.5, "reports on", rotation=90, ha="center", va="center",
           fontsize=5.5, color=MUTED)
    m.text(-0.75, 5.5, "reports off", rotation=90, ha="center", va="center",
           fontsize=5.5, color=MUTED)
    null_line(a, x=0)
    for y, rowkey, (pt, lo, hi) in zip(ys, FIG4_ROWS, D["fig4a"]):
        if rowkey == (0, 0, 0):
            a.plot(0, y, "o", ms=2.6, mfc="white", mec=MUTED, mew=0.6, zorder=5)
            continue
        errpt(a, pt, y, lo, hi, SHALLOW, filled=lo > 0 or hi < 0)
    a.set_xlim(-1.5, 7.6); a.set_xticks([0, 2, 4, 6])
    a.spines["left"].set_visible(False)
    a.tick_params(axis="y", length=0, labelleft=False)
    hl_xlabel(a, "<Shallow> − <typed> failure probability (10$^{-4}$)", [SHALLOW, INK])
    direction_cue(a, 0, "", "typed better", y=1.06)
    letter(fig, 0.02, "a")

    ey = np.array([5.0, 4.0, 3.0, 1.5, 0.5, -0.5])
    null_line(b, x=0)
    for i, (y, (pt, lo, hi)) in enumerate(zip(ey, D["fig4b"])):
        errpt(b, pt, y, lo, hi, SHALLOW if i < 3 else SHALLOW_DARK,
              filled=lo > 0 or hi < 0, ms=3.8)
    b.set_yticks(ey); b.set_yticklabels(FIG4B_LABELS)
    b.tick_params(axis="y", length=0)
    b.spines["left"].set_visible(False)
    b.set_ylim(-1.0, 5.6); b.set_xlim(-1.0, 5.3); b.set_xticks([0, 2, 4])
    b.axhline(2.25, color=HAIR, lw=0.5)
    b.text(5.3, 5.45, "main effects", ha="right", va="center", fontsize=5.5,
           color=MUTED, style="italic")
    b.text(5.3, 2.0, "interactions", ha="right", va="center", fontsize=5.5,
           color=MUTED, style="italic")
    b.set_xlabel("Factorial effect (10$^{-4}$)")
    letter(fig, 4.15, "b")
    return fig


# ---------------------------------------------------------------------------
CAPS = [2, 5, 10, 20]
NOCAP_X = 42.0


def cap_panel(ax, data, title, title_colors, ylim, yticks, show_x, legend):
    null_line(ax, y=0)
    k = 1.075
    for key, col, mk in (("background", INK, "o"), ("oracle", ORACLE, "s")):
        f = 1 / k if key == "background" else k
        pts = data[key]
        xs = np.array(CAPS + [NOCAP_X], float) * f
        ax.plot(xs[:4], [q[0] for q in pts[:4]], color=col, lw=0.5, alpha=0.45, zorder=3)
        for x, (p, lo, hi) in zip(xs, pts):
            errpt(ax, x, p, lo, hi, col, filled=lo > 0 or hi < 0, vertical=True,
                  marker=mk, ms=3.4 if mk == "o" else 3.0)
    ax.axvline(30, color=HAIR, lw=0.6, zorder=0)
    ax.set_xscale("log"); ax.set_xlim(1.6, 56)
    ax.set_xticks(CAPS + [NOCAP_X])
    ax.set_xticklabels([str(c) for c in CAPS] + ["No cap"] if show_x else [])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(*ylim); ax.set_yticks(yticks)
    header(ax, title, title_colors, x=0.02, y=0.88)
    if legend:
        ax.plot([0.035], [0.70], "o", ms=3.4, color=INK, transform=ax.transAxes)
        ax.text(0.065, 0.70, "Background matching", transform=ax.transAxes, fontsize=6,
                va="center")
        ax.plot([0.035], [0.57], "s", ms=3.0, mfc="white", mec=ORACLE, mew=0.8,
                transform=ax.transAxes)
        ax.text(0.065, 0.57, "Event-informed oracle", transform=ax.transAxes, fontsize=6,
                va="center", color=ORACLE)


def fig5():
    fig = plt.figure(figsize=(W, 3.15))
    a = ax_in(fig, 1.62, 0.40, 2.05, 2.35)
    for y, (p, lo, hi) in zip(ROW_Y, D["fig5a"]):
        errpt(a, p, y, lo, hi, DEEP, filled=lo > 1 or hi < 1)
    a.set_xscale("log"); a.set_xlim(0.4, 9.0)
    a.set_xticks([0.5, 1, 2, 4, 8]); a.set_xticklabels(["0.5", "1", "2", "4", "8"])
    a.xaxis.set_minor_formatter(NullFormatter())
    null_line(a, x=1.0)
    decoder_axis(a)
    hl_xlabel(a, "<Deep> / <typed> excess risk, $d$ = 13", [DEEP, INK])
    direction_cue(a, 1.0, "deep better", "typed better", y=1.035)
    letter(fig, 0.02, "a")

    b = ax_in(fig, 4.62, 1.80, 2.45, 1.00)
    c = ax_in(fig, 4.62, 0.40, 2.45, 1.00)
    cap_panel(b, D["fig5b"], "<Shallow> − <typed>", [SHALLOW, INK], (-2.4, 9.5),
              [0, 4, 8], show_x=False, legend=True)
    cap_panel(c, D["fig5c"], "<Deep> − <typed>", [DEEP, INK], (-1.6, 4.2),
              [-1, 0, 1, 2, 3, 4], show_x=True, legend=False)
    c.set_xlabel("Burst-duration cap (rounds of 1.1 µs)")
    fig.text(4.20 / W, 1.60 / 3.15, "Excess-risk difference (10$^{-4}$)", rotation=90,
             ha="center", va="center", fontsize=7)
    letter(fig, 4.05, "b")
    letter(fig, 4.05, "c", y_in=1.42)
    return fig


if __name__ == "__main__":
    out = os.path.join(HERE, "out")
    os.makedirs(out, exist_ok=True)
    jobs = (("fig02_placement_resources", lambda: fig2(sys.argv[1])),
            ("fig03_boundary_paths", fig3),
            ("fig04_channel_interventions", fig4),
            ("fig05_decoder_duration", fig5))
    for name, fn in jobs:
        f = fn()
        pdf = os.path.join(out, f"{name}.pdf")
        f.savefig(pdf)
        f.savefig(os.path.join(out, f"{name}.png"), dpi=600)
        svg = os.path.join(out, f"{name}.svg")
        f.savefig(svg)
        S._postprocess(pdf, svg, name)
        plt.close(f)
    print("done")
