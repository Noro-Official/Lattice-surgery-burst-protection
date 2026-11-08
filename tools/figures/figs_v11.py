"""PACER manuscript, redesigned Figs. 2 and 5 (v11).

What changed relative to make_manuscript_figures.py
  Fig. 2a  log ratio axis shared with Fig. 5a; decoder rows grouped by the
           information the decoder receives; exact values printed beside each
           row; a secondary top axis converts the ratio eta into the risk
           reduction 1 - 1/eta, so the 38% headline is read directly.
  Fig. 2b  the crossing is drawn where it happens, on the shallow reference
           line, and the protected-ancilla saving is drawn as a dimension arrow
           from the crossing to the shallow band at 168 sites.  The reference
           interval is two hairlines instead of a filled box.
  Fig. 5a  same row grammar as Fig. 2a; the distance-9 ratios from Fig. 2a are
           overlaid as grey ghost markers so the change with distance is visible.
  Fig. 5b,c  per-panel axis labels (no shared label colliding with the panel
           letter), a top axis in microseconds, an axis break before the
           no-cap reference, a marker key that states the filled/open rule, and
           the five-round reversal is annotated.

Open-source components
  SciencePlots          github.com/garrettj403/SciencePlots       base 'nature' style
  cmcrameri             github.com/callumrollo/cmcrameri          perceptual colour ramps
  highlight_text        github.com/znstrider/highlight_text       colour-coded label words
  colorspacious         github.com/njsmith/colorspacious          CVD simulation of the render
  pikepdf / scour       github.com/pikepdf/pikepdf, scour-project font audit, SVG cleanup
  URW Nimbus Sans       github.com/ArtifexSoftware/urw-base35-fonts  Helvetica metrics
"""
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator, FixedLocator, FixedFormatter
import matplotlib.transforms as mtransforms

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pacer_style as S                                               # noqa: E402
from pacer_style import (ax_in, letter, errpt, null_line, direction_cue,  # noqa: E402
                         hl_xlabel, INK, MUTED, FAINT, HAIR, DEEP, SHALLOW,
                         SHALLOW_DARK, ORACLE)
from highlight_text import ax_text                                    # noqa: E402


W = 518.4 / 72
# Plotted values: out/figure_data.json from tools/build_figure_inputs.py (override with FIGURE_DATA).
REPO = os.path.dirname(os.path.dirname(HERE))
FIGURE_DATA = os.environ.get("FIGURE_DATA", os.path.join(REPO, "out", "figure_data.json"))
D = json.load(open(FIGURE_DATA))
GHOST = "#a3a3a3"

# Decoder rows: (group header, row label, primary?)  -- order matches figure_data
ROWS = [("Background noise model", "Matching", False),
        ("Background noise model", "+ correlations", True),
        ("Event-informed oracle", "Matching", False),
        ("Event-informed oracle", "+ correlations", False)]
ROW_Y = np.array([3.0, 2.15, 0.85, 0.0])
GROUP_Y = {"Background noise model": 3.72, "Event-informed oracle": 1.50}
RULE_Y = 1.84


def decoder_rows(ax, gutter_in=1.24):
    ax.set_ylim(-0.55, 4.15)
    ax.set_yticks(ROW_Y)
    ax.set_yticklabels([r[1] for r in ROWS])
    for tl, (_, _, prim) in zip(ax.get_yticklabels(), ROWS):
        if prim:
            tl.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=5)
    ax.spines["left"].set_visible(False)
    w_in = ax.get_position().width * ax.figure.get_size_inches()[0]
    x0 = -gutter_in / w_in
    tr = ax.get_yaxis_transform()
    for g, y in GROUP_Y.items():
        ax.text(x0, y, g, transform=tr, ha="left", va="center", fontsize=6,
                color=MUTED, style="italic")
    ax.plot([x0, 1.0], [RULE_Y, RULE_Y], transform=tr, color=HAIR, lw=0.5,
            clip_on=False, zorder=0)


def log_ratio_axis(ax, lo, hi, ticks):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())


# ---------------------------------------------------------------------------
def fig2(results_path):
    R = json.load(open(results_path))
    V = R["decoders"]["D0C"]
    counts = np.array(R["protected_counts"], float)
    fig = plt.figure(figsize=(W, 2.55))
    B, H = 0.42, 1.55

    # ---- a: placement ratio at d = 9 ------------------------------------
    a = ax_in(fig, 1.30, B, 1.78, H)
    for y, (p, lo, hi), (_, _, prim) in zip(ROW_Y, D["fig2a"], ROWS):
        errpt(a, p, y, lo, hi, DEEP, filled=lo > 1, ms=4.4 if prim else 3.6,
              lw=1.15 if prim else 0.8)
        a.text(1.04, y, f"{p:.2f} ({lo:.2f}–{hi:.2f})", transform=a.get_yaxis_transform(),
               ha="left", va="center", fontsize=6, color=INK,
               fontweight="bold" if prim else "normal")
    a.text(1.04, 4.02, "$\\eta$ (95% CI)", transform=a.get_yaxis_transform(), ha="left",
           va="center", fontsize=6, color=MUTED)
    log_ratio_axis(a, 0.94, 2.0, [1.0, 1.2, 1.4, 1.6, 1.8])
    null_line(a, x=1.0)
    decoder_rows(a)
    hl_xlabel(a, "<Deep> / <typed> excess risk, $\\eta$", [DEEP, INK])
    # top axis: risk reduction for typed, 1 - 1/eta (ticks placed at the exact eta)
    top = a.twiny()
    top.set_xscale("log"); top.set_xlim(a.get_xlim())
    red = [0, 20, 30, 40]
    top.xaxis.set_major_locator(FixedLocator([1 / (1 - r / 100) for r in red]))
    top.xaxis.set_major_formatter(FixedFormatter(["0", "20", "30", "40%"]))
    top.xaxis.set_minor_locator(NullLocator())
    for sp in ("left", "right", "bottom"):
        top.spines[sp].set_visible(False)
    top.spines["top"].set_visible(True); top.spines["top"].set_color(MUTED)
    top.tick_params(axis="x", labelsize=6, colors=MUTED, length=2, pad=1.2)
    top.set_xlabel("Risk reduction for typed, $1-1/\\eta$", fontsize=6, color=MUTED,
                   labelpad=2.5)
    letter(fig, 0.02, "a")

    # ---- b: protected-count sweep ----------------------------------------
    b = ax_in(fig, 4.62, B, 2.45, H)
    sc = 1e4
    t = V["target"] * sc
    tlo, thi = (v * sc for v in V["target_ci"])
    xmin, xmax = 109, 213
    b.set_yscale("log")
    b.set_ylim(0.62, 24)
    # shallow reference and its interval (hairlines, no filled box)
    b.axhline(t, color=SHALLOW, lw=0.9, ls=(0, (3.6, 1.8)), zorder=2)
    for v in (tlo, thi):
        b.axhline(v, color=SHALLOW, lw=0.45, ls=(0, (0.6, 1.4)), alpha=0.9, zorder=1)
    b.plot(168, t, "s", ms=4.0, mfc="white", mec=SHALLOW_DARK, mew=0.9, zorder=6)
    ax_text(xmax - 1, thi * 1.07, "<shallow> band, 168 sites", ax=b, ha="right",
            va="bottom", fontsize=6, color=MUTED,
            highlight_textprops=[dict(color=SHALLOW_DARK, fontweight="bold")])
    # infeasible region
    b.axvspan(xmin, 117, facecolor="white", edgecolor=HAIR, hatch="//////", lw=0, zorder=0)
    b.axvline(117, color=FAINT, lw=0.5, zorder=1)
    b.text((xmin + 117) / 2, 1.35, "not possible", rotation=90, ha="center", va="center",
           fontsize=5.3, color=MUTED, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    # typed sweep
    ys = np.array(V["typed_excess"]) * sc
    cis = np.array(V["typed_excess_ci"]) * sc
    b.plot(counts, ys, color=INK, lw=0.85, zorder=3, label="typed")
    for x, y, (lo, hi) in zip(counts, ys, cis):
        errpt(b, x, y, lo, hi, INK, filled=True, vertical=True, ms=3.5, lw=0.8)
    b.text(186, 1.86, "typed", ha="center", va="bottom", fontsize=6, fontweight="bold",
           color=INK)
    # crossing on the reference line, with its interval along the line
    n, (nlo, nhi) = V["n_star"], V["n_star_ci"]
    b.plot([nlo, nhi], [t, t], color=INK, lw=1.7, solid_capstyle="butt", zorder=5)
    b.plot(n, t, "D", ms=4.2, mfc="white", mec=INK, mew=1.0, zorder=7)
    # saving as a dimension arrow from the crossing to the shallow band
    ya = 15.0
    for xv in (n, 168):
        b.plot([xv, xv], [t * 1.12, ya * 0.93], color=MUTED, lw=0.4, zorder=2)
    b.annotate("", xy=(n, ya), xytext=(168, ya),
               arrowprops=dict(arrowstyle="<|-|>,head_length=0.32,head_width=0.16",
                               lw=0.7, color=INK, shrinkA=0, shrinkB=0))
    s, (slo, shi) = V["saving"], V["saving_ci"]
    b.text((n + 168) / 2, ya * 1.1, f"$\\Delta N$ = {s:.0f} ({slo:.0f}–{shi:.0f}) fewer",
           ha="center", va="bottom", fontsize=6, color=INK)
    b.set_xlim(xmin, xmax)
    b.set_xticks([117, 132, 150, 168, 204]); b.set_xticks([123], minor=True)
    b.set_yticks([1, 2, 5, 10]); b.set_yticklabels(["1", "2", "5", "10"])
    b.yaxis.set_minor_formatter(NullFormatter())
    hl_xlabel(b, "Protected ancillas, <typed> band", [INK])
    b.set_ylabel("Excess risk (10$^{-4}$)")
    letter(fig, 4.08, "b")
    return fig


# ---------------------------------------------------------------------------
CAPS = [2, 5, 10, 20]
NOCAP_X = 44.0
ROUND_US = 1.1


def axis_break(ax, x, size=0.022):
    tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for dx in (-0.55, 0.55):
        xs = x * 10 ** (dx * 0.03)
        ax.plot([xs * 10 ** (-size), xs * 10 ** size], [-0.035, 0.035], transform=tr,
                color=INK, lw=0.6, clip_on=False, zorder=10)
    ax.plot([x * 10 ** -0.02, x * 10 ** 0.02], [0, 0], transform=tr, color="white",
            lw=1.6, clip_on=False, zorder=9)


def cap_panel(ax, data, ylim, yticks, show_x):
    null_line(ax, y=0)
    k = 1.08
    for key, col, mk in (("background", INK, "o"), ("oracle", ORACLE, "s")):
        f = 1 / k if key == "background" else k
        xs = np.array(CAPS + [NOCAP_X], float) * f
        pts = data[key]
        ax.plot(xs[:4], [q[0] for q in pts[:4]], color=col, lw=0.55, alpha=0.5, zorder=3)
        for x, (p, lo, hi) in zip(xs, pts):
            errpt(ax, x, p, lo, hi, col, filled=lo > 0 or hi < 0, vertical=True,
                  marker=mk, ms=3.6 if mk == "o" else 3.2, lw=0.8)
    ax.set_xscale("log"); ax.set_xlim(1.6, 58)
    ax.xaxis.set_major_locator(FixedLocator(CAPS + [NOCAP_X]))
    ax.xaxis.set_major_formatter(FixedFormatter([str(c) for c in CAPS] + ["No cap"]
                                                if show_x else [""] * 5))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(*ylim); ax.set_yticks(yticks)
    axis_break(ax, 30)


def fig5():
    fig = plt.figure(figsize=(W, 3.25))

    # ---- a: d = 13 ratios, d = 9 ghosts ----------------------------------
    a = ax_in(fig, 1.30, 0.42, 2.20, 2.42)
    for y, (p9, _, _) in zip(ROW_Y, D["fig2a"]):
        a.plot(p9, y, "D", ms=3.0, mfc=GHOST, mec="white", mew=0.4, zorder=4.6)
    for y, (p, lo, hi), (_, _, prim) in zip(ROW_Y, D["fig5a"], ROWS):
        errpt(a, p, y, lo, hi, DEEP, filled=lo > 1 or hi < 1, ms=4.4 if prim else 3.6,
              lw=1.15 if prim else 0.8)
    log_ratio_axis(a, 0.4, 9.0, [0.5, 1, 2, 4, 8])
    null_line(a, x=1.0)
    decoder_rows(a)
    hl_xlabel(a, "<Deep> / <typed> excess risk, $\\eta$, $d$ = 13", [DEEP, INK])
    direction_cue(a, 1.0, "deep better", "typed better", y=1.05)
    # key
    a.plot([0.66], [0.06], "D", ms=2.9, mfc=GHOST, mec="none", transform=a.transAxes)
    a.text(0.69, 0.06, "$d$ = 9 (Fig. 2a)", transform=a.transAxes, fontsize=5.6,
           color=MUTED, va="center")
    letter(fig, 0.02, "a")

    # ---- b, c: burst-duration caps ---------------------------------------
    b = ax_in(fig, 4.70, 1.86, 2.37, 0.98)
    c = ax_in(fig, 4.70, 0.42, 2.37, 0.98)
    cap_panel(b, D["fig5b"], (-2.7, 9.6), [0, 4, 8], show_x=False)
    cap_panel(c, D["fig5c"], (-1.6, 4.3), [-1, 0, 1, 2, 3, 4], show_x=True)
    for ax, txt, cols in ((b, "<Shallow> − <typed>", [SHALLOW, INK]),
                          (c, "<Deep> − <typed>", [DEEP, INK])):
        ax_text(0.02, 0.97, txt, ax=ax, ha="left", va="top", fontsize=6.5, color=INK,
                highlight_textprops=[dict(color=cols[0], fontweight="bold"),
                                     dict(color=cols[1], fontweight="bold")],
                annotationbbox_kw=dict(xycoords=ax.transAxes))
        ax.set_ylabel("$\\Delta$ excess risk (10$^{-4}$)", fontsize=6.5)
    c.set_xlabel("Burst-duration cap (rounds)")
    # top axis in microseconds
    top = b.secondary_xaxis("top", functions=(lambda r: r * ROUND_US, lambda u: u / ROUND_US))
    top.xaxis.set_major_locator(FixedLocator([c_ * ROUND_US for c_ in CAPS]))
    top.xaxis.set_major_formatter(FixedFormatter(["2.2", "5.5", "11", "22 µs"]))
    top.xaxis.set_minor_locator(NullLocator())
    top.tick_params(axis="x", labelsize=6, colors=MUTED, length=2, pad=1.2)
    top.spines["top"].set_color(MUTED)
    # marker key (states the filled/open rule)
    kx, ky = 0.035, 0.66
    b.plot([kx], [ky], "o", ms=3.6, color=INK, transform=b.transAxes)
    b.text(kx + 0.03, ky, "Background matching", transform=b.transAxes, fontsize=5.8,
           va="center")
    b.plot([kx], [ky - 0.13], "s", ms=3.2, color=ORACLE, transform=b.transAxes)
    b.text(kx + 0.03, ky - 0.13, "Event-informed oracle", transform=b.transAxes,
           fontsize=5.8, va="center", color=ORACLE)
    b.plot([kx - 0.008, kx + 0.008], [ky - 0.26] * 2, transform=b.transAxes)
    b.plot([kx], [ky - 0.26], "o", ms=3.4, mfc="white", mec=MUTED, mew=0.7,
           transform=b.transAxes)
    b.text(kx + 0.03, ky - 0.26, "open: 95% CI includes 0", transform=b.transAxes,
           fontsize=5.6, va="center", color=MUTED)
    # five-round reversal
    p5 = D["fig5b"]["background"][1][0]
    b.annotate("reversal", xy=(5 / 1.08, p5), xytext=(2.35, -1.95), fontsize=5.8,
               color=INK, va="center", ha="left",
               arrowprops=dict(arrowstyle="-", lw=0.45, color=MUTED, shrinkA=0.5,
                               shrinkB=2.5))
    letter(fig, 4.12, "b")
    letter(fig, 4.12, "c", y_in=1.47)
    return fig


def cvd_preview(png, out_png):
    """Simulate 100% deuteranomaly and greyscale on the rendered figure."""
    from colorspacious import cspace_convert
    import matplotlib.image as mpimg
    img = mpimg.imread(png)[..., :3]
    deut = np.clip(cspace_convert(img, {"name": "sRGB1+CVD", "cvd_type": "deuteranomaly",
                                        "severity": 100}, "sRGB1"), 0, 1)
    grey = np.repeat(np.dot(img, [0.2126, 0.7152, 0.0722])[..., None], 3, axis=2)
    stack = np.concatenate([img, deut, grey], axis=0)
    plt.imsave(out_png, stack)


if __name__ == "__main__":
    out = os.environ.get("FIGURE_OUT_DIR", os.path.join(REPO, "out", "figures"))
    os.makedirs(out, exist_ok=True)
    jobs = (("fig02_placement_resources", lambda: fig2(sys.argv[1])),
            ("fig05_decoder_duration", fig5))
    for name, fn in jobs:
        f = fn()
        pdf = os.path.join(out, f"{name}.pdf")
        f.savefig(pdf)
        f.savefig(os.path.join(out, f"{name}.png"), dpi=300)
        svg = os.path.join(out, f"{name}.svg")
        f.savefig(svg)
        S._postprocess(pdf, svg, name)
        plt.close(f)
        cvd_preview(os.path.join(out, f"{name}.png"), os.path.join(out, f"{name}_cvd.png"))
    print("done")
