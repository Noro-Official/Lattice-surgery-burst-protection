"""
PACER manuscript, Figures 2-5. Publication styling after Google Quantum AI
Nature / Science papers.

Layout system
- All figures are 7.0 in wide (double column, revtex / Nature 183 mm).
- Axes are placed in absolute inches (ax_in), so panels share left edges,
  baselines and panel-letter positions across figures when stacked in the PDF.
- Forest plots (2a, 2b, 5a) use one row grid: decoder family (Background /
  Oracle) x matching (plain / correlated). Fig 2b reuses 2a's rows directly.

Encoding (constant across the paper)
- deep = violet, shallow = green, typed = black reference line at the null
- filled marker = 95% interval excludes the null; open = includes it
- oracle rows sit on a pale band: they need the true simulated impact

DATA PROVENANCE
EXACT     = taken from the manuscript text / tables
DIGITIZED = read off the submitted figures (about +/-0.1 units)
Replace every DIGITIZED block with the frozen summaries before submission.
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle
from matplotlib.ticker import NullFormatter, NullLocator
import matplotlib.patheffects as pe

# Open-source building blocks
#   SciencePlots  (github.com/garrettj403/SciencePlots)   - 'nature' base style
#   cmcrameri     (github.com/callumrollo/cmcrameri)      - Crameri perceptually
#                                                           uniform, CVD-safe maps
#   labellines    (github.com/cphyc/matplotlib-label-lines) - inline curve labels
import scienceplots  # noqa: F401  (registers styles)
import cmcrameri.cm as cmc
from labellines import labelLines
#   highlight_text (github.com/znstrider/highlight_text)   - colour-coded words
#                                                           in axis labels
#   URW base35 Nimbus Sans (github.com/ArtifexSoftware/urw-base35-fonts)
#                                                         - metric Helvetica clone
#   colorspacious (github.com/njsmith/colorspacious)      - CVD check (see below)
from highlight_text import ax_text
from matplotlib import font_manager

plt.style.use(["nature"])

_FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts")
if os.path.isdir(_FONT_DIR):
    for _f in os.listdir(_FONT_DIR):
        if _f.endswith(".ttf"):
            font_manager.fontManager.addfont(os.path.join(_FONT_DIR, _f))

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans",
                        "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.35, "ytick.minor.width": 0.35,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "xtick.minor.size": 1.4, "ytick.minor.size": 1.4,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#262626", "xtick.color": "#262626",
    "ytick.color": "#262626", "axes.labelcolor": "#1a1a1a",
    "text.color": "#1a1a1a", "axes.labelpad": 2.5,
    "lines.linewidth": 0.9, "lines.solid_capstyle": "round",
    "legend.frameon": False, "hatch.linewidth": 0.35,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "sans",
    "mathtext.it": "sans:italic", "mathtext.bf": "sans:bold",
    "svg.fonttype": "none",
})

INK, MUTED, FAINT, HAIR = "#1a1a1a", "#707070", "#b8b8b8", "#d9d9d9"
BAND = "#f3f3f1"
DEEP, DEEP_GHOST = "#5B2C91", "#c9b6e2"
SHALLOW, SHALLOW_DARK, SHALLOW_TINT = "#23894D", "#16603A", "#eaf5ee"
ORACLE = "#2C67C4"
D13 = "#9A78C9"            # d = 13 companion to DEEP
def _ramp(cmap, lo, hi, n):
    return [mpl.colors.to_hex(cmap(v)) for v in np.linspace(lo, hi, n)]


ZZ = _ramp(cmc.oslo, 0.66, 0.24, 3)        # light -> dark blue
XX = _ramp(cmc.lajolla, 0.74, 0.22, 4)     # amber -> deep brown-red
DATA_Q = "#c9a227"                          # data-qubit gold, as in device schematics

W = 183 / 25.4    # Nature double column, 7.20 in
LETTER_SIZE = 8


# Palette check (colorspacious, CAM02-UCS, 100% severity): worst distinct-pair
# dE = 13.0 deuteranomaly, 10.6 protanomaly, 10.0 tritanomaly  (all > 8).


def hl_xlabel(ax, text, colors, pad_in=0.19, size=7):
    """Axis label with colour-coded words: '<Deep> / typed excess risk'."""
    fh = ax.figure.get_size_inches()[1]
    h = ax.get_position().height * fh
    props = [dict(color=c, fontweight="bold") for c in colors]
    ax_text(0.5, -pad_in / h, text, ax=ax, ha="center", va="top",
            highlight_textprops=props, fontsize=size, color=INK,
            annotationbbox_kw=dict(xycoords=ax.transAxes))


def ax_in(fig, left, bottom, width, height, **kw):
    fw, fh = fig.get_size_inches()
    return fig.add_axes([left / fw, bottom / fh, width / fw, height / fh], **kw)


def letter(fig, x_in, ch, y_in=None):
    fw, fh = fig.get_size_inches()
    y = (1 - 0.04 / fh) if y_in is None else y_in / fh
    fig.text(x_in / fw, y, ch, fontsize=LETTER_SIZE,
             fontweight="bold", va="top", ha="left")


def errpt(ax, x, y, lo, hi, color, filled=True, vertical=False,
          marker="o", ms=3.8, lw=0.85, z=4):
    if vertical:
        ax.plot([x, x], [lo, hi], color=color, lw=lw, zorder=z, solid_capstyle="butt")
    else:
        ax.plot([lo, hi], [y, y], color=color, lw=lw, zorder=z, solid_capstyle="butt")
    ax.plot(x, y, marker=marker, ms=ms, mfc=color if filled else "white",
            mec=color, mew=0.8, ls="none", zorder=z + 1)


def null_line(ax, x=None, y=None):
    if x is not None:
        ax.axvline(x, color=INK, lw=0.6, zorder=2)
    else:
        ax.axhline(y, color=INK, lw=0.6, zorder=2)


def direction_cue(ax, null, left, right, y=1.035):
    """Small 'which side is better' cue above an axis."""
    tr = ax.get_xaxis_transform()
    kw = dict(arrowstyle="-|>,head_length=0.35,head_width=0.18", lw=0.45,
              color=MUTED, shrinkA=0, shrinkB=0)
    x0, x1 = ax.get_xlim()
    gap = (np.log10(x1 / x0) if ax.get_xscale() == "log" else x1 - x0) * 0.012

    def shift(v, d):
        return v * 10 ** d if ax.get_xscale() == "log" else v + d
    if left:
        ax.annotate("", xy=(shift(null, -gap * 9), y), xytext=(shift(null, -gap), y),
                    xycoords=tr, textcoords=tr, arrowprops=kw)
        ax.text(shift(null, -gap * 10), y, left + " ", transform=tr, ha="right",
                va="center", fontsize=6, color=MUTED)
    if right:
        ax.annotate("", xy=(shift(null, gap * 9), y), xytext=(shift(null, gap), y),
                    xycoords=tr, textcoords=tr, arrowprops=kw)
        ax.text(shift(null, gap * 10), y, " " + right, transform=tr, ha="left",
                va="center", fontsize=6, color=MUTED)


# ---------------------------------------------------------------------------
# Layout pictograms (schematic; equal strong count = 16 in each)
# ---------------------------------------------------------------------------
X_WEAK, Z_WEAK = "#d6604d", "#4a7fc1"
X_STRONG, Z_STRONG = "#f7dcd6", "#dbe6f5"


def _strong(layout, is_x, dist):
    if layout == "deep":
        return dist <= 2
    if layout == "typed":
        return is_x
    if layout == "shallow":
        return dist <= 1 or dist == 4
    raise ValueError(layout)


def pictogram(fig, left, bottom, layout, w=0.46, h=0.19):
    ax = ax_in(fig, left, bottom, w, h)
    rows, cols = 3, 8
    for r in range(rows):
        for c in range(cols):
            is_x = (r + c) % 2 == 0
            dist = int(abs(c - 3.5) + 0.5)
            strong = _strong(layout, is_x, dist)
            col = (X_STRONG if is_x else Z_STRONG) if strong else (X_WEAK if is_x else Z_WEAK)
            ax.add_patch(Rectangle((c + 0.08, r + 0.08), 0.84, 0.84, facecolor=col,
                                   edgecolor="none"))
    ax.plot([4, 4], [-0.2, rows + 0.2], color=INK, lw=0.8)
    ax.set_xlim(0, cols); ax.set_ylim(-0.2, rows + 0.2)
    ax.set_aspect("auto"); ax.axis("off")
    return ax


LAYOUT_COLOR = {"typed": INK, "deep": DEEP, "shallow": SHALLOW}


def layout_key(fig, y_in, x_start, layouts=("typed", "deep", "shallow"), gap=1.45):
    fw, fh = fig.get_size_inches()
    for i, lay in enumerate(layouts):
        x = x_start + i * gap
        pictogram(fig, x, y_in, lay)
        fig.text((x + 0.53) / fw, (y_in + 0.145) / fh, lay.capitalize(),
                 fontsize=6.5, fontweight="bold", color=LAYOUT_COLOR[lay],
                 va="center")
        fig.text((x + 0.53) / fw, (y_in + 0.04) / fh,
                 {"typed": "measured type only", "deep": "both types, wide",
                  "shallow": "both types, narrow"}[lay],
                 fontsize=5.5, color=MUTED, va="center")
    # weak / strong swatches
    xs = x_start + len(layouts) * gap
    for j, (lab, cx, cz) in enumerate((("weak", X_WEAK, Z_WEAK),
                                       ("strong", X_STRONG, Z_STRONG))):
        yy = y_in + 0.145 - j * 0.105
        for k, cc in enumerate((cx, cz)):
            fig.patches.append(Rectangle(((xs + k * 0.075) / fw, (yy - 0.035) / fh),
                                         0.065 / fw, 0.07 / fh, facecolor=cc,
                                         edgecolor="none", transform=fig.transFigure))
        fig.text((xs + 0.18) / fw, yy / fh, lab, fontsize=5.5, color=MUTED,
                 va="center")


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
# Row order for all decoder forest plots (top to bottom)
ROWS = [("Background", "plain"), ("Background", "correlated"),
        ("Oracle", "plain"), ("Oracle", "correlated")]
ROW_Y = np.array([3.0, 2.0, 0.8, -0.2])     # small gap between families

FIG2A = [(1.36, 1.28, 1.45), (1.62, 1.49, 1.77),       # EXACT
         (1.28, 1.21, 1.36), (1.43, 1.32, 1.56)]
FIG2B = [(15, -1, 36, False), (36, 28, 36, True),       # plain rows DIGITIZED
         (11, -11, 36, False), (36, 9, 36, True)]       # corr. rows EXACT (>=36; lower 28, 9)
FIG2B_RESOLVED = [False, True, False, True]

FIG3A = {"5": [-7.6, -2.0, 0, 0, 0, 0, 0, 0],                       # DIGITIZED
         "9": [-16.04, -10.50, -4.96, 0, 0, 0, 0, 0],               # EXACT
         "13": [-16.9, -11.3, -6.2, -1.4, 0, 0, 0, 0]}              # DIGITIZED
FIG3B = {7: [-16.5, -16.5, -16.1, -10.6, -5.1, 0],                  # DIGITIZED
         9: [-20.72, -20.72, -20.32, -20.32, -14.78, -9.24, -3.71, 0],  # EXACT
         11: [-25.0, -25.0, -24.6, -24.6, -24.5, -19.0, -13.5, -8.0, -2.5, 0],
         13: [-29.2, -29.2, -28.8, -28.8, -28.8, -28.5, -22.5, -17.3, -11.8,
              -6.3, -0.5, 0]}
FIG3C = {"All channels": [-14.9, -9.2, -3.9, -1.1, -1.1, -1.1, -1.1, -1.1],
         "Reports only": [-12.3, -7.3, -1.5, 0, 0, 0, 0, 0],
         "Pauli + reset": [-11.5, -6.8, -1.3, -1.1, -1.1, -1.1, -1.1, -1.1]}
# depth-1 values EXACT; remainder DIGITIZED

FIG4A = [(0, 0, 0, 0.0, 0.0, 0.0), (0, 0, 1, -0.2, -0.7, 0.3),     # DIGITIZED
         (0, 1, 0, 0.05, -0.7, 0.8), (0, 1, 1, -0.2, -0.8, 0.4),
         (1, 0, 0, 2.2, 1.2, 3.3), (1, 0, 1, 2.0, 0.8, 3.5),
         (1, 1, 0, 4.1, 2.3, 6.3), (1, 1, 1, 4.6, 2.5, 7.0)]
FIG4B = [("Report", 3.33, 1.9, 5.0), ("Ancilla Pauli", 1.15, 0.55, 1.8),
         ("Reset", -0.05, -0.65, 0.65),
         ("Report × ancilla", 1.12, 0.71, 1.60),
         ("Report × reset", 0.20, -0.05, 0.50),
         ("Ancilla × reset", 0.12, -0.45, 0.75)]
# Report / Ancilla main effects and Report x ancilla EXACT; other CIs DIGITIZED

FIG5A = [(1.70, 1.38, 2.25), (3.30, 2.20, 5.60),       # background rows DIGITIZED
         (0.96, 0.62, 1.56), (1.12, 0.48, 2.65)]       # oracle rows EXACT
CAPS = np.array([2, 5, 10, 20], float)                 # all DIGITIZED
FIG5B_BG = [(0.3, -0.6, 0.9), (-0.9, -1.7, -0.1), (0.6, -0.4, 1.6), (2.5, 0.4, 4.6)]
FIG5B_OR = [(0.3, -0.35, 0.8), (-0.8, -1.6, 0.15), (0.3, -0.6, 1.1), (1.2, -0.1, 2.5)]
ROUND_US = 1.1


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------
def decoder_rows(ax, labels=True):
    ax.set_ylim(-0.75, 3.55)
    ax.set_yticks(ROW_Y)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    if labels:
        ax.set_yticklabels([m for _, m in ROWS], color=INK)
        fig = ax.figure
        fw = fig.get_size_inches()[0]
        x0 = ax.get_position().x0 - 1.38 / fw       # family column
        import matplotlib.transforms as mtr
        tr = mtr.blended_transform_factory(fig.transFigure, ax.transData)
        for fam, yc in (("Background", 2.5), ("Oracle", 0.55)):
            ax.text(x0, yc, fam, transform=tr, ha="left", va="center",
                    fontsize=7, fontweight="bold")
        ax.text(x0, -0.05, "uses true impact", transform=tr, ha="left",
                va="center", fontsize=5.5, color=MUTED, style="italic")
    else:
        ax.set_yticklabels([])


def value_column(fig, ax, data, x_in, header):
    fw, fh = fig.get_size_inches()
    for y, (p, lo, hi) in zip(ROW_Y, data):
        res = lo > 1 or hi < 1
        disp = ax.transData.transform((1, y))
        yf = fig.transFigure.inverted().transform(disp)[1]
        fig.text(x_in / fw, yf, f"{p:.2f}", ha="right", va="center", fontsize=6.5,
                 color=INK if res else MUTED)
        fig.text((x_in + 0.04) / fw, yf, f"[{lo:.2f}, {hi:.2f}]", ha="left",
                 va="center", fontsize=6, color=MUTED)
    top = fig.transFigure.inverted().transform(ax.transAxes.transform((0, 1)))[1]
    fig.text((x_in + 0.04) / fw, top, header, ha="center", va="bottom",
             fontsize=6, color=MUTED)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
# Stim (github.com/quantumlib/Stim) supplies the rotated-surface-code geometry:
# qubit coordinates and which measure qubits are X-type (Hadamard-sandwiched).
import stim
from matplotlib.patches import Polygon, Wedge


def _stim_patch(d=5):
    c = stim.Circuit.generated("surface_code:rotated_memory_z", distance=d, rounds=1)
    co = c.get_final_qubit_coordinates()
    mr, hq, mq = set(), set(), set()
    for ins in c.flattened():
        tg = [t.value for t in ins.targets_copy() if t.is_qubit_target]
        if ins.name == "MR":
            mr |= set(tg)
        elif ins.name == "H":
            hq |= set(tg)
        elif ins.name == "M":
            mq |= set(tg)
    data = [tuple(co[q]) for q in mq]
    meas = [(tuple(co[q]), q in hq) for q in mr]
    return data, meas


_PATCH = _stim_patch(5)


def big_pictogram(ax, layout):
    """Two d=5 rotated surface-code patches (geometry from Stim) about a merge
    seam. Weak = saturated plaquette, strong = pale. Schematic budget."""
    data, meas = _PATCH
    seam = 11.0
    for off, side in ((0.0, -1), (12.0, +1)):
        dset = {(x + off, y) for x, y in data}
        cx = 5.0 + off
        for (mx, my), is_x in meas:
            mx += off
            nb = [(mx + a, my + b) for a in (-1, 1) for b in (-1, 1)
                  if (mx + a, my + b) in dset]
            dist = int(abs(mx - seam) // 2) + 1
            if layout == "deep":
                strong = dist <= 2
            elif layout == "typed":
                strong = is_x and dist <= 4
            else:
                strong = dist <= 1
            fc = (X_STRONG if is_x else Z_STRONG) if strong else (X_WEAK if is_x else Z_WEAK)
            if len(nb) == 4:
                ax.add_patch(Polygon([(mx - 1, my - 1), (mx + 1, my - 1),
                                      (mx + 1, my + 1), (mx - 1, my + 1)],
                                     closed=True, facecolor=fc, edgecolor="white",
                                     lw=0.35, zorder=1))
            elif len(nb) == 2:
                (x1, y1), (x2, y2) = nb
                ccx, ccy = (x1 + x2) / 2, (y1 + y2) / 2
                if x1 == x2:                      # vertical pair -> bulge in x
                    th = 90 if mx > ccx else -90
                    th1, th2 = (-90, 90) if mx > ccx else (90, 270)
                else:                             # horizontal pair -> bulge in y
                    th1, th2 = (0, 180) if my > ccy else (180, 360)
                ax.add_patch(Wedge((ccx, ccy), 1.0, th1, th2, facecolor=fc,
                                   edgecolor="white", lw=0.35, zorder=1))
        for x, y in dset:
            ax.add_patch(Circle((x, y), 0.32, facecolor=DATA_Q, edgecolor="white",
                                lw=0.3, zorder=3))
    ax.plot([seam, seam], [-0.6, 10.6], color=INK, lw=1.0, zorder=4,
            solid_capstyle="round")
    ax.set_xlim(-1.2, 23.2); ax.set_ylim(-1.2, 11.2)
    ax.set_aspect("equal"); ax.axis("off")


def fig2():
    """Placement advantage and its limits (merges former Figs 2 and 5)."""
    fig = plt.figure(figsize=(W, 3.75))
    fw, fh = fig.get_size_inches()

    # ---- a: layout schematic -------------------------------------------
    top = 3.30
    for i, lay in enumerate(("typed", "deep", "shallow")):
        y = top - i * 0.50
        ax = ax_in(fig, 0.06, y - 0.05, 0.90, 0.44)
        big_pictogram(ax, lay)
        fig.text(1.00 / fw, (y + 0.22) / fh, lay.capitalize(), fontsize=7,
                 fontweight="bold", color=LAYOUT_COLOR[lay], va="center")
        fig.text(1.00 / fw, (y + 0.09) / fh,
                 {"typed": "measured type only",
                  "deep": "both types, widest band",
                  "shallow": "both types, fixed depth"}[lay],
                 fontsize=6, color=MUTED, va="center")
    # swatches
    ky = top - 1.5 + 0.18
    for j, (lab, cx, cz) in enumerate((("weak", X_WEAK, Z_WEAK),
                                       ("strong", X_STRONG, Z_STRONG))):
        x0 = 0.12 + j * 0.62
        for k, cc in enumerate((cx, cz)):
            fig.patches.append(Rectangle(((x0 + k * 0.085) / fw, (ky - 0.04) / fh),
                                         0.075 / fw, 0.08 / fh, facecolor=cc,
                                         edgecolor="none", transform=fig.transFigure))
        fig.text((x0 + 0.20) / fw, ky / fh, lab, fontsize=6, color=MUTED, va="center")
    fig.text(1.36 / fw, ky / fh, "│ merge seam", fontsize=6, color=MUTED, va="center")
    letter(fig, 0.02, "a")

    # ---- b: deep/typed at d = 9 and d = 13 -----------------------------
    b = ax_in(fig, 3.95, 2.18, 3.10, 1.25)
    dy = 0.17
    for y, (p, lo, hi) in zip(ROW_Y, FIG2A):
        errpt(b, p, y + dy, lo, hi, DEEP, filled=lo > 1, ms=3.6)
    for y, (p, lo, hi) in zip(ROW_Y, FIG5A):
        errpt(b, p, y - dy, lo, hi, D13, filled=lo > 1 or hi < 1,
              marker="D", ms=3.0)
    b.set_xscale("log"); b.set_xlim(0.42, 7.0)
    b.set_xticks([0.5, 1, 2, 4]); b.set_xticklabels(["0.5", "1", "2", "4"])
    b.xaxis.set_minor_formatter(NullFormatter())
    null_line(b, x=1.0)
    decoder_rows(b)
    hl_xlabel(b, "<Deep> / <typed> excess risk", [DEEP, INK])
    direction_cue(b, 1.0, "deep better", "typed better", y=1.08)
    b.plot([4.0], [3.40], "o", ms=3.4, color=DEEP, clip_on=False)
    b.text(4.45, 3.40, "$d$ = 9", fontsize=6, va="center")
    b.plot([4.0], [2.98], "D", ms=2.8, color=D13, clip_on=False)
    b.text(4.45, 2.98, "$d$ = 13", fontsize=6, va="center", color=D13)
    letter(fig, 2.52, "b")

    # ---- c: protection-count saving ------------------------------------
    c = ax_in(fig, 1.48, 0.38, 2.05, 1.18)
    c.add_patch(Rectangle((36, -0.75), 14, 4.3, facecolor="white",
                          edgecolor=HAIR, hatch="//////", lw=0, zorder=0.6))
    c.axvline(36, color=FAINT, lw=0.5, ls=(0, (2, 1.5)), zorder=1)
    for y, (p, lo, hi, cens), res in zip(ROW_Y, FIG2B, FIG2B_RESOLVED):
        errpt(c, p, y, lo, hi, SHALLOW, filled=res)
        if cens:
            c.annotate("", xy=(46, y), xytext=(37.2, y),
                       arrowprops=dict(arrowstyle="-|>,head_length=0.4,head_width=0.2",
                                       lw=0.85, color=SHALLOW, shrinkA=0, shrinkB=0),
                       zorder=5)
            c.text(41.5, y + 0.2, "≥36", ha="center", va="bottom", fontsize=6)
    null_line(c, x=0)
    c.set_xlim(-15, 50); c.set_xticks([-12, 0, 12, 24, 36])
    decoder_rows(c)
    c.text(43, 3.5, "beyond\ngrid", ha="center", va="top", fontsize=5.5,
           color=MUTED, linespacing=0.95,
           bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    hl_xlabel(c, "Protected ancillas saved, <typed> vs <shallow>", [INK, SHALLOW])
    direction_cue(c, 0, "", "typed needs fewer", y=1.08)
    letter(fig, 0.02, "c", y_in=1.86)

    # ---- d: duration cap -------------------------------------------------
    d = ax_in(fig, 4.55, 0.38, 2.50, 1.18)
    null_line(d, y=0)
    k = 1.075
    for ser, col, mk, f in ((FIG5B_BG, INK, "o", 1 / k), (FIG5B_OR, ORACLE, "s", k)):
        xs = CAPS * f
        d.plot(xs, [q[0] for q in ser], color=col, lw=0.5, alpha=0.45, zorder=3)
        for x, (p, lo, hi) in zip(xs, ser):
            errpt(d, x, p, lo, hi, col, filled=lo > 0 or hi < 0, vertical=True,
                  marker=mk, ms=3.4 if mk == "o" else 3.0)
    d.set_xscale("log"); d.set_xlim(1.6, 30)
    d.set_xticks(CAPS); d.set_xticklabels([f"{int(v)}" for v in CAPS])
    d.xaxis.set_minor_locator(NullLocator())
    d.set_ylim(-2.3, 5.6); d.set_yticks([-2, 0, 2, 4])
    d.plot([1.85], [4.45], "o", ms=3.4, color=INK)
    d.text(2.05, 4.45, "Background", fontsize=6, va="center")
    d.plot([1.85], [3.75], "s", ms=3.0, mfc="white", mec=ORACLE, mew=0.8)
    d.text(2.05, 3.75, "Oracle", fontsize=6, va="center", color=ORACLE)
    d.text(28.5, 5.45, "typed better $\\uparrow$", fontsize=5.5, color=MUTED,
           va="top", ha="right")
    d.text(28.5, -2.15, "shallow better $\\downarrow$", fontsize=5.5, color=MUTED,
           va="bottom", ha="right")
    d.set_xlabel("Burst-duration cap (rounds of 1.1 µs)")
    d.set_ylabel("Shallow − typed (10$^{-4}$)")
    letter(fig, 3.95, "d", y_in=1.86)
    return fig


def fig3():
    """Mechanism: schedule schematic + failing-path costs."""
    fig = plt.figure(figsize=(W, 2.75))
    fw, fh = fig.get_size_inches()

    # ---- a: CNOT schedule (d = 9) ----------------------------------------
    t = ax_in(fig, 0.62, 2.02, 6.43, 0.52)
    pre, dur, post = 12, 9, 12
    zz0, zz1 = pre, pre + dur
    xx0, xx1 = zz1 + dur, zz1 + 2 * dur
    end = xx1 + post
    rows = {"C": 2, "A": 1, "T": 0}
    for name, yv in rows.items():
        stop = xx1 if name == "A" else end
        t.plot([0, stop], [yv, yv], color="#cfcfcf", lw=2.4, solid_capstyle="round",
               zorder=1)
        t.text(-0.8, yv, name, ha="right", va="center", fontsize=6.5,
               fontweight="bold")
    t.add_patch(FancyBboxPatch((zz0, 0.72), dur, 1.56,
                               boxstyle="round,pad=0,rounding_size=0.9",
                               mutation_aspect=0.12, facecolor=Z_STRONG,
                               edgecolor=ZZ[1], lw=0.6, zorder=2))
    t.text((zz0 + zz1) / 2, 1.5, "$ZZ$ merge", ha="center", va="center",
           fontsize=6, color=ZZ[2], zorder=3)
    t.add_patch(FancyBboxPatch((xx0, -0.28), dur, 1.56,
                               boxstyle="round,pad=0,rounding_size=0.9",
                               mutation_aspect=0.12, facecolor=X_STRONG,
                               edgecolor=XX[1], lw=0.6, zorder=2))
    t.text((xx0 + xx1) / 2, 0.5, "$XX$ merge", ha="center", va="center",
           fontsize=6, color=XX[3], zorder=3)
    # A readout = temporal boundary
    t.plot([xx1, xx1], [0.62, 1.38], color=INK, lw=1.1, zorder=4)
    t.plot(xx1, 1, marker="s", ms=4, mfc="white", mec=INK, mew=0.8, zorder=5)
    t.annotate("A read out in $Z$: no later round for its $X$ checks",
               xy=(xx1 + 0.4, 1.0), xytext=(xx1 - 6, 2.75), fontsize=5.5,
               color=MUTED, va="center", ha="left", linespacing=0.95,
               arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED,
                               shrinkA=0, shrinkB=1))
    for x0, x1, lab in ((0, zz0, "12 pre-rounds"), (zz1, xx0, "$d$ rounds"),
                        (xx1, end, "12 post-rounds")):
        t.text((x0 + x1) / 2, -0.75, lab, ha="center", va="center", fontsize=5.5,
               color=MUTED)
    t.annotate("", xy=(end + 1.5, -0.75), xytext=(end - 0.2, -0.75),
               arrowprops=dict(arrowstyle="-|>,head_length=0.3,head_width=0.15",
                               lw=0.45, color=MUTED))
    t.text(end + 1.8, -0.75, "time", fontsize=5.5, color=MUTED, va="center")
    t.set_xlim(-1.5, end + 5); t.set_ylim(-1.05, 3.05)
    t.axis("off")
    letter(fig, 0.02, "a")

    # ---- b-d: path costs -------------------------------------------------
    B, H, Wd = 0.38, 1.18, 1.72
    lefts = [0.62, 2.90, 5.15]
    a, b, c = (ax_in(fig, l, B, Wd, H) for l in lefts)

    def quiet(ax, label=False):
        ax.axhline(0, color=FAINT, lw=0.5, ls=(0, (2, 1.5)), zorder=0)
        if label:
            ax.text(ax.get_xlim()[1], 0.6, "quiet", ha="right", va="bottom",
                    fontsize=5.5, color=MUTED)

    def header(ax, txt):
        ax.text(0.0, 1.05, txt, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=6.5)

    for (kk, ys), col, mk in zip(FIG3A.items(), ZZ, ["^", "o", "s"]):
        a.plot(np.arange(1, len(ys) + 1), ys, color=col, marker=mk, ms=2.6, mew=0,
               lw=0.9, zorder=3)
    a.text(0.75, -7.6, "5", color=ZZ[0], fontsize=6, ha="right", va="center")
    a.text(0.75, -15.4, "9", color=ZZ[1], fontsize=6, ha="right", va="center")
    a.text(0.75, -17.4, "13", color=ZZ[2], fontsize=6, ha="right", va="center")
    a.text(0.15, -3.2, "merge\nrounds", fontsize=5.5, color=MUTED, ha="center",
           va="center", linespacing=0.95)
    a.set_xlim(-0.6, 8.4); a.set_xticks([1, 2, 4, 6, 8])
    a.set_ylim(-18.8, 2); a.set_yticks([-15, -10, -5, 0])
    quiet(a, label=True)
    header(a, "$ZZ$ merge, $d$ = 9")
    a.set_ylabel("Path-cost change, $W - W_{\\rm quiet}$")
    a.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 0.02, "b", y_in=1.80)

    xx_lines = []
    for (dd, ys), col in zip(FIG3B.items(), XX):
        ln, = b.plot(np.arange(1, len(ys) + 1), ys, color=col, marker="o", ms=2.2,
                     mew=0, lw=0.9, zorder=3, label=f"$d$ = {dd}")
        xx_lines.append(ln)
    labelLines(xx_lines, xvals=[4.5, 5.5, 7.5, 8.5], fontsize=5.5,
               align=False, outline_width=0, zorder=6,
               bbox=dict(facecolor="white", edgecolor="none", pad=0.4))
    b.annotate("plateau set by\nA readout", xy=(5.2, -28.6), xytext=(6.8, -30.2),
               fontsize=5.5, color=MUTED, va="center", ha="left", linespacing=0.95,
               arrowprops=dict(arrowstyle="-", lw=0.4, color=MUTED, shrinkA=1,
                               shrinkB=1.5))
    b.set_xlim(0.3, 12.6); b.set_xticks([1, 4, 8, 12])
    b.set_ylim(-33, 2.2); b.set_yticks([-30, -20, -10, 0])
    quiet(b)
    header(b, "$XX$ merge")
    b.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 2.40, "c", y_in=1.80)

    sty = {"All channels": (ZZ[2], "-", "o"),
           "Reports only": (ZZ[1], (0, (3.2, 1.4)), "s"),
           "Pauli + reset": (ZZ[0], (0, (1, 1.1)), "^")}
    for kk, ys in FIG3C.items():
        col, ls, mk = sty[kk]
        c.plot(np.arange(1, len(ys) + 1), ys, color=col, ls=ls, marker=mk, ms=2.6,
               mew=0, lw=0.95, label=kk, zorder=3)
    c.legend(loc="lower right", handlelength=2.6, borderaxespad=0.1,
             labelspacing=0.3, handletextpad=0.5, fontsize=6)
    c.set_xlim(0.4, 8.4); c.set_xticks([1, 2, 4, 6, 8])
    c.set_ylim(-18.8, 2); c.set_yticks([-15, -10, -5, 0])
    quiet(c)
    header(c, "Helper $Z$ checks, reference dose")
    c.set_xlabel("Same-type depth $\\ell$")
    letter(fig, 4.65, "d", y_in=1.80)
    return fig


def fig4():
    """Channel interventions (compact)."""
    fig = plt.figure(figsize=(W, 1.95))
    B, H = 0.36, 1.25
    ys = np.arange(8)[::-1].astype(float)
    m = ax_in(fig, 0.42, B, 0.50, H)
    a = ax_in(fig, 1.00, B, 2.95, H, sharey=m)
    b = ax_in(fig, 5.10, B, 1.95, H)
    for ax in (m, a):
        ax.axhline(3.5, color=HAIR, lw=0.6, zorder=0)
    for y, (r, an, p, *_r) in zip(ys, FIG4A):
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
    for y, (r, an, p, pt, lo, hi) in zip(ys, FIG4A):
        if (r, an, p) == (0, 0, 0):
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
    for y, (lab, pt, lo, hi), i in zip(ey, FIG4B, range(6)):
        main = i < 3
        errpt(b, pt, y, lo, hi, SHALLOW if main else SHALLOW_DARK,
              filled=lo > 0 or hi < 0, marker="o", ms=3.8)
    b.set_yticks(ey); b.set_yticklabels([f[0] for f in FIG4B])
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


def _postprocess(pdf_path, svg_path, name):
    """Non-matplotlib finishing steps.
    scour  (github.com/scour-project/scour): lossless SVG clean-up, keeps text
    pikepdf (github.com/pikepdf/pikepdf, qpdf): PDF metadata + font audit
    """
    try:
        from scour import scour as _sc
        opts = _sc.sanitizeOptions()
        opts.enable_viewboxing = True
        opts.strip_comments = True
        opts.remove_metadata = True
        opts.indent_type = "none"
        with open(svg_path, encoding="utf-8") as fh:
            cleaned = _sc.scourString(fh.read(), opts)
        with open(svg_path, "w", encoding="utf-8") as fh:
            fh.write(cleaned)
    except ImportError:
        pass
    try:
        import pikepdf
        with pikepdf.open(pdf_path, allow_overwriting_input=True) as pdf:
            with pdf.open_metadata() as meta:
                meta["dc:title"] = f"PACER {name}"
                meta["dc:creator"] = ["Eric Kim"]
            for page in pdf.pages:
                for _, font in page.Resources.get("/Font", {}).items():
                    desc = font.get("/DescendantFonts", [font])[0]
                    fd = desc.get("/FontDescriptor", {})
                    if not any(k in fd for k in ("/FontFile", "/FontFile2", "/FontFile3")):
                        print(f"WARNING {name}: font {desc.get('/BaseFont')} not embedded")
            pdf.save(pdf_path)
    except ImportError:
        pass


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "out")
    os.makedirs(out, exist_ok=True)
    for name, fn in (("fig2", fig2), ("fig3", fig3), ("fig4", fig4)):
        f = fn()
        f.savefig(os.path.join(out, f"{name}.pdf"))
        f.savefig(os.path.join(out, f"{name}.png"), dpi=600)
        svg = os.path.join(out, f"{name}.svg")
        f.savefig(svg)
        _postprocess(os.path.join(out, f"{name}.pdf"), svg, name)
        plt.close(f)
    print("done")
