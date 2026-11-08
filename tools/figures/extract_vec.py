"""Recover plotted values from the manuscript's vector figure PDFs.

For each axes frame, tick-label text positions define the pixel->data map
(linear or log). Small closed paths are markers; straight segments are
whiskers. Output is JSON per figure for the restyled plotting script.
"""
import json
import math
import re
import sys
import warnings

warnings.filterwarnings("ignore")
import pymupdf as fitz

NUM = re.compile(r"^[−\-]?\d+(\.\d+)?$")


def num(s):
    return float(s.replace("−", "-"))


def axes_frames(pg, dr):
    """Axes = rect frames, or spine pairs. Use bounding boxes of long spines."""
    frames = []
    for d in dr:
        r = d["rect"]
        it = d["items"]
        if len(it) == 1 and it[0][0] == "re" and 60 < r.width < 500 and 40 < r.height < 300:
            frames.append(fitz.Rect(r))
    # spines: long vertical line + long horizontal line sharing a corner
    vs = [d["rect"] for d in dr if len(d["items"]) == 1 and d["items"][0][0] == "l"
          and d["rect"].width < 0.5 and d["rect"].height > 40]
    hs = [d["rect"] for d in dr if len(d["items"]) == 1 and d["items"][0][0] == "l"
          and d["rect"].height < 0.5 and d["rect"].width > 60]
    for v in vs:
        for h in hs:
            if abs(v.x0 - h.x0) < 1.5 and abs(v.y1 - h.y0) < 1.5:
                frames.append(fitz.Rect(v.x0, v.y0, h.x1, v.y1))
    out = []
    for f in frames:
        if not any(abs(f.x0 - g.x0) < 2 and abs(f.y1 - g.y1) < 2 for g in out):
            out.append(f)
    return out


def tick_map(words, fr, axis):
    pts = []
    for w in words:
        x0, y0, x1, y1, t = w[:5]
        if not NUM.match(t):
            continue
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        if axis == "x" and fr.x0 - 8 <= cx <= fr.x1 + 8 and fr.y1 < cy < fr.y1 + 16:
            pts.append((cx, num(t)))
        if axis == "y" and fr.x0 - 34 < cx < fr.x0 and fr.y0 - 6 <= cy <= fr.y1 + 6:
            pts.append((cy, num(t)))
    pts = sorted(set(pts))
    if len(pts) < 2:
        return None
    vals = [v for _, v in pts]
    log = all(v > 0 for v in vals) and len(vals) >= 3 and \
        abs((vals[2] / vals[1]) - (vals[1] / vals[0])) < 1e-6 and vals[1] / vals[0] != 1 and \
        abs((vals[2] - vals[1]) - (vals[1] - vals[0])) > 1e-6
    tv = [math.log10(v) if log else v for v in vals]
    ps = [p for p, _ in pts]
    n = len(ps)
    mp, mv = sum(ps) / n, sum(tv) / n
    a = sum((p - mp) * (v - mv) for p, v in zip(ps, tv)) / sum((p - mp) ** 2 for p in ps)
    b = mv - a * mp
    resid = max(abs(a * p + b - v) for p, v in zip(ps, tv))
    return dict(a=a, b=b, log=log, resid=resid, ticks=vals)


def apply(m, p):
    v = m["a"] * p + m["b"]
    return 10 ** v if m["log"] else v


def extract(path):
    doc = fitz.open(path)
    pg = doc[0]
    dr = pg.get_drawings()
    words = pg.get_text("words")
    res = []
    for fr in sorted(axes_frames(pg, dr), key=lambda r: (round(r.y0 / 20), r.x0)):
        mx, my = tick_map(words, fr, "x"), tick_map(words, fr, "y")
        ax = dict(frame=[fr.x0, fr.y0, fr.x1, fr.y1], xmap=mx, ymap=my,
                  markers=[], segments=[])
        for d in dr:
            r = d["rect"]
            if not fr.intersects(r) and not fr.contains(r):
                continue
            kinds = [i[0] for i in d["items"]]
            if 2 <= r.width <= 8 and 2 <= r.height <= 8 and (
                    "c" in kinds or kinds == ["re"] or len(kinds) >= 3):
                cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
                shape = "circle" if "c" in kinds else ("square" if kinds == ["re"] or len(kinds) == 4 else "tri")
                if len(kinds) == 4 and "c" not in kinds:
                    # diamond vs square
                    p = [i[1] for i in d["items"]]
                    shape = "diamond" if abs(p[0].x - cx) < 0.3 or abs(p[0].y - cy) < 0.3 else "square"
                ax["markers"].append(dict(
                    px=cx, py=cy,
                    x=apply(mx, cx) if mx else None, y=apply(my, cy) if my else None,
                    fill=d.get("fill"), stroke=d.get("color"), shape=shape,
                    w=r.width))
            elif kinds == ["l"] and (r.width > 3 or r.height > 3) and \
                    fr.x0 - 1 <= r.x0 and r.x1 <= fr.x1 + 1 and fr.y0 - 1 <= r.y0 and r.y1 <= fr.y1 + 1:
                p0, p1 = d["items"][0][1], d["items"][0][2]
                ax["segments"].append(dict(
                    x0=apply(mx, p0.x) if mx else None, y0=apply(my, p0.y) if my else None,
                    x1=apply(mx, p1.x) if mx else None, y1=apply(my, p1.y) if my else None,
                    px0=p0.x, py0=p0.y, px1=p1.x, py1=p1.y,
                    color=d.get("color"), width=d.get("width")))
        res.append(ax)
    return res


if __name__ == "__main__":
    out = {}
    for p in sys.argv[2:]:
        out[p.rsplit("/", 1)[-1]] = extract(p)
    json.dump(out, open(sys.argv[1], "w"), indent=1, default=str)
    for k, axes in out.items():
        print("==", k)
        for i, a in enumerate(axes):
            print(f"  ax{i} frame={[round(v) for v in a['frame']]} "
                  f"x={a['xmap'] and (a['xmap']['ticks'], a['xmap']['log'], round(a['xmap']['resid'], 4))} "
                  f"y={a['ymap'] and (a['ymap']['ticks'], a['ymap']['log'], round(a['ymap']['resid'], 4))} "
                  f"markers={len(a['markers'])} segs={len(a['segments'])}")
