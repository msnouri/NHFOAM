"""Reading data from vector figures in a PDF (PyMuPDF): the plot frames,
their axes from the tick marks and tick labels, and the markers and lines
drawn in them.  Used by digitize_zijlema2011.py.

A frame is a stroked rectangle; its ticks are the short line segments that
start on its bottom and left edges, and the labels are the numbers printed
next to them.  The map from page coordinates to data values is a straight
line fitted to the labelled ticks, and the largest misfit is reported, so a
wrong pairing shows.
"""
import numpy as np
import pymupdf


def _num(s):
    s = s.replace("−", "-").replace("–", "-")
    try:
        return float(s)
    except ValueError:
        return None


class Panel:
    def __init__(self, page, rect, drawings, words, xscale=1.0, yscale=1.0, log=None):
        self.rect = rect
        x0, y0, x1, y1 = rect
        tx, ty = [], []
        for d in drawings:
            if len(d["items"]) != 1 or d["items"][0][0] != "l":
                continue
            p, q = d["items"][0][1], d["items"][0][2]
            if abs(p.x - q.x) < 0.01 and 0.5 < abs(p.y - q.y) < 5 and abs(max(p.y, q.y) - y1) < 0.05:
                tx.append(p.x)
            if abs(p.y - q.y) < 0.01 and 0.5 < abs(p.x - q.x) < 5 and abs(min(p.x, q.x) - x0) < 0.05:
                ty.append(p.y)
        tx, ty = np.unique(np.round(tx, 2)), np.unique(np.round(ty, 2))
        px, vx, py, vy = [], [], [], []
        for w in words:
            v = _num(w[4])
            if v is None:
                continue
            cx, cy = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
            if y1 < cy < y1 + 14 and x0 - 8 < cx < x1 + 8 and len(tx):
                i = np.argmin(np.abs(tx - cx))
                if abs(tx[i] - cx) < 4:
                    px.append(tx[i]); vx.append(v * xscale)
            if x0 - 30 < w[2] < x0 and y0 - 5 < cy < y1 + 5 and len(ty):
                i = np.argmin(np.abs(ty - cy))
                if abs(ty[i] - cy) < 4:
                    py.append(ty[i]); vy.append(v * yscale)
        if len(py) < 2:
            raise RuntimeError(f"frame {rect}: {len(px)} x labels, {len(py)} y labels")
        # stacked frames carry the x labels on the lowest one only: ax is
        # filled in by panels() from the frame below
        self.ax = np.polyfit(px, vx, 1) if len(px) > 1 else None
        self.ay = np.polyfit(py, vy, 1)
        self.misfit = (np.abs(np.polyval(self.ax, px) - vx).max() / np.ptp(vx)
                       if len(px) > 1 else 0.0,
                       np.abs(np.polyval(self.ay, py) - vy).max() / np.ptp(vy))
        self.xticks, self.yticks = sorted(zip(vx, px)), sorted(zip(vy, py))
        eps = 3.0          # markers on the frame stick out by their radius
        self.drawings = [d for d in drawings
                         if d["rect"].x0 >= x0 - eps and d["rect"].x1 <= x1 + eps
                         and d["rect"].y0 >= y0 - eps and d["rect"].y1 <= y1 + eps]

    def X(self, px):
        return np.polyval(self.ax, np.asarray(px, float))

    def Y(self, py):
        return np.polyval(self.ay, np.asarray(py, float))

    def _centres(self, sel):
        pts = [(0.5 * (d["rect"].x0 + d["rect"].x1), 0.5 * (d["rect"].y0 + d["rect"].y1))
               for d in self.drawings if sel(d)]
        pts = np.array(sorted(set((round(a, 2), round(b, 2)) for a, b in pts)))
        if not len(pts):
            return np.empty((0, 2))
        return np.c_[self.X(pts[:, 0]), self.Y(pts[:, 1])]

    def circles(self):
        """Centres of the circles (four Bezier arcs each)."""
        return self._centres(lambda d: len(d["items"]) == 4 and all(i[0] == "c" for i in d["items"])
                             and abs(d["rect"].width - d["rect"].height) < 0.3)

    def crosses(self):
        """Centres of the crosses (two diagonal strokes each)."""
        def diag(d):
            if len(d["items"]) != 1 or d["items"][0][0] != "l":
                return False
            p, q = d["items"][0][1], d["items"][0][2]
            return 1 < abs(p.x - q.x) < 6 and abs(abs(p.x - q.x) - abs(p.y - q.y)) < 0.3
        return self._centres(diag)

    def lines(self, min_points=3, min_width=0.0, max_width=1e9):
        """Polylines drawn in the frame: list of (line width, x, y)."""
        out = []
        for d in self.drawings:
            w = d.get("width") or 0.0
            if not (min_width <= w <= max_width) or d.get("dashes") not in ("[] 0", None):
                continue
            if len(d["items"]) == 1 and d["items"][0][0] == "l":
                # the frame and its ticks: axis-parallel, starting on an edge
                p, q = d["items"][0][1], d["items"][0][2]
                x0, y0, x1, y1 = self.rect
                on_edge = min(abs(p.x - x0), abs(p.x - x1), abs(p.y - y0), abs(p.y - y1)) < 0.05
                if (abs(p.x - q.x) < 0.01 or abs(p.y - q.y) < 0.01) and on_edge:
                    continue
            pts = []
            for it in d["items"]:
                if it[0] == "l":
                    if not pts:
                        pts.append((it[1].x, it[1].y))
                    pts.append((it[2].x, it[2].y))
                elif it[0] == "c":
                    if not pts:
                        pts.append((it[1].x, it[1].y))
                    pts.append((it[4].x, it[4].y))
            if len(pts) >= min_points:
                pts = np.array(pts)
                if np.ptp(pts[:, 0]) < 8 and np.ptp(pts[:, 1]) < 8 \
                        and np.hypot(*(pts[0] - pts[-1])) < 0.6 * max(np.ptp(pts[:, 0]), 0.1):
                    continue                       # a closed marker
                out.append((w, self.X(pts[:, 0]), self.Y(pts[:, 1])))
        return out

    def vticks(self, min_len=5.0):
        """x of the free vertical segments in the frame (gauge marks)."""
        xs = []
        x0, y0, x1, y1 = self.rect
        for d in self.drawings:
            if len(d["items"]) != 1 or d["items"][0][0] != "l":
                continue
            p, q = d["items"][0][1], d["items"][0][2]
            if abs(p.x - q.x) < 0.01 and min_len < abs(p.y - q.y) < 0.5 * (y1 - y0) \
                    and min(p.y, q.y) > y0 + 0.5 and max(p.y, q.y) < y1 - 0.5:
                xs.append(p.x)
        return self.X(np.unique(np.round(xs, 2)))


def panels(pdf, page_no, clip=None, min_size=40.0, xscale=1.0, yscale=1.0):
    """The plot frames of a page (1-based number), top to bottom and left to
    right, within clip = (x0, y0, x1, y1) in points.  yscale: factor of the
    y labels, one number or one for each frame (the "x 10^-3" that MATLAB
    prints above an axis)."""
    doc = pymupdf.open(pdf) if isinstance(pdf, str) else pdf
    page = doc[page_no - 1]
    drawings = page.get_drawings()
    words = page.get_text("words")
    rects = set()
    for d in drawings:
        for it in d["items"]:
            if it[0] == "re":
                r = it[1]
                if r.width > min_size and r.height > min_size:
                    if clip is None or (r.x0 >= clip[0] and r.y0 >= clip[1] and r.x1 <= clip[2] and r.y1 <= clip[3]):
                        rects.add((round(r.x0, 2), round(r.y0, 2), round(r.x1, 2), round(r.y1, 2)))
    out = []
    order = sorted(rects, key=lambda r: (round(r[1] / 20), r[0]))
    ys = list(yscale) if np.ndim(yscale) else [yscale] * len(order)
    for r, sc in zip(order, ys):
        try:
            out.append(Panel(page, r, drawings, words, xscale, sc))
        except RuntimeError as e:
            print("skipped:", e)
    for p in out:
        if p.ax is None:
            same = [q for q in out if q.ax is not None and abs(q.rect[0] - p.rect[0]) < 0.1
                    and abs(q.rect[2] - p.rect[2]) < 0.1 and q.rect[1] > p.rect[1]]
            if not same:
                raise RuntimeError(f"frame {p.rect}: no x labels and no frame below it")
            p.ax = min(same, key=lambda q: q.rect[1]).ax
    return out


def joined(lines, width, tol=0.05):
    """The polylines of one line width as one curve, ordered along x."""
    parts = [np.c_[x, y] for w, x, y in lines if abs(w - width) < tol]
    a = np.vstack(sorted(parts, key=lambda q: q[0, 0]))
    return a[:, 0], a[:, 1]


def match(x, *sets, tol=0.02):
    """Rows (x, y of set 1, y of set 2, ...) for the points of the sets that
    share an x (within tol); NaN where a set has none."""
    out = []
    for xi in x:
        row = [xi]
        for s in sets:
            i = np.argmin(np.abs(s[:, 0] - xi)) if len(s) else 0
            row.append(s[i, 1] if len(s) and abs(s[i, 0] - xi) < tol else np.nan)
        out.append(row)
    return np.array(out)
