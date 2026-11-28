#!/usr/bin/env python3
"""The laboratory data that Zijlema, Stelling & Smit (2011, Coastal
Engineering 58, 992-1012) compare SWASH with, read from the figures of the
paper.  The figures are vector graphics, so the markers and lines come out
at the coordinates they were drawn at (the axes are fitted to the tick
labels, misfit below 0.1% of the axis); the values are as good as the
plotting of the paper, not as the original measurements.

    digitize_zijlema2011.py <paper.pdf> <outdir>

writes (x in m along the flume, as in the paper)

  bj78_profile.txt        bed of Battjes & Janssen's (1978) bar-trough beach,
                          z relative to the still water of Fig. 14
  bj78_gauges.txt         x of the gauges of Fig. 14
  bj78_run13.txt, bj78_run15.txt
                          x, Hm0 measured, Hm0 of SWASH (Fig. 15)
  boers_profile.txt       bed of Boers's (1996) flume (Fig. 17)
  boers_gauges.txt        x of the eight stations of Fig. 17
  boers1b_hm0.txt, boers1b_tm02.txt          measured (Fig. 19, circles)
  boers1b_hm0_swash.txt, boers1b_tm02_swash.txt   SWASH (Fig. 19, lines)
  boers1b_spec_<n>.txt, boers1b_spec_<n>_swash.txt
                          f [Hz], E [m2/Hz] at station n: observed, SWASH
                          (Fig. 18)
  haller_B_H_<y>.txt, haller_B_setup_<y>.txt
                          x, measured wave height / set-up [cm] on the
                          sections y = 11.23 and 13.68 m (Figs. 24-25), and
                          ..._swash.txt the lines of SWASH
  haller_B_U_<x>.txt, haller_B_V_<x>.txt
                          y, measured mean cross-shore / longshore current
                          [cm/s] on five sections (Figs. 26-27), and
                          ..._swash.txt
  haller_B_bathy.txt      the bed of the basin as SWASH used it (Fig. 23):
                          z [m] on 88 x 122 points, 0.15 m apart, with x and
                          y in the header (see haller_bathymetry)

The paper is not part of this repository and neither are these files.
"""
import os
import sys
import numpy as np
import figdig


def save(out, name, a, header):
    a = np.atleast_2d(np.asarray(a, float))
    np.savetxt(os.path.join(out, name), a[np.argsort(a[:, 0], kind="stable")], fmt="%.5g",
               header=header + "\nread from the figures of Zijlema et al. (2011), Coast. Eng. 58")
    print(f"{name}: {len(a)} points")


def longest(lines):
    return max(lines, key=lambda l: len(l[1]))


def haller_bathymetry(pdf, page_no=18, ymax=320.0):
    """The bed of Fig. 23, a MATLAB mesh plot.  Every face of the mesh is in
    the file as two triangles, so all vertices are there, also the hidden
    ones.  The view is a parallel projection: the axes box (x 2..18 m, y
    0..20 m, z -0.4..0.1 m, from the tick labels) gives the map from (x, y,
    z) to the page; the mesh edges along x and along y differ in their
    horizontal step on the page, which orders the vertices on a lattice
    (i, j); the horizontal position then gives the spacing and one
    combination of the two offsets, and the vertical position gives z.
    The other offset is taken so that y starts at 0, which puts the
    seaward edge at x = 3.95 m, where the paper has its boundary.
    Returns x, y, z[i, j]."""
    import collections
    import pymupdf
    page = pymupdf.open(pdf)[page_no - 1]
    dr = [d for d in page.get_drawings() if d["rect"].y1 < ymax]
    box = [d for d in dr if d["type"] == "f" and len(d["items"]) == 3]
    floor = [box[0]["items"][0][1]] + [it[2] for it in box[0]["items"]]
    wall = [box[1]["items"][0][1]] + [it[2] for it in box[1]["items"]]
    o, py, px = floor[0], floor[1], floor[3]          # (2, 0), (2, 20), (18, 0) at z = -0.4
    assert abs(wall[0].x - py.x) < 0.01 and abs(wall[1].x - py.x) < 0.01
    ax, ay = (px.x - o.x) / 16.0, (py.x - o.x) / 20.0
    bx, by = (px.y - o.y) / 16.0, (py.y - o.y) / 20.0
    bz = (wall[1].y - wall[0].y) / 0.5
    tris = []
    for d in dr:
        its = d["items"]
        if len(its) == 2 and its[0][0] == "l" and its[1][0] == "l" \
                and d.get("dashes") in ("[] 0", None):
            tris.append([(round(q.x, 2), round(q.y, 2)) for q in (its[0][1], its[0][2], its[1][2])])
    V = sorted(set(v for t in tris for v in t))
    idx = {v: i for i, v in enumerate(V)}
    P = np.array(V)
    sx_x, sx_y = abs(ax) * 0.15, abs(ay) * 0.15       # steps of a 0.15 m mesh, to tell edges apart
    adj = collections.defaultdict(list)
    for t in tris:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            ia, ib = idx[a], idx[b]
            d = P[ib, 0] - P[ia, 0]
            if abs(abs(d) - sx_x) < 0.12:
                st = (1, 0) if d * ax > 0 else (-1, 0)
            elif abs(abs(d) - sx_y) < 0.12:
                st = (0, 1) if d * ay > 0 else (0, -1)
            else:
                continue                               # the diagonal of a face
            adj[ia].append((ib, st)); adj[ib].append((ia, (-st[0], -st[1])))
    best, seen = {}, set()
    for start in adj:
        if start in seen:
            continue
        ij = {start: (0, 0)}
        stack = [start]
        while stack:
            a = stack.pop()
            for b, st in adj[a]:
                c = (ij[a][0] + st[0], ij[a][1] + st[1])
                if b not in ij:
                    ij[b] = c; stack.append(b)
                else:
                    assert ij[b] == c, "the mesh does not order on a lattice"
        seen.update(ij)
        if len(ij) > len(best):
            best = ij
    K = np.array(list(best))
    I = np.array([best[k] for k in K])
    I -= I.min(axis=0)
    c = np.linalg.lstsq(np.c_[I, np.ones(len(I))], P[K, 0], rcond=None)[0]
    dx, dy = c[0] / ax, c[1] / ay
    y = dy * np.arange(I[:, 1].max() + 1)
    x0 = 2.0 + (c[2] - o.x) / ax                       # with y = 0 at j = 0
    x = x0 + dx * np.arange(I[:, 0].max() + 1)
    zv = (P[K, 1] - o.y - bx * (x[I[:, 0]] - 2.0) - by * y[I[:, 1]]) / bz - 0.4
    S = np.zeros((len(x), len(y))); N = np.zeros_like(S); Q = np.zeros_like(S)
    np.add.at(S, (I[:, 0], I[:, 1]), zv); np.add.at(N, (I[:, 0], I[:, 1]), 1)
    np.add.at(Q, (I[:, 0], I[:, 1]), zv * zv)
    assert N.min() > 0, "vertices missing"
    Z = S / N
    spread = np.sqrt(np.maximum(Q / N - Z * Z, 0)).max()
    print(f"haller bathymetry: {len(x)} x {len(y)} points, spacing {dx:.4f} x {dy:.4f} m, "
          f"x from {x[0]:.3f} m, copies of a vertex agree within {1000 * spread:.2f} mm")
    return x, y, Z


def main(pdf, out):
    os.makedirs(out, exist_ok=True)
    worst = 0.0

    # Battjes & Janssen: Fig. 14 (set-up), Fig. 15 (Hm0, runs 13 and 15)
    ps = figdig.panels(pdf, 14, clip=(0, 40, 600, 310))
    worst = max([worst] + [max(p.misfit) for p in ps])
    w, x, z = longest(ps[0].lines(min_width=1.0))
    save(out, "bj78_profile.txt", np.c_[x, z], "x [m], z [m] (still water of Fig. 14 at z = 0)")
    save(out, "bj78_gauges.txt", np.c_[ps[0].vticks(), np.arange(1, 14)], "x [m], gauge")
    for p, run in zip(ps[1:], (13, 15)):
        lab, sw = p.crosses(), p.circles()
        save(out, f"bj78_run{run}.txt", figdig.match(lab[:, 0], lab, sw),
             "x [m], Hm0 measured [m], Hm0 SWASH [m]")

    # Boers 1B: Fig. 17 (bed), Fig. 18 (spectra), Fig. 19 (Hm0, Tm02)
    ps = figdig.panels(pdf, 15, clip=(0, 40, 600, 200))
    worst = max([worst] + [max(p.misfit) for p in ps])
    w, x, z = longest(ps[0].lines(min_width=1.0))
    save(out, "boers_profile.txt", np.c_[x, z], "x [m], z [m]")
    save(out, "boers_gauges.txt", np.c_[ps[0].vticks(), np.arange(1, 9)], "x [m], station")
    ps = figdig.panels(pdf, 15, clip=(0, 300, 600, 720),
                       yscale=[1, 1, 1, 1, 1e-3, 1e-3, 1e-3, 1e-3])
    worst = max([worst] + [max(p.misfit) for p in ps])
    for n, p in enumerate(ps, 1):
        ls = p.lines(min_points=2)
        widths = sorted(set(round(l[0], 2) for l in ls))
        assert len(widths) == 2, widths
        save(out, f"boers1b_spec_{n}.txt", np.c_[figdig.joined(ls, widths[1])],
             "f [Hz], E observed [m2/Hz]")
        save(out, f"boers1b_spec_{n}_swash.txt", np.c_[figdig.joined(ls, widths[0])],
             "f [Hz], E SWASH [m2/Hz]")
    ps = figdig.panels(pdf, 16, clip=(0, 40, 600, 160))
    worst = max([worst] + [max(p.misfit) for p in ps])
    for p, q, unit in zip(ps, ("hm0", "tm02"), ("Hm0 [m]", "Tm02 [s]")):
        save(out, f"boers1b_{q}.txt", p.circles(), f"x [m], {unit} measured")
        ls = p.lines(min_points=2)
        save(out, f"boers1b_{q}_swash.txt", np.c_[figdig.joined(ls, longest(ls)[0])],
             f"x [m], {unit} SWASH")

    # Haller et al. (2002), test B: Figs. 24-27
    def sections(page, clip, q, names, unit, coord):
        nonlocal worst
        ps = figdig.panels(pdf, page, clip=clip, min_size=25.0)
        worst = max([worst] + [max(p.misfit) for p in ps])
        assert len(ps) == len(names), (page, q, len(ps))
        for p, n in zip(ps, names):
            save(out, f"haller_B_{q}_{n}.txt", p.circles(), f"{coord} [m], {unit} measured")
            ls = p.lines(min_points=2)
            save(out, f"haller_B_{q}_{n}_swash.txt", np.c_[figdig.joined(ls, longest(ls)[0])],
                 f"{coord} [m], {unit} SWASH")
    ys = ("11.23", "13.68")
    xs = ("14", "13", "12.2", "11.2", "10")
    sections(18, (0, 440, 600, 720), "H", ys, "H [cm]", "x")
    sections(19, (0, 40, 600, 310), "setup", ys, "set-up [cm]", "x")
    sections(19, (0, 440, 600, 720), "U", xs, "U [cm/s]", "y")
    sections(20, (0, 40, 600, 400), "V", xs, "V [cm/s]", "y")
    x, y, Z = haller_bathymetry(pdf)
    np.savetxt(os.path.join(out, "haller_B_bathy.txt"), Z, fmt="%.4f",
               header=f"z [m] of the bed, rows x = {x[0]:.4f} + {x[1] - x[0]:.5f} i (i = 0..{len(x) - 1}), "
                      f"columns y = {y[1] - y[0]:.5f} j (j = 0..{len(y) - 1})\n"
                      "read from Fig. 23 of Zijlema et al. (2011), Coast. Eng. 58")
    print(f"largest misfit of an axis: {100 * worst:.3f}% of its range")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
