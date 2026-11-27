#!/usr/bin/env python3
"""Monai valley on meshes refined around the valley (monai_refined.py)
against the uniform grids (monai_valley.py) and the measured run-up.

usage: plot_monai_refined.py <figure prefix> <mesh run> <run> [<run> ...]
    run = <path>:<label>; <mesh run>: the path of the refined run whose
    mesh is drawn
writes <prefix>.png/.pdf and prints the numbers of the README/report table.
"""
import os, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import monai_valley as mv
import monai_refined as mr
import plot_monai_valley as pmv

STYLES = [("-", "#1f5fa8", 1.3), (":", "#4a90d9", 1.6), ("--", "#c0392b", 1.3),
          ((0, (5, 1.5, 1, 1.5)), "#e67e22", 1.4), ("-.", "#7f8c8d", 1.1)]
VALLEY = (4.6, 5.4)


def unstructured(case):
    return "mesh" in mv.params(case)


def polygons(case):
    """Plan-view polygons of the columns, from the bottom faces of mesh.msh
    (the flat mesh before it is laid on the bed)."""
    t = open(os.path.join(case, "mesh.msh")).read()
    nd = t.split("$Nodes")[1].split("$EndNodes")[0].strip().split("\n")[1:]
    N = {}
    for l in nd:
        v = l.split()
        N[int(v[0])] = (float(v[1]), float(v[2]), float(v[3]))
    out = []
    for l in t.split("$Elements")[1].split("$EndElements")[0].strip().split("\n")[1:]:
        v = l.split()
        if int(v[1]) in (2, 3):
            p = [N[int(i)] for i in v[3 + int(v[2]):]]
            if all(abs(q[2] + mv.D) < 1e-9 for q in p):
                out.append([(q[0], q[1]) for q in p])
    return out


def flooded(case):
    """(x, y, zb, flooded) of the columns of a run of either kind."""
    if unstructured(case):
        x, y, z, e, fl = mr.emax(case)
        return x, y, z, fl
    pr = mv.params(case)
    e, fl, z = mv.emax(case)
    xc, yc = pmv.cell_xy(pr["dx"])
    X, Y = np.meshgrid(xc, yc)
    return X.ravel(), Y.ravel(), z.ravel(), fl.ravel()


def runup_profile(case):
    """Highest flooded ground of the valley (VALLEY in x) in bins of y: the
    rows of a uniform grid; bins of 1.4 cm over the refined box otherwise."""
    x, y, z, fl = flooded(case)
    if unstructured(case):
        bins = np.arange(mr.BOX[2], mr.BOX[3] + 1e-9, mv.DX0)
    else:
        dy = mv.grid(mv.params(case)["dx"])[2]
        bins = np.arange(0.0, mv.LY + 1e-9, dy)
        bins = bins[(bins > 1.43) & (bins < 2.47)]
    s = fl & (z > 0.0) & (x > VALLEY[0]) & (x < VALLEY[1])
    k = np.digitize(y[s], bins) - 1
    R = np.full(len(bins) - 1, np.nan)
    for i in range(len(R)):
        m = k == i
        if m.any():
            R[i] = z[s][m].max()
    ok = np.isfinite(R)
    return (0.5 * (bins[1:] + bins[:-1]))[ok], R[ok]


def numbers(case):
    """Run-up at the three lab points, highest in the valley, columns,
    time steps and run time."""
    if unstructured(case):
        r = [a["R"] for a in mr.runup(case)]
        ncol = int(mv.params(case)["columns"])
    else:
        r = [a["R"] for a in mv.runup(case)]
        nx, ny, _ = mv.grid(mv.params(case)["dx"])
        ncol = nx * ny
    x, y, z, fl = flooded(case)
    s = fl & (z > 0.0) & (y > 1.75) & (y < 2.05) & (x > VALLEY[0]) & (x < VALLEY[1])
    steps, clock = 0, np.nan
    for l in open(os.path.join(case, "log.NHFOAM"), errors="replace"):
        if l.startswith("Time ="):
            steps += 1
        elif l.startswith("ExecutionTime"):
            clock = float(l.split()[2])
    return r, z[s].max(), ncol, steps, clock


def main():
    prefix, meshrun = sys.argv[1], sys.argv[2]
    runs = [tuple(a.split(":", 1)) for a in sys.argv[3:]]
    obs = mv.read_runup()
    tl, gl = mv.read_lab()

    fig = plt.figure(figsize=(12.5, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15])
    axm = fig.add_subplot(gs[0, :])
    axz = fig.add_subplot(gs[1, 0])
    axr = fig.add_subplot(gs[1, 1])

    P = polygons(meshrun)
    xn, yn, Zn = mv.nodes(mv.DX0)
    axm.add_collection(PolyCollection(P, facecolors="none", edgecolors="0.35",
                                      linewidths=0.15, rasterized=True))
    axm.contour(xn, yn, 100 * Zn, levels=[0.0], colors="#c0392b", linewidths=0.9)
    b = mr.BOX
    axm.plot([b[0], b[1], b[1], b[0], b[0]], [b[2], b[2], b[3], b[3], b[2]],
             color="#1f5fa8", lw=0.9)
    for g, (x, y) in mv.GAUGES.items():
        if x is not None:
            axm.plot(x, y, "o", color="#1f5fa8", ms=3)
    axm.set_xlim(0, mv.LX)
    axm.set_ylim(0, mv.LY)
    axm.set_aspect("equal")
    axm.set_xlabel(r"$x$ [m]")
    axm.set_ylabel(r"$y$ [m]")
    axm.set_title(rf"(a) mesh of {len(P)} columns (still-water shoreline red, "
                  r"refined box and gauges 5, 7, 9 blue)", loc="left")

    pmv.valley_map(axz, box=(4.7, 5.3, 1.6, 2.45))
    axz.add_collection(PolyCollection(P, facecolors="none", edgecolors="0.3",
                                      linewidths=0.12, rasterized=True))
    print(f"{'run':30s} {'head':>6s} {'N1':>6s} {'N2':>6s} {'valley':>7s} {'columns':>8s} "
          f"{'steps':>6s} {'time [min]':>10s}   crest g5/g7/g9 [%], bore [s], rms [cm]")
    for k, (path, label) in enumerate(runs):
        st = STYLES[min(k, len(STYLES) - 1)]
        x, y, z, fl = flooded(path)
        s = (x > 4.6) & (x < 5.4) & (y > 1.5) & (y < 2.55)
        axz.tricontour(x[s], y[s], fl[s].astype(float), levels=[0.5], colors=[st[1]],
                       linestyles=[st[0]], linewidths=st[2])
        yb, R = runup_profile(path)
        axr.plot(yb, 100 * R, ls=st[0], color=st[1], lw=st[2], label=label)
        r, rv, ncol, steps, clock = numbers(path)
        tg, vg = mv.model_gauges(path)
        txt = []
        for g in ("ch5", "ch7", "ch9"):
            HL, _ = pmv.crest(tl, gl[g], *pmv.WIN)
            Hm, _ = pmv.crest(tg[g], vg[g], *pmv.WIN)
            wl = (tl >= pmv.WIN[0]) & (tl <= pmv.WIN[1])
            rms = np.sqrt(np.mean((np.interp(tl[wl], tg[g], vg[g]) - gl[g][wl]) ** 2))
            txt.append(f"{100 * (Hm - HL) / HL:+.0f}% {mv.arrival(tg[g], vg[g]) - mv.arrival(tl, gl[g]):+.2f} "
                       f"{100 * rms:.2f}")
        print(f"{label:30s} {100 * r[0]:6.2f} {100 * r[1]:6.2f} {100 * r[2]:6.2f} "
              f"{100 * rv:7.2f} {ncol:8d} {steps:6d} {clock / 60:10.1f}   " + " | ".join(txt))
    for k, (x, y, o) in enumerate(obs):
        axz.plot(x, y, "ko", ms=4)
        axr.errorbar(y, 100 * o.mean(), yerr=[[100 * (o.mean() - o.min())],
                                               [100 * (o.max() - o.mean())]],
                     fmt="ko", ms=4, capsize=3, label="lab (6 repeats)" if k == 0 else None)
    axz.set_xlabel(r"$x$ [m]")
    axz.set_ylabel(r"$y$ [m]")
    axz.set_title("(b) the mesh in the valley and the highest inundation", loc="left")
    axr.set_xlabel(r"$y$ [m]")
    axr.set_ylabel("run-up [cm]")
    axr.set_xlim(1.45, 2.45)
    axr.grid(alpha=0.3)
    axr.set_title(r"(c) highest flooded ground, $4.6 < x < 5.4$ m", loc="left")
    axr.legend(loc="upper right", fontsize=9, frameon=False)
    fig.tight_layout()
    save(fig, prefix + ".png")


if __name__ == "__main__":
    main()
