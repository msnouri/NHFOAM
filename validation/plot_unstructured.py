#!/usr/bin/env python3
"""Unstructured meshes (unstructured.py): figures and numbers of the report.

usage: plot_unstructured.py runup <out.png> <structuredCase> <labRunup.txt> <labProfDir> <case>:<label> ...
       plot_unstructured.py cylinder <out.png> <case>
(also writes out.pdf)
"""
import os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from latex_style import use_latex, save
import wetdry
from longshore import cell_centres
import unstructured as un

use_latex()
G = 9.81


def read_scalar(case, time, name):
    txt = open(os.path.join(case, time, name)).read()
    i = txt.index("internalField")
    m = re.search(r"nonuniform\s+List<scalar>\s*(\d+)\s*\(", txt[i:])
    if m is None:
        v = float(re.search(r"uniform\s+([-0-9.eE+]+)", txt[i:]).group(1))
        return None, v
    n = int(m.group(1))
    a = i + m.end()
    return np.array(txt[a:txt.index(")", a)].split(), float)[:n], None


def times(case):
    out = []
    for t in os.listdir(case):
        try:
            out.append((float(t), t))
        except ValueError:
            pass
    return sorted(out)


def shoreline(case):
    s = np.loadtxt(os.path.join(case, "postProcessing/wettingDrying/0/shoreline.dat"),
                   comments="#")
    return s[:, 0], s[:, 2]


def runup(argv):
    out, ref, labrunup, labprof = argv[:4]
    runs = [(a.split(":")[0], a.split(":")[1]) for a in argv[4:]]
    d = wetdry.RU_D
    H = 0.0588
    ts = np.sqrt(d / G)
    xc0, xshore = wetdry.runup_geometry(H, d)
    film = 0.005
    TST = [15, 20, 25, 30]
    allruns = [(ref, r"structured, $\Delta x = 2$ cm")] + runs
    cols = ["k", "#1f5fa8", "#e67e22", "#2e8b57"]
    lss = ["-", "--", "-.", ":"]
    fig = plt.figure(figsize=(7.5, 9.0))
    gs = fig.add_gridspec(len(TST) + 1, 1, height_ratios=[1] * len(TST) + [1.3])
    axes = [fig.add_subplot(gs[i]) for i in range(len(TST))]
    axs = fig.add_subplot(gs[-1])
    xb = np.linspace(-6, 30, 400)
    prof = {}
    for case, label in allruns:
        C = cell_centres(case, "0")
        x = C[:, 0]
        h = wetdry.runup_depth(x, d)
        p = {}
        for tv, tn in times(case):
            tsr = round(tv / ts)
            if tsr not in TST:
                continue
            e, u = read_scalar(case, tn, "eta")
            if e is None:
                continue
            D = e + h
            # cross-strip mean per 2 cm bin of x
            bins = np.arange(0, wetdry.RU_LEN + 0.02, 0.02)
            ib = np.digitize(x, bins)
            xm = np.array([x[ib == i].mean() for i in np.unique(ib)])
            em = np.array([e[ib == i].mean() for i in np.unique(ib)])
            Dm = np.array([D[ib == i].mean() for i in np.unique(ib)])
            em = np.where(Dm > 2 * film, em, np.nan)
            p[tsr] = (xm, em)
        prof[label] = p
    for ax, tsr in zip(axes, TST):
        xx = xshore - xb * d
        ax.fill_between(xb, -0.2 * np.ones_like(xb), -wetdry.runup_depth(xx, d) / d,
                        color="0.85", lw=0)
        f = os.path.join(labprof, f"Case0_3.t={tsr}")
        if os.path.exists(f):
            L = np.loadtxt(f)
            ax.plot(L[:, 0], L[:, 1], "o", ms=3, mfc="none", color="#c0392b",
                    label="Synolakis (1987)")
        for (case, label), col, ls in zip(allruns, cols, lss):
            if tsr in prof[label]:
                xm, em = prof[label][tsr]
                ax.plot((xshore - xm) / d, em / d, ls, color=col, lw=1.3, label=label)
        ax.set_ylim(-0.1, 0.65)
        ax.set_xlim(-6, 30)
        ax.invert_xaxis()
        ax.set_ylabel(r"$\eta/d$")
        ax.text(0.98, 0.82, rf"$t\sqrt{{g/d}} = {tsr}$", transform=ax.transAxes, ha="right")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel(r"$x/d$ (offshore of the initial shoreline)")
    lab = np.loadtxt(labrunup, comments="#")
    near = np.abs(lab[:, 0] - H / d) < 0.02
    axs.axhspan(lab[near, 1].min(), lab[near, 1].max(), color="#c0392b", alpha=0.15, lw=0,
                label=f"lab maximum, $H/d = {H/d - 0.02:.2f}$--${H/d + 0.02:.2f}$")
    print(f"lab maximum run-up R/d = {lab[near, 1].min():.3f}-{lab[near, 1].max():.3f}")
    for (case, label), col, ls in zip(allruns, cols, lss):
        t, R = shoreline(case)
        axs.plot(t / ts, R / d, ls, color=col, lw=1.3, label=label)
        nc = len(cell_centres(case, "0"))
        logf = os.path.join(case, "log.NHFOAM")
        cpu = re.findall(r"ExecutionTime = ([0-9.]+)", open(logf).read())
        print(f"{label}: {nc} cells, maximum run-up R/d = {R.max()/d:.3f} at t* = "
              f"{t[np.argmax(R)]/ts:.1f}, CPU {float(cpu[-1]):.0f} s")
        for tsr in TST:
            if tsr in prof[label] and label != allruns[0][1]:
                xm, em = prof[label][tsr]
                xr, er = prof[allruns[0][1]][tsr]
                ei = np.interp(xr, xm, np.nan_to_num(em, nan=-1))
                ok = np.isfinite(er) & (ei > -0.5)
                print(f"    t* = {tsr}: rms difference to the structured run "
                      f"{100*np.sqrt(np.mean((ei[ok]-er[ok])**2))/d:.2f}% of d")
    axs.set_xlabel(r"$t\sqrt{g/d}$")
    axs.set_ylabel(r"shoreline $z/d$")
    axs.set_xlim(0, 60)
    axs.grid(alpha=0.3)
    h_, l_ = [], []
    for ax in axes + [axs]:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l_:
                h_.append(hh); l_.append(ll)
    fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save(fig, out)


def cylinder(argv):
    out, case = argv[:2]
    p = dict(un.CYL)
    for l in open(os.path.join(case, "params.txt")):
        k, v = l.split()
        try:
            p[k] = float(v)
        except ValueError:
            pass
    h, T, H, a = p["h"], p["T"], p["H"], p["a"]
    w = 2 * np.pi / T
    k = w * w / G
    for _ in range(60):
        k -= (G * k * np.tanh(k * h) - w * w) / (G * np.tanh(k * h) + G * k * h / np.cosh(k * h) ** 2)
    fn = os.path.join(case, "postProcessing/gauges/0/eta.dat")
    d = np.loadtxt(fn, comments="#")
    t = d[:, 0]
    names, xy = [], []
    for l in open(fn):
        if l.startswith("# "):
            m = re.match(r"# (\S+) \(([-0-9.eE+]+) ([-0-9.eE+]+)", l)
            if m:
                names.append(m.group(1)); xy.append((float(m.group(2)), float(m.group(3))))
    xy = np.array(xy)
    # analysis windows of 4 periods: the last one before the side-wall and
    # paddle reflections return, and the end of the run
    res = {}
    for t2 in (12.0, t[-1]):
        m = (t > t2 - 4 * T) & (t <= t2)
        A = 2 * np.mean(d[m, 1:] * np.exp(1j * w * t[m])[:, None], axis=0)
        res[t2] = np.abs(A) / (H / 2)
    mf = np.abs(un.maccamy_fuchs(xy[:, 0], xy[:, 1], k, a))
    kinds = np.array([n[0] for n in names])
    fig = plt.figure(figsize=(7.5, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    # mesh near the cylinder (top surface triangles)
    lines = open(os.path.join(case, "mesh.msh")).read().split("\n")
    i = lines.index("$Nodes"); n = int(lines[i + 1])
    nodes = {int(l.split()[0]): list(map(float, l.split()[1:4])) for l in lines[i + 2:i + 2 + n]}
    i = lines.index("$Elements"); n = int(lines[i + 1])
    polys = []
    for l in lines[i + 2:i + 2 + n]:
        s = list(map(int, l.split()))
        if s[1] == 2:
            v = s[3 + s[2]:]
            P = [nodes[q] for q in v]
            if abs(P[0][2]) < 1e-6 and max(abs(q[0]) for q in P) < 1.2 and max(q[1] for q in P) < 1.2:
                polys.append([q[:2] for q in P])
    ax0.add_collection(PolyCollection(polys, facecolors="none", edgecolors="0.35", linewidths=0.3))
    ax0.add_patch(plt.Circle((0, 0), a, color="0.6"))
    c = kinds == "c"
    ax0.plot(xy[c, 0], xy[c, 1], "o", ms=3, color="#c0392b")
    ax0.set_xlim(-1.2, 1.2); ax0.set_ylim(0, 1.2)
    ax0.set_aspect("equal")
    ax0.set_xlabel(r"$x$ [m]"); ax0.set_ylabel(r"$y$ [m]")
    ax0.set_title(r"(a) mesh near the cylinder", fontsize=10)
    ax1 = fig.add_subplot(gs[0, 1])
    th = np.degrees(np.arctan2(xy[c, 1], xy[c, 0]))
    tt = np.linspace(0, np.pi, 181)
    rr = a + 0.02
    ax1.plot(np.degrees(tt), np.abs(un.maccamy_fuchs(rr * np.cos(tt), rr * np.sin(tt), k, a)),
             "k-", lw=1.0, label="MacCamy \\& Fuchs")
    ax1.plot(th, res[12.0][c], "o", ms=4, color="#1f5fa8", label=r"NHFOAM, $t = 8$--12 s")
    ax1.plot(th, res[t[-1]][c], "x", ms=4, color="#e67e22", label=rf"NHFOAM, $t = {t[-1]-4*T:g}$--{t[-1]:g} s")
    ax1.set_xlabel(r"$\theta$ [deg] (0: lee side)")
    ax1.set_ylabel(r"$|\eta|/a_I$ at $r = a + 2$ cm")
    ax1.set_title(r"(b) around the cylinder", fontsize=10)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)
    ax2 = fig.add_subplot(gs[1, :])
    for kind, lab, col in (("f", "axis, up-wave", "#1f5fa8"), ("b", "axis, lee", "#2e8b57"),
                           ("s", r"across ($x = 0$)", "#8e44ad")):
        s = kinds == kind
        r = np.hypot(xy[s, 0], xy[s, 1])
        o = np.argsort(r)
        ax2.plot(r[o] / a, mf[s][o], "-", color=col, lw=1.0)
        ax2.plot(r[o] / a, res[12.0][s][o], "o", ms=3, color=col, label=lab)
    ax2.set_xlabel(r"$r/a$")
    ax2.set_ylabel(r"$|\eta|/a_I$")
    ax2.set_title(r"(c) along lines from the cylinder: NHFOAM ($t = 8$--12 s, symbols)"
                  r" and MacCamy \& Fuchs (lines)", fontsize=10)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, ncol=3)
    fig.tight_layout()
    save(fig, out)
    for t2, r in res.items():
        err = r - mf
        print(f"window ending at {t2:g} s: around the cylinder rms error {np.sqrt(np.mean(err[c]**2)):.3f},"
              f" max |error| {np.abs(err[c]).max():.3f} (of a_I); all gauges rms {np.sqrt(np.mean(err**2)):.3f}")
    print(f"kh = {k*h:.3f}, ka = {k*a:.3f}; largest |eta|/a_I: model {res[12.0][c].max():.3f},"
          f" theory {mf[c].max():.3f} (theta = 180)")


if __name__ == "__main__":
    {"runup": runup, "cylinder": cylinder}[sys.argv[1]](sys.argv[2:])
