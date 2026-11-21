#!/usr/bin/env python3
"""Solitary wave of H/h = 0.5 in constant depth (Ma et al. 2012, section 4.2):
surface elevation and velocities at mid-depth (sigma = 0.5) at t = 0, 10, 20
and 40 s, NHFOAM and NHWAVE against the exact solution.

usage: plot_ma2012_solitary.py [--old <nhfoamCase>] <nhfoamCase> <nhwaveCase> <out.png> [nLayers]
    --old: a second NHFOAM run (e.g. before a fix), drawn dotted
(also writes out.pdf)
Top: elevation over the flume; bottom: eta, u and w at mid-depth around the
crest at t = 40 s.
"""
import os, re, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import ma2012
from solitary_exact import SolitaryWave

args = sys.argv[1:]
old_case = None
if "--old" in args:
    i = args.index("--old")
    old_case = args[i + 1]
    del args[i:i + 2]
nf_case, nw_case, out = args[0:3]
nL = int(args[3]) if len(args) > 3 else 3
h = ma2012.SOL_h
sw = SolitaryWave(ma2012.SOL_H, d=h)
Nx = int(round(ma2012.SOL_LEN / ma2012.SOL_DX))
x = (np.arange(Nx) + 0.5) * ma2012.SOL_DX


def foam(case, t, name):
    s = open(os.path.join(case, f"{t:g}", name)).read().split("internalField")[1].split(";")[0]
    if name == "U":
        v = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", s), float)
        return v[:, 0].reshape(nL, Nx), v[:, 2].reshape(nL, Nx)
    return np.array(re.findall(r"[-0-9.eE+]+", s.split("(", 1)[1]), float).reshape(nL, Nx)[0]


def nhwave(case, i, name):
    return np.loadtxt(os.path.join(case, "results", f"{name}_{i:05d}"))


tw = np.loadtxt(os.path.join(nw_case, "results", "time"))
times = [0.0, 10.0, 20.0, 40.0]
mid = nL // 2

fig = plt.figure(figsize=(7.5, 6.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.15])
ax_all = fig.add_subplot(gs[0, :])
zoom = [fig.add_subplot(gs[1, j]) for j in range(3)]
c_nf, c_nw, c_old = "#1f5fa8", "#e67e22", "#7f8c8d"
print("t [s]  model   eta_max   u_max   w_max   crest x   (rel. err. of peaks vs exact)")
for t in times:
    # NHFOAM
    e_nf = foam(nf_case, t, "eta") if t > 0 else sw.eta(x - ma2012.SOL_X0)
    if t > 0:
        u_nf, w_nf = foam(nf_case, t, "U")
        u_nf, w_nf = u_nf[mid], w_nf[mid]
    # NHWAVE (outputs slightly after the nominal times)
    if t > 0:
        i = int(np.argmin(np.abs(tw - t))) + 1
        t_nw = tw[i - 1]
        e_nw = nhwave(nw_case, i, "eta")
        u_nw = nhwave(nw_case, i, "u")[mid]
        w_nw = nhwave(nw_case, i, "w")[mid]
    # exact
    for tt, tag in ([(t, "nf")] + ([(t_nw, "nw")] if t > 0 else [])):
        xs = x - ma2012.SOL_X0 - sw.c * tt
        ee = sw.eta(xs)
        zm = -h + (np.arange(nL)[mid] + 0.5) / nL * (h + ee)
        ue, we = sw.velocity(xs, zm)
        if tag == "nf":
            e_ex, u_ex, w_ex = ee, ue, we
        if t > 0:
            if tag == "nf":
                e, u, w = e_nf, u_nf, w_nf
            else:
                e, u, w = e_nw, u_nw, w_nw
            xc = x[np.argmax(e)]
            print(f"{tt:6.3f} {('NHFOAM' if tag == 'nf' else 'NHWAVE'):7s} "
                  f"{e.max():.4f}  {u.max():.4f}  {w.max():.4f}  {xc:7.2f}   "
                  f"({100*(e.max()/ee.max()-1):+.1f}%, {100*(u.max()/ue.max()-1):+.1f}%, "
                  f"{100*(w.max()/we.max()-1):+.1f}%; crest {xc - (ma2012.SOL_X0 + sw.c*tt):+.2f} m)")
    ax_all.plot(x, e_ex, "o", ms=2.4, mfc="white", mew=0.6, color="0.25",
                label="exact" if t == 0 else None, markevery=2)
    if t > 0:
        ax_all.plot(x, e_nw, "--", color=c_nw, lw=1.3, label="NHWAVE" if t == 10 else None)
        ax_all.plot(x, e_nf, "-", color=c_nf, lw=1.3, label="NHFOAM" if t == 10 else None)
    if t == times[-1]:
        zex = (e_ex, u_ex, w_ex)
        znw = (e_nw, u_nw, w_nw)
        znf = (e_nf, u_nf, w_nf)
        xc_ex = ma2012.SOL_X0 + sw.c * t
if old_case is not None:
    e_old = foam(old_case, times[-1], "eta")
    u_old, w_old = foam(old_case, times[-1], "U")
    zold = (e_old, u_old[mid], w_old[mid])
    xo = x[np.argmax(e_old)]
    print(f"old run at t = {times[-1]:g} s: eta_max {e_old.max():.4f}, crest x {xo:.2f}")
m = np.abs(x - xc_ex) < 9
labels = [r"$\eta$ [m]", r"$u$ at mid-depth [m/s]", r"$w$ at mid-depth [m/s]"]
for j, ax in enumerate(zoom):
    ax.plot(x[m], zex[j][m], "o", ms=2.6, mfc="white", mew=0.7, color="0.25", label="exact")
    if old_case is not None:
        ax.plot(x[m], zold[j][m], ":", color=c_old, lw=1.6, label="NHFOAM, first order in time")
    ax.plot(x[m], znw[j][m], "--", color=c_nw, lw=1.4, label="NHWAVE")
    ax.plot(x[m], znf[j][m], "-", color=c_nf, lw=1.4, label="NHFOAM")
    ax.set_title(labels[j], fontsize=11)
    ax.set_xlabel(r"$x$ [m]")
    ax.grid(alpha=0.3)
ax_all.set_ylabel(r"$\eta$ [m]")
ax_all.set_xlabel(r"$x$ [m]")
ax_all.set_xlim(0, 200)
ax_all.grid(alpha=0.3)
ax_all.set_title(r"$t$ = 0, 10, 20, 40 s", fontsize=11)
h, l = zoom[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=10, frameon=False)
fig.suptitle(rf"Solitary wave, $H/h = 0.5$, {nL} layers, $\Delta x = 0.1$ m; bottom: $t = {times[-1]:g}$ s")
fig.tight_layout(rect=(0, 0.05, 1, 1))
save(fig, out)
