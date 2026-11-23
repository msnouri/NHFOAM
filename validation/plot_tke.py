#!/usr/bin/env python3
"""Ting & Kirby (1994) spilling breakers: time-mean turbulent kinetic energy
of RANS runs (averaging fields (k nut): kMean) at the measuring stations.

usage: plot_tke.py <fig7.dat> <out.png> <case>[:label] ...
    fig7.dat: station, x_TK [m], (x - xb)/hb, h/hb, (z - mean level)/h,
              (k/(g h))^(1/2) (Ting & Kirby's Fig. 7a, digitised; not
              distributed with NHFOAM)
The models are averaged at fixed sigma and plotted at the mean height of
each level, as in plot_undertow.py.  Printed: the depth-mean of
(k/(g h))^(1/2) below the trough level at each station, model and lab.
(also writes out.pdf)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import longshore as ls
import wetdry

G = 9.81
labfile, out = sys.argv[1], sys.argv[2]
runs = []
for a in sys.argv[3:]:
    p = a.split(":", 1)
    runs.append((p[0], p[1] if len(p) > 1 else os.path.basename(p[0])))
lab = np.loadtxt(labfile, comments="#")
xst = sorted(set(lab[:, 1]))


def depth(x):
    return np.where(x < 0, wetdry.TK_D, wetdry.TK_D - wetdry.TK_SLOPE * x)


def profiles(case):
    t = ls.latest_time(case)
    C = ls.cell_centres(case, t)
    km = ls._read_field(os.path.join(case, t, "kMean"))
    em = ls._read_field(os.path.join(case, t, "etaMean"))
    e = ls._read_field(os.path.join(case, t, "eta"))
    xr = np.round(C[:, 0], 6)
    xs = np.unique(xr)
    prof = {}
    for xs0 in xst:
        x = xs[np.argmin(abs(xs - xs0))]
        m = np.where(xr == x)[0]
        m = m[np.argsort(C[m, 2])]
        h = float(depth(np.array(x)))
        zb = -h
        sig = (C[m, 2] - zb) / (e[m[0]] - zb)
        zbar = zb + sig * (em[m[0]] - zb)
        prof[xs0] = ((zbar - em[m[0]]) / h, np.sqrt(np.maximum(km[m], 0) / (G * h)))
    return prof


data = [(label, profiles(case)) for case, label in runs]
n = len(xst)
kmax = max([lab[:, 5].max()] + [p[x][1].max() for _, p in data for x in xst])
fig, axes = plt.subplots(1, n, figsize=(7.5, 3.6), sharey=True)
cols = ["#1f5fa8", "#e67e22", "#2e8b57", "#8e44ad"]
for ax, x in zip(np.atleast_1d(axes), xst):
    L = lab[lab[:, 1] == x]
    ax.plot(L[:, 5], L[:, 4], "o", ms=3.5, mfc="none", color="#c0392b", label="lab")
    line = f"x = {x:.2f} m: lab {L[:, 5].mean():.3f}"
    for (label, prof), col in zip(data, cols):
        zz, kk = prof[x]
        ax.plot(kk, zz, "-", color=col, lw=1.3, label=label)
        below = (zz > -1) & (zz < L[:, 4].max())
        line += f", {label} {kk[below].mean():.3f}"
    print(line)
    ax.set_title(rf"$x = {x:.2f}$ m", fontsize=9)
    ax.set_xlim(0, 1.05 * kmax)
    ax.grid(alpha=0.3)
    ax.set_xlabel(r"$(k/gh)^{1/2}$")
np.atleast_1d(axes)[0].set_ylabel(r"$(z - \bar\eta)/h$")
np.atleast_1d(axes)[0].set_ylim(-1.0, 0.1)
h_, l_ = np.atleast_1d(axes)[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=len(l_), fontsize=9, frameon=False)
fig.tight_layout(rect=(0, 0.08, 1, 1))
save(fig, out)
