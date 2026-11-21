#!/usr/bin/env python3
"""Circular basin: unstructured surface mesh and gauge time series vs linear theory.

usage: plot_basin.py <case> <out.png>   (case = a run of tutorials/circularBasin;
also writes out.pdf)
"""
import sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.special import jv
from latex_style import use_latex, save

use_latex()
case, out = sys.argv[1], sys.argv[2]
R, h, a = 1.0, 0.5, 0.005
k = 1.841183781
w = math.sqrt(9.81 * k * math.tanh(k * h))
T = 2 * math.pi / w

# --- top-surface triangles from the gmsh file (physical group 'surface')
lines = open(f"{case}/basin.msh").read().split("\n")
i = lines.index("$PhysicalNames"); n = int(lines[i + 1])
phys = {l.split()[2].strip('"'): int(l.split()[1]) for l in lines[i + 2:i + 2 + n]}
i = lines.index("$Nodes"); n = int(lines[i + 1])
nodes = {int(l.split()[0]): list(map(float, l.split()[1:4])) for l in lines[i + 2:i + 2 + n]}
i = lines.index("$Elements"); n = int(lines[i + 1])
tris = []
for l in lines[i + 2:i + 2 + n]:
    s = list(map(int, l.split()))
    if s[1] == 2 and s[3] == phys["surface"]:
        nt = s[2]
        tris.append(s[3 + nt:])
polys = [[nodes[v][:2] for v in t] for t in tris]
cent = np.array([np.mean(p, axis=0) for p in polys])
r = np.hypot(cent[:, 0], cent[:, 1])
th = np.arctan2(cent[:, 1], cent[:, 0])
eta0 = a * jv(1, k * r) * np.cos(th)

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(7.5, 4.0),
                               gridspec_kw={"width_ratios": [1, 1.6]})
pc = PolyCollection(polys, array=100 * eta0, cmap="RdBu_r",
                    edgecolors="0.35", linewidths=0.25)
ax0.add_collection(pc)
ax0.set_xlim(-1.05, 1.05); ax0.set_ylim(-1.05, 1.05); ax0.set_aspect("equal")
ax0.set_xlabel(r"$x$ [m]"); ax0.set_ylabel(r"$y$ [m]")
ax0.set_title(f"{len(tris)} triangles," + "\n" + r"initial $\eta$ [cm]")
fig.colorbar(pc, ax=ax0, shrink=0.8)

d = np.loadtxt(f"{case}/postProcessing/gauges/0/eta.dat", comments="#")
t = d[:, 0]
gauges = [(r"$r = 0.95$ m, $\theta = 0$", "#3969ac"),
          (r"$r = 0.5$ m, $\theta = 0$", "#11a579"),
          (r"$r = 0.95$ m, $\theta = 45^\circ$", "#e73f74")]
# r = 0.5 m and theta = 45 deg have almost the same amplitude: draw the first
# as a wider band underneath so both stay visible
widths = [1.4, 3.2, 1.2]
for j, (name, c) in enumerate(gauges):
    e = d[:, j + 1]
    ax1.plot(t, 100 * e, color=c, lw=widths[j], label="NHFOAM, " + name,
             alpha=0.55 if j == 1 else 1.0)
    ax1.plot(t, 100 * e[0] * np.cos(w * t), "k--", lw=0.8)
ax1.plot([], [], "k--", lw=0.8, label=rf"linear theory, $T = {T:.4f}$ s")
ax1.set_xlabel(r"$t$ [s]"); ax1.set_ylabel(r"$\eta$ [cm]")
ax1.grid(alpha=0.3)
h, l = ax1.get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=2, fontsize=10, frameon=False)
ax1.set_title(r"First sloshing mode, $h = 0.5$ m, 8 layers")
fig.tight_layout(rect=(0, 0.14, 1, 1))
save(fig, out)
