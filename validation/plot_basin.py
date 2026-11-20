#!/usr/bin/env python3
"""Circular basin: unstructured surface mesh and gauge time series vs linear theory."""
import sys, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.special import jv
case, out = sys.argv[1], sys.argv[2]
R, h, a = 1.0, 0.5, 0.005
k = 1.841183781; w = math.sqrt(9.81 * k * math.tanh(k * h)); T = 2 * math.pi / w
# --- read top-surface triangles from the gmsh file (physical group 'surface')
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
        nt = s[2]; tris.append(s[3 + nt:])
polys = [[nodes[v][:2] for v in t] for t in tris]
cent = np.array([np.mean(p, axis=0) for p in polys])
r = np.hypot(cent[:, 0], cent[:, 1]); th = np.arctan2(cent[:, 1], cent[:, 0])
eta0 = a * jv(1, k * r) * np.cos(th)
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1, 1.6]})
pc = PolyCollection(polys, array=100 * eta0, cmap="RdBu_r", edgecolors="0.35", linewidths=0.25)
ax0.add_collection(pc); ax0.set_xlim(-1.05, 1.05); ax0.set_ylim(-1.05, 1.05); ax0.set_aspect("equal")
ax0.set_title(f"Unstructured surface mesh ({len(tris)} triangles)\ninitial η [cm]", fontsize=10)
fig.colorbar(pc, ax=ax0, shrink=0.8)
d = np.loadtxt(f"{case}/postProcessing/gauges/0/eta.dat", comments="#")
t = d[:, 0]
for j, (name, c) in enumerate([("r = 0.95, θ = 0", "#3969ac"), ("r = 0.5, θ = 0", "#11a579"), ("r = 0.95, θ = 45°", "#e73f74")]):
    e = d[:, j + 1]
    ax1.plot(t, 100 * e, color=c, lw=1.5, label=f"NHFoam, {name}")
    ax1.plot(t, 100 * e[0] * np.cos(w * t), "k--", lw=0.8)
ax1.plot([], [], "k--", lw=0.8, label=f"linear theory, T = {T:.4f} s")
ax1.set_xlabel("t [s]"); ax1.set_ylabel("η [cm]"); ax1.grid(alpha=0.3)
ax1.legend(fontsize=8, loc="lower left")
ax1.set_title("First sloshing mode, h = 0.5 m, 8 layers (prisms)", fontsize=10)
fig.tight_layout(); fig.savefig(out, dpi=140)
