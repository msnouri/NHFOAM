#!/usr/bin/env python3
"""Solitary wave (H/h = 0.1) profiles vs the initial sech^2 shape moved at c = sqrt(g(h+H))."""
import sys, os, re, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
case, out = sys.argv[1], sys.argv[2]
H, h, x0, L, dx = 0.1, 1.0, 10.0, 60.0, 0.05
Nx = int(round(L / dx)); x = (np.arange(Nx) + 0.5) * dx
K = math.sqrt(3 * H / (4 * h ** 3)); c = math.sqrt(9.81 * (h + H))
times = sorted([d for d in os.listdir(case) if re.match(r"^[0-9.]+$", d) and float(d) > 0], key=float)
fig, ax = plt.subplots(figsize=(9, 3.6))
for tn in ["0"] + times:
    if tn == "0":
        e = H / np.cosh(K * (x - x0)) ** 2
    else:
        txt = open(f"{case}/{tn}/eta").read().split("internalField")[1].split(";")[0]
        vals = np.array(re.findall(r"[-0-9.eE+]+", txt.split("(", 1)[1]), float)[:Nx]
        e = vals
    ax.plot(x, 100 * e, lw=1.6, label=f"t = {float(tn):g} s")
    ax.plot(x, 100 * H / np.cosh(K * (x - x0 - c * float(tn))) ** 2, "k:", lw=0.9)
ax.plot([], [], "k:", label="sech² moved at √(g(h+H))")
ax.set_xlim(0, 50); ax.set_xlabel("x [m]"); ax.set_ylabel("η [cm]"); ax.grid(alpha=0.3)
ax.set_title("Solitary wave, H/h = 0.1, 10 layers, Δx = 0.05 m")
ax.legend(fontsize=8, ncol=3)
fig.tight_layout(); fig.savefig(out, dpi=140)
