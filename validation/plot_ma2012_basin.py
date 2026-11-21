#!/usr/bin/env python3
"""Standing wave in a closed basin (Ma et al. 2012, section 4.1): normalised
rms error at x = 17.5 m against the number of layers, for kD = pi, 2pi, 3pi.

usage: plot_ma2012_basin.py <runDir> <out.png> [prefix]
    runDir holds nf_kD<m>_L<n> (NHFOAM, limitedLinearV),
    <prefix>_kD<m>_L<n> (NHFOAM, filteredLinear2V; prefix nfc by default) and
    nw_kD<m>_L<n> (NHWAVE) for m = 1, 2, 3 and n = 3, 5, 8, 10.
(also writes out.pdf)
"""
import os, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import ma2012

run, out = sys.argv[1:3]
pre = sys.argv[3] if len(sys.argv) > 3 else "nfc"
layers = [3, 5, 8, 10]
depths = {1: 10.0, 2: 20.0, 3: 30.0}
colors = {1: "#1f5fa8", 2: "#e67e22", 3: "#2e8b57"}


def nf(case, D):
    d = np.loadtxt(os.path.join(run, case, "postProcessing/gauges/0/eta.dat"), comments="#")
    return ma2012.basin_error(d[:, 0], d[:, 2], D)


def nw(case, D):
    d = np.loadtxt(os.path.join(run, case, "results/probe_0002"))
    return ma2012.basin_error(d[:, 0], d[:, 1], D)


fig, ax = plt.subplots(figsize=(6.4, 4.6))
print("kD   layers  paper   NHWAVE  NHFOAM  NHFOAM(central)  [%]")
for m, D in depths.items():
    c = colors[m]
    e_nf = [nf(f"nf_kD{m}_L{n}", D) for n in layers]
    e_nfc = [nf(f"{pre}_kD{m}_L{n}", D) for n in layers]
    e_nw = [nw(f"nw_kD{m}_L{n}", D) for n in layers]
    e_pa = [ma2012.BASIN_PAPER[m][n] for n in layers]
    lab = {1: r"\pi", 2: r"2\pi", 3: r"3\pi"}[m]
    ax.plot(layers, e_pa, "o", color=c, mfc="white", mew=1.3, ms=7)
    ax.plot(layers, e_nw, "--", color=c, lw=1.4)
    ax.plot(layers, e_nf, "-", color=c, lw=1.8, marker="s", ms=4, label=f"$kD = {lab}$")
    ax.plot(layers, e_nfc, ":", color=c, lw=1.8)
    for n, a, b, d1, d2 in zip(layers, e_pa, e_nw, e_nf, e_nfc):
        print(f"{m}pi  {n:3d}    {100*a:5.2f}   {100*b:5.2f}   {100*d1:5.2f}   {100*d2:5.2f}")
ax.plot([], [], "o", color="0.3", mfc="white", mew=1.3, label="Ma et al. (2012), NHWAVE")
ax.plot([], [], "--", color="0.3", label="NHWAVE, rerun here")
ax.plot([], [], "-", color="0.3", label=r"NHFOAM (\texttt{limitedLinearV})")
ax.plot([], [], ":", color="0.3", lw=1.8, label=r"NHFOAM (\texttt{filteredLinear2V})")
ax.set_xlabel("Number of layers")
ax.set_ylabel(r"$\sqrt{\overline{(\eta_a-\eta)^2}}\,/\,H$ at $x = 17.5$ m")
ax.set_xlim(2.5, 10.5)
ax.set_ylim(0, 0.052)
ax.set_xticks(layers)
ax.grid(alpha=0.3)
ax.legend(fontsize=10, ncol=2, loc="upper right")
ax.set_title(r"Standing wave in a closed basin, $L = 20$ m, $a = 0.1$ m, 30 s")
fig.tight_layout()
save(fig, out)
