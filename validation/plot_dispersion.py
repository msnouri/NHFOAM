#!/usr/bin/env python3
"""Relative period error vs kh for different numbers of layers (standing waves).

usage: plot_dispersion.py <dispersion.csv> <out.png> [title]   (also writes out.pdf)

The magnitude of the error is plotted; open symbols mark periods shorter than
linear theory.
"""
import sys, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex()
rows = list(csv.DictReader(open(sys.argv[1])))
out = sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else "Linear dispersion: standing waves, 64 cells per wavelength"
layers = sorted({int(r["layers"]) for r in rows})
cols = ["#7f3c8d", "#11a579", "#3969ac", "#f2b701", "#e73f74", "#80ba5a"]

fig, ax = plt.subplots(figsize=(6.0, 4.0))
for c, n in zip(cols, layers):
    rr = sorted([r for r in rows if int(r["layers"]) == n], key=lambda r: float(r["kh"]))
    kh = [float(r["kh"]) for r in rr]
    err = [100 * (float(r["ratio"]) - 1) for r in rr]
    ax.semilogy(kh, [abs(e) for e in err], "-", color=c, lw=1.5,
                label=f"{n} layer" + ("s" if n > 1 else ""))
    for k, e in zip(kh, err):
        ax.semilogy([k], [abs(e)], "o", color=c, ms=4.5,
                    mfc=(c if e > 0 else "white"), mew=1.2)
ax.axhline(1, color="0.5", lw=0.8, ls="--")
ax.set_xlabel(r"$kh$")
ax.set_ylabel(r"Period error, $|T_{\mathrm{NH}}/T_{\mathrm{lin}} - 1|$ $[\%]$")
ax.set_title(title)
ax.plot([], [], "o", color="0.3", mfc="white", mew=1.2, ms=4.5, label=r"$T_{\mathrm{NH}} < T_{\mathrm{lin}}$")
ax.grid(alpha=0.3, which="both")
ax.legend(ncol=2)
fig.tight_layout()
save(fig, out)
