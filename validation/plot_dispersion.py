#!/usr/bin/env python3
"""Relative period error vs kh for different numbers of layers (standing waves).

usage: plot_dispersion.py <dispersion.csv> <out.png>   (also writes out.pdf)
"""
import sys, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex(10)
rows = list(csv.DictReader(open(sys.argv[1])))
out = sys.argv[2]
layers = sorted({int(r["layers"]) for r in rows})
cols = ["#7f3c8d", "#11a579", "#3969ac", "#f2b701", "#e73f74", "#80ba5a"]

fig, ax = plt.subplots(figsize=(6.0, 4.0))
for c, n in zip(cols, layers):
    rr = sorted([r for r in rows if int(r["layers"]) == n], key=lambda r: float(r["kh"]))
    kh = [float(r["kh"]) for r in rr]
    err = [100 * (float(r["ratio"]) - 1) for r in rr]
    ax.semilogy(kh, err, "o-", color=c, lw=1.5, ms=4,
                label=f"{n} layer" + ("s" if n > 1 else ""))
ax.axhline(1, color="0.5", lw=0.8, ls="--")
ax.set_xlabel(r"$kh$")
ax.set_ylabel(r"Period error, $T_{\mathrm{NH}}/T_{\mathrm{lin}} - 1$ $[\%]$")
ax.set_title("Linear dispersion: standing waves, 64 cells per wavelength")
ax.grid(alpha=0.3, which="both")
ax.legend(ncol=2)
fig.tight_layout()
save(fig, out)
