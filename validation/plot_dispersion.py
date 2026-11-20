#!/usr/bin/env python3
"""Relative period error vs kh for different numbers of layers (standing waves)."""
import sys, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open(sys.argv[1])))
out = sys.argv[2]
layers = sorted({int(r["layers"]) for r in rows})
cols = ["#7f3c8d", "#11a579", "#3969ac", "#f2b701", "#e73f74", "#80ba5a"]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for c, n in zip(cols, layers):
    rr = sorted([r for r in rows if int(r["layers"]) == n], key=lambda r: float(r["kh"]))
    kh = [float(r["kh"]) for r in rr]
    err = [100 * (float(r["ratio"]) - 1) for r in rr]
    ax.semilogy(kh, err, "o-", color=c, lw=1.6, ms=4, label=f"{n} layer{'s' if n > 1 else ''}")
ax.axhline(1, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("kh")
ax.set_ylabel("period error  (T_NHFoam / T_linear - 1)  [%]")
ax.set_title("Linear dispersion: standing waves, 64 cells per wavelength")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(out, dpi=140)
