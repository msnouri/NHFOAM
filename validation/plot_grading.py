#!/usr/bin/env python3
"""Period error of standing waves vs kh for uniform and surface-refined layers.

usage: plot_grading.py <grading.csv> <out.png> [title]   (from
run_grading_study.sh; also writes out.pdf)
"""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex(10)
d = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
cols = {1.0: "#1f5fa8", 0.5: "#e67e22", 0.25: "#2e8b57"}
for ax, nl in zip(axes, (3, 5)):
    for gr, col in cols.items():
        m = (d["layers"] == nl) & np.isclose(d["grading"], gr)
        kh, r = d["kh"][m], d["ratio"][m]
        o = np.argsort(kh)
        lbl = "uniform" if gr == 1 else rf"top/bottom = {gr:g}"
        ax.plot(kh[o], 100 * (r[o] - 1), "o-", color=col, ms=4, label=lbl)
    ax.axhline(1, color="0.5", lw=0.8, ls="--")
    ax.set_title(f"{nl} layers")
    ax.set_xlabel(r"$kh$")
    ax.grid(alpha=0.3)
axes[0].set_ylabel(r"period error [\%]")
axes[0].legend()
fig.suptitle(sys.argv[3] if len(sys.argv) > 3 else
             r"Standing waves: uniform vs surface-refined $\sigma$-layers")
fig.tight_layout()
save(fig, sys.argv[2], dpi=150)
