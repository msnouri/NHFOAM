#!/usr/bin/env python3
"""Progressive waves in a flat flume: wave height along the flume for
different convection schemes (numerical damping of div(phi,U)).

usage: plot_flume.py <T> <out.png> <case>:<label> [...]   (also writes out.pdf)
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex(10)
T, out = float(sys.argv[1]), sys.argv[2]
cases = [(c.split(":")[0], c.split(":")[1]) for c in sys.argv[3:]]
cols = ["#c0392b", "#e67e22", "#2e8b57", "#1f5fa8", "#7f8c8d"]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
for (case, label), col in zip(cases, cols):
    fn = os.path.join(case, "postProcessing/gauges/0/eta.dat")
    d = np.loadtxt(fn, comments="#")
    t = d[:, 0]
    xg = [float(l.split()[1][1:]) for l in open(fn) if l.startswith("# G")]
    m = t > t[-1] - 5 * T
    H = [d[m, i + 1].max() - d[m, i + 1].min() for i in range(len(xg))]
    ax.plot(xg, 100 * np.array(H), "o-", color=col, ms=4, label=label)
    print(f"{label}: H = " + " ".join(f"{100*h:.2f}" for h in H) + " cm")
ax.set_xlabel(r"$x$ [m] (wave-making zone $0$--$8$ m)")
ax.set_ylabel(r"$H$ [cm]")
ax.set_title(r"Flat flume, $h = 0.4$ m, $T = 1.01$ s, target $H = 4.1$ cm, 10 layers, $\Delta x = 2$ cm")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
save(fig, out, dpi=150)
