#!/usr/bin/env python3
"""Submerged bar (Luth et al. 1994 case A): NHFoam gauges vs lab data.

usage: plot_bar.py <labDataDir> <out.png> <case>[:label] [<case>[:label] ...]
(also writes out.pdf)

Lab data: fengyanshi/BENCHMARK_FUNWAVE, car_luth_AC/lab_data/CaseA.
For each case one time shift (fitted at the first gauge, x = 22 m) is applied
to all gauges; the window starts 6 periods before the end of the run.
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex(10)

labdir, out = sys.argv[1:3]
cases = [(c.split(":")[0], c.split(":")[1] if ":" in c else os.path.basename(c)) for c in sys.argv[3:]]
T = 2.02
files = {22: "Exp_22m.txt", 24: "Exp_24m.txt", 30.5: "Exp_30.5m.txt", 32.5: "Exp_32.5m.txt",
         33.5: "Exp_33.5m.txt", 34.5: "Exp_34.5m.txt", 35.7: "Exp_35.7m.txt", 37.3: "Exp_37.3m.txt",
         39.0: "Exp_39.0m.txt", 41.0: "Exp_41m.txt"}
lab = {x: np.loadtxt(os.path.join(labdir, f)) for x, f in files.items()}
xs = list(files)
cols = ["#1f5fa8", "#e67e22", "#2e8b57"]

fig, axes = plt.subplots(5, 2, figsize=(11, 12.5), sharex=True)
report = []
for (case, label), col in zip(cases, cols):
    fn = os.path.join(case, "postProcessing/gauges/0/eta.dat")
    d = np.loadtxt(fn, comments="#")
    t = d[:, 0]
    names = [l.split()[1] for l in open(fn) if l.startswith("# G")]
    gx = [float(n[1:]) for n in names]
    win0 = t[-1] - 6 * T
    L0 = lab[22]
    i0 = gx.index(22.0)
    best, bestErr = 0.0, np.inf
    for s in np.linspace(0, T, 404, endpoint=False):
        ei = np.interp(L0[:, 0], t - win0 - s, d[:, i0 + 1])
        err = np.mean((ei - L0[:, 1]) ** 2)
        if err < bestErr:
            best, bestErr = s, err
    rms = []
    for k, x in enumerate(xs):
        ax = axes.flat[k]
        i = gx.index(x)
        tm = t - win0 - best
        m = (tm >= 0) & (tm <= 2 * T + 0.1)
        ax.plot(tm[m], 100 * d[m, i + 1], "-", color=col, lw=1.5, label=label)
        ei = np.interp(lab[x][:, 0], tm, d[:, i + 1])
        rms.append(100 * np.sqrt(np.mean((ei - lab[x][:, 1]) ** 2)))
    report.append((label, rms))
for k, x in enumerate(xs):
    ax = axes.flat[k]
    ax.plot(lab[x][:, 0], 100 * lab[x][:, 1], "o", ms=3.5, mfc="none", color="#c0392b",
            label="Luth et al. (1994)" if k == 0 else None)
    ax.set_title(rf"$x = {x:g}$ m")
    ax.set_ylim(-1.7, 3.2)
    ax.grid(alpha=0.3)
    if k % 2 == 0:
        ax.set_ylabel(r"$\eta$ [cm]")
for ax in axes[-1]:
    ax.set_xlabel(r"$t$ [s] (one time shift per run)")
axes.flat[0].legend(loc="upper right")
fig.suptitle(r"Submerged bar, case A ($H = 2$ cm, $T = 2.02$ s)")
fig.tight_layout()
save(fig, out, dpi=130)
print("rms error [cm] per gauge:")
print("x      " + "  ".join(f"{lbl:>18s}" for lbl, _ in report))
for k, x in enumerate(xs):
    print(f"{x:5.1f}  " + "  ".join(f"{r[k]:18.3f}" for _, r in report))
