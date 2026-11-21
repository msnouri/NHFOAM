#!/usr/bin/env python3
"""Submerged bar (Luth et al. 1994, cases A and C): gauges of NHFOAM (and of
NHWAVE / REEF3D::NHFLOW runs of the same set-up) vs lab data.

usage: plot_bar.py [--caseC] [--tend <t>] <labDataDir> <out.png> [kind:]<case>[:label] ...
    kind = nhfoam (default) | nhwave | reef3d   (see other_models.py)
    --tend: end every run at time t (same window for all models; the slow
            long-wave motion of the flume changes the phase over the bar)
(also writes out.pdf)

Lab data: fengyanshi/BENCHMARK_FUNWAVE, car_luth_AC/lab_data/CaseA (CaseC).
For each run one time shift (fitted at the first gauge, x = 22 m) is applied
to all gauges; the comparison window is the last 4-5 s of the run.
Printed: rms error with that shift, and the wave height (max - min over the
lab record) error per gauge.
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save
import other_models as om

use_latex(10)

args = sys.argv[1:]
T = 2.02
title = r"Submerged bar, case A ($H = 2$ cm, $T = 2.02$ s)"
if "--caseC" in args:
    args.remove("--caseC")
    T = 1.01
    title = r"Submerged bar, case C ($H = 4.1$ cm, $T = 1.01$ s)"
tEnd = None
if "--tend" in args:
    i = args.index("--tend")
    tEnd = float(args[i + 1])
    del args[i:i + 2]
labdir, out = args[0:2]
runs = [om.parse_case(a) for a in args[2:]]
files = {22: "Exp_22m.txt", 24: "Exp_24m.txt", 30.5: "Exp_30.5m.txt", 32.5: "Exp_32.5m.txt",
         33.5: "Exp_33.5m.txt", 34.5: "Exp_34.5m.txt", 35.7: "Exp_35.7m.txt", 37.3: "Exp_37.3m.txt",
         39.0: "Exp_39.0m.txt", 41.0: "Exp_41m.txt"}
lab = {x: np.loadtxt(os.path.join(labdir, f)) for x, f in files.items()}
xs = list(files)
cols = ["#1f5fa8", "#e67e22", "#2e8b57", "#8e44ad", "#7f8c8d"]
styles = ["-", "--", "-.", ":", (0, (1, 1))]

allLab = np.concatenate([v[:, 1] for v in lab.values()])
labTmax = max(v[:, 0].max() for v in lab.values())
ylo, yhi = 100 * 1.6 * allLab.min(), 100 * 1.35 * allLab.max()
fig, axes = plt.subplots(5, 2, figsize=(11, 12.5), sharex=True)
report = []
for (kind, path, label), col, ls in zip(runs, cols, styles):
    t, eta = om.read_any(kind, path, xs)
    if tEnd is not None:
        keep = t <= tEnd + 1e-9
        t = t[keep]
        eta = {x: e[keep] for x, e in eta.items()}
    # window at the end of the run (steady state also at the last gauges:
    # in case C the wave front needs about 40 s to reach x = 41 m)
    win0 = t[-1] - labTmax - 0.5 - T
    L0 = lab[22]
    best, bestErr = 0.0, np.inf
    for s in np.linspace(0, T, 404, endpoint=False):
        ei = np.interp(L0[:, 0], t - win0 - s, eta[22])
        err = np.mean((ei - L0[:, 1]) ** 2)
        if err < bestErr:
            best, bestErr = s, err
    rms, dH = [], []
    for k, x in enumerate(xs):
        ax = axes.flat[k]
        tm = t - win0 - best
        m = (tm >= 0) & (tm <= lab[x][:, 0].max() + 0.1)
        ax.plot(tm[m], 100 * eta[x][m], ls, color=col, lw=1.4, label=label)
        ei = np.interp(lab[x][:, 0], tm, eta[x])
        rms.append(100 * np.sqrt(np.mean((ei - lab[x][:, 1]) ** 2)))
        Hl = lab[x][:, 1].max() - lab[x][:, 1].min()
        dH.append(((eta[x][m].max() - eta[x][m].min()) / Hl - 1) * 100)
    report.append((label, rms, dH))
for k, x in enumerate(xs):
    ax = axes.flat[k]
    ax.plot(lab[x][:, 0], 100 * lab[x][:, 1], "o", ms=3.5, mfc="none", color="#c0392b",
            label="Luth et al. (1994)" if k == 0 else None)
    ax.set_title(rf"$x = {x:g}$ m")
    ax.set_ylim(ylo, yhi)
    ax.grid(alpha=0.3)
    if k % 2 == 0:
        ax.set_ylabel(r"$\eta$ [cm]")
for ax in axes[-1]:
    ax.set_xlabel(r"$t$ [s] (one time shift per run)")
axes.flat[0].legend(loc="upper right", fontsize=8)
fig.suptitle(title)
fig.tight_layout()
save(fig, out, dpi=130)
w = 22
print("rms error [cm] / wave height error [%] per gauge:")
print("x      " + "  ".join(f"{lbl[:w]:>{w}s}" for lbl, _, _ in report))
for k, x in enumerate(xs):
    print(f"{x:5.1f}  " + "  ".join(f"{r[k]:9.3f} / {h[k]:+7.1f}%   " for _, r, h in report))
print("mean   " + "  ".join(f"{np.mean(r):9.3f} / {np.mean(np.abs(h)):7.1f}%   " for _, r, h in report))
