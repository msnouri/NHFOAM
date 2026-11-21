#!/usr/bin/env python3
"""Elliptic shoal (Berkhoff et al. 1982; Ma et al. 2012, section 4.4): wave
height along the eight measurement sections, NHFOAM and NHWAVE against the
lab data, and the normalised rms error of Ma et al. per section.

usage: plot_ma2012_berkhoff.py <labDir> <out.png> [kind:]<case>[:label] ...
    kind = nhfoam (default) | nhwave; the gauges/stations are those of
    ma2012.berk_points() (ma2012.berkhoff_nhfoam / berkhoff_nhwave)
    labDir: NHWAVE examples/berkhoff/postprocessing (section1-5.dat: y and
            amplitude [mm]; section678.dat: -x and the amplitudes of 6-8)
(also writes out.pdf)

Wave height: mean of max - min over each of the last five periods.
"""
import os, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import ma2012
import other_models as om

labdir, out = sys.argv[1:3]
runs = [om.parse_case(a) for a in sys.argv[3:]]
pts = ma2012.berk_points()
T = ma2012.BERK_T


def heights(t, eta, nwave=5):
    t1 = t[-1]
    H = np.zeros(eta.shape[1])
    for m in range(nwave):
        sel = (t > t1 - (m + 1) * T) & (t <= t1 - m * T)
        H += eta[sel].max(axis=0) - eta[sel].min(axis=0)
    return H / nwave


def read(kind, path):
    if kind == "nhwave":
        cols = [np.loadtxt(os.path.join(path, "results", f"probe_{i + 1:04d}")) for i in range(len(pts))]
        n = min(len(c) for c in cols)
        t = cols[0][:n, 0]
        return t, np.column_stack([c[:n, 1] for c in cols])
    for sub in ("", "processor0"):
        f = os.path.join(path, sub, "postProcessing/gauges/0/eta.dat")
        if os.path.exists(f):
            d = np.loadtxt(f, comments="#")
            return d[:, 0], d[:, 1:]
    raise FileNotFoundError(path)


lab = {k: np.loadtxt(os.path.join(labdir, f"section{k}.dat")) for k in range(1, 6)}
l678 = np.loadtxt(os.path.join(labdir, "section678.dat"))
for j, k in enumerate((6, 7, 8)):
    lab[k] = np.column_stack([-l678[:, 0], l678[:, j + 1]])
for k in lab:
    lab[k] = lab[k][np.argsort(lab[k][:, 0])]
    lab[k][:, 1] *= 2e-3            # amplitude [mm] -> height [m]

sec = np.array([p[0] for p in pts])
pos = np.array([p[2] if p[0] <= 5 else p[1] for p in pts])   # y (1-5) or x (6-8)
H0 = ma2012.BERK_H

fig, axes = plt.subplots(4, 2, figsize=(7.5, 9.6), sharey=False)
styles = [("-", "#1f5fa8"), ("--", "#e67e22"), ("-.", "#7f8c8d")]
report = []
for (kind, path, label), (ls, col) in zip(runs, styles):
    t, eta = read(kind, path)
    H = heights(t, eta)
    errs = []
    for k in range(1, 9):
        m = sec == k
        ax = axes.flat[k - 1]
        ax.plot(pos[m], H[m] / H0, ls, color=col, lw=1.5, label=label)
        Hi = np.interp(lab[k][:, 0], pos[m], H[m])
        errs.append(np.sqrt(np.mean((Hi - lab[k][:, 1]) ** 2)) / lab[k][:, 1].mean())
    report.append((label, errs, t[-1]))
for k in range(1, 9):
    ax = axes.flat[k - 1]
    ax.plot(lab[k][:, 0], lab[k][:, 1] / H0, "o", ms=3.5, mfc="none", color="#c0392b",
            label="Berkhoff et al. (1982)" if k == 1 else None)
    if k <= 5:
        ax.set_title(rf"section {k}: $x = {ma2012.BERK_SEC_X[k - 1]:g}$ m", fontsize=11)
        ax.set_xlabel(r"$y$ [m]")
    else:
        ax.set_title(rf"section {k}: $y = {ma2012.BERK_SEC_Y[k - 6]:g}$ m", fontsize=11)
        ax.set_xlabel(r"$x$ [m]")
    ax.set_ylabel(r"$H/H_0$")
    ax.grid(alpha=0.3)
h, l = axes.flat[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=10, frameon=False)
fig.suptitle(r"Elliptic shoal (Berkhoff et al. 1982), $T = 1$ s, $H_0 = 4.64$ cm")
fig.tight_layout(rect=(0, 0.035, 1, 1))
save(fig, out)

print("normalised rms error of H per section (Ma et al. 2012 definition):")
print("run" + " " * 22 + "  ".join(f"  {k}  " for k in range(1, 9)) + "   mean   (t_end)")
for label, e, te in report:
    print(f"{label[:24]:24s} " + "  ".join(f"{x:.3f}" for x in e) + f"   {np.mean(e):.3f}   ({te:.1f} s)")
