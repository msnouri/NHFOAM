#!/usr/bin/env python3
"""Ting & Kirby (1994) spilling breakers: crest and trough envelopes, wave
height and mean water level of NHFOAM (and NHWAVE) against the measurements.

usage: plot_tingkirby.py <fig4.dat> <out.png> [--periods N] [--nhwave-x0 X]
                         [kind:]<case>[:label] ...
    fig4.dat: x_TK [m], mean level [cm], crest - mean [cm], trough - mean
              [cm] (the measurements of Ting & Kirby's Fig. 4, digitised;
              not distributed with NHFOAM)
    kind = nhfoam (gauges of wetdry.tingkirby, x = x_TK) or nhwave (surface
           fields results/eta_*, x_TK = x - X, default X = 5.6 m: the toe of
           the NHWAVE example's bed)
    statistics over the last N wave periods (default 10)
(also writes out.pdf)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import other_models as om
import wetdry

args = sys.argv[1:]
nper, xnw = 10, 5.6
if "--periods" in args:
    i = args.index("--periods"); nper = int(args[i + 1]); del args[i:i + 2]
if "--nhwave-x0" in args:
    i = args.index("--nhwave-x0"); xnw = float(args[i + 1]); del args[i:i + 2]
labfile, out = args[0], args[1]
runs = [om.parse_case(a) for a in args[2:]]
T = wetdry.TK_T


def stats(t, e):
    """Per-period crest and trough over the last nper periods, mean level."""
    tend = t[-1]
    cr, tr = [], []
    for p in range(nper):
        m = (t > tend - (p + 1) * T) & (t <= tend - p * T)
        cr.append(e[m].max(axis=0)); tr.append(e[m].min(axis=0))
    sel = t > tend - nper * T
    return np.mean(cr, axis=0), np.mean(tr, axis=0), e[sel].mean(axis=0)


def nhfoam(case):
    f = os.path.join(case, "postProcessing", "gauges", "0", "eta.dat")
    xs = []
    for line in open(f):
        if line.startswith("# G"):
            xs.append(float(line.split("(")[1].split()[0]))
    d = np.loadtxt(f, comments="#")
    return np.array(xs), d[:, 0], d[:, 1:]


def nhwave(case):
    res = os.path.join(case, "results")
    files = sorted(f for f in os.listdir(res) if f.startswith("eta_"))
    if os.path.exists(os.path.join(res, "time")):
        tt = np.atleast_1d(np.loadtxt(os.path.join(res, "time")))[:len(files)]
    else:                                           # PLOT_START + n PLOT_INTV
        keys = dict(l.split("=", 1) for l in open(os.path.join(case, "input.txt"))
                    if "=" in l and not l.lstrip().startswith("!"))
        t0 = float(keys["PLOT_START "].split()[0]) if "PLOT_START " in keys else 0.0
        dtp = float(keys["PLOT_INTV "].split()[0])
        tt = t0 + dtp * np.arange(1, len(files) + 1)
    e = np.array([np.loadtxt(os.path.join(res, f)).ravel() for f in files])
    f = os.path.join(res, "depth")
    dep = np.loadtxt(f if os.path.exists(f) else os.path.join(case, "depth.txt")).ravel()
    x = (np.arange(e.shape[1]) + 0.5) * 0.025 - xnw
    e = np.where(e + dep > 0.0201, e, np.nan)       # dry: MinDep 0.02 m
    return x, tt, e


lab = np.loadtxt(labfile, comments="#")
lx, lmean, lcr, ltr = lab[:, 0], lab[:, 1] / 100, lab[:, 2] / 100, lab[:, 3] / 100
lH = lcr - ltr

data = []
for kind, path, label in runs:
    x, t, e = nhwave(path) if kind == "nhwave" else nhfoam(path)
    cr, tr, mn = stats(t, e)
    data.append((label, x, cr, tr, mn))

fig, axes = plt.subplots(3, 1, figsize=(7.5, 9.0), sharex=True,
                         gridspec_kw=dict(height_ratios=[1.3, 1, 1]))
styles = [("-", "#1f5fa8", 1.5), ("--", "#e67e22", 1.4), ("-.", "#7f8c8d", 1.4),
          (":", "#2e8b57", 1.6)]
ax = axes[0]
ax.plot(lx, lcr, "o", ms=4, mfc="none", color="#c0392b", label="Ting \\& Kirby (1994)"
        if plt.rcParams["text.usetex"] else "Ting & Kirby (1994)")
ax.plot(lx, ltr, "o", ms=4, mfc="none", color="#c0392b")
for (label, x, cr, tr, mn), (ls, c, lw) in zip(data, styles):
    ax.plot(x, cr - mn, ls=ls, color=c, lw=lw, label=label)
    ax.plot(x, tr - mn, ls=ls, color=c, lw=lw)
ax.set_ylabel(r"$\eta_{\max} - \bar\eta$, $\eta_{\min} - \bar\eta$ [m]")
ax = axes[1]
ax.plot(lx, lH, "o", ms=4, mfc="none", color="#c0392b")
for (label, x, cr, tr, mn), (ls, c, lw) in zip(data, styles):
    ax.plot(x, cr - tr, ls=ls, color=c, lw=lw)
ax.set_ylabel(r"$H$ [m]")
h, l = axes[0].get_legend_handles_labels()
ax.legend(h, l, fontsize=10, frameon=False, ncol=1, loc="lower left")
ax = axes[2]
ax.plot(lx, lmean, "o", ms=4, mfc="none", color="#c0392b")
for (label, x, cr, tr, mn), (ls, c, lw) in zip(data, styles):
    ax.plot(x, mn, ls=ls, color=c, lw=lw)
ax.set_ylabel(r"$\bar\eta$ [m]")
ax.set_xlabel(r"$x$ [m] (from the toe of the 1:35 slope)")
ax.set_xlim(-2, 13.5)
for a in axes:
    a.grid(alpha=0.3)
fig.tight_layout()
save(fig, out)

# summary
print("lab: largest H %.4f m at x = %.2f m" % (lH.max(), lx[np.argmax(lH)]))
pre, post = lx < 6.0, lx > 7.0
for label, x, cr, tr, mn in data:
    H = cr - tr
    ok = np.isfinite(H) & (x > -2) & (x < 13)
    i = np.nanargmax(np.where(ok, H, np.nan))
    Hl = np.interp(lx, x[ok], H[ok])
    ml = np.interp(lx, x[ok], mn[ok])
    cl = np.interp(lx, x[ok], (cr - mn)[ok])
    tl = np.interp(lx, x[ok], (tr - mn)[ok])
    print("%s: largest H %.4f m at x = %.2f m; H rms error %.2f cm before x = 6 m,"
          " %.2f cm after x = 7 m (bias %+.2f cm); crest %+.2f cm, trough %+.2f cm"
          " (mean errors); mean level rms %.2f cm (lab %.2f..%.2f cm)"
          % (label, H[i], x[i], 100 * np.sqrt(np.mean((Hl - lH)[pre] ** 2)),
             100 * np.sqrt(np.mean((Hl - lH)[post] ** 2)), 100 * np.mean((Hl - lH)[post]),
             100 * np.mean(cl - lcr), 100 * np.mean(tl - ltr),
             100 * np.sqrt(np.mean((ml - lmean) ** 2)), 100 * lmean.min(), 100 * lmean.max()))
    for xx in [0.0, 4.0, 6.0, 6.4, 7.0, 8.0, 9.0, 10.0, 11.0]:
        j = np.argmin(abs(lx - xx))
        print("   x %5.2f: H %.4f (lab %.4f at %.2f), mean %+.4f (lab %+.4f)"
              % (lx[j], Hl[j], lH[j], lx[j], ml[j], lmean[j]))
