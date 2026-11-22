#!/usr/bin/env python3
"""Solitary wave on a plane beach (Synolakis 1987; NTHMP benchmark 4; Ma et
al. 2012, 4.5): free-surface profiles of NHFOAM and NHWAVE against the lab
profiles, and the run-up of NHFOAM over time.

usage: plot_runup.py <labProfDir> <labRunup.txt> <out.png> <H> [kind:]<case>[:label] ...
    kind = nhfoam (default) | nhwave; cases made by wetdry.runup_nhfoam /
    wetdry.runup_nhwave (fields every 5 t* or more often, t* = t sqrt(g/d));
    also plots the run-up over time (out_shoreline.png): NHFOAM from
    postProcessing/wettingDrying, NHWAVE from its surface fields (only
    cases with 20 or more of them)
    labProfDir: NTHMP BP04 profs (files Case0_3.t=<t*>: x/d offshore of the
    initial shoreline, eta/d)
(also writes out.pdf)
"""
import os, re, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import wetdry
import other_models as om
import make_case as mc

labdir, labrunup, out, H = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
runs = [om.parse_case(a) for a in sys.argv[5:]]
d = wetdry.RU_D
ts = np.sqrt(d / wetdry.G)
xc0, xshore = wetdry.runup_geometry(H, d)
TSTARS = [10, 15, 20, 25, 30, 50]
film = 0.005


def nhfoam_profiles(case, nz=3):
    xs, hs, zt = mc.beach_profile([0.0, wetdry.RU_TOE, wetdry.RU_LEN],
                                  [d, d, d - (wetdry.RU_LEN - wetdry.RU_TOE) * wetdry.RU_SLOPE], film)
    xcell = wetdry.profile_cell_x(xs, wetdry.RU_DX, nz)
    first = np.zeros(len(xcell), bool)        # top-layer cells: one per column
    pos = 0
    for i in range(len(xs) - 1):
        nx = max(1, int(round((xs[i + 1] - xs[i]) / wetdry.RU_DX)))
        first[pos:pos + nx] = True
        pos += nx * nz
    x = xcell[first]
    out = {}
    for t in os.listdir(case):
        try:
            tv = float(t)
        except ValueError:
            continue
        if tv == 0:
            continue
        s = open(os.path.join(case, t, "eta")).read().split("internalField")[1]
        e = np.array(s.split("(", 1)[1].split(")", 1)[0].split(), float)[first]
        h = wetdry.runup_depth(x, d)
        e = np.where(e + h > 2 * film, e, np.nan)
        out[round(tv / ts)] = (x, e)
    return out


def nhwave_profiles(case):
    tt = np.atleast_1d(np.loadtxt(os.path.join(case, "results", "time")))
    out = {}
    for n, tv in enumerate(tt, start=1):
        e = np.loadtxt(os.path.join(case, "results", f"eta_{n:05d}")).ravel()
        x = (np.arange(len(e)) + 0.5) * wetdry.RU_DX
        D = e + wetdry.runup_depth(x, d)
        e = np.where(D > 2 * film, e, np.nan)
        out[round(tv / ts)] = (x, e)
    return out


data = []
for kind, path, label in runs:
    prof = nhwave_profiles(path) if kind == "nhwave" else nhfoam_profiles(path)
    data.append((label, prof))

fig, axes = plt.subplots(len(TSTARS), 1, figsize=(7.5, 10.5), sharex=True)
styles = [("-", "#1f5fa8"), ("--", "#e67e22"), ("-.", "#7f8c8d")]
xb = np.linspace(-6, 30, 400)
for ax, tsr in zip(axes, TSTARS):
    # bed
    xx = xshore - xb * d
    ax.fill_between(xb, -0.2 * np.ones_like(xb), -wetdry.runup_depth(xx, d) / d,
                    color="0.85", lw=0)
    f = os.path.join(labdir, f"Case0_3.t={tsr}")
    if os.path.exists(f):
        L = np.loadtxt(f)
        ax.plot(L[:, 0], L[:, 1], "o", ms=3, mfc="none", color="#c0392b",
                label="Synolakis (1987)")
    for (label, prof), (ls, col) in zip(data, styles):
        if tsr in prof:
            x, e = prof[tsr]
            ax.plot((xshore - x) / d, e / d, ls, color=col, lw=1.4, label=label)
    ax.set_ylim(-0.1, 0.65)
    ax.set_ylabel(r"$\eta/d$")
    ax.text(0.98, 0.85, rf"$t\sqrt{{g/d}} = {tsr}$", transform=ax.transAxes, ha="right")
    ax.grid(alpha=0.3)
axes[-1].set_xlim(-6, 30)
axes[-1].invert_xaxis()
axes[-1].set_xlabel(r"$x/d$ (offshore of the initial shoreline)")
h, l = [], []
for ax in axes:
    for hh, ll in zip(*ax.get_legend_handles_labels()):
        if ll not in l:
            h.append(hh); l.append(ll)
fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=10, frameon=False)
fig.suptitle(rf"Solitary wave on a 1:20 beach, $H/d = {H/d:.2f}$")
fig.tight_layout(rect=(0, 0.04, 1, 1))
save(fig, out)

# shoreline (highest wet bed) over time
lab = np.loadtxt(labrunup, comments="#")
near = np.abs(lab[:, 0] - H / d) < 0.02
print(f"lab, H/d = {H/d - 0.02:.2f}-{H/d + 0.02:.2f}: R/d = " + ", ".join(f"{r:.3f}" for r in lab[near, 1]))


def nhwave_shoreline(case):
    res = os.path.join(case, "results")
    tt = np.atleast_1d(np.loadtxt(os.path.join(res, "time")))
    dep = np.loadtxt(os.path.join(res, "depth")).ravel()
    R = []
    for n in range(1, len(tt) + 1):
        f = os.path.join(res, f"eta_{n:05d}")
        if not os.path.exists(f):
            break
        e = np.loadtxt(f).ravel()
        R.append(-dep[e + dep > film + 1e-5].min())
    return tt[:len(R)], np.array(R)


fig, ax = plt.subplots(figsize=(7.5, 3.8))
ax.axhspan(lab[near, 1].min(), lab[near, 1].max(), color="#c0392b", alpha=0.15, lw=0,
           label=f"lab maximum, $H/d = {H/d - 0.02:.2f}$--${H/d + 0.02:.2f}$")
for (kind, path, label), (ls, col) in zip(runs, styles + [(":", "#2e8b57"), ((0, (5, 1)), "#8e44ad")]):
    if kind == "nhwave":
        t, R = nhwave_shoreline(path)
        if len(t) < 20:
            continue
        print(f"{label}: maximum run-up R/d = {R.max() / d:.3f} at t* = {t[np.argmax(R)] / ts:.1f}")
    else:
        f = os.path.join(path, "postProcessing/wettingDrying/0/shoreline.dat")
        if not os.path.exists(f):
            continue
        sh = np.loadtxt(f, comments="#")
        t, R = sh[:, 0], sh[:, 2]
        print(f"{label}: maximum run-up R/d = {R.max() / d:.3f} at t* = {t[np.argmax(R)] / ts:.1f}"
              f" (highest surface on land {sh[:, 3].max() / d:.3f})")
    ax.plot(t / ts, R / d, ls=ls, color=col, lw=1.4, label=label)
ax.set_xlim(15, 60)
ax.set_ylim(-0.2, 0.65)
ax.set_xlabel(r"$t\sqrt{g/d}$")
ax.set_ylabel(r"$R/d$ (highest wet bed)")
ax.grid(alpha=0.3)
ax.legend(fontsize=9, frameon=False, ncol=2, loc="lower center")
fig.tight_layout()
save(fig, out.rsplit(".", 1)[0] + "_shoreline.png")
