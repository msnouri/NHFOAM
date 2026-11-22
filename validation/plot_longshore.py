#!/usr/bin/env python3
"""Longshore current on a plane beach (Visser 1991, case 4; Ma et al. 2012,
Section 4.7): wave height, mean water level and depth-averaged longshore
current against the distance from the still-water shoreline.

usage: plot_longshore.py <out.png> [--snapshot case] [kind:]<case>[:label] ...
    kind = nhfoam (validation/longshore.visser) or nhwave
           (longshore.nhwave_visser)
    --snapshot case: also draws the surface elevation of that NHFOAM case at
           its last time in plan view (out_snapshot.png)
(also writes the .pdf files)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import other_models as om
import longshore as ls

args = sys.argv[1:]
snap = None
if "--snapshot" in args:
    i = args.index("--snapshot"); snap = args[i + 1]; del args[i:i + 2]
out = args[0]
runs = [om.parse_case(a) for a in args[1:]]

lab = ls.VISSER4
fig, axs = plt.subplots(3, 1, figsize=(6.5, 7.6), sharex=True)
kw = dict(ls="none", marker="o", mfc="none", mec="C3", ms=5,
          label="Visser (1991), case 4")
axs[0].plot(lab[:, 0], lab[:, 1], **kw)
axs[1].plot(lab[:, 0], lab[:, 3], **kw)
axs[2].plot(lab[:, 0], lab[:, 2], **kw)

styles = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
cols = ["C0", "C1", "C2", "C4", "C5"]


def metrics(label, xh, H, mwl, xv, V):
    """Largest H and V with their distances from the shoreline, mean level
    at the shoreline, rms errors and biases at the measuring stations."""
    o, ov = np.argsort(xh), np.argsort(xv)
    out = {}
    for name, col, xs_, ys_ in (("H", 1, xh[o], H[o]), ("V", 2, xv[ov], V[ov]),
                                ("mwl", 3, xh[o], mwl[o])):
        ok = ~np.isnan(lab[:, col])
        err = np.interp(lab[ok, 0], xs_, ys_) - lab[ok, col]
        out[name] = (np.sqrt(np.mean(err ** 2)), err.mean())
    m = xh > 0.1
    print("%s: largest H %.4f m at %.2f m; mean level at the shoreline %.4f m;"
          " largest V %.3f m/s at %.2f m from the shoreline"
          % (label, H[m].max(), xh[m][np.argmax(H[m])],
             np.interp(0.0, xh[o], mwl[o]), V.max(), xv[np.argmax(V)]))
    for k, (r, b) in out.items():
        print("   %s: rms error %.4f, bias %+.4f" % (k, r, b))


for n, r in enumerate(runs):
    kind, case, label = r
    st = dict(color=cols[n % 5], ls=styles[n % 5], lw=1.4, label=label)
    if kind == "nhwave":
        xo, H, setup, U, V = ls.nhwave_profiles(case)
        m = xo > -0.5
        # NHWAVE's wave height needs zero crossings: none in the swash
        mh = m & (H > 0)
        axs[0].plot(xo[mh], H[mh], **st)
        axs[1].plot(xo[m], setup[m], **st)
        axs[2].plot(xo[m], V[m], **st)
        metrics(label, xo[mh], H[mh], setup[mh], xo[m], V[m])
        continue
    x, H, mwl = ls.wave_stats(case)
    xo = ls.V_XSHORE - x
    m = xo > -0.5
    axs[0].plot(xo[m], H[m], **st)
    axs[1].plot(xo[m], mwl[m], **st)
    xc, ux, uy, sy = ls.mean_current(case)
    xv = ls.V_XSHORE - xc
    mv = xv > -0.5
    axs[2].plot(xv[mv], uy[mv], **st)
    metrics(label, xo[m], H[m], mwl[m], xv[mv], uy[mv])
    print("   alongshore spread of the column means of V up to %.3f m/s;"
          " mean cross-shore velocity %.3f..%.3f m/s"
          % (sy[mv].max(), ux[mv].min(), ux[mv].max()))

for a in axs:
    a.axvline(ls.V_XBR, color="0.6", lw=0.8, ls=":")
    a.grid(alpha=0.3)
axs[2].text(ls.V_XBR + 0.05, 0.37, "breaker line (Visser)", fontsize=8,
            color="0.4")
axs[0].set_ylabel(r"$H$ [m]")
axs[1].set_ylabel(r"$\bar\eta$ [m]")
axs[2].set_ylabel(r"$\bar V$ [m/s]")
axs[2].set_xlabel("distance from the still-water shoreline [m]")
axs[0].set_ylim(0, None)
axs[2].set_xlim(-0.5, 6.0)
axs[0].legend(fontsize=8, loc="lower right", framealpha=0.9)
fig.tight_layout()
save(fig, out)

if snap:
    C = ls.cell_centres(snap, ls.latest_time(snap))
    eta = ls._read_field(os.path.join(snap, ls.latest_time(snap), "eta"))
    p = dict(l.split() for l in open(os.path.join(snap, "params.txt")))
    nz = int(p["nz"])
    # top cells of the columns: highest z per (x, y)
    key = np.round(C[:, 0], 6) + 1j * np.round(C[:, 1], 6)
    o = np.lexsort((C[:, 2], key))
    k = key[o]
    last = np.r_[k[1:] != k[:-1], True]
    xs, ys, es = C[o][last, 0], C[o][last, 1], eta[o][last]
    ux, uy = np.unique(np.round(xs, 6)), np.unique(np.round(ys, 6))
    E = np.full((len(uy), len(ux)), np.nan)
    ix = np.searchsorted(ux, np.round(xs, 6))
    iy = np.searchsorted(uy, np.round(ys, 6))
    E[iy, ix] = es
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    xo = ls.V_XSHORE - ux
    pc = ax.pcolormesh(xo, uy, E, shading="nearest", cmap="RdBu_r",
                       vmin=-0.06, vmax=0.06)
    ax.set_xlim(xo.max(), -0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("distance from the still-water shoreline [m]")
    ax.set_ylabel(r"$y$ [m]")
    fig.colorbar(pc, ax=ax, label=r"$\eta$ [m]", shrink=0.8)
    fig.tight_layout()
    base = os.path.splitext(out)[0]
    save(fig, base + "_snapshot.png")
