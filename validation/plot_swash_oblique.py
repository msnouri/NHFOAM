#!/usr/bin/env python3
"""SWASH's tutorial 4 in NHFOAM (swash_tutorials.nhfoam_oblique):
short-crested waves at an angle to a plane beach, periodic alongshore.
Significant wave height, mean water level, undertow and longshore current
across the beach, and the free surface at the end of the run.  No
measurements exist for this case; the references are linear theory without
dissipation for the wave components that the wavemaker makes (wave height)
and the return of the mass flux of linear waves with the model's wave
heights (cross-shore mean velocity).

usage: plot_swash_oblique.py <out.png> <case>
    (also writes out.pdf and prints the numbers)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import swash_tutorials as st

out, case = sys.argv[1:3]
pr = st.params(case)
comps = np.loadtxt(os.path.join(case, "components.txt")).T
free = st.swash_components(pr["Hs"], pr["Tp"], pr["direction"], pr["dspr"], fit=False)
xs0 = L4 = st.L4 * st.D4 / (st.D4 + st.ZTOP4)          # still-water shoreline
g = st.oblique_stats(case)
xm, em, um, vm, vsd = st.oblique_means(case)
xl = np.arange(0.5, 900.0, 5.0)
hl, ul, _, sxy = st.linear_oblique(xl, comps)
hl0, _, _, sxy0 = st.linear_oblique(xl, free)

fig = plt.figure(figsize=(8.4, 10.0))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 0.8, 1, 0.75], hspace=0.32)
a0, a1, a2 = (fig.add_subplot(gs[i]) for i in range(3))
a3 = fig.add_subplot(gs[3])

x, dep, hs, hss, hig, su = g[0]
a0.plot(x, hs, color="#1f5fa8", lw=1.4, label=r"NHFOAM, $y = 102$ m")
a0.plot(g[1][0], g[1][2], color="#1f5fa8", lw=1.0, ls="--", label=r"NHFOAM, $y = 50$ m")
a0.plot(x, hss, color="#e67e22", lw=1.0, label=r"$y = 102$ m, above $0.5 f_p$")
a0.plot(x, hig, color="#2e8b57", lw=1.0, label=r"$y = 102$ m, below $0.5 f_p$")
a0.plot(xl[xl < 830], hl[xl < 830], color="k", lw=0.8, ls=(0, (1, 1)),
        label="linear, no dissipation")
a0.set_ylabel(r"$H_s$ [m]")
a0.set_ylim(0, 3.2)
a0.legend(frameon=False, fontsize=9, ncol=2, loc="upper left")

wet = dep > 3 * pr["film"]
a1.plot(x[wet], 100 * su[wet], color="#1f5fa8", lw=1.4, label=r"gauges, $y = 102$ m")
mw = st.depth4(xm) > 3 * pr["film"]
a1.plot(xm[mw], 100 * em[mw], color="#c0392b", lw=1.0, ls="--", label="alongshore mean")
a1.set_ylabel("mean water level [cm]")
a1.legend(frameon=False, fontsize=9, loc="upper left")

wm = st.depth4(xm) > 0.0
a2.fill_between(xm[wm], (vm - vsd)[wm], (vm + vsd)[wm], color="#c0392b", alpha=0.2, lw=0)
a2.plot(xm[wm], vm[wm], color="#c0392b", lw=1.4, label=r"alongshore, $\overline{V}$")
a2.plot(xm[wm], um[wm], color="#1f5fa8", lw=1.4, label=r"cross-shore, $\overline{U}$")
# the current over the last interval between two written times
times = sorted((d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0),
               key=float)
vlast = None
if len(times) > 1:
    t1, t2 = float(times[-2]), float(times[-1])
    v1 = st.oblique_means(case, times[-2])[3]
    vlast = ((t2 - pr["tStat"]) * vm - (t1 - pr["tStat"]) * v1) / (t2 - t1)
    a2.plot(xm[wm], vlast[wm], color="#c0392b", lw=0.9, ls="--",
            label=rf"$\overline{{V}}$, last {t2 - t1:g} s")
# return of the mass flux of linear waves that have the model's Hs
sel = xl < 895
ur = ul * (np.interp(xl, x, hs) / hl) ** 2
a2.plot(xl[sel], ur[sel], color="k", lw=0.8, ls=(0, (1, 1)), label="return of the wave mass flux")
a2.axhline(0, color="0.5", lw=0.5)
a2.set_ylabel("mean velocity [m/s]")
a2.set_xlabel(r"$x$ [m]")
a2.legend(frameon=False, fontsize=9, loc="upper left")
for a in (a0, a1, a2):
    a.grid(alpha=0.3)
    a.axvline(xs0, color="0.6", lw=0.6)
    a.set_xlim(0, 960)

sdir = os.path.join(case, "postProcessing", "surface")
snaps = [f[:-4] for f in os.listdir(sdir) if f.endswith(".dat")] if os.path.isdir(sdir) else []
a3.set_xlim(0, 960); a3.set_ylim(0, st.W4)
a3.set_xlabel(r"$x$ [m]"); a3.set_ylabel(r"$y$ [m]")
if snaps:
    tlast = max(snaps, key=float)
    S = np.loadtxt(os.path.join(sdir, tlast + ".dat"))
    X = np.unique(np.round(S[:, 0], 6)); Y = np.unique(np.round(S[:, 1], 6))
    E = np.full((len(Y), len(X)), np.nan)
    ix = np.searchsorted(X, np.round(S[:, 0], 6)); iy = np.searchsorted(Y, np.round(S[:, 1], 6))
    E[iy, ix] = np.where(S[:, 3] - S[:, 2] > 2 * pr["film"], S[:, 3], np.nan)
    pc = a3.pcolormesh(X, Y, E, cmap="RdBu_r", vmin=-1.5, vmax=1.5, shading="nearest",
                       rasterized=True)
    cb = fig.colorbar(pc, ax=a3, orientation="vertical", pad=0.015, fraction=0.03)
    a3.text(0.012, 0.9, rf"$\eta$ [m], $t = {float(tlast):g}$ s", transform=a3.transAxes,
            ha="left", va="center", bbox=dict(fc="w", ec="none", alpha=0.8, pad=1.5))
for a, l in zip((a0, a1, a2, a3), "abcd"):
    a.text(0.985, 0.9, f"({l})", transform=a.transAxes, ha="right", va="center",
           bbox=dict(fc="w", ec="none", alpha=0.8, pad=1.5) if a is a3 else None)
fig.subplots_adjust(left=0.1, right=0.97, top=0.985, bottom=0.055)
# the colour bar takes width from panel (d) only: give the others the same
p3 = a3.get_position()
for a in (a0, a1, a2):
    p = a.get_position()
    a.set_position([p.x0, p.y0, p3.x1 - p.x0, p.height])
save(fig, out)

f, a, th, _ = comps
w = a * a
md = lambda th, w: np.degrees(np.arctan2((w * np.sin(np.radians(th))).sum(),
                                         (w * np.cos(np.radians(th))).sum()))
print(f"components: {len(f)}; shore-normal share of the energy "
      f"{w[np.abs(th) < 0.01].sum() / w.sum():.2f}; mean direction {md(th, w):.1f} deg "
      f"(before the fit {md(free[2], free[1] ** 2):.1f}); Sxy/rho {sxy:.3f} m3/s2 "
      f"(before the fit {sxy0:.3f})")
for X0 in (2.5, 200, 400, 600, 700, 760, 800, 840, 860, 880, 900, 910, 920):
    i = int(np.argmin(np.abs(x - X0)))
    j = int(np.argmin(np.abs(xl - x[i])))
    k = int(np.argmin(np.abs(g[1][0] - X0)))
    print(f"x = {x[i]:6.1f} m, h = {dep[i]:5.2f} m: Hs {hs[i]:.3f} m (y = 50 m: {g[1][2][k]:.3f}; "
          f"linear {hl[j]:.3f}, {100 * (hs[i] / hl[j] - 1):+.1f}%), above 0.5 fp {hss[i]:.3f}, "
          f"below {hig[i]:.3f}, Hs/(h + setup) {hs[i] / max(dep[i] + su[i], 1e-9):.2f}, "
          f"mean level {100 * su[i]:+.1f} cm")
i = int(np.nanargmax(hs))
print(f"largest Hs {hs[i]:.3f} m at x = {x[i]:.0f} m (h = {dep[i]:.2f} m)")
iv = int(np.argmax(np.where(wm, vm, -9)))
print(f"largest longshore current {vm[iv]:.3f} m/s at x = {xm[iv]:.1f} m (h = {st.depth4(xm[iv]):.2f} m), "
      f"alongshore std there {vsd[iv]:.3f} m/s")
if vlast is not None:
    il = int(np.argmax(np.where(wm, vlast, -9)))
    print(f"   over the last interval: {vlast[il]:.3f} m/s at x = {xm[il]:.1f} m")
half = xm[wm & (vm > 0.5 * vm[iv])]
print(f"   above half its maximum from x = {half.min():.0f} to {half.max():.0f} m")
iu = int(np.argmin(np.where(wm, um, 9)))
print(f"strongest offshore mean flow {um[iu]:.3f} m/s at x = {xm[iu]:.1f} m (h = {st.depth4(xm[iu]):.2f} m)")
for X0 in (100, 300, 500, 700, 800, 850, 880, 900, 910):
    i = int(np.argmin(np.abs(xm - X0))); j = int(np.argmin(np.abs(xl - X0)))
    last = f", last interval {vlast[i]:+.4f}" if vlast is not None else ""
    print(f"x = {xm[i]:6.1f} m: U {um[i]:+.4f} (flux return {ur[j]:+.4f}; with the target "
          f"heights {ul[j]:+.4f}), V {vm[i]:+.4f} (+-{vsd[i]:.4f}{last}) m/s")
ish = int(np.argmin(np.abs(xm - xs0)))
print(f"mean level at the still-water shoreline {100 * em[ish]:.1f} cm; lowest {100 * em[mw].min():.1f} cm "
      f"at x = {xm[mw][np.argmin(em[mw])]:.0f} m")
