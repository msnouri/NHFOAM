#!/usr/bin/env python3
"""SWASH's tutorials 1-3 in NHFOAM (swash_tutorials.py): significant wave
height and set-up along the beach, against the linear shoaling of the
spectrum (no measurements exist for these cases).

usage: plot_swash_tutorials.py <out.png> <run> [<run> ...]
    run = <path>:<label>   (also writes out.pdf and prints the numbers)
"""
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import swash_tutorials as st

COL = ["#1f5fa8", "#c0392b", "#e67e22", "#7f8c8d", "#2e8b57"]
LS = ["-", "--", "-.", ":", (0, (5, 1.5, 1, 1.5))]

runs = [tuple(a.split(":", 1)) for a in sys.argv[2:]]
fig, axs = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
seen = set()
for k, (path, label) in enumerate(runs):
    pr = st.params(path)
    x, dep, hs, su = st.nhfoam_stats(path)
    axs[0].plot(x, hs, ls=LS[k % 5], color=COL[k % 5], lw=1.3, label=label)
    wet = dep > 3 * pr["film"]
    axs[1].plot(x[wet], 100 * su[wet], ls=LS[k % 5], color=COL[k % 5], lw=1.3)
    key = (pr["Hs"], pr["Tp"])
    if key not in seen:
        seen.add(key)
        xs = np.arange(0.0, 760.0, 10.0)
        axs[0].plot(xs, st.linear_shoaling(xs, pr["Hs"], pr["Tp"]), color="k",
                    lw=0.8, ls=(0, (1, 1)),
                    label="linear shoaling" if len(seen) == 1 else None)
    xr = [2.5, 100.0, 300.0, 500.0, 650.0]
    lin = st.linear_shoaling(xr, pr["Hs"], pr["Tp"])
    b = st.band_heights(path, xr)
    i0 = int(np.nanargmax(hs))
    print(f"{label}: largest Hs {hs[i0]:.3f} m at x = {x[i0]:.0f} m (depth {dep[i0]:.2f} m); "
          f"set-up at the still-water shoreline "
          f"{100 * su[np.argmin(np.abs(dep))]:.1f} cm; lowest mean level "
          f"{100 * np.nanmin(su[wet]):.2f} cm")
    for j, xx in enumerate(xr):
        print(f"    x = {xx:5.1f} m: Hs {b[j, 0]:.3f} m (linear {lin[j]:.3f}, "
              f"{100 * (b[j, 0] / lin[j] - 1):+.1f}%); 0.5-1.5 fp {b[j, 1]:.3f}, "
              f"above 1.5 fp {b[j, 2]:.3f}")
axs[0].set_ylabel(r"$H_s$ [m]")
axs[1].set_ylabel("mean water level [cm]")
axs[1].set_xlabel(r"$x$ [m]")
axs[0].legend(frameon=False, fontsize=9, loc="lower left")
for a in axs:
    a.grid(alpha=0.3)
    a.axvline(1000.0 * st.D0 / (st.D0 + st.ZTOP), color="0.6", lw=0.6)
axs[0].set_xlim(0, 900)
fig.tight_layout()
save(fig, sys.argv[1])
