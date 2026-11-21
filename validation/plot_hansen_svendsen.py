#!/usr/bin/env python3
"""Hansen & Svendsen (1979): wave height and mean water level across the
shoaling and surf zones, NHFOAM vs lab data.

usage: plot_hansen_svendsen.py <labFile> <T> <out.png> [kind:]<case>[:label] [...]
    kind = nhfoam (default) | nhwave | reef3d   (see other_models.py)
(also writes out.pdf)

Lab data: fengyanshi/BENCHMARK_FUNWAVE, car_hansen_svendsen/postprocessing/
ht_setup_061071.txt (spilling, T = 1.667 s) or ht_setup_031041.txt
(plunging, T = 3.333 s); columns x [m] from the toe of the 1:34.26 slope,
H [m], setup [m]. The model is averaged over the last 10 periods
(wave height: mean of max - min over each period).
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save
import other_models as om

use_latex(10)
labfile, T, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
runs = [om.parse_case(c) for c in sys.argv[4:]]
lab = np.loadtxt(labfile)
cols = ["#1f5fa8", "#e67e22", "#2e8b57", "#8e44ad", "#7f8c8d"]
nPer = 10


def wave_stats(t, e, T):
    """Mean wave height (max - min in each period of the last nPer periods)
    and mean water level."""
    t0 = t[-1] - nPer * T
    m = t > t0
    H = []
    for i in range(nPer):
        w = (t > t0 + i * T) & (t <= t0 + (i + 1) * T)
        H.append(e[w].max() - e[w].min())
    return np.mean(H), e[m].mean()


fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.0), sharex=True,
                         gridspec_kw=dict(height_ratios=[3, 2, 1.2]))
report = []
styles = ["-", "--", "-.", ":", (0, (1, 1))]
for (kind, case, label), col, ls in zip(runs, cols, styles):
    if kind == "nhfoam":
        t, eta = om.read_nhfoam(case)
        xg = np.array(sorted(eta))
    else:
        xg = np.loadtxt(os.path.join(case, "gauges_x.txt"))
        t, eta = om.read_any(kind, case, list(xg))
    st = np.array([wave_stats(t, eta[x], T) for x in xg])
    axes[0].plot(xg, 100 * st[:, 0], ls, color=col, lw=1.6, label=label)
    axes[1].plot(xg, 1000 * st[:, 1], ls, color=col, lw=1.6, label=label)
    Hm = np.interp(lab[:, 0], xg, st[:, 0])
    Sm = np.interp(lab[:, 0], xg, st[:, 1])
    report.append((label, 100 * np.sqrt(np.mean((Hm - lab[:, 1]) ** 2)),
                   100 * xg[np.nanargmax(st[:, 0])], 100 * np.nanmax(st[:, 0]),
                   1000 * np.sqrt(np.mean((Sm - lab[:, 2]) ** 2))))
axes[0].plot(lab[:, 0], 100 * lab[:, 1], "o", ms=4, mfc="none", color="#c0392b",
             label="Hansen \\& Svendsen (1979)")
axes[1].plot(lab[:, 0], 1000 * lab[:, 2], "o", ms=4, mfc="none", color="#c0392b")
axes[0].set_ylabel(r"$H$ [cm]")
axes[1].set_ylabel(r"$\overline{\eta}$ [mm]")
xb = np.linspace(-1, 12, 200)
hb = np.where(xb < 0, 0.36, 0.36 - xb / 34.26)
axes[2].fill_between(xb, -hb, -0.4, color="0.8")
axes[2].plot(xb, -hb, "k-", lw=1)
axes[2].axhline(0, color="#1f5fa8", lw=0.8)
axes[2].set_ylim(-0.4, 0.05)
axes[2].set_ylabel(r"$z$ [m]")
axes[2].set_xlabel(r"$x$ [m] from the toe of the 1:34.26 slope")
for ax in axes:
    ax.grid(alpha=0.3)
    ax.set_xlim(-1, 11.5)
axes[0].legend(loc="upper left")
fig.suptitle(rf"Hansen \& Svendsen (1979), $T = {T:.3g}$ s")
fig.tight_layout()
save(fig, out, dpi=130)
print("label: rms(H) [cm], x of max H [cm], max H [cm]; lab max H = "
      f"{100*lab[:,1].max():.2f} cm at x = {lab[np.argmax(lab[:,1]),0]:.2f} m")
for r in report:
    print(f"{r[0]}: rms {r[1]:.3f} cm, Hmax {r[3]:.2f} cm at x = {r[2]/100:.2f} m, "
          f"rms(setup) {r[4]:.2f} mm")
