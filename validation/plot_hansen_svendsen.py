#!/usr/bin/env python3
"""Hansen & Svendsen (1979): wave height and mean water level across the
shoaling and surf zones, NHFoam vs lab data.

usage: plot_hansen_svendsen.py <labFile> <T> <out.png> <case>[:label] [...]
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

use_latex(10)
labfile, T, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
cases = [(c.split(":")[0], c.split(":")[1] if ":" in c else os.path.basename(c))
         for c in sys.argv[4:]]
lab = np.loadtxt(labfile)
cols = ["#1f5fa8", "#e67e22", "#2e8b57"]
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
for (case, label), col in zip(cases, cols):
    fn = os.path.join(case, "postProcessing/gauges/0/eta.dat")
    d = np.loadtxt(fn, comments="#")
    hdr = [l.split() for l in open(fn) if l.startswith("# G")]
    xg = np.array([float(h[2].strip("(")) for h in hdr])
    st = np.array([wave_stats(d[:, 0], d[:, i + 1], T) for i in range(len(xg))])
    axes[0].plot(xg, 100 * st[:, 0], "-", color=col, lw=1.6, label=label)
    axes[1].plot(xg, 1000 * st[:, 1], "-", color=col, lw=1.6, label=label)
    Hm = np.interp(lab[:, 0], xg, st[:, 0])
    report.append((label, 100 * np.sqrt(np.mean((Hm - lab[:, 1]) ** 2)),
                   100 * xg[np.nanargmax(st[:, 0])], 100 * np.nanmax(st[:, 0])))
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
    print(f"{r[0]}: rms {r[1]:.3f} cm, Hmax {r[3]:.2f} cm at x = {r[2]/100:.2f} m")
