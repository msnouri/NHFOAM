#!/usr/bin/env python3
"""The flume tests of swash_flumes.py against the measurements and against
SWASH's results as published by Zijlema et al. (2011), both read from the
figures of that paper (digitize_zijlema2011.py).

usage: plot_swash_flumes.py bj <out.png> <dir> [<suffix>:<label> ...]
       plot_swash_flumes.py boers <out.png> <dir> [<suffix>:<label> ...]
    <dir> holds the runs bj13_<suffix>, bj15_<suffix> / b1b_<suffix>
    (default: L1, L1_wv, L2_wv, L2_hyb; runs that are not there are left
    out); also writes the pdf and prints the numbers
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import swash_flumes as sf

kind, out, root = sys.argv[1:4]
runs = [a.split(":", 1) for a in sys.argv[4:]] or [
    ["L1", "NHFOAM, 1 layer, front only"],
    ["L1_wv", "NHFOAM, 1 layer, whole wave + mixing length"],
    ["L2_wv", "NHFOAM, 2 layers, whole wave + mixing length"],
    ["L2_hyb", "NHFOAM, 2 layers, Kennedy et al."]]
COL = ["#1f5fa8", "#c0392b", "#2e8b57", "#8e44ad", "#e67e22"]
LS = ["-", "--", "-.", ":", (0, (5, 1))]


def errors(xl, yl, x, y):
    """rms and mean of y(x) - yl at the measured points that are wet in the run."""
    m = np.isfinite(y)
    d = np.interp(xl, x[m], y[m]) - yl
    ok = (xl >= x[m].min()) & (xl <= x[m].max())
    return np.sqrt(np.mean(d[ok] ** 2)), np.mean(d[ok]), d


if kind == "bj":
    fig, axs = plt.subplots(2, 2, figsize=(8.4, 5.2), sharex="col",
                            gridspec_kw=dict(height_ratios=[1.5, 1], hspace=0.08, wspace=0.22))
    for c, run in enumerate((13, 15)):
        a, b = axs[0, c], axs[1, c]
        lab = sf.load(f"bj78_run{run}.txt")
        xb, db = sf.bj_profile(run)
        a.plot(lab[:, 0], lab[:, 1], "kx", ms=6, mew=1.0, label="measured", zorder=5)
        a.plot(lab[:, 0], lab[:, 2], "o", mfc="none", mec="0.35", ms=5, mew=0.8,
               label="SWASH (Zijlema et al.)", zorder=4)
        r, m, _ = errors(lab[:, 0], lab[:, 1], lab[:, 0], lab[:, 2])
        print(f"run {run}: SWASH as published: rms {100 * r:.2f} cm, mean {100 * m:+.2f} cm "
              f"({len(lab)} gauges, Hm0 at the first {lab[0, 1]:.3f} m)")
        b.fill_between(xb, -db, -1.0, color="0.85", lw=0)
        b.plot(xb, -db, "k", lw=0.8)
        for k, (suf, label) in enumerate(runs):
            case = os.path.join(root, f"bj{run}_{suf}")
            if not os.path.exists(os.path.join(case, "ok")):
                continue
            s = sf.finish(case)
            a.plot(s["x"], s["hm0"], color=COL[k], ls=LS[k], lw=1.2, label=label)
            wet = np.isfinite(s["setup"]) & (s["depth"] > 0)
            b.plot(s["x"][wet], 10 * s["setup"][wet], color=COL[k], ls=LS[k], lw=1.2)
            r, m, d = errors(lab[:, 0], lab[:, 1], s["x"], s["hm0"])
            i1 = int(np.nanargmin(np.abs(s["x"] - lab[0, 0])))
            sel = s["x"] <= lab[-1, 0]
            print(f"   {label}: rms {100 * r:.2f} cm, mean {100 * m:+.2f} cm; at the gauges "
                  + " ".join(f"{100 * v:+.1f}" for v in d) + f" cm; Hm0 at the first gauge "
                  f"{s['hm0'][i1]:.3f} m; highest mean level {100 * np.nanmax(s['setup'][wet]):.2f} cm, "
                  f"lowest {100 * np.nanmin(s['setup'][wet]):.2f} cm")
        a.set_ylim(0.04, 0.25)
        a.set_xlim(0, xb[-2] + 0.5)
        a.grid(alpha=0.3)
        a.text(0.03, 0.08, f"({'ab'[c]}) run {run}", transform=a.transAxes)
        b.axhline(0, color="0.4", lw=0.5, ls=(0, (4, 3)))
        b.set_ylim(-0.8, 0.25)
        b.set_xlabel(r"$x$ [m]")
        b.grid(alpha=0.3)
        b.text(0.03, 0.82, f"({'cd'[c]})", transform=b.transAxes)
    axs[0, 0].set_ylabel(r"$H_{m0}$ [m]")
    axs[1, 0].set_ylabel(r"$z$ [m]; $10\,\overline{\eta}$ [m]")
    h, l = axs[0, 1].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc="upper center", ncol=3)
    fig.subplots_adjust(left=0.09, right=0.985, top=0.9, bottom=0.1)
    save(fig, out)

else:
    lab = sf.load("boers1b_hm0.txt")
    labT = sf.load("boers1b_tm02.txt")
    swH = sf.load("boers1b_hm0_swash.txt")
    swT = sf.load("boers1b_tm02_swash.txt")
    st = sf.load("boers_gauges.txt")[:, 0]
    xb, db = sf.boers_profile()
    # the measured spectra over the band of the figure: Hm0 and Tm02 there
    obs = []
    for n in range(1, 9):
        sp = sf.load(f"boers1b_spec_{n}.txt")
        f, S = sp[:, 0], sp[:, 1]
        m0 = np.sum(0.5 * (S[1:] + S[:-1]) * np.diff(f))
        m2 = np.sum(0.5 * (S[1:] * f[1:] ** 2 + S[:-1] * f[:-1] ** 2) * np.diff(f))
        obs.append((4 * np.sqrt(m0), np.sqrt(m0 / m2)))
    obs = np.array(obs)
    fp = 1.0 / sf.BOERS_1B["Tp"]

    def low(f, S):
        m = f <= 0.5 * fp
        return 4 * np.sqrt(np.sum(0.5 * (S[m][1:] + S[m][:-1]) * np.diff(f[m])))
    obs_low = np.array([low(*sf.load(f"boers1b_spec_{n}.txt").T) for n in range(1, 9)])
    fig = plt.figure(figsize=(8.4, 7.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.25, 0.7, 1.0], hspace=0.32, wspace=0.42)
    aH = fig.add_subplot(gs[0, :2]); aT = fig.add_subplot(gs[0, 2:])
    aB = fig.add_subplot(gs[1, :2]); aL = fig.add_subplot(gs[1, 2:])
    shown = (1, 4, 6, 8)
    aS = [fig.add_subplot(gs[2, k]) for k in range(4)]
    aH.plot(lab[:, 0], lab[:, 1], "o", mfc="none", mec="k", ms=4.5, mew=0.8, label="measured", zorder=5)
    aH.plot(swH[:, 0], swH[:, 1], color="0.45", lw=1.0, label="SWASH (Zijlema et al.)")
    aT.plot(labT[:, 0], labT[:, 1], "o", mfc="none", mec="k", ms=4.5, mew=0.8, zorder=5,
            label="measured")
    aT.plot(st, obs[:, 1], "s", mfc="k", mec="k", ms=4, zorder=6,
            label="measured spectra, 0--2 Hz")
    aT.plot(swT[:, 0], swT[:, 1], color="0.45", lw=1.0)
    r, m, _ = errors(lab[:, 0], lab[:, 1], swH[:, 0], swH[:, 1])
    print(f"Boers 1B: SWASH as published: Hm0 rms {100 * r:.2f} cm, mean {100 * m:+.2f} cm "
          f"({len(lab)} points)")
    print("   Hm0 of the measured spectra (0-2 Hz) / of the circles at the stations: "
          + " ".join(f"{a:.3f}/{b:.3f}" for a, b in zip(obs[:, 0], np.interp(st, lab[:, 0], lab[:, 1]))))
    aB.fill_between(xb, -db, -1.0, color="0.85", lw=0)
    aB.plot(xb, -db, "k", lw=0.8)
    aB.axhline(0, color="0.4", lw=0.5, ls=(0, (4, 3)))
    for n, xs in enumerate(st, 1):
        aB.plot([xs, xs], [0.0, 0.06], "k", lw=0.6)
        aB.text(xs, 0.08, str(n), ha="center", fontsize=7)
    for k, n in enumerate(shown):
        sp = sf.load(f"boers1b_spec_{n}.txt")
        aS[k].plot(sp[:, 0], 1e3 * sp[:, 1], "k", lw=1.3, label="measured")
    for k, (suf, label) in enumerate(runs):
        case = os.path.join(root, f"b1b_{suf}")
        if not os.path.exists(os.path.join(case, "ok")):
            continue
        s = sf.finish(case)
        aH.plot(s["x"], s["hm0"], color=COL[k], ls=LS[k], lw=1.2, label=label)
        aT.plot(s["x"], s["tm02"], color=COL[k], ls=LS[k], lw=1.2)
        wet = np.isfinite(s["setup"]) & (s["depth"] > 0)
        aL.plot(s["x"][wet], 100 * s["setup"][wet], color=COL[k], ls=LS[k], lw=1.2)
        aL.plot(s["x"][wet], 100 * s["hlow"][wet], color=COL[k], ls=LS[k], lw=0.7)
        # the spectra at the resolution of the measured ones (2 Hz / 86 points)
        nb = max(int(round(0.023 / (s["spec_f"][1] - s["spec_f"][0]))), 1)
        for kk, n in enumerate(shown):
            aS[kk].plot(s["spec_f"], 1e3 * np.convolve(s[f"spec_{n}"], np.ones(nb) / nb, "same"),
                        color=COL[k], ls=LS[k], lw=0.9)
        r, m, d = errors(lab[:, 0], lab[:, 1], s["x"], s["hm0"])
        ist = [int(np.nanargmin(np.abs(s["x"] - v))) for v in st]
        print(f"   {label}: Hm0 rms {100 * r:.2f} cm, mean {100 * m:+.2f} cm; at the stations "
              + " ".join(f"{100 * (s['hm0'][i] - o):+.1f}" for i, o in zip(ist, obs[:, 0]))
              + " cm; Tm02 (0-2 Hz) model/measured spectra "
              + " ".join(f"{s['tm02'][i]:.2f}/{o:.2f}" for i, o in zip(ist, obs[:, 1]))
              + f"; highest mean level {100 * np.nanmax(s['setup'][wet]):.2f} cm, lowest "
              f"{100 * np.nanmin(s['setup'][wet]):.2f} cm; Hm0 below fp/2 at the stations "
              + " ".join(f"{100 * s['hlow'][i]:.1f}" for i in ist) + " cm")
    aH.set_ylabel(r"$H_{m0}$ [m]"); aH.set_ylim(0.04, 0.25)
    aT.set_ylabel(r"$T_{m02}$ [s]"); aT.set_ylim(0.8, 2.0)
    aB.set_ylabel(r"$z$ [m]"); aB.set_ylim(-0.8, 0.2)
    aL.plot(st, 100 * obs_low, "s", mfc="k", mec="k", ms=4, zorder=6)
    aL.set_ylabel("[cm]")
    aL.text(0.03, 0.9, r"thin: $H_{m0}$ below $f_p/2$ (squares: measured, incoming)",
            transform=aL.transAxes, fontsize=8, va="center")
    aL.text(0.03, 0.2, r"thick: mean level", transform=aL.transAxes, fontsize=8, va="center")
    print("   measured Hm0 below fp/2 (shoreward-propagating part) at the stations: "
          + " ".join(f"{100 * v:.1f}" for v in obs_low) + " cm")
    for a, l in zip((aH, aT, aB, aL), "abcd"):
        a.set_xlim(0, 31); a.grid(alpha=0.3)
        a.text(0.975, 0.9 if l == "a" else 0.1, f"({l})", transform=a.transAxes, ha="right",
               va="center")
    aB.set_xlabel(r"$x$ [m]"); aL.set_xlabel(r"$x$ [m]"); aL.set_ylim(-1, 9.8)
    for k, n in enumerate(shown):
        aS[k].set_xlim(0, 2); aS[k].set_ylim(0, None); aS[k].grid(alpha=0.3)
        aS[k].set_xlabel(r"$f$ [Hz]")
        aS[k].set_title(f"({'efgh'[k]}) station {n}", fontsize=9, pad=3)
    aS[0].set_ylabel(r"$E$ [$10^{-3}$ m$^2$/Hz]")
    aH.legend(frameon=False, fontsize=8, loc="lower left")
    aT.legend(frameon=False, fontsize=8, loc="upper left")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.985, bottom=0.065)
    save(fig, out)
