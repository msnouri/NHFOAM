#!/usr/bin/env python3
"""Solitary waves on a conical island (Briggs et al. 1995; NTHMP benchmark
6): gauges and run-up of NHFOAM and NHWAVE against the lab data.

usage: plot_conical_island.py <gauges.png> <runup.png> <run> [<run> ...]
    run = <case>:<kind>:<path>:<dx>:<label>
    case = A | B | C, kind = nhfoam | nhwave, cases made by
    conical_island.nhfoam / conical_island.nhwave with the grid spacing dx
(also writes the PDFs, and prints peak and rms errors)

Lab records: the mean of each gauge before the wave (20-24 s) is removed;
each run is shifted in time so that its crest at the incident gauges meets
the lab crest (mean of gauges 1-4).  Model gauge g3 is the mirror image of
lab gauge 2 (and g4 of gauge 1).
"""
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import conical_island as ci

GAUGE_PANELS = [("g3", ["g2", "g3"], "gauges 2, 3"), ("g6", ["g6"], "gauge 6"),
                ("g9", ["g9"], "gauge 9"), ("g16", ["g16"], "gauge 16"),
                ("g22", ["g22"], "gauge 22")]
NHFOAM_STYLES = [("-", "#1f5fa8", 1.3), (":", "#4a90d9", 1.5), ("-.", "#7f8c8d", 1.1)]
NHWAVE_STYLE = ("--", "#e67e22", 1.3)


def styles(runs):
    """NHFOAM runs in blues (the first one solid), NHWAVE dashed orange."""
    out, n = {}, 0
    for case, kind, path, dx, label in runs:
        if label in out:
            continue
        if kind == "nhwave":
            out[label] = NHWAVE_STYLE
        else:
            out[label] = NHFOAM_STYLES[min(n, len(NHFOAM_STYLES) - 1)]
            n += 1
    return out


def crest(t, e):
    """Time and height of the maximum (parabolic fit)."""
    i = int(np.argmax(e))
    if 0 < i < len(e) - 1:
        a, b, c = e[i - 1], e[i], e[i + 1]
        den = a - 2 * b + c
        off = 0.5 * (a - c) / den if den != 0 else 0.0
        return t[i] + off * (t[i + 1] - t[i]), b - 0.25 * (a - c) * off
    return t[i], e[i]


def lab(case):
    t, g = ci.read_lab(case)
    pre = (t >= 20.0) & (t < 24.0)
    g = {k: v - v[pre].mean() for k, v in g.items()}
    tc = np.mean([crest(t, g[k])[0] for k in ("g1", "g2", "g3", "g4")])
    return t, g, tc


def model(case, kind, path, dx, tlab):
    if kind == "nhwave":
        t, g = ci.read_nhwave_gauges(path, case, dx)
        emax, _ = ci.nhwave_hmax(path, case, dx)
        flooded = ci.flooded_nhwave
    else:
        t, g = ci.read_nhfoam_gauges(path)
        emax, _ = ci.nhfoam_etamax(path, case, dx=dx)
        flooded = ci.flooded_nhfoam
    tm = np.mean([crest(t, g[k])[0] for k in ("g3", "g4")])
    return t + (tlab - tm), g, (emax, flooded)


def main():
    gpng, rpng = sys.argv[1], sys.argv[2]
    runs = []
    for arg in sys.argv[3:]:
        case, kind, path, dx, label = arg.split(":", 4)
        runs.append((case, kind, path, float(dx), label))
    cases = sorted(set(r[0] for r in runs))
    sty = styles(runs)

    figg, axg = plt.subplots(len(cases), len(GAUGE_PANELS),
                             figsize=(13.5, 2.5 * len(cases) + 0.6), squeeze=False)
    figr, axr = plt.subplots(1, len(cases), figsize=(4.6 * len(cases), 3.6),
                             squeeze=False)
    print(f"{'case':4s} {'run':24s} {'gauge':6s} {'peak model':>10s} {'peak lab':>9s} "
          f"{'err %':>6s} {'rms/H':>6s}")
    for r, case in enumerate(cases):
        H = ci.CASES[case]["Hd"] * ci.D
        tl, gl, tcl = lab(case)
        t0, t1 = tcl - 2.0, tcl + 12.0
        wl = (tl >= t0) & (tl <= t1)
        mods = [(kind, path, dx, label) + model(case, kind, path, dx, tcl)
                for c, kind, path, dx, label in runs if c == case]
        for c, (mg, lgs, title) in enumerate(GAUGE_PANELS):
            ax = axg[r, c]
            for k, lgn in enumerate(lgs):
                ax.plot(tl[wl], 100 * gl[lgn][wl], "k-" if k == 0 else "k:",
                        lw=1.0 if k == 0 else 0.8, label="lab" if k == 0 else None)
            tpL, epL = crest(tl[wl], gl[lgs[0]][wl])
            for kind, path, dx, label, tm, gm, _ in mods:
                st = sty[label]
                w = (tm >= t0) & (tm <= t1)
                ax.plot(tm[w], 100 * gm[mg][w], st[0], color=st[1], lw=st[2], label=label)
                tp, ep = crest(tm[w], gm[mg][w])
                em = np.interp(tl[wl], tm, gm[mg])
                rms = np.sqrt(np.mean((em - gl[lgs[0]][wl]) ** 2))
                print(f"{case:4s} {label:24s} {mg:6s} {100*ep:9.2f}  {100*epL:8.2f} "
                      f"{100*(ep-epL)/epL:6.1f} {rms/H:6.3f}")
            ax.set_xlim(t0, t1)
            ax.grid(alpha=0.3)
            if r == 0:
                ax.set_title(title)
            if c == 0:
                ax.set_ylabel(rf"case {case}, $H/d = {ci.CASES[case]['Hd']:g}$" + "\n" + r"$\eta$ [cm]")
            if r == len(cases) - 1:
                ax.set_xlabel(r"$t$ [s]")
        # run-up
        ax = axr[0, r]
        thL, RL = ci.read_lab_runup(case)
        ax.plot(thL, RL / ci.D, "ko", ms=4, label="lab")
        ang = np.arange(0.0, 360.0 + 1e-9, 2.5)
        for kind, path, dx, label, tm, gm, (emax, flooded) in mods:
            st = sty[label]
            Rm = ci.runup(emax, case, flooded, ang, dx=dx)
            ax.plot(ang, Rm / ci.D, st[0], color=st[1], lw=st[2], label=label)
            RmL = ci.runup(emax, case, flooded, thL, dx=dx)
            ok = np.isfinite(RmL)
            rms = np.sqrt(np.mean((RmL[ok] - RL[ok]) ** 2))
            print(f"{case:4s} {label:24s} run-up: mean model/lab {np.mean(RmL[ok]/RL[ok]):.3f}, "
                  f"rms {100*rms:.2f} cm (R/d {rms/ci.D:.3f}), "
                  f"max {100*np.nanmax(RmL):.2f} vs {100*RL.max():.2f} cm")
        ax.set_xlim(0, 360)
        ax.set_xticks(np.arange(0, 361, 90))
        ax.set_xlabel(r"angle around the island [$^\circ$]")
        if r == 0:
            ax.set_ylabel(r"run-up $R/d$")
        ax.set_title(rf"case {case}, $H/d = {ci.CASES[case]['Hd']:g}$")
        ax.grid(alpha=0.3)
    axg[0, 0].legend(loc="upper right", fontsize=9, frameon=False)
    axr[0, 0].legend(loc="lower right", fontsize=9, frameon=False)
    figg.tight_layout()
    figr.tight_layout()
    save(figg, gpng)
    save(figr, rpng)


if __name__ == "__main__":
    main()
