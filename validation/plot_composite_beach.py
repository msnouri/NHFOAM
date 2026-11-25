#!/usr/bin/env python3
"""Solitary waves on a composite beach with a vertical wall (NTHMP
benchmark 5): gauges of NHFOAM (three ways of making the wave) and NHWAVE
against the lab, and the run-up on the wall.

usage: plot_composite_beach.py <gauges.png> <dir>
    dir holds the runs of composite_beach.py: on the 1 cm grid
    nf_<A|B|C>_<gauge4|gauge4c|paddle|solitary> (gauge4c: celerity="pulse")
    and nw_A_01, on the 2 cm grid nf_<A|B|C>_solitary_02 and nw_<A|B|C>
(also writes the PDF and prints the crest errors, rms differences and the
run-up on the wall)

Times are those of the lab records; no run is shifted.
"""
import os, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import composite_beach as cb

RUNS = [("nf_{}_gauge4", "NHFOAM, gauge 4 record", ("-", "#1f5fa8", 1.3)),
        ("nf_{}_paddle", "NHFOAM, paddle record", (":", "#4a90d9", 1.6)),
        ("nf_{}_solitary_02", "NHFOAM, solitary wave, 2 cm", ("-.", "#7f8c8d", 1.1)),
        ("nw_{}_inv", "NHWAVE, solitary wave, 2 cm", ("--", "#e67e22", 1.3)),
        # in the table only: NHWAVE with its k-epsilon model, the gauge 4
        # record with c = sqrt(g (d + H)), the solitary wave on the 1 cm
        # grid, NHWAVE on the 1 cm grid
        ("nw_{}", "NHWAVE k-eps, 2 cm", None),
        ("nf_{}_gauge4c", "NHFOAM, gauge 4, c of pulse", None),
        ("nf_{}_solitary", "NHFOAM, solitary, 1 cm", None),
        ("nw_{}_01", "NHWAVE k-eps, 1 cm", None)]
PANELS = [("g4", "gauge 4"), ("g5", "gauge 5"), ("g7", "gauge 7"),
          ("g9", "gauge 9"), ("g10", "gauge 10")]


def main():
    out, root = sys.argv[1], sys.argv[2]
    names = ["A", "B", "C"]
    fig, ax = plt.subplots(len(names), len(PANELS), figsize=(13.5, 2.5 * len(names) + 0.9),
                           squeeze=False)
    Rlab = cb.read_runup()
    print(f"{'case':4s} {'run':24s} {'rms [cm]':>8s} {'wall R [cm]':>11s} {'lab':>6s}  peaks (model/lab - 1)")
    for r, name in enumerate(names):
        tl, gl = cb.read_lab(name)
        t0 = tl[np.argmax(gl["g4"])] - 2.0
        t1 = {"A": cb.T_END, "B": 290.0, "C": 287.0}[name]
        wl = (tl >= t0) & (tl <= t1)
        H4 = gl["g4"].max()
        for c, (g, title) in enumerate(PANELS):
            ax[r, c].plot(tl[wl], 100 * gl[g][wl], "k-", lw=1.0, label="lab")
        for pat, label, st in RUNS:
            case = os.path.join(root, pat.format(name))
            try:
                tm, gm = cb.read_model(case)
            except (OSError, ValueError):
                print(f"{name:4s} {label:24s} no results in {case}")
                continue
            w = (tm >= t0) & (tm <= t1)
            for c, (g, title) in enumerate(PANELS):
                if st:
                    ax[r, c].plot(tm[w], 100 * gm[g][w], st[0], color=st[1], lw=st[2],
                                  label=label)
            try:
                pk = cb.crests(case, name)
                print(f"{name:4s} {label:24s} {100 * cb.rms(case, name):8.2f} "
                      f"{100 * cb.wall_runup(case):11.2f} {100 * Rlab[name]:6.2f}  "
                      + " ".join(f"{k}:{100 * v[0]:+.0f}%" for k, v in pk.items()))
            except (OSError, ValueError, IndexError):
                print(f"{name:4s} {label:24s} unfinished")
        for c, (g, title) in enumerate(PANELS):
            a = ax[r, c]
            a.set_xlim(t0, t1)
            a.grid(alpha=0.3)
            if r == 0:
                a.set_title(title)
            if c == 0:
                a.set_ylabel(rf"case {name}, $H/d = {H4 / cb.D:.3g}$" + "\n" + r"$\eta$ [cm]")
            if r == len(names) - 1:
                a.set_xlabel(r"$t$ [s]")
    handles, labels = ax[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), fontsize=10,
               frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, out)


if __name__ == "__main__":
    main()
