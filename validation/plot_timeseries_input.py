#!/usr/bin/env python3
"""Time-series wave input on a flat tank (timeseries_input.py): a solitary
wave made at x = 0 from its record, against the exact wave, at x = 12 m.

usage: plot_timeseries_input.py <out.png> <dir>
    dir holds the runs ts_<mode>_<10|30|60> (modes exact, eta, etaC, paddle)
(also writes the PDF and prints the table of timeseries_input.analyse)
"""
import os, sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import timeseries_input as ti
import other_models as om

RUNS = [("exact", "exact initial condition", ("-", "#e67e22", 1.1)),
        ("eta", r"record, $c = \sqrt{g(h+\eta)}$", ("-", "#1f5fa8", 1.3)),
        ("etaC", r"record, $c = c_{\mathrm{exact}}$", (":", "#4a90d9", 1.6)),
        ("paddle", "paddle record", ("-.", "#7f8c8d", 1.1))]
X = 12.0


def main():
    out, root = sys.argv[1], sys.argv[2]
    fig, ax = plt.subplots(2, 3, figsize=(14, 5.6), sharex=True,
                           gridspec_kw=dict(height_ratios=[1.6, 1]))
    for c, Hh in enumerate((0.1, 0.3, 0.6)):
        sw = ti.wave(Hh)
        H = Hh * ti.H0
        tl = None
        for mode, label, st in RUNS:
            case = os.path.join(root, f"ts_{mode}_{int(round(100 * Hh)):02d}")
            p = dict(l.split() for l in open(os.path.join(case, "params.txt")))
            t, g = om.read_nhfoam(case)
            ex = ti.expected(sw, mode, X, t, float(p["t0"]))
            tp = t[np.argmax(ex)]
            w = np.abs(t - tp) < 2.5
            if tl is None:
                tt = np.linspace(-2.5, 2.5, 501)
                ax[0, c].plot(tt, sw.eta(-sw.c * tt) / H, "k-", lw=0.8, label="exact wave")
                tl = True
            ax[0, c].plot(t[w] - tp, g[X][w] / H, st[0], color=st[1], lw=st[2], label=label)
            ax[1, c].plot(t[w] - tp, 100 * (g[X][w] - ex[w]) / H, st[0], color=st[1], lw=st[2])
        ax[0, c].set_title(rf"$H/h = {Hh:g}$, $x = {X:g}$ m from the inflow")
        ax[1, c].set_xlabel(r"$t - t_{\mathrm{crest}}$ [s]")
        for a in ax[:, c]:
            a.grid(alpha=0.3)
            a.set_xlim(-2.5, 2.5)
    ax[0, 0].set_ylabel(r"$\eta/H$")
    ax[1, 0].set_ylabel(r"$(\eta - \eta_{\mathrm{exact}})/H$ [\%]")
    ax[0, 0].legend(loc="upper left", fontsize=9, frameon=False)
    fig.tight_layout()
    save(fig, out)
    cases = [os.path.join(root, f"ts_{m}_{h}") for h in ("10", "30", "60") for m, _, _ in RUNS]
    for r in ti.analyse(cases):
        print("%-14s %-6s H/h=%.2f x=%5.1f  peak/H=%.4f  dt_peak=%+.4f s  rms/H=%.4f  trail/H=%.4f"
              % ((os.path.basename(r[0]),) + r[1:]))


if __name__ == "__main__":
    main()
