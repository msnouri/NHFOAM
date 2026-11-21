#!/usr/bin/env python3
"""NTHMP landslide benchmark 2 (Enet & Grilli 2007, d = 61 mm): NHFOAM gauges
vs lab data.

usage: plot_enet_grilli.py <labFile> <out.png> [kind:]<case>[:label] ...
    kind = nhfoam (default) | nhwave   (see other_models.py)
(also writes out.pdf)

Lab data: JimKirby/NHWAVE, examples/landslide_benchmarks/BM2_enet_grilli/d61
(columns t [s], g2, g3, g4 [mm]; the slide is released at t = 0).
No time shift is applied.
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save
import other_models as om

use_latex(10)

labfile, out = sys.argv[1:3]
runs = [om.parse_case(c) for c in sys.argv[3:]]
lab = np.loadtxt(labfile)
b, ut, a0 = 0.395, 1.70, 1.12
t0 = ut / a0
cols = ["#1f5fa8", "#e67e22", "#2e8b57", "#8e44ad", "#7f8c8d"]

# 5-point running mean of the lab records (as in the NHWAVE benchmark script)
def smooth(y, n=5):
    return np.convolve(y, np.ones(n) / n, mode="same")

panels = [("g2", 1, r"$g_2$: $(x, y) = (1.469, 0.35)$ m"),
          ("g3", 2, r"$g_3$: $(x, y) = (1.929, 0)$ m"),
          ("g4", 3, r"$g_4$: $(x, y) = (1.919, 0.5)$ m")]
fig, axes = plt.subplots(len(panels) + 1, 1, figsize=(7.5, 9.5), sharex=True)
report = []
styles = ["-", "--", "-.", ":", (0, (1, 1))]


def gauges(kind, case):
    """Time and {g1..g4: eta} for an NHFOAM or NHWAVE run."""
    if kind == "nhfoam":
        fn = os.path.join(case, "postProcessing/gauges/0/eta.dat")
        d = np.loadtxt(fn, comments="#")
        names = [l.split()[1] for l in open(fn) if l.startswith("# g")]
        return d[:, 0], {n: d[:, i + 1] for i, n in enumerate(names)}
    t, e = om.read_nhwave(case, ["g1", "g2", "g3", "g4"])
    return t, e


for (kind, case, label), col, ls in zip(runs, cols, styles):
    t, eg = gauges(kind, case)
    axes[0].plot(t / t0, eg["g1"] / b, ls, color=col, lw=1.5, label=label)
    row = [label]
    for ax, (g, k, title) in zip(axes[1:], panels):
        e = eg[g]
        ax.plot(t / t0, e / b, ls, color=col, lw=1.5, label=label)
        m = (lab[:, 0] > 0) & (lab[:, 0] <= t[-1])
        el = np.interp(lab[m, 0], t, e)
        yl = smooth(lab[:, k] / 1000)[m]
        # extremes of the first wave train (t < 2.5 s)
        w = lab[m, 0] < 2.5
        row += [1000 * yl[w].min(), 1000 * el[w].min(), 1000 * yl[w].max(), 1000 * el[w].max()]
    report.append(row)
axes[0].set_title(r"$g_1$: $(x, y) = (x_0, 0)$, above the initial slide position (no lab record)")
for ax, (g, k, title) in zip(axes[1:], panels):
    ax.plot(lab[:, 0] / t0, smooth(lab[:, k] / 1000) / b, "--", color="k", lw=1.2,
            label="Enet \\& Grilli (2007)")
    ax.set_title(title)
for ax in axes:
    ax.set_ylabel(r"$\eta / b$")
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 2.5)
axes[1].legend(loc="lower right")
axes[0].legend(loc="lower right")
axes[-1].set_xlabel(r"$t / t_0$, $t_0 = u_t / a_0 = 1.52$ s")
fig.suptitle(r"3D rigid slide, Enet \& Grilli (2007), $d = 61$ mm")
fig.tight_layout()
save(fig, out, dpi=130)
print("first wave train (t < 2.5 s), mm:  lab min / model min / lab max / model max")
for r in report:
    print(r[0])
    for (g, k, title), i in zip(panels, range(1, len(r), 4)):
        print(f"  {g}: " + " ".join(f"{v:8.2f}" for v in r[i:i + 4]))
