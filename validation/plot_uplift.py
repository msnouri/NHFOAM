#!/usr/bin/env python3
"""Bottom uplift (Hammack 1973 type): NHFOAM gauges vs linear potential theory.

Bottom: zeta(x, t) = zeta0 (1 - exp(-alpha t)) for |x| < b, alpha = 1.11/tc.
Linear theory (constant depth h):
    eta(x, t) = (1/pi) int_0^inf etah(k, t) cos(k x) dk
    etah = 2 zeta0 sin(k b)/(k cosh(k h)) * alpha/(alpha^2 + omega^2)
           * [omega sin(omega t) + alpha cos(omega t) - alpha exp(-alpha t)]

usage: plot_uplift.py <case> <out.png>   (also writes out.pdf)
"""
import sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex(10)
case, out = sys.argv[1], sys.argv[2]
p = dict(l.split() for l in open(f"{case}/params.txt"))
h, b, z0, tc, dx = (float(p[k]) for k in ("h", "b", "zeta0", "tc", "dx"))
g = 9.81
alpha = 1.11 / tc

fn = f"{case}/postProcessing/gauges/0/eta.dat"
d = np.loadtxt(fn, comments="#")
t = d[:, 0]
gx = [float(l.split("(")[1].split()[0]) for l in open(fn) if l.startswith("# G")]

# Linear theory at the column centres the gauges sample
k = np.linspace(1e-6, 40.0 / h, 200001)
dk = k[1] - k[0]
w = np.sqrt(g * k * np.tanh(k * h))
amp = 2 * z0 * np.sin(k * b) / (k * np.cosh(k * h)) * alpha / (alpha**2 + w**2)


def eta_lin(x, tt):
    res = np.empty_like(tt)
    ck = np.cos(k * x) * amp
    for i, ti in enumerate(tt):
        f = ck * (w * np.sin(w * ti) + alpha * np.cos(w * ti) - alpha * math.exp(-alpha * ti))
        res[i] = (f.sum() - 0.5 * (f[0] + f[-1])) * dk / np.pi
    return res


tt = t[::4]
# gauges up to x/b = 20 (further out the record is cut by the end of the run)
sel = [i for i, x in enumerate(gx) if x / b <= 20.5]
fig, axes = plt.subplots(len(sel), 1, figsize=(7.0, 1.55 * len(sel)), sharex=True)
for i, ax in zip(sel, axes):
    x = gx[i]
    xc = (math.floor(x / dx) + 0.5) * dx
    el = eta_lin(xc, tt)
    ax.plot(t, d[:, i + 1] / z0, color="#1f5fa8", lw=1.4, label="NHFOAM, 5 layers")
    ax.plot(tt, el / z0, "k--", lw=1.0, label="linear potential theory")
    ax.set_ylabel(r"$\eta/\zeta_0$")
    ax.text(0.99, 0.82, rf"$x/b = {x / b:.0f}$", transform=ax.transAxes, ha="right")
    ax.grid(alpha=0.3)
    err = np.max(np.abs(np.interp(tt, t, d[:, i + 1]) - el)) / z0
    print(f"x/b = {x/b:4.1f}: max |eta - eta_lin| / zeta0 = {err:.4f}, peak NHFOAM {d[:, i+1].max()/z0:.4f}, linear {el.max()/z0:.4f}")
axes[0].legend(loc="lower right", ncol=2)
axes[-1].set_xlabel(r"$t$ [s]")
fig.suptitle(rf"Bottom uplift: $h = {h:g}$ m, $b = {b:g}$ m, $\zeta_0/h = {z0 / h:g}$, $t_c = {tc:g}$ s")
fig.tight_layout()
save(fig, out)
