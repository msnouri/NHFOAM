#!/usr/bin/env python3
"""Ting & Kirby (1994) spilling breakers: time-mean horizontal velocity
(undertow) at the measuring stations.

usage: plot_undertow.py <fig5.dat> <out.png> [kind:]<case>[:label] ...
    fig5.dat: panel, x_TK [m], (x - xb)/hb, h/hb, (z - mean level)/h,
              u/sqrt(g h) (Ting & Kirby's Fig. 5, digitised; not
              distributed with NHFOAM)
    kind = nhfoam (wetdry.tingkirby with averaging: UMean and etaMean at the
           last time) or nhwave (results/u_* and eta_* of its example, x_TK =
           x - 5.6 m)
The models are averaged at fixed sigma; the mean height of a sigma level is
z_b + sigma (mean depth).  Below the trough this is close to the Eulerian
mean at fixed z measured in the flume.
(also writes out.pdf)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import other_models as om
import longshore as ls
import wetdry

G = 9.81
args = sys.argv[1:]
labfile, out = args[0], args[1]
runs = [om.parse_case(a) for a in args[2:]]
lab = np.loadtxt(labfile, comments="#")
panels = sorted(set(lab[:, 0].astype(int)))


def depth(x):
    return np.where(x < 0, wetdry.TK_D, wetdry.TK_D - wetdry.TK_SLOPE * x)


def nhfoam(case):
    t = ls.latest_time(case)
    C = ls.cell_centres(case, t)
    Um = ls._read_field(os.path.join(case, t, "UMean"))
    em = ls._read_field(os.path.join(case, t, "etaMean"))
    e = ls._read_field(os.path.join(case, t, "eta"))
    xr = np.round(C[:, 0], 6)
    xs = np.unique(xr)
    prof = {}
    for xst in sorted(set(lab[:, 1])):
        x = xs[np.argmin(abs(xs - xst))]
        m = np.where(xr == x)[0]
        m = m[np.argsort(C[m, 2])]
        h = float(depth(np.array(x)))
        zb = -h
        sig = (C[m, 2] - zb) / (e[m[0]] - zb)
        zbar = zb + sig * (em[m[0]] - zb)
        prof[xst] = ((zbar - em[m[0]]) / h, Um[m, 0] / np.sqrt(G * h))
    return prof


def nhwave(case, x0=5.6):
    p = {}
    for line in open(os.path.join(case, "input.txt")):
        if "=" in line and not line.startswith("!"):
            k, v = line.split("=", 1)
            p[k.strip()] = v.split()[0] if v.split() else ""
    Nx, Nz, dx = int(p["Mglob"]), int(p["Kglob"]), float(p["DX"])
    res = os.path.join(case, "results")
    uf = sorted(f for f in os.listdir(res) if f.startswith("u_"))
    ef = sorted(f for f in os.listdir(res) if f.startswith("eta_"))
    U = np.mean([np.loadtxt(os.path.join(res, f)).reshape(Nz, Nx) for f in uf], axis=0)
    E = np.mean([np.loadtxt(os.path.join(res, f)).reshape(-1)[:Nx] for f in ef], axis=0)
    hb = np.loadtxt(os.path.join(case, "depth.txt")).reshape(-1)[:Nx]
    x = (np.arange(Nx) + 0.5) * dx - x0
    sig = (np.arange(Nz) + 0.5) / Nz
    prof = {}
    for xst in sorted(set(lab[:, 1])):
        i = np.argmin(abs(x - xst))
        h = float(depth(np.array(x[i])))
        zbar = -hb[i] + sig * (hb[i] + E[i])
        prof[xst] = ((zbar - E[i]) / h, U[:, i] / np.sqrt(G * h))
    return prof


models = []
for kind, case, label in runs:
    models.append((label, nhwave(case) if kind == "nhwave" else nhfoam(case)))

fig, axs = plt.subplots(2, 4, figsize=(9.0, 6.2), sharey=True)
axs = axs.ravel()
styles = ["-", "--", "-.", ":"]
colors = ["C0", "C1", "C2", "C4"]
for ip, pn in enumerate(panels):
    ax = axs[ip]
    d = lab[lab[:, 0] == pn]
    xst = d[0, 1]
    ax.plot(d[:, 5], d[:, 4], "o", mfc="none", mec="C3", ms=4,
            label="Ting \\& Kirby (1994)" if ip == 0 else None)
    for n, (label, prof) in enumerate(models):
        z, u = prof[xst]
        ax.plot(u, z, color=colors[n % 4], ls=styles[n % 4], lw=1.3,
                label=label if ip == 0 else None)
    ax.axvline(0, color="0.6", lw=0.6)
    ax.set_title(f"$x = {xst:.2f}$ m", fontsize=10)
    ax.set_ylim(-1.0, 0.1)
    ax.grid(alpha=0.3)
    if ip % 4 == 0:
        ax.set_ylabel(r"$(z-\bar\eta)/h$")
    if ip >= 4:
        ax.set_xlabel(r"$\bar u/\sqrt{gh}$")
fig.legend(loc="lower center", ncol=3, fontsize=8,
           bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=(0, 0.09, 1, 1))
save(fig, out)
# rms difference below the trough ((z - mean level)/h < -0.25), inside the
# surf zone (x > 6.5 m); and the largest offshore velocity in each profile
for label, prof in models:
    errs, peak = [], []
    for pn in panels:
        d = lab[lab[:, 0] == pn]
        xst = d[0, 1]
        if xst < 6.5:
            continue
        d = d[d[:, 4] < -0.25]
        z, u = prof[xst]
        o = np.argsort(z)
        errs += list(np.interp(d[:, 4], z[o], u[o]) - d[:, 5])
        m = z < -0.25
        peak.append((u[m].min(), d[:, 5].min()))
    errs = np.array(errs)
    peak = np.array(peak)
    print("%s: undertow x > 6.5 m below the trough, rms difference %.4f, bias"
          " %+.4f (u/sqrt(gh), %d points); strongest undertow model/lab %.2f"
          % (label, np.sqrt(np.mean(errs ** 2)), errs.mean(), len(errs),
             np.mean(peak[:, 0] / peak[:, 1])))
