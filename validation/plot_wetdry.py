#!/usr/bin/env python3
"""Wetting-and-drying verification figures (cases made by wetdry.py).

usage:
    plot_wetdry.py thacker  <case> <out.png>
        Thacker (1981) planar surface in a parabolic basin: surface profiles
        and shorelines against the exact solution; prints the errors.
    plot_wetdry.py dambreak <out.png> [kind:]<case>[:label] ...
        Dam break on a dry bed (wetdry.dambreak; kind nhfoam or nhwave):
        depth and velocity at t = 1, 2, 3 s against Ritter's solution and
        Stoker's solution for a wet bed as deep as the film.
(each also writes out.pdf)
"""
import os
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import wetdry
import make_case as mc
import other_models as om

G = wetdry.G


def foam_list(path, vector=False):
    s = open(path).read().split("internalField")[1]
    first = s.split("\n")[0].split()
    if first and first[0] == "uniform":
        return None
    n = int(s.split("(", 1)[0].split()[-1])
    body = s.split("(", 1)[1]
    if vector:
        return np.array(body.replace("(", " ").replace(")", " ").split()[:3 * n],
                        float).reshape(n, 3)
    return np.array(body.split(")", 1)[0].split(), float)


def column_index(xs, dx, nz):
    """x of each column, and the column and layer of every cell of a
    blockmesh_profile mesh."""
    xc = wetdry.profile_cell_x(xs, dx, nz)
    col = np.zeros(len(xc), int)
    lay = np.zeros(len(xc), int)
    pos = c0 = 0
    for i in range(len(xs) - 1):
        nx = max(1, int(round((xs[i + 1] - xs[i]) / dx)))
        for k in range(nz):
            col[pos + k * nx:pos + (k + 1) * nx] = c0 + np.arange(nx)
            lay[pos + k * nx:pos + (k + 1) * nx] = k
        pos += nx * nz
        c0 += nx
    x = np.zeros(c0)
    x[col] = xc
    return x, col, lay


def times(case):
    out = []
    for t in os.listdir(case):
        try:
            if float(t) > 0:
                out.append(t)
        except ValueError:
            pass
    return sorted(out, key=float)


# --------------------------------------------------------------------------
def thacker(case, out, a=100.0, h0=0.5, B=0.5, L=400.0, dx=1.0, nz=3, film=0.001):
    xs0 = np.linspace(0, L, 81)
    zb0 = h0 * (((xs0 - L / 2) / a) ** 2 - 1)
    xs, hs, _ = mc.beach_profile(list(xs0), list(-zb0), film, dx=dx)
    x, col, lay = column_index(xs, dx, nz)
    zb = -np.interp(x, xs, hs)
    _, _, _, T = wetdry.thacker_exact(x, 0.0, a, h0, B, L)
    tt, xl, xr, rows = [], [], [], []
    snaps = {}
    for t in times(case):
        e = np.zeros(len(x))
        e[col] = foam_list(os.path.join(case, t, "eta"))
        tv = float(t)
        ex, _, shore, _ = wetdry.thacker_exact(x, tv, a, h0, B, L)
        wet = e - zb > film + 1e-6
        both = wet & np.isfinite(ex)
        err = e[both] - ex[both]
        tt.append(tv / T); xl.append(x[wet].min()); xr.append(x[wet].max())
        rows.append((tv / T, np.sqrt(np.mean(err ** 2)), np.abs(err).max(),
                     x[wet].min() - shore[0], x[wet].max() - shore[1]))
        snaps[round(tv / T, 3)] = np.where(wet, e, np.nan)
    print("t/T     rms [mm]  max [mm]  left shore - exact [m]  right [m]")
    for r in rows:
        print("%5.3f  %8.2f  %8.2f  %10.2f  %10.2f" % (r[0], 1e3 * r[1], 1e3 * r[2], r[3], r[4]))
    r = np.array(rows)
    print("over all: rms %.2f mm (max of the rms %.2f), max %.2f mm, shoreline error"
          " max %.2f m (mean |.| %.2f m)" % (1e3 * np.sqrt(np.mean(r[:, 1] ** 2)),
                                             1e3 * r[:, 1].max(), 1e3 * r[:, 2].max(),
                                             np.abs(r[:, 3:5]).max(), np.abs(r[:, 3:5]).mean()))

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.5),
                             gridspec_kw=dict(height_ratios=[1.4, 1]))
    ax = axes[0]
    xf = np.linspace(0, L, 801)
    zf = h0 * (((xf - L / 2) / a) ** 2 - 1)
    ax.fill_between(xf, -0.6, zf, color="0.85", lw=0)
    cols = ["#1f5fa8", "#e67e22", "#2e8b57", "#8e44ad"]
    for k, (ts0, c) in enumerate(zip([1.25, 1.5, 1.75, 2.0], cols)):
        ts = min(snaps, key=lambda v: abs(v - ts0))
        ex, _, _, _ = wetdry.thacker_exact(xf, ts * T, a, h0, B, L)
        ax.plot(xf, ex, "-", color=c, lw=0.9)
        sel = slice(None, None, 6)
        ax.plot(x[sel], snaps[ts][sel], "o", ms=3, mfc="none", color=c,
                label=rf"$t/T = {ts:.2f}$")
    ax.set_xlim(60, 340)
    ax.set_ylim(-0.12, 0.12)
    ax.set_ylabel(r"$\eta$ [m]")
    ax.set_xlabel(r"$x$ [m]")
    ax.legend(ncol=4, fontsize=10, loc="upper center", frameon=False)
    ax.set_title("surface: exact (lines) and \\textsc{NHFOAM} (circles)" if
                 plt.rcParams["text.usetex"] else "surface: exact (lines) and NHFOAM (circles)",
                 fontsize=11)
    ax = axes[1]
    tf = np.linspace(0, max(tt), 400)
    w = np.sqrt(2 * G * h0) / a
    beta = B / (a * w)
    cc = np.cos(w * tf * T)
    ax.plot(tf, L / 2 + a * (-beta * cc - 1), "k-", lw=0.9, label="exact")
    ax.plot(tf, L / 2 + a * (-beta * cc + 1), "k-", lw=0.9)
    ax.plot(tt, xl, "o", ms=4, mfc="none", color="#1f5fa8", label="NHFOAM, left")
    ax.plot(tt, xr, "s", ms=4, mfc="none", color="#e67e22", label="NHFOAM, right")
    ax.set_xlabel(r"$t/T$")
    ax.set_ylabel("shoreline $x$ [m]")
    ax.legend(ncol=3, fontsize=10, frameon=False, loc="center right")
    fig.tight_layout()
    save(fig, out)


# --------------------------------------------------------------------------
def stoker(x, t, h0, hR, x0):
    """Dam break on a wet bed of depth hR (Stoker 1957)."""
    from scipy.optimize import brentq
    c0 = np.sqrt(G * h0)
    f = lambda hm: 2 * (c0 - np.sqrt(G * hm)) - (hm - hR) * np.sqrt(G * (hm + hR) / (2 * hm * hR))
    hm = brentq(f, hR * 1.0001, h0)
    um = 2 * (c0 - np.sqrt(G * hm))
    s = hm * um / (hm - hR)
    xi = (x - x0) / t
    h = np.where(xi < -c0, h0, (2 * c0 - xi) ** 2 / (9 * G))
    u = np.where(xi < -c0, 0.0, 2.0 / 3.0 * (xi + c0))
    mid = (xi >= um - np.sqrt(G * hm)) & (xi < s)
    h = np.where(mid, hm, h); u = np.where(mid, um, u)
    h = np.where(xi >= s, hR, h); u = np.where(xi >= s, 0.0, u)
    return h, u, s


def dambreak(out, runs, h0=0.1, x0=5.0, dx=0.02, nz=3, film=0.001, L=15.0):
    xc = (np.arange(int(round(L / dx))) + 0.5) * dx
    n = len(xc)
    tlist = [1.0, 2.0, 3.0]
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.0), sharex=True)
    styles = [("-", "#1f5fa8", 1.4), ("--", "#e67e22", 1.3), ("-.", "#7f8c8d", 1.3),
              (":", "#2e8b57", 1.6)]
    fronts = {}
    for t in tlist:
        h, u = wetdry.ritter_exact(xc, t, h0, x0)
        hs, us, s = stoker(xc, t, h0, film, x0)
        lr = axes[0].plot(xc, h, "k-", lw=0.7, label="Ritter (dry bed)" if t == 1 else None)
        axes[1].plot(xc, u, "k-", lw=0.7)
        axes[0].plot(xc, hs, "k:", lw=1.0, label="Stoker, bed as deep as the film" if t == 1 else None)
        axes[1].plot(xc, us, "k:", lw=1.0)
        for (kind, path, label), (ls, c, lw) in zip(runs, styles):
            if kind == "nhwave":
                k = int(round(t / 0.5))
                e = np.loadtxt(os.path.join(path, "results", f"eta_{k:05d}")).ravel()
                uu = np.loadtxt(os.path.join(path, "results", f"u_{k:05d}")).reshape(-1, n).mean(axis=0)
                D = e + h0
            else:
                e = foam_list(os.path.join(path, f"{t:g}", "eta"))
                U = foam_list(os.path.join(path, f"{t:g}", "U"), vector=True)
                D = e[:n] + h0
                uu = U[:, 0].reshape(nz, n).mean(axis=0)
            wet = D > 1.1 * film
            axes[0].plot(xc, np.where(wet, D, np.nan), ls=ls, color=c, lw=lw,
                         label=label if t == 1 else None)
            axes[1].plot(xc, np.where(wet, uu, np.nan), ls=ls, color=c, lw=lw)
            fronts.setdefault(label, []).append(xc[wet].max())
        fronts.setdefault("Ritter", []).append(x0 + 2 * np.sqrt(G * h0) * t)
        fronts.setdefault("Stoker (film)", []).append(x0 + s * t)
        axes[0].text(x0 + 2 * np.sqrt(G * h0) * t, 0.006, rf"$t = {t:g}$ s",
                     ha="center", fontsize=10)
    print("front position [m] at t = " + ", ".join("%g" % t for t in tlist))
    for k, v in fronts.items():
        print("  %-28s " % k + "  ".join("%.2f" % x for x in v))
    axes[0].set_ylabel(r"$h$ [m]")
    axes[1].set_ylabel(r"$u$ [m/s]")
    axes[1].set_xlabel(r"$x$ [m]")
    axes[0].set_xlim(1.5, 11.5)
    axes[0].legend(fontsize=10, frameon=False, loc="upper right")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, out)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "thacker":
        thacker(sys.argv[2], sys.argv[3])
    elif what == "dambreak":
        dambreak(sys.argv[2], [om.parse_case(a) for a in sys.argv[3:]])
    else:
        sys.exit(__doc__)
