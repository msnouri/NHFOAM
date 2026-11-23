#!/usr/bin/env python3
"""Wave-input tests of wavemaker.py: figures and the numbers of the report.

usage: plot_wavemaker.py <root> <outDir>
    root: the cases written by wavemaker.py (run)
Writes <outDir>/wavemaker_regular.{png,pdf}, wavemaker_absorber.{png,pdf},
wavemaker_irregular.{png,pdf} and prints the tables.
"""
import os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save

use_latex()
G, h = 9.81, 0.4
root, outdir = sys.argv[1], sys.argv[2]


def gauges(case):
    fn = os.path.join(root, case, "postProcessing/gauges/0/eta.dat")
    d = np.loadtxt(fn, comments="#")
    xs = [float(l.split()[2][1:]) for l in open(fn) if l.startswith("# G")]
    return d[:, 0], np.array(xs), d[:, 1:]


def wavenumber(w, h=h):
    k = w * w / G
    for _ in range(60):
        k -= (G * k * np.tanh(k * h) - w * w) / (G * np.tanh(k * h) + G * k * h / np.cosh(k * h) ** 2)
    return k


def harmonic(t, e, w, t1, t2):
    """Complex amplitudes A (eta = Re{A exp(-i w t)}) over [t1, t2]."""
    m = (t > t1) & (t <= t2)
    return 2 * np.mean(e[m] * np.exp(1j * w * t[m])[:, None], axis=0)


def reflection(x, A, k):
    """Least squares A_j = A_I exp(i k x_j) + A_R exp(-i k x_j), with k
    adjusted (+-5%) to the model's own wave number."""
    best = None
    for kk in np.linspace(0.95 * k, 1.05 * k, 401):
        M = np.stack([np.exp(1j * kk * x), np.exp(-1j * kk * x)], 1)
        c = np.linalg.lstsq(M, A, rcond=None)[0]
        r = np.linalg.norm(M @ c - A)
        if best is None or r < best[0]:
            best = (r, kk, c)
    r, kk, c = best
    return abs(c[0]), abs(c[1]), kk


def volume(case):
    t, v, tc = [], [], None
    for l in open(os.path.join(root, case, "log.NHFOAM")):
        if l.startswith("Time = "):
            tc = float(l.split()[2])
        elif "water volume change" in l and tc is not None:
            t.append(tc); v.append(float(l.split()[-1]))
    return np.array(t), np.array(v)


def done(case):
    f = os.path.join(root, case, "log.NHFOAM")
    return os.path.exists(f) and "\nEnd\n" in open(f).read()


# ---------------------------------------------------------------- regular
T, H = 1.5, 0.02
w = 2 * np.pi / T
k = wavenumber(w)
fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.4))
ax = axes[0]
print("Regular waves, H = 2 cm, T = 1.5 s, kh = %.2f" % (k * h))
for case, off, lab, col in (("reg_bc", 0.0, "boundary wavemaker", "#1f5fa8"),
                            ("reg_zone", 6.0, "relaxation zone (6 m)", "#e67e22")):
    if not done(case):
        continue
    t, x, e = gauges(case)
    A = harmonic(t, e, w, t[-1] - 20 * T, t[-1])
    aI, aR, kk = reflection(x, A, k)
    Hx = 2 * np.abs(A)
    ax.plot(x - off, Hx / H, "-", color=col, lw=1.4, label=lab)
    print(f"  {lab}: H/H_target = {Hx.mean()/H:.4f} (min {Hx.min()/H:.4f}, max {Hx.max()/H:.4f}),"
          f" 2a_I/H = {2*aI/H:.4f}, K_r = {aR/aI:.4f}, k/k_lin = {kk/k:.4f}")
ax.set_xlabel(r"distance from the wavemaker (or the inner end of the zone) [m]")
ax.set_ylabel(r"$H/H_\mathrm{target}$")
ax.set_ylim(0.96, 1.02)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="lower left")
ax.set_title(r"(a) Progressive waves, $H = 2$ cm, $T = 1.5$ s, $h = 0.4$ m ($kh = 0.96$)",
             fontsize=10)

ax1, ax2 = axes[1], axes[2]
print("Closed 10 m flume, boundary wavemaker at x = 0, wall at x = 10 m")
for case, lab, col, ls in (("wall_abs", "active absorption", "#1f5fa8", "-"),
                           ("wall_noabs", "no absorption", "#c0392b", "--")):
    if not done(case):
        continue
    t, x, e = gauges(case)
    # incident and reflected amplitude over sliding 4-period windows
    nper = int(round(t[-1] / T))
    tw, aIw, aRw = [], [], []
    for i in range(4, nper + 1):
        A = harmonic(t, e, w, (i - 4) * T, i * T)
        aI, aR, kk = reflection(x, A, k)
        tw.append(i * T); aIw.append(aI); aRw.append(aR)
    tw, aIw, aRw = map(np.array, (tw, aIw, aRw))
    ax1.plot(tw, aIw / (H / 2), ls, color=col, lw=1.4, label=lab)
    print(f"  {lab}: a_I/a over the last 4 periods = {aIw[-1]/(H/2):.3f} (largest"
          f" {aIw.max()/(H/2):.3f}, smallest after 30 s {aIw[tw > 30].min()/(H/2):.3f}),"
          f" a_R/a_I = {aRw[-1]/aIw[-1]:.3f}")
ax1.axhline(1.0, color="k", lw=0.8, ls="-.")
ax1.set_ylabel(r"incident amplitude $a_I/a$")
ax1.set_xlabel(r"$t$ [s]")
ax1.grid(alpha=0.3)
ax1.legend(fontsize=9, loc="lower left")
ax1.set_title(r"(b) Closed flume: wavemaker at $x=0$, wall at $x = 10$ m", fontsize=10)

print("Mean water level of closed flumes (mean over the gauges, one-period running mean)")
for case, lab, col, ls in (("zone_wall_norc", "relaxation zone, no return current", "#c0392b", "--"),
                           ("zone_wall_rc", "relaxation zone, return current", "#2e8b57", "-"),
                           ("wall_abs", "boundary wavemaker, absorption", "#1f5fa8", "-.")):
    if not done(case):
        continue
    t, x, e = gauges(case)
    em = e.mean(axis=1)
    n = int(round(T / np.mean(np.diff(t))))
    rm = np.convolve(em, np.ones(n) / n, mode="valid")
    tr = t[n - 1:]
    ax2.plot(tr, 1000 * rm, ls, color=col, lw=1.3, label=lab)
    late = tr > tr[-1] - 10 * T
    print(f"  {lab}: mean level over the last 10 periods {1000*rm[late].mean():+.2f} mm,"
          f" range after 20 s {1000*rm[tr > 20].min():+.2f} to {1000*rm[tr > 20].max():+.2f} mm")
M = G * (H / 2) ** 2 / (2 * w / k)
ax2.set_ylabel(r"mean level [mm]")
ax2.set_xlabel(r"$t$ [s]")
ax2.grid(alpha=0.3)
ax2.legend(fontsize=9, loc="upper left")
ax2.set_title(r"(c) Mean water level of closed flumes (zone: 0--6 m, wall at 16 m;"
              r" wavemaker: as in (b))", fontsize=10)
fig.tight_layout()
save(fig, os.path.join(outdir, "wavemaker_regular.png"))

# ---------------------------------------------------------------- absorber
print("Absorbing boundary (nhWaveVelocity, zero wave height) at x = 14 m, generation in 0-6 m")
rows = []
for Tn in (1.0, 1.5, 2.5):
    tag = f"{Tn:g}".replace(".", "p")
    wn = 2 * np.pi / Tn
    kn = wavenumber(wn)
    row = [kn * h]
    for case in (f"abs_T{tag}_shallow", f"abs_T{tag}_wave", f"sponge_T{tag}"):
        if not done(case):
            row.append(np.nan); continue
        t, x, e = gauges(case)
        A = harmonic(t, e, wn, t[-1] - 10 * Tn, t[-1])
        aI, aR, kk = reflection(x, A, kn)
        row.append(aR / aI)
    rows.append(row)
    print(f"  T = {Tn:g} s, kh = {kn*h:.2f}: K_r shallow = {row[1]:.3f}, wave = {row[2]:.3f},"
          f" sponge (8 m) = {row[3]:.3f}")
rows = np.array(rows)
fig, ax = plt.subplots(figsize=(6.0, 3.4))
khs = np.linspace(0.3, 2.0, 100)
# shallow-water absorption of a linear wave: depth-averaged velocity c/h eta
# against the imposed sqrt(g/h) eta gives K_r = |1 - c/sqrt(gh)|/(1 + c/sqrt(gh))
ch = np.sqrt(np.tanh(khs) / khs)
ax.plot(khs, (1 - ch) / (1 + ch), "k:", lw=0.9,
        label=r"$(1 - c/\sqrt{gh})/(1 + c/\sqrt{gh})$")
ax.plot(rows[:, 0], rows[:, 1], "o-", color="#c0392b", label=r"absorbing boundary, $\sqrt{gh}$")
ax.plot(rows[:, 0], rows[:, 2], "s-", color="#1f5fa8", label=r"absorbing boundary, $c$ of the wave")
ax.plot(rows[:, 0], rows[:, 3], "^-", color="#2e8b57", label=r"sponge layer, 8 m")
ax.set_xlabel(r"$kh$")
ax.set_ylabel(r"reflection coefficient $K_r$")
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
save(fig, os.path.join(outdir, "wavemaker_absorber.png"))

# ---------------------------------------------------------------- irregular
Hs, Tp, gam = 0.03, 1.5, 3.3
fp = 1 / Tp


def jonswap(f):
    s = np.where(f <= fp, 0.07, 0.09)
    return f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gam ** np.exp(-(f - fp) ** 2 / (2 * (s * fp) ** 2))


fmin, fmax, nf = 0.5 * fp, 3 * fp, 100
df = (fmax - fmin) / nf
fc = fmin + (np.arange(nf) + 0.5) * df
S0 = jonswap(fc)
S0 *= (Hs / 4) ** 2 / (S0.sum() * df)
fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4))
print("Irregular waves: JONSWAP Hs = 3 cm, Tp = 1.5 s, gamma = 3.3, 100 components")
from scipy.signal import welch
ax, axh = axes
ax.plot(fc, 1e4 * S0, "k-", lw=2.2, alpha=0.3, label="target")
for case, off, lab, col in (("irr_zone", 6.0, "zone, TVD", "#e67e22"),
                            ("irr_bc", 0.0, "wavemaker, TVD", "#7f8c8d"),
                            ("irr_bc_flt", 0.0, "wavemaker, filtered central", "#1f5fa8")):
    if not done(case):
        continue
    t, x, e = gauges(case)
    m = t > 30.0
    fs = 1 / np.mean(np.diff(t))
    hm0, tp = [], []
    for j in range(len(x)):
        f, S = welch(e[m, j] - e[m, j].mean(), fs=fs, nperseg=int(round(40.0 * fs)))
        band = (f > fmin) & (f < fmax)
        hm0.append(4 * np.sqrt(np.trapezoid(S[band], f[band])))
        tp.append(1 / f[np.argmax(S)])
        if case != "irr_bc" and j in (0, len(x) - 1):
            ax.plot(f, 1e4 * S, "-" if j == 0 else "--", color=col, lw=1.0,
                    label=f"{lab}, {x[j]-off:g} m")
    hm0 = np.array(hm0)
    axh.plot(x - off, hm0 / Hs, "o-", color=col, ms=4, lw=1.2, label=lab)
    print(f"  {lab}: Hm0/Hs at x - x0 = " + ", ".join(f"{xx-off:g} m: {hh/Hs:.3f}" for xx, hh in zip(x, hm0))
          + f"; Tp = {np.mean(tp):.3f} s")
ax.set_xlim(0, 2.2)
ax.set_xlabel(r"$f$ [Hz]")
ax.set_ylabel(r"$S$ [cm$^2$/Hz]")
ax.set_title(r"(a) spectra 2 and 12 m from the wavemaker", fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=7)
axh.axhline(1.0, color="k", lw=0.8, ls="-.")
axh.set_xlabel(r"distance from the wavemaker (or zone) [m]")
axh.set_ylabel(r"$H_{m0}/H_s$")
axh.set_title(r"(b) significant wave height", fontsize=10)
axh.grid(alpha=0.3)
axh.legend(fontsize=8)
fig.tight_layout()
save(fig, os.path.join(outdir, "wavemaker_irregular.png"))
