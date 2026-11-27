#!/usr/bin/env python3
"""Run-up onto a complex 3D beach, the Monai valley (NTHMP benchmark 7):
gauges, run-up in the valley and shorelines of NHFOAM and NHWAVE against
the lab data.

usage: plot_monai_valley.py <figure prefix> <run> [<run> ...]
    run = <path>:<label>, cases made by monai_valley.nhfoam / .nhwave
writes <prefix>_gauges, <prefix>_runup and <prefix>_frames (.png and .pdf)
and prints the crest, rms and run-up numbers of the README/report tables.

Lab gauges as given (they read +2 to +5 mm before the wave arrives); no
run is shifted in time.  Panel (a): the surface of the first column
against the input record (the lag of each model's wave there, by
monai_valley.lag, is printed).
"""
import sys
import numpy as np
from latex_style import use_latex, save
use_latex()
import matplotlib.pyplot as plt
import monai_valley as mv

NHFOAM_STYLES = [("-", "#1f5fa8", 1.3), (":", "#4a90d9", 1.6), ("-.", "#7f8c8d", 1.1)]
NHWAVE_STYLES = [("--", "#e67e22", 1.3), ((0, (5, 1.5, 1, 1.5)), "#c0392b", 1.1)]
T0, T1 = 0.0, mv.T_END
WIN = (13.0, 23.0)              # main wave at the gauges: crest and rms window


def styles(runs):
    out, nf, nw = {}, 0, 0
    for path, label in runs:
        if mv.params(path)["model"] == "nhwave":
            out[label] = NHWAVE_STYLES[min(nw, len(NHWAVE_STYLES) - 1)]
            nw += 1
        else:
            out[label] = NHFOAM_STYLES[min(nf, len(NHFOAM_STYLES) - 1)]
            nf += 1
    return out


def crest(t, e, t0, t1):
    """Height and time of the highest value in [t0, t1] (parabolic fit)."""
    s = np.where((t >= t0) & (t <= t1))[0]
    i = s[np.argmax(e[s])]
    if 0 < i < len(e) - 1:
        a, b, c = e[i - 1], e[i], e[i + 1]
        den = a - 2 * b + c
        off = 0.5 * (a - c) / den if den != 0 else 0.0
        return b - 0.25 * (a - c) * off, t[i] + off * (t[i + 1] - t[i])
    return e[i], t[i]


def gauges_figure(runs, sty, prefix):
    tl, gl = mv.read_lab()
    tr, er = mv.read_input()
    fig, axs = plt.subplots(4, 1, figsize=(8.0, 9.6), sharex=True)
    ax = axs[0]
    ax.plot(tr, 100 * er, "k-", lw=1.0, label="input record")
    print(f"{'run':26s} {'lag at x=0 [s]':>14s} {'rms 0-12 s [mm]':>16s} "
          f"{'max 9-12 s [mm]':>16s} (record {1000 * er[(tr >= 9) & (tr <= 12)].max():.2f} mm)")
    mods = []
    for path, label in runs:
        tg, vg = mv.model_gauges(path)
        st = sty[label]
        traw = tg["x0"]
        ax.plot(traw, 100 * vg["x0"], ls=st[0], color=st[1], lw=st[2], label=label)
        Lm = mv.lag(traw, vg["x0"])
        w = tr <= 12.0
        rms0 = np.sqrt(np.mean((np.interp(tr[w], traw - Lm, vg["x0"]) - er[w]) ** 2))
        mx = vg["x0"][(traw - Lm >= 9.0) & (traw - Lm <= 12.0)].max()
        print(f"{label:26s} {Lm:14.3f} {1000 * rms0:16.3f} {1000 * mx:16.2f}")
        mods.append((label, tg, vg))
    ax.set_ylabel(r"$\eta$ [cm]")
    ax.set_title(r"(a) first column, $x = 0$ to $\Delta x$, $y = 1.70$ m", loc="left")
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    print(f"\n{'run':26s} {'gauge':5s} {'crest [cm]':>10s} {'lab':>6s} {'err %':>6s} "
          f"{'t [s]':>6s} {'lab':>6s} {'bore [s]':>8s} {'lab':>6s} {'rms [cm]':>9s}")
    for k, g in enumerate(("ch5", "ch7", "ch9")):
        ax = axs[k + 1]
        ax.plot(tl, 100 * gl[g], "k-", lw=1.0, label="lab")
        HL, tL = crest(tl, gl[g], *WIN)
        aL = mv.arrival(tl, gl[g])
        wl = (tl >= WIN[0]) & (tl <= WIN[1])
        for label, tg, vg in mods:
            st = sty[label]
            ax.plot(tg[g], 100 * vg[g], ls=st[0], color=st[1], lw=st[2], label=label)
            Hm, tm = crest(tg[g], vg[g], *WIN)
            em = np.interp(tl[wl], tg[g], vg[g])
            rms = np.sqrt(np.mean((em - gl[g][wl]) ** 2))
            print(f"{label:26s} {g:5s} {100 * Hm:10.2f} {100 * HL:6.2f} "
                  f"{100 * (Hm - HL) / HL:6.1f} {tm:6.2f} {tL:6.2f} "
                  f"{mv.arrival(tg[g], vg[g]):8.2f} {aL:6.2f} {100 * rms:9.3f}")
        x, y = mv.GAUGES[g]
        ax.set_title(rf"({'bcd'[k]}) gauge {g[2:]}, $x = {x:g}$ m, $y = {y:g}$ m", loc="left")
        ax.set_ylabel(r"$\eta$ [cm]")
        if k == 0:
            ax.legend(loc="upper left", fontsize=9, frameon=False)
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.set_xlim(T0, T1)
    axs[-1].set_xlabel(r"$t$ [s]")
    fig.tight_layout()
    save(fig, prefix + "_gauges.png")


def valley_map(ax, dx=mv.DX0, box=(4.6, 5.35, 1.45, 2.45)):
    """Ground of the valley shaded by its height above the still water
    (contours every 2 cm), the still-water shoreline thick."""
    xn, yn, Zn = mv.nodes(dx)
    sx = (xn >= box[0] - 0.05) & (xn <= box[1] + 0.05)
    sy = (yn >= box[2] - 0.05) & (yn <= box[3] + 0.05)
    Z = 100 * Zn[np.ix_(sy, sx)]
    ax.contourf(xn[sx], yn[sy], Z, levels=[-20, 0], colors=["#e8f1fa"])
    ax.contourf(xn[sx], yn[sy], Z, levels=np.arange(0, 14.1, 1.0), cmap="Greys",
                vmin=-4, vmax=40, extend="max")
    cs = ax.contour(xn[sx], yn[sy], Z, levels=np.arange(2, 13, 2), colors="0.55",
                    linewidths=0.45)
    ax.clabel(cs, fmt="%d", fontsize=7, inline_spacing=1)
    ax.contour(xn[sx], yn[sy], Z, levels=[0.0], colors="0.25", linewidths=0.9)
    ax.set_xlim(box[0], box[1])
    ax.set_ylim(box[2], box[3])
    ax.set_aspect("equal")


def cell_xy(dx):
    nx, ny, dy = mv.grid(dx)
    return ((np.arange(nx) + 0.5) * (mv.LX / nx), (np.arange(ny) + 0.5) * dy)


def runup_figure(runs, sty, prefix):
    obs = mv.read_runup()
    fig, (axm, axr) = plt.subplots(1, 2, figsize=(11.0, 5.0),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    valley_map(axm)
    print(f"\n{'run':26s} {'point':16s} {'R [cm]':>7s} {'x_in [m]':>9s} {'lab [cm]':>12s}")
    ys = None
    for path, label in runs:
        st = sty[label]
        pr = mv.params(path)
        e, fl, z = mv.emax(path)
        xc, yc = cell_xy(pr["dx"])
        axm.contour(xc, yc, fl.astype(float), levels=[0.5], colors=[st[1]],
                    linestyles=[st[0]], linewidths=st[2])
        axm.plot([], [], ls=st[0], color=st[1], lw=st[2], label=label)
        # run-up along the valley: highest flooded bed of each row
        X, Y = np.meshgrid(xc, yc)
        land = fl & (z > 0.0) & (X > 4.6) & (X < 5.4)
        R = np.where(land, z, -np.inf).max(axis=1)
        sel = (yc >= 1.45) & (yc <= 2.45)
        axr.plot(yc[sel], 100 * np.where(np.isfinite(R[sel]), R[sel], np.nan),
                 ls=st[0], color=st[1], lw=st[2], label=label)
        for r in mv.runup(path):
            print(f"{label:26s} ({r['x']:.3f}, {r['y']:.3f}) {100 * r['R']:7.2f} "
                  f"{r['x_inland']:9.3f} {100 * r['obs'].min():5.2f}-{100 * r['obs'].max():5.2f}")
        v = (yc >= 1.75) & (yc <= 2.05) & np.isfinite(R)
        k = np.argmax(np.where(v, R, -np.inf))
        area = (fl & (z > 0.0)).sum() * (mv.LX / len(xc)) * (yc[1] - yc[0])
        print(f"{label:26s} highest flooded ground in the valley (1.75 < y < 2.05): "
              f"{100 * R[k]:.2f} cm at y = {yc[k]:.3f}; flooded land {area:.3f} m^2")
    for k, (x, y, o) in enumerate(obs):
        axm.plot(x, y, "ko", ms=4, label="lab run-up" if k == 0 else None)
        axr.errorbar(y, 100 * o.mean(), yerr=[[100 * (o.mean() - o.min())],
                                               [100 * (o.max() - o.mean())]],
                     fmt="ko", ms=4, capsize=3, label="lab (6 repeats)" if k == 0 else None)
    axm.set_xlabel(r"$x$ [m]")
    axm.set_ylabel(r"$y$ [m]")
    axm.set_title("(a) highest inundation in the valley", loc="left")
    axm.legend(loc="lower left", fontsize=8.5, frameon=True, framealpha=0.85)
    axr.set_xlabel(r"$y$ [m]")
    axr.set_ylabel(r"run-up [cm]")
    axr.set_title(r"(b) highest flooded ground, $4.6 < x < 5.4$ m", loc="left")
    axr.set_xlim(1.45, 2.45)
    axr.grid(alpha=0.3)
    axr.legend(loc="upper right", fontsize=9, frameon=False)
    fig.tight_layout()
    save(fig, prefix + "_runup.png")


def frames_figure(runs, sty, prefix, times=mv.FRAMES):
    box = (4.7, 5.25, 1.5, 2.25)
    fig, axs = plt.subplots(1, len(times), figsize=(13.5, 4.3), sharey=True)
    for k, t in enumerate(times):
        ax = axs[k]
        valley_map(ax, box=box)
        for path, label in runs:
            st = sty[label]
            pr = mv.params(path)
            e, z, wet, ts = mv.surface(path, t)
            xc, yc = cell_xy(pr["dx"])
            ax.contour(xc, yc, wet.astype(float), levels=[0.5], colors=[st[1]],
                       linestyles=[st[0]], linewidths=st[2])
            if k == 0:
                ax.plot([], [], ls=st[0], color=st[1], lw=st[2], label=label)
        ax.set_title(rf"$t = {t:g}$ s (frame {10 + 15 * k})")
        ax.set_xlabel(r"$x$ [m]")
        ax.set_xticks([4.8, 5.0, 5.2])
    axs[0].set_ylabel(r"$y$ [m]")
    fig.legend(loc="lower center", ncol=len(runs), frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, prefix + "_frames.png")


def main():
    prefix = sys.argv[1]
    runs = [tuple(a.split(":", 1)) for a in sys.argv[2:]]
    sty = styles(runs)
    gauges_figure(runs, sty, prefix)
    runup_figure(runs, sty, prefix)
    frames_figure(runs, sty, prefix)


if __name__ == "__main__":
    main()
