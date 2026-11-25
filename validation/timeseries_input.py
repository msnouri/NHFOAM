#!/usr/bin/env python3
"""Time-series wave input (regularWave, theory timeSeries) on a flat 2DV
tank: a solitary wave made at x = 0 from its own record, against the exact
wave of Babenko's equation travelling at its celerity c.

modes
  exact    the exact wave as the initial condition, crest at x = XC0 (the
           model's own propagation error, for reference)
  eta      record of the elevation at x = 0 with the long-wave velocity
           u = eta sqrt(g/(h + eta)) of the NTHMP benchmarks, active
           absorption on
  etaC     the same with u = c eta/(h + eta), c the celerity of the wave:
           the exact depth-averaged velocity of a wave of permanent form
  paddle   the piston displacement X(t), the time integral of the etaC
           velocity, without absorption (a laboratory paddle that stays at
           its rest position)

The record is sampled every RDT = 0.02 s (a typical gauge rate) and passes
through the input as a table (linear interpolation).

usage (python):
    import timeseries_input as ti
    ti.flat_case("ts_eta_10", "eta", Hh=0.1)     # then ./Allrun
    ti.analyse(["ts_exact_10", "ts_eta_10", ...], Hh=0.1)
"""
import os, math
import numpy as np
import make_case as mc

G = 9.81
H0 = 0.218                       # depth [m] (the composite beach, BP5)
L, SPONGE = 18.0, 3.0
GAUGES = [0.5, 4.0, 8.0, 12.0]
XC0 = 3.5                        # crest of the initial wave (mode exact)
RDT = 0.02


def wave(Hh, h=H0):
    from solitary_exact import SolitaryWave
    return SolitaryWave(Hh * h, d=h)


def crest_time(sw, h=H0, level=1e-3):
    """Time at which the crest passes x = 0 for a record starting at
    eta = level*H (first-order decay rate)."""
    gamma = math.sqrt(0.75 * sw.eta(np.array([0.0]))[0] / h)
    return h * math.acosh(math.sqrt(1.0 / level)) / gamma / sw.c


def records(sw, h=H0, rdt=RDT, taper=0.3):
    """Record of the elevation at x = 0 and the paddle displacement whose
    velocity is c eta/(h + eta), both from 0 to twice the crest time,
    tapered to zero at both ends."""
    t0 = crest_time(sw, h)
    t = np.arange(0.0, 2 * t0 + rdt / 2, rdt)
    eta = sw.eta(-sw.c * (t - t0))
    w = np.ones_like(t)
    a = t < taper
    w[a] = np.sin(0.5 * np.pi * t[a] / taper) ** 2
    b = t > t[-1] - taper
    w[b] = np.sin(0.5 * np.pi * (t[-1] - t[b]) / taper) ** 2
    eta = eta * w
    u = sw.c * eta / (h + eta)
    X = np.concatenate([[0.0], np.cumsum(0.5 * (u[1:] + u[:-1]) * rdt)])
    return t, eta, X, t0


def write_table(path, t, v):
    with open(path, "w") as f:
        f.write(f"{len(t)}\n(\n")
        f.writelines(f"({a:.6f} {b:.10e})\n" for a, b in zip(t, v))
        f.write(")\n")


def flat_case(case, mode, Hh=0.1, h=H0, dx=0.02, nz=5, dt=0.005, endTime=None):
    sw = wave(Hh, h)
    t, eta, X, t0 = records(sw, h)
    endTime = endTime or (t0 + GAUGES[-1] / sw.c + t0 + 1.0)
    Nx = int(round(L / dx))
    dy = dx
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in GAUGES)
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({L - SPONGE:g} 0 0);
        end         ({L:g} 0 0);
        sigmaMax    20;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps,
                    nCorr=3, pTol=1e-10)
    mc.blockmesh_flat(case, L, h, dy, Nx, nz)
    mc.common_fields(case)
    mc.allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)
    if mode == "exact":
        import ma2012
        i = np.arange(Nx)
        k = np.arange(nz)
        xc = np.tile((i + 0.5) * dx, nz)
        sigma = np.repeat((k + 0.5) / nz, Nx)
        e = sw.eta(xc - XC0)
        zc = -h + sigma * (h + e)
        u, w, prgh = ma2012.solitary_fields(sw, xc - XC0, zc, h)
        ma2012._set_internal(os.path.join(case, "0", "eta"), ma2012._foam_list(e))
        ma2012._set_internal(os.path.join(case, "0", "U"),
                             ma2012._foam_list(list(zip(u, w)), vector=True))
        ma2012._set_internal(os.path.join(case, "0", "p_rgh"), ma2012._foam_list(prgh))
    else:
        rec = "paddle" if mode == "paddle" else "eta"
        write_table(os.path.join(case, "constant", "record.dat"), t,
                    X if mode == "paddle" else eta)
        cel = f"celerity    {sw.c:.10g};" if mode == "etaC" else ""
        text = f"""theory      timeSeries;
        depth       {h:g};
        direction   (1 0 0);
        {rec}
        {{
            type    tableFile;
            file    "$FOAM_CASE/constant/record.dat";
        }}
        {cel}"""
        mc.set_patch(case, "U", "left", mc.wave_bc("left", text, absorb=(mode != "paddle")))
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"mode {mode}\nHh {Hh}\nh {h}\nc {sw.c:.10g}\nt0 {t0:.10g}\n")
    return dict(c=sw.c, t0=t0, endTime=endTime)


def expected(sw, mode, x, t, t0):
    """Exact elevation at x and times t for the mode's start."""
    if mode == "exact":
        return sw.eta(x - XC0 - sw.c * t)
    return sw.eta(x - sw.c * (t - t0))


def analyse(cases):
    """Peak height relative to H, arrival of the peak against the exact
    wave [s], rms error over the pulse relative to H, and the largest
    elevation behind the pulse (trailing waves) relative to H."""
    import other_models as om
    rows = []
    for case in cases:
        p = dict(l.split() for l in open(os.path.join(case, "params.txt")))
        mode, Hh, h = p["mode"], float(p["Hh"]), float(p["h"])
        sw = wave(Hh, h)
        H = Hh * h
        t0 = float(p["t0"])
        t, g = om.read_nhfoam(case)
        for x in sorted(g):
            if mode == "exact" and x < XC0 + 0.4:
                continue
            e = g[x]
            ex = expected(sw, mode, x, t, t0)
            tp_ex = t[np.argmax(ex)]
            i = np.argmax(e)
            # parabola through the three samples around the peak
            if 0 < i < len(e) - 1:
                a, b, c3 = e[i - 1], e[i], e[i + 1]
                den = a - 2 * b + c3
                s = 0.5 * (a - c3) / den if den != 0 else 0.0
                tp = t[i] + s * (t[1] - t[0])
                peak = b - 0.25 * (a - c3) * s
            else:
                tp, peak = t[i], e[i]
            win = np.abs(t - tp_ex) < 2.0 * crest_time(sw, h, 0.01)
            rms = math.sqrt(np.mean((e[win] - ex[win]) ** 2)) / H
            # behind the pulse, before waves from the sponge come back
            behind = ((t > tp_ex + 1.2 * crest_time(sw, h))
                      & (t < tp_ex + 2.0 * crest_time(sw, h) + 1.0))
            trail = np.max(np.abs(e[behind])) / H if behind.any() else float("nan")
            rows.append((case, mode, Hh, x, peak / H, tp - tp_ex, rms, trail))
    return rows


if __name__ == "__main__":
    import sys
    for r in analyse(sys.argv[1:]):
        print("%-14s %-6s H/h=%.2f x=%5.1f  peak/H=%.4f  dt_peak=%+.4f s  rms/H=%.4f  trail/H=%.4f" % r)
