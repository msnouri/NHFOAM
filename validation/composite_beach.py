#!/usr/bin/env python3
"""Solitary waves on a composite beach with a vertical wall (Revere Beach
model of the Coastal and Hydraulics Laboratory; Briggs et al. 1995; NTHMP
benchmark problem 5), NHFOAM with three ways of making the wave, and NHWAVE.

Lab: a flume of depth d = 0.218 m from the wavemaker (X = 0, its rest
position) to X = 15.04 m, then slopes 1:53 (4.36 m), 1:150 (2.93 m) and
1:13 (0.90 m) up to a vertical wall at X = 23.23 m, in 4.7 cm of water.
Gauge 4 (the incident gauge) at X = 12.64 / 14.06 / 14.40 m for cases
A / B / C, gauges 5-10 at 15.04, 17.22, 19.40, 20.86, 22.33 and 22.80 m.
Target H/d = 0.05, 0.30, 0.70, measured at gauge 4 0.038, 0.26 and 0.68.
Data (tsunami3_runup.zip of the NTHMP benchmark set): ts3a/b/c.txt, the
gauges from 265 to 295 s; fdbk3abc.txt, the paddle from 258 to 268 s;
run3abc.txt, the highest run-up on the wall.  The gauges are taken
relative to their mean over the first 2 s (this removes the 1.5 mm offset
of gauge 4 in case C that the benchmark notes, and the smaller ones).

NHFOAM modes
  gauge4    domain from gauge 4 to the wall; the direct pulse of gauge 4
            (to 275 s) as a time series at the inflow, with the active
            absorption: the way the benchmark recommends, with its velocity
            u = eta sqrt(g/(d + eta)).  celerity="pulse" uses
            u = c eta/(d + eta), c = sqrt(g (d + H)) of the measured pulse
            height H instead (better up to H/d = 0.3, worse above:
            timeseries_input.py)
  paddle    the whole flume, from the middle of the paddle stroke; the
            paddle trajectory, without absorption (a paddle at rest
            reflects, as in the lab)
  solitary  the whole flume; the exact solitary wave of Babenko's equation
            with the height measured at gauge 4, its crest timed to the
            peak there (also NHWAVE's input, nhwave())

Both models: 3 sigma layers, log-law friction with roughness 0.1 mm
(NHWAVE's Zob in its examples of the inundation benchmarks), film /
minimum depth 5 mm; breaking: NHFOAM's hybrid model, NHWAVE's k-epsilon
model.  NHFOAM's three inputs run on dx = 1 cm (the grid of NHWAVE's
benchmark report).  On that grid NHWAVE grows short-wave noise in this 2DV
flume (crests 3-8% high, with or without its k-epsilon model; a flat tank
shows it too), so the two models are compared on dx = 2 cm, where NHWAVE
carries the solitary wave cleanly, both with the 'solitary' input.  Its
k-epsilon model (on in the template, as in NHWAVE's inundation examples)
also damps the steep solitary waves of cases B and C on the flat bed (a
flat tank: crest 0.87 H after 3.5 m at H/d = 0.26, 1.02 H without it), so
NHWAVE runs with and without it (kepsilon=False).

usage (python):
    import composite_beach as cb
    cb.nhfoam("nf_A_gauge4", "A", "gauge4")   # then ./Allrun
    cb.nhfoam("nf_A_solitary_02", "A", "solitary", dx=0.02, dt=0.004)
    cb.nhwave("nw_A", "A", dx=0.02, px=1)     # then mpirun -np 1 nhwave
    cb.nhwave("nw_A_inv", "A", dx=0.02, px=1, kepsilon=False)
"""
import os, io, math, zipfile
import numpy as np
import make_case as mc

G = 9.81
D = 0.218
XW = 23.23
XB = [0.0, 15.04, 19.40, 22.33, XW]          # breaks in slope
SLOPES = [0.0, 1 / 53, 1 / 150, 1 / 13]
GAUGE_X = {"g5": 15.04, "g6": 17.22, "g7": 19.40, "g8": 20.86,
           "g9": 22.33, "g10": 22.80}
DX, NZ, KS, FILM = 0.01, 3, 1e-4, 0.005
DATA = ("/home/claude/ext/nthmp/nthmp-benchmark-problems/"
        "BP05-ElenaT-Solitary_wave_on_composite_beach_laboratory/tsunami3_runup.zip")
NHWAVE_TEMPLATE = ("/home/claude/ext/nhw_forks/NHWAVE_NHWAVE/examples/inundation_benchmarks/"
                   "Solitary_wave_on_a_simple_beach_laboratory/H_d_0.30/input/input.txt")

CASES = {"A": dict(x4=12.64, ts="ts3a.txt", col=0),
         "B": dict(x4=14.06, ts="ts3b.txt", col=1),
         "C": dict(x4=14.40, ts="ts3c.txt", col=2)}
T_LAB0 = 265.0          # start of the gauge records
T_PULSE = 275.0         # end of the direct pulse at gauge 4 (benchmark)
T_PADDLE0 = 258.0       # start of the paddle records
T_END = 295.0


def depth(x):
    x = np.asarray(x, float)
    h = np.full_like(x, D)
    for i in range(1, len(XB) - 1):
        h -= SLOPES[i] * np.clip(x - XB[i], 0.0, XB[i + 1] - XB[i])
    return h


def _lines(member):
    with zipfile.ZipFile(DATA) as z:
        return io.TextIOWrapper(z.open(member), encoding="latin-1").read().splitlines()


def _numbers(lines, n):
    rows = []
    for l in lines:
        try:
            v = [float(x) for x in l.split()]
        except ValueError:
            continue
        if len(v) == n:
            rows.append(v)
    return np.array(rows)


def read_lab(name):
    """Lab time [s] and gauges g4..g10 [m], relative to their mean over the
    first 2 s."""
    a = _numbers(_lines(CASES[name]["ts"]), 8)
    t = a[:, 0]
    pre = t < t[0] + 2.0
    names = ["g4", "g5", "g6", "g7", "g8", "g9", "g10"]
    return t, {n: a[:, k + 1] - a[pre, k + 1].mean() for k, n in enumerate(names)}


def read_paddle(name):
    """Lab time [s] and paddle position [m] (from the rest position)."""
    a = _numbers(_lines("fdbk3abc.txt"), 4)
    return a[:, 0], a[:, 1 + CASES[name]["col"]] / 100.0


def read_runup():
    """Highest run-up on the wall [m] of cases A, B, C."""
    out = {}
    for l in _lines("run3abc.txt"):
        v = l.split()
        if len(v) == 5 and v[0] in CASES:
            out[v[0]] = float(v[3]) / 100.0
    return out


def incident(name):
    """Direct pulse of gauge 4 (to 275 s, the last 0.5 s tapered to zero)
    and its peak height and time."""
    t, g = read_lab(name)
    e = g["g4"].copy()
    e[t >= T_PULSE] = 0.0
    tap = (t > T_PULSE - 0.5) & (t < T_PULSE)
    e[tap] *= np.sin(0.5 * np.pi * (T_PULSE - t[tap]) / 0.5) ** 2
    i = np.argmax(e)
    return t, e, e[i], t[i]


def write_table(path, t, v):
    with open(path, "w") as f:
        f.write(f"{len(t)}\n(\n")
        f.writelines(f"({a:.6f} {b:.8e})\n" for a, b in zip(t, v))
        f.write(")\n")


def nhfoam(case, name, mode, dx=DX, nz=NZ, dt=0.002, ks=KS, film=FILM,
           breaking=True, conv="Gauss limitedLinearV 1", celerity=None):
    c = CASES[name]
    x4 = c["x4"]
    tg, eg, Hg, tpk = incident(name)
    if mode == "gauge4":
        x0, t0 = x4, T_LAB0
    elif mode == "paddle":
        tp, xp = read_paddle(name)
        x0, t0 = 0.5 * (xp.min() + xp.max()), T_PADDLE0
    elif mode == "solitary":
        from solitary_exact import SolitaryWave
        sw = SolitaryWave(Hg, d=D)
        gamma = math.sqrt(0.75 * Hg / D)
        xc0 = D * math.acosh(math.sqrt(1000.0)) / gamma        # back tail H/1000 at the paddle
        x0, t0 = 0.0, tpk - (x4 - xc0) / sw.c
    else:
        raise ValueError(mode)
    endTime = T_END - t0
    xs = [x0] + [x for x in XB[1:]]
    hs = list(depth(xs))
    gx = dict(g4=x4, **GAUGE_X, wall=XW - 0.5 * dx)
    gtxt = "\n".join(f"    {g}    ({x:.4f} {dx/2:g} 0);" for g, x in gx.items())
    brk = """
breaking
{
    model       hybrid;
    alpha       0.6;
    beta        0.3;
    betaFinal   0.15;
    Tstar       5;
    delta       1.2;
}
""" if breaking else ""
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film:g};
}}

bottomFriction
{{
    ks          {ks:g};             // log law, as NHWAVE's Zob
}}
{brk}
etaMax          yes;

gauges
{{
{gtxt}
}}
""")
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=max(1, int(round(5.0 / dt))), nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dx, dx, nz)
    mc.common_fields(case)
    mc.allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)
    if mode == "gauge4":
        sel = tg < T_PULSE + 1.0
        write_table(os.path.join(case, "constant", "gauge4.dat"), tg[sel], eg[sel])
        cel = ""
        if celerity == "pulse":
            cel = f"celerity    {math.sqrt(G * (D + Hg)):.6g};   // sqrt(g (d + H))"
        elif celerity:
            cel = f"celerity    {celerity:g};"
        text = f"""theory      timeSeries;
        depth       {D:g};
        direction   (1 0 0);
        origin      ({x4:g} 0 0);
        eta
        {{
            type    tableFile;
            file    "$FOAM_CASE/constant/gauge4.dat";
        }}
        {cel}
        timeShift   {T_LAB0:g};"""
        mc.set_patch(case, "U", "left", mc.wave_bc("left", text, absorb=True))
    elif mode == "paddle":
        write_table(os.path.join(case, "constant", "paddle.dat"), tp, xp * 100.0)
        text = f"""theory      timeSeries;
        depth       {D:g};
        direction   (1 0 0);
        origin      ({x0:.4f} 0 0);
        paddle
        {{
            type    tableFile;
            file    "$FOAM_CASE/constant/paddle.dat";   // [cm]
        }}
        scale       0.01;
        timeShift   {T_PADDLE0:g};"""
        mc.set_patch(case, "U", "left", mc.wave_bc("left", text, absorb=False))
    else:
        _solitary_fields(case, xs, hs, dx, nz, sw, xc0)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"name {name}\nmode {mode}\nt0 {t0:.10g}\nx0 {x0:.10g}\n"
                f"dx {dx:g}\nnz {nz}\nH4 {Hg:.6g}\ntpeak4 {tpk:.6g}\n")
    return dict(x0=x0, t0=t0, endTime=endTime, H4=Hg, nsteps=nsteps)


def _solitary_fields(case, xs, hs, dx, nz, sw, xc0):
    """Exact solitary wave on the cell centres; cells of blockMesh block by
    block, x fastest, then the layers.  The wave lies in the first (flat)
    block; the other blocks start at rest."""
    import ma2012
    eta, U, W, P = [], [], [], []
    for b in range(len(xs) - 1):
        nx = max(1, int(round((xs[b + 1] - xs[b]) / dx)))
        xc = xs[b] + (np.arange(nx) + 0.5) * (xs[b + 1] - xs[b]) / nx
        for k in range(nz):
            if b == 0:
                e = sw.eta(xc - xc0)
                zc = -D + (k + 0.5) / nz * (D + e)
                u, w, p = ma2012.solitary_fields(sw, xc - xc0, zc, D)
            else:
                e = u = w = p = np.zeros(nx)
            eta.append(e); U.append(u); W.append(w); P.append(p)
    eta = np.concatenate(eta); U = np.concatenate(U); W = np.concatenate(W)
    ma2012._set_internal(os.path.join(case, "0", "eta"), ma2012._foam_list(eta))
    ma2012._set_internal(os.path.join(case, "0", "U"),
                         ma2012._foam_list(list(zip(U, W)), vector=True))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"),
                         ma2012._foam_list(np.concatenate(P)))


def nhwave(case, name, template=NHWAVE_TEMPLATE, dx=DX, nz=NZ, ks=KS, film=FILM,
           px=2, cfl=0.5, kepsilon=True):
    """NHWAVE on the whole flume with the exact solitary wave of the
    'solitary' mode as its initial condition (no wavemaker), walls at both
    ends, with its k-epsilon model (the template's settings) or without it
    (kepsilon=False: VISCOUS_FLOW F); gauges as NHFOAM's."""
    from solitary_exact import SolitaryWave
    import other_models as om
    c = CASES[name]
    x4 = c["x4"]
    _, _, Hg, tpk = incident(name)
    sw = SolitaryWave(Hg, d=D)
    gamma = math.sqrt(0.75 * Hg / D)
    xc0 = D * math.acosh(math.sqrt(1000.0)) / gamma
    t0 = tpk - (x4 - xc0) / sw.c
    M = int(round(XW / dx))
    x = (np.arange(M) + 0.5) * dx
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    np.savetxt(os.path.join(case, "depth.txt"), depth(x)[None, :], fmt="%.6f")
    gx = dict(g4=x4, **GAUGE_X, wall=XW - 0.5 * dx)
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for g, xg in gx.items():
            f.write(f"{xg:.4f} {0.5 * dx:.4f} 0.0\n")
    txt = open(template, encoding="latin-1").read()
    txt = om.set_keys(txt, {
        "Mglob": M, "Nglob": 1, "Kglob": nz, "PX": px, "PY": 1,
        "DX": dx, "DY": dx, "TOTAL_TIME": T_END - t0, "PLOT_START": 0.0,
        "PLOT_INTV": T_END - t0, "SCREEN_INTV": 1.0, "PLOT_INTV_STAT": 0.01,
        "NSTAT": len(gx), "INITIAL_EUVW": "T", "WAVEMAKER": "LEF_SOL",
        "AMP": 0.0, "BC_X0": 1, "BC_Xn": 1, "Ibot": 2, "Zob": ks,
        "MinDep": film, "CFL": cfl, "SPONGE_ON": "F",
        "VISCOUS_FLOW": "T" if kepsilon else "F",
        "Sponge_West_Width": 0.0, "Sponge_East_Width": 0.0,
        "OUT_H": "F", "OUT_E": "T", "OUT_M": "T", "OUT_U": "F", "OUT_W": "F"})
    open(os.path.join(case, "input.txt"), "w", encoding="latin-1").write(txt)
    eta = sw.eta(x - xc0)
    U = np.zeros((nz, M)); W = np.zeros((nz, M))
    for k in range(nz):
        zk = -depth(x) + (k + 0.5) / nz * (depth(x) + eta)
        flat = x < XB[1]
        U[k][flat], W[k][flat] = sw.velocity((x - xc0)[flat], zk[flat])
    eta = np.where(x < XB[1], eta, 0.0)
    np.savetxt(os.path.join(case, "eta0.txt"), eta[None, :], fmt="%.10e")
    np.savetxt(os.path.join(case, "uvw0.txt"), np.vstack([U, np.zeros_like(U), W]),
               fmt="%.10e")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"name {name}\nmode nhwave\nt0 {t0:.10g}\nx0 0\n"
                f"dx {dx:g}\nnz {nz}\nH4 {Hg:.6g}\ntpeak4 {tpk:.6g}\n"
                f"closure {'kepsilon' if kepsilon else 'none'}\n")
    return dict(t0=t0, M=M, H4=Hg)


# --------------------------------------------------------------------------
# Results

GAUGE_NAMES = ["g4", "g5", "g6", "g7", "g8", "g9", "g10", "wall"]


def params(case):
    return {l.split()[0]: l.split()[1] for l in open(os.path.join(case, "params.txt"))}


def read_model(case):
    """Lab time [s] and gauges of a NHFOAM or NHWAVE run."""
    p = params(case)
    t0 = float(p["t0"])
    if p["mode"] == "nhwave":
        res = os.path.join(case, "results")
        files = os.listdir(res)
        out, t = {}, None
        for k, g in enumerate(GAUGE_NAMES):
            stem = f"probe_{k + 1:04d}"
            rows = [np.loadtxt(os.path.join(res, f), ndmin=2)[:, :2]
                    for f in files if f.startswith(stem)
                    and os.path.getsize(os.path.join(res, f)) > 0]
            rows = [r for r in rows if r.size]
            if not rows:
                # NHWAVE leaves a probe in the last half cell empty (the
                # wall gauge): its run-up comes from hmax
                continue
            d = np.vstack(rows)
            d = d[np.argsort(d[:, 0], kind="stable")]
            keep = np.r_[True, np.diff(d[:, 0]) > 0]
            t, out[g] = d[keep, 0], d[keep, 1]
        return t + t0, out
    import conical_island as ci
    t, g = ci.read_nhfoam_gauges(case)
    return t + t0, g


def wall_runup(case):
    """Highest free surface in the column at the wall [m]: NHFOAM's etaMax
    (the last cell of blockMesh's last block is in that column), NHWAVE's
    hmax (which stays zero in its last cell: the one before it)."""
    p = params(case)
    if p["mode"] == "nhwave":
        res = os.path.join(case, "results")
        f = sorted(x for x in os.listdir(res) if x.startswith("hmax_"))[-1]
        return float(np.atleast_2d(np.loadtxt(os.path.join(res, f)))[0, -2:].max())
    import conical_island as ci
    times = [d for d in os.listdir(case) if d.replace(".", "", 1).isdigit() and d != "0"]
    v = ci._field_values(os.path.join(case, max(times, key=float), "etaMax"))
    return float(v[-1])


def peaks(case, name, t1=T_END):
    """Highest crest at each gauge, model / lab - 1, for the incident
    wave (before the reflection comes back) and over the whole record."""
    tl, gl = read_lab(name)
    tm, gm = read_model(case)
    out = {}
    for g in GAUGE_NAMES[:-1]:
        lab = gl[g][tl <= t1].max()
        mod = gm[g][tm <= t1].max()
        out[g] = (mod / lab - 1.0, mod, lab)
    return out


def crests(case, name):
    """Incident and reflected crests, model / lab - 1, at gauges 4-8
    (split at the time of the measured crest at gauge 10, which the
    incident crest passes before and the reflected one after), and the
    highest crest at gauges 9 and 10, where the two overlap."""
    tl, gl = read_lab(name)
    tm, gm = read_model(case)
    ts = tl[np.argmax(gl["g10"])]
    out = {}
    for g in GAUGE_NAMES[:-1]:
        if g in ("g9", "g10"):
            parts = [("max", tl <= T_END, tm <= T_END)]
        else:
            parts = [("inc", tl < ts, tm < ts), ("ref", tl >= ts, tm >= ts)]
        for k, sl, sm in parts:
            lab = gl[g][sl].max()
            mod = gm[g][sm].max() if sm.any() else np.nan     # unfinished run
            out[f"{g}_{k}"] = (mod / lab - 1.0, mod, lab)
    return out


def arrivals(case, name, gauges=("g5", "g7", "g8")):
    """Arrival of the incident crest, model - lab [s], at the gauges (the
    split of crests())."""
    tl, gl = read_lab(name)
    tm, gm = read_model(case)
    ts = tl[np.argmax(gl["g10"])]
    out = {}
    for g in gauges:
        sl, sm = tl < ts, tm < ts
        out[g] = float(tm[sm][np.argmax(gm[g][sm])] - tl[sl][np.argmax(gl[g][sl])])
    return out


def rms(case, name, t1=T_END):
    """rms difference model - lab [m] over gauges 5-10 on the lab times
    from the arrival of the wave at gauge 5 - 2 s to t1."""
    tl, gl = read_lab(name)
    tm, gm = read_model(case)
    ta = tl[np.argmax(gl["g5"])] - 2.0
    sel = (tl >= max(ta, tm[0])) & (tl <= min(t1, tm[-1]))
    err = [np.interp(tl[sel], tm, gm[g]) - gl[g][sel] for g in GAUGE_NAMES[1:-1]]
    return float(np.sqrt(np.mean(np.concatenate(err) ** 2)))


if __name__ == "__main__":
    import sys
    lab = read_runup()
    for case in sys.argv[1:]:
        p = params(case)
        name = p["name"]
        pk = peaks(case, name)
        print(f"{case:16s} {p['mode']:8s} rms {100 * rms(case, name):.2f} cm  wall {100 * wall_runup(case):.2f} cm"
              f" (lab {100 * lab[name]:.2f})  "
              + " ".join(f"{g}:{100 * v[0]:+.0f}%" for g, v in pk.items()))
