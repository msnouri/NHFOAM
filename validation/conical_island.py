#!/usr/bin/env python3
"""Solitary waves on a conical island (Briggs et al. 1995; NTHMP benchmark
problem 6), NHFOAM and NHWAVE on the same grid.

Lab: still-water depth d = 0.32 m; the island is a truncated cone with its
toe radius 3.6 m, crest radius 1.1 m and 1:4 faces, centred at
(12.96, 13.80) m.  Incident gauges 1-4 on x = 5.76 / 6.82 / 7.56 m for cases
A / B / C, gauges 6 and 9 in front of the island, 16 on its side, 22 behind
it.  Data (NTHMP benchmark files): ts2a.txt, ts2b.txt, ts2cnew1.txt (gauges)
and run2a.txt, run2b.txt, run2c.txt (run-up around the island).

Both models use the grid of NHWAVE's own example
(inundation_benchmarks/Solitary_wave_on_a_conical_island): dx = dy = 5 cm,
3 uniform sigma layers, quadratic bottom friction Cd = 0.01.  To save time:

* only the half y <= 13.80 m is computed: the island centre lies on the
  symmetry line, a slip wall in both models (gauges 1 and 2 are the mirror
  images of 4 and 3);
* the wave starts as the exact solitary wave of Babenko's equation with its
  crest L = d acosh(sqrt 20)/gamma upstream of the incident gauges
  (eta = H/20 there), gamma = sqrt(3H/4d), as in the run-up test (6.6),
  instead of NHWAVE's boundary wavemaker; the basin starts far enough
  upstream to hold the tail of the wave down to H/200;
* the basin ends at x = 22 m with a 2 m sponge (NHWAVE's example: 2 m).

Each model keeps its own breaking closure (NHFOAM: hybrid breaking model;
NHWAVE: the k-epsilon model of its example) and the same film/minimum depth.
Run-up: the highest free surface of each column while wet (NHFOAM etaMax,
NHWAVE hmax) on the edge of the flooded area, in 5-degree sectors around the
island (NHWAVE's script for this benchmark uses the same edge cells).

usage (python):
    import conical_island as ci
    ci.nhfoam("nf_A", "A")              # then ./Allrun (2 processors)
    ci.nhwave("nw_A", "A", template)    # then mpirun -np 2 nhwave
"""
import os, re, math, subprocess
import numpy as np
import make_case as mc

G = 9.81
D = 0.32                         # still-water depth [m]
XC, YC = 12.96, 13.80            # island centre [m]
R_TOE, R_TOP, SLOPE = 3.6, 1.1, 0.25
Z_TOP = -D + (R_TOE - R_TOP) * SLOPE     # crest of the island, 0.305 m
X1, SPONGE = 22.0, 2.0
DX, NZ = 0.05, 3
FILM, CD = 0.005, 0.01
DATA = "/home/claude/ext/nthmp/nthmp-benchmark-problems/BP06-FrankG-Solitary_wave_on_a_conical_island"

# H/d as measured at the incident gauges (NHWAVE's example: 0.045, 0.091,
# 0.181); x of gauges 1-4; west end of the basin
CASES = {
    "A": dict(Hd=0.045, xg=5.76, x0=-4.0, ts="ts2a.txt", ru="run2a.txt"),
    "B": dict(Hd=0.091, xg=6.82, x0=0.0, ts="ts2b.txt", ru="run2b.txt"),
    "C": dict(Hd=0.181, xg=7.56, x0=2.4, ts="ts2cnew1.txt", ru="run2c.txt"),
}

# lab gauges (x, y); gauges 1-4 at x = xg
GAUGES_LAB = {"g1": (None, 16.05), "g2": (None, 14.55), "g3": (None, 13.05),
              "g4": (None, 11.55), "g6": (9.36, 13.80), "g9": (10.36, 13.80),
              "g16": (12.96, 11.22), "g22": (15.56, 13.80)}
MODEL_GAUGES = ["g4", "g3", "g6", "g9", "g16", "g22"]   # g1 = g4, g2 = g3


def zb(x, y):
    """Bed elevation [m] (still water level z = 0)."""
    r = np.hypot(np.asarray(x) - XC, np.asarray(y) - YC)
    return -D + np.clip((R_TOE - r) * SLOPE, 0.0, (R_TOE - R_TOP) * SLOPE)


def crest_start(Hd, xg, d=D):
    """Crest position at t = 0: eta = H/20 at the incident gauges."""
    gamma = math.sqrt(0.75 * Hd)
    return xg - d * math.acosh(math.sqrt(20.0)) / gamma


def grid(case_name, dx=DX):
    c = CASES[case_name]
    x0 = c["x0"]
    nx = int(round((X1 - x0) / dx))
    ny = int(round(YC / dx))
    return x0, nx, ny


def gauge_cells(case_name, dx=DX):
    """Model gauges: the cell containing the (mirrored) lab position."""
    x0, nx, ny = grid(case_name, dx)
    out = {}
    for g in MODEL_GAUGES:
        x, y = GAUGES_LAB[g]
        if x is None:
            x = CASES[case_name]["xg"]
        y = min(y, 2 * YC - y)              # mirror into y <= YC
        i = min(int(math.floor((x - x0) / dx)), nx - 1)
        j = min(int(math.floor(y / dx)), ny - 1)
        out[g] = (i, j, x0 + (i + 0.5) * dx, (j + 0.5) * dx)
    return out


def end_time(case_name):
    """Crest at the incident gauges after about 2 s (A), 1.4 s (B), 0.9 s
    (C); the lab records show the main peaks within 8 s after that."""
    return {"A": 20.0, "B": 18.0, "C": 16.0}[case_name]


# --------------------------------------------------------------------------
# NHFOAM

def nhfoam(case, name, dx=DX, nz=NZ, dt=0.01, endTime=None, film=FILM, Cd=CD,
           nproc=2, conv="Gauss limitedLinearV 1", breaking=True,
           writeEvery=4.0):
    from solitary_exact import SolitaryWave
    import ma2012
    c = CASES[name]
    x0, nx, ny = grid(name, dx)
    endTime = endTime or end_time(name)
    H = c["Hd"] * D
    xc0 = crest_start(c["Hd"], c["xg"])
    mc.write(case, "system/blockMeshDict", "dictionary", f"""
convertToMeters 1;
vertices
(
    ({x0} 0 {-D}) ({X1} 0 {-D}) ({X1} {YC} {-D}) ({x0} {YC} {-D})
    ({x0} 0 0) ({X1} 0 0) ({X1} {YC} 0) ({x0} {YC} 0)
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );
edges ();
boundary
(
    surface {{ type patch; faces ((4 5 6 7)); }}
    bottom  {{ type wall;  faces ((0 3 2 1)); }}
    left    {{ type wall;  faces ((0 4 7 3)); }}
    right   {{ type wall;  faces ((1 2 6 5)); }}
    south   {{ type wall;  faces ((0 1 5 4)); }}
    north   {{ type wall;  faces ((3 7 6 2)); }}
);
mergePatchPairs ();
""")
    gtxt = "\n".join(f"    {g}    ({gx:.4f} {gy:.4f} 0);"
                     for g, (i, j, gx, gy) in gauge_cells(name, dx).items())
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
    Cd          {Cd:g};
}}
{brk}
etaMax          yes;

gauges
{{
{gtxt}
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({X1 - SPONGE:g} 0 0);
        end         ({X1:g} 0 0);
        sigmaMax    20;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(writeEvery / dt)), nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    mc.common_fields(case, walls=("left", "right", "south", "north"), empty=False)
    if nproc > 1:
        mc.write(case, "system/decomposeParDict", "dictionary", f"""numberOfSubdomains {nproc};
method          structured;     // whole columns on one processor
structuredCoeffs
{{
    method      simple;
    patches     (bottom);
    simpleCoeffs
    {{
        n       ({nproc} 1 1);
    }}
}}
""")
        mc.allrun(case, """
runApplication decomposePar
runParallel NHFOAM
runApplication reconstructPar
""")
    else:
        mc.allrun(case, "runApplication NHFOAM\n")

    # mesh: flat box, then the points onto the bed (dry columns: the film)
    subprocess.run(["blockMesh", "-case", case], check=True,
                   stdout=open(os.path.join(case, "log.blockMesh"), "w"))
    pf = os.path.join(case, "constant/polyMesh/points")
    txt = open(pf).read()
    head, rest = txt.split("(", 1)
    n = int(head.strip().split()[-1])
    body, tail = rest.rsplit(")", 1)
    P = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", body), float)
    assert len(P) == n
    zbp = zb(P[:, 0], P[:, 1])
    ztop = np.maximum(0.0, zbp + film)
    sig = (P[:, 2] + D) / D
    P[:, 2] = zbp + sig * (ztop - zbp)
    with open(pf, "w") as f:
        f.write(head + "(\n")
        f.writelines(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n" for p in P)
        f.write(")" + tail)

    # initial fields: exact solitary wave on the cell centres (blockMesh
    # order: x fastest, then y, then the layers)
    sw = SolitaryWave(H, d=D)
    xi = x0 + np.arange(nx + 1) * dx
    yj = np.arange(ny + 1) * dx
    XP, YP = np.meshgrid(xi, yj)                 # (ny+1, nx+1)
    zp = zb(XP, YP)
    zbc = 0.25 * (zp[:-1, :-1] + zp[:-1, 1:] + zp[1:, :-1] + zp[1:, 1:])  # (ny, nx)
    xc = x0 + (np.arange(nx) + 0.5) * dx
    eta1 = sw.eta(xc - xc0)
    eta2 = np.broadcast_to(eta1, (ny, nx))
    wet = eta2 - zbc > film
    eta2 = np.where(wet, eta2, zbc + film)
    # dry columns: below the bed, so that NHFOAM puts the film on its own bed
    # (the centre of a warped bottom face is not the mean of its corners)
    eta_out = np.where(wet, eta2, -10.0)
    etaC, uC, wC, pC = [], [], [], []
    for k in range(nz):
        Dk = eta2 - zbc
        zk = zbc + (k + 0.5) / nz * Dk
        X2 = np.broadcast_to(xc, (ny, nx))
        u, w = sw.velocity((X2 - xc0)[wet], zk[wet])
        U = np.zeros((ny, nx)); W = np.zeros((ny, nx))
        U[wet], W[wet] = u, w
        p = np.where(wet, 0.5 * (sw.c ** 2 - (U - sw.c) ** 2 - W ** 2), G * eta2)
        etaC.append(eta_out.ravel()); uC.append(U.ravel()); wC.append(W.ravel())
        pC.append(p.ravel())
    eta_all = np.concatenate(etaC)
    u_all = np.concatenate(uC); w_all = np.concatenate(wC)
    ma2012._set_internal(os.path.join(case, "0", "eta"), ma2012._foam_list(eta_all))
    ma2012._set_internal(os.path.join(case, "0", "U"),
                         ma2012._foam_list(list(zip(u_all, w_all)), vector=True))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"),
                         ma2012._foam_list(np.concatenate(pC)))
    return dict(x0=x0, nx=nx, ny=ny, xc0=xc0, c=sw.c, H=H, nsteps=nsteps)


# --------------------------------------------------------------------------
# NHWAVE

NHWAVE_EXAMPLE = ("/home/claude/ext/nhw_forks/NHWAVE_NHWAVE/examples/inundation_benchmarks/"
                  "Solitary_wave_on_a_conical_island/1_Hd0.045/input/input.txt")


def nhwave(case, name, template=NHWAVE_EXAMPLE, dx=DX, nz=NZ, endTime=None,
           film=FILM, Cd=CD, px=2, cfl=0.5):
    """NHWAVE's example on the same half basin, exact initial condition
    instead of its boundary wavemaker, output every 4 s (eta, hmax)."""
    from solitary_exact import SolitaryWave
    import other_models as om
    c = CASES[name]
    x0, nx, ny = grid(name, dx)
    endTime = endTime or end_time(name)
    H = c["Hd"] * D
    xc0 = crest_start(c["Hd"], c["xg"])
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    xc = x0 + (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xc, yc)                   # (ny, nx)
    z = zb(X, Y)
    np.savetxt(os.path.join(case, "depth.txt"), -z, fmt="%.6f")
    gc = gauge_cells(name, dx)
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for g, (i, j, gx, gy) in gc.items():
            f.write(f"{gx - x0:.4f} {gy:.4f} 0.0\n")
    txt = open(template, encoding="latin-1").read()
    txt = om.set_keys(txt, {
        "Mglob": nx, "Nglob": ny, "Kglob": nz, "PX": px, "PY": 1,
        "DX": dx, "DY": dx, "TOTAL_TIME": endTime, "PLOT_START": 0.0,
        "PLOT_INTV": 4.0, "SCREEN_INTV": 1.0, "PLOT_INTV_STAT": 0.02,
        "NSTAT": len(gc), "INITIAL_EUVW": "T", "WAVEMAKER": "LEF_SOL",
        "AMP": 0.0, "BC_X0": 1, "Cd0": Cd, "MinDep": film, "CFL": cfl,
        "Sponge_West_Width": 0.0, "Sponge_East_Width": SPONGE,
        "OUT_H": "F", "OUT_E": "T", "OUT_M": "T"})
    open(os.path.join(case, "input.txt"), "w", encoding="latin-1").write(txt)
    sw = SolitaryWave(H, d=D)
    eta = np.broadcast_to(sw.eta(xc - xc0), (ny, nx)).copy()
    wet = eta - z > 0
    eta = np.where(wet, eta, z)
    U = np.zeros((nz, ny, nx)); W = np.zeros((nz, ny, nx))
    for k in range(nz):
        zk = z + (k + 0.5) / nz * np.maximum(eta - z, 0.0)
        u, w = sw.velocity((X - xc0)[wet], zk[wet])
        U[k][wet], W[k][wet] = u, w
    np.savetxt(os.path.join(case, "eta0.txt"), eta, fmt="%.10e")
    with open(os.path.join(case, "uvw0.txt"), "w") as f:
        for A in (U, np.zeros_like(U), W):
            for k in range(nz):
                np.savetxt(f, A[k], fmt="%.10e")
    return dict(x0=x0, nx=nx, ny=ny, xc0=xc0, c=sw.c, H=H)


# --------------------------------------------------------------------------
# Results

def read_lab(name):
    rows = []
    for l in open(os.path.join(DATA, CASES[name]["ts"])):
        try:
            v = [float(x) for x in l.split()]
        except ValueError:
            continue
        if len(v) == 9:
            rows.append(v)
    d = np.array(rows)
    names = ["g1", "g2", "g3", "g4", "g6", "g9", "g16", "g22"]
    return d[:, 0], {n: d[:, k + 1] for k, n in enumerate(names)}


def read_lab_runup(name):
    out = []
    for l in open(os.path.join(DATA, CASES[name]["ru"])):
        v = l.split()
        if len(v) == 4:
            try:
                out.append([float(x) for x in v])
            except ValueError:
                pass
    a = np.array(out)
    return a[:, 1], a[:, 2] / 100.0          # degrees, run-up [m]


def read_nhfoam_gauges(case):
    """Gauge records; a restarted run writes postProcessing/gauges/<t0>, and
    its records replace those of the earlier runs from t0 on."""
    root = os.path.join(case, "postProcessing", "gauges")
    starts = sorted((d for d in os.listdir(root)), key=float)
    parts, names = [], None
    for k, t0 in enumerate(starts):
        fn = os.path.join(root, t0, "eta.dat")
        d = np.loadtxt(fn, comments="#", ndmin=2)
        if names is None:
            names = [l.split()[1] for l in open(fn) if l.startswith("# ") and "(" in l]
        if k + 1 < len(starts):
            d = d[d[:, 0] <= float(starts[k + 1]) + 1e-9]
        parts.append(d)
    d = np.vstack(parts)
    keep = np.r_[True, np.diff(d[:, 0]) > 0]
    d = d[keep]
    return d[:, 0], {n: d[:, k + 1] for k, n in enumerate(names)}


def read_nhwave_gauges(case, name, dx=DX):
    """NHWAVE's probes builds its file names from an 80-character buffer of
    which only 4 characters are set, so the records of one probe can end up
    in several files (probe_0001 followed by stray bytes): all are read and
    merged by time."""
    res = os.path.join(case, "results")
    files = os.listdir(res)
    out, t = {}, None
    for k, g in enumerate(gauge_cells(name, dx)):
        stem = f"probe_{k + 1:04d}"
        rows = []
        for f in files:
            fn = os.path.join(res, f)
            if f.startswith(stem) and os.path.getsize(fn) > 0:
                a = np.loadtxt(fn, ndmin=2)
                if a.size:
                    rows.append(a[:, :2])
        d = np.vstack(rows)
        d = d[np.argsort(d[:, 0], kind="stable")]
        keep = np.r_[True, np.diff(d[:, 0]) > 0]
        t = d[keep, 0]
        out[g] = d[keep, 1]
    return t, out


def _field_values(path):
    s = open(path).read()
    i = s.index("internalField")
    j = s.index(";", i)
    body = s[i:j]
    if "uniform" in body.split("\n")[0] and "nonuniform" not in body:
        raise ValueError("uniform field")
    nums = body.split("(", 1)[1].rsplit(")", 1)[0].split()
    return np.array(nums, float)


def nhfoam_etamax(case, name, time=None, dx=DX):
    """etaMax per column (ny, nx) of the last (or the given) time."""
    x0, nx, ny = grid(name, dx)
    if time is None:
        times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"]
        time = max(times, key=float)
    v = _field_values(os.path.join(case, time, "etaMax"))
    return v[: nx * ny].reshape(ny, nx), time


def nhwave_hmax(case, name, dx=DX):
    files = sorted(f for f in os.listdir(os.path.join(case, "results"))
                   if re.fullmatch(r"hmax_[0-9]{5}", f))
    a = np.loadtxt(os.path.join(case, "results", files[-1]))
    return a, files[-1]


def shoreline(emax, name, flooded_fn, dx=DX):
    """Cells on the edge of the flooded area (wet at some time, next to a
    column that never was), as NHWAVE's own script for this benchmark does
    with hmax; returns the mask, the angle of each cell around the island
    and emax.  Angles as in the lab files: 0 deg towards -y (gauge 16),
    90 deg behind the island (gauge 22), 270 deg in front (gauges 6, 9)."""
    x0, nx, ny = grid(name, dx)
    xc = x0 + (np.arange(nx) + 0.5) * dx
    yc = (np.arange(ny) + 0.5) * dx
    X, Y = np.meshgrid(xc, yc)
    z = zb(X, Y)
    fl = (z < 0.0) | flooded_fn(emax, z)
    dry_nb = np.zeros_like(fl)
    dry_nb[:, 1:] |= ~fl[:, :-1]
    dry_nb[:, :-1] |= ~fl[:, 1:]
    dry_nb[1:, :] |= ~fl[:-1, :]
    dry_nb[:-1, :] |= ~fl[1:, :]
    th = np.degrees(np.arctan2(X - XC, -(Y - YC))) % 360.0
    return fl & dry_nb, th, emax


def runup(emax, name, flooded_fn, angles, dx=DX, width=5.0):
    """Run-up at the given angles: the highest free surface of the edge
    cells within +-width/2.  The computed half covers 270..360 and 0..90
    deg; other angles are mirrored (theta -> 180 - theta)."""
    edge, th, e = shoreline(emax, name, flooded_fn, dx)
    out = []
    for a in np.atleast_1d(angles):
        a = a % 360.0
        if 90.0 < a < 270.0:
            a = (180.0 - a) % 360.0
        d = (th - a + 180.0) % 360.0 - 180.0
        sel = edge & (np.abs(d) <= 0.5 * width)
        out.append(e[sel].max() if sel.any() else np.nan)
    return np.array(out)


def flooded_nhfoam(emax, z, film=FILM):
    """Columns wet at some time: eta > z_b + film (never-wet columns keep
    their bed in etaMax)."""
    return emax > z + 0.5 * film


def flooded_nhwave(hmax, z):
    """NHWAVE's hmax starts at zero and grows only while the cell is wet;
    cells below the still-water level count as flooded (shoreline())."""
    return hmax > 0.0


# --------------------------------------------------------------------------
# All runs of the README/report (section 6.9):
#     python3 conical_island.py <dir>     (then <dir>/runall.sh, about 6 h on
#                                          2 cores, NHWAVE 4.5 h of it)
#     plot_conical_island.py with the runs listed in validation/make_figures.sh

if __name__ == "__main__":
    import sys
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    os.chdir(out)
    for name in ("A", "B", "C"):
        nhfoam(f"nf_{name}_05", name, dx=0.05, dt=0.01)
        nhfoam(f"nf_{name}_10", name, dx=0.1, dt=0.02)
        nhwave(f"nw_{name}_10", name, dx=0.1)
    with open("runall.sh", "w") as f:
        f.write("#!/bin/sh\n# NHFOAM (Allrun: 2 processors) and NHWAVE (2 processors)\n"
                "NHWAVE=${NHWAVE:-nhwave}\n")
        for name in ("A", "B", "C"):
            for c in (f"nf_{name}_10", f"nf_{name}_05"):
                f.write(f"(cd {c} && ./Allrun)\n")
            f.write(f"(cd nw_{name}_10 && mpirun -np 2 $NHWAVE > log.nhwave 2>&1)\n")
    os.chmod("runall.sh", 0o755)
