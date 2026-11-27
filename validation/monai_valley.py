#!/usr/bin/env python3
"""Tsunami run-up onto a complex three-dimensional beach: the Monai valley
(Matsuyama & Tanaka 2001; Liu et al. 2008; NTHMP benchmark problem 7),
NHFOAM and NHWAVE.

Lab: a 1:400 model of the coast near Monai (Okushiri Island, 1993
Hokkaido Nansei-Oki tsunami) in the large tank of CRIEPI.  The benchmark
covers 5.488 x 3.402 m with the bathymetry and topography on a 1.4 cm grid
(393 x 244 nodes, depth 13.535 cm at x = 0, land up to 12.5 cm above the
still water), walls at y = 0, y = 3.402 m and x = 5.488 m.  The incident
wave is given at x = 0 for 22.5 s (a depression of -2.9 mm, then a crest
of 16.1 mm, H/d = 0.12); after that the boundary should let the waves
out.  Gauges 5, 7 and 9 at x = 4.521 m, y = 1.196, 1.696 and 2.196 m (in
1.2, 0.3 and 0.6 cm of water by the bathymetry file); in six repeats of
the experiment the run-up in the narrow valley at (5.158, 1.880) m was
8.0-10.0 cm, at two points on its north side (5.030, 2.206) and
(4.998, 2.320) m 5.25-6.75 and 5.25-5.75 cm.  Data
(the NTHMP benchmark files, read from DATA): Benchmark_2_Bathymetry.txt,
Benchmark_2_input.txt, output_ch5-7-9.xls, OBS_RUNUP.txt.

NHFOAM: the incident wave at x = 0 as a time series (theory timeSeries)
with the active absorption, so that the boundary lets the reflected waves
out, also after 22.5 s (the record tapered to zero over its last 0.5 s);
depth-uniform velocity u = c eta/(d + eta), c = sqrt(g (d + eta)).

NHWAVE: its own example (inundation_benchmarks/Runup_onto_a_complex_3-D_
beach_Monai_Valley): the domain extended 1.4 m offshore at d = 0.135 m,
the internal wavemaker INT_IRR between x = -0.14 and 0 with a Fourier fit
of the input wave (irr2d.txt, 399 components), a 1 m sponge at the west
end, 4 sigma layers, the k-epsilon model (RNG), log-law friction with
Zob = 0.1 mm and MinDep = 2 mm.  Its wave at x = 0 keeps time with the
record (lag() over the first 12 s: 0.001 s), so no run is shifted in time
(NHWAVE's example shifts its results by -0.4 s for its comparisons).

Both models: 4 sigma layers, log-law friction with roughness 0.1 mm, film /
minimum depth 2 mm (NHWAVE's example); breaking: NHFOAM's hybrid model,
NHWAVE's k-epsilon model.  NHFOAM runs on the benchmark grid (dx = 1.4 cm,
392 x 243 cells) and on the grid coarsened twice (dx = 2.8 cm,
dy = 3.402/122 = 2.79 cm).  NHWAVE takes about 10 s per time step on the
benchmark grid on two cores (about 16 h for 25 s), so it runs on the coarse grid
only: bed by bilinear interpolation of the benchmark nodes; NHWAVE's cell
depths are the means of the four corners, which is how its example's
depth.txt is made from the benchmark nodes on the fine grid.

Results: gauges (NHFOAM gauges, NHWAVE probes, at the cell centres around
the lab positions), the highest surface while wet (NHFOAM etaMax, NHWAVE
hmax) and the run-up (runup(): the highest flooded ground, which is how the
lab run-up was taken: the ground at the three lab points is 8.9, 6.2 and
5.5 cm high), the free surface at the times of the video frames (NHFOAM
surfaceSnapshots every 0.1 s, NHWAVE's eta output every 0.1 s).

usage (python):
    import monai_valley as mv
    mv.nhfoam("nf_14", dx=0.014, dt=0.004)   # then ./Allrun (2 processors)
    mv.nhfoam("nf_28", dx=0.028, dt=0.008)
    mv.nhwave("nw_28", dx=0.028)             # then mpirun -np 2 nhwave
    mv.nhwave("nw_28_inv", dx=0.028, kepsilon=False)
"""
import os, re, math, subprocess
import numpy as np
import make_case as mc

G = 9.81
D = 0.135                      # depth of the first column and of NHWAVE's extension
LX, LY, DX0 = 5.488, 3.402, 0.014
NZ, FILM, KS = 4, 0.002, 1e-4
T_REC = 22.5                   # end of the input record
T_END = 25.0
FRAMES = [15.3, 15.8, 16.3, 16.8, 17.3]   # video frames 10, 25, 40, 55, 70
GAUGES = {"x0": (None, 1.701), "ch5": (4.521, 1.196), "ch7": (4.521, 1.696),
          "ch9": (4.521, 2.196)}
DATA = ("/home/claude/ext/nthmp/nthmp-benchmark-problems/"
        "BP07-DmitryN-Monai_valley_beach")
NHWAVE_EXAMPLE = ("/home/claude/ext/nhw_forks/NHWAVE_NHWAVE/examples/inundation_benchmarks/"
                  "Runup_onto_a_complex_3-D_beach_Monai_Valley/input")
X_EXT = 1.4                    # NHWAVE's offshore extension [m]


# --------------------------------------------------------------------------
# Benchmark data

_BATHY = None


def bathymetry():
    """Depth at the benchmark nodes, Z[j, i] at x = 0.014 i, y = 0.014 j
    (positive below the still water level)."""
    global _BATHY
    if _BATHY is None:
        d = np.loadtxt(os.path.join(DATA, "Benchmark_2_Bathymetry.txt"), skiprows=1)
        i = np.rint(d[:, 0] / DX0).astype(int)
        j = np.rint(d[:, 1] / DX0).astype(int)
        Z = np.full((j.max() + 1, i.max() + 1), np.nan)
        Z[j, i] = d[:, 2]
        assert not np.isnan(Z).any()
        _BATHY = Z
    return _BATHY


def zb(x, y):
    """Bed elevation [m] (still water z = 0), bilinear between the nodes."""
    Z = bathymetry()
    x = np.clip(np.asarray(x, float), 0.0, LX) / DX0
    y = np.clip(np.asarray(y, float), 0.0, LY) / DX0
    i = np.clip(np.floor(x).astype(int), 0, Z.shape[1] - 2)
    j = np.clip(np.floor(y).astype(int), 0, Z.shape[0] - 2)
    tx, ty = x - i, y - j
    z = ((1 - tx) * (1 - ty) * Z[j, i] + tx * (1 - ty) * Z[j, i + 1]
         + (1 - tx) * ty * Z[j + 1, i] + tx * ty * Z[j + 1, i + 1])
    return -z


def read_input():
    d = np.loadtxt(os.path.join(DATA, "Benchmark_2_input.txt"), skiprows=1)
    return d[:, 0], d[:, 1]


def incident():
    """The input record with its last 0.5 s tapered to zero (sin^2), then
    zero: after 22.5 s the boundary only absorbs."""
    t, e = read_input()
    e = e.copy()
    tap = t > T_REC - 0.5
    e[tap] *= np.sin(0.5 * np.pi * np.clip(T_REC - t[tap], 0.0, None) / 0.5) ** 2
    t = np.r_[t, T_REC + 0.05, 1000.0]
    e = np.r_[e, 0.0, 0.0]
    return t, e


def read_lab():
    """Gauges 5, 7, 9 [m] as given (they read 2-5 mm before the wave);
    pandas needs xlrd or python-calamine for the .xls file."""
    import pandas as pd
    a = pd.read_excel(os.path.join(DATA, "output_ch5-7-9.xls"), header=0).to_numpy(float)
    return a[:, 0], {"ch5": a[:, 1] / 100.0, "ch7": a[:, 2] / 100.0,
                     "ch9": a[:, 3] / 100.0}


def read_runup():
    """[(x, y, run-up [m] of the six repeats)] from OBS_RUNUP.txt."""
    out = []
    for l in open(os.path.join(DATA, "OBS_RUNUP.txt")):
        v = l.split()
        if len(v) == 8:
            try:
                a = [float(s) for s in v]
            except ValueError:
                continue
            out.append((a[0], a[1], np.array(a[2:])))
    return out


# --------------------------------------------------------------------------
# Grids

def grid(dx):
    """Cells in x and y, and dy: the benchmark grid (dx = 1.4 cm) or a
    coarser one with dx a multiple of 1.4 cm and ny = round(3.402/dx)."""
    nx = int(round(LX / dx))
    ny = int(round(LY / dx))
    return nx, ny, LY / ny


def nodes(dx):
    nx, ny, dy = grid(dx)
    xn = np.arange(nx + 1) * (LX / nx)
    yn = np.arange(ny + 1) * dy
    X, Y = np.meshgrid(xn, yn)
    return xn, yn, zb(X, Y)


def cell_bed(dx):
    """Bed elevation of the cells: mean of the four corner nodes."""
    _, _, Zn = nodes(dx)
    return 0.25 * (Zn[:-1, :-1] + Zn[:-1, 1:] + Zn[1:, :-1] + Zn[1:, 1:])


def gauge_cells(dx):
    nx, ny, dy = grid(dx)
    ddx = LX / nx
    out = {}
    for g, (x, y) in GAUGES.items():
        if x is None:
            x = 0.5 * ddx
        i = min(int(math.floor(x / ddx)), nx - 1)
        j = min(int(math.floor(y / dy)), ny - 1)
        out[g] = (i, j, (i + 0.5) * ddx, (j + 0.5) * dy)
    return out


# --------------------------------------------------------------------------
# NHFOAM

def nhfoam(case, dx=DX0, nz=NZ, dt=0.004, endTime=T_END, film=FILM, ks=KS,
           breaking=True, nproc=2, conv="Gauss limitedLinearV 1",
           frames=(15.0, 18.0, 0.1), writeEvery=5.0):
    """The benchmark domain with the incident wave at x = 0 (nhWaveVelocity,
    timeSeries, active absorption); free-surface snapshots every frames[2] s
    from frames[0] to frames[1] (around the times of the video frames)."""
    nx, ny, dy = grid(dx)
    mc.write(case, "system/blockMeshDict", "dictionary", f"""
convertToMeters 1;
vertices
(
    (0 0 {-D}) ({LX} 0 {-D}) ({LX} {LY} {-D}) (0 {LY} {-D})
    (0 0 0) ({LX} 0 0) ({LX} {LY} 0) (0 {LY} 0)
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
                     for g, (i, j, gx, gy) in gauge_cells(dx).items())
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
    snaps = " ".join(f"{t:g}" for t in
                     np.round(np.arange(frames[0], frames[1] + 1e-9, frames[2]), 6)) if frames else ""
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

surfaceSnapshots ({snaps});
""")
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(writeEvery / dt)), nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    mc.common_fields(case, walls=("left", "right", "south", "north"), empty=False)
    t, e = incident()
    with open(os.path.join(case, "constant", "incident.dat"), "w") as f:
        f.write(f"{len(t)}\n(\n")
        f.writelines(f"({a:.6f} {b:.8e})\n" for a, b in zip(t, e))
        f.write(")\n")
    text = f"""theory      timeSeries;
        depth       {D:g};
        direction   (1 0 0);
        origin      (0 0 0);
        eta
        {{
            type    tableFile;
            file    "$FOAM_CASE/constant/incident.dat";
        }}"""
    mc.set_patch(case, "U", "left", mc.wave_bc("left", text, absorb=True))
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

    # still water; dry columns below their bed, so that NHFOAM puts the film
    # on them (cells: x fastest, then y, then the layers)
    import ma2012
    zbc = cell_bed(dx)
    wet = -zbc > film
    eta = np.where(wet, 0.0, -10.0)
    p = np.where(wet, 0.0, G * (zbc + film))
    ma2012._set_internal(os.path.join(case, "0", "eta"),
                         ma2012._foam_list(np.tile(eta.ravel(), nz)))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"),
                         ma2012._foam_list(np.tile(p.ravel(), nz)))
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\ndx {dx:g}\ndy {dy:.10g}\nnz {nz}\ndt {dt:g}\n"
                f"film {film:g}\nks {ks:g}\nbreaking {breaking}\n")
    return dict(nx=nx, ny=ny, dy=dy, nsteps=int(round(endTime / dt)))


# --------------------------------------------------------------------------
# NHWAVE

def nhwave(case, dx=0.028, nz=NZ, endTime=T_END + 0.5, example=NHWAVE_EXAMPLE,
           px=2, cfl=0.5, kepsilon=True, plot=(15.0, 0.1)):
    """NHWAVE's example on the grid of dx (the lab part as NHFOAM's, plus
    the 1.4 m extension offshore), its wavemaker, sponge and settings;
    eta and hmax every plot[1] s from plot[0] s; probes at the gauges and
    at x = 0 (lab coordinates; NHWAVE's x = lab x + 1.4 m)."""
    import other_models as om
    nx, ny, dy = grid(dx)
    ddx = LX / nx
    next_ = int(round(X_EXT / ddx))
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    dep = np.hstack([np.full((ny, next_), D), -cell_bed(dx)])
    np.savetxt(os.path.join(case, "depth.txt"), dep, fmt="%.6f")
    for f in ("irr2d.txt",):
        with open(os.path.join(example, f)) as a, open(os.path.join(case, f), "w") as b:
            b.write(a.read())
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for g, (i, j, gx, gy) in gauge_cells(dx).items():
            f.write(f"{gx + next_ * ddx:.4f} {gy:.4f} 0.0\n")   # x y z
    txt = open(os.path.join(example, "input.txt"), encoding="latin-1").read()
    txt = om.set_keys(txt, {
        "Mglob": nx + next_, "Nglob": ny, "Kglob": nz, "PX": px, "PY": 1,
        "DX": ddx, "DY": f"{dy:.10g}", "TOTAL_TIME": endTime,
        "PLOT_START": plot[0], "PLOT_INTV": plot[1], "SCREEN_INTV": 1.0,
        "PLOT_INTV_STAT": 0.01, "NSTAT": len(GAUGES), "CFL": cfl,
        "VISCOUS_FLOW": "T" if kepsilon else "F",
        "OUT_H": "F", "OUT_E": "T", "OUT_M": "T"})
    open(os.path.join(case, "input.txt"), "w", encoding="latin-1").write(txt)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhwave\ndx {ddx:g}\ndy {dy:.10g}\nnz {nz}\n"
                f"next {next_}\nfilm {FILM:g}\n"
                f"closure {'kepsilon' if kepsilon else 'none'}\n")
    return dict(nx=nx + next_, ny=ny, dy=dy)


# --------------------------------------------------------------------------
# Results

def params(case):
    out = {}
    for l in open(os.path.join(case, "params.txt")):
        k, v = l.split(None, 1)
        v = v.strip()
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def read_nhfoam(case):
    import conical_island as ci
    return ci.read_nhfoam_gauges(case)


def read_nhwave(case):
    """Probes (NHWAVE's times); see conical_island.read_nhwave_gauges for
    the stray bytes in the file names.  An empty probe is skipped."""
    res = os.path.join(case, "results")
    files = os.listdir(res)
    out, t = {}, {}
    for k, g in enumerate(GAUGES):
        stem = f"probe_{k + 1:04d}"
        rows = []
        for f in files:
            fn = os.path.join(res, f)
            if f.startswith(stem) and os.path.getsize(fn) > 0:
                a = np.loadtxt(fn, ndmin=2)
                if a.size:
                    rows.append(a[:, :2])
        if not rows:
            continue
        d = np.vstack(rows)
        d = d[np.argsort(d[:, 0], kind="stable")]
        keep = np.r_[True, np.diff(d[:, 0]) > 0]
        t[g] = d[keep, 0]
        out[g] = d[keep, 1]
    return t, out


def lag(t, eta, t1=12.0, max_lag=1.5):
    """Lag [s] of a model's surface at x = 0 behind the input record over
    0..t1 s (the rise of the main crest, before the waves reflected from
    the island and the coast come back; NHFOAM's boundary column follows
    the record to 12.3 s), by the maximum of the cross-correlation on a
    0.01 s grid (parabolic peak)."""
    tr, er = read_input()
    tt = np.arange(0.0, t1, 0.01)
    r = np.interp(tt, tr, er)
    m = np.interp(tt, t, eta, left=0.0)
    lags = np.arange(-int(max_lag / 0.01), int(max_lag / 0.01) + 1)
    c = np.array([np.sum(r[max(0, -L):len(tt) - max(0, L)] *
                         m[max(0, L):len(tt) - max(0, -L)]) for L in lags])
    k = int(np.argmax(c))
    if 0 < k < len(c) - 1:
        a, b, cc = c[k - 1], c[k], c[k + 1]
        dk = 0.5 * (a - cc) / (a - 2 * b + cc)
    else:
        dk = 0.0
    return (lags[k] + dk) * 0.01


def model_gauges(case):
    """(t, {gauge: eta}) of a run (NHFOAM gauges or NHWAVE probes); t is a
    dict of NHWAVE's probe times, one array per probe."""
    pr = params(case)
    if pr["model"] == "nhfoam":
        t, v = read_nhfoam(case)
        return {g: t for g in v}, v
    return read_nhwave(case)


def arrival(t, eta, level=0.005, t0=14.0):
    """First time after t0 at which eta rises through level (the bore)."""
    i = np.where((t[:-1] >= t0) & (eta[:-1] < level) & (eta[1:] >= level))[0]
    if not len(i):
        return np.nan
    i = i[0]
    return t[i] + (level - eta[i]) * (t[i + 1] - t[i]) / (eta[i + 1] - eta[i])


def _times(case):
    return sorted((d for d in os.listdir(case)
                   if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"), key=float)


def nhwave_outputs(case):
    """{file number: NHWAVE time} of the field outputs (log.txt)."""
    out = {}
    for l in open(os.path.join(case, "log.txt"), errors="replace"):
        m = re.search(r"Printing file No\.\s*(\d+)\s+TIME/TOTAL:\s*([-0-9.Ee+]+)", l)
        if m:
            out[int(m.group(1))] = float(m.group(2))
    return out


def _snapshot(case, name=None):
    """Column bottoms and free surface (ny, nx) of an NHFOAM snapshot
    (postProcessing/surface/<name>; default: the first one)."""
    pr = params(case)
    nx, ny, dy = grid(pr["dx"])
    root = os.path.join(case, "postProcessing", "surface")
    if name is None:
        name = min(os.listdir(root), key=lambda s: float(s[:-4]))
    a = np.loadtxt(os.path.join(root, name))
    i = np.floor(a[:, 0] / (LX / nx)).astype(int)
    j = np.floor(a[:, 1] / dy).astype(int)
    zc = np.full((ny, nx), np.nan)
    e = np.full((ny, nx), np.nan)
    zc[j, i] = a[:, 2]
    e[j, i] = a[:, 3]
    return zc, e


def surface(case, time):
    """Free surface of the columns (ny, nx) on the lab part of the grid at
    a time (NHFOAM: the nearest snapshot; NHWAVE: the nearest output), the
    bed, the wet mask and the time of the snapshot."""
    pr = params(case)
    nx, ny, dy = grid(pr["dx"])
    z = cell_bed(pr["dx"])
    if pr["model"] == "nhfoam":
        root = os.path.join(case, "postProcessing", "surface")
        tt = min(os.listdir(root), key=lambda s: abs(float(s[:-4]) - time))
        zc, e = _snapshot(case, tt)
        return e, z, e - zc > 1.5 * pr["film"], float(tt[:-4])
    outs = nhwave_outputs(case)
    k = min(outs, key=lambda k: abs(outs[k] - time))
    e = np.loadtxt(os.path.join(case, "results", f"eta_{k:05d}"))[:, int(pr["next"]):]
    return e, z, e - z > 1.5 * pr["film"], outs[k]


def emax(case):
    """Highest surface of each column while wet, the flooded mask and the
    bed (lab part of the grid): NHFOAM etaMax of the last time (never-wet
    columns keep their bottom), NHWAVE hmax of the last output (zero where
    never wet)."""
    import conical_island as ci
    pr = params(case)
    nx, ny, dy = grid(pr["dx"])
    z = cell_bed(pr["dx"])
    if pr["model"] == "nhfoam":
        tt = _times(case)[-1]
        e = ci._field_values(os.path.join(case, tt, "etaMax"))[: nx * ny].reshape(ny, nx)
        zc = _snapshot(case)[0]          # NHFOAM's own column bottoms
        return e, (z < 0.0) | (e > zc + 0.5 * pr["film"]), z
    files = sorted(f for f in os.listdir(os.path.join(case, "results"))
                   if re.fullmatch(r"hmax_[0-9]{5}", f))
    h = np.loadtxt(os.path.join(case, "results", files[-1]))[:, int(pr["next"]):]
    return h, (z < 0.0) | (h > 0.0), z


def runup(case, half_width=None):
    """Run-up at the lab points: the highest bed among the flooded columns
    in the strip |y - y_lab| <= half_width (default: half a cell) on land
    within 0.3 m of the lab x; also the most inland x reached there and the
    highest surface of the edge columns (as for the conical island)."""
    pr = params(case)
    e, fl, z = emax(case)
    nx, ny, dy = grid(pr["dx"])
    ddx = LX / nx
    xc = (np.arange(nx) + 0.5) * ddx
    yc = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(xc, yc)
    hw = half_width if half_width is not None else 0.5 * dy
    dry_nb = np.zeros_like(fl)
    dry_nb[:, 1:] |= ~fl[:, :-1]
    dry_nb[:, :-1] |= ~fl[:, 1:]
    dry_nb[1:, :] |= ~fl[:-1, :]
    dry_nb[:-1, :] |= ~fl[1:, :]
    edge = fl & dry_nb & (z > 0.0)
    out = []
    for x, y, obs in read_runup():
        strip = (np.abs(Y - y) <= hw + 1e-9) & (np.abs(X - x) <= 0.3)
        land = strip & fl & (z > 0.0)
        R = z[land].max() if land.any() else np.nan
        xr = X[land].max() if land.any() else np.nan
        sel = strip & edge
        Re = e[sel].max() if sel.any() else np.nan
        out.append(dict(x=x, y=y, obs=obs, R=R, x_inland=xr, R_edge=Re))
    return out


# --------------------------------------------------------------------------
# All runs of the README/report (section 6.11):
#     python3 monai_valley.py <dir>     (then <dir>/runall.sh, about 11 h on
#                                        2 cores: NHFOAM 1.4 cm 6.2 h,
#                                        NHWAVE 2.8 cm 2 h each, NHFOAM
#                                        2.8 cm 35 min each)

if __name__ == "__main__":
    import sys
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    os.chdir(out)
    nhfoam("nf_14", dx=0.014, dt=0.004)
    nhfoam("nf_28", dx=0.028, dt=0.008)
    nhfoam("nf_28_nobrk", dx=0.028, dt=0.008, breaking=False)
    nhwave("nw_28", dx=0.028)
    nhwave("nw_28_inv", dx=0.028, kepsilon=False)
    with open("runall.sh", "w") as f:
        f.write("#!/bin/sh\n# NHFOAM (Allrun: 2 processors) and NHWAVE (2 processors)\n"
                "NHWAVE=${NHWAVE:-nhwave}\n")
        for c in ("nf_14", "nf_28", "nf_28_nobrk"):
            f.write(f"(cd {c} && ./Allrun)\n")
        for c in ("nw_28", "nw_28_inv"):
            f.write(f"(cd {c} && mpirun -np 2 $NHWAVE > log.nhwave 2>&1)\n")
    os.chmod("runall.sh", 0o755)
