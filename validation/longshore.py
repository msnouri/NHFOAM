"""Wave-induced longshore current on a plane beach: Visser (1991) case 4, the
test of Ma, Shi & Kirby (2012), Section 4.7.

Regular waves, H = 7.8 cm (amplitude 3.9 cm), T = 1.02 s, 15.4 degrees from
the shore normal at h = 0.35 m, break on a smooth 1:20 concrete slope and
drive a longshore current.  The basin is periodic alongshore over one
alongshore wavelength, Ly = 2 pi/(k sin(theta)).

x is measured from the toe of the slope; the still-water shoreline is at
x = 7 m.  NHFOAM makes stream-function waves in a relaxation zone on a flat
bottom ahead of the toe (Ma et al. had 0.6 m of flat bottom and the waves
from the boundary).

Measurements (Visser 1991, as tabulated by Kraus & Larson 1991, NMLONG,
Table B3): distance offshore of the still-water shoreline [m], wave height
[m], longshore current [m/s] (depth averaged, (Vs + 2 Vm + Vb)/4 of the dye
velocities at the surface, mid-depth and near the bottom), mean water level
[m].  NaN: not measured.
"""
import os
import re
import numpy as np
import make_case as mc
from stream_function import fenton, foam_entries

V_H, V_T, V_D, V_SLOPE, V_ANGLE = 0.078, 1.02, 0.35, 1.0 / 20.0, 15.4
V_XSHORE = V_D / V_SLOPE          # still-water shoreline, from the toe

nan = np.nan
VISSER4 = np.array([
    # x_off   H      V      mwl
    [0.01,   nan,  0.000,  0.016],
    [0.21,   nan,  0.114,  0.015],
    [0.41,   nan,  0.275,  0.014],
    [0.81,  0.016, 0.332,  0.011],
    [1.21,  0.026, 0.389,  0.008],
    [1.61,  0.034, 0.375,  0.003],
    [2.01,  0.056, 0.301, -0.001],
    [2.41,  0.087, 0.218, -0.001],
    [2.81,  0.084, 0.139, -0.001],
    [3.21,  0.078, 0.099, -0.001],
    [3.61,  0.080, 0.063, -0.001],
    [4.01,  0.079, 0.034, -0.001],
    [4.41,  0.076, 0.016,  nan],
    [4.81,  0.076, 0.007,  nan],
    [5.21,  0.075, 0.000,  nan],
    [5.61,  0.072, -0.006, nan],
])
# Visser (1991), experiment 4: breaker line 2.51 m and plunge line 1.92 m
# from the shoreline, largest longshore current 0.404 m/s
V_XBR, V_XPL, V_VMAX = 2.51, 1.92, 0.404

AGE = "\n    betaFinal   0.15;\n    Tstar       5;"


def wave(N=24):
    sol = fenton(V_H, V_T, V_D, N=N)
    th = np.radians(V_ANGLE)
    Ly = 2 * np.pi / (sol["k"] * np.sin(th))
    return sol, th, Ly


def visser(case, dx=0.04, ny=40, nz=5, dt=0.0102, tEnd=100.0, avgPeriods=5,
           zone=1.5, flat=0.5, xLand=0.6, film=0.003, ks=0.0008,
           brk_model="hybrid", brk_extra=AGE + "\n    delta       1.2;",
           alpha=0.6, beta=0.3, les_Ck=0.094, N=24, gauge_dx=0.05,
           writeEvery=10.0, conv="Gauss limitedLinearV 1"):
    """NHFOAM case of Visser's case 4.  les_Ck: Smagorinsky subgrid model of
    the OpenFOAM library (None: laminar).  brk_model None: no breaking
    model.  The mean velocity UMean is averaged over the last avgPeriods
    periods; eta is recorded along y = Ly/2 and y = Ly/4."""
    sol, th, Ly = wave(N)
    dy = Ly / ny
    x0 = -(zone + flat)
    xEnd = V_XSHORE + xLand
    xs, hs, zt = mc.beach_profile([x0, 0.0, xEnd],
                                  [V_D, V_D, V_D - V_SLOPE * xEnd], film, dx=dx)
    gx = np.arange(-1.0, V_XSHORE + 0.3 + 1e-9, gauge_dx)
    gl = []
    for j, y in enumerate((0.5 * Ly, 0.25 * Ly)):
        gl += [f"    G{j}_{i:03d}    ({x:.4f} {y:.6g} 0);" for i, x in enumerate(gx)]
    brk = "" if brk_model is None else f"""
breaking
{{
    model       {brk_model};
    alpha       {alpha:g};
    beta        {beta:g};{brk_extra}
}}
"""
    d = np.cos(th), np.sin(th)
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{chr(10).join(gl)}
}}

forcingZones
{{
    inlet
    {{
        type        wave;
        start       ({x0 + zone:g} 0 0);
        end         ({x0:g} 0 0);
        sigmaMax    20;
        exponent    2;
{foam_entries(sol, indent="        ")}        direction   ({d[0]:.10g} {d[1]:.10g} 0);
        origin      ({x0:g} 0 0);
        rampTime    {3 * V_T:g};
    }}
}}

wettingDrying
{{
    minDepth    {film};
}}

bottomFriction
{{
    ks          {ks};
}}
{brk}""")
    nsteps = int(round(tEnd / dt))
    mc.system_files(case, endTime=nsteps * dt, deltaT=dt,
                    writeInterval=int(round(writeEvery / dt)), nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    tAvg = nsteps * dt - avgPeriods * V_T
    nh = open(os.path.join(case, "constant/NHFOAMProperties")).read()
    nh = nh.replace("\ngauges\n", f"""
averaging
{{
    timeStart   {tAvg:.10g};
}}

gauges
""", 1)
    open(os.path.join(case, "constant/NHFOAMProperties"), "w").write(nh)
    mc.blockmesh_profile(case, xs, hs, Ly, dx, nz, Ny=ny, ztop=zt, cyclic=True)
    mc.common_fields(case, walls=("left", "right"), empty=False,
                     constraint_types=True)
    if les_Ck:
        mc.les_files(case, Ck=les_Ck, patches={
            "surface": "zeroGradient", "bottom": "zeroGradient",
            "left": "zeroGradient", "right": "zeroGradient",
            "front": "cyclic", "back": "cyclic"})
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"T {V_T}\nH {V_H}\nLy {Ly:.10g}\nny {ny}\nnz {nz}\ndx {dx}\n"
                f"k {sol['k']:.10g}\ntAvg {tAvg:.10g}\n")
    return Ly


# --------------------------------------------------------------------------
# Post-processing

def _read_field(path):
    """Internal field of an ASCII OpenFOAM vol field (scalar or vector)."""
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*\n?(\d+)\s*\(", txt)
    if not m:
        m = re.search(r"internalField\s+uniform\s+(\S+);", txt)
        raise ValueError("uniform field " + path)
    n = int(m.group(2))
    body = txt[m.end():]
    if m.group(1) == "vector":
        vals = re.findall(r"\(([^()]+)\)", body[: body.find("\n)")])
        return np.array([[float(v) for v in s.split()] for s in vals[:n]])
    vals = body[: body.find("\n)")].split()
    return np.array([float(v) for v in vals[:n]])


def _read_list(path):
    """Body of an ASCII OpenFOAM list file (points, faces, owner, ...)."""
    txt = open(path).read()
    txt = txt[txt.index("FoamFile"):]
    txt = txt[txt.index("}") + 1:]
    m = re.search(r"\n(\d+)\s*\n\(", txt)
    return int(m.group(1)), txt[m.end():txt.rfind(")")]


def cell_centres(case, time):
    """Cell centres as the mean of the cell's points (exact in x and y for
    the extruded sigma columns), from the mesh of that time (or the initial
    mesh).  The OpenFOAM function objects (writeCellCentres) do not run with
    every OpenFOAM build, hence this."""
    pdir = os.path.join(case, time, "polyMesh")
    if not os.path.exists(os.path.join(pdir, "points")):
        pdir = os.path.join(case, "constant", "polyMesh")
    cdir = os.path.join(case, "constant", "polyMesh")
    n, body = _read_list(os.path.join(pdir, "points"))
    pts = np.array(body.replace("(", " ").replace(")", " ").split(),
                   dtype=float).reshape(n, 3)
    n, body = _read_list(os.path.join(cdir, "faces"))
    toks = re.findall(r"(\d+)\(([^)]*)\)", body)
    faces = [np.array(t[1].split(), dtype=int) for t in toks]
    n, body = _read_list(os.path.join(cdir, "owner"))
    own = np.array(body.split(), dtype=int)
    n, body = _read_list(os.path.join(cdir, "neighbour"))
    nei = np.array(body.split(), dtype=int)
    nc = max(own.max(), nei.max()) + 1
    cellpts = [set() for _ in range(nc)]
    for f, c in enumerate(own):
        cellpts[c].update(faces[f])
    for f, c in enumerate(nei):
        cellpts[c].update(faces[f])
    return np.array([pts[list(cp)].mean(axis=0) for cp in cellpts])


def latest_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            ts.append((float(d), d))
        except ValueError:
            pass
    return max(ts)[1]


def mean_current(case, time=None, field="UMean"):
    """Cross-shore profile of the time-mean depth-averaged velocity: the
    layer average of UMean in each column (uniform sigma layers), averaged
    alongshore.  Returns x, U_x, U_y and the alongshore spread (std) of U_y."""
    time = time or latest_time(case)
    C = cell_centres(case, time)
    U = _read_field(os.path.join(case, time, field))
    xr = np.round(C[:, 0], 6)
    xs = np.unique(xr)
    ux, uy, sy = [], [], []
    yr = np.round(C[:, 1], 6)
    for x in xs:
        m = xr == x
        ux.append(U[m, 0].mean())
        uy.append(U[m, 1].mean())
        # depth averages of the columns at this x
        cols = [U[m & (yr == y), 1].mean() for y in np.unique(yr[m])]
        sy.append(np.std(cols))
    return xs, np.array(ux), np.array(uy), np.array(sy)


def gauges(case):
    f = os.path.join(case, "postProcessing", "gauges", "0", "eta.dat")
    names, pts = [], []
    for line in open(f):
        if line.startswith("# G"):
            names.append(line.split()[1])
            pts.append([float(v) for v in line.split("(")[1].split(")")[0].split()])
    d = np.loadtxt(f, comments="#")
    return names, np.array(pts), d[:, 0], d[:, 1:]


def wave_stats(case, nper=5, line=0):
    """Wave height (mean of the per-period crest-to-trough heights) and mean
    level over the last nper periods along gauge line 0 (y = Ly/2) or 1."""
    names, pts, t, e = gauges(case)
    sel = np.array([n.startswith(f"G{line}_") for n in names])
    x = pts[sel, 0]
    e = e[:, sel]
    tend = t[-1]
    Hs = []
    for p in range(nper):
        m = (t > tend - (p + 1) * V_T) & (t <= tend - p * V_T)
        Hs.append(e[m].max(axis=0) - e[m].min(axis=0))
    m = t > tend - nper * V_T
    return x, np.mean(Hs, axis=0), e[m].mean(axis=0)


# --------------------------------------------------------------------------
# NHWAVE on Ma et al.'s domain (8.4 m: 0.6 m flat, the slope, dry land to a
# bed 4 cm above the still water level), its linear left-boundary wavemaker
# (AMP is the wave height), periodic alongshore (Nglob a power of two; NHWAVE
# adjusts the angle so that the alongshore wavelength fits), 3D Smagorinsky
# (IVTURB = 20, Cvs = Cs), log-law friction (Ibot = 2, Zob = ks), and its
# wave averages (wave height, set-up, depth-averaged mean velocity) over the
# last avgPeriods periods.

def linear_k(omega, h, g=9.81):
    k = omega ** 2 / g
    for _ in range(50):
        k -= (g * k * np.tanh(k * h) - omega ** 2) / (
            g * np.tanh(k * h) + g * k * h / np.cosh(k * h) ** 2)
    return k


def nhwave_visser(case, template_dir, dx=0.04, ny=32, nz=5, tEnd=100.0,
                  avgPeriods=5, ks=0.0008, Cs=0.168, mindep=0.005, cfl=0.5):
    import other_models as om
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    flat, xLand = 0.6, 0.8
    Lx = flat + V_XSHORE + xLand
    Nx = int(round(Lx / dx))
    # alongshore wavelength of linear theory (NHWAVE's wavemaker)
    k = linear_k(2 * np.pi / V_T, V_D)
    Ly = 2 * np.pi / (k * np.sin(np.radians(V_ANGLE)))
    dy = Ly / ny
    x = (np.arange(Nx) + 0.5) * dx
    h = np.where(x < flat, V_D, V_D - (x - flat) * V_SLOPE)
    np.savetxt(os.path.join(case, "depth.txt"), np.tile(h, (ny, 1)), fmt="%.6f")
    txt = open(os.path.join(template_dir, "input.txt")).read()
    txt = om.set_keys(txt, {
        "Mglob": Nx, "Nglob": ny, "Kglob": nz, "DX": dx, "DY": dy,
        "PX": 1, "PY": 1, "TOTAL_TIME": tEnd, "PLOT_START": tEnd - 10 * V_T,
        "PLOT_INTV": V_T, "SCREEN_INTV": 1.0, "DEPTH_TYPE": "CELL_CENTER",
        "Ibot": 2, "Zob": ks, "VISCOUS_FLOW": "T", "IVTURB": 20, "IHTURB": 10,
        "Cvs": Cs, "MinDep": mindep, "CFL": cfl, "PERIODIC_Y": "T",
        "WAVEMAKER": "LEF_LIN", "AMP": V_H, "PER": V_T, "DEP": V_D,
        "THETA": V_ANGLE, "SPONGE_ON": "F", "WAVE_AVERAGE_ON": "T",
        "WAVE_AVERAGE_START": tEnd - avgPeriods * V_T,
        "WAVE_AVERAGE_END": tEnd, "WaveheightID": 1, "NSTAT": 0,
        "OUT_E": "T", "OUT_U": "F", "OUT_V": "F", "OUT_H": "T",
        "INITIAL_EUVW": "F"})
    open(os.path.join(case, "input.txt"), "w").write(txt)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"Nx {Nx}\nNy {ny}\ndx {dx}\ndy {dy:.10g}\nLy {Ly:.10g}\n"
                f"flat {flat}\nxshore {flat + V_XSHORE}\n")
    return Ly


def nhwave_profiles(case):
    """Alongshore averages of NHWAVE's wave averages: x from the still-water
    shoreline (offshore positive), H, set-up, U, V (depth averaged)."""
    p = dict(l.split() for l in open(os.path.join(case, "params.txt")))
    Nx, Ny, dx = int(p["Nx"]), int(p["Ny"]), float(p["dx"])
    xoff = float(p["xshore"]) - (np.arange(Nx) + 0.5) * dx
    out = [xoff]
    for f in ("waveheight", "setup", "umean", "vmean"):
        a = np.loadtxt(os.path.join(case, "results", f)).reshape(Ny, Nx)
        out.append(a.mean(axis=0))
    return out
