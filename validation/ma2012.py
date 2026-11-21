#!/usr/bin/env python3
"""Test cases of Ma, Shi & Kirby (2012), "Shock-capturing non-hydrostatic
model for fully dispersive surface wave processes", Ocean Modelling 43-44,
22-35 (the NHWAVE paper), set up identically for NHFOAM and NHWAVE.

    4.1  standing wave in a closed basin      basin_nhfoam / basin_nhwave
    4.2  solitary wave in constant depth      solitary_nhfoam / solitary_nhwave
    4.3  periodic waves over a submerged bar  bar_nhfoam / bar_nhwave
    4.4  elliptic shoal (Berkhoff et al.)     berkhoff_nhfoam / berkhoff_nhwave
    (4.6, the landslide of Enet & Grilli, uses the cases of section 5.5 of the
    README; 4.5 and 4.7 need wetting and drying)

usage (python):
    import ma2012
    ma2012.basin_nhfoam("nf_basin_kD1_L3", D=10.0, nLayers=3)
    ma2012.basin_nhwave("nw_basin_kD1_L3", template, D=10.0, nLayers=3)
"""
import os, math, re
import numpy as np
import make_case as mc
import other_models as om

g = 9.81

# --------------------------------------------------------------------------
# 4.1 Standing wave in a closed basin: L = 20 m, eta = a cos(pi x/10),
# a = 0.1 m, depth D = 10, 20, 30 m (kD = pi, 2pi, 3pi), dx = 0.2 m,
# 30 s, gauges at x = 0.1 m and 17.5 m.

BASIN_L = 20.0
BASIN_A = 0.1
BASIN_DX = 0.2
BASIN_TEND = 30.0
BASIN_GAUGES = [0.1, 17.5]

# Fig. 3 of the paper (read off the figure): normalised rms error at
# x = 17.5 m against the number of layers
BASIN_PAPER = {
    1: {3: 0.0163, 5: 0.0090, 8: 0.0061, 10: 0.0053},
    2: {3: 0.0307, 5: 0.0150, 8: 0.0081, 10: 0.0064},
    3: {3: 0.0490, 5: 0.0250, 8: 0.0135, 10: 0.0100},
}


def basin_period(D):
    k = 2 * math.pi / BASIN_L
    return 2 * math.pi / math.sqrt(g * k * math.tanh(k * D))


def basin_nhfoam(case, D, nLayers, dt=0.005):
    Nx = int(round(BASIN_L / BASIN_DX))
    dy = BASIN_DX
    gauges = "".join(f"    G{i}      ({x:.6g} {dy/2:.6g} 0);\n"
                     for i, x in enumerate(BASIN_GAUGES))
    mc.constant_files(case, nh_extra=f"\ngauges\n{{\n{gauges}}}\n")
    nSteps = int(round(BASIN_TEND / dt))
    mc.system_files(case, endTime=BASIN_TEND, deltaT=dt, writeInterval=nSteps)
    mc.blockmesh_flat(case, BASIN_L, D, dy, Nx, nLayers)
    mc.common_fields(case)
    mc.write(case, "system/nhSetWavesDict", "dictionary", f"""
type            standing;
amplitude       {BASIN_A:.10g};
wavelength      {BASIN_L:.10g};
origin          (0 0 0);
direction       (1 0 0);
""")
    mc.allrun(case, """
    runApplication blockMesh
    runApplication nhSetWaves
    runApplication NHFOAM
    """)
    return dict(T=basin_period(D))


def basin_nhwave(case, template, D, nLayers, cfl=0.5):
    Nx = int(round(BASIN_L / BASIN_DX))
    dx = BASIN_DX
    x = (np.arange(Nx) + 0.5) * dx
    T = basin_period(D)
    keys = om.nhwave_common(Nx, 1, nLayers, dx, dx, BASIN_TEND, px=1)
    keys.update({"WAVEMAKER": "LEF_LIN", "AMP": 0.0, "PER": T, "DEP": D,
                 "BC_X0": 1, "BC_Xn": 1, "SPONGE_ON": "F",
                 "Sponge_West_Width": 0.0, "Sponge_East_Width": 0.0,
                 "VISCOUS_FLOW": "F", "INITIAL_EUVW": "T", "CFL": cfl,
                 "DT_INI": 0.001, "DT_MAX": 1.0, "PLOT_INTV_STAT": 0.005})
    om.nhwave_write(case, template, keys, np.full((1, Nx), D),
                    [(xg, 0.5 * dx) for xg in BASIN_GAUGES])
    np.savetxt(os.path.join(case, "eta0.txt"),
               (BASIN_A * np.cos(2 * math.pi * x / BASIN_L))[None, :], fmt="%.12e")
    np.savetxt(os.path.join(case, "uvw0.txt"), np.zeros((3 * nLayers, Nx)), fmt="%.1f")
    return dict(T=T)


def basin_error(t, eta, D, x=17.5):
    """Normalised rms error of Ma et al. (2012): sqrt(mean((eta_a - eta)^2))/H
    with the linear standing wave eta_a = a cos(kx) cos(sigma t), H = 2a."""
    k = 2 * math.pi / BASIN_L
    sig = math.sqrt(g * k * math.tanh(k * D))
    ea = BASIN_A * math.cos(k * x) * np.cos(sig * t)
    return math.sqrt(np.mean((ea - eta) ** 2)) / (2 * BASIN_A)


# --------------------------------------------------------------------------
# 4.2 Solitary wave in constant depth: h = 1 m, H/h = 0.5, 200 m long,
# dx = 0.1 m, 3 layers, crest at x = 8 m at t = 0, snapshots at 10, 20, 40 s.
# Initial condition: the exact (fully nonlinear) solitary wave.

SOL_H = 0.5
SOL_h = 1.0
SOL_LEN = 200.0
SOL_DX = 0.1
SOL_X0 = 8.0
SOL_TEND = 40.0
SOL_TIMES = [10.0, 20.0, 40.0]


def _foam_list(values, vector=False):
    n = len(values)
    if vector:
        body = "\n".join(f"({u:.10g} 0 {w:.10g})" for u, w in values)
        return f"nonuniform List<vector>\n{n}\n(\n{body}\n)"
    body = "\n".join(f"{v:.10g}" for v in values)
    return f"nonuniform List<scalar>\n{n}\n(\n{body}\n)"


def _set_internal(path, text):
    s = open(path).read()
    i = s.index("internalField")
    j = s.index(";", i)
    open(path, "w").write(s[:i] + "internalField   " + text + s[j:])


def solitary_fields(sw, x, z, h):
    """eta per column, and u, w, p_rgh at points (x, z) (fixed frame), for
    a wave with its crest at x = 0."""
    u, w = sw.velocity(x, z)
    c = sw.c
    # Bernoulli in the frame of the wave: p_rgh = p/rho + g z = (c^2 - |v|^2)/2
    p_rgh = 0.5 * (c * c - (u - c) ** 2 - w ** 2)
    return u, w, p_rgh


def solitary_nhfoam(case, nLayers=3, dt=0.005, H=SOL_H, conv=None):
    """Writes the case and its initial fields (exact solitary wave of
    Babenko's equation; the cell centres of the flat single-block mesh are
    computed here, so no OpenFOAM utility is run)."""
    from solitary_exact import SolitaryWave
    Nx = int(round(SOL_LEN / SOL_DX))
    mc.constant_files(case)
    mc.system_files(case, endTime=SOL_TEND, deltaT=dt,
                    writeInterval=int(round(10.0 / dt)), conv=conv)
    mc.blockmesh_flat(case, SOL_LEN, SOL_h, SOL_DX, Nx, nLayers)
    mc.common_fields(case)
    mc.allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)
    # cell centres of the single-block flat mesh: x fastest, then z
    i = np.arange(Nx)
    k = np.arange(nLayers)
    xc = np.tile((i + 0.5) * SOL_DX, nLayers)
    sigma = np.repeat((k + 0.5) / nLayers, Nx)
    sw = SolitaryWave(H, d=SOL_h)
    eta = sw.eta(xc - SOL_X0)
    zc = -SOL_h + sigma * (SOL_h + eta)
    u, w, prgh = solitary_fields(sw, xc - SOL_X0, zc, SOL_h)
    _set_internal(os.path.join(case, "0", "eta"), _foam_list(eta))
    _set_internal(os.path.join(case, "0", "U"), _foam_list(list(zip(u, w)), vector=True))
    _set_internal(os.path.join(case, "0", "p_rgh"), _foam_list(prgh))
    return dict(c=sw.c, F=sw.F)


def _read_vectors(path):
    import re
    s = open(path).read().split("internalField")[1].split(";")[0]
    return np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", s), float)


def solitary_nhwave(case, template, nLayers=3, cfl=0.5, H=SOL_H):
    from solitary_exact import SolitaryWave
    Nx = int(round(SOL_LEN / SOL_DX))
    dx = SOL_DX
    x = (np.arange(Nx) + 0.5) * dx
    keys = om.nhwave_common(Nx, 1, nLayers, dx, dx, SOL_TEND, px=1)
    keys.update({"WAVEMAKER": "LEF_LIN", "AMP": 0.0, "PER": 1.0, "DEP": SOL_h,
                 "BC_X0": 1, "BC_Xn": 1, "SPONGE_ON": "F",
                 "Sponge_West_Width": 0.0, "Sponge_East_Width": 0.0,
                 "VISCOUS_FLOW": "F", "INITIAL_EUVW": "T", "CFL": cfl,
                 "DT_INI": 0.001, "DT_MAX": 1.0, "PLOT_INTV": 10.0,
                 "PLOT_INTV_STAT": 0.01, "OUT_E": "T", "OUT_U": "T", "OUT_W": "T"})
    om.nhwave_write(case, template, keys, np.full((1, Nx), SOL_h),
                    [(SOL_X0, 0.5 * dx)])
    sw = SolitaryWave(H, d=SOL_h)
    eta = sw.eta(x - SOL_X0)
    sig = (np.arange(nLayers) + 0.5) / nLayers
    U = np.zeros((nLayers, Nx)); W = np.zeros((nLayers, Nx))
    for k, s in enumerate(sig):
        zc = -SOL_h + s * (SOL_h + eta)
        U[k], W[k] = sw.velocity(x - SOL_X0, zc)
    np.savetxt(os.path.join(case, "eta0.txt"), eta[None, :], fmt="%.12e")
    np.savetxt(os.path.join(case, "uvw0.txt"),
               np.vstack([U, np.zeros_like(U), W]), fmt="%.12e")
    return dict(c=sw.c, F=sw.F)


# --------------------------------------------------------------------------
# 4.3 Periodic waves over a submerged bar (Beji & Battjes 1993), NHWAVE's
# examples/submerged_bar: T = 2.02 s, H = 2 cm, 3 layers, dx = 2 cm, 40 s.
# NHWAVE generates the waves at its left boundary (x = 0 in Beji & Battjes'
# coordinates); NHFOAM has its wave-making zone just upstream of that point,
# so the waves cross the same distance to the bar.  Both are written in the
# coordinates of the FUNWAVE benchmark data (Beji & Battjes x + 20 m).

BAR_XOFF = 20.0
BAR_GAUGES_BB = [2.0, 4.0, 10.5, 12.5, 13.5, 14.5, 15.7, 17.3, 19.0, 21.0]


def bar_nhfoam(case, nLayers=3, dx=0.02, dt=0.005, endTime=40.0, H=0.02, T=2.02,
               conv=None):
    xs = [12, 26, 32, 34, 37, 55]
    hs = [0.4, 0.4, 0.1, 0.1, 0.4, 0.4]
    dy = 0.02
    gtxt = "\n".join(f"    G{g + BAR_XOFF:g}    ({g + BAR_XOFF:g} {dy/2:g} 0);"
                     for g in BAR_GAUGES_BB)
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

forcingZones
{{
    inlet
    {{
        type        wave;
        start       (20 0 0);
        end         (12 0 0);
        sigmaMax    20;
        exponent    2;
        waveHeight  {H:g};
        wavePeriod  {T:g};
        depth       0.4;
        direction   (1 0 0);
        rampTime    {2*T:g};
    }}
    outlet
    {{
        type        sponge;
        start       (45 0 0);
        end         (55 0 0);
        sigmaMax    20;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps,
                    nCorr=3, nNonOrth=1, pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dy, dx, nLayers)
    mc.common_fields(case)
    mc.allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)


def bar_nhwave(case, example_dir, nLayers=3, px=2, total=40.0):
    """NHWAVE's own example, with all ten gauges and px MPI ranks."""
    import shutil
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    shutil.copy(os.path.join(example_dir, "depth.txt"), case)
    txt = open(os.path.join(example_dir, "input.txt")).read()
    txt = om.set_keys(txt, {"PX": px, "PY": 1, "NSTAT": len(BAR_GAUGES_BB),
                            "Kglob": nLayers, "TOTAL_TIME": total})
    open(os.path.join(case, "input.txt"), "w").write(txt)
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for g in BAR_GAUGES_BB:
            f.write(f"{g:g} 0.01 0.0\n")


# --------------------------------------------------------------------------
# 4.4 Elliptic shoal on a sloping bottom (Berkhoff et al. 1982), NHWAVE's
# examples/berkhoff: T = 1 s, H = 4.64 cm, walls at y = 0 and 20 m, sponge in
# the last 3 m. Coordinates of NHWAVE's domain (x along the waves, 0 - 30 m);
# the bathymetry of its depth.txt is matched to 2 mm by the formulas of Ma et
# al. (2012), eqs. (49)-(51), with the shoal centre at (12.025, 10.025).
# The measurement sections are relative to the shoal centre: 1-5 across the
# flume at x = 1, 3, 5, 7, 9 m, 6-8 along it at y = -2, 0, 2 m.

BERK_C = (12.025, 10.025)
BERK_H, BERK_T, BERK_h0 = 0.0464, 1.0, 0.45
BERK_SEC_X = [1.0, 3.0, 5.0, 7.0, 9.0]
BERK_SEC_Y = [-2.0, 0.0, 2.0]


def berk_depth(x, y):
    a = np.deg2rad(20.0)
    xr = (x - BERK_C[0]) * np.cos(a) - (y - BERK_C[1]) * np.sin(a)
    yr = (x - BERK_C[0]) * np.sin(a) + (y - BERK_C[1]) * np.cos(a)
    h = np.where(xr < -5.84, 0.45, np.maximum(0.07, 0.45 - 0.02 * (5.84 + xr)))
    inside = (xr / 3) ** 2 + (yr / 4) ** 2 < 1
    th = -0.3 + 0.5 * np.sqrt(np.clip(1 - (xr / 3.75) ** 2 - (yr / 5) ** 2, 0, None))
    return np.where(inside, h - th, h)


def berk_points():
    """[(section, xs, ys)] in the shoal frame (x along the flume), every
    0.5 m (174 points: NHWAVE takes at most 200 stations)."""
    pts = []
    ys = np.arange(-5.0, 5.0 + 1e-9, 0.5)
    for k, xs in enumerate(BERK_SEC_X):
        pts += [(k + 1, xs, y) for y in ys]
    xs = np.arange(0.0, 11.0 + 1e-9, 0.5)
    for k, ys0 in enumerate(BERK_SEC_Y):
        pts += [(k + 6, x, ys0) for x in xs]
    return pts


def berkhoff_nhfoam(case, dx=0.05, dy=0.1, nLayers=3, dt=0.01, endTime=30.0,
                    conv=None, x0=-3.0, x1=30.0, W=20.0):
    """Flat blockMesh, then the points are moved onto the bathymetry
    (z = -h + sigma h); wave zone in x0..0, sponge in x1-3..x1."""
    import subprocess
    nx, ny = int(round((x1 - x0) / dx)), int(round(W / dy))
    h0 = BERK_h0
    mc.write(case, "system/blockMeshDict", "dictionary", f"""
convertToMeters 1;
vertices
(
    ({x0} 0 {-h0}) ({x1} 0 {-h0}) ({x1} {W} {-h0}) ({x0} {W} {-h0})
    ({x0} 0 0) ({x1} 0 0) ({x1} {W} 0) ({x0} {W} 0)
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nLayers}) simpleGrading (1 1 1) );
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
    gtxt = "\n".join(f"    S{s}_{i}    ({BERK_C[0] + xs:.4f} {BERK_C[1] + ys:.4f} 0);"
                     for i, (s, xs, ys) in enumerate(berk_points()))
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

forcingZones
{{
    inlet
    {{
        type        wave;
        start       (0 0 0);
        end         ({x0} 0 0);
        sigmaMax    20;
        exponent    2;
        waveHeight  {BERK_H:g};
        wavePeriod  {BERK_T:g};
        depth       {h0:g};
        direction   (1 0 0);
        rampTime    {2*BERK_T:g};
    }}
    outlet
    {{
        type        sponge;
        start       ({x1 - 3:g} 0 0);
        end         ({x1:g} 0 0);
        sigmaMax    20;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps,
                    nCorr=3, nNonOrth=1, pTol=1e-9, conv=conv)
    mc.common_fields(case, walls=("left", "right", "south", "north"), empty=False)
    mc.allrun(case, """
    runApplication NHFOAM
    """)
    # mesh: flat box, then the points onto the bathymetry
    subprocess.run(["blockMesh", "-case", case], check=True,
                   stdout=open(os.path.join(case, "log.blockMesh"), "w"))
    pf = os.path.join(case, "constant/polyMesh/points")
    txt = open(pf).read()
    head, rest = txt.split("(", 1)
    n = int(head.strip().split()[-1])
    body, tail = rest.rsplit(")", 1)
    P = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", body), float)
    assert len(P) == n
    h = berk_depth(P[:, 0], P[:, 1])
    sig = (P[:, 2] + h0) / h0
    P[:, 2] = -h + sig * h
    with open(pf, "w") as f:
        f.write(head + "(\n")
        f.writelines(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n" for p in P)
        f.write(")" + tail)


def berkhoff_nhwave(case, example_input, dx=0.05, dy=0.1, nLayers=3, total=30.0):
    """NHWAVE's example on the same (coarser) grid, depth from berk_depth,
    stations at the measurement sections."""
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    M, N = int(round(30.0 / dx)), int(round(20.0 / dy))
    x = (np.arange(M) + 0.5) * dx
    y = (np.arange(N) + 0.5) * dy
    X, Y = np.meshgrid(x, y)
    np.savetxt(os.path.join(case, "depth.txt"), berk_depth(X, Y), fmt="%.5f")
    pts = berk_points()
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for s, xs, ys in pts:
            f.write(f"{BERK_C[0] + xs:.4f} {BERK_C[1] + ys:.4f} 0.0\n")
    txt = open(example_input).read()
    txt = om.set_keys(txt, {"Mglob": M, "Nglob": N, "Kglob": nLayers, "DX": dx, "DY": dy,
                            "PX": 1, "PY": 1, "TOTAL_TIME": total, "NSTAT": len(pts),
                            "PLOT_INTV": 10.0, "PLOT_INTV_STAT": 0.02,
                            "SCREEN_INTV": 1.0})
    open(os.path.join(case, "input.txt"), "w").write(txt)
