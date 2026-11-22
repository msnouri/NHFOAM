"""Wetting-and-drying cases for NHFOAM (and NHWAVE).

    lake_at_rest(case)          still water against a dry beach (well-balanced)
    thacker(case, ...)          planar surface oscillating in a parabolic basin
                                (Thacker 1981; nonlinear shallow-water
                                solution with moving shorelines)
    thacker_exact(x, t, ...)    its analytical surface and shoreline
    dambreak(case, ...)         dam break on a dry (or wet) flat bed
    ritter_exact(x, t, ...)     Ritter's dry-bed solution
    runup_nhfoam(case, ...)     solitary wave on a 1:20 beach (Synolakis 1987,
    runup_nhwave(case, ...)     NTHMP benchmark 4; Ma et al. 2012, 4.5), on
                                NHWAVE's example solitary_wave_breaking
    tingkirby(case, ...)        spilling breakers on a 1:35 beach (Ting &
                                Kirby 1994), stream-function waves

Options of the wettingDrying dictionary can be passed as wd_extra, e.g.
wd_extra="\n    faceDepth upwind;".  Example:

    import wetdry
    wetdry.runup_nhfoam("ru", Cd=0.0)        # then ./Allrun in ru
"""
import os
import numpy as np
import make_case as mc

G = 9.81


def profile_cell_x(xs, dx, nz):
    """x of the cell centres of a blockmesh_profile mesh, in mesh order
    (block by block, x fastest, then the layers)."""
    out = []
    for i in range(len(xs) - 1):
        nx = max(1, int(round((xs[i + 1] - xs[i]) / dx)))
        xc = xs[i] + (np.arange(nx) + 0.5) * (xs[i + 1] - xs[i]) / nx
        out.append(np.tile(xc, nz))
    return np.concatenate(out)


def write_scalar_field(case, name, values, dims="[0 1 0 0 0 0 0]", patches=None):
    vals = "\n".join(f"{v:.10g}" for v in values)
    bf = patches or "\".*\"  { type zeroGradient; }\n    frontAndBack { type empty; }"
    mc.write(case, f"0/{name}", "volScalarField", f"""
dimensions      {dims};
internalField   nonuniform List<scalar>
{len(values)}
(
{vals}
);
boundaryField
{{
    {bf}
}}
""")


# --------------------------------------------------------------------------
def lake_at_rest(case, film=0.005, dx=0.05, nz=3, endTime=10.0):
    xs, hs, zt = mc.beach_profile([0, 5, 15], [0.4, 0.4, -0.1], film)
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film};
}}
""")
    mc.system_files(case, endTime=endTime, deltaT=0.01, writeInterval=500,
                    nCorr=3, nNonOrth=1, pTol=1e-10, conv="Gauss limitedLinearV 1")
    mc.blockmesh_profile(case, xs, hs, dx, dx, nz, ztop=zt)
    mc.common_fields(case)
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")


# --------------------------------------------------------------------------
# Thacker (1981), 1D: bed z = h0 ((x - L/2)^2/a^2 - 1); the free surface stays
# planar, eta = -h0 beta (2 (x - L/2)/a + beta cos(w t)) cos(w t), with
# w = sqrt(2 g h0)/a, beta = B/(a w) and u = B sin(w t) in the wet region,
# whose edges are x - L/2 = a (-beta cos(w t) -+ 1).

def thacker_exact(x, t, a=100.0, h0=0.5, B=0.5, L=400.0):
    w = np.sqrt(2 * G * h0) / a
    beta = B / (a * w)
    xi = (x - L / 2) / a
    c = np.cos(w * t)
    zb = h0 * (xi ** 2 - 1)
    eta = -h0 * beta * c * (2 * xi + beta * c)
    wet = eta > zb
    shore = (L / 2 + a * (-beta * c - 1), L / 2 + a * (-beta * c + 1))
    return np.where(wet, eta, np.nan), zb, shore, 2 * np.pi / w


def thacker(case, a=100.0, h0=0.5, B=0.5, L=400.0, dx=1.0, nz=3, dt=0.25,
            periods=2.0, film=0.001, nseg=80, conv="Gauss limitedLinearV 1",
            wd_extra=""):
    """Parabolic basin approximated by nseg linear segments."""
    xs0 = np.linspace(0, L, nseg + 1)
    zb0 = h0 * (((xs0 - L / 2) / a) ** 2 - 1)
    xs, hs, zt = mc.beach_profile(list(xs0), list(-zb0), film, dx=dx)
    T = 2 * np.pi * a / np.sqrt(2 * G * h0)
    endTime = periods * T
    nsteps = int(round(endTime / dt))
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film};{wd_extra}
}}
""")
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(T / 8 / dt)), nCorr=3, nNonOrth=1,
                    pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dx, dx, nz, ztop=zt)
    mc.common_fields(case)
    xc = profile_cell_x(xs, dx, nz)
    eta0, _, _, _ = thacker_exact(xc, 0.0, a, h0, B, L)
    write_scalar_field(case, "eta", np.nan_to_num(eta0, nan=-10.0))
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    return T


# --------------------------------------------------------------------------
# Dam break on a dry, flat, frictionless bed (Ritter 1892): for |x - x0| in
# [-c0 t, 2 c0 t], c0 = sqrt(g h0): h = (2 c0 - (x - x0)/t)^2/(9 g),
# u = 2/3 ((x - x0)/t + c0); the front runs at 2 c0.

def ritter_exact(x, t, h0=0.1, x0=5.0):
    c0 = np.sqrt(G * h0)
    xi = (x - x0) / max(t, 1e-12)
    h = np.where(xi < -c0, h0, (2 * c0 - xi) ** 2 / (9 * G))
    u = np.where(xi < -c0, 0.0, 2.0 / 3.0 * (xi + c0))
    dry = xi > 2 * c0
    return np.where(dry, 0.0, h), np.where(dry, 0.0, u)


def dambreak(case, h0=0.1, L=15.0, x0=5.0, dx=0.02, nz=3, dt=0.004, tend=3.0,
             film=0.001, hydrostatic=True, conv="Gauss limitedLinearV 1",
             wd_extra="", writeEvery=0.5, hR=None):
    """hR: depth downstream (default: dry, the film)"""
    xs, hs, zt = [0.0, L], [h0, h0], [0.0, 0.0]
    hyd = "\n    hydrostaticDepth 100;" if hydrostatic else ""
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film};{hyd}{wd_extra}
}}
""")
    mc.system_files(case, endTime=tend, deltaT=dt,
                    writeInterval=int(round(writeEvery / dt)), nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dx, dx, nz, ztop=zt)
    mc.common_fields(case)
    xc = profile_cell_x(xs, dx, nz)
    write_scalar_field(case, "eta", np.where(xc < x0, 0.0,
                                             -10.0 if hR is None else hR - h0))
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")


# --------------------------------------------------------------------------
# Ting & Kirby (1994), spilling breakers: T = 2 s, H = 0.125 m in h = 0.4 m,
# 1:35 slope.  x = x_TK, measured from the toe (the slope starts at the origin,
# which reproduces the measured breaking point; the still-water shoreline is
# at x = 14 m).  Stream-function waves (Fenton 1988, zero mass flux) are made
# in a relaxation zone one wavelength long; the beach runs up onto dry land.

TK_H, TK_T, TK_D, TK_SLOPE = 0.125, 2.0, 0.4, 1.0 / 35.0
TK_STATIONS = [6.665, 7.275, 7.885, 8.495, 9.110, 9.725, 10.335, 10.945]


def tk_depth(x):
    return np.where(x < 0, TK_D, TK_D - TK_SLOPE * np.asarray(x, float))


def tingkirby(case, nz=10, dx=0.025, dt=0.005, nPeriods=30, zone=3.8,
              flat=2.2, xEnd=17.0, film=0.005, Cd=0.006, breaking=True,
              conv="Gauss limitedLinearV 1", wd_extra="", N=24,
              theory="streamFunction", alpha=0.6, beta=0.3, gauge_dx=0.1,
              brk_model="HFA", brk_extra="", les_Ck=None):
    from stream_function import fenton, foam_entries
    x0 = -(zone + flat)
    xs0 = [x0, 0.0, xEnd]
    hs0 = [TK_D, TK_D, TK_D - TK_SLOPE * xEnd]
    xs, hs, zt = mc.beach_profile(xs0, hs0, film, dx=dx)
    dy = dx
    gx = sorted(set([round(x, 3) for x in np.arange(-3.0, 13.0 + 1e-9, gauge_dx)]
                    + TK_STATIONS))
    gtxt = "\n".join(f"    G{i:03d}    ({x:g} {dy/2:g} 0);" for i, x in enumerate(gx))
    if theory == "streamFunction":
        sol = fenton(TK_H, TK_T, TK_D, N=N)
        wave = foam_entries(sol, indent="        ")
    else:
        wave = (f"        waveHeight  {TK_H};\n        wavePeriod  {TK_T};\n"
                f"        depth       {TK_D};\n")
    brk = f"""
breaking
{{
    model       {brk_model};
    alpha       {alpha:g};
    beta        {beta:g};{brk_extra}
}}
""" if breaking else ""
    fric = f"""
bottomFriction
{{
    Cd          {Cd};
}}
""" if Cd > 0 else ""
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
        start       ({x0 + zone:g} 0 0);
        end         ({x0:g} 0 0);
        sigmaMax    20;
        exponent    2;
{wave}        direction   (1 0 0);
        origin      ({x0:g} 0 0);
        rampTime    {2 * TK_T:g};
    }}
}}

wettingDrying
{{
    minDepth    {film};{wd_extra}
}}
{fric}{brk}""")
    endTime = nPeriods * TK_T
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(endTime / dt)), nCorr=3, nNonOrth=1,
                    pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dy, dx, nz, ztop=zt)
    mc.common_fields(case)
    if les_Ck:
        mc.les_files(case, Ck=les_Ck)
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"T {TK_T}\nH {TK_H}\nnz {nz}\ndx {dx}\n")
    return gx


# --------------------------------------------------------------------------
# Solitary wave on a plane beach (Synolakis 1987; NTHMP benchmark 4; Ma et
# al. 2012, 4.5), on the set-up of NHWAVE's examples/solitary_wave_breaking:
# d = 0.21 m, toe at x = 4 m, slope 1:20 to x = 11 m (bed 0.14 m above the
# still water level), dx = 2 cm, 3 layers, MinDep 5 mm, Cd 0.006.  Both
# models start from the exact solitary wave (Babenko's equation) with its
# crest L = d acosh(sqrt 20)/gamma from the toe, gamma = sqrt(3H/4d), which
# is t = 0 of the benchmark.

RU_D, RU_TOE, RU_SLOPE, RU_LEN, RU_DX = 0.21, 4.0, 1.0 / 20.0, 11.0, 0.02


def runup_geometry(H, d=RU_D):
    gamma = np.sqrt(0.75 * H / d)
    L = d * np.arccosh(np.sqrt(20.0)) / gamma
    xshore = RU_TOE + d / RU_SLOPE
    return RU_TOE - L, xshore


def runup_depth(x, d=RU_D):
    return np.where(x < RU_TOE, d, d - (x - RU_TOE) * RU_SLOPE)


def runup_nhfoam(case, H=0.0588, d=RU_D, nz=3, dt=0.002, tstar_end=60.0,
                 film=0.005, Cd=0.006, breaking=True,
                 conv="Gauss limitedLinearV 1", wd_extra=""):
    from solitary_exact import SolitaryWave
    import ma2012
    xc0, xshore = runup_geometry(H, d)
    xs, hs, zt = mc.beach_profile([0.0, RU_TOE, RU_LEN],
                                  [d, d, d - (RU_LEN - RU_TOE) * RU_SLOPE], film)
    tscale = np.sqrt(d / G)
    brk = """
breaking
{
    model       HFA;
    alpha       0.6;
    beta        0.3;
}
""" if breaking else ""
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film};{wd_extra}
}}

bottomFriction
{{
    Cd          {Cd};
}}
{brk}""")
    mc.system_files(case, endTime=tstar_end * tscale, deltaT=dt, writeInterval=1,
                    nCorr=3, nNonOrth=1, pTol=1e-10, conv=conv)
    cd = open(os.path.join(case, "system/controlDict")).read()
    import re
    cd = re.sub(r"writeControl\s+\w+;", "writeControl    runTime;", cd)
    cd = re.sub(r"writeInterval\s+[^;]+;", f"writeInterval   {5 * tscale:.8g};", cd)
    open(os.path.join(case, "system/controlDict"), "w").write(cd)
    mc.blockmesh_profile(case, xs, hs, RU_DX, RU_DX, nz, ztop=zt)
    mc.common_fields(case)
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    # initial fields: exact solitary wave, crest at xc0, on the cell centres
    sw = SolitaryWave(H, d=d)
    x = profile_cell_x(xs, RU_DX, nz)
    ncell = len(x)
    # layer index of each cell (blocks: x fastest, then the layers)
    k = np.zeros(ncell, int)
    pos = 0
    for i in range(len(xs) - 1):
        nx = max(1, int(round((xs[i + 1] - xs[i]) / RU_DX)))
        k[pos:pos + nx * nz] = np.repeat(np.arange(nz), nx)
        pos += nx * nz
    h = runup_depth(x, d)
    eta = sw.eta(x - xc0)
    wet = eta + h > film
    D = np.where(wet, h + eta, film)
    zc = -h + (k + 0.5) / nz * D
    u, w = sw.velocity(x - xc0, zc)
    u = np.where(wet, u, 0.0); w = np.where(wet, w, 0.0)
    p_rgh = np.where(wet, 0.5 * (sw.c ** 2 - (u - sw.c) ** 2 - w ** 2), G * (film - h))
    etaw = np.where(wet, eta, -h + film)
    ma2012._set_internal(os.path.join(case, "0", "eta"), ma2012._foam_list(etaw))
    ma2012._set_internal(os.path.join(case, "0", "U"), ma2012._foam_list(list(zip(u, w)), vector=True))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"), ma2012._foam_list(p_rgh))
    return dict(xc0=xc0, xshore=xshore, tscale=tscale, c=sw.c)


def runup_nhwave(case, template_dir, H=0.0588, d=RU_D, nz=3, tstar_end=60.0,
                 cfl=0.5):
    """NHWAVE's example with the exact initial condition instead of its
    boundary wavemaker (as in Ma et al. 2012), output every 5 t*."""
    from solitary_exact import SolitaryWave
    import other_models as om
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    xc0, xshore = runup_geometry(H, d)
    tscale = np.sqrt(d / G)
    Nx = int(round(RU_LEN / RU_DX))
    x = (np.arange(Nx) + 0.5) * RU_DX
    h = runup_depth(x, d)
    np.savetxt(os.path.join(case, "depth.txt"), h[None, :], fmt="%.6f")
    txt = open(os.path.join(template_dir, "input.txt")).read()
    txt = om.set_keys(txt, {"Mglob": Nx, "Nglob": 1, "Kglob": nz, "DX": RU_DX,
                            "PX": 1, "PY": 1, "TOTAL_TIME": tstar_end * tscale,
                            "PLOT_INTV": 5 * tscale, "WAVEMAKER": "LEF_SOL",
                            "AMP": 0.0, "BC_X0": 1, "INITIAL_EUVW": "T",
                            "CFL": cfl, "NSTAT": 0, "OUT_E": "T", "OUT_U": "T",
                            "OUT_D": "T", "SCREEN_INTV": 1.0})
    open(os.path.join(case, "input.txt"), "w").write(txt)
    sw = SolitaryWave(H, d=d)
    eta = sw.eta(x - xc0)
    eta = np.where(eta + h > 0, eta, -h)
    sig = (np.arange(nz) + 0.5) / nz
    U = np.zeros((nz, Nx)); W = np.zeros((nz, Nx))
    for kk, s in enumerate(sig):
        zc = -h + s * np.maximum(h + eta, 0)
        U[kk], W[kk] = sw.velocity(x - xc0, zc)
    wet = eta + h > 0
    U[:, ~wet] = 0; W[:, ~wet] = 0
    np.savetxt(os.path.join(case, "eta0.txt"), eta[None, :], fmt="%.12e")
    np.savetxt(os.path.join(case, "uvw0.txt"),
               np.vstack([U, np.zeros_like(U), W]), fmt="%.12e")
    return dict(xc0=xc0, xshore=xshore, tscale=tscale)
