#!/usr/bin/env python3
"""The tutorials of SWASH (https://delftwaves.github.io/swash-docs, tutorials
1-4) in NHFOAM, with the SWASH input files next to them.

Tutorials 1-2 ("wavbrk"): JONSWAP waves, Hs = 0.5 m, Tp = 10 s, on a 1 km
long beach that rises from 10 m below to 2 m above the still water
(1:83.3), dx = 1 m, one layer, 30 min, dt = 0.05 s; output: the significant
wave height and the set-up over the last 20 min.  SWASH: weakly reflective
boundary with second-order bound waves (ADDBoundwave), BREAK (the
hydrostatic front approximation, alpha 0.6, beta 0.3), Manning n = 0.02,
upwind momentum.  Tutorial 3: Tp = 5 s, two layers, 15 min, statistics over
the last 10 min.

NHFOAM: the same profile, grid, layers and time step; the boundary
wavemaker nhWaveVelocity with a JONSWAP spectrum (linear: no bound waves)
and active absorption, the HFA breaking model with the same alpha and beta,
Manning n = 0.02, a film of 1 cm on the dry beach.  The statistics come
from gauges every 5 m.

SWASH's own command file of tutorials 1-3 sits in its Docker image; the one
written here is tutorial 4's listing (given in the documentation) brought
back to one dimension as tutorials 2-3 describe it.

Tutorial 4 ("obqwav"): short-crested waves at an angle on a beach that is
periodic alongshore (nhfoam_oblique, swash_oblique; see there), with the
wave components made as SWASH's source makes them (swash_components) and
linear theory for them (linear_oblique).

usage (python):
    import swash_tutorials as st
    st.nhfoam_beach("nf_t1")                         # then ./Allrun
    st.nhfoam_beach("nf_t3", Tp=5.0, nz=2, endTime=900.0, tStat=300.0)
    st.swash_beach("sw_t1")                          # then swashrun -input wavbrk
    st.nhfoam_oblique("nf_t4")                       # 4 h on two cores
"""
import os
import numpy as np
import make_case as mc

G = 9.81
L, D0, ZTOP = 1000.0, 10.0, 2.0          # length, offshore depth, height of the beach top


def depth(x, d0=D0):
    return d0 - (d0 + ZTOP) * np.asarray(x, float) / L


def nhfoam_beach(case, Hs=0.5, Tp=10.0, nz=1, dx=1.0, dt=0.05, endTime=1800.0,
                 tStat=600.0, film=0.01, manning=0.02, brk_model="HFA", alpha=0.6,
                 beta=0.3, nFrequencies=300, seed=1, gauge_dx=5.0,
                 conv="Gauss limitedLinearV 1", fMax=None, d0=D0, brk_extra="",
                 pTol=1e-10, layerAverage=False):
    """Tutorials 1-3 in NHFOAM (2DV); fMax: highest frequency of the
    spectrum (default 3/Tp); d0: depth at the wavemaker (23: the profile of
    tutorial 4); brk_extra: further entries of the breaking dictionary,
    e.g. "persistence crest;"; layerAverage: the wavemaker takes the
    average of the wave's velocity profile over each layer."""
    xs, hs, zt = mc.beach_profile([0.0, L], [d0, -ZTOP], film, dx=dx)
    extra = "\n    " + brk_extra if brk_extra else ""
    dy = dx
    gx = np.arange(0.5 * dx, L, gauge_dx)
    gtxt = "\n".join(f"    G{i:03d}    ({x:g} {dy/2:g} 0);" for i, x in enumerate(gx))
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

averaging
{{
    timeStart   {tStat:g};
}}

wettingDrying
{{
    minDepth    {film:g};
}}

bottomFriction
{{
    manning     {manning:g};
}}

breaking
{{
    model       {brk_model};
    alpha       {alpha:g};
    beta        {beta:g};{extra}
}}
""")
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(endTime / dt)), nCorr=3, nNonOrth=1,
                    pTol=pTol, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dy, dx, nz, ztop=zt)
    mc.common_fields(case)
    fmax = f"\n        fMax        {fMax:g};" if fMax else ""
    mc.set_patch(case, "U", "left", mc.wave_bc("left", f"""theory      irregular;
        Hs          {Hs:g};
        Tp          {Tp:g};
        gamma       3.3;
        nFrequencies {nFrequencies};{fmax}
        seed        {seed};
        depth       {d0:g};
        direction   (1 0 0);
        rampTime    {2 * Tp:g};""", absorb=True, layerAverage=layerAverage))
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\nHs {Hs:g}\nTp {Tp:g}\nnz {nz}\ndx {dx:g}\ndt {dt:g}\n"
                f"endTime {endTime:g}\ntStat {tStat:g}\nfilm {film:g}\nd0 {d0:g}\n")
    return gx


def swash_beach(case, Hs=0.5, Tp=10.0, nz=1, dx=1.0, dt=0.05, endTime=1800.0,
                tStat=600.0, cycle=None, addb=True):
    """The input files of tutorials 1-3 for SWASH: bathy.txt and bc.txt as
    the tutorials give them, wavbrk.sws from tutorial 4's listing."""
    os.makedirs(case, exist_ok=True)
    cyc = cycle if cycle else (endTime - tStat) / 60.0
    hh, mm = divmod(int(round(endTime / 60.0)), 60)
    with open(os.path.join(case, "bathy.txt"), "w") as f:
        f.write(f"{D0:g}. {-ZTOP:g}.\n")
    with open(os.path.join(case, "bc.txt"), "w") as f:
        f.write("BOUND SHAPE JONSWAP SIG PEAK\n"
                f"BOUND SIDE West BTYPE WEAK {'ADDBoundwave ' if addb else ''}CON SPECTrum "
                f"{Hs:g} {Tp:g}. cycle {cyc:g} min\n")
    n = int(round(L / dx))
    with open(os.path.join(case, "wavbrk.sws"), "w") as f:
        f.write(f"""!************* HEADING ********************************************
!
!  tutorial: breaking waves over a sloping beach (1D)
!
!*********** MODEL INPUT ******************************************

MODE DYN ONED

CGRID  0. 0. 0. {L:g}. 0. {n} 0
{'VERT ' + str(nz) if nz > 1 else '!VERT 1'}

INPGRID BOTTOM 0. 0. 0. 1 0 {L:g}. 1.
READINPUT BOTTOM 1. 'bathy.txt' 1 0 FREE

NONHYDROSTATIC
BREAK
FRIC MANNING 0.02

include 'bc.txt'

DISCRET UPW  MOM
DISCRET UPW WMOM H NONE

!************* OUTPUT REQUESTS ************************************

QUANT DIST HEXP 10.
QUANT HS DUR {(endTime - tStat) / 60.0:g} MIN

TABLE 'COMPGRID' NOHEAD 'output.txt' DIST BOTLEV HS SETUP

COMPUTE 000000.000 {dt:g} SEC {hh:02d}{mm:02d}00.000

STOP
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model swash\nHs {Hs:g}\nTp {Tp:g}\nnz {nz}\ndx {dx:g}\ndt {dt:g}\n"
                f"endTime {endTime:g}\ntStat {tStat:g}\n")


# --------------------------------------------------------------------------
# Results

def params(case):
    out = {}
    for l in open(os.path.join(case, "params.txt")):
        k, v = l.split(None, 1)
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v.strip()
    return out


def nhfoam_stats(case):
    """x, bed level (positive down, as SWASH's BOTLEV), Hs = 4 std(eta) and
    set-up = mean(eta) over tStat..endTime at the gauges; dry gauges (the
    surface never above bed + 2 films) get NaN."""
    import conical_island as ci
    pr = params(case)
    t, v = ci.read_nhfoam_gauges(case)
    names = sorted(v)
    x = np.array([float(l.split("(")[1].split()[0])
                  for l in open(os.path.join(case, "postProcessing", "gauges", "0", "eta.dat"))
                  if l.startswith("# G")])
    s = t >= pr["tStat"]
    hs = np.full(len(x), np.nan)
    su = np.full(len(x), np.nan)
    dep = depth(x, pr.get("d0", D0))
    for i, n in enumerate(names):
        e = v[n][s]
        if e.max() > -dep[i] + 2 * pr["film"]:
            # time-weighted (uniform steps here)
            hs[i] = 4.0 * e.std()
            su[i] = e.mean()
    return x, dep, hs, su


def swash_stats(case):
    """x, bed level, Hs and set-up of SWASH's output.txt (exception values
    -9 and -99 -> NaN)."""
    a = np.loadtxt(os.path.join(case, "output.txt"))
    hs = np.where(a[:, 2] < -8.9, np.nan, a[:, 2])
    su = np.where(a[:, 3] < -8.9, np.nan, a[:, 3])
    return a[:, 0], a[:, 1], hs, su


# --------------------------------------------------------------------------
# Linear shoaling of the spectrum (energy flux of each component conserved)

def _cg(om, h):
    k = om * om / G
    for _ in range(100):
        k = om * om / (G * np.tanh(k * h))
    return 0.5 * om / k * (1.0 + 2.0 * k * h / np.sinh(2.0 * k * h))


def jonswap(f, Tp, gamma=3.3):
    fp = 1.0 / Tp
    sig = np.where(f <= fp, 0.07, 0.09)
    return f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(
        -(f - fp) ** 2 / (2.0 * sig ** 2 * fp ** 2))


def linear_shoaling(x, Hs, Tp, fMin=None, fMax=None, n=300, d0=D0):
    """Hs(x) of the JONSWAP spectrum that NHFOAM's wavemaker makes at
    x = 0 (n equal bins from fMin = 0.5/Tp to fMax = 3/Tp), shoaled by
    linear theory without dissipation."""
    fe = np.linspace(fMin or 0.5 / Tp, fMax or 3.0 / Tp, n + 1)
    f = 0.5 * (fe[1:] + fe[:-1])
    S = jonswap(f, Tp)
    S *= (Hs / 4.0) ** 2 / S.sum()
    out = []
    for h in depth(np.atleast_1d(x), d0):
        out.append(4.0 * np.sqrt(sum(S[j] * _cg(2 * np.pi * f[j], d0) / _cg(2 * np.pi * f[j], h)
                                     for j in range(n))) if h > 0.05 else np.nan)
    return np.array(out)


def band_heights(case, xs, bands=((0.5, 1.5), (1.5, 1e9))):
    """4 std of eta at the gauges nearest to xs, in total and in the
    frequency bands (in units of the peak frequency)."""
    import conical_island as ci
    pr = params(case)
    t, v = ci.read_nhfoam_gauges(case)
    names = sorted(v)
    gx = np.array([float(l.split("(")[1].split()[0])
                   for l in open(os.path.join(case, "postProcessing", "gauges", "0", "eta.dat"))
                   if l.startswith("# G")])
    s = t >= pr["tStat"]
    dt = t[1] - t[0]
    out = []
    for x in xs:
        e = v[names[int(np.argmin(np.abs(gx - x)))]][s]
        e = e - e.mean()
        F = np.fft.rfftfreq(len(e), dt) * pr["Tp"]
        P = np.abs(np.fft.rfft(e)) ** 2
        tot = P[1:].sum()
        out.append([4 * e.std()] + [4 * e.std() * np.sqrt(P[(F >= a) & (F < b)].sum() / tot)
                                    for a, b in bands])
    return np.array(out)


# --------------------------------------------------------------------------
# Tutorial 4 ("obqwav"): short-crested waves at an angle to a plane beach,
# periodic alongshore; undertow and longshore current.  1000 m x 200 m, the
# depth falls from 23 m to -2 m (1:40, shoreline at x = 920 m), dx = 1 m,
# dy = 4 m, one layer, JONSWAP Hs = 2 m, Tp = 10 s, mean direction 15 deg,
# directional spreading 12 deg, Manning 0.02, Smagorinsky, 30 min with
# dt = 0.05 s, statistics over the last 20 min.

L4, W4, D4, ZTOP4 = 1000.0, 200.0, 23.0, 2.0


def depth4(x):
    return D4 - (D4 + ZTOP4) * np.asarray(x, float) / L4


def _k(om, h):
    k = om * om / G / np.sqrt(np.tanh(om * om * h / G))
    for _ in range(50):
        t = np.tanh(k * h)
        k -= (G * k * t - om * om) / (G * t + G * k * h * (1 - t * t))
    return k


def swash_components(Hs=2.0, Tp=10.0, direction=15.0, dspr=12.0, cycle=1200.0,
                     width=W4, depth=D4, seed=123456789, ndir=361, fit=True):
    """The wave components as SWASH's spectrum boundary makes them
    (SwashBCspectrum, SwashBCboundwave): frequencies every 1/cycle from
    fp/2 to 3 fp, the amplitudes of the JONSWAP spectrum (gamma 3.3), one
    direction for each frequency drawn from cos^m with m = 1.01/sigma^2 - 2
    (sigma = dspr in radians, below 17 deg) within a sector of 90 deg (m
    between 6 and 30), and each direction moved so that the alongshore wave
    number is a multiple of 2 pi/width (fit).  Returns frequency [Hz],
    amplitude [m], direction [deg from x, counter-clockwise] and phase
    [rad]; the random numbers are numpy's, not SWASH's."""
    df = 1.0 / cycle
    nfp = int(round(1.0 / Tp / df))
    fp = nfp * df
    fmin = (nfp // 2) * df
    n = int(round((3.0 * fp - fmin) / df))
    f = fmin + df * np.arange(n)
    S = jonswap(f, 1.0 / fp)
    a = np.sqrt(2.0 * S * df)
    a *= Hs / (4.0 * np.sqrt(0.5 * (a ** 2).sum()))
    sig = np.radians(dspr)
    m = max((1.2 if dspr > 23 else 1.096 if dspr > 17 else 1.01) / sig ** 2 - 2.0, 1.0)
    sector = np.pi * (1.0 if m <= 6 else 0.5 if m <= 30 else 0.25 if m <= 100 else 0.125)
    th = sector * (np.arange(ndir) / (ndir - 1.0) - 0.5)
    cdf = np.cumsum(np.maximum(np.cos(th), 0.0) ** m)
    cdf /= cdf[-1]
    rng = np.random.default_rng(seed)
    theta = np.radians(direction) + np.interp(rng.random(n), cdf, th)
    phase = 2.0 * np.pi * rng.random(n)
    if fit:
        k = _k(2 * np.pi * f, depth)
        dk = 2.0 * np.pi / width / k
        r = np.round(np.sin(theta) / dk) * dk
        r = np.where(r > 1.0, r - dk, np.where(r < -1.0, r + dk, r))
        theta = np.arcsin(r)
    return f, a, np.degrees(theta), phase


def linear_oblique(x, comps, depth0=D4):
    """Linear theory without dissipation for the components (f, a, theta,
    phase) made at the depth depth0: Hs(x) (energy flux towards the shore
    conserved, the alongshore wave number conserved), the mean
    depth-averaged velocities that return the mass flux of the waves,
    U = -sum a^2 c cos(theta)/(2 h^2) and the same with sin(theta) for V,
    and the radiation stress Sxy/rho [m3/s2] (constant)."""
    f, a, th0, _ = comps
    om = 2 * np.pi * f
    k0 = _k(om, depth0)
    ky = k0 * np.sin(np.radians(th0))
    cg0 = 0.5 * om / k0 * (1 + 2 * k0 * depth0 / np.sinh(2 * k0 * depth0))
    c0 = np.cos(np.radians(th0))
    hs, U, V = [], [], []
    for h in depth4(np.atleast_1d(x)):
        if h < 0.05:
            hs.append(np.nan); U.append(np.nan); V.append(np.nan)
            continue
        k = _k(om, h)
        s = np.clip(ky / k, -1, 1)
        c = np.sqrt(1 - s * s)
        cg = 0.5 * om / k * (1 + 2 * k * h / np.sinh(2 * k * h))
        a2 = a * a * cg0 * c0 / (cg * c)
        hs.append(4 * np.sqrt(0.5 * a2.sum()))
        U.append(-(a2 * om / k * c).sum() / (2 * h * h))
        V.append(-(a2 * om / k * s).sum() / (2 * h * h))
    sxy = (0.5 * G * a * a * cg0 / (om / k0) * np.sin(np.radians(th0)) * c0).sum()
    return np.array(hs), np.array(U), np.array(V), sxy


def nhfoam_oblique(case, Hs=2.0, Tp=10.0, direction=15.0, dspr=12.0, nz=1, dx=1.0,
                   ny=50, dt=0.05, endTime=1800.0, tStat=600.0, film=0.01,
                   manning=0.02, brk_model="HFA", alpha=0.6, beta=0.3,
                   persistence="crest", les_Ck=0.094,
                   seed=123456789, gauge_dx=5.0, writeEvery=200.0, nproc=2,
                   frames=(1500.0, 1800.0), conv="Gauss limitedLinearV 1", pTol=1e-8):
    """Tutorial 4 in NHFOAM: the grid, layer, time step, friction and
    breaking parameters of SWASH's listing; the wave components of
    swash_components at the boundary wavemaker (linear, with active
    absorption); breaking columns stay hydrostatic until the crest has
    passed, as in SWASH (persistence crest); the Smagorinsky model of the
    OpenFOAM library (les_Ck None: none); gauges along y = W/2 + dy/2
    (line 0) and y = W/4 (line 1).
    The run restarts from the last written time (every writeEvery seconds)
    when ./Allrun is called again."""
    dy = W4 / ny
    xs, hs, zt = mc.beach_profile([0.0, L4], [D4, -ZTOP4], film, dx=dx)
    comps = swash_components(Hs, Tp, direction, dspr, seed=seed)
    gx = np.arange(0.5 * dx, L4, gauge_dx)
    gl = []
    for j, (y, step) in enumerate(((0.5 * W4 + 0.5 * dy, 1), (0.25 * W4, 2))):
        gl += [f"    G{j}_{i:03d}    ({x:g} {y:g} 0);" for i, x in enumerate(gx[::step])]
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{chr(10).join(gl)}
}}

averaging
{{
    timeStart   {tStat:g};
}}

surfaceSnapshots ({" ".join(f"{t:g}" for t in frames)});

wettingDrying
{{
    minDepth    {film:g};
}}

bottomFriction
{{
    manning     {manning:g};
}}

breaking
{{
    model       {brk_model};
    alpha       {alpha:g};
    beta        {beta:g};
    persistence {persistence};
}}
""")
    np.savetxt(os.path.join(case, "components.txt"), np.c_[comps],
               header="frequency [Hz], amplitude [m], direction [deg], phase [rad]")
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(writeEvery / dt)), nCorr=3, nNonOrth=1,
                    pTol=pTol, conv=conv)
    cd = os.path.join(case, "system/controlDict")
    txt = open(cd).read().replace("startFrom       startTime;", "startFrom       latestTime;")
    txt = txt.replace("purgeWrite      0;", "purgeWrite      2;")
    open(cd, "w").write(txt)
    mc.blockmesh_profile(case, xs, hs, W4, dx, nz, Ny=ny, ztop=zt, cyclic=True)
    mc.common_fields(case, walls=("left", "right"), empty=False, constraint_types=True)
    if les_Ck:
        mc.les_files(case, Ck=les_Ck, patches={
            "surface": "zeroGradient", "bottom": "zeroGradient",
            "left": "zeroGradient", "right": "zeroGradient",
            "front": "cyclic", "back": "cyclic"})
    ctxt = "\n".join(f"            ({a:.8g} {1 / f:.10g} {th:.10g} {ph:.8g})"
                     for f, a, th, ph in zip(*comps))
    mc.set_patch(case, "U", "left", mc.wave_bc("left", f"""theory      irregular;
        depth       {D4:g};
        direction   (1 0 0);
        rampTime    {2 * Tp:g};
        components
        (
{ctxt}
        );""", absorb=True))
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
    mc.allrun(case, f"""
if [ ! -d processor0 ]
then
    runApplication blockMesh
    runApplication decomposePar
fi
mpirun -n {nproc} NHFOAM -parallel >> log.NHFOAM 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\nHs {Hs:g}\nTp {Tp:g}\ndirection {direction:g}\n"
                f"dspr {dspr:g}\nnz {nz}\ndx {dx:g}\ndy {dy:g}\ndt {dt:g}\n"
                f"endTime {endTime:g}\ntStat {tStat:g}\nfilm {film:g}\n")
    return comps


def swash_oblique(case):
    """Tutorial 4's own input (the listing of the documentation)."""
    os.makedirs(case, exist_ok=True)
    open(os.path.join(case, "bathy.txt"), "w").write("23. -2.\n23. -2.\n")
    open(os.path.join(case, "obqwav.sws"), "w").write("""PROJECT 'obqwav' '1'

MODE DYNANIC TWODIMENSIONAL

CGRID  0. 0. 0. 1000. 200. 1000 50 REPEAT Y

INPGRID BOTTOM 0. 0. 0. 1 1 1000. 200.
READINPUT BOTTOM 1. 'bathy.txt' 1 0 FREE

NONHYDROSTATIC
BREAK
FRIC MANNING 0.02
VISC SMAG

SET SEED 123456789
BOUND SHAPE JONSWAP SIG PEAK DSPR DEGRees
BOUND SIDE West CCW BTYPE WEAK ADDB CON SPECTrum h=2. per=10. dir=15. dd=12. cycle 20 min

DISCRET UPW  MOM
DISCRET UPW WMOM H NONE

QUANT DIST HEXP 10.
QUANT HS MVEL DUR 20 MIN

CURVE 'transect' 0. 100. 1000 1000. 100.
TABLE 'transect' NOHEAD 'output.txt' DIST BOTLEV HS SETUP MVEL
BLOCK 'COMPGRID' NOHEAD 'output.mat' LAY 3 XP YP BOTLEV MVEL

TEST 1 0
COMPUTE 000000.000 0.05 SEC 003000.000
STOP
""")


def oblique_stats(case, fSplit=0.5):
    """Along the gauge lines 0 (y = W/2 + dy/2) and 1 (y = W/4), over
    tStat..endTime: {line: (x, depth, Hs, Hs of the band above fSplit fp,
    Hs of the band below it, set-up)}, Hs = 4 std(eta); dry gauges (the
    surface never above bed + 2 films) get NaN."""
    import conical_island as ci
    pr = params(case)
    t, v = ci.read_nhfoam_gauges(case)
    root = os.path.join(case, "postProcessing", "gauges")
    fn = os.path.join(root, sorted(os.listdir(root), key=float)[0], "eta.dat")
    pos = {l.split()[1]: float(l.split("(")[1].split()[0])
           for l in open(fn) if l.startswith("# G")}
    s = t >= pr["tStat"]
    dt = t[1] - t[0]
    F = np.fft.rfftfreq(int(s.sum()), dt) * pr["Tp"]
    out = {}
    for line in (0, 1):
        names = sorted(n for n in v if n.startswith(f"G{line}_"))
        x = np.array([pos[n] for n in names])
        dep = depth4(x)
        res = np.full((4, len(x)), np.nan)
        for i, n in enumerate(names):
            e = v[n][s]
            if e.max() > -dep[i] + 2 * pr["film"]:
                P = np.abs(np.fft.rfft(e - e.mean())) ** 2
                tot = max(P[1:].sum(), 1e-300)
                hs = 4.0 * e.std()
                res[:, i] = (hs, hs * np.sqrt(P[F >= fSplit].sum() / tot),
                             hs * np.sqrt(P[(F > 0) & (F < fSplit)].sum() / tot), e.mean())
        out[line] = (x, dep) + tuple(res)
    return out


def oblique_means(case, time=None):
    """Alongshore averages of the time means: x, etaMean, U (cross-shore),
    V (alongshore) and the alongshore standard deviation of V."""
    import longshore as ls
    time = time or ls.latest_time(case)
    C = ls.cell_centres(case, time)
    U = ls._read_field(os.path.join(case, time, "UMean"))
    e = ls._read_field(os.path.join(case, time, "etaMean"))
    xr = np.round(C[:, 0], 6)
    xs, inv = np.unique(xr, return_inverse=True)
    cnt = np.bincount(inv)
    mean = lambda q: np.bincount(inv, q) / cnt
    v = mean(U[:, 1])
    sd = np.sqrt(np.maximum(mean(U[:, 1] ** 2) - v * v, 0))
    return xs, mean(e), mean(U[:, 0]), v, sd


if __name__ == "__main__":
    import sys
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    os.chdir(out)
    nhfoam_beach("nf_t1")
    nhfoam_beach("nf_t3", Tp=5.0, nz=2, endTime=900.0, tStat=300.0)
    # tutorial 3: grid and time step halved; no components above 1.5 fp;
    # five layers
    nhfoam_beach("nf_t3_fine", Tp=5.0, nz=2, endTime=900.0, tStat=300.0, dx=0.5, dt=0.025)
    nhfoam_beach("nf_t3_f15", Tp=5.0, nz=2, endTime=900.0, tStat=300.0, fMax=0.3,
                 nFrequencies=120)
    nhfoam_beach("nf_t3_L5", Tp=5.0, nz=5, endTime=900.0, tStat=300.0)
    swash_beach("sw_t1")
    swash_beach("sw_t3", Tp=5.0, nz=2, endTime=900.0, tStat=300.0)
    nhfoam_oblique("nf_t4")
    swash_oblique("sw_t4")
    # the profile and waves of tutorial 4 in 2DV, 10 min: the breaking model
    # as in NHFOAM's default, with SWASH's persistence, and all columns
    # hydrostatic
    kw = dict(Hs=2.0, Tp=10.0, d0=D4, endTime=600.0, tStat=200.0, pTol=1e-8)
    nhfoam_beach("t4_2dv", **kw)
    nhfoam_beach("t4_2dv_crest", brk_extra="persistence crest;", **kw)
    nhfoam_beach("t4_2dv_hyd", alpha=-1000, beta=-1000, **kw)
