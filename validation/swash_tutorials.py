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

usage (python):
    import swash_tutorials as st
    st.nhfoam_beach("nf_t1")                         # then ./Allrun
    st.nhfoam_beach("nf_t3", Tp=5.0, nz=2, endTime=900.0, tStat=300.0)
    st.swash_beach("sw_t1")                          # then swashrun -input wavbrk
"""
import os
import numpy as np
import make_case as mc

G = 9.81
L, D0, ZTOP = 1000.0, 10.0, 2.0          # length, offshore depth, height of the beach top


def depth(x):
    return D0 - (D0 + ZTOP) * np.asarray(x, float) / L


def nhfoam_beach(case, Hs=0.5, Tp=10.0, nz=1, dx=1.0, dt=0.05, endTime=1800.0,
                 tStat=600.0, film=0.01, manning=0.02, brk_model="HFA", alpha=0.6,
                 beta=0.3, nFrequencies=300, seed=1, gauge_dx=5.0,
                 conv="Gauss limitedLinearV 1", fMax=None):
    """Tutorials 1-3 in NHFOAM (2DV); fMax: highest frequency of the
    spectrum (default 3/Tp)."""
    xs, hs, zt = mc.beach_profile([0.0, L], [D0, -ZTOP], film, dx=dx)
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
    beta        {beta:g};
}}
""")
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(endTime / dt)), nCorr=3, nNonOrth=1,
                    pTol=1e-10, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dy, dx, nz, ztop=zt)
    mc.common_fields(case)
    fmax = f"\n        fMax        {fMax:g};" if fMax else ""
    mc.set_patch(case, "U", "left", mc.wave_bc("left", f"""theory      irregular;
        Hs          {Hs:g};
        Tp          {Tp:g};
        gamma       3.3;
        nFrequencies {nFrequencies};{fmax}
        seed        {seed};
        depth       {D0:g};
        direction   (1 0 0);
        rampTime    {2 * Tp:g};""", absorb=True))
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\nHs {Hs:g}\nTp {Tp:g}\nnz {nz}\ndx {dx:g}\ndt {dt:g}\n"
                f"endTime {endTime:g}\ntStat {tStat:g}\nfilm {film:g}\n")
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
    dep = depth(x)
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


def linear_shoaling(x, Hs, Tp, fMin=None, fMax=None, n=300):
    """Hs(x) of the JONSWAP spectrum that NHFOAM's wavemaker makes at
    x = 0 (n equal bins from fMin = 0.5/Tp to fMax = 3/Tp), shoaled by
    linear theory without dissipation."""
    fe = np.linspace(fMin or 0.5 / Tp, fMax or 3.0 / Tp, n + 1)
    f = 0.5 * (fe[1:] + fe[:-1])
    S = jonswap(f, Tp)
    S *= (Hs / 4.0) ** 2 / S.sum()
    out = []
    for h in depth(np.atleast_1d(x)):
        out.append(4.0 * np.sqrt(sum(S[j] * _cg(2 * np.pi * f[j], D0) / _cg(2 * np.pi * f[j], h)
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
