#!/usr/bin/env python3
"""Two flume tests with irregular breaking waves that Zijlema et al. (2011)
compare SWASH with, in NHFOAM: the bar-trough beach of Battjes & Janssen
(1978), runs 13 and 15, and the barred beach of Boers (1996), case 1B.

The measurements, the bed profiles and SWASH's own results are read from
the figures of the paper by digitize_zijlema2011.py into a directory
(DATA below, or the environment variable ZIJLEMA2011_DATA); nothing of
that is in the repository.  The still-water depths of Battjes & Janssen's
runs are from their Table 1: the bed of Fig. 14 of Zijlema et al. is
drawn for 0.70 m at the wavemaker (their run 11), and runs 13 and 15 had
0.762 and 0.616 m, which gives the depths over the bar crest that they
list (0.267 and 0.120 m).

SWASH's set-up in the paper: depth-averaged, dx = 0.05 m (Battjes &
Janssen) and 0.02 m (Boers), no bottom friction, and a mixing length
calibrated for each case (0.15, 0.20 and 0.13 m) instead of a breaking
criterion.  NHFOAM has no mixing-length viscosity: the runs here use the
hydrostatic front approximation with SWASH's later defaults
(brk="hfa": alpha 0.6, beta 0.3, persistence crest) or NHFOAM's breaking
eddy viscosity (brk="hybrid", the settings of the Ting & Kirby test).

    import swash_flumes as sf
    sf.battjes_janssen("bj15_L1", 15)                    # then ./Allrun
    sf.boers("b1b_L2", nz=2)
    sf.stats("bj15_L1")        # x, depth, Hm0, Tm02, mean level; cached
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_case as mc

DATA = os.environ.get("ZIJLEMA2011_DATA", "/home/claude/data_z11")
G = 9.81

# Battjes & Janssen (1978), Table 1: still-water depth at the wavemaker;
# Zijlema et al. (2011): Hm0 and Tp of the incident JONSWAP spectrum
BJ = {13: dict(d1=0.762, Hm0=0.147, Tp=2.01), 15: dict(d1=0.616, Hm0=0.202, Tp=1.89)}
BJ_SKETCH_DEPTH = 0.70
BOERS_1B = dict(Hm0=0.206, Tp=2.03)

BREAKING = {
    "hfa": "model       HFA;\n    alpha       0.6;\n    beta        0.3;\n"
           "    persistence crest;",
    "hybrid": "model       hybrid;\n    alpha       0.6;\n    beta        0.3;\n"
              "    betaFinal   0.15;\n    Tstar       5;\n    delta       1.2;",
    # SWASH's BREAK command as its source has it: the front approximation
    # and, in the breaking columns, a mixing length equal to the depth
    "swash": "model       hybrid;\n    alpha       0.6;\n    beta        0.3;\n"
             "    persistence crest;\n    viscosity   depth;\n    nufac       1;",
    "swash_limit": "model       hybrid;\n    alpha       0.6;\n    beta        0.3;\n"
                   "    persistence crest;\n    viscosity   depth;\n    nufac       1;\n"
                   "    explicitLimit yes;",
    # the whole broken wave hydrostatic, without and with the mixing length
    "wave": "model       HFA;\n    alpha       0.6;\n    beta        0.3;\n"
            "    persistence wave;",
    "wave_visc": "model       hybrid;\n    alpha       0.6;\n    beta        0.3;\n"
                 "    persistence wave;\n    viscosity   depth;\n    nufac       1;",
    "none": None,
}


def load(name):
    return np.atleast_2d(np.loadtxt(os.path.join(DATA, name)))


def _k(om, h):
    k = om * om / G / np.sqrt(np.tanh(om * om * h / G))
    for _ in range(50):
        t = np.tanh(k * h)
        k -= (G * k * t - om * om) / (G * t + G * k * h * (1 - t * t))
    return k


def _cg(om, h):
    k = _k(om, h)
    return 0.5 * om / k * (1 + 2 * k * h / np.sinh(2 * k * h))


def bj_profile(run, above=0.2):
    """x and still-water depth of the bar-trough beach for a run, the last
    slope extended to `above` over the still water."""
    p = load("bj78_profile.txt")
    x, z = p[:, 0].copy(), p[:, 1] - (BJ[run]["d1"] - BJ_SKETCH_DEPTH)
    x[0] = 0.0
    s = (z[-1] - z[-2]) / (x[-1] - x[-2])
    if z[-1] < above:
        x = np.append(x, x[-1] + (above - z[-1]) / s)
        z = np.append(z, above)
    return x, -z


def boers_profile():
    p = load("boers_profile.txt")
    x = p[:, 0].copy()
    x[0] = 0.0
    return x, -p[:, 1]


def spectrum_components(f, S, Hm0, depth_from, depth_to, T, fmin, fmax, seed=1):
    """Wave components (amplitude, period, direction, phase) with the shape
    of the spectrum S(f) given at the depth depth_from, taken to depth_to
    with linear shoaling, scaled to Hm0, at the frequencies n/T between
    fmin and fmax (so that the record repeats after T), random phases."""
    df = 1.0 / T
    fi = df * np.arange(int(np.ceil(fmin / df)), int(np.floor(fmax / df)) + 1)
    Si = np.maximum(np.interp(fi, f, S), 0.0)
    om = 2 * np.pi * fi
    a = np.sqrt(2 * Si * df * _cg(om, depth_from) / _cg(om, depth_to))
    a *= Hm0 / (4 * np.sqrt(0.5 * (a * a).sum()))
    ph = 2 * np.pi * np.random.default_rng(seed).random(len(fi))
    return np.c_[a, 1.0 / fi, np.zeros(len(fi)), ph]


def flume(case, x, depth, wave, d0, nz=1, dx=0.05, dt=0.01, endTime=1260.0, tStat=60.0,
          film=0.005, brk="hfa", Cd=None, stations=(), gauge_dx=0.5,
          layerAverage=True, conv="Gauss limitedLinearV 1", pTol=1e-10, meta=""):
    """A 2DV flume with the bed (x, depth), a boundary wavemaker with active
    absorption at x = 0 (entries `wave` of the target wave) and gauges every
    gauge_dx and at the stations."""
    xs, hs, zt = mc.beach_profile(list(x), list(depth), film, dx=dx)
    dy = dx
    xw = np.interp(film, depth[::-1], x[::-1]) if depth[-1] < film else x[-1]
    gx = sorted(set(np.round(np.r_[np.arange(0.25, min(xw + 1.0, x[-1]), gauge_dx),
                                   np.asarray(stations, float)], 3)))
    gtxt = "\n".join(f"    G{i:03d}    ({v:g} {dy/2:g} 0);" for i, v in enumerate(gx))
    extra = f"""
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
"""
    if Cd:
        extra += f"\nbottomFriction\n{{\n    Cd          {Cd:g};\n}}\n"
    if BREAKING[brk]:
        extra += f"\nbreaking\n{{\n    {BREAKING[brk]}\n}}\n"
    mc.constant_files(case, nh_extra=extra)
    mc.system_files(case, endTime=endTime, deltaT=dt,
                    writeInterval=int(round(endTime / dt)), nCorr=3, nNonOrth=1,
                    pTol=pTol, conv=conv)
    mc.blockmesh_profile(case, xs, hs, dy, dx, nz, ztop=zt)
    mc.common_fields(case)
    mc.set_patch(case, "U", "left", mc.wave_bc("left", wave + f"""
        depth       {d0:g};
        direction   (1 0 0);""", absorb=True, layerAverage=layerAverage))
    mc.allrun(case, """
runApplication blockMesh
runApplication NHFOAM
""")
    np.savetxt(os.path.join(case, "bed.txt"), np.c_[x, depth], header="x [m], still-water depth [m]")
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\nnz {nz}\ndx {dx:g}\ndt {dt:g}\nendTime {endTime:g}\n"
                f"tStat {tStat:g}\nfilm {film:g}\nbrk {brk}\n{meta}")


def battjes_janssen(case, run, nz=1, dx=0.05, dt=0.01, endTime=1260.0, tStat=60.0,
                    brk="hfa", nFrequencies=1500, seed=1, **kw):
    """Run 13 (mildly breaking) or 15 (strongly breaking): JONSWAP waves
    (gamma 3.3) with the Hm0 and Tp of Zijlema et al.; 1500 frequencies
    from 0.5 to 3 fp, so the record repeats after 1200 s."""
    r = BJ[run]
    x, depth = bj_profile(run)
    g = load("bj78_gauges.txt")[:, 0]
    lab = load(f"bj78_run{run}.txt")[:, 0]
    wave = f"""theory      irregular;
        Hs          {r['Hm0']:g};
        Tp          {r['Tp']:g};
        gamma       3.3;
        nFrequencies {nFrequencies};
        seed        {seed};
        rampTime    {2 * r['Tp']:g};"""
    flume(case, x, depth, wave, r["d1"], nz=nz, dx=dx, dt=dt, endTime=endTime, tStat=tStat,
          brk=brk, stations=np.r_[g, lab], meta=f"test bj{run}\nTp {r['Tp']:g}\n", **kw)


def boers(case, nz=1, dx=0.04, dt=0.008, endTime=1680.0, tStat=60.0, brk="hfa",
          fmin=0.2, fmax=1.5, seed=1, **kw):
    """Case 1B: the spectrum observed at station 1 (x = 5 m, Fig. 18 of
    Zijlema et al.), taken back to the depth at x = 0 with linear shoaling
    and scaled to the Hm0 measured there (0.206 m), as components that
    repeat after endTime - tStat.  SWASH: dx = 0.02 m, 1680 s."""
    x, depth = boers_profile()
    st = load("boers_gauges.txt")[:, 0]
    sp = load("boers1b_spec_1.txt")
    d1 = float(np.interp(st[0], x, depth))
    comps = spectrum_components(sp[:, 0], sp[:, 1], BOERS_1B["Hm0"], d1, depth[0],
                                endTime - tStat, fmin, fmax, seed)
    ctxt = "\n".join(f"            ({a:.8g} {T:.10g} 0 {ph:.8g})" for a, T, _, ph in comps)
    wave = f"""theory      irregular;
        rampTime    {2 * BOERS_1B['Tp']:g};
        components
        (
{ctxt}
        );"""
    flume(case, x, depth, wave, depth[0], nz=nz, dx=dx, dt=dt, endTime=endTime, tStat=tStat,
          brk=brk, stations=st,
          meta=f"test boers1b\nTp {BOERS_1B['Tp']:g}\n", **kw)
    np.savetxt(os.path.join(case, "components.txt"), comps,
               header="amplitude [m], period [s], direction [deg], phase [rad]")


def params(case):
    out = {}
    for line in open(os.path.join(case, "params.txt")):
        k, v = line.split(None, 1)
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v.strip()
    return out


def moments(t, e, fmin, fmax):
    """m0 and m2 of a record over fmin..fmax (one-sided periodogram)."""
    dt = t[1] - t[0]
    P = 2.0 * np.abs(np.fft.rfft(e - e.mean())) ** 2 / len(e) ** 2
    f = np.fft.rfftfreq(len(e), dt)
    s = (f >= fmin) & (f <= fmax) & (f > 0)
    return P[s].sum(), (P[s] * f[s] ** 2).sum()


def stats(case, band=(0.0, 1e9), spectra_at=(), nseg=16):
    """Over tStat..endTime at every gauge: x, still-water depth, Hm0 =
    4 sqrt(m0), Tm02 = sqrt(m0/m2) over the band [Hz], mean level, and the
    part of Hm0 below half the peak frequency; dry gauges NaN.  With
    spectra_at, also the spectra there (nseg half-overlapping Hann
    segments).  The result is kept in <case>/stats.npz, so the gauge
    record can be deleted afterwards."""
    key = f"{band[0]:g}_{band[1]:g}"
    cache = os.path.join(case, f"stats_{key}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    import conical_island as ci
    pr = params(case)
    t, v = ci.read_nhfoam_gauges(case)
    root = os.path.join(case, "postProcessing", "gauges")
    fn = os.path.join(root, sorted(os.listdir(root), key=float)[0], "eta.dat")
    pos = {l.split()[1]: float(l.split("(")[1].split()[0])
           for l in open(fn) if l.startswith("# G")}
    names = sorted(v, key=lambda n: pos[n])
    x = np.array([pos[n] for n in names])
    bed = np.loadtxt(os.path.join(case, "bed.txt"))
    dep = np.interp(x, bed[:, 0], bed[:, 1])
    s = t >= pr["tStat"]
    ts = t[s]
    n = len(x)
    hm0, tm02, eta, hlow = (np.full(n, np.nan) for _ in range(4))
    fp = 1.0 / pr["Tp"]
    for i, name in enumerate(names):
        e = v[name][s]
        if e.max() > -dep[i] + 2 * pr["film"]:
            m0, m2 = moments(ts, e, *band)
            hm0[i], tm02[i], eta[i] = 4 * np.sqrt(m0), np.sqrt(m0 / m2), e.mean()
            hlow[i] = 4 * np.sqrt(moments(ts, e, 0.0, 0.5 * fp)[0])
    out = dict(x=x, depth=dep, hm0=hm0, tm02=tm02, setup=eta, hlow=hlow)
    if len(spectra_at):
        from scipy.signal import welch
        dt = ts[1] - ts[0]
        nper = int(len(ts) / (nseg / 2 + 0.5))
        for k, xs in enumerate(spectra_at):
            i = int(np.argmin(np.abs(x - xs)))
            f, S = welch(v[names[i]][s], fs=1.0 / dt, nperseg=nper, window="hann")
            out[f"spec_f"] = f[f <= 2.5]
            out[f"spec_{k + 1}"] = S[f <= 2.5]
    np.savez(cache, **out)
    return out


BAND = (0.0, 2.0)       # the range of the spectra in the paper [Hz]


def finish(case, keep=False):
    """Statistics of a finished run (kept in the case), with the spectra at
    the stations of the experiment; then the gauge record is deleted."""
    import shutil
    pr = params(case)
    st = load("boers_gauges.txt" if pr["test"] == "boers1b" else "bj78_gauges.txt")[:, 0]
    out = stats(case, BAND, spectra_at=st)
    if not keep:
        shutil.rmtree(os.path.join(case, "postProcessing"), ignore_errors=True)
    return out


RUNS = {
    # case: (builder, arguments); L1: one layer with SWASH's front rule
    "bj13_L1": ("bj", dict(run=13)),
    "bj13_L1_wv": ("bj", dict(run=13, brk="wave_visc")),
    "bj13_L2_wv": ("bj", dict(run=13, nz=2, brk="wave_visc")),
    "bj15_L1": ("bj", dict(run=15)),
    "bj15_L1_wv": ("bj", dict(run=15, brk="wave_visc")),
    "bj15_L2_wv": ("bj", dict(run=15, nz=2, brk="wave_visc")),
    "bj15_L2_hyb": ("bj", dict(run=15, nz=2, brk="hybrid")),
    "b1b_L1": ("boers", dict()),
    "b1b_L1_wv": ("boers", dict(brk="wave_visc")),
    "b1b_L2_wv": ("boers", dict(nz=2, brk="wave_visc")),
}


if __name__ == "__main__":
    if sys.argv[1:2] == ["finish"]:
        for name in sys.argv[2:]:
            finish(name)
    else:
        for name in sys.argv[1:] or RUNS:
            kind, kw = RUNS[name]
            (battjes_janssen if kind == "bj" else boers)(name, **kw)
