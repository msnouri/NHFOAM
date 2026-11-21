#!/usr/bin/env python3
"""Input files for NHWAVE and REEF3D::NHFLOW on the NHFOAM validation cases,
with the same geometry, gauges and (where possible) resolution, so that their
results can be plotted with ours.

NHWAVE:  github.com/JimKirby/NHWAVE (built with -DPARALLEL -DDOUBLE_PRECISION
         -DFROUDE_CAP, plus -DLANDSLIDE for the slide); the input template is
         the NHWAVE example examples/submerged_bar/input.txt.
REEF3D:  REEF3D::NHFLOW 26.08 (conda-forge packages reef3d, divemesh).

usage (python):
    import other_models as om
    om.nhwave_bar("nw_barA", template, H=0.02, T=2.02)
    om.reef3d_bar("r3_barA", H=0.02, T=2.02)
"""
import os, re, math
import numpy as np

BAR_GAUGES = [22, 24, 30.5, 32.5, 33.5, 34.5, 35.7, 37.3, 39.0, 41.0]


def bar_depth(x):
    """Luth et al. (1994) bar in the NHFOAM/FUNWAVE coordinates."""
    xs = [0, 26, 32, 34, 37, 54]
    hs = [0.4, 0.4, 0.1, 0.1, 0.4, 0.4]
    return np.interp(x, xs, hs)


# --------------------------------------------------------------------------
# NHWAVE

def set_keys(text, values):
    """Replace 'KEY = value' lines of an NHWAVE input file."""
    for key, val in values.items():
        pat = re.compile(rf"^{re.escape(key)}\s*=.*$", re.M)
        if not pat.search(text):
            raise KeyError(key)
        text = pat.sub(f"{key} = {val}", text, count=1)
    return text


def nhwave_write(case, template, keys, depth, stations):
    """depth: array (Nglob, Mglob) of still-water depths at cell centres
    (negative on land); stations: list of (x, y)."""
    os.makedirs(os.path.join(case, "results"), exist_ok=True)
    keys = dict(keys)
    keys.setdefault("NSTAT", len(stations))
    keys.setdefault("DEPTH_TYPE", "CELL_CENTER")
    txt = set_keys(open(template).read(), keys)
    open(os.path.join(case, "input.txt"), "w").write(txt)
    depth = np.atleast_2d(depth)
    np.savetxt(os.path.join(case, "depth.txt"), depth, fmt="%.6f")
    with open(os.path.join(case, "stat.txt"), "w") as f:
        for x, y in stations:
            f.write(f"{x:.6g} {y:.6g} 0.0\n")


def nhwave_common(Mglob, Nglob, Kglob, dx, dy, total, px=2):
    return {
        "Mglob": Mglob, "Nglob": Nglob, "Kglob": Kglob, "PX": px, "PY": 1,
        "DX": dx, "DY": dy, "IVGRD": 1, "TOTAL_TIME": total,
        "PLOT_INTV": 1000.0, "SCREEN_INTV": 1.0, "PLOT_INTV_STAT": 0.01,
        "HIGH_ORDER": "SECOND", "TIME_ORDER": "SECOND", "CONVECTION": "HLPA",
        "NON_HYDRO": "T", "CFL": 0.5, "HOTSTART": "F",
        "OUT_H": "F", "OUT_E": "F",
    }


def nhwave_bar(case, template, H, T, nLayers=10, dx=0.02, total=50.0, px=2):
    """Submerged bar: linear waves from the left boundary (x = 0), 8 m
    sponge at the right end, as the NHFOAM case."""
    M = int(round(54.0 / dx))
    x = (np.arange(M) + 0.5) * dx
    keys = nhwave_common(M, 1, nLayers, dx, dx, total, px)
    keys.update({"WAVEMAKER": "LEF_LIN", "AMP": H, "PER": T, "DEP": 0.4,
                 "BC_X0": 3, "BC_Xn": 1, "SPONGE_ON": "T",
                 "Sponge_West_Width": 0.0, "Sponge_East_Width": 8.0,
                 "VISCOUS_FLOW": "F"})
    nhwave_write(case, template, keys, bar_depth(x)[None, :],
                 [(g, 0.5 * dx) for g in BAR_GAUGES])


def hs_depth(x, xToe=0.0, h0=0.36, slope=1 / 34.26):
    """Hansen & Svendsen (1979): flat 0.36 m, then 1:34.26 up the beach."""
    return np.where(x < xToe, h0, h0 - slope * (x - xToe))


def nhwave_hs(case, template, H=0.0686, T=1.667, nLayers=10, dx=0.025, x0=-4.0,
              x1=13.5, total=66.68, px=2, viscous=True, gauge_dx=0.1):
    """Hansen & Svendsen spilling breakers with the full beach (NHWAVE has
    wetting and drying): linear waves from the left boundary at x = x0."""
    M = int(round((x1 - x0) / dx))
    x = x0 + (np.arange(M) + 0.5) * dx
    keys = nhwave_common(M, 1, nLayers, dx, dx, total, px)
    keys.update({"WAVEMAKER": "LEF_LIN", "AMP": H, "PER": T, "DEP": 0.36,
                 "BC_X0": 3, "BC_Xn": 1, "SPONGE_ON": "F",
                 "Sponge_West_Width": 0.0, "Sponge_East_Width": 0.0,
                 "VISCOUS_FLOW": "T" if viscous else "F", "MinDep": 0.005})
    gx = [round(-1.0 + gauge_dx * i, 3) for i in range(int(round((11.0 + 1.0) / gauge_dx)) + 1)]
    nhwave_write(case, template, keys, hs_depth(x)[None, :],
                 [(g - x0, 0.5 * dx) for g in gx])
    with open(os.path.join(case, "gauges_x.txt"), "w") as f:
        f.write("\n".join(f"{g:g}" for g in gx) + "\n")


# --------------------------------------------------------------------------
# REEF3D::NHFLOW

def reef3d_write(case, control, ctrl):
    os.makedirs(case, exist_ok=True)
    open(os.path.join(case, "control.txt"), "w").write(control.strip() + "\n")
    open(os.path.join(case, "ctrl.txt"), "w").write(ctrl.strip() + "\n")


def reef3d_bar(case, H, T, nLayers=10, dx=0.02, total=50.0, nproc=2):
    """Submerged bar: 8 m relaxation zones for generation and absorption,
    linear waves, as the NHFOAM case. The bar is built from solids."""
    nx = int(round(54.0 / dx))
    control = f"""
C 11 6
C 12 3
C 13 3
C 14 7
C 15 21
C 16 3

B 2 {nx} 1 {nLayers}
B 10 0.0 54.0 0.0 {dx:g} 0.0 1.0

S 61 26.0 32.0 0.0 {dx:g} 0.0 0.3
S 10 32.0 34.0 0.0 {dx:g} 0.0 0.3
S 61 34.0 37.0 0.0 {dx:g} 0.3 0.0

M 10 {nproc}
"""
    gauges = "\n".join(f"P 51 {g:g} {0.5 * dx:g}" for g in BAR_GAUGES)
    ctrl = f"""
A 10 5

B 90 1
B 92 2
B 93 {H:g} {T:g}
B 96 8.0 8.0
B 98 2
B 99 2

F 60 0.4

N 41 {total:g}
N 47 1.0
M 10 {nproc}

P 10 0

W 22 -9.81

{gauges}
"""
    reef3d_write(case, control, ctrl)


def reef3d_hs(case, H=0.0686, T=1.667, nLayers=10, dx=0.025, total=66.68,
              zone=8.0, flat=6.0, nproc=2, breaking=True):
    """Hansen & Svendsen spilling breakers with the full beach (NHFLOW has
    wetting and drying). x = 0 at the toe in the NHFOAM/lab coordinates;
    the REEF3D domain starts at x0 = -(zone + flat)."""
    x0 = -(zone + flat)
    x1 = 13.5
    L = x1 - x0
    nx = int(round(L / dx))
    xToe = -x0                         # toe in the REEF3D coordinates
    # beach as a wedge solid from the toe to the end of the domain
    zEnd = (x1 - 0.0) / 34.26          # bed rise at the end of the domain
    control = f"""
C 11 6
C 12 3
C 13 3
C 14 21
C 15 21
C 16 3

B 2 {nx} 1 {nLayers}
B 10 0.0 {L:g} 0.0 {dx:g} 0.0 1.0

S 61 {xToe:g} {L:g} 0.0 {dx:g} 0.0 {zEnd:.6g}

M 10 {nproc}
"""
    gx = [round(-1.0 + 0.1 * i, 3) for i in range(int(round((11.0 + 1.0) / 0.1)) + 1)]
    gauges = "\n".join(f"P 51 {g - x0:g} {0.5 * dx:g}" for g in gx)
    brk = "A 550 1\nA 551 1\nA 554 0.6\n" if breaking else ""
    ctrl = f"""
A 10 5
{brk}
B 90 1
B 92 2
B 93 {H:g} {T:g}
B 96 {zone:g} 0.0
B 98 2
B 99 0

F 60 0.36

N 41 {total:g}
N 47 1.0
M 10 {nproc}

P 10 0

W 22 -9.81

{gauges}
"""
    reef3d_write(case, control, ctrl)
    with open(os.path.join(case, "gauges_x.txt"), "w") as f:
        f.write("\n".join(f"{g:g}" for g in gx) + "\n")


# --------------------------------------------------------------------------
# Readers: gauge time series as (t, {x: eta})

def read_nhfoam(case, xs=None):
    """NHFOAM gauges; keys are the gauge x positions (or names if xs is
    given as a list of names)."""
    fn = os.path.join(case, "postProcessing", "gauges", "0", "eta.dat")
    d = np.loadtxt(fn, comments="#")
    hdr = [l.split() for l in open(fn) if l.startswith("# ") and "(" in l]
    out = {}
    for k, h in enumerate(hdr):
        x = float(h[2].strip("("))
        out[x] = d[:, k + 1]
    return d[:, 0], out


def parse_case(arg):
    """'[nhfoam|nhwave|reef3d:]path[:label]' -> (kind, path, label)."""
    parts = arg.split(":")
    kind = "nhfoam"
    if parts[0] in ("nhfoam", "nhwave", "reef3d"):
        kind = parts.pop(0)
    path = parts[0]
    label = parts[1] if len(parts) > 1 else os.path.basename(path.rstrip("/"))
    return kind, path, label


def read_any(kind, path, xs):
    """Gauge series at positions xs (in the NHFOAM coordinates) for any of
    the three models. NHWAVE and REEF3D gauges are in the order of xs."""
    if kind == "nhfoam":
        t, d = read_nhfoam(path)
        return t, {x: d[min(d, key=lambda g: abs(g - x))] for x in xs}
    if kind == "nhwave":
        return read_nhwave(path, xs)
    return read_reef3d(path, xs)

def read_nhwave(case, xs):
    out = {}
    t = None
    for i, x in enumerate(xs):
        d = np.loadtxt(os.path.join(case, "results", f"probe_{i + 1:04d}"))
        t = d[:, 0]
        out[x] = d[:, 1]
    return t, out


def read_reef3d(case, xs):
    fn = os.path.join(case, "REEF3D_NHFLOW_WSF", "REEF3D-NHFLOW-WSF-HG.dat")
    lines = open(fn).read().splitlines()
    i0 = [k for k, l in enumerate(lines) if l.startswith("time")][0]
    d = np.array([list(map(float, l.split())) for l in lines[i0 + 1:] if l.strip()])
    return d[:, 0], {x: d[:, k + 1] for k, x in enumerate(xs)}


def nhwave_standing(case, template, kh, nLayers, h=1.0, Nx=64, a_over_h=1e-3,
                    nPeriods=5, cfl=0.3):
    """Standing wave in a closed basin one wavelength long, as the NHFOAM
    dispersion study: eta = a cos(k x), U = 0 (eta0.txt, uvw0.txt)."""
    g = 9.81
    k = kh / h
    L = 2 * math.pi / k
    T = 2 * math.pi / math.sqrt(g * k * math.tanh(kh))
    dx = L / Nx
    x = (np.arange(Nx) + 0.5) * dx
    keys = nhwave_common(Nx, 1, nLayers, dx, dx, nPeriods * T, px=1)
    keys.update({"WAVEMAKER": "LEF_LIN", "AMP": 0.0, "PER": T, "DEP": h,
                 "BC_X0": 1, "BC_Xn": 1, "SPONGE_ON": "F",
                 "Sponge_West_Width": 0.0, "Sponge_East_Width": 0.0,
                 "VISCOUS_FLOW": "F", "INITIAL_EUVW": "T", "CFL": cfl,
                 "DT_MAX": T / 200, "DT_INI": T / 200, "PLOT_INTV_STAT": T / 200})
    nhwave_write(case, template, keys, np.full((1, Nx), h), [(0.5 * dx, 0.5 * dx)])
    np.savetxt(os.path.join(case, "eta0.txt"), (a_over_h * h * np.cos(k * x))[None, :], fmt="%.10e")
    np.savetxt(os.path.join(case, "uvw0.txt"), np.zeros((3 * nLayers, Nx)), fmt="%.1f")
    return dict(T=T, L=L)


def standing_period(t, e):
    """Mean period from the upward zero crossings (linear interpolation)."""
    s = np.sign(e - e.mean())
    i = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    tc = t[i] - e[i] * (t[i + 1] - t[i]) / (e[i + 1] - e[i])
    return np.mean(np.diff(tc)) if len(tc) > 1 else np.nan
