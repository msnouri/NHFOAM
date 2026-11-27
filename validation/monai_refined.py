#!/usr/bin/env python3
"""Monai valley (monai_valley.py) on unstructured meshes refined around the
valley: gmsh quadrilaterals (or triangles) extruded into sigma columns.

The run-up at the head of the valley is low on uniform grids (6.6 cm at
2.8 cm, 7.4 cm at 1.4 cm against 8.0-10.0 cm measured): the valley narrows
to 4-6 cm there.  A uniform grid pays for that resolution everywhere (and
NHWAVE has no other choice).  Here the cell size follows the depth,

    lc = lcDeep sqrt(h/d)   between lcShore and lcDeep in the water,
    lc = lcShore            on the land that the wave can reach (z < zmax),
    lc = lcDeep             on the land above,
    lc = lcValley           in the box of the valley and of the three lab
                            run-up points (BOX),
    lc = lcShelf            (option) on the shelf in front of the valley,
                            with the gauges (SHELF),

with its gradient limited to `grad` (size change per distance).  The bed is
the bilinear interpolation of the benchmark nodes (1.4 cm), so cells
smaller than that add numerical resolution, not bathymetric detail.

Everything else as monai_valley.nhfoam (incident wave at x = 0 as a time
series with active absorption, 4 layers, log-law friction, 2 mm film,
hybrid breaking model).  The time step follows the Courant number
(adjustTimeStep, maxCo) below maxDeltaT, so the small cells cost time
steps only while water moves through them.

The run-up at the head goes 6.6 (uniform 2.8 cm), 7.0-7.4 (1.4 cm),
8.1 (0.7 cm in the valley, 22 455 columns, 22 min on 2 cores against 6.2 h
for the uniform 1.4 cm grid) and 8.5 cm (0.35 cm at the head); README 6.12.

usage (python):
    import monai_refined as mr
    mr.nhfoam("nr_14", lcValley=0.014)     # then ./Allrun (2 processors)
    mr.nhfoam("nr_07", lcValley=0.007)
    mr.nhfoam("nr_s14_07", lcShelf=0.014, lcValley=0.007)   # shelf at 1.4 cm
    mr.nhfoam("nr_035", lcShelf=0.007, shelf=mr.BOX,        # nested boxes
              lcValley=0.0035, box=mr.HEAD)
    mr.nhfoam("nr_07_tri", lcValley=0.007, quads=False)     # triangles
    mr.nhfoam("mv_L1_n012", lcShelf=0.014, lcValley=0.007,  # depth-averaged,
              nz=1, manning=0.012)                          # Manning (README 6.13)
"""
import os, re
import numpy as np
import make_case as mc
import monai_valley as mv
import unstructured as us
from longshore import cell_centres

G = 9.81
BOX = (4.80, 5.25, 1.70, 2.40)       # valley and lab run-up points: x0, x1, y0, y1
SHELF = (4.00, 5.30, 1.20, 2.70)     # the shelf in front of the valley, with the gauges
HEAD = (4.90, 5.22, 1.80, 2.02)      # the head of the valley


# --------------------------------------------------------------------------
# Mesh

def size_field(lcDeep=0.056, lcShore=0.028, lcValley=0.014, box=BOX, grad=0.2,
               zmax=0.08, lcShelf=None, shelf=SHELF):
    """Cell size at the benchmark nodes, lc[j, i] (see the header); with
    lcShelf also that size in the box `shelf` in front of the valley."""
    Z = mv.bathymetry()
    ny, nx = Z.shape
    x = np.arange(nx) * mv.DX0
    y = np.arange(ny) * mv.DX0
    X, Y = np.meshgrid(x, y)
    lc = np.clip(lcDeep * np.sqrt(np.maximum(Z, 0.0) / mv.D), lcShore, lcDeep)
    lc[(Z < 0) & (-Z <= zmax)] = lcShore
    lc[-Z > zmax] = lcDeep
    inbox = (X >= box[0]) & (X <= box[1]) & (Y >= box[2]) & (Y <= box[3])
    lc[inbox] = np.minimum(lc[inbox], lcValley)
    if lcShelf:
        inshelf = ((X >= shelf[0]) & (X <= shelf[1]) & (Y >= shelf[2]) & (Y <= shelf[3])
                   & (-Z <= zmax))
        lc[inshelf] = np.minimum(lc[inshelf], lcShelf)
    # gradient limit: lc(a) <= lc(b) + grad |a - b| (sweeps over the 8 neighbours)
    for _ in range(5000):
        old = lc.copy()
        for sj, si in ((0, 1), (1, 0), (1, 1), (1, -1)):
            d = grad * mv.DX0 * np.hypot(sj, si)
            for sgn in (1, -1):
                a, b = sgn * sj, sgn * si
                dst = lc[max(a, 0):ny + min(a, 0), max(b, 0):nx + min(b, 0)]
                src = lc[max(-a, 0):ny + min(-a, 0), max(-b, 0):nx + min(-b, 0)]
                np.minimum(dst, src + d, out=dst)
        if np.abs(lc - old).max() < 1e-9:
            break
    return x, y, lc


def _write_size_field(fn, x, y, lc):
    """gmsh Structured field (text): origin, spacing, counts, values with x
    outermost; two identical levels in z."""
    with open(fn, "w") as f:
        f.write(f"0 0 -1\n{mv.DX0} {mv.DX0} 2\n{len(x)} {len(y)} 2\n")
        for i in range(len(x)):
            f.writelines(f"{v:.6f} {v:.6f}\n" for v in lc[:, i])


def mesh(case, nz=mv.NZ, quads=True, film=mv.FILM, **size):
    """gmsh mesh of the benchmark domain, extruded into nz layers, with the
    points on the bed (dry columns: the film); patches surface, bottom,
    left (x = 0), right, south, north."""
    os.makedirs(case, exist_ok=True)
    x, y, lc = size_field(**size)
    _write_size_field(os.path.join(case, "lc.dat"), x, y, lc)
    alg = ("Mesh.Algorithm = 8;\nMesh.RecombinationAlgorithm = 3;\nRecombine Surface{1};"
           if quads else "Mesh.Algorithm = 6;")
    us._gmsh(case, f"""// Monai valley: cell size from lc.dat, {nz} layers
Point(1) = {{0, 0, {-mv.D}}};
Point(2) = {{{mv.LX}, 0, {-mv.D}}};
Point(3) = {{{mv.LX}, {mv.LY}, {-mv.D}}};
Point(4) = {{0, {mv.LY}, {-mv.D}}};
Line(1) = {{1, 2}};
Line(2) = {{2, 3}};
Line(3) = {{3, 4}};
Line(4) = {{4, 1}};
Curve Loop(1) = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};
Field[1] = Structured;
Field[1].FileName = "lc.dat";
Field[1].TextFormat = 1;
Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
{alg}
out[] = Extrude {{0, 0, {mv.D}}} {{ Surface{{1}}; Layers{{{nz}}}; Recombine; }};
Physical Surface("surface") = {{out[0]}};
Physical Surface("bottom") = {{1}};
Physical Surface("south") = {{out[2]}};
Physical Surface("right") = {{out[3]}};
Physical Surface("north") = {{out[4]}};
Physical Surface("left") = {{out[5]}};
Physical Volume("water") = {{out[1]}};
""")
    pf, head, P, tail = us._points(case)
    zbp = mv.zb(P[:, 0], P[:, 1])
    ztop = np.maximum(0.0, zbp + film)
    P[:, 2] = zbp + (P[:, 2] + mv.D) / mv.D * (ztop - zbp)
    us._write_points(pf, head, P, tail)


# --------------------------------------------------------------------------
# Case

def nhfoam(case, lcDeep=0.056, lcShore=0.028, lcValley=0.014, quads=True,
           nz=mv.NZ, dt=0.008, maxCo=1.0, endTime=mv.T_END, film=mv.FILM,
           ks=mv.KS, breaking=True, nproc=2, conv="Gauss limitedLinearV 1",
           frames=(15.0, 18.0, 0.1), box=BOX, grad=0.2, lcShelf=None, shelf=SHELF,
           manning=None):
    """monai_valley.nhfoam on the refined mesh; dt is the largest time step
    (maxCo = None: fixed time step dt); manning: Manning's n instead of the
    log law with ks."""
    os.makedirs(case, exist_ok=True)
    mc.system_files(case, endTime=endTime, deltaT=dt, writeInterval=1, nCorr=3,
                    nNonOrth=1, pTol=1e-10, conv=conv)
    cd = os.path.join(case, "system/controlDict")
    s = open(cd).read()
    s = re.sub(r"writeControl\s+\w+;", "writeControl    adjustableRunTime;", s)
    s = re.sub(r"writeInterval\s+[^;]+;", f"writeInterval   {endTime:g};", s)
    if maxCo:
        s = re.sub(r"adjustTimeStep\s+\w+;", "adjustTimeStep  yes;", s)
        s = re.sub(r"maxCo\s+[^;]+;", f"maxCo           {maxCo:g};", s)
        s = re.sub(r"maxDeltaT\s+[^;]+;", f"maxDeltaT       {dt:g};", s)
    open(cd, "w").write(s)
    mesh(case, nz=nz, quads=quads, film=film, lcDeep=lcDeep, lcShore=lcShore,
         lcValley=lcValley, box=box, grad=grad, lcShelf=lcShelf, shelf=shelf)

    gtxt = "\n".join(f"    {g}    ({0.5 * lcDeep if x is None else x:.4f} {y:.4f} 0);"
                     for g, (x, y) in mv.GAUGES.items())
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
                     np.round(np.arange(frames[0], frames[1] + 1e-9, frames[2]), 6))
    fric = (f"manning     {manning:g};" if manning
            else f"ks          {ks:g};             // log law, as NHWAVE's Zob")
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film:g};
}}

bottomFriction
{{
    {fric}
}}
{brk}
etaMax          yes;

gauges
{{
{gtxt}
}}

surfaceSnapshots ({snaps});
""")
    mc.common_fields(case, walls=("left", "right", "south", "north"), empty=False)
    t, e = mv.incident()
    with open(os.path.join(case, "constant", "incident.dat"), "w") as f:
        f.write(f"{len(t)}\n(\n")
        f.writelines(f"({a:.6f} {b:.8e})\n" for a, b in zip(t, e))
        f.write(")\n")
    text = f"""theory      timeSeries;
        depth       {mv.D:g};
        direction   (1 0 0);
        origin      (0 0 0);
        eta
        {{
            type    tableFile;
            file    "$FOAM_CASE/constant/incident.dat";
        }}"""
    mc.set_patch(case, "U", "left", mc.wave_bc("left", text, absorb=True))
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
""" if nproc > 1 else "runApplication NHFOAM\n")

    # still water; dry columns below their bed, so that NHFOAM puts the film
    # on them
    import ma2012
    C = cell_centres(case, "0")
    np.save(os.path.join(case, "cellCentres.npy"), C[:, :2])
    zbc = mv.zb(C[:, 0], C[:, 1])
    wet = -zbc > film
    ma2012._set_internal(os.path.join(case, "0", "eta"),
                         ma2012._foam_list(np.where(wet, 0.0, -10.0)))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"),
                         ma2012._foam_list(np.where(wet, 0.0, G * (zbc + film))))
    ncol = len(C) // nz
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"model nhfoam\nmesh {'quads' if quads else 'triangles'}\n"
                f"lcDeep {lcDeep:g}\nlcShore {lcShore:g}\nlcValley {lcValley:g}\n"
                f"lcShelf {lcShelf if lcShelf else lcShore:g}\n"
                f"columns {ncol}\nnz {nz}\ndt {dt:g}\nmaxCo {maxCo}\nfilm {film:g}\n"
                f"ks {ks:g}\nbreaking {breaking}\n")
    return ncol


# --------------------------------------------------------------------------
# Results (columns as scattered points)

def columns(case, name=None):
    """x, y, bottom and free surface of the columns of a snapshot
    (postProcessing/surface/<name>; default: the first one)."""
    root = os.path.join(case, "postProcessing", "surface")
    if name is None:
        name = min(os.listdir(root), key=lambda s: float(s[:-4]))
    a = np.loadtxt(os.path.join(root, name))
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def surface(case, time):
    """(x, y, zb, eta, wet, t) of the snapshot nearest to a time."""
    root = os.path.join(case, "postProcessing", "surface")
    name = min(os.listdir(root), key=lambda s: abs(float(s[:-4]) - time))
    x, y, z, e = columns(case, name)
    film = mv.params(case)["film"]
    return x, y, z, e, e - z > 1.5 * film, float(name[:-4])


def emax(case):
    """(x, y, zb, etaMax, flooded) of the columns: etaMax of the last time
    (never-wet columns keep their bottom)."""
    from scipy.spatial import cKDTree
    import conical_island as ci
    pr = mv.params(case)
    x, y, z, _ = columns(case)
    tt = mv._times(case)[-1]
    em = ci._field_values(os.path.join(case, tt, "etaMax"))
    C = np.load(os.path.join(case, "cellCentres.npy"))
    # the snapshot has the face centres of the columns, cellCentres.npy the
    # means of the cell vertices: each cell goes to the nearest column, and
    # every column must get its nz cells
    _, k = cKDTree(np.c_[x, y]).query(C)
    assert (np.bincount(k, minlength=len(x)) == int(pr["nz"])).all()
    e = np.empty(len(x))
    e[k] = em
    return x, y, z, e, (z < 0.0) | (e > z + 0.5 * pr["film"])


def runup(case, half_width=0.5 * mv.DX0):
    """As monai_valley.runup: the highest flooded ground in the strip
    |y - y_lab| <= half_width (half a benchmark cell) within 0.3 m of the
    lab x, and the most inland x reached there."""
    x, y, z, e, fl = emax(case)
    out = []
    for xl, yl, obs in mv.read_runup():
        s = (np.abs(y - yl) <= half_width) & (np.abs(x - xl) <= 0.3) & fl & (z > 0.0)
        out.append(dict(x=xl, y=yl, obs=obs, R=z[s].max() if s.any() else np.nan,
                        x_inland=x[s].max() if s.any() else np.nan))
    return out


# --------------------------------------------------------------------------
# All runs of the README/report (section 6.12):
#     python3 monai_refined.py <dir>    (then <dir>/runall.sh: 16-28 min each
#                                        on 2 cores, the triangles longer)

RUNS = {
    "nr_14": dict(lcValley=0.014),
    "nr_s14": dict(lcShelf=0.014, lcValley=0.014),
    "nr_07": dict(lcValley=0.007),
    "nr_s14_07": dict(lcShelf=0.014, lcValley=0.007),
    "nr_035": dict(lcShelf=0.007, shelf=BOX, lcValley=0.0035, box=HEAD),
    "nr_07_tri": dict(lcValley=0.007, quads=False),
}

if __name__ == "__main__":
    import sys
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    os.chdir(out)
    for name, kw in RUNS.items():
        print(name, nhfoam(name, **kw), "columns")
    # the film thickness on the uniform 2.8 cm grid
    mv.nhfoam("nf_28_f1", dx=0.028, dt=0.008, film=0.001)
    mv.nhfoam("nf_28_f05", dx=0.028, dt=0.008, film=0.0005)
    with open("runall.sh", "w") as f:
        f.write("#!/bin/sh\n# NHFOAM (Allrun: 2 processors)\n")
        for c in list(RUNS) + ["nf_28_f1", "nf_28_f05"]:
            f.write(f"(cd {c} && ./Allrun)\n")
    os.chmod("runall.sh", 0o755)
