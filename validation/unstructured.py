#!/usr/bin/env python3
"""Cases on unstructured (triangulated) meshes: gmsh triangles extruded into
sigma columns of prisms.  gmsh and the OpenFOAM environment (gmshToFoam) are
needed to write them.

    runup_strip(case, ...)   solitary wave on a 1:20 beach (Synolakis 1987,
                             H/d = 0.28, breaking, wetting and drying) on a
                             triangulated strip of width W, the set-up of
                             wetdry.runup_nhfoam (the structured run of the
                             report) otherwise
    cylinder(case, ...)      linear diffraction of regular waves by a vertical
                             circular cylinder (MacCamy and Fuchs 1954) on a
                             mesh refined around the cylinder, waves from the
                             boundary wavemaker (half domain, symmetry wall
                             on the axis y = 0)

usage: unstructured.py runup <case> [W lc] | cylinder <case>
"""
import os, re, subprocess, sys
import numpy as np
import make_case as mc
import wetdry
import ma2012
from longshore import cell_centres

G = 9.81


def _gmsh(case, geo):
    with open(os.path.join(case, "mesh.geo"), "w") as f:
        f.write(geo)
    subprocess.run(["gmsh", "-3", "mesh.geo", "-format", "msh2", "-o", "mesh.msh"],
                   cwd=case, check=True, stdout=open(os.path.join(case, "log.gmsh"), "w"),
                   stderr=subprocess.STDOUT)
    subprocess.run(["gmshToFoam", "mesh.msh"], cwd=case, check=True,
                   stdout=open(os.path.join(case, "log.gmshToFoam"), "w"),
                   stderr=subprocess.STDOUT)
    # all patches but the free surface are walls
    bf = os.path.join(case, "constant/polyMesh/boundary")
    txt = open(bf).read()
    txt = re.sub(r"(\n\s*(\w+)\s*\{\s*type\s+)patch;",
                 lambda m: m.group(1) + ("patch;" if m.group(2) == "surface" else "wall;"), txt)
    txt = re.sub(r"\s*physicalType\s+\w+;", "", txt)
    open(bf, "w").write(txt)


def _points(case):
    pf = os.path.join(case, "constant/polyMesh/points")
    txt = open(pf).read()
    i = txt.index("FoamFile")
    i = txt.index("}", i) + 1
    m = re.search(r"\n(\d+)\s*\n\(", txt[i:])
    n = int(m.group(1))
    a, b = i + m.end(), txt.rfind(")")
    P = np.array(txt[a:b].replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)
    return pf, txt[:a], P, txt[b:]


def _write_points(pf, head, P, tail):
    with open(pf, "w") as f:
        f.write(head + "\n")
        f.writelines(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n" for p in P)
        f.write(tail)


# --------------------------------------------------------------------------
# Solitary wave run-up on a triangulated strip

def runup_strip(case, W=0.1, lc=0.02, nz=3, H=0.0588, dt=0.002, tstar_end=60.0,
                film=0.005, Cd=0.0, conv="Gauss limitedLinearV 1", lcShore=None):
    """As wetdry.runup_nhfoam (d = 0.21 m, toe at 4 m, slope 1:20 to 11 m),
    on gmsh triangles of size lc (lcShore from x = 7 m on, if given),
    extruded into nz layers; slip walls at y = 0 and y = W."""
    from solitary_exact import SolitaryWave
    d, toe, slope, L = wetdry.RU_D, wetdry.RU_TOE, wetdry.RU_SLOPE, wetdry.RU_LEN
    os.makedirs(case, exist_ok=True)
    xc0, xshore = wetdry.runup_geometry(H, d)
    tscale = np.sqrt(d / G)
    field = ""
    if lcShore:
        field = f"""
Field[1] = MathEval;
Field[1].F = "{lc} + ({lcShore} - {lc})*(0.5 + 0.5*Tanh((x - 7.0)/0.3))";
Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
"""
    mc.system_files(case, endTime=tstar_end * tscale, deltaT=dt, writeInterval=1,
                    nCorr=3, nNonOrth=1, pTol=1e-10, conv=conv)
    cd = open(os.path.join(case, "system/controlDict")).read()
    cd = re.sub(r"writeControl\s+\w+;", "writeControl    runTime;", cd)
    cd = re.sub(r"writeInterval\s+[^;]+;", f"writeInterval   {5 * tscale:.8g};", cd)
    open(os.path.join(case, "system/controlDict"), "w").write(cd)
    _gmsh(case, f"""// solitary wave run-up: triangulated strip, {nz} layers
lc = {lc};
Point(1) = {{0, 0, {-d}, lc}};
Point(2) = {{{L}, 0, {-d}, lc}};
Point(3) = {{{L}, {W}, {-d}, lc}};
Point(4) = {{0, {W}, {-d}, lc}};
Line(1) = {{1, 2}};
Line(2) = {{2, 3}};
Line(3) = {{3, 4}};
Line(4) = {{4, 1}};
Curve Loop(1) = {{1, 2, 3, 4}};
Plane Surface(1) = {{1}};
Mesh.Algorithm = 5;
{field}
out[] = Extrude {{0, 0, {d}}} {{ Surface{{1}}; Layers{{{nz}}}; Recombine; }};
Physical Surface("surface") = {{out[0]}};
Physical Surface("bottom") = {{1}};
Physical Surface("sides") = {{out[2], out[4]}};
Physical Surface("beach") = {{out[3]}};
Physical Surface("offshore") = {{out[5]}};
Physical Volume("water") = {{out[1]}};
""")
    # points onto the beach: bed z_b(x), top max(0, z_b + film)
    pf, head, P, tail = _points(case)
    zb = -wetdry.runup_depth(P[:, 0], d)
    zt = np.maximum(0.0, zb + film)
    P[:, 2] = zb + (P[:, 2] + d) / d * (zt - zb)
    _write_points(pf, head, P, tail)

    brk = """
breaking
{
    model       HFA;
    alpha       0.6;
    beta        0.3;
}
"""
    mc.constant_files(case, nh_extra=f"""
wettingDrying
{{
    minDepth    {film};
}}

bottomFriction
{{
    Cd          {Cd};
}}
{brk}""")
    mc.common_fields(case, walls=("offshore", "beach", "sides"), empty=False)
    mc.allrun(case, """
runApplication NHFOAM
""")
    # initial fields: exact solitary wave at the cell centres, on the heights
    # the sigma mesh gets when NHFOAM places it on the new free surface
    C = cell_centres(case, "0")
    x = C[:, 0]
    h = wetdry.runup_depth(x, d)
    zbc = -h
    ztc = np.maximum(0.0, zbc + film)
    k = np.clip(np.round((C[:, 2] - zbc) / (ztc - zbc) * nz - 0.5), 0, nz - 1)
    sw = SolitaryWave(H, d=d)
    eta = sw.eta(x - xc0)
    wet = eta + h > film
    D = np.where(wet, h + eta, film)
    zc = -h + (k + 0.5) / nz * D
    u, w = sw.velocity(x - xc0, zc)
    u = np.where(wet, u, 0.0); w = np.where(wet, w, 0.0)
    p_rgh = np.where(wet, 0.5 * (sw.c ** 2 - (u - sw.c) ** 2 - w ** 2), G * (film - h))
    etaw = np.where(wet, eta, -h + film)
    ma2012._set_internal(os.path.join(case, "0", "eta"), ma2012._foam_list(etaw))
    ma2012._set_internal(os.path.join(case, "0", "U"),
                         ma2012._foam_list(list(zip(u, w)), vector=True))
    ma2012._set_internal(os.path.join(case, "0", "p_rgh"), ma2012._foam_list(p_rgh))
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"W {W}\nlc {lc}\nlcShore {lcShore}\nnz {nz}\nH {H}\nnCells {len(x)}\n")
    return len(x)


# --------------------------------------------------------------------------
# Diffraction by a vertical cylinder (MacCamy and Fuchs 1954)

CYL = dict(h=0.3, T=1.0, H=0.01, a=0.2, x0=-3.0, x1=5.0, W=4.0, sponge=2.6)


def cylinder_gauges(p=CYL):
    """(name, x, y): around the cylinder at r = a + 2 cm (theta from the
    down-wave axis), along the axis in front and behind, and across."""
    a = p["a"]
    g = []
    r = a + 0.02
    for i, th in enumerate(np.arange(0, 181, 15)):
        t = np.radians(th)
        g.append((f"c{th:03d}", r * np.cos(t), r * np.sin(t)))
    for x in np.arange(-2.0, -a - 0.04, 0.1):
        g.append((f"f{len(g):03d}", x, 0.0))
    for x in np.arange(a + 0.1, 2.31, 0.1):
        g.append((f"b{len(g):03d}", x, 0.0))
    for y in np.arange(a + 0.1, 2.01, 0.1):
        g.append((f"s{len(g):03d}", 0.0, y))
    return [(n, round(x, 4), round(y, 4)) for n, x, y in g]


def cylinder(case, lcFar=0.06, lcNear=0.025, nz=4, dt=0.01, endTime=16.0, p=CYL):
    """Half domain x0..x1, 0..W, cylinder of radius a at the origin; boundary
    wavemaker at x = x0 (active absorption), sponge over the last 'sponge'
    metres; slip walls elsewhere (the axis y = 0 is the symmetry plane)."""
    os.makedirs(case, exist_ok=True)
    h, a, x0, x1, W = p["h"], p["a"], p["x0"], p["x1"], p["W"]
    nsteps = int(round(endTime / dt))
    mc.system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps, nCorr=3,
                    nNonOrth=1, pTol=1e-10, libs=True)
    _gmsh(case, f"""// diffraction by a vertical cylinder: half domain, {nz} layers
SetFactory("OpenCASCADE");
Rectangle(1) = {{{x0}, 0, {-h}, {x1 - x0}, {W}}};
Disk(2) = {{0, 0, {-h}, {a}, {a}}};
BooleanDifference(3) = {{ Surface{{1}}; Delete; }}{{ Surface{{2}}; Delete; }};
Field[1] = MathEval;
Field[1].F = "{lcNear} + ({lcFar} - {lcNear})*Min(1, Max(0, (Sqrt(x*x + y*y) - {a})/0.6))";
Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 5;
out[] = Extrude {{0, 0, {h}}} {{ Surface{{3}}; Layers{{{nz}}}; Recombine; }};
eps = 1e-4;
top[] = Surface In BoundingBox{{{x0} - eps, -eps, -eps, {x1} + eps, {W} + eps, eps}};
bot[] = Surface In BoundingBox{{{x0} - eps, -eps, {-h} - eps, {x1} + eps, {W} + eps, {-h} + eps}};
inl[] = Surface In BoundingBox{{{x0} - eps, -eps, {-h} - eps, {x0} + eps, {W} + eps, eps}};
outl[] = Surface In BoundingBox{{{x1} - eps, -eps, {-h} - eps, {x1} + eps, {W} + eps, eps}};
axs[] = Surface In BoundingBox{{{x0} - eps, -eps, {-h} - eps, {x1} + eps, eps, eps}};
sid[] = Surface In BoundingBox{{{x0} - eps, {W} - eps, {-h} - eps, {x1} + eps, {W} + eps, eps}};
cyl[] = Surface In BoundingBox{{{-a} - eps, -eps, {-h} - eps, {a} + eps, {a} + eps, eps}};
cyl[] -= {{top[], bot[], axs[]}};
Physical Surface("surface") = {{top[]}};
Physical Surface("bottom") = {{bot[]}};
Physical Surface("inlet") = {{inl[]}};
Physical Surface("outlet") = {{outl[]}};
Physical Surface("axis") = {{axs[]}};
Physical Surface("side") = {{sid[]}};
Physical Surface("cylinder") = {{cyl[]}};
Physical Volume("water") = {{out[1]}};
""")
    gtxt = "\n".join(f"    {n}    ({x:g} {y:g} 0);" for n, x, y in cylinder_gauges(p))
    mc.constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({x1 - p['sponge']:g} 0 0);
        end         ({x1:g} 0 0);
        sigmaMax    20;
        exponent    2;
    }}
}}
""")
    mc.common_fields(case, walls=("inlet", "outlet", "axis", "side", "cylinder"), empty=False)
    mc.set_patch(case, "U", "inlet", mc.wave_bc("inlet", f"""theory      linear;
        waveHeight  {p['H']:g};
        wavePeriod  {p['T']:g};
        depth       {h:g};
        direction   (1 0 0);
        origin      (0 0 0);
        rampTime    {2 * p['T']:g};"""))
    mc.allrun(case, """
runApplication NHFOAM
""")
    nc = len(cell_centres(case, "0"))
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write("".join(f"{k} {v}\n" for k, v in p.items()))
        f.write(f"lcFar {lcFar}\nlcNear {lcNear}\nnz {nz}\nnCells {nc}\n")
    return nc


def maccamy_fuchs(x, y, k, a, nmax=60):
    """Complex amplitude of eta/a_I around a cylinder of radius a at the
    origin, for the incident wave Re{a_I exp(i(k x - omega t))}."""
    from scipy.special import jv, jvp, hankel1, h1vp
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    out = np.zeros(np.shape(r), complex)
    for m in range(nmax):
        em = 1.0 if m == 0 else 2.0
        out += em * 1j ** m * (jv(m, k * r) - jvp(m, k * a) / h1vp(m, k * a)
                               * hankel1(m, k * r)) * np.cos(m * th)
    return out


if __name__ == "__main__":
    what, case = sys.argv[1], sys.argv[2]
    if what == "runup":
        kw = {}
        if len(sys.argv) > 3:
            kw["W"] = float(sys.argv[3])
        if len(sys.argv) > 4:
            kw["lc"] = float(sys.argv[4])
        if len(sys.argv) > 5:
            kw["lcShore"] = float(sys.argv[5])
        print("cells:", runup_strip(case, **kw))
    elif what == "cylinder":
        print("cells:", cylinder(case))
