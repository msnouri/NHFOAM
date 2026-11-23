#!/usr/bin/env python3
"""Generate 2DV NHFOAM cases (standing wave, solitary wave, Beji-Battjes bar)."""
import os, sys, math, textwrap

HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | NHFOAM
   \\\\    /   O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
FTR = "\n// ************************************************************************* //\n"

def write(case, rel, cls, body):
    path = os.path.join(case, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    loc, obj = os.path.split(rel)
    with open(path, "w") as f:
        f.write(HDR.format(cls=cls, loc=loc, obj=obj) + textwrap.dedent(body).strip() + "\n" + FTR)

def common_fields(case, walls=("left", "right"), extra_U=None, extra_p=None, empty=True,
                  constraint_types=False):
    wallU = "".join(f"    {w}\n    {{\n        type            slip;\n    }}\n" for w in walls)
    wallP = "".join(f"    {w}\n    {{\n        type            fixedFluxPressure;\n        value           uniform 0;\n    }}\n" for w in walls)
    fbE = "    frontAndBack\n    {\n        type            empty;\n    }\n" if empty else ""
    wallE = "".join(f"    {w}\n    {{\n        type            zeroGradient;\n    }}\n" for w in walls)
    write(case, "0/U", "volVectorField", f"""
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{{
    surface
    {{
        type            zeroGradient;
    }}
    bottom
    {{
        type            slip;
    }}
{wallU}{fbE}}}
""")
    write(case, "0/p_rgh", "volScalarField", f"""
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    // Free surface: coefficients are set by NHFOAM every corrector
    surface
    {{
        type            mixed;
        refValue        uniform 0;
        refGradient     uniform 0;
        valueFraction   uniform 1;
        value           uniform 0;
    }}
    bottom
    {{
        type            fixedFluxPressure;
        value           uniform 0;
    }}
{wallP}{fbE}}}
""")
    write(case, "0/eta", "volScalarField", f"""
dimensions      [0 1 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    surface
    {{
        type            zeroGradient;
    }}
    bottom
    {{
        type            zeroGradient;
    }}
{wallE}{fbE}}}
""")
    if constraint_types:
        # entries for the constraint patches (cyclic, processor, ...)
        for f in os.listdir(os.path.join(case, "0")):
            path = os.path.join(case, "0", f)
            txt = open(path).read()
            txt = txt.replace("boundaryField\n{\n",
                              "boundaryField\n{\n    #includeEtc \"caseDicts/setConstraintTypes\"\n", 1)
            open(path, "w").write(txt)


# Vertical location of the pressure unknowns (pressureLocation in
# NHFOAMProperties): "interface" (layer interfaces, NHWAVE layout; default)
# or "cell"
PRESSURE = "interface"

def constant_files(case, nu=1e-6, nh_extra=""):
    write(case, "constant/g", "uniformDimensionedVectorField", """
dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);
""")
    write(case, "constant/transportProperties", "dictionary", f"""
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {nu:g};
""")
    write(case, "constant/turbulenceProperties", "dictionary", """
simulationType  laminar;
""")
    write(case, "constant/NHFOAMProperties", "dictionary", f"""
// Free-surface (top) and bottom patches of the sigma mesh
surfacePatch    surface;
bottomPatch     bottom;

// Smallest water depth used to place the mesh points [m]
minDepth        1e-4;

stillWaterLevel 0;

pressureLocation {PRESSURE};
{nh_extra}
""")

def les_files(case, Ck=0.094, patches=None):
    """Smagorinsky subgrid model of the OpenFOAM library (LES, filter width
    the cube root of the cell volume, as in NHWAVE): nu_t = (Cs Delta)^2
    sqrt(2 S:S) with Cs^2 = Ck sqrt(Ck/Ce), Ce = 1.048 (Ck = 0.094: Cs =
    0.168).  patches: {name: type} for nut (default: zeroGradient on the
    surface, bottom and walls, empty frontAndBack)."""
    write(case, "constant/turbulenceProperties", "dictionary", f"""
simulationType  LES;

LES
{{
    LESModel        Smagorinsky;
    turbulence      on;
    printCoeffs     on;
    delta           cubeRootVol;

    SmagorinskyCoeffs
    {{
        Ck              {Ck:g};
        Ce              1.048;
    }}

    cubeRootVolCoeffs
    {{
        deltaCoeff      1;
    }}
}}
""")
    if patches is None:
        patches = {"surface": "zeroGradient", "bottom": "zeroGradient",
                   "left": "zeroGradient", "right": "zeroGradient",
                   "frontAndBack": "empty"}
    bf = "".join(f"    {n}\n    {{\n        type            {t};\n    }}\n"
                 for n, t in patches.items())
    write(case, "0/nut", "volScalarField", f"""
dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
{bf}}}
""")


def rans_files(case, model="kOmegaLF", k0=1e-7, omega0=1.0, coeffs="",
               walls=("left", "right"), empty=True, noslip_bottom=True,
               conv_turb="Gauss limitedLinear 1"):
    """RAS turbulence (default: the stabilised k-omega model kOmegaLF of
    libNHFOAMturbulenceModels).  The bottom gets wall functions and, with
    noslip_bottom, a no-slip velocity (then switch the quadratic bottom
    friction off): the wall function supplies the bed shear stress."""
    write(case, "constant/turbulenceProperties", "dictionary", f"""
simulationType  RAS;

RAS
{{
    RASModel        {model};
    turbulence      on;
    printCoeffs     on;
{coeffs}
}}
""")
    fbE = ("    frontAndBack\n    {\n        type            empty;\n    }\n"
           if empty else "")

    def field(name, dims, value, bottom, surface, wall):
        wl = "".join(f"    {w}\n    {{\n        {wall}\n    }}\n" for w in walls)
        write(case, f"0/{name}", "volScalarField", f"""
dimensions      {dims};

internalField   uniform {value:g};

boundaryField
{{
    surface
    {{
        {surface}
    }}
    bottom
    {{
        {bottom}
    }}
{wl}{fbE}}}
""")

    field("k", "[0 2 -2 0 0 0 0]", k0,
          f"type kqRWallFunction; value uniform {k0:g};",
          "type zeroGradient;", "type zeroGradient;")
    field("omega", "[0 0 -1 0 0 0 0]", omega0,
          f"type omegaWallFunction; value uniform {omega0:g};",
          "type zeroGradient;", "type zeroGradient;")
    field("nut", "[0 2 -1 0 0 0 0]", 0.0,
          "type nutkWallFunction; value uniform 0;",
          "type calculated; value uniform 0;",
          "type calculated; value uniform 0;")
    if noslip_bottom:
        u = open(os.path.join(case, "0/U")).read()
        i = u.index("    bottom")
        j = u.index("}", i)
        u = (u[:i] + "    bottom\n    {\n        type            noSlip;\n    "
             + u[j:])
        open(os.path.join(case, "0/U"), "w").write(u)
    sch = open(os.path.join(case, "system/fvSchemes")).read()
    sch = sch.replace(
        "    div((nuEff*dev2(T(grad(U))))) Gauss linear;",
        "    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n"
        f"    div(phi,k)      {conv_turb};\n"
        f"    div(phi,omega)  {conv_turb};", 1)
    sch = sch.rstrip()
    k = sch.rfind("\n//")
    sch = (sch[:k] + "\n\nwallDist\n{\n    method          meshWave;\n}\n"
           + sch[k:] + "\n")
    open(os.path.join(case, "system/fvSchemes"), "w").write(sch)
    sol = open(os.path.join(case, "system/fvSolution")).read()
    sol = sol.replace(
        '    "(U|UFinal)"',
        '    "(k|omega)(Final)?"\n    {\n        solver          smoothSolver;\n'
        '        smoother        symGaussSeidel;\n        tolerance       1e-12;\n'
        '        relTol          0;\n    }\n\n    "(U|UFinal)"', 1)
    open(os.path.join(case, "system/fvSolution"), "w").write(sol)


# Convection of momentum: a TVD scheme. linearUpwind loses ~12% of H = 4 cm
# over 8 wavelengths, limitedLinearV ~3%; central differencing (Gauss linear)
# loses least but lets grid-scale noise grow in long nonlinear runs (the
# bar cases failed after 50-60 s), and bores need a TVD scheme anyway.
CONV = "Gauss limitedLinearV 1"
CONV_BREAKING = "Gauss limitedLinearV 1"


# Pressure solvers: the interface matrix is solved with PCG and the
# sigmaColumn (line) preconditioner
PSOLVER = {
    "cell": "solver          PCG;\n        preconditioner  DIC;",
    "interface": "solver          PCG;\n        preconditioner  sigmaColumn;",
}

def system_files(case, endTime, deltaT, writeInterval, ddt="backward 1", nCorr=3, nNonOrth=0, pTol=1e-12, extraSolver="", libs=False, conv=None):
    libsLine = 'libs            ("libNHFOAM.so");   // NHFOAM boundary conditions for utilities\n\n' if libs else ""
    write(case, "system/controlDict", "dictionary", f"""
{libsLine}application     NHFOAM;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {endTime:.10g};
deltaT          {deltaT:.10g};

writeControl    timeStep;
writeInterval   {writeInterval};
purgeWrite      0;

writeFormat     ascii;
writePrecision  10;
writeCompression off;
timeFormat      general;
timePrecision   8;

runTimeModifiable true;

adjustTimeStep  no;
maxCo           0.5;
maxDeltaT       1;
""")
    write(case, "system/fvSchemes", "dictionary", f"""
ddtSchemes
{{
    default         {ddt};
}}

gradSchemes
{{
    default         Gauss linear;
}}

divSchemes
{{
    default         none;
    div(phi,U)      {conv or CONV};
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         corrected;
}}
""")
    write(case, "system/fvSolution", "dictionary", f"""
solvers
{{
    p_rgh
    {{
        {PSOLVER[PRESSURE]}
        tolerance       {pTol:g};
        relTol          0.01;
    }}

    p_rghFinal
    {{
        $p_rgh;
        relTol          0;
    }}

    "(U|UFinal)"
    {{
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-12;
        relTol          0;
    }}
{extraSolver}
}}

PIMPLE
{{
    momentumPredictor         yes;
    nOuterCorrectors          1;
    nCorrectors               {nCorr};
    nNonOrthogonalCorrectors  {nNonOrth};
    moveMeshOuterCorrectors   false;
}}

relaxationFactors
{{
    equations
    {{
        ".*"            1;
    }}
}}
""")

def blockmesh_flat(case, L, h, dy, Nx, Nz, zgrade=1.0):
    write(case, "system/blockMeshDict", "dictionary", f"""
scale   1;

vertices
(
    (0 0 {-h})
    ({L} 0 {-h})
    ({L} {dy} {-h})
    (0 {dy} {-h})
    (0 0 0)
    ({L} 0 0)
    ({L} {dy} 0)
    (0 {dy} 0)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({Nx} 1 {Nz}) simpleGrading (1 1 {zgrade})
);

boundary
(
    surface
    {{
        type patch;
        faces ((4 5 6 7));
    }}
    bottom
    {{
        type wall;
        faces ((0 3 2 1));
    }}
    left
    {{
        type wall;
        faces ((0 4 7 3));
    }}
    right
    {{
        type wall;
        faces ((1 2 6 5));
    }}
    frontAndBack
    {{
        type empty;
        faces ((0 1 5 4) (3 7 6 2));
    }}
);
""")

def beach_profile(xs, hs, film, dx=None):
    """Profile with dry ground (hs < 0: bed above the still water level):
    inserts the point where the depth equals the film and returns
    (xs, hs, ztop) with the mesh top at max(0, z_b + film), so that the
    dry columns start with the film (wettingDrying).  With dx given, the
    point is not inserted closer than dx to a vertex (no narrow block; the
    solver puts the surface points in place at the start anyway)."""
    X, H = [xs[0]], [hs[0]]
    for i in range(1, len(xs)):
        h0, h1 = hs[i - 1], hs[i]
        if (h0 - film) * (h1 - film) < 0:
            t = (h0 - film) / (h0 - h1)
            xi = xs[i - 1] + t * (xs[i] - xs[i - 1])
            if dx is None or min(xi - xs[i - 1], xs[i] - xi) >= dx:
                X.append(xi)
                H.append(film)
        X.append(xs[i]); H.append(h1)
    ztop = [max(0.0, -h + film) for h in H]
    return X, H, ztop


def blockmesh_profile(case, xs, hs, dy, dx, Nz, zgrade=1.0, Ny=1, ygrade=1.0,
                      ztop=None, cyclic=False):
    """Piecewise-linear bottom: depth hs[i] at xs[i], one block per segment.
    Ny = 1: 2D (empty frontAndBack, width dy). Ny > 1: 3D channel of width
    dy with walls 'axis' (y = 0) and 'side' (y = dy), or with cyclic: the
    periodic pair 'front' (y = 0) and 'back' (y = dy). ztop: top of the
    mesh per point (default 0; see beach_profile for dry ground)."""
    if ztop is None:
        ztop = [0.0] * len(xs)
    verts = []
    for i, (x, h, zt) in enumerate(zip(xs, hs, ztop)):
        verts += [f"({x:.10g} 0 {-h:.10g})", f"({x:.10g} {dy} {-h:.10g})",
                  f"({x:.10g} 0 {zt:.10g})", f"({x:.10g} {dy} {zt:.10g})"]
    blocks, surf, bot = [], [], []
    for i in range(len(xs) - 1):
        a = 4 * i; b = 4 * (i + 1)
        # hex: bottom (a0, b0, b1, a1), top (a2, b2, b3, a3)
        nx = max(1, int(round((xs[i + 1] - xs[i]) / dx)))
        blocks.append(f"hex ({a} {b} {b+1} {a+1} {a+2} {b+2} {b+3} {a+3}) ({nx} {Ny} {Nz}) simpleGrading (1 {ygrade:g} {zgrade})")
        surf.append(f"({a+2} {b+2} {b+3} {a+3})")
        bot.append(f"({a} {a+1} {b+1} {b})")
    n = len(xs) - 1
    left = "(0 2 3 1)"
    right = f"({4*n} {4*n+1} {4*n+3} {4*n+2})"
    front, back = [], []
    for i in range(n):
        a = 4 * i; b = 4 * (i + 1)
        front.append(f"({a} {b} {b+2} {a+2})")
        back.append(f"({a+1} {a+3} {b+3} {b+1})")
    if Ny == 1:
        sides = f"""    frontAndBack
    {{
        type empty;
        faces ({' '.join(front + back)});
    }}"""
    elif cyclic:
        sides = f"""    front
    {{
        type cyclic;
        neighbourPatch back;
        faces ({' '.join(front)});
    }}
    back
    {{
        type cyclic;
        neighbourPatch front;
        faces ({' '.join(back)});
    }}"""
    else:
        sides = f"""    axis
    {{
        type wall;
        faces ({' '.join(front)});
    }}
    side
    {{
        type wall;
        faces ({' '.join(back)});
    }}"""
    nl = "\n    "
    write(case, "system/blockMeshDict", "dictionary", f"""
scale   1;

vertices
(
    {nl.join(verts)}
);

blocks
(
    {nl.join(blocks)}
);

boundary
(
    surface
    {{
        type patch;
        faces ({' '.join(surf)});
    }}
    bottom
    {{
        type wall;
        faces ({' '.join(bot)});
    }}
    left
    {{
        type wall;
        faces ({left});
    }}
    right
    {{
        type wall;
        faces ({right});
    }}
{sides}
);
""")

def allrun(case, body):
    p = os.path.join(case, "Allrun")
    with open(p, "w") as f:
        f.write("#!/bin/sh\ncd \"${0%/*}\" || exit 1\n. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions\n\n" + textwrap.dedent(body).strip() + "\n")
    os.chmod(p, 0o755)
    p = os.path.join(case, "Allclean")
    with open(p, "w") as f:
        f.write("#!/bin/sh\ncd \"${0%/*}\" || exit 1\n. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions\n\ncleanCase0 2>/dev/null || cleanCase\nrm -rf postProcessing\n")
    os.chmod(p, 0o755)

def g_omega(k, h, g=9.81):
    return math.sqrt(g * k * math.tanh(k * h))

def standing_case(case, kh, nLayers, Nx=64, a_over_h=1e-3, h=1.0, nPeriods=5, stepsPerPeriod=200, nu=1e-6, ddt="backward 1", zgrade=1.0):
    k = kh / h
    L = 2 * math.pi / k
    T = 2 * math.pi / g_omega(k, h)
    dt = T / stepsPerPeriod
    dy = L / Nx
    constant_files(case, nu=nu, nh_extra=f"""
gauges
{{
    G0      (0 {dy/2:.6g} 0);
    Gq      ({L/4:.6g} {dy/2:.6g} 0);
}}
""")
    system_files(case, endTime=nPeriods * T, deltaT=dt, writeInterval=stepsPerPeriod * nPeriods, ddt=ddt)
    blockmesh_flat(case, L, h, dy, Nx, nLayers, zgrade)
    common_fields(case)
    write(case, "system/nhSetWavesDict", "dictionary", f"""
type            standing;
amplitude       {a_over_h * h:.10g};
wavelength      {L:.10g};
origin          (0 0 0);
direction       (1 0 0);
""")
    allrun(case, """
    runApplication blockMesh
    runApplication nhSetWaves
    runApplication NHFOAM
    """)
    return dict(k=k, L=L, T=T, dt=dt)

if __name__ == "__main__":
    case = sys.argv[1]
    kh = float(sys.argv[2]); nl = int(sys.argv[3])
    print(standing_case(case, kh, nl))


def wave_bc(patch, wave, absorb=True, celerity=None):
    """U entry of a boundary wavemaker (nhWaveVelocity) with the given entries
    of the target wave (see regularWave.H)."""
    cel = f"        absorptionCelerity {celerity};\n" if celerity else ""
    body = "\n".join("            " + l.strip() for l in wave.strip().splitlines() if l.strip())
    return f"""    {patch}
    {{
        type            nhWaveVelocity;
        activeAbsorption {"yes" if absorb else "no"};
{cel}        wave
        {{
{body}
        }}
        value           uniform (0 0 0);
    }}
"""


def set_patch(case, field, patch, text, old_type="slip"):
    """Replace the entry of a patch (of type old_type) in 0/<field>."""
    p = os.path.join(case, "0", field)
    s = open(p).read()
    old = f"    {patch}\n    {{\n        type            {old_type};\n    }}\n"
    assert old in s, (case, field, patch)
    open(p, "w").write(s.replace(old, text))


def bar_case(case, nLayers=10, dx=0.02, dt=0.005, endTime=50.0, H=0.02, T=2.02,
             waves=True, sigmaMax=20.0, writeEvery=None, breaking=None, conv=None,
             zgrade=1.0, wavemaker="zone", returnCurrent=False, gauges=None):
    """Luth et al. (1994) / Beji & Battjes (1993) submerged bar, case A
    (flume coordinates of the FUNWAVE benchmark: bar from x = 26 to 37 m).
    wavemaker: "zone" (relaxation zone in 0..8 m, as in the report's runs,
    made before the return current existed) or "bc" (nhWaveVelocity at
    x = 0 with active absorption).  gauges: x positions (default: the
    measuring stations)."""
    xs = [0, 26, 32, 34, 37, 54]
    hs = [0.4, 0.4, 0.1, 0.1, 0.4, 0.4]
    dy = 0.02
    if gauges is None:
        gauges = [22, 24, 30.5, 32.5, 33.5, 34.5, 35.7, 37.3, 39.0, 41.0]
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in gauges)
    rc = "yes" if returnCurrent else "no"
    wave = f"""waveHeight  {H:g};
        wavePeriod  {T:g};
        depth       0.4;
        direction   (1 0 0);
        rampTime    {2*T:g};
        returnCurrent {rc};"""
    zones = ""
    if waves:
        inlet = ""
        if wavemaker == "zone":
            inlet = f"""    inlet
    {{
        type        wave;
        start       (8 0 0);
        end         (0 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    2;
        {wave}
    }}
"""
        zones = f"""
forcingZones
{{
{inlet}    outlet
    {{
        type        sponge;
        start       (46 0 0);
        end         (54 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    2;
    }}
}}
"""
    brk = ""
    if breaking:
        brk = f"""
breaking
{{
    model       HFA;
    alpha       {breaking.get('alpha', 0.6):g};
    beta        {breaking.get('beta', 0.3):g};
}}
"""
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}
{zones}{brk}""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt,
                 writeInterval=writeEvery or nsteps, nCorr=3, nNonOrth=1, pTol=1e-10,
                 conv=conv or (CONV_BREAKING if breaking else CONV))
    blockmesh_profile(case, xs, hs, dy, dx, nLayers, zgrade=zgrade)
    common_fields(case)
    if waves and wavemaker == "bc":
        set_patch(case, "U", "left", wave_bc("left", wave))
    allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)


def solitary_case(case, H=0.1, h=1.0, nLayers=10, dx=0.05, dt=0.01, L=60.0, x0=10.0, endTime=10.0):
    """Solitary wave propagating over a flat bottom between closed walls."""
    Nx = int(round(L / dx))
    dy = dx
    gauges = [10, 20, 30, 40]
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in gauges)
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}
""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps // 4, nCorr=3, pTol=1e-10)
    blockmesh_flat(case, L, h, dy, Nx, nLayers)
    common_fields(case)
    write(case, "system/nhSetWavesDict", "dictionary", f"""
type            solitary;
waveHeight      {H:g};
depth           {h:g};
origin          ({x0:g} 0 0);
direction       (1 0 0);
""")
    allrun(case, """
    runApplication blockMesh
    runApplication nhSetWaves
    runApplication NHFOAM
    """)



def uplift_case(case, h=0.05, b=0.305, zeta0=0.0005, tc=0.1, history="exponential",
                L=15.0, dx=0.005, nLayers=5, dt=0.005, endTime=12.0, sponge=3.0):
    """Tsunami generation by a rectangular bottom uplift (Hammack 1973 type),
    half domain with a symmetry wall at x = 0."""
    Nx = int(round(L / dx))
    dy = dx
    gx = [0.0, b, 5 * b, 10 * b, 20 * b, 30 * b]
    gtxt = "\n".join(f"    G{i}    ({x:.6g} {dy/2:g} 0);" for i, x in enumerate(gx))
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

bottomMotion
{{
    type        uplift;
    min         ({-b:g} -1 0);
    max         ({b:g} 1 0);
    amplitude   {zeta0:g};
    history     {history};
    riseTime    {tc:g};
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({L - sponge:g} 0 0);
        end         ({L:g} 0 0);
        sigmaMax    5;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps, nCorr=3, pTol=1e-12, libs=True)
    blockmesh_flat(case, L, h, dy, Nx, nLayers)
    common_fields(case)
    p = os.path.join(case, "0/U")
    s = open(p).read().replace("""    bottom
    {
        type            slip;
    }""", """    bottom
    {
        type            nhMovingSlipVelocity;
        value           uniform (0 0 0);
    }""")
    open(p, "w").write(s)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"h {h}\nb {b}\nzeta0 {zeta0}\ntc {tc}\ndx {dx}\n")
    allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)


def slide_case(case, nLayers=10, dx=0.01, dt=0.005, endTime=4.0, slope=15.0,
               hTop=0.1, hDeep=1.5, L=15.0, xc0=0.8):
    """2D rigid solid slide on a plane slope (Grilli & Watts 2005 / Enet &
    Grilli 2007 shape and law of motion), demonstration case."""
    tanb = math.tan(math.radians(slope))
    xToe = (hDeep - hTop) / tanb
    xs = [0, xToe, L]
    hs = [hTop, hDeep, hDeep]
    dy = dx
    gauges = [1.5, 2.5, 3.5, 5.0, 8.0]
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in gauges)
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}
{wd}
bottomMotion
{{
    type        slide;
    centre      ({xc0:g} 0 0);
    direction   (1 0 0);
    slopeAngle  {slope:g};
    thickness   0.082;
    length      0.395;
    width       0;                  // 2D slide
    terminalVelocity     1.70;
    initialAcceleration  1.21;
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({L - 3:g} 0 0);
        end         ({L:g} 0 0);
        sigmaMax    10;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=max(1, nsteps // 4),
                 nCorr=3, nNonOrth=1, pTol=1e-10, libs=True)
    blockmesh_profile(case, xs, hs, dy, dx, nLayers)
    common_fields(case)
    p = os.path.join(case, "0/U")
    s = open(p).read().replace("""    bottom
    {
        type            slip;
    }""", """    bottom
    {
        type            nhMovingSlipVelocity;
        value           uniform (0 0 0);
    }""")
    open(p, "w").write(s)
    allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)


def flume_case(case, H=0.041, T=1.01, h=0.4, L=30.0, dx=0.02, dt=0.01, nLayers=10,
               zone=8.0, sigmaMax=20.0, exponent=2.0, endTime=30.0, conv=None):
    """Flat flume: wave-making zone at the left, sponge at the right; for
    checking the generated wave height."""
    Nx = int(round(L / dx))
    dy = dx
    gauges = [zone + 1 + 2 * i for i in range(7)]
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in gauges)
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}

forcingZones
{{
    inlet
    {{
        type        wave;
        start       ({zone:g} 0 0);
        end         (0 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    {exponent:g};
        waveHeight  {H:g};
        wavePeriod  {T:g};
        depth       {h:g};
        direction   (1 0 0);
        rampTime    {2*T:g};
    }}
    outlet
    {{
        type        sponge;
        start       ({L - 8:g} 0 0);
        end         ({L:g} 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    {exponent:g};
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps, nCorr=3, pTol=1e-10, conv=conv)
    blockmesh_flat(case, L, h, dy, Nx, nLayers)
    common_fields(case)


def enet_grilli_case(case, d=61, dx=0.025, nLayers=5, dt=0.005, endTime=4.0,
                     xWall=0.2, L=10.0, Ly=1.8, ygrade=1.0, sponge=3.0, conv=None,
                     zgrade=1.0, film=0.005):
    """NTHMP landslide benchmark 2: 3D rigid slide of Enet & Grilli (2007) on
    a 15 deg slope, half tank with a symmetry wall on the axis (y = 0).
    x is measured from the still-water shoreline; the domain starts at
    x = xWall and the slope meets the 1.5 m deep tank bottom at x = 5.6 m.
    xWall > 0: a wall in the water (no wetting and drying, the runs of the
    report's sections 4 and 5); xWall < 0: the beach up to the bed height
    -xWall tan(15 deg) above the still water, with wetting and drying
    (minDepth film)."""
    # Initial position of the slide centre and fitted law of motion
    # (Enet & Grilli 2007; values as in the NHWAVE benchmark input)
    slides = {61: dict(x0=0.551, ut=1.70, a0=1.12)}
    s = slides[d]
    tanb = math.tan(math.radians(15.0))
    h0 = 1.5
    xs = [xWall, h0 / tanb, L]
    hs = [xWall * tanb, h0, h0]
    ztop = None
    wd = ""
    if xWall < 0:
        xs, hs, ztop = beach_profile(xs, hs, film, dx=dx)
        wd = f"""
wettingDrying
{{
    minDepth    {film:g};
}}
"""
    Ny = int(round(Ly / dx))
    gauges = {"g1": (s["x0"], 0.0), "g2": (1.469, 0.35), "g3": (1.929, 0.0),
              "g4": (1.919, 0.5)}
    gtxt = "\n".join(f"    {n}    ({x:g} {y:g} 0);" for n, (x, y) in gauges.items())
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}
{wd}
bottomMotion
{{
    type        slide;
    centre      ({s['x0']:g} 0 0);
    direction   (1 0 0);
    slopeAngle  15;
    thickness   0.082;
    length      0.395;
    width       0.680;
    epsilon     0.717;
    terminalVelocity     {s['ut']:g};
    initialAcceleration  {s['a0']:g};
}}

forcingZones
{{
    outlet
    {{
        type        sponge;
        start       ({L - sponge:g} 0 0);
        end         ({L:g} 0 0);
        sigmaMax    5;
        exponent    2;
    }}
}}
""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=max(1, nsteps // 8),
                 nCorr=3, nNonOrth=1, pTol=1e-10, libs=True, conv=conv)
    blockmesh_profile(case, xs, hs, Ly, dx, nLayers, zgrade=zgrade, Ny=Ny, ygrade=ygrade,
                      ztop=ztop)
    common_fields(case, walls=("left", "right", "axis", "side"), empty=False)
    p = os.path.join(case, "0/U")
    txt = open(p).read().replace("""    bottom
    {
        type            slip;
    }""", """    bottom
    {
        type            nhMovingSlipVelocity;
        value           uniform (0 0 0);
    }""")
    open(p, "w").write(txt)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"d {d}\nx0 {s['x0']}\nut {s['ut']}\na0 {s['a0']}\nb 0.395\n")
    allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)


def hansen_svendsen_case(case, test="061071", nLayers=10, dx=0.02, dt=0.01, nPeriods=30,
                         breaking=None, conv=None, H=None, zone=8.0, flat=6.0,
                         xEnd=11.0, shelf=3.0, sigmaMax=20.0, zgrade=1.0):
    """Hansen & Svendsen (1979): regular waves shoaling and breaking on a
    1:34.26 slope from h0 = 0.36 m. x is measured from the toe of the slope
    (as the lab data); the slope is cut at x = xEnd (no wetting and drying)
    and continues as a flat shelf with a sponge that absorbs the bores."""
    tests = {"061071": dict(T=1.667, H=0.0686), "031041": dict(T=3.333, H=0.0411)}
    T = tests[test]["T"]
    H = H or tests[test]["H"]
    h0, s = 0.36, 1 / 34.26
    hEnd = h0 - s * xEnd
    x0 = -(zone + flat)
    xs = [x0, 0.0, xEnd, xEnd + shelf]
    hs = [h0, h0, hEnd, hEnd]
    dy = dx
    gx = [round(-1.0 + 0.1 * i, 3) for i in range(int(round((xEnd - 0.2 + 1.0) / 0.1)) + 1)]
    gtxt = "\n".join(f"    G{i:03d}    ({x:g} {dy/2:g} 0);" for i, x in enumerate(gx))
    brk = ""
    if breaking:
        brk = f"""
breaking
{{
    model       HFA;
    alpha       {breaking.get('alpha', 0.6):g};
    beta        {breaking.get('beta', 0.3):g};
}}
"""
    constant_files(case, nh_extra=f"""
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
        sigmaMax    {sigmaMax:g};
        exponent    2;
        waveHeight  {H:g};
        wavePeriod  {T:g};
        depth       {h0:g};
        direction   (1 0 0);
        rampTime    {2*T:g};
    }}
    outlet                   // damps the bores, keeps the water volume
    {{
        type        sponge;
        start       ({xEnd + 0.2:g} 0 0);
        end         ({xEnd + shelf:g} 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    2;
        relaxEta    false;
    }}
}}
{brk}""")
    endTime = nPeriods * T
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt, writeInterval=nsteps, nCorr=3,
                 nNonOrth=1, pTol=1e-10, conv=conv or (CONV_BREAKING if breaking else CONV))
    blockmesh_profile(case, xs, hs, dy, dx, nLayers, zgrade=zgrade)
    common_fields(case)
    with open(os.path.join(case, "params.txt"), "w") as f:
        f.write(f"test {test}\nT {T}\nH {H}\n")
    allrun(case, """
    runApplication blockMesh
    runApplication NHFOAM
    """)
