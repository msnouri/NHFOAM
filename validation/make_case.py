#!/usr/bin/env python3
"""Generate 2DV NHFoam cases (standing wave, solitary wave, Beji-Battjes bar)."""
import os, sys, math, textwrap

HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | NHFoam
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

def common_fields(case, walls=("left", "right"), extra_U=None, extra_p=None):
    wallU = "".join(f"    {w}\n    {{\n        type            slip;\n    }}\n" for w in walls)
    wallP = "".join(f"    {w}\n    {{\n        type            fixedFluxPressure;\n        value           uniform 0;\n    }}\n" for w in walls)
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
{wallU}    frontAndBack
    {{
        type            empty;
    }}
}}
""")
    write(case, "0/p_rgh", "volScalarField", f"""
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    // Free surface: coefficients are set by NHFoam every corrector
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
{wallP}    frontAndBack
    {{
        type            empty;
    }}
}}
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
{wallE}    frontAndBack
    {{
        type            empty;
    }}
}}
""")

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
    write(case, "constant/NHFoamProperties", "dictionary", f"""
// Free-surface (top) and bottom patches of the sigma mesh
surfacePatch    surface;
bottomPatch     bottom;

// Smallest water depth used to place the mesh points [m]
minDepth        1e-4;

stillWaterLevel 0;
{nh_extra}
""")

def system_files(case, endTime, deltaT, writeInterval, ddt="backward 1", nCorr=3, nNonOrth=0, pTol=1e-12, extraSolver=""):
    write(case, "system/controlDict", "dictionary", f"""
application     NHFoam;

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
    div(phi,U)      Gauss linearUpwind grad(U);
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
        solver          PCG;
        preconditioner  DIC;
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

def blockmesh_profile(case, xs, hs, dy, dx, Nz, zgrade=1.0):
    """Piecewise-linear bottom: depth hs[i] at xs[i], one block per segment."""
    verts = []
    for i, (x, h) in enumerate(zip(xs, hs)):
        verts += [f"({x:.10g} 0 {-h:.10g})", f"({x:.10g} {dy} {-h:.10g})",
                  f"({x:.10g} 0 0)", f"({x:.10g} {dy} 0)"]
    blocks, surf, bot = [], [], []
    for i in range(len(xs) - 1):
        a = 4 * i; b = 4 * (i + 1)
        # hex: bottom (a0, b0, b1, a1), top (a2, b2, b3, a3)
        nx = max(1, int(round((xs[i + 1] - xs[i]) / dx)))
        blocks.append(f"hex ({a} {b} {b+1} {a+1} {a+2} {b+2} {b+3} {a+3}) ({nx} 1 {Nz}) simpleGrading (1 1 {zgrade})")
        surf.append(f"({a+2} {b+2} {b+3} {a+3})")
        bot.append(f"({a} {a+1} {b+1} {b})")
    n = len(xs) - 1
    left = "(0 2 3 1)"
    right = f"({4*n} {4*n+1} {4*n+3} {4*n+2})"
    fb = []
    for i in range(n):
        a = 4 * i; b = 4 * (i + 1)
        fb.append(f"({a} {b} {b+2} {a+2})")
        fb.append(f"({a+1} {a+3} {b+3} {b+1})")
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
    frontAndBack
    {{
        type empty;
        faces ({' '.join(fb)});
    }}
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
    runApplication NHFoam
    """)
    return dict(k=k, L=L, T=T, dt=dt)

if __name__ == "__main__":
    case = sys.argv[1]
    kh = float(sys.argv[2]); nl = int(sys.argv[3])
    print(standing_case(case, kh, nl))


def bar_case(case, nLayers=10, dx=0.02, dt=0.01, endTime=50.0, H=0.02, T=2.02,
             waves=True, sigmaMax=20.0, writeEvery=None):
    """Luth et al. (1994) / Beji & Battjes (1993) submerged bar, case A
    (flume coordinates of the FUNWAVE benchmark: bar from x = 26 to 37 m)."""
    xs = [0, 26, 32, 34, 37, 54]
    hs = [0.4, 0.4, 0.1, 0.1, 0.4, 0.4]
    dy = 0.02
    gauges = [22, 24, 30.5, 32.5, 33.5, 34.5, 35.7, 37.3, 39.0, 41.0]
    gtxt = "\n".join(f"    G{g:g}    ({g:g} {dy/2:g} 0);" for g in gauges)
    zones = ""
    if waves:
        zones = f"""
forcingZones
{{
    inlet
    {{
        type        wave;
        start       (8 0 0);
        end         (0 0 0);
        sigmaMax    {sigmaMax:g};
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
        start       (46 0 0);
        end         (54 0 0);
        sigmaMax    {sigmaMax:g};
        exponent    2;
    }}
}}
"""
    constant_files(case, nh_extra=f"""
gauges
{{
{gtxt}
}}
{zones}""")
    nsteps = int(round(endTime / dt))
    system_files(case, endTime=endTime, deltaT=dt,
                 writeInterval=writeEvery or nsteps, nCorr=3, nNonOrth=1, pTol=1e-10)
    blockmesh_profile(case, xs, hs, dy, dx, nLayers)
    common_fields(case)
    allrun(case, """
    runApplication blockMesh
    runApplication NHFoam
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
    runApplication NHFoam
    """)
