Info<< "Reading NHFoamProperties\n" << endl;

IOdictionary nhProperties
(
    IOobject
    (
        "NHFoamProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

#include "readGravitationalAcceleration.H"

const dimensionedScalar magg(mag(g));
const scalar gMag = magg.value();

if (gMag < SMALL || (g.value() & vector(0, 0, -1)) < 0.9999*gMag)
{
    FatalErrorInFunction
        << "Gravity must point in the -z direction, g = " << g.value()
        << exit(FatalError);
}

// Column structure of the sigma mesh
sigmaColumns columns(mesh, nhProperties);

const label surfaceID = columns.surfacePatchID();

// Smallest water depth used to place the mesh points
const scalar minDepth
(
    nhProperties.getOrDefault<scalar>("minDepth", 1e-4)
);


Info<< "Reading field eta\n" << endl;

volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Free-surface elevation per water column (the value of the surface cell)
scalarField etaCol(columns.columnValues(eta.primitiveField()));

// Time levels and coefficients of the free-surface time integration
scalarField etaN(etaCol);       // eta at time level n
scalarField etaNm1(etaCol);     // eta at time level n-1
scalarField etaHat(etaCol);     // history part of the time derivative
scalarField etaPred(etaCol);    // extrapolated eta used to place the mesh
scalar aEff = runTime.deltaTValue();   // dt/coefficient of eta^{n+1}
label nEtaLevels = 0;

columns.setCellValues(etaCol, eta.primitiveFieldRef());
eta.correctBoundaryConditions();


// Move the mesh with the free surface (false: fixed mesh, for linear
// small-amplitude studies only)
const bool moveMesh(nhProperties.getOrDefault<bool>("moveMesh", true));

// Place the mesh on the initial free surface
if (moveMesh)
{
    tmp<pointField> tnewPoints(columns.newPoints(etaCol, minDepth));

    const scalar maxShift =
        gMax(mag(tnewPoints() - mesh.points())());

    if (maxShift > SMALL)
    {
        Info<< "Moving the sigma mesh onto the initial free surface"
            << " (largest point shift " << maxShift << " m)\n" << endl;

        mesh.movePoints(tnewPoints());
    }
}


// Free-surface time integration: follows the ddt scheme of U
bool bdf2 = false;
{
    ITstream& is = mesh.ddtScheme("ddt(U)");
    const word ddtName(is);

    // Optional coefficient of the ddt flux correction (ddtPhiCoeff)
    scalar ddtPhiCoeff = -1;
    if (is.good() && !is.eof())
    {
        is >> ddtPhiCoeff;
    }

    if (ddtName == "backward")
    {
        bdf2 = true;
    }
    else if (ddtName != "Euler")
    {
        WarningInFunction
            << "ddt scheme " << ddtName << " is not supported for the free"
            << " surface: using Euler for eta (use Euler or backward)"
            << endl;
    }

    if (ddtPhiCoeff < 0)
    {
        WarningInFunction
            << "No ddtPhiCoeff given for ddt(U): OpenFOAM then limits the"
            << " ddt flux correction, which damps waves noticeably." << nl
            << "    Use e.g. 'default backward 1;' in ddtSchemes." << endl;
    }

    Info<< "Free-surface time integration: "
        << (bdf2 ? "backward (BDF2)" : "Euler") << nl << endl;
}


Info<< "Reading field p_rgh\n" << endl;

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

if (!isA<mixedFvPatchScalarField>(p_rgh.boundaryField()[surfaceID]))
{
    FatalErrorInFunction
        << "p_rgh must use the 'mixed' condition on the free-surface patch "
        << mesh.boundaryMesh()[surfaceID].name() << " (its coefficients are"
        << " set by the solver)" << exit(FatalError);
}

mesh.setFluxRequired(p_rgh.name());


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Calculating face flux field phi\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

// Face velocity, used to rebuild the flux after the mesh has moved
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(U)
);


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"


// Derived pressures (kinematic, m2/s2)
//   q: non-hydrostatic part, p: total pressure p/rho
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh - magg*eta
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + (g & mesh.C())
);


// Sponge and wave-making zones
forcingZones zones(nhProperties, gMag);

volScalarField::Internal zoneSigma
(
    IOobject
    (
        "zoneSigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0)
);

volVectorField::Internal zoneSigmaU
(
    IOobject
    (
        "zoneSigmaU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity/dimTime, Zero)
);

scalarField colSigma(columns.nColumns(), 0);
scalarField colSigmaEta(columns.nColumns(), 0);

const scalar volume0 = gSum(columns.columnArea()*etaCol);
