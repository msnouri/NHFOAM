Info<< "Reading NHFOAMProperties\n" << endl;

IOdictionary nhProperties
(
    IOobject
    (
        "NHFOAMProperties",
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

// Prescribed bottom motion (uplift, solid slide)
bottomMotion bedMotion(nhProperties);

if (bedMotion.active())
{
    if (!moveMesh)
    {
        FatalErrorInFunction
            << "bottomMotion needs a moving mesh (moveMesh true)"
            << exit(FatalError);
    }

    columns.setLineBottomZ
    (
        columns.lineBottomZ0()
      + bedMotion.displacement(columns.lineCentre(), runTime.value())
    );
}

// Place the mesh on the initial free surface (and bottom)
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
        columns.updateGeometry();
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

// Water volume at the start (changes only through relaxation zones)
const scalar volume0 =
    gSum(columns.columnArea()*(etaCol - columns.columnBottomZ()));


// Wave breaking: hydrostatic front approximation (Smit et al. 2013)
//     breaking  { model HFA; alpha 0.6; beta 0.3; hydrostaticFactor 1000; }
// A column starts breaking when d(eta)/dt > alpha sqrt(g D) and keeps (or
// spreads to a neighbouring column) while d(eta)/dt > beta sqrt(g D).
// In breaking columns the non-hydrostatic pressure is suppressed by coupling
// the cells of the column (and the surface face) with a large vertical
// mobility in the pressure equation, so p_rgh stays uniform in the vertical
// (hydrostatic) and the fluxes through the layer faces follow from
// continuity.
const dictionary breakingDict(nhProperties.subOrEmptyDict("breaking"));
const word breakingModel
(
    breakingDict.getOrDefault<word>("model", "none")
);
const bool hfa(breakingModel == "HFA");

if (breakingModel != "none" && !hfa)
{
    FatalIOErrorInFunction(nhProperties)
        << "Unknown breaking model " << breakingModel
        << ". Valid: none HFA" << exit(FatalIOError);
}

const scalar hfaAlpha(breakingDict.getOrDefault<scalar>("alpha", 0.6));
const scalar hfaBeta(breakingDict.getOrDefault<scalar>("beta", 0.3));
const scalar hfaFactor
(
    breakingDict.getOrDefault<scalar>("hydrostaticFactor", 1000)
);

if (hfa)
{
    Info<< "Wave breaking: hydrostatic front approximation, alpha "
        << hfaAlpha << ", beta " << hfaBeta << nl << endl;
}

boolList colBreaking(columns.nColumns(), false);

volScalarField breaking
(
    IOobject
    (
        "breaking",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        hfa ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);
