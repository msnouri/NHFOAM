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

// Wetting and drying (NHFOAMProperties/wettingDrying; off if absent)
wettingDrying wetDry(mesh, columns, nhProperties);

// Smallest water depth used to place the mesh points: the film of the dry
// columns with wetting and drying
const scalar minDepth
(
    wetDry.active()
  ? wetDry.minDepth()
  : nhProperties.getOrDefault<scalar>("minDepth", 1e-4)
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

// Free-surface elevation per water column (the value of the surface cell);
// dry columns get the film
scalarField etaCol(columns.columnValues(eta.primitiveField()));
wetDry.initialise(etaCol);

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

// Face velocities of the last two time levels, n and n-1: the convecting
// flux of the first outer corrector is extrapolated from them
// (fluxExtrapolation)
surfaceVectorField UfN("UfN", Uf);
surfaceVectorField UfNm1("UfNm1", Uf);


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


// Vertical location of the pressure unknowns
//     interface : layer interfaces (sigma levels) with q = 0 exactly on the
//                 free surface, as in NHWAVE (see interfacePressure.H);
//                 the default
//     cell      : cell centres, surface condition on the top faces
const word pressureLocation
(
    nhProperties.getOrDefault<word>("pressureLocation", "interface")
);

if (pressureLocation != "cell" && pressureLocation != "interface")
{
    FatalIOErrorInFunction(nhProperties)
        << "Unknown pressureLocation " << pressureLocation
        << ". Valid: cell interface" << exit(FatalIOError);
}

const bool interfaceP(pressureLocation == "interface");

Info<< "Pressure unknowns: "
    << (
           interfaceP
         ? "layer interfaces, q = 0 on the free surface"
         : "cell centres"
       )
    << nl << endl;

autoPtr<interfacePressure> pIntPtr;   // interface layout
autoPtr<volScalarField> pIPtr;        // p_rgh on the bottom interface of
                                      // each cell
autoPtr<volVectorField> gradPPtr;     // cell gradient of p_rgh
scalarField psCol(gMag*etaCol);       // surface pressure g eta per column

if (interfaceP)
{
    forAll(p_rgh.boundaryField(), patchi)
    {
        const fvPatchScalarField& pp = p_rgh.boundaryField()[patchi];

        if (patchi != surfaceID && !pp.coupled() && pp.fixesValue())
        {
            FatalErrorInFunction
                << "pressureLocation interface: patch " << pp.patch().name()
                << " fixes p_rgh (" << pp.type() << ")." << nl
                << "    Only flux conditions (fixedFluxPressure,"
                << " zeroGradient) are supported on walls and bottom."
                << exit(FatalError);
        }
    }

    pIntPtr.reset(new interfacePressure(mesh, columns));

    IOobject pIHeader
    (
        "pI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (pIHeader.typeHeaderOk<volScalarField>(true))
    {
        Info<< "Reading field pI\n" << endl;
        pIPtr.reset(new volScalarField(pIHeader, mesh));
    }
    else
    {
        Info<< "Field pI (p_rgh on the bottom interface of each cell)"
            << " initialised from p_rgh\n" << endl;
        pIHeader.readOpt() = IOobject::NO_READ;
        pIPtr.reset
        (
            new volScalarField
            (
                pIHeader,
                p_rgh,
                calculatedFvPatchScalarField::typeName
            )
        );
    }

    pIntPtr->cellPressure(pIPtr(), psCol, p_rgh);

    gradPPtr.reset
    (
        new volVectorField
        (
            "gradP",
            pIntPtr->grad(pIPtr(), psCol, p_rgh)
        )
    );

    if (wetDry.active())
    {
        pIntPtr->setFaceMask(wetDry.faceOpen(), wetDry.patchFaceOpen());
        gradPPtr() = pIntPtr->grad(pIPtr(), psCol, p_rgh);
    }
}

if (wetDry.active() && !interfaceP)
{
    FatalErrorInFunction
        << "wettingDrying needs pressureLocation interface"
        << exit(FatalError);
}

// Columns that changed state (flooded or dried) at the start of the step
boolList wdChanged(columns.nColumns(), false);

// Wetting and drying: the side fluxes carry the upwind (TVD) water depth
// through scale factors (wettingDrying::faceFactors), evaluated on the
// surface of the current mesh
scalarField etaMesh(etaCol);
scalarField faceFactor(mesh.nInternalFaces(), 1.0);
List<scalarField> patchFaceFactor(mesh.boundary().size());
forAll(patchFaceFactor, patchi)
{
    patchFaceFactor[patchi].setSize(mesh.boundary()[patchi].size(), 1.0);
}

if (wetDry.scaleFaceDepth())
{
    pIntPtr->setFaceFactor(faceFactor, patchFaceFactor);
}

// Dry columns: the velocity is held at zero by a large implicit sink
volScalarField::Internal wdSink
(
    IOobject
    (
        "wdSink",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0)
);

scalarField dryFlag(mesh.nCells(), 0);

// Quadratic bottom friction on the bottom cells,
//     tau_b/rho = Cd |u_h| u   (Cd constant, or g n^2/D^(1/3) with Manning n)
// applied as an implicit sink Cd |u_h| A_b/V in the bottom cells
const dictionary* botFrictionDictPtr =
    nhProperties.findDict("bottomFriction");
scalar botCd = 0;
scalar botManning = 0;

if (botFrictionDictPtr)
{
    botCd = botFrictionDictPtr->getOrDefault<scalar>("Cd", 0);
    botManning = botFrictionDictPtr->getOrDefault<scalar>("manning", 0);

    Info<< "Bottom friction: "
        << (botManning > 0 ? "Manning n = " : "Cd = ")
        << (botManning > 0 ? botManning : botCd) << nl << endl;
}

volScalarField::Internal botDrag
(
    IOobject
    (
        "botDrag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0)
);

// Shoreline record (wetting and drying): time, wet columns, highest wet
// bed, highest free surface over ground above the still-water level
autoPtr<OFstream> shorelineFilePtr;
const scalarField zbInitial(columns.columnBottomZ());
const scalar wdStillWaterLevel
(
    nhProperties.getOrDefault<scalar>("stillWaterLevel", 0)
);

if (wetDry.active() && Pstream::master())
{
    const fileName dir
    (
        runTime.globalPath()/"postProcessing"/"wettingDrying"
       /runTime.timeName()
    );
    mkDir(dir);
    shorelineFilePtr.reset(new OFstream(dir/"shoreline.dat"));
    shorelineFilePtr()
        << "# time  wetColumns  maxWetBed  maxEtaOnLand" << endl;
}
