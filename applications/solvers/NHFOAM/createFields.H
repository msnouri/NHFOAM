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

// Restart state, written with the fields (endTimeStep.H): the number of
// time levels of eta and the averaging time in uniform/NHFOAMState, eta of
// the previous level (etaOld), the face velocities of the last two levels
// (Uf, UfN), the surface pressure unknowns (pSurf), and the breaking flags
// and ages.  With it a restart continues as if the run had not stopped
// (bitwise with binary output); without it the run starts afresh: Euler
// first step, mesh placed on eta.
IOdictionary nhState
(
    IOobject
    (
        "NHFOAMState",
        runTime.timeName(),
        "uniform",
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    )
);
const bool restartState(nhState.found("nEtaLevels"));

volScalarField etaOld
(
    IOobject
    (
        "etaOld",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    eta
);

if (restartState)
{
    nEtaLevels = nhState.get<label>("nEtaLevels");
    etaN = columns.columnValues(etaOld.primitiveField());

    Info<< "Restart with the saved state: " << nEtaLevels
        << " time levels, eta of the previous level and the face"
        << " velocities read; the mesh is kept\n" << endl;
}


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

    // On a restart the mesh holds the displaced bottom: recover the
    // reference bottom of the displacement
    if (restartState)
    {
        columns.setLineBottomZ0
        (
            columns.lineBottomZ0()
          - bedMotion.displacement(columns.lineCentre(), runTime.value())
        );
    }

    columns.setLineBottomZ
    (
        columns.lineBottomZ0()
      + bedMotion.displacement(columns.lineCentre(), runTime.value())
    );
}

// Place the mesh on the initial free surface (and bottom); a restart keeps
// the mesh it was written with
if (moveMesh && !restartState)
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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(U)()
);

// Face velocities of the last two time levels, n and n-1: the convecting
// flux of the first outer corrector is extrapolated from them
// (fluxExtrapolation)
surfaceVectorField UfN
(
    IOobject
    (
        "UfN",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uf
);
surfaceVectorField UfNm1("UfNm1", UfN);


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


// Wave breaking (updateBreaking.H)
//     breaking
//     {
//         model       HFA;    // none | HFA | eddyViscosity | hybrid
//         onset       etaT;   // etaT | kinematic
//         alpha       0.6;    // onset: d(eta)/dt > alpha sqrt(g D)
//         Bkin        0.85;   // kinematic onset: u_s.n > Bkin C, crest
//         crestSlope  0.01;   //   speed C from the front faces steeper
//                             //   than crestSlope (updateBreaking.H)
//         beta        0.3;    // persistence: d(eta)/dt > beta sqrt(g D)
//         persistence threshold;  // threshold: a breaking column stops
//                             //   below beta sqrt(g D); crest: it stops
//                             //   when the surface no longer rises (the
//                             //   crest has passed), and beta only lets
//                             //   neighbours join
//         betaFinal   0.15;   // optional: persistence threshold falling
//         Tstar       5;      //   from beta to betaFinal over the breaking
//                             //   age T* = Tstar sqrt(D/g)
//         delta       1.2;    // eddy viscosity mixing length
//         verticalMixing 0;   // eddy viscosity between the layers: a
//                             // fraction of nu_b, or yes (1) / no (0)
//         hydrostaticFactor 1000;
//     }
// A column starts breaking at the onset criterion and keeps breaking while
// d(eta)/dt stays above the persistence threshold; a neighbouring column
// joins (with the breaker's age) while d(eta)/dt > beta sqrt(g D).
//   HFA            hydrostatic front approximation (Smit et al. 2013): in
//                  breaking columns the non-hydrostatic pressure is
//                  suppressed by coupling the cells of the column (and the
//                  surface face) with a large vertical mobility in the
//                  pressure equation, so p_rgh stays uniform in the vertical
//                  and the fluxes through the layer faces follow from
//                  continuity.
//   eddyViscosity  breaking eddy viscosity after Kennedy et al. (2000),
//                  nu_b = B delta^2 D d(eta)/dt in the breaking columns,
//                  B = min(max(d(eta)/dt/eta_t* - 1, 0), 1) with eta_t* the
//                  persistence threshold, as momentum diffusion through the
//                  side faces (and the layer faces with verticalMixing); the
//                  non-hydrostatic pressure is kept.
//   hybrid         both.
const dictionary breakingDict(nhProperties.subOrEmptyDict("breaking"));
const word breakingModel
(
    breakingDict.getOrDefault<word>("model", "none")
);

if
(
    breakingModel != "none" && breakingModel != "HFA"
 && breakingModel != "eddyViscosity" && breakingModel != "hybrid"
)
{
    FatalIOErrorInFunction(nhProperties)
        << "Unknown breaking model " << breakingModel
        << ". Valid: none HFA eddyViscosity hybrid" << exit(FatalIOError);
}

const bool breakingOn(breakingModel != "none");
const bool hfa(breakingModel == "HFA" || breakingModel == "hybrid");
const bool breakingVisc
(
    breakingModel == "eddyViscosity" || breakingModel == "hybrid"
);

const word breakingOnset
(
    breakingDict.getOrDefault<word>("onset", "etaT")
);

if (breakingOnset != "etaT" && breakingOnset != "kinematic")
{
    FatalIOErrorInFunction(nhProperties)
        << "Unknown breaking onset " << breakingOnset
        << ". Valid: etaT kinematic" << exit(FatalIOError);
}

const bool kinematicOnset(breakingOnset == "kinematic");
const scalar hfaAlpha(breakingDict.getOrDefault<scalar>("alpha", 0.6));
const scalar brkBkin(breakingDict.getOrDefault<scalar>("Bkin", 0.85));
const scalar brkCrestSlope
(
    breakingDict.getOrDefault<scalar>("crestSlope", 0.01)
);
const scalar hfaBeta(breakingDict.getOrDefault<scalar>("beta", 0.3));
const word brkPersistence
(
    breakingDict.getOrDefault<word>("persistence", "threshold")
);

if (brkPersistence != "threshold" && brkPersistence != "crest")
{
    FatalIOErrorInFunction(nhProperties)
        << "Unknown breaking persistence " << brkPersistence
        << ". Valid: threshold crest" << exit(FatalIOError);
}

const bool brkToCrest(brkPersistence == "crest");
const bool brkAgeing(breakingDict.found("betaFinal"));
const scalar brkBetaFinal
(
    breakingDict.getOrDefault<scalar>("betaFinal", hfaBeta)
);
const scalar brkTstar(breakingDict.getOrDefault<scalar>("Tstar", 5));
const scalar brkDelta(breakingDict.getOrDefault<scalar>("delta", 1.2));
// verticalMixing: fraction of nu_b on the layer faces (yes = 1, no = 0)
scalar brkVertical = 0;
if (breakingDict.found("verticalMixing"))
{
    ITstream& is = breakingDict.lookup("verticalMixing");
    const token t(is);

    if (t.isNumber())
    {
        brkVertical = t.number();
    }
    else if (t.isWord())
    {
        brkVertical = (Switch(t.wordToken()) ? 1 : 0);
    }
    else
    {
        FatalIOErrorInFunction(breakingDict)
            << "verticalMixing: a switch (yes/no) or a fraction expected"
            << exit(FatalIOError);
    }
}
const scalar hfaFactor
(
    breakingDict.getOrDefault<scalar>("hydrostaticFactor", 1000)
);

if (breakingOn)
{
    Info<< "Wave breaking: " << breakingModel << ", onset ";
    if (kinematicOnset)
    {
        Info<< "u_s.n > " << brkBkin << " C (crest speed)";
    }
    else
    {
        Info<< "d(eta)/dt > " << hfaAlpha << " sqrt(g D)";
    }
    Info<< ", persistence d(eta)/dt > ";
    if (brkAgeing)
    {
        Info<< "(" << hfaBeta << " to " << brkBetaFinal << " over "
            << brkTstar << " sqrt(D/g)) sqrt(g D)";
    }
    else
    {
        Info<< hfaBeta << " sqrt(g D)";
    }
    if (brkToCrest)
    {
        Info<< " to join, breaking columns until the crest has passed";
    }
    if (breakingVisc)
    {
        Info<< ", eddy viscosity delta " << brkDelta
            << (brkVertical > 0 ? " (horizontal, vertical x " : " (horizontal");
        if (brkVertical > 0)
        {
            Info<< brkVertical;
        }
        Info<< ")";
    }
    Info<< nl << endl;
}

boolList colBreaking(columns.nColumns(), false);
scalarField colAge(columns.nColumns(), 0);
scalarField colNuB(columns.nColumns(), 0);

volScalarField breaking
(
    IOobject
    (
        "breaking",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        breakingOn ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField breakingAge
(
    IOobject
    (
        "breakingAge",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        breakingOn ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimTime, 0),
    zeroGradientFvPatchScalarField::typeName
);

if (breakingOn && restartState)
{
    IOobject bHeader
    (
        "breaking",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );
    IOobject aHeader
    (
        "breakingAge",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if
    (
        bHeader.typeHeaderOk<volScalarField>(true)
     && aHeader.typeHeaderOk<volScalarField>(true)
    )
    {
        const volScalarField b(bHeader, mesh);
        const volScalarField a(aHeader, mesh);
        const scalarField bCol(columns.columnValues(b.primitiveField()));

        colAge = columns.columnValues(a.primitiveField());
        forAll(colBreaking, coli)
        {
            colBreaking[coli] = (bCol[coli] > 0.5);
        }
        breaking.primitiveFieldRef() = b.primitiveField();
        breaking.correctBoundaryConditions();
    }
}

// Breaking ages exchanged across processor boundaries (updateBreaking.H)
volScalarField brkExchange
(
    IOobject
    (
        "brkExchange",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimensionedScalar("minusOne", dimless, -1),
    zeroGradientFvPatchScalarField::typeName
);

// Breaking eddy viscosity in the cells, and the face weights of its
// diffusion: 1 on the side faces (including coupled patches), 1 or 0 on the
// layer faces, 0 on the other boundaries
volScalarField nuBreaking
(
    IOobject
    (
        "nuBreaking",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        breakingVisc ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimViscosity, 0),
    zeroGradientFvPatchScalarField::typeName
);

surfaceScalarField brkFaceWeight
(
    IOobject
    (
        "brkFaceWeight",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

if (breakingVisc)
{
    scalarField& w = brkFaceWeight.primitiveFieldRef();
    w = 1;

    const labelListList& lf = columns.columnLayerFaces();
    forAll(lf, coli)
    {
        forAll(lf[coli], i)
        {
            w[lf[coli][i]] = brkVertical;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        brkFaceWeight.boundaryFieldRef()[patchi] =
            (mesh.boundary()[patchi].coupled() ? 1 : 0);
    }
}


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

// Surface pressure unknowns, written for restarts
volScalarField pSurf
(
    IOobject
    (
        "pSurf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        interfaceP ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(p_rgh.dimensions(), Zero),
    zeroGradientFvPatchScalarField::typeName
);

if (interfaceP && restartState && pSurf.headerOk())
{
    psCol = columns.columnValues(pSurf.primitiveField());
}

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
//     tau_b/rho = Cd |u_h| u
// with Cd constant (Cd), g n^2/D^(1/3) (Manning n, manning) or the log law
// [kappa/ln(30 z1/ks)]^2 at the bottom cell centre z1 (roughness height ks,
// as NHWAVE's Zob), applied as an implicit sink Cd |u_h| A_b/V in the bottom
// cells
const dictionary* botFrictionDictPtr =
    nhProperties.findDict("bottomFriction");
scalar botCd = 0;
scalar botManning = 0;
scalar botKs = 0;

if (botFrictionDictPtr)
{
    botCd = botFrictionDictPtr->getOrDefault<scalar>("Cd", 0);
    botManning = botFrictionDictPtr->getOrDefault<scalar>("manning", 0);
    botKs = botFrictionDictPtr->getOrDefault<scalar>("ks", 0);

    Info<< "Bottom friction: ";
    if (botManning > 0)
    {
        Info<< "Manning n = " << botManning;
    }
    else if (botKs > 0)
    {
        Info<< "log law, ks = " << botKs;
    }
    else
    {
        Info<< "Cd = " << botCd;
    }
    Info<< nl << endl;
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

// Time averages from timeStart on (endTimeStep.H), written with the fields:
//     averaging { timeStart 95; fields (k); }
// UMean (cells; with uniform sigma layers the layer average of a column is
// the mean depth-averaged velocity), etaMean, and <field>Mean of the listed
// scalar fields
const dictionary* avgDictPtr = nhProperties.findDict("averaging");
const scalar avgStart
(
    avgDictPtr ? avgDictPtr->get<scalar>("timeStart") : GREAT
);
scalar avgTime = nhState.getOrDefault<scalar>("avgTime", 0);

autoPtr<volVectorField> UMeanPtr;
autoPtr<volScalarField> etaMeanPtr;

if (avgDictPtr)
{
    UMeanPtr.reset
    (
        new volVectorField
        (
            IOobject
            (
                "UMean",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector(dimVelocity, Zero)
        )
    );
    etaMeanPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "etaMean",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimLength, Zero)
        )
    );

    Info<< "Time averages UMean, etaMean from t = " << avgStart << nl
        << endl;
}

// Time averages of further scalar fields, e.g. averaging { fields (k nut); }:
// <name>Mean
const wordList avgScalarNames
(
    avgDictPtr
  ? avgDictPtr->getOrDefault<wordList>("fields", wordList())
  : wordList()
);
PtrList<volScalarField> avgScalarMeans(avgScalarNames.size());

forAll(avgScalarNames, i)
{
    const volScalarField& f =
        mesh.lookupObject<volScalarField>(avgScalarNames[i]);

    avgScalarMeans.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                avgScalarNames[i] + "Mean",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(f.dimensions(), Zero)
        )
    );
}

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

// Highest free surface of each column over the time steps in which it is wet
// (maximum inundation, run-up), written as the field etaMax:
//     etaMax yes;
// Columns that were never wet keep their bed elevation, so etaMax - z_b > 0
// marks the flooded area.  A restart continues the record.
const bool etaMaxOn(nhProperties.getOrDefault<Switch>("etaMax", false));
autoPtr<volScalarField> etaMaxPtr;
scalarField etaMaxCol;

if (etaMaxOn)
{
    etaMaxPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "etaMax",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("etaMax", dimLength, -GREAT),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    scalarField start(etaCol);
    if (wetDry.active())
    {
        const boolList& wet = wetDry.colWet();
        const scalarField& zb = columns.columnBottomZ();
        forAll(start, coli)
        {
            if (!wet[coli])
            {
                start[coli] = zb[coli];
            }
        }
    }
    etaMaxCol = max(columns.columnValues(etaMaxPtr->primitiveField()), start);
    columns.setCellValues(etaMaxCol, etaMaxPtr->primitiveFieldRef());
    etaMaxPtr->correctBoundaryConditions();

    Info<< "Recording the highest free surface (etaMax)" << nl << endl;
}
