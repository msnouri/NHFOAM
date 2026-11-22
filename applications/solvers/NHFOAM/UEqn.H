// Momentum predictor.  Gravity is absorbed in p_rgh = g eta + q.

MRF.correctBoundaryVelocity(U);

if (zones.active())
{
    zones.cellTargets
    (
        mesh.C().primitiveField(),
        runTime.value(),
        zoneSigma.field(),
        zoneSigmaU.field()
    );
}

// Quadratic bottom friction on the bottom cells (lagged |u_h|)
if (botFrictionDictPtr)
{
    const fvPatch& bp = mesh.boundary()[columns.bottomPatchID()];
    const labelUList& fc = bp.faceCells();
    const scalarField& magSb = bp.magSf();
    const scalarField& Vc = mesh.V();
    const labelList& cellCol = columns.cellColumn();
    const scalarField& zbCol = columns.columnBottomZ();

    botDrag.field() = 0;

    forAll(fc, i)
    {
        const label c = fc[i];
        const scalar uh = Foam::sqrt(sqr(U[c].x()) + sqr(U[c].y()));
        const scalar D = max(etaCol[cellCol[c]] - zbCol[cellCol[c]], minDepth);
        const scalar Cd =
            botManning > 0
          ? gMag*sqr(botManning)/Foam::cbrt(D)
          : botCd;

        botDrag[c] = Cd*uh*magSb[i]/Vc[c];
    }
}

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  + fvm::Sp(zoneSigma, U)
  + fvm::Sp(wdSink, U)
  + fvm::Sp(botDrag, U)
 ==
    fvOptions(U)
  + zoneSigmaU
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    if (interfaceP)
    {
        solve(UEqn == -gradPPtr());
    }
    else
    {
        solve(UEqn == -fvc::grad(p_rgh));
    }

    fvOptions.correct(U);
}
