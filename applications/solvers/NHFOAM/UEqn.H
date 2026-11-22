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
        scalar Cd = botCd;

        if (botManning > 0)
        {
            Cd = gMag*sqr(botManning)/Foam::cbrt(D);
        }
        else if (botKs > 0)
        {
            // Log law at the bottom cell centre, z1 half its height
            const scalar z1 = 0.5*Vc[c]/magSb[i];
            Cd = sqr(0.41/max(Foam::log(30*z1/botKs), scalar(1)));
        }

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

// Breaking eddy viscosity: momentum diffusion through the side faces (and
// the layer faces with verticalMixing), not through closed faces of dry
// columns
if (breakingVisc)
{
    surfaceScalarField nuBreakingf
    (
        "nuBreakingf",
        fvc::interpolate(nuBreaking)*brkFaceWeight
    );

    if (wetDry.active())
    {
        const boolList& open = wetDry.faceOpen();
        scalarField& nf = nuBreakingf.primitiveFieldRef();

        forAll(nf, facei)
        {
            if (!open[facei])
            {
                nf[facei] = 0;
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                const boolList& open = wetDry.patchFaceOpen()[patchi];
                scalarField& pnf = nuBreakingf.boundaryFieldRef()[patchi];

                forAll(pnf, i)
                {
                    if (!open[i])
                    {
                        pnf[i] = 0;
                    }
                }
            }
        }
    }

    UEqn -= fvm::laplacian(nuBreakingf, U);
}

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
