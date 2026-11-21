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

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  + fvm::Sp(zoneSigma, U)
 ==
    fvOptions(U)
  + zoneSigmaU
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_rgh));

    fvOptions.correct(U);
}
