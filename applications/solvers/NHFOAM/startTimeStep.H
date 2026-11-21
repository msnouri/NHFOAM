// Shift the free-surface time levels and set the coefficients of
//     d(eta)/dt = (eta^{n+1} - etaHat)/aEff
// Euler on the first step (and after a restart), then backward (BDF2) if
// ddt(U) is backward, with variable time-step coefficients.
{
    etaNm1 = etaN;
    etaN = etaCol;
    ++nEtaLevels;

    const scalar dt = runTime.deltaTValue();
    const scalar dt0 = runTime.deltaT0Value();

    if (bdf2 && nEtaLevels >= 2)
    {
        const scalar coefft = 1 + dt/(dt + dt0);
        const scalar coefft00 = dt*dt/(dt0*(dt + dt0));
        const scalar coefft0 = coefft + coefft00;

        aEff = dt/coefft;
        etaHat = (coefft0*etaN - coefft00*etaNm1)/coefft;
        etaPred = etaN + (dt/dt0)*(etaN - etaNm1);
    }
    else
    {
        aEff = dt;
        etaHat = etaN;
        etaPred = etaN;
    }

    // Bottom position at the new time
    if (bedMotion.active())
    {
        columns.setLineBottomZ
        (
            columns.lineBottomZ0()
          + bedMotion.displacement(columns.lineCentre(), runTime.value())
        );
    }

    // Breaking columns (hydrostatic front approximation)
    if (hfa)
    {
        #include "updateBreaking.H"
    }

    // Relaxation zones on the free surface at the new time
    if (zones.active())
    {
        zones.columnTargets
        (
            columns.columnCentre(),
            runTime.value(),
            colSigma,
            colSigmaEta
        );
    }
}
