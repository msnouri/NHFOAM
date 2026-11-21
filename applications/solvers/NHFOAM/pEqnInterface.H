// Pressure equation with the pressure on the layer interfaces (NHWAVE
// layout, see interfacePressure.H): p_rgh = g eta + q on the sigma levels,
// q = 0 on the free surface, eta^{n+1} solved together with it.

{
    volScalarField& pI = pIPtr();
    volVectorField& gradP = gradPPtr();
    interfacePressure& pInt = pIntPtr();

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    if (ddtCorr)
    {
        phiHbyA += MRF.zeroFilter(rAUf*fvc::ddtCorr(U, Uf));
    }

    MRF.makeRelative(phiHbyA);

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    // Boundary gradients of the (derived) cell pressure on flux boundaries
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Free-surface coefficients; breaking columns get a large vertical
    // mobility (hydrostatic front approximation)
    const scalarField etaHatEff(etaHat + aEff*colSigmaEta);
    const scalarField colS(1 + aEff*colSigma);
    scalarField vFactor(columns.nColumns(), 1);

    if (hfa)
    {
        forAll(colBreaking, coli)
        {
            if (colBreaking[coli])
            {
                vFactor[coli] = hfaFactor;
            }
        }
    }

    // Lagged terms are updated in every non-orthogonal corrector
    while (pimple.correctNonOrthogonal())
    {
        pInt.solve
        (
            pI, psCol, p_rgh, gradP,
            rAU, rAUf, HbyA, phiHbyA, U,
            etaHatEff, colS, aEff, gMag, vFactor,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        pInt.cellPressure(pI, psCol, p_rgh);

        if (pimple.finalNonOrthogonalIter())
        {
            pInt.fluxes(p_rgh, phiHbyA, U, phi);
        }

        gradP = pInt.grad(pI, psCol, p_rgh);
    }

    // New free-surface elevation from the (absolute) surface flux
    {
        const scalarField& phis = phi.boundaryField()[surfaceID];
        const scalarField& Ah = columns.columnArea();

        forAll(etaCol, coli)
        {
            etaCol[coli] =
            (
                etaHat[coli]
              + aEff*(phis[coli]/Ah[coli] + colSigmaEta[coli])
            )/(1 + aEff*colSigma[coli]);
        }
    }

    #include "continuityErrs.H"

    U = HbyA - rAU*gradP;
    U.correctBoundaryConditions();
    fvOptions.correct(U);

    // Face velocity (absolute) for rebuilding the flux after mesh motion
    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf += n*(phi/mesh.magSf() - (n & Uf));
    }

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi, U);
}
