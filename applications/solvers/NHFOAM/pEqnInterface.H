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

    // Wetting and drying: supercritical side faces take the predicted flux
    // from the upwind cell (as an HLL flux with both wave speeds positive)
    wetDry.upwindSupercritical(phiHbyA, HbyA, U, etaN, gMag);

    // Wetting and drying: faces opened by flooding carry no predicted flux
    // out of the newly wet column
    wetDry.limitOpenedFaces(phiHbyA, wdChanged, etaN, gMag);

    // Wetting and drying: side fluxes with the upwind (TVD) water depth
    // (the implicit coefficients are scaled in interfacePressure)
    if (wetDry.scaleFaceDepth())
    {
        wetDry.faceFactors
        (
            max(etaMesh, wetDry.etaMin()),
            phiHbyA,
            faceFactor,
            patchFaceFactor
        );

        phiHbyA.primitiveFieldRef() *= faceFactor;

        forAll(phiHbyA.boundaryField(), patchi)
        {
            if (phiHbyA.boundaryField()[patchi].coupled())
            {
                phiHbyA.boundaryFieldRef()[patchi] *= patchFaceFactor[patchi];
            }
        }
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

    // Wetting and drying: shallow columns hydrostatic (optional)
    if (wetDry.active() && wetDry.hydrostaticDepth() > 0)
    {
        const scalarField& zb = columns.columnBottomZ();

        forAll(vFactor, coli)
        {
            if (etaN[coli] - zb[coli] < wetDry.hydrostaticDepth())
            {
                vFactor[coli] = max(vFactor[coli], hfaFactor);
            }
        }
    }

    // Lowest admissible surface (positivity limiter)
    const scalarField etaMinWD
    (
        wetDry.active() ? wetDry.etaMin() : tmp<scalarField>(new scalarField())
    );

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
            if (wetDry.active())
            {
                pInt.fluxes
                (
                    p_rgh, phiHbyA, U, phi,
                    &etaMinWD, &etaHatEff, &colS, aEff
                );
            }
            else
            {
                pInt.fluxes(p_rgh, phiHbyA, U, phi);
            }
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

    // Wetting and drying: never below the film (the limiter keeps the
    // columns there up to round-off); dry columns stay at rest
    if (wetDry.active())
    {
        forAll(etaCol, coli)
        {
            etaCol[coli] = max(etaCol[coli], etaMinWD[coli]);
        }
    }

    U = HbyA - rAU*gradP;

    if (wetDry.active())
    {
        forAll(dryFlag, celli)
        {
            if (dryFlag[celli] > 0.5)
            {
                U[celli] = Zero;
            }
        }
    }

    U.correctBoundaryConditions();
    fvOptions.correct(U);

    // Face velocity (absolute) for rebuilding the flux after mesh motion.
    // With the upwind face depth it is the flux over the face depth, not
    // over the mesh face.
    {
        Uf = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());

        if (wetDry.scaleFaceDepth())
        {
            surfaceScalarField un(phi/mesh.magSf());
            scalarField& uni = un.primitiveFieldRef();
            const scalarField nUfi(n.primitiveField() & Uf.primitiveField());

            forAll(uni, facei)
            {
                const scalar ff = faceFactor[facei];
                uni[facei] =
                    ff > 1e-3 ? uni[facei]/ff : (ff > 0 ? nUfi[facei] : 0);
            }

            forAll(un.boundaryField(), patchi)
            {
                if (un.boundaryField()[patchi].coupled())
                {
                    scalarField& unp = un.boundaryFieldRef()[patchi];
                    const scalarField nUfp
                    (
                        n.boundaryField()[patchi] & Uf.boundaryField()[patchi]
                    );
                    const scalarField& ffp = patchFaceFactor[patchi];

                    forAll(unp, i)
                    {
                        unp[i] =
                        (
                            ffp[i] > 1e-3 ? unp[i]/ffp[i]
                          : (ffp[i] > 0 ? nUfp[i] : 0)
                        );
                    }
                }
            }

            Uf += n*(un - (n & Uf));
        }
        else
        {
            Uf += n*(phi/mesh.magSf() - (n & Uf));
        }
    }

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi, U);
}
