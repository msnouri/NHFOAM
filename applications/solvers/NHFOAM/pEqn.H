// Pressure equation for p_rgh = g eta + q with the implicit free surface.
//
// On the surface patch the kinematic condition (with relaxation zones)
//     eta^{n+1} = (etaHat + aEff*(phi_s/A_h + sigma*eta_t))/(1 + aEff*sigma)
// and the dynamic condition p_rgh = g eta^{n+1}, with the surface flux
//     phi_s = phiHbyA_s - rAUf |Sf| deltaCoeff (p_rgh_s - p_rgh_c),
// combine into a Robin condition, written as a mixed condition:
//     p_rgh_s = f refValue + (1 - f) p_rgh_c
//     refValue = g (etaHat + aEff*sigma*eta_t + aEff*phiHbyA_s/A_h)/(1 + s)
//     f        = (1 + s)/(1 + s + beta),   s = aEff*sigma
//     beta     = g aEff rAUf |Sf| deltaCoeff/A_h

{
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

    // Zero-flux walls and bottom (fixedFluxPressure)
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    // Breaking columns: large vertical mobility of the surface face (the
    // layer faces are handled in the pressure matrix below)
    bool anyBreaking = false;

    if (hfa)
    {
        scalarField& rAUfs = rAUf.boundaryFieldRef()[surfaceID];

        forAll(colBreaking, coli)
        {
            if (colBreaking[coli])
            {
                anyBreaking = true;
                rAUfs[coli] *= hfaFactor;
            }
        }
    }

    // Implicit free-surface condition
    {
        mixedFvPatchScalarField& pSurf =
            refCast<mixedFvPatchScalarField>
            (
                p_rgh.boundaryFieldRef()[surfaceID]
            );

        const scalarField& phiHs = phiHbyA.boundaryField()[surfaceID];
        const scalarField& rAUs = rAUf.boundaryField()[surfaceID];
        const scalarField& magSfs = mesh.magSf().boundaryField()[surfaceID];
        const scalarField& deltas = pSurf.patch().deltaCoeffs();
        const scalarField& Ah = columns.columnArea();

        forAll(pSurf, coli)
        {
            const scalar s = aEff*colSigma[coli];
            const scalar beta =
                gMag*aEff*rAUs[coli]*magSfs[coli]*deltas[coli]/Ah[coli];

            pSurf.refValue()[coli] =
                gMag
               *(
                    etaHat[coli]
                  + aEff*colSigmaEta[coli]
                  + aEff*phiHs[coli]/Ah[coli]
                )/(1 + s);

            pSurf.valueFraction()[coli] = (1 + s)/(1 + s + beta);
            pSurf.refGrad()[coli] = 0;
        }
    }

    // Non-orthogonal pressure corrector loop
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        // Breaking columns: couple the cells of a column through a large
        // extra mobility on the layer faces, so that p_rgh is uniform in the
        // vertical (hydrostatic). Only the difference between the two cell
        // values is used: scaling rAUf itself would also scale the explicit
        // non-orthogonal correction, i.e. the horizontal pressure gradient
        // on sloping layer faces, and drive spurious vertical fluxes.
        if (anyBreaking)
        {
            scalarField& upper = p_rghEqn.upper();
            scalarField& diag = p_rghEqn.diag();
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            const scalarField& magSf = mesh.magSf().primitiveField();
            const scalarField& dc = mesh.nonOrthDeltaCoeffs().primitiveField();
            const scalarField& rAUfi = rAUf.primitiveField();
            const labelListList& layerFaces = columns.columnLayerFaces();

            forAll(colBreaking, coli)
            {
                if (colBreaking[coli])
                {
                    const labelList& lf = layerFaces[coli];
                    forAll(lf, i)
                    {
                        const label facei = lf[i];
                        const scalar c =
                            (hfaFactor - 1)*rAUfi[facei]*magSf[facei]*dc[facei];

                        upper[facei] += c;
                        diag[own[facei]] -= c;
                        diag[nei[facei]] -= c;
                    }
                }
            }
        }

        p_rghEqn.solve
        (
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();
        }
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

    // Explicitly relax pressure for the momentum corrector
    p_rgh.relax();

    U = HbyA - rAU*fvc::grad(p_rgh);
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
