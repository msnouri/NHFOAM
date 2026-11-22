// Move the sigma mesh onto the free surface: the extrapolated surface on the
// first outer corrector, the latest iterate on the following ones.
// The face flux is rebuilt on the moved faces from the face velocity and
// made relative to the mesh motion for the momentum equation. On the first
// outer corrector the face velocity is extrapolated from levels n and n-1
// (fluxExtrapolation): convecting with the flux of level n alone makes the
// scheme first order in time.
{
    etaMesh = (pimple.firstIter() ? etaPred : etaCol);

    mesh.movePoints(columns.newPoints(etaMesh, minDepth));

    columns.updateGeometry();

    if (pimple.firstIter() && fluxExtrapolation && bdf2 && nEtaLevels >= 2)
    {
        const scalar r = runTime.deltaTValue()/runTime.deltaT0Value();
        phi = mesh.Sf() & ((1 + r)*UfN - r*UfNm1);
    }
    else
    {
        phi = mesh.Sf() & Uf;
    }

    // Wetting and drying: the face velocity times the upwind (TVD) water
    // depth at the side faces
    if (wetDry.scaleFaceDepth())
    {
        wetDry.faceFactors
        (
            max(etaMesh, wetDry.etaMin()),
            phi,
            faceFactor,
            patchFaceFactor
        );

        phi.primitiveFieldRef() *= faceFactor;

        forAll(phi.boundaryField(), patchi)
        {
            if (phi.boundaryField()[patchi].coupled())
            {
                phi.boundaryFieldRef()[patchi] *= patchFaceFactor[patchi];
            }
        }
    }

    fvc::makeRelative(phi, U);

    // Wetting and drying: no flux through closed faces
    if (wetDry.active())
    {
        const boolList& open = wetDry.faceOpen();
        scalarField& phii = phi.primitiveFieldRef();

        forAll(open, facei)
        {
            if (!open[facei])
            {
                phii[facei] = 0;
            }
        }

        forAll(phi.boundaryField(), patchi)
        {
            if (phi.boundaryField()[patchi].coupled())
            {
                const boolList& openp = wetDry.patchFaceOpen()[patchi];
                scalarField& phip = phi.boundaryFieldRef()[patchi];

                forAll(openp, i)
                {
                    if (!openp[i])
                    {
                        phip[i] = 0;
                    }
                }
            }
        }
    }
}
