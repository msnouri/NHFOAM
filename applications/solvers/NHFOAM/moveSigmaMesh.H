// Move the sigma mesh onto the free surface: the extrapolated surface on the
// first outer corrector, the latest iterate on the following ones.
// The face flux is rebuilt on the moved faces from the face velocity and
// made relative to the mesh motion for the momentum equation. On the first
// outer corrector the face velocity is extrapolated from levels n and n-1
// (fluxExtrapolation): convecting with the flux of level n alone makes the
// scheme first order in time.
{
    mesh.movePoints
    (
        columns.newPoints
        (
            pimple.firstIter() ? etaPred : etaCol,
            minDepth
        )
    );

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
    fvc::makeRelative(phi, U);
}
