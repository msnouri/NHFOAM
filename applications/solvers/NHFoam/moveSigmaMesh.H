// Move the sigma mesh onto the free surface: the extrapolated surface on the
// first outer corrector, the latest iterate on the following ones.
// The face flux is rebuilt on the moved faces from the face velocity and
// made relative to the mesh motion for the momentum equation.
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

    phi = mesh.Sf() & Uf;
    fvc::makeRelative(phi, U);
}
