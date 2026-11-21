// Gravity-wave Courant number sqrt(g D) dt/dx (for information only: the
// free surface is integrated implicitly and this does not limit stability)
{
    scalar surfaceCoNum = 0;

    const scalarField& Ah = columns.columnArea();
    const scalarField& zb = columns.columnBottomZ();

    forAll(etaCol, coli)
    {
        const scalar D = max(etaCol[coli] - zb[coli], 0);

        surfaceCoNum = max
        (
            surfaceCoNum,
            Foam::sqrt(gMag*D)*runTime.deltaTValue()/Foam::sqrt(Ah[coli])
        );
    }

    reduce(surfaceCoNum, maxOp<scalar>());

    Info<< "Gravity-wave Courant Number max: " << surfaceCoNum << endl;
}
