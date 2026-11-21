// Hydrostatic front approximation (Smit, Zijlema & Stelling 2013):
//   onset       d(eta)/dt > alpha sqrt(g D)
//   persistence d(eta)/dt > beta  sqrt(g D), for columns that were breaking
//               or neighbour a breaking column
// d(eta)/dt is taken from the last two time levels.
{
    const boolList wasBreaking(colBreaking);
    colBreaking = false;

    if (nEtaLevels >= 2)
    {
        const scalar dt0 = runTime.deltaT0Value();
        const scalarField& zb = columns.columnBottomZ();
        const labelListList& nbrs = columns.columnNeighbours();

        boolList persist(colBreaking.size(), false);

        forAll(colBreaking, coli)
        {
            const scalar dEtaDt = (etaN[coli] - etaNm1[coli])/dt0;
            const scalar c = Foam::sqrt(gMag*max(etaN[coli] - zb[coli], 0));

            persist[coli] = (dEtaDt > hfaBeta*c);
            colBreaking[coli] =
                (dEtaDt > hfaAlpha*c) || (wasBreaking[coli] && persist[coli]);
        }

        // Spread to neighbouring columns on the wave front
        bool changed = true;
        label sweep = 0;

        while (changed && sweep++ < colBreaking.size())
        {
            changed = false;

            forAll(colBreaking, coli)
            {
                if (!colBreaking[coli] && persist[coli])
                {
                    const labelList& nb = nbrs[coli];
                    forAll(nb, i)
                    {
                        if (colBreaking[nb[i]])
                        {
                            colBreaking[coli] = true;
                            changed = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    label nBreaking = 0;
    boundBox bb;
    forAll(colBreaking, coli)
    {
        if (colBreaking[coli])
        {
            ++nBreaking;
            bb.add(columns.columnCentre()[coli]);
        }
    }

    reduce(nBreaking, sumOp<label>());
    bb.reduce();

    Info<< "Breaking columns: " << nBreaking;
    if (nBreaking)
    {
        Info<< " in " << bb.min() << " - " << bb.max();
    }
    Info<< endl;

    // Breaking indicator in the cells
    scalarField colB(colBreaking.size(), 0);
    forAll(colBreaking, coli)
    {
        if (colBreaking[coli])
        {
            colB[coli] = 1;
        }
    }
    columns.setCellValues(colB, breaking.primitiveFieldRef());
    breaking.correctBoundaryConditions();
}
