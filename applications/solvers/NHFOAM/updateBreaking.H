// Wave breaking: the breaking columns, their breaking age and the breaking
// eddy viscosity.
//   onset           d(eta)/dt > alpha sqrt(g D)                 (etaT)
//                   |u_s| > Bkin sqrt(g D)                      (kinematic)
//                   u_s the horizontal velocity of the surface cell and
//                   sqrt(g D) the long-wave speed at the local total depth
//                   (an estimate of the crest speed)
//   persistence     d(eta)/dt > eta_t*, for columns that were breaking;
//                   a column next to a breaking column joins above
//                   beta sqrt(g D)
//                     eta_t* = beta sqrt(g D)                   (default)
//                     eta_t* = [beta + (betaFinal - beta) min(age/T*, 1)]
//                              sqrt(g D),  T* = Tstar sqrt(D/g)
//                   (with betaFinal; Kennedy et al. 2000 start at the onset
//                   threshold, beta = alpha).  A column next to an older
//                   breaker takes over its age, so the age is that of the
//                   breaking event and travels with the wave.
//   eddy viscosity  nu_b = B delta^2 D d(eta)/dt,
//                   B = min(max(d(eta)/dt/eta_t* - 1, 0), 1)
// d(eta)/dt is taken from the last two time levels.  The same rules with
// the defaults (etaT, fixed beta) are the hydrostatic front approximation of
// Smit, Zijlema & Stelling (2013).
{
    const boolList wasBreaking(colBreaking);
    const scalarField age0(colAge);
    const label nCol = colBreaking.size();

    colBreaking = false;
    colAge = 0;
    colNuB = 0;

    if (nEtaLevels >= 2)
    {
        const scalar dt0 = runTime.deltaT0Value();
        const scalarField& zb = columns.columnBottomZ();
        const labelListList& nbrs = columns.columnNeighbours();
        const boolList noWetMask;
        const boolList& wet = (wetDry.active() ? wetDry.colWet() : noWetMask);

        scalarField etaT(nCol), Dc(nCol), cc(nCol);
        forAll(etaT, coli)
        {
            etaT[coli] = (etaN[coli] - etaNm1[coli])/dt0;
            Dc[coli] = max(etaN[coli] - zb[coli], scalar(0));
            cc[coli] = Foam::sqrt(gMag*Dc[coli]);
        }

        // Persistence threshold eta_t* of a column for breaking age a
        auto threshold = [&](const label coli, const scalar a)
        {
            if (!brkAgeing)
            {
                return hfaBeta*cc[coli];
            }

            const scalar Ts = brkTstar*Foam::sqrt(Dc[coli]/gMag);
            const scalar s = (Ts > VSMALL ? min(a/Ts, scalar(1)) : scalar(1));

            return (hfaBeta + (brkBetaFinal - hfaBeta)*s)*cc[coli];
        };

        // Onset
        boolList onset(nCol, false);

        if (kinematicOnset)
        {
            scalarField uh(mesh.nCells());
            forAll(uh, celli)
            {
                uh[celli] = Foam::sqrt(sqr(U[celli].x()) + sqr(U[celli].y()));
            }
            const scalarField uS(columns.columnValues(uh));

            forAll(onset, coli)
            {
                onset[coli] = (uS[coli] > brkBkin*cc[coli]);
            }
        }
        else
        {
            forAll(onset, coli)
            {
                onset[coli] = (etaT[coli] > hfaAlpha*cc[coli]);
            }
        }

        // Breaking columns and their age
        forAll(colBreaking, coli)
        {
            if (wet.size() && !wet[coli])
            {
                continue;
            }

            if (wasBreaking[coli])
            {
                const scalar a = age0[coli] + dt0;

                if (onset[coli] || etaT[coli] > threshold(coli, a))
                {
                    colBreaking[coli] = true;
                    colAge[coli] = a;
                }
            }
            else if (onset[coli])
            {
                colBreaking[coli] = true;
            }
        }

        // Spread to neighbouring columns on the wave front.  A column joins
        // a breaking neighbour above the threshold of a new breaker,
        // beta sqrt(g D): with the lower threshold of an old breaker the
        // breaking region would creep along the crests of oblique waves
        // into deep water.  A breaking or joining column takes over the age
        // of an older breaking neighbour, so the age is that of the breaking
        // event and travels with the wave.
        bool changed = true;
        label sweep = 0;

        while (changed && sweep++ < nCol)
        {
            changed = false;

            forAll(colBreaking, coli)
            {
                if (wet.size() && !wet[coli])
                {
                    continue;
                }

                const labelList& nb = nbrs[coli];
                scalar a = (colBreaking[coli] ? colAge[coli] : scalar(-1));
                bool join = false;

                forAll(nb, i)
                {
                    const label cn = nb[i];

                    if
                    (
                        colBreaking[cn]
                     && colAge[cn] > a
                     && (colBreaking[coli] || etaT[coli] > threshold(coli, 0))
                    )
                    {
                        a = colAge[cn];
                        join = true;
                    }
                }

                if (join)
                {
                    colBreaking[coli] = true;
                    colAge[coli] = a;
                    changed = true;
                }
            }
        }

        // Breaking eddy viscosity (Kennedy et al. 2000)
        if (breakingVisc)
        {
            forAll(colNuB, coli)
            {
                if (colBreaking[coli] && etaT[coli] > 0)
                {
                    const scalar ts = threshold(coli, colAge[coli]);
                    const scalar B =
                        (
                            ts > VSMALL
                          ? min(max(etaT[coli]/ts - 1, scalar(0)), scalar(1))
                          : scalar(0)
                        );

                    colNuB[coli] = B*sqr(brkDelta)*Dc[coli]*etaT[coli];
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
    if (breakingVisc)
    {
        Info<< ", max nu_b " << gMax(colNuB);
    }
    Info<< endl;

    // Breaking indicator and eddy viscosity in the cells
    scalarField colB(nCol, 0);
    forAll(colBreaking, coli)
    {
        if (colBreaking[coli])
        {
            colB[coli] = 1;
        }
    }
    columns.setCellValues(colB, breaking.primitiveFieldRef());
    breaking.correctBoundaryConditions();

    if (breakingVisc)
    {
        columns.setCellValues(colNuB, nuBreaking.primitiveFieldRef());
        nuBreaking.correctBoundaryConditions();
    }
}
