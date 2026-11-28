// Wave breaking: the breaking columns, their breaking age and the breaking
// eddy viscosity.
//   onset           d(eta)/dt > alpha sqrt(g D)                 (etaT)
//                   u_s.n > Bkin C                              (kinematic)
//                   u_s the horizontal velocity of the surface cell, C the
//                   crest speed and n the direction of propagation, both
//                   from the front faces around the column (C = d(eta)/dt
//                   over |grad(eta)|), or C = sqrt(g D) if there are none
//   persistence     for columns that were breaking: d(eta)/dt > eta_t*
//                   (persistence threshold), or d(eta)/dt > 0, i.e. until
//                   the crest has passed (persistence crest; SWASH), or
//                   through the whole wave (persistence wave): over the
//                   crest, down the back and under the next front, until a
//                   crest passes whose front never rose faster than
//                   beta sqrt(g D); a column next to a breaking column
//                   joins above beta sqrt(g D)
//                     eta_t* = beta sqrt(g D)                   (default)
//                     eta_t* = [beta + (betaFinal - beta) min(age/T*, 1)]
//                              sqrt(g D),  T* = Tstar sqrt(D/g)
//                   (with betaFinal; Kennedy et al. 2000 start at the onset
//                   threshold, beta = alpha).  A column next to an older
//                   breaker takes over its age, so the age is that of the
//                   breaking event and travels with the wave.
//   eddy viscosity  nu_b = B delta^2 D d(eta)/dt,
//                   B = min(max(d(eta)/dt/eta_t* - 1, 0), 1)   (kennedy)
//                   nu_b = nufac D^2 |S_h| in the breaking columns, S_h the
//                   horizontal strain rate (depth; SWASH's BREAK command)
// d(eta)/dt is taken from the last two time levels.  The same rules with
// the defaults (etaT, fixed beta) are the hydrostatic front approximation of
// Smit, Zijlema & Stelling (2013).
{
    const boolList wasBreaking(colBreaking);
    const labelList wasPhase(colPhase);
    colPhase = 0;
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
            // Horizontal velocity of the surface cells
            const scalarField ux(columns.columnValues(U.primitiveField().component(0)));
            const scalarField uy(columns.columnValues(U.primitiveField().component(1)));

            // Surface slope of each column
            const volVectorField gradEta(fvc::grad(eta));
            const scalarField gx
            (
                columns.columnValues(gradEta.primitiveField().component(0))
            );
            const scalarField gy
            (
                columns.columnValues(gradEta.primitiveField().component(1))
            );

            // Crest speed and direction from the wave faces around each
            // column (the column and two rings of neighbours): for a wave of
            // permanent form d(eta)/dt = C |grad(eta)| on the front face, so
            //     C = sum(d(eta)/dt) / sum(|grad(eta)|),
            //     n = -sum(grad(eta)) / |sum(grad(eta))|
            // over the front-face columns (d(eta)/dt > 0, slope above
            // crestSlope).  Without such columns C = sqrt(g D).
            forAll(onset, coli)
            {
                labelHashSet ring(16);
                ring.insert(coli);
                for (const label n1 : nbrs[coli])
                {
                    ring.insert(n1);
                    for (const label n2 : nbrs[n1])
                    {
                        ring.insert(n2);
                    }
                }

                scalar sumEtaT = 0, sumSlope = 0;
                vector dir(Zero);

                for (const label cj : ring)
                {
                    const scalar slope = Foam::sqrt(sqr(gx[cj]) + sqr(gy[cj]));

                    if (etaT[cj] > 0 && slope > brkCrestSlope)
                    {
                        sumEtaT += etaT[cj];
                        sumSlope += slope;
                        dir -= vector(gx[cj], gy[cj], 0);
                    }
                }

                scalar C = cc[coli];
                scalar un = Foam::sqrt(sqr(ux[coli]) + sqr(uy[coli]));

                if (sumSlope > VSMALL && mag(dir) > VSMALL)
                {
                    C = sumEtaT/sumSlope;
                    dir /= mag(dir);
                    un = ux[coli]*dir.x() + uy[coli]*dir.y();
                }

                onset[coli] = (un > brkBkin*C);
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
                bool keep;

                if (brkWholeWave)
                {
                    // The column stays breaking over the crest, down the
                    // back of the wave and under the next front, and is
                    // released at the next crest if that front never rose
                    // faster than beta sqrt(g D) (the wave has reformed).
                    // Phase: 1 rising, no steep front yet; 2 rising with a
                    // steep front; 3 behind the crest.
                    const bool steep = etaT[coli] > hfaBeta*cc[coli];

                    // behind the crest, a rise slower than brkRise times
                    // the joining threshold is not a new front
                    if
                    (
                        etaT[coli]
                      > (wasPhase[coli] == 3 ? brkRise*hfaBeta*cc[coli] : 0)
                    )
                    {
                        keep = true;
                        colPhase[coli] =
                        (
                            (steep || wasPhase[coli] == 2) ? 2 : 1
                        );
                    }
                    else
                    {
                        keep = (wasPhase[coli] != 1);
                        colPhase[coli] = (keep ? 3 : 0);
                    }
                }
                else
                {
                    keep =
                        etaT[coli]
                      > (brkToCrest ? scalar(0) : threshold(coli, a));
                }

                if (onset[coli] || keep)
                {
                    colBreaking[coli] = true;
                    colAge[coli] = a;

                    if (onset[coli])
                    {
                        colPhase[coli] = 2;
                    }
                }
            }
            else if (onset[coli])
            {
                colBreaking[coli] = true;
                colPhase[coli] = 2;
            }
        }

        // Spread to neighbouring columns on the wave front.  A column joins
        // a breaking neighbour above the threshold of a new breaker,
        // beta sqrt(g D): with the lower threshold of an old breaker the
        // breaking region would creep along the crests of oblique waves
        // into deep water.  A breaking or joining column takes over the age
        // of an older breaking neighbour, so the age is that of the breaking
        // event and travels with the wave.
        // Join (or take the age of) a breaking neighbour of age aNb
        // With persistence wave the breaking columns outlast the front, so
        // a new column joins only ahead of a steep front: next to a column
        // of phase 2 with a higher surface (frontNb, etaNb).  Otherwise the
        // breaking region would grow seawards, to where the fronts reach
        // beta instead of alpha.
        auto join = [&]
        (
            const label coli,
            const scalar aNb,
            const bool frontNb,
            const scalar etaNb
        )
        {
            if (aNb < 0 || (wet.size() && !wet[coli]))
            {
                return false;
            }

            const scalar a = (colBreaking[coli] ? colAge[coli] : scalar(-1));

            if
            (
                aNb > a
             && (
                    colBreaking[coli]
                 || (
                        etaT[coli] > threshold(coli, 0)
                     && (
                            !brkWholeWave
                         || (frontNb && etaN[coli] < etaNb)
                        )
                    )
                )
            )
            {
                if (!colBreaking[coli])
                {
                    colPhase[coli] = 2;     // joins on a steep front
                }
                colBreaking[coli] = true;
                colAge[coli] = aNb;
                return true;
            }

            return false;
        };

        // Sweeps over the local columns until nothing changes
        auto spreadLocally = [&]()
        {
            bool changed = true;
            label sweep = 0;

            while (changed && sweep++ < nCol)
            {
                changed = false;

                forAll(colBreaking, coli)
                {
                    const labelList& nb = nbrs[coli];

                    forAll(nb, i)
                    {
                        const label cn = nb[i];

                        if
                        (
                            colBreaking[cn]
                         && join(coli, colAge[cn], colPhase[cn] == 2, etaN[cn])
                        )
                        {
                            changed = true;
                        }
                    }
                }
            }
        };

        spreadLocally();

        // Across processor boundaries: exchange the age of the breaking
        // columns (-1: not breaking) through the coupled patches, join, and
        // repeat the local sweeps until nothing changes anywhere
        if (Pstream::parRun())
        {
            const labelList& cellCol = columns.cellColumn();

            for (label iter = 0; iter < 1000; ++iter)
            {
                // the surface of the neighbours, then the phase and the age
                scalarField colA(etaN);
                columns.setCellValues(colA, brkExchange.primitiveFieldRef());
                brkExchange.correctBoundaryConditions();
                PtrList<scalarField> nbEta(mesh.boundary().size());
                PtrList<scalarField> nbPhase(mesh.boundary().size());
                forAll(mesh.boundary(), patchi)
                {
                    if (mesh.boundary()[patchi].coupled())
                    {
                        nbEta.set
                        (
                            patchi,
                            brkExchange.boundaryField()[patchi]
                           .patchNeighbourField()
                        );
                    }
                }

                forAll(colA, coli)
                {
                    colA[coli] = colPhase[coli];
                }
                columns.setCellValues(colA, brkExchange.primitiveFieldRef());
                brkExchange.correctBoundaryConditions();
                forAll(mesh.boundary(), patchi)
                {
                    if (mesh.boundary()[patchi].coupled())
                    {
                        nbPhase.set
                        (
                            patchi,
                            brkExchange.boundaryField()[patchi]
                           .patchNeighbourField()
                        );
                    }
                }

                colA = -1;
                forAll(colBreaking, coli)
                {
                    if (colBreaking[coli])
                    {
                        colA[coli] = colAge[coli];
                    }
                }
                columns.setCellValues(colA, brkExchange.primitiveFieldRef());
                brkExchange.correctBoundaryConditions();

                bool changed = false;

                forAll(mesh.boundary(), patchi)
                {
                    const fvPatch& fvp = mesh.boundary()[patchi];

                    if (!fvp.coupled())
                    {
                        continue;
                    }

                    const scalarField nbA
                    (
                        brkExchange.boundaryField()[patchi]
                       .patchNeighbourField()
                    );
                    const labelUList& fc = fvp.faceCells();

                    forAll(fc, i)
                    {
                        if
                        (
                            join
                            (
                                cellCol[fc[i]],
                                nbA[i],
                                nbPhase[patchi][i] > 1.5
                             && nbPhase[patchi][i] < 2.5,
                                nbEta[patchi][i]
                            )
                        )
                        {
                            changed = true;
                        }
                    }
                }

                if (changed)
                {
                    spreadLocally();
                }

                if (!returnReduce(changed, orOp<bool>()))
                {
                    break;
                }
            }
        }

        // Breaking eddy viscosity: the mixing length of SWASH's BREAK
        // command, nu_b = nufac D^2 |S_h| with the horizontal strain rate of
        // the column, S_h^2 = c (u_x^2 + v_y^2) + (u_y + v_x)^2, c = 2 (c = 1
        // in two dimensions, where SWASH has |u_x|)
        if (breakingVisc && brkDepthVisc)
        {
            const volTensorField gradU(fvc::grad(U));
            const scalarField& Vc = mesh.V();
            const labelListList& colCells = columns.columnCells();
            const scalarField& colA = columns.columnArea();
            const scalar cN = (mesh.nSolutionD() == 3 ? 2 : 1);
            const scalar dt = runTime.deltaTValue();

            forAll(colNuB, coli)
            {
                if (!colBreaking[coli])
                {
                    continue;
                }

                scalar s = 0, v = 0;
                for (const label celli : colCells[coli])
                {
                    const tensor& g = gradU[celli];
                    s +=
                        Vc[celli]
                       *Foam::sqrt
                        (
                            cN*(sqr(g.xx()) + sqr(g.yy()))
                          + sqr(g.xy() + g.yx())
                        );
                    v += Vc[celli];
                }

                colNuB[coli] = brkNuFac*sqr(Dc[coli])*s/max(v, VSMALL);

                if (brkNuLimit)
                {
                    // square columns: dx^2 = dy^2 = column area
                    colNuB[coli] = min
                    (
                        colNuB[coli],
                        colA[coli]/(2*dt*(mesh.nSolutionD() - 1))
                    );
                }
            }
        }

        // Breaking eddy viscosity (Kennedy et al. 2000)
        if (breakingVisc && !brkDepthVisc)
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
            colB[coli] = (brkWholeWave ? max(colPhase[coli], 1) : 1);
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
