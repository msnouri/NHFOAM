// Free-surface and derived fields at the end of the time step
{
    columns.setCellValues(etaCol, eta.primitiveFieldRef());
    eta.correctBoundaryConditions();

    q = p_rgh - magg*eta;
    p = p_rgh + (g & mesh.C());

    const scalarField depth(etaCol - columns.columnBottomZ());
    const scalar volume = gSum(columns.columnArea()*depth);

    Info<< "Free surface: eta min/max = " << gMin(etaCol) << " / "
        << gMax(etaCol) << ", min depth = " << gMin(depth)
        << ", water volume change = " << volume - volume0 << endl;

    if (!wetDry.active() && gMin(depth) < minDepth)
    {
        WarningInFunction
            << "Water depth below minDepth (" << minDepth << "): wetting and"
            << " drying is off (NHFOAMProperties/wettingDrying)" << endl;
    }

    // Time averages
    if (avgDictPtr && runTime.value() > avgStart + 0.5*runTime.deltaTValue())
    {
        const scalar dt = runTime.deltaTValue();
        const scalar w = dt/(avgTime + dt);

        UMeanPtr->primitiveFieldRef() +=
            w*(U.primitiveField() - UMeanPtr->primitiveField());
        etaMeanPtr->primitiveFieldRef() +=
            w*(eta.primitiveField() - etaMeanPtr->primitiveField());

        forAll(avgScalarNames, i)
        {
            const volScalarField& f =
                mesh.lookupObject<volScalarField>(avgScalarNames[i]);

            avgScalarMeans[i].primitiveFieldRef() +=
                w*(f.primitiveField() - avgScalarMeans[i].primitiveField());
        }

        avgTime += dt;
    }

    // Restart state (createFields.H)
    if (runTime.writeTime())
    {
        columns.setCellValues(etaN, etaOld.primitiveFieldRef());
        etaOld.correctBoundaryConditions();

        if (breakingOn)
        {
            columns.setCellValues(colAge, breakingAge.primitiveFieldRef());
            breakingAge.correctBoundaryConditions();
        }

        if (interfaceP)
        {
            columns.setCellValues(psCol, pSurf.primitiveFieldRef());
            pSurf.correctBoundaryConditions();
        }

        nhState.set("nEtaLevels", nEtaLevels);
        nhState.set("avgTime", avgTime);
        nhState.instance() = runTime.timeName();
        nhState.regIOobject::write();
    }

    // Wetting and drying: limiter and shoreline record
    if (wetDry.active())
    {
        const boolList& wet = wetDry.colWet();
        const scalarField& zb = columns.columnBottomZ();
        label nWet = 0;
        scalar maxWetBed = -GREAT;
        scalar maxEtaLand = -GREAT;

        forAll(wet, coli)
        {
            if (wet[coli])
            {
                ++nWet;
                maxWetBed = max(maxWetBed, zb[coli]);

                if (zbInitial[coli] > wdStillWaterLevel)
                {
                    maxEtaLand = max(maxEtaLand, etaCol[coli]);
                }
            }
        }
        reduce(nWet, sumOp<label>());
        reduce(maxWetBed, maxOp<scalar>());
        reduce(maxEtaLand, maxOp<scalar>());

        Info<< "Wetting and drying: limited faces "
            << pIntPtr->nLimited() << ", highest wet bed " << maxWetBed
            << endl;

        if (Pstream::master())
        {
            shorelineFilePtr()
                << runTime.value() << tab << nWet << tab << maxWetBed
                << tab << (maxEtaLand > -GREAT ? maxEtaLand : maxWetBed)
                << endl;
        }
    }

    // Highest free surface while wet (createFields.H)
    if (etaMaxOn)
    {
        if (wetDry.active())
        {
            const boolList& wet = wetDry.colWet();
            forAll(etaMaxCol, coli)
            {
                if (wet[coli])
                {
                    etaMaxCol[coli] = max(etaMaxCol[coli], etaCol[coli]);
                }
            }
        }
        else
        {
            etaMaxCol = max(etaMaxCol, etaCol);
        }

        if (runTime.writeTime())
        {
            columns.setCellValues(etaMaxCol, etaMaxPtr->primitiveFieldRef());
            etaMaxPtr->correctBoundaryConditions();
        }
    }

    // Wave gauges
    if (gaugePoints.size())
    {
        scalarField gaugeEta(gaugePoints.size(), 0);

        forAll(gaugeColumn, gaugei)
        {
            if (gaugeColumn[gaugei] >= 0)
            {
                gaugeEta[gaugei] = etaCol[gaugeColumn[gaugei]];
            }
            reduce(gaugeEta[gaugei], sumOp<scalar>());
        }

        if (Pstream::master())
        {
            OFstream& os = gaugeFilePtr();

            os  << runTime.value();
            forAll(gaugeEta, gaugei)
            {
                os  << tab << gaugeEta[gaugei];
            }
            os  << endl;
        }
    }

    // Free-surface snapshots (createGauges.H)
    if
    (
        nextSnapshot < snapshotTimes.size()
     && runTime.value()
      > snapshotTimes[nextSnapshot] - 0.5*runTime.deltaTValue()
    )
    {
        while
        (
            nextSnapshot < snapshotTimes.size()
         && runTime.value()
          > snapshotTimes[nextSnapshot] - 0.5*runTime.deltaTValue()
        )
        {
            ++nextSnapshot;
        }

        List<vectorField> allC(Pstream::nProcs());
        List<scalarField> allZb(Pstream::nProcs());
        List<scalarField> allEta(Pstream::nProcs());
        allC[Pstream::myProcNo()] = columns.columnCentre();
        allZb[Pstream::myProcNo()] = columns.columnBottomZ();
        allEta[Pstream::myProcNo()] = etaCol;
        Pstream::gatherList(allC);
        Pstream::gatherList(allZb);
        Pstream::gatherList(allEta);

        if (Pstream::master())
        {
            OFstream os(snapshotDir/(runTime.timeName() + ".dat"));
            os.precision(8);
            os  << "# Free surface of the columns at t = " << runTime.value()
                << " s" << nl << "# x y zb eta" << nl;
            forAll(allC, proci)
            {
                forAll(allC[proci], coli)
                {
                    os  << allC[proci][coli].x() << ' '
                        << allC[proci][coli].y() << ' '
                        << allZb[proci][coli] << ' '
                        << allEta[proci][coli] << nl;
                }
            }
        }
        Info<< "Free-surface snapshot at t = " << runTime.value() << endl;
    }
}
