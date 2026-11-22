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
}
