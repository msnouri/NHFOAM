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

    if (gMin(depth) < minDepth)
    {
        WarningInFunction
            << "Water depth below minDepth (" << minDepth << "): wetting and"
            << " drying is not supported" << endl;
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
