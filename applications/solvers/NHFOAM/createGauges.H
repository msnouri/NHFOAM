// Wave gauges: free-surface elevation of the column nearest to each point,
// written every time step to postProcessing/gauges/<startTime>/eta.dat
//
//     gauges
//     {
//         G1  (5.7 0 0);
//         G2  (10.5 0 0);
//     }

wordList gaugeNames;
List<point> gaugePoints;

if (nhProperties.found("gauges"))
{
    const dictionary& gaugeDict = nhProperties.subDict("gauges");

    gaugeNames = gaugeDict.toc();
    gaugePoints.setSize(gaugeNames.size());

    forAll(gaugeNames, gaugei)
    {
        gaugePoints[gaugei] = gaugeDict.get<point>(gaugeNames[gaugei]);
    }
}

// Local column of each gauge on its owning processor, -1 elsewhere
labelList gaugeColumn(gaugePoints.size(), -1);

forAll(gaugePoints, gaugei)
{
    scalar dist = GREAT;
    const label coli = columns.findNearestColumn(gaugePoints[gaugei], dist);
    const scalar minDist = returnReduce(dist, minOp<scalar>());

    // Ties (within round-off) go to the lowest processor
    label owner =
    (
        (coli >= 0 && dist <= minDist*(1 + 1e-9) + 1e-14)
      ? Pstream::myProcNo()
      : labelMax
    );
    reduce(owner, minOp<label>());

    if (owner == Pstream::myProcNo())
    {
        gaugeColumn[gaugei] = coli;
    }
}

autoPtr<OFstream> gaugeFilePtr;

if (gaugePoints.size() && Pstream::master())
{
    const fileName caseDir
    (
        Pstream::parRun() ? runTime.path()/".." : runTime.path()
    );
    const fileName gaugeDir
    (
        caseDir/"postProcessing"/"gauges"/runTime.timeName()
    );

    mkDir(gaugeDir);

    gaugeFilePtr.reset(new OFstream(gaugeDir/"eta.dat"));

    OFstream& os = gaugeFilePtr();
    os.precision(10);

    os  << "# Free-surface elevation eta [m] at the nearest column" << nl;
    forAll(gaugeNames, gaugei)
    {
        os  << "# " << gaugeNames[gaugei] << " " << gaugePoints[gaugei] << nl;
    }
    os  << "# Time";
    forAll(gaugeNames, gaugei)
    {
        os  << tab << gaugeNames[gaugei];
    }
    os  << endl;
}


// Free-surface snapshots: x, y, bed and eta of every column at the end of the
// time step nearest to each listed time, written to
// postProcessing/surface/<time>.dat by the master (the free surface of
// NHWAVE's field output, without writing the full fields)
//
//     surfaceSnapshots (15.3 15.8 16.3);

scalarList snapshotTimes
(
    nhProperties.getOrDefault<scalarList>("surfaceSnapshots", scalarList())
);
Foam::sort(snapshotTimes);

label nextSnapshot = 0;
while
(
    nextSnapshot < snapshotTimes.size()
 && snapshotTimes[nextSnapshot] < runTime.value()
)
{
    ++nextSnapshot;
}

fileName snapshotDir;

if (snapshotTimes.size())
{
    snapshotDir =
    (
        Pstream::parRun() ? runTime.path()/".." : runTime.path()
    )/"postProcessing"/"surface";

    if (Pstream::master())
    {
        mkDir(snapshotDir);
    }

    Info<< "Free-surface snapshots at " << snapshotTimes.size() - nextSnapshot
        << " times in postProcessing/surface" << nl << endl;
}
