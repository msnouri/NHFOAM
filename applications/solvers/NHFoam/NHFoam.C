/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | NHFoam: non-hydrostatic free-surface wave solver
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of NHFoam, built on OpenFOAM.

    NHFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Application
    NHFoam

Group
    grpIncompressibleSolvers

Description
    Non-hydrostatic free-surface wave solver on a surface- and
    terrain-following (sigma) moving mesh, in the spirit of NHWAVE and SWASH,
    for unstructured horizontal meshes.

    Model
      - Water only.  The free surface z = eta(x, y, t) is single valued and
        is the top boundary (patch 'surface') of the mesh.
      - The mesh is a vertical extrusion of a 2D mesh (any cell shape) into
        columns.  Points move vertically only and keep the sigma
        distribution of the initial mesh: z = zb + sigma*(eta - zb).
      - Incompressible Navier-Stokes (laminar, RANS or inviscid) in
        physical coordinates on the moving mesh (arbitrary
        Lagrangian-Eulerian), with the kinematic pressure written as
            p_rgh = p/rho + g z = g eta + q
        (hydrostatic part g*eta, non-hydrostatic part q).  Gravity does not
        appear explicitly and still water stays at rest on any mesh.
      - Semi-implicit free surface: the kinematic condition
            d(eta)/dt = phi_surface/A_h
        is inserted in the surface condition p_rgh = g eta^{n+1}, giving a
        Robin (mixed) condition for p_rgh on the surface patch.  Each
        pressure solve therefore yields q and eta^{n+1} together.  No
        gravity-wave Courant limit; Euler or backward (BDF2) in time,
        following ddt(U).
      - Sponge and wave-making relaxation zones, wave gauges.

    Not (yet) included: wetting and drying, wave breaking models.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "mixedFvPatchFields.H"
#include "OFstream.H"
#include "sigmaColumns.H"
#include "forcingZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Non-hydrostatic free-surface wave solver on a sigma"
        " (surface- and terrain-following) moving mesh"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createPimpleControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createGauges.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readNHControls.H"
        #include "CourantNo.H"
        #include "surfaceCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "startTimeStep.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (moveMesh && (pimple.firstIter() || moveMeshOuterCorrectors))
            {
                #include "moveSigmaMesh.H"
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        #include "endTimeStep.H"

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
