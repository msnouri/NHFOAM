/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | NHFOAM: non-hydrostatic free-surface wave solver
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of NHFOAM, built on OpenFOAM.

    NHFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Description
    Incompressible turbulence models of NHFOAM (RAS kOmegaLF).

\*---------------------------------------------------------------------------*/

#include "turbulentTransportModels.H"

#include "kOmegaLF.H"
makeRASModel(kOmegaLF);

// ************************************************************************* //
