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

\*---------------------------------------------------------------------------*/

#include "nhMovingSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nhMovingSlipVelocityFvPatchVectorField::
nhMovingSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::nhMovingSlipVelocityFvPatchVectorField::
nhMovingSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::nhMovingSlipVelocityFvPatchVectorField::
nhMovingSlipVelocityFvPatchVectorField
(
    const nhMovingSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::nhMovingSlipVelocityFvPatchVectorField::
nhMovingSlipVelocityFvPatchVectorField
(
    const nhMovingSlipVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


Foam::nhMovingSlipVelocityFvPatchVectorField::
nhMovingSlipVelocityFvPatchVectorField
(
    const nhMovingSlipVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nhMovingSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();
    const vectorField n(p.nf());
    const vectorField Uc(patchInternalField());

    // Tangential velocity of the adjacent cell (slip)
    vectorField Up(Uc - n*(n & Uc));

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        // Normal velocity of the moving wall from the mesh-motion flux
        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());

        tmp<surfaceScalarField> tmeshPhi(fvc::meshPhi(U));
        const scalarField phip(tmeshPhi().boundaryField()[p.index()]);

        Up += n*phip/(p.magSf() + VSMALL);
    }

    vectorField::operator=(Up);

    fixedValueFvPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        nhMovingSlipVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
