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

Class
    Foam::nhMovingSlipVelocityFvPatchVectorField

Description
    Slip condition on a moving wall: the tangential velocity is taken from
    the adjacent cell and the normal velocity follows the wall motion, so
    that the face flux equals the mesh-motion flux (kinematic condition of
    a moving, frictionless bottom). On a static mesh it reduces to slip.

    Usage
    \verbatim
    bottom
    {
        type            nhMovingSlipVelocity;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    nhMovingSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef nhMovingSlipVelocityFvPatchVectorField_H
#define nhMovingSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class nhMovingSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class nhMovingSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("nhMovingSlipVelocity");


    // Constructors

        //- Construct from patch and internal field
        nhMovingSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nhMovingSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        nhMovingSlipVelocityFvPatchVectorField
        (
            const nhMovingSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        nhMovingSlipVelocityFvPatchVectorField
        (
            const nhMovingSlipVelocityFvPatchVectorField&
        );

        //- Copy construct setting internal field reference
        nhMovingSlipVelocityFvPatchVectorField
        (
            const nhMovingSlipVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new nhMovingSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new nhMovingSlipVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
