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
    Foam::nhWaveVelocityFvPatchVectorField

Description
    Boundary wavemaker: the velocity of a target wave (see regularWave.H:
    linear, irregular or stream-function theory, with the return current)
    at the face centres of a side patch, which follow the free surface of
    the sigma mesh.  The face flux U.Sf drives the free surface of the
    boundary columns, so the patch acts as a wave paddle whose motion is
    the wave's own velocity profile.

    With activeAbsorption the waves that reach the patch from the inside
    are let out by the shallow-water correction (Schaffer and Klopman
    2000), uniform over the depth,
        U = U_target + (c_abs/h) (eta - swl - eta_target) n,
    with n the outward normal, eta the surface elevation of the boundary
    columns and eta_target the target at their centres, both at the old
    time level.  The celerity c_abs is sqrt(g h) (absorptionCelerity
    shallow, the default), the phase speed omega/k of the target wave
    (absorptionCelerity wave: the depth-averaged velocity of a linear wave
    of that period is matched exactly) or a given value.  Only the
    deviation from the target is let out, so a patch with zero wave height
    (waveHeight 0) is an absorbing boundary.

    The pressure (p_rgh) on the patch is fixedFluxPressure, as on a wall.

    Usage
    \verbatim
    inlet
    {
        type             nhWaveVelocity;
        activeAbsorption yes;       // optional (default yes)
        absorptionCelerity shallow; // shallow | wave | <value> [m/s]
        stillWaterLevel  0;         // optional (default: NHFOAMProperties)
        eta              eta;       // optional: name of the eta field
        wave
        {
            theory      linear;
            waveHeight  0.02;
            wavePeriod  2.02;
            depth       0.4;
            direction   (1 0 0);
            rampTime    4;
        }
        value            uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    nhWaveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef nhWaveVelocityFvPatchVectorField_H
#define nhWaveVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "regularWave.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class nhWaveVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class nhWaveVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Target wave entries
        dictionary waveDict_;

        //- Let out the waves reaching the patch from the inside
        bool activeAbsorption_;

        //- Celerity of the absorption: shallow, wave or value
        word absorptionMode_;
        scalar absorptionCelerity_;

        //- Still-water level (from NHFOAMProperties unless given)
        bool swlGiven_;
        mutable scalar swl_;

        //- Name of the surface elevation field
        word etaName_;

        //- Target wave, constructed on first use
        mutable autoPtr<regularWave> wavePtr_;
        mutable scalar gMag_;


    // Private Member Functions

        //- Target wave (constructed on first use, when g is known)
        const regularWave& wave() const;


public:

    //- Runtime type information
    TypeName("nhWaveVelocity");


    // Constructors

        //- Construct from patch and internal field
        nhWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nhWaveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        nhWaveVelocityFvPatchVectorField
        (
            const nhWaveVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        nhWaveVelocityFvPatchVectorField
        (
            const nhWaveVelocityFvPatchVectorField&
        );

        //- Copy construct setting internal field reference
        nhWaveVelocityFvPatchVectorField
        (
            const nhWaveVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new nhWaveVelocityFvPatchVectorField(*this)
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
                new nhWaveVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
