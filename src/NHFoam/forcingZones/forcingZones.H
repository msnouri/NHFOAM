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

Class
    Foam::forcingZone
    Foam::forcingZones

Description
    Sponge (absorbing) and wave-making relaxation zones.

    Inside a zone the solution is relaxed towards a target state with a rate
    sigma [1/s]:
        dU/dt   = ... + sigma (U_target - U)
        deta/dt = ... + sigma (eta_target - eta)
    Both terms are treated implicitly (the eta term inside the free-surface
    condition).  Targets are the still water (type sponge) or a linear wave
    (type wave), so a wave zone generates waves and absorbs reflected ones.

    The strength grows from 0 at 'start' to sigmaMax at 'end' (horizontal
    distance along end - start):
        sigma = sigmaMax*chi^exponent,  chi in [0, 1]

    Example (in constant/NHFoamProperties)
    \verbatim
    stillWaterLevel 0;

    forcingZones
    {
        inlet
        {
            type        wave;
            start       (4 0 0);
            end         (0 0 0);
            sigmaMax    10;
            exponent    2;
            waveHeight  0.02;
            wavePeriod  2.02;
            depth       0.4;
            direction   (1 0 0);
            rampTime    4;
        }
        outlet
        {
            type        sponge;
            start       (25 0 0);
            end         (35 0 0);
            sigmaMax    10;
        }
    }
    \endverbatim

SourceFiles
    forcingZones.C

\*---------------------------------------------------------------------------*/

#ifndef forcingZones_H
#define forcingZones_H

#include "linearWave.H"
#include "PtrList.H"
#include "autoPtr.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class forcingZone Declaration
\*---------------------------------------------------------------------------*/

class forcingZone
{
    // Private data

        word name_;
        bool isWave_;
        bool relaxEta_;
        point start_;
        vector dir_;
        scalar length_;
        scalar sigmaMax_;
        scalar exponent_;
        scalar swl_;
        autoPtr<linearWave> wave_;


public:

    //- Construct from dictionary
    forcingZone
    (
        const word& name,
        const dictionary& dict,
        const scalar g,
        const scalar swl
    );


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        //- Whether the zone also relaxes the surface elevation (false: it
        //  damps the velocity only and conserves the water volume)
        bool relaxEta() const
        {
            return relaxEta_;
        }

        //- Relaxation rate at a position (0 outside the zone)
        scalar sigma(const point& x) const;

        //- Target surface elevation
        scalar etaTarget(const point& x, const scalar t) const;

        //- Target velocity
        vector UTarget(const point& x, const scalar t) const;
};


/*---------------------------------------------------------------------------*\
                        Class forcingZones Declaration
\*---------------------------------------------------------------------------*/

class forcingZones
{
    // Private data

        PtrList<forcingZone> zones_;


public:

    //- Construct from the NHFoam dictionary (optional sub-dictionary
    //  forcingZones) and gravity magnitude
    forcingZones(const dictionary& dict, const scalar g);


    // Member Functions

        bool active() const
        {
            return zones_.size() > 0;
        }

        //- Relaxation rate and rate*target velocity at cell centres
        void cellTargets
        (
            const vectorField& C,
            const scalar t,
            scalarField& sigma,
            vectorField& sigmaU
        ) const;

        //- Relaxation rate and rate*target elevation at column centres
        void columnTargets
        (
            const vectorField& colCentre,
            const scalar t,
            scalarField& sigma,
            scalarField& sigmaEta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
