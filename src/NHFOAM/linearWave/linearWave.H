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
    Foam::linearWave

Description
    Linear (Airy) progressive wave on constant depth: surface elevation and
    velocity, with an optional smooth start-up ramp.

    Dictionary entries
    \verbatim
        waveHeight      0.02;       // H [m]
        wavePeriod      2.02;       // T [s]
        depth           0.4;        // still-water depth h [m]
        direction       (1 0 0);    // propagation direction (horizontal)
        phase           0;          // [rad], optional
        rampTime        4;          // [s], optional (default 0: no ramp)
        origin          (0 0 0);    // phase origin, optional
    \endverbatim

    Elevations are relative to the still-water level swl:
        eta = (H/2) cos(k (x - x0).d - omega t + phase)
    Velocities use z measured from swl (z - swl).

SourceFiles
    linearWave.C

\*---------------------------------------------------------------------------*/

#ifndef linearWave_H
#define linearWave_H

#include "dictionary.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearWave Declaration
\*---------------------------------------------------------------------------*/

class linearWave
{
    // Private data

        scalar H_;
        scalar T_;
        scalar h_;
        scalar phase_;
        scalar rampTime_;
        scalar swl_;
        scalar g_;
        vector direction_;
        point origin_;
        scalar omega_;
        scalar k_;


public:

    // Constructors

        //- Construct from dictionary, gravity magnitude and still-water level
        linearWave(const dictionary& dict, const scalar g, const scalar swl);


    // Static Functions

        //- Wave number from the linear dispersion relation
        static scalar waveNumber
        (
            const scalar omega,
            const scalar h,
            const scalar g
        );


    // Member Functions

        scalar k() const
        {
            return k_;
        }

        scalar omega() const
        {
            return omega_;
        }

        scalar depth() const
        {
            return h_;
        }

        //- Start-up ramp factor (0 to 1)
        scalar ramp(const scalar t) const;

        //- Surface elevation above the still-water level
        scalar eta(const point& x, const scalar t) const;

        //- Velocity at position x (z of x is used for the vertical profile)
        vector velocity(const point& x, const scalar t) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
