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
    Foam::regularWave

Description
    Target wave on constant depth for the wave-making zones, the boundary
    wavemaker (nhWaveVelocity) and the initial conditions (nhSetWaves):
    surface elevation and velocity, with an optional smooth start-up ramp.
    Three theories:

    linear (default): the Airy wave,
        eta = (H/2) cos(theta),  theta = k (x - x0).d - omega t + phase,
    with k from the linear dispersion relation.

    irregular: a sum of linear components (amplitude a_i, frequency
    omega_i, direction d_i, phase phi_i),
        eta = sum_i a_i cos(theta_i),
        u   = sum_i a_i omega_i cosh(k_i z_b)/sinh(k_i h) cos(theta_i) d_i,
        w   = sum_i a_i omega_i sinh(k_i z_b)/sinh(k_i h) sin(theta_i),
    from a JONSWAP spectrum (random phases, optional cos^{2s} directional
    spreading) or from an explicit list of components.  The spectrum is
    scaled so that the discrete components have the requested
    Hm0 = 4 sqrt(m0).

    streamFunction: a steady wave of any height from the Fourier
    approximation of Fenton (1988), whose coefficients are computed
    beforehand (validation/stream_function.py writes these entries):
        eta = sum_{j=0}^{N} E_j cos(j theta)
        u   = uMean + sum_j j k B_j cosh(j k z_b)/cosh(j k h) cos(j theta)
        w   =         sum_j j k B_j sinh(j k z_b)/cosh(j k h) sin(j theta)
    with z_b the height above the bottom and uMean the Eulerian mean current
    (the return current of a closed flume for zero mass flux).

    Return current: a linear wave applied up to the moving free surface
    carries the mass flux M = g a^2/(2 c) per component (the Stokes
    transport between trough and crest).  With returnCurrent (default yes,
    as in the boundary wavemakers of NHWAVE) the depth-uniform current
    -M/h along each component's direction is added, so that the target has
    no net mass flux, as in a closed flume.  It grows with the square of
    the ramp.

    Dictionary entries
    \verbatim
        theory          linear;     // linear | irregular | streamFunction
        waveHeight      0.02;       // H [m]
        wavePeriod      2.02;       // T [s]
        depth           0.4;        // still-water depth h [m]
        direction       (1 0 0);    // propagation direction (horizontal)
        phase           0;          // [rad], optional
        rampTime        4;          // [s], optional (default 0: no ramp)
        origin          (0 0 0);    // phase origin, optional
        returnCurrent   yes;        // linear and irregular, optional

        // irregular: JONSWAP spectrum ...
        Hs              0.05;       // significant wave height Hm0 [m]
        Tp              1.5;        // peak period [s]
        gamma           3.3;        // peak enhancement (1: Pierson-Moskowitz)
        nFrequencies    100;        // equally spaced in [fMin, fMax]
        fMin            0.333;      // [Hz], optional (default 0.5/Tp)
        fMax            2.0;        // [Hz], optional (default 3/Tp)
        seed            1;          // random phases, optional
        nDirections     1;          // optional (default 1: unidirectional)
        spreading       10;         // s of cos^{2s}((theta - theta0)/2)
        maxAngle        60;         // half-width of the directions [deg]
        periodicWidth   2;          // optional: width of a periodic domain;
                                    // the directions are moved to the
                                    // nearest ones that fit it

        // ... or explicit components (instead of the spectrum)
        components
        (
            // amplitude [m]  period [s]  direction [deg]  phase [rad]
            (0.01  1.5  0   0)
            (0.005 1.2  20  1.3)
        );

        // streamFunction only
        waveNumber      1.659;      // k [1/m]
        uMean           -0.024;     // Eulerian mean current [m/s]
        etaCoeffs       (...);      // E_0 ... E_N [m]
        psiCoeffs       (...);      // B_1 ... B_N [m^2/s]
    \endverbatim

    Directions of the components are measured from 'direction',
    anticlockwise seen from above.  Elevations are relative to the
    still-water level swl; velocities use z measured from swl (z - swl).

SourceFiles
    regularWave.C

\*---------------------------------------------------------------------------*/

#ifndef regularWave_H
#define regularWave_H

#include "dictionary.H"
#include "point.H"
#include "scalarField.H"
#include "vectorList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regularWave Declaration
\*---------------------------------------------------------------------------*/

class regularWave
{
    // Private data

        //- Representative height and period (Hs, Tp for irregular waves)
        scalar H_;
        scalar T_;
        scalar h_;
        scalar rampTime_;
        scalar swl_;
        scalar g_;
        vector direction_;
        point origin_;
        scalar omega_;
        scalar k_;

        //- Linear components (one for theory linear)
        scalarList amp_;
        scalarList omegas_;
        scalarList ks_;
        scalarList phases_;
        vectorList dirs_;

        //- Return current of the linear components (at full ramp)
        vector uReturn_;

        //- Stream-function wave (Fenton 1988)
        bool streamFunction_;
        scalar uMean_;
        scalar phase_;
        scalarList E_;
        scalarList B_;


    // Private Member Functions

        //- Phase of the stream-function wave at x and t
        scalar theta(const point& x, const scalar t) const;


public:

    // Constructors

        //- Construct from dictionary, gravity magnitude and still-water level
        regularWave(const dictionary& dict, const scalar g, const scalar swl);


    // Static Functions

        //- Wave number from the linear dispersion relation
        static scalar waveNumber
        (
            const scalar omega,
            const scalar h,
            const scalar g
        );

        //- JONSWAP spectral shape (not normalised) at frequency f [Hz]
        static scalar jonswapShape
        (
            const scalar f,
            const scalar Tp,
            const scalar gamma
        );


    // Member Functions

        //- Wave number of the (peak) wave
        scalar k() const
        {
            return k_;
        }

        //- Angular frequency of the (peak) wave
        scalar omega() const
        {
            return omega_;
        }

        scalar depth() const
        {
            return h_;
        }

        //- Number of linear components (0 for the stream-function wave)
        label nComponents() const
        {
            return amp_.size();
        }

        //- Return current at full ramp
        const vector& returnCurrent() const
        {
            return uReturn_;
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
