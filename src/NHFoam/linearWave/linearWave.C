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

\*---------------------------------------------------------------------------*/

#include "linearWave.H"
#include "mathematicalConstants.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearWave::linearWave
(
    const dictionary& dict,
    const scalar g,
    const scalar swl
)
:
    H_(dict.get<scalar>("waveHeight")),
    T_(dict.get<scalar>("wavePeriod")),
    h_(dict.get<scalar>("depth")),
    phase_(dict.getOrDefault<scalar>("phase", 0)),
    rampTime_(dict.getOrDefault<scalar>("rampTime", 0)),
    swl_(swl),
    g_(g),
    direction_(dict.getOrDefault<vector>("direction", vector(1, 0, 0))),
    origin_(dict.getOrDefault<point>("origin", point::zero)),
    omega_(constant::mathematical::twoPi/T_),
    k_(waveNumber(omega_, h_, g_))
{
    direction_.z() = 0;

    if (mag(direction_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Wave direction must have a horizontal component"
            << exit(FatalIOError);
    }

    direction_ /= mag(direction_);
}


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::linearWave::waveNumber
(
    const scalar omega,
    const scalar h,
    const scalar g
)
{
    const scalar omega2 = sqr(omega);

    // Eckart's approximation as the initial guess, then Newton
    scalar k = omega2/(g*Foam::sqrt(Foam::tanh(omega2*h/g)));

    for (label iter = 0; iter < 100; ++iter)
    {
        const scalar th = Foam::tanh(k*h);
        const scalar f = g*k*th - omega2;
        const scalar df = g*th + g*k*h*(1 - sqr(th));
        const scalar dk = f/df;

        k -= dk;

        if (mag(dk) < 1e-14*k)
        {
            break;
        }
    }

    return k;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::linearWave::ramp(const scalar t) const
{
    if (rampTime_ <= 0 || t >= rampTime_)
    {
        return 1;
    }
    if (t <= 0)
    {
        return 0;
    }

    return 0.5*(1 - Foam::cos(constant::mathematical::pi*t/rampTime_));
}


Foam::scalar Foam::linearWave::eta(const point& x, const scalar t) const
{
    const scalar arg = k_*((x - origin_) & direction_) - omega_*t + phase_;

    return ramp(t)*0.5*H_*Foam::cos(arg);
}


Foam::vector Foam::linearWave::velocity(const point& x, const scalar t) const
{
    const scalar arg = k_*((x - origin_) & direction_) - omega_*t + phase_;

    // Height above the bottom, limited to the water column
    const scalar zb = max(x.z() - swl_ + h_, 0);
    const scalar a = 0.5*H_*omega_/Foam::sinh(k_*h_);

    const scalar u = a*Foam::cosh(k_*zb)*Foam::cos(arg);
    const scalar w = a*Foam::sinh(k_*zb)*Foam::sin(arg);

    return ramp(t)*(u*direction_ + vector(0, 0, w));
}


// ************************************************************************* //
