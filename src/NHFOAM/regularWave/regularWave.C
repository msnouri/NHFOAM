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

#include "regularWave.H"
#include "mathematicalConstants.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regularWave::regularWave
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
    k_(waveNumber(omega_, h_, g_)),
    streamFunction_(false),
    uMean_(0)
{
    direction_.z() = 0;

    if (mag(direction_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Wave direction must have a horizontal component"
            << exit(FatalIOError);
    }

    direction_ /= mag(direction_);

    const word theory(dict.getOrDefault<word>("theory", "linear"));

    if (theory == "streamFunction")
    {
        streamFunction_ = true;
        k_ = dict.get<scalar>("waveNumber");
        uMean_ = dict.get<scalar>("uMean");
        E_ = dict.get<scalarList>("etaCoeffs");
        B_ = dict.get<scalarList>("psiCoeffs");

        if (E_.size() != B_.size() + 1)
        {
            FatalIOErrorInFunction(dict)
                << "etaCoeffs needs one entry (E_0) more than psiCoeffs"
                << exit(FatalIOError);
        }
    }
    else if (theory != "linear")
    {
        FatalIOErrorInFunction(dict)
            << "theory " << theory << " unknown; use linear or streamFunction"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::regularWave::waveNumber
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


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::scalar Foam::regularWave::theta(const point& x, const scalar t) const
{
    return k_*((x - origin_) & direction_) - omega_*t + phase_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::regularWave::ramp(const scalar t) const
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


Foam::scalar Foam::regularWave::eta(const point& x, const scalar t) const
{
    const scalar th = theta(x, t);

    if (streamFunction_)
    {
        scalar e = 0;
        forAll(E_, j)
        {
            e += E_[j]*Foam::cos(j*th);
        }
        return ramp(t)*e;
    }

    return ramp(t)*0.5*H_*Foam::cos(th);
}


Foam::vector Foam::regularWave::velocity(const point& x, const scalar t) const
{
    const scalar th = theta(x, t);

    // Height above the bottom, limited to the water column
    const scalar zb = max(x.z() - swl_ + h_, 0);

    scalar u = 0, w = 0;

    if (streamFunction_)
    {
        u = uMean_;

        forAll(B_, i)
        {
            const label j = i + 1;
            const scalar jk = j*k_;
            const scalar ch = Foam::cosh(jk*h_);

            u += jk*B_[i]*Foam::cosh(jk*zb)/ch*Foam::cos(j*th);
            w += jk*B_[i]*Foam::sinh(jk*zb)/ch*Foam::sin(j*th);
        }
    }
    else
    {
        const scalar a = 0.5*H_*omega_/Foam::sinh(k_*h_);

        u = a*Foam::cosh(k_*zb)*Foam::cos(th);
        w = a*Foam::sinh(k_*zb)*Foam::sin(th);
    }

    return ramp(t)*(u*direction_ + vector(0, 0, w));
}


// ************************************************************************* //
