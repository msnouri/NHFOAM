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

#include "bottomMotion.H"
#include "mathematicalConstants.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bottomMotion::bottomMotion(const dictionary& nhDict)
:
    type_(NONE),
    min_(Zero),
    max_(Zero),
    edgeWidth_(0),
    amplitude_(0),
    history_("exponential"),
    riseTime_(1),
    centre0_(Zero),
    dir_(1, 0, 0),
    theta_(0),
    T_(0),
    b_(1),
    w_(0),
    eps_(0.717),
    ut_(1),
    a0_(1)
{
    if (!nhDict.found("bottomMotion"))
    {
        return;
    }

    const dictionary& dict = nhDict.subDict("bottomMotion");
    const word type(dict.get<word>("type"));

    if (type == "none")
    {
        type_ = NONE;
    }
    else if (type == "uplift")
    {
        type_ = UPLIFT;
        min_ = dict.get<point>("min");
        max_ = dict.get<point>("max");
        edgeWidth_ = dict.getOrDefault<scalar>("edgeWidth", 0);
        amplitude_ = dict.get<scalar>("amplitude");
        history_ = dict.getOrDefault<word>("history", "exponential");
        riseTime_ = dict.get<scalar>("riseTime");

        if
        (
            history_ != "exponential"
         && history_ != "halfCosine"
         && history_ != "linear"
        )
        {
            FatalIOErrorInFunction(dict)
                << "Unknown history " << history_
                << ". Valid: exponential halfCosine linear"
                << exit(FatalIOError);
        }

        Info<< "Bottom motion: uplift of " << amplitude_ << " m, "
            << history_ << " history, rise time " << riseTime_ << " s"
            << nl << endl;
    }
    else if (type == "slide")
    {
        type_ = SLIDE;
        centre0_ = dict.get<point>("centre");
        centre0_.z() = 0;
        dir_ = dict.getOrDefault<vector>("direction", vector(1, 0, 0));
        dir_.z() = 0;
        dir_ /= mag(dir_);
        theta_ =
            dict.get<scalar>("slopeAngle")*constant::mathematical::pi/180;
        T_ = dict.get<scalar>("thickness");
        b_ = dict.get<scalar>("length");
        w_ = dict.getOrDefault<scalar>("width", 0);
        eps_ = dict.getOrDefault<scalar>("epsilon", 0.717);
        ut_ = dict.get<scalar>("terminalVelocity");
        a0_ = dict.get<scalar>("initialAcceleration");

        Info<< "Bottom motion: solid slide T = " << T_ << " m, b = " << b_
            << " m, w = " << w_ << " m, slope " << dict.get<scalar>("slopeAngle")
            << " deg, ut = " << ut_ << " m/s, a0 = " << a0_ << " m/s2"
            << nl << endl;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown bottomMotion type " << type
            << ". Valid: none uplift slide" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::bottomMotion::upliftHistory(const scalar t) const
{
    if (t <= 0)
    {
        return 0;
    }

    if (history_ == "exponential")
    {
        // Hammack (1973): alpha = 1.11/tc, tc = time of 2/3 of the motion
        return 1 - Foam::exp(-1.11*t/riseTime_);
    }

    if (t >= riseTime_)
    {
        return 1;
    }

    if (history_ == "halfCosine")
    {
        return 0.5*(1 - Foam::cos(constant::mathematical::pi*t/riseTime_));
    }

    return t/riseTime_;
}


Foam::scalar Foam::bottomMotion::upliftShape(const point& x) const
{
    scalar f = 1;

    for (direction cmpt = 0; cmpt < 2; ++cmpt)
    {
        const scalar lo = min_[cmpt];
        const scalar hi = max_[cmpt];
        const scalar xi = x[cmpt];

        if (edgeWidth_ <= 0)
        {
            if (xi < lo || xi > hi)
            {
                return 0;
            }
        }
        else
        {
            // Cosine taper of width edgeWidth inside the box edges
            const scalar d = min(xi - lo, hi - xi);

            if (d <= 0)
            {
                return 0;
            }
            if (d < edgeWidth_)
            {
                f *= 0.5*(1 - Foam::cos(constant::mathematical::pi*d/edgeWidth_));
            }
        }
    }

    return f;
}


Foam::scalar Foam::bottomMotion::slideDistance(const scalar t) const
{
    if (t <= 0)
    {
        return 0;
    }

    const scalar s0 = sqr(ut_)/a0_;
    const scalar t0 = ut_/a0_;

    return s0*Foam::log(Foam::cosh(t/t0));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::bottomMotion::displacement
(
    const point& x,
    const scalar t
) const
{
    if (type_ == UPLIFT)
    {
        return amplitude_*upliftHistory(t)*upliftShape(x);
    }

    if (type_ == SLIDE)
    {
        const scalar cosTheta = Foam::cos(theta_);
        const scalar C = Foam::acosh(1/eps_);
        const scalar kb = 2*C/b_;

        // Horizontal position of the slide centre
        const point centre = centre0_ + slideDistance(t)*cosTheta*dir_;

        vector d = x - centre;
        d.z() = 0;

        // Along-slope and cross-slope coordinates
        const scalar xi = (d & dir_)/cosTheta;
        const scalar chi = mag(d - (d & dir_)*dir_);

        scalar shape = 1/Foam::cosh(kb*xi);

        if (w_ > 0)
        {
            const scalar kw = 2*C/w_;
            shape /= Foam::cosh(kw*chi);
        }

        const scalar thickness = T_*max(shape - eps_, 0)/(1 - eps_);

        // Thickness measured normal to the slope -> vertical
        return thickness/cosTheta;
    }

    return 0;
}


Foam::tmp<Foam::scalarField> Foam::bottomMotion::displacement
(
    const vectorField& x,
    const scalar t
) const
{
    tmp<scalarField> tdz(new scalarField(x.size(), 0));
    scalarField& dz = tdz.ref();

    if (active())
    {
        forAll(x, i)
        {
            dz[i] = displacement(x[i], t);
        }
    }

    return tdz;
}


// ************************************************************************* //
