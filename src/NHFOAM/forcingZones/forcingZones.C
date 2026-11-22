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

#include "forcingZones.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forcingZone::forcingZone
(
    const word& name,
    const dictionary& dict,
    const scalar g,
    const scalar swl
)
:
    name_(name),
    isWave_(false),
    relaxEta_(dict.getOrDefault<bool>("relaxEta", true)),
    start_(dict.get<point>("start")),
    dir_(Zero),
    length_(0),
    sigmaMax_(dict.get<scalar>("sigmaMax")),
    exponent_(dict.getOrDefault<scalar>("exponent", 2)),
    swl_(swl),
    wave_(nullptr)
{
    const word type(dict.get<word>("type"));

    if (type == "wave")
    {
        isWave_ = true;
        wave_.reset(new regularWave(dict, g, swl));
    }
    else if (type != "sponge")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown forcing zone type " << type
            << " (zone " << name << "). Valid types: sponge wave"
            << exit(FatalIOError);
    }

    const point end(dict.get<point>("end"));

    dir_ = end - start_;
    dir_.z() = 0;
    length_ = mag(dir_);

    if (length_ < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Zone " << name << ": start and end must be separated"
            << " horizontally" << exit(FatalIOError);
    }

    dir_ /= length_;

    Info<< "    forcing zone " << name_ << " (" << type << "): from "
        << start_ << " to " << end << ", sigmaMax " << sigmaMax_;
    if (!relaxEta_)
    {
        Info<< ", velocity only";
    }
    Info<< endl;
}


Foam::forcingZones::forcingZones(const dictionary& dict, const scalar g)
:
    zones_()
{
    if (!dict.found("forcingZones"))
    {
        return;
    }

    const scalar swl = dict.getOrDefault<scalar>("stillWaterLevel", 0);
    const dictionary& zonesDict = dict.subDict("forcingZones");
    const wordList names(zonesDict.toc());

    label n = 0;
    forAll(names, i)
    {
        if (zonesDict.isDict(names[i]))
        {
            ++n;
        }
    }

    zones_.setSize(n);

    n = 0;
    forAll(names, i)
    {
        if (zonesDict.isDict(names[i]))
        {
            zones_.set
            (
                n++,
                new forcingZone(names[i], zonesDict.subDict(names[i]), g, swl)
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::forcingZone::sigma(const point& x) const
{
    vector d = x - start_;
    d.z() = 0;

    const scalar chi = (d & dir_)/length_;

    if (chi < 0 || chi > 1)
    {
        return 0;
    }

    return sigmaMax_*Foam::pow(chi, exponent_);
}


Foam::scalar Foam::forcingZone::etaTarget
(
    const point& x,
    const scalar t
) const
{
    if (isWave_)
    {
        return swl_ + wave_->eta(x, t);
    }

    return swl_;
}


Foam::vector Foam::forcingZone::UTarget
(
    const point& x,
    const scalar t
) const
{
    if (isWave_)
    {
        return wave_->velocity(x, t);
    }

    return Zero;
}


void Foam::forcingZones::cellTargets
(
    const vectorField& C,
    const scalar t,
    scalarField& sigma,
    vectorField& sigmaU
) const
{
    sigma = 0;
    sigmaU = Zero;

    forAll(zones_, zonei)
    {
        const forcingZone& zone = zones_[zonei];

        forAll(C, celli)
        {
            const scalar s = zone.sigma(C[celli]);

            if (s > 0)
            {
                sigma[celli] += s;
                sigmaU[celli] += s*zone.UTarget(C[celli], t);
            }
        }
    }
}


void Foam::forcingZones::columnTargets
(
    const vectorField& colCentre,
    const scalar t,
    scalarField& sigma,
    scalarField& sigmaEta
) const
{
    sigma = 0;
    sigmaEta = 0;

    forAll(zones_, zonei)
    {
        const forcingZone& zone = zones_[zonei];

        if (!zone.relaxEta())
        {
            continue;
        }

        forAll(colCentre, coli)
        {
            const scalar s = zone.sigma(colCentre[coli]);

            if (s > 0)
            {
                sigma[coli] += s;
                sigmaEta[coli] += s*zone.etaTarget(colCentre[coli], t);
            }
        }
    }
}


// ************************************************************************* //
