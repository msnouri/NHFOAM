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
#include "Random.H"
#include "FixedList.H"
#include "unitConversion.H"
#include "error.H"
#include "TableBase.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regularWave::regularWave
(
    const dictionary& dict,
    const scalar g,
    const scalar swl
)
:
    H_(0),
    T_(0),
    h_(dict.get<scalar>("depth")),
    rampTime_(dict.getOrDefault<scalar>("rampTime", 0)),
    swl_(swl),
    g_(g),
    direction_(dict.getOrDefault<vector>("direction", vector(1, 0, 0))),
    origin_(dict.getOrDefault<point>("origin", point::zero)),
    omega_(0),
    k_(0),
    uReturn_(Zero),
    streamFunction_(false),
    uMean_(0),
    phase_(dict.getOrDefault<scalar>("phase", 0)),
    timeSeries_(false),
    paddle_(false),
    series_(nullptr),
    linearRelation_(false),
    celerity_(0),
    timeShift_(0),
    scale_(1),
    dtDiff_(1e-3),
    verticalVelocity_(true)
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

    // Periodic domain: transverse wave numbers must fit the width
    const scalar periodicWidth
    (
        dict.getOrDefault<scalar>("periodicWidth", 0)
    );
    vector periodicAxis
    (
        dict.getOrDefault<vector>("periodicAxis", vector(0, 1, 0))
    );
    periodicAxis.z() = 0;
    periodicAxis /= max(mag(periodicAxis), SMALL);
    const vector normalAxis(periodicAxis.y(), -periodicAxis.x(), 0);

    DynamicList<scalar> amp, omegas, ks, phases;
    DynamicList<vector> dirs;

    // Linear component at an angle [rad] from the main direction
    auto addComponent = [&]
    (
        const scalar a,
        const scalar omega,
        const scalar angle,
        const scalar phase
    )
    {
        const scalar k = waveNumber(omega, h_, g_);

        const scalar ca = Foam::cos(angle), sa = Foam::sin(angle);
        vector d
        (
            ca*direction_.x() - sa*direction_.y(),
            sa*direction_.x() + ca*direction_.y(),
            0
        );

        if (periodicWidth > 0)
        {
            const scalar dp = d & periodicAxis;
            const scalar dn = d & normalAxis;
            const scalar dk = twoPi/(k*periodicWidth);
            const scalar s =
                max(min(dk*std::round(dp/dk), 1 - SMALL), -1 + SMALL);

            d = s*periodicAxis
              + (dn < 0 ? -1 : 1)*Foam::sqrt(1 - sqr(s))*normalAxis;
        }

        amp.append(a);
        omegas.append(omega);
        ks.append(k);
        phases.append(phase);
        dirs.append(d);
    };

    if (theory == "linear" || theory == "streamFunction")
    {
        H_ = dict.get<scalar>("waveHeight");
        T_ = dict.get<scalar>("wavePeriod");
        omega_ = twoPi/T_;
        k_ = waveNumber(omega_, h_, g_);
    }

    if (theory == "linear")
    {
        addComponent(0.5*H_, omega_, 0, phase_);

        if (periodicWidth > 0)
        {
            Info<< "    linear wave direction for the periodic width "
                << periodicWidth << ": " << dirs[0] << endl;
        }
    }
    else if (theory == "irregular")
    {
        if (dict.found("components"))
        {
            const List<FixedList<scalar, 4>> comps
            (
                dict.get<List<FixedList<scalar, 4>>>("components")
            );

            forAll(comps, i)
            {
                addComponent
                (
                    comps[i][0],
                    twoPi/comps[i][1],
                    degToRad(comps[i][2]),
                    comps[i][3]
                );
            }
        }
        else
        {
            const scalar Hs = dict.get<scalar>("Hs");
            const scalar Tp = dict.get<scalar>("Tp");
            const scalar gamma = dict.getOrDefault<scalar>("gamma", 3.3);
            const label nF = dict.getOrDefault<label>("nFrequencies", 100);
            const scalar fMin = dict.getOrDefault<scalar>("fMin", 0.5/Tp);
            const scalar fMax = dict.getOrDefault<scalar>("fMax", 3/Tp);
            const label seed = dict.getOrDefault<label>("seed", 1);
            const label nD = dict.getOrDefault<label>("nDirections", 1);
            const scalar s = dict.getOrDefault<scalar>("spreading", 10);
            const scalar maxAngle =
                degToRad(dict.getOrDefault<scalar>("maxAngle", 60));

            if (nF < 1 || fMax <= fMin || fMin <= 0 || nD < 1)
            {
                FatalIOErrorInFunction(dict)
                    << "Irregular wave needs nFrequencies > 0,"
                    << " 0 < fMin < fMax and nDirections > 0"
                    << exit(FatalIOError);
            }

            // Frequencies: midpoints of nF equal bins
            const scalar df = (fMax - fMin)/nF;
            scalarField S(nF);
            scalar m0 = 0;
            forAll(S, i)
            {
                S[i] = jonswapShape(fMin + (i + 0.5)*df, Tp, gamma);
                m0 += S[i]*df;
            }

            // Directions: midpoints of nD equal bins in [-maxAngle, maxAngle]
            // with weights cos^{2s}(angle/2) summing to one
            scalarField angles(nD, Zero), D(nD, scalar(1));
            if (nD > 1)
            {
                const scalar dth = 2*maxAngle/nD;
                scalar sumD = 0;
                forAll(angles, j)
                {
                    angles[j] = -maxAngle + (j + 0.5)*dth;
                    D[j] = Foam::pow(Foam::cos(0.5*angles[j]), 2*s);
                    sumD += D[j];
                }
                D /= sumD;
            }

            // Scale to Hm0 = 4 sqrt(m0) = Hs, random phases
            const scalar scale = sqr(0.25*Hs)/max(m0, VSMALL);
            Random rng(seed);

            forAll(S, i)
            {
                forAll(angles, j)
                {
                    addComponent
                    (
                        Foam::sqrt(2*scale*S[i]*df*D[j]),
                        twoPi*(fMin + (i + 0.5)*df),
                        angles[j],
                        twoPi*rng.sample01<scalar>()
                    );
                }
            }
        }

        if (amp.empty())
        {
            FatalIOErrorInFunction(dict)
                << "Irregular wave without components"
                << exit(FatalIOError);
        }

        // Representative height Hm0 and the period of the largest component
        scalar m0 = 0;
        label iMax = 0;
        forAll(amp, i)
        {
            m0 += 0.5*sqr(amp[i]);
            if (amp[i] > amp[iMax])
            {
                iMax = i;
            }
        }
        H_ = 4*Foam::sqrt(m0);
        omega_ = omegas[iMax];
        T_ = twoPi/omega_;
        k_ = ks[iMax];
    }
    else if (theory == "streamFunction")
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
    else if (theory == "timeSeries")
    {
        timeSeries_ = true;

        const bool hasEta = dict.found("eta");
        paddle_ = dict.found("paddle");

        if (hasEta == paddle_)
        {
            FatalIOErrorInFunction(dict)
                << "timeSeries needs either eta (a record of the surface"
                << " elevation) or paddle (a record of the paddle"
                << " displacement)" << exit(FatalIOError);
        }

        series_ = Function1<scalar>::New(paddle_ ? "paddle" : "eta", dict);

        const word relation(dict.getOrDefault<word>("velocity", "longWave"));
        if (relation == "linear")
        {
            linearRelation_ = true;
        }
        else if (relation != "longWave")
        {
            FatalIOErrorInFunction(dict)
                << "velocity " << relation << " unknown; use longWave or"
                << " linear" << exit(FatalIOError);
        }

        celerity_ = dict.getOrDefault<scalar>("celerity", 0);
        timeShift_ = dict.getOrDefault<scalar>("timeShift", 0);
        scale_ = dict.getOrDefault<scalar>("scale", 1);
        verticalVelocity_ = dict.getOrDefault<bool>("verticalVelocity", true);

        // Default step of the time derivatives: half the smallest sampling
        // interval of a table (the slope of its linear interpolation,
        // averaged over one interval), otherwise 1 ms
        const auto* tablePtr =
            dynamic_cast<const Function1Types::TableBase<scalar>*>
            (
                series_.get()
            );

        scalar t0 = 0, t1 = 0;
        if (tablePtr)
        {
            const scalarField ts(tablePtr->x());
            scalar dtMin = GREAT;
            for (label i = 1; i < ts.size(); ++i)
            {
                if (ts[i] - ts[i-1] > SMALL)
                {
                    dtMin = min(dtMin, ts[i] - ts[i-1]);
                }
            }
            if (dtMin < GREAT)
            {
                dtDiff_ = 0.5*dtMin;
            }
            if (ts.size())
            {
                t0 = ts.first();
                t1 = ts.last();
            }
        }
        dtDiff_ = dict.getOrDefault<scalar>("differenceStep", dtDiff_);

        if (dtDiff_ <= 0 || (celerity_ < 0))
        {
            FatalIOErrorInFunction(dict)
                << "timeSeries needs differenceStep > 0 and celerity >= 0"
                << exit(FatalIOError);
        }

        Info<< "    time series of the "
            << (paddle_ ? "paddle displacement" : "surface elevation")
            << " (" << series_->type() << ")";
        if (tablePtr)
        {
            Info<< ", record from " << t0 << " to " << t1 << " s";
        }
        Info<< ", time shift " << timeShift_ << " s, scale " << scale_
            << ", velocity " << relation << ", difference step " << dtDiff_
            << " s" << endl;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "theory " << theory << " unknown; use linear, irregular,"
            << " streamFunction or timeSeries" << exit(FatalIOError);
    }

    amp_.transfer(amp);
    omegas_.transfer(omegas);
    ks_.transfer(ks);
    phases_.transfer(phases);
    dirs_.transfer(dirs);

    // Return current: minus the Stokes transport g a^2/(2c) over the depth
    if (!streamFunction_ && dict.getOrDefault<bool>("returnCurrent", true))
    {
        forAll(amp_, i)
        {
            const scalar c = omegas_[i]/ks_[i];
            uReturn_ -= g_*sqr(amp_[i])/(2*c*h_)*dirs_[i];
        }
    }

    if (theory == "irregular")
    {
        Info<< "    irregular wave: " << amp_.size() << " components, Hm0 "
            << H_ << ", period of the largest component " << T_ << endl;
    }
    if (mag(uReturn_) > 0)
    {
        Info<< "    return current " << uReturn_ << endl;
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


Foam::scalar Foam::regularWave::jonswapShape
(
    const scalar f,
    const scalar Tp,
    const scalar gamma
)
{
    const scalar fp = 1/Tp;
    const scalar sigma = (f <= fp ? 0.07 : 0.09);
    const scalar r = Foam::exp(-sqr(f - fp)/(2*sqr(sigma*fp)));

    return
        Foam::pow(f, -5)*Foam::exp(-1.25*Foam::pow(fp/f, 4))
       *Foam::pow(gamma, r);
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::scalar Foam::regularWave::theta(const point& x, const scalar t) const
{
    return k_*((x - origin_) & direction_) - omega_*t + phase_;
}


Foam::scalar Foam::regularWave::recordTime
(
    const point& x,
    const scalar t
) const
{
    // The record holds at the origin; elsewhere it arrives after the
    // travel time of a long wave along the direction
    const scalar c0 = celerity_ > 0 ? celerity_ : Foam::sqrt(g_*h_);

    return t + timeShift_ - ((x - origin_) & direction_)/c0;
}


Foam::scalar Foam::regularWave::velocityFromEta(const scalar eta) const
{
    if (linearRelation_)
    {
        return eta*Foam::sqrt(g_/h_);
    }

    // Flux c eta through the depth h + eta (kept positive)
    const scalar D = max(h_ + eta, 0.01*h_);
    const scalar c = celerity_ > 0 ? celerity_ : Foam::sqrt(g_*D);

    return c*eta/D;
}


Foam::scalar Foam::regularWave::etaFromVelocity(const scalar u) const
{
    if (linearRelation_)
    {
        return u*Foam::sqrt(h_/g_);
    }

    if (celerity_ > 0)
    {
        // u = c eta/(h + eta)
        return u*h_/max(celerity_ - u, 0.01*celerity_);
    }

    // u = eta sqrt(g/(h + eta)): g eta^2 - u^2 eta - u^2 h = 0, the root
    // with the sign of u
    const scalar u2 = sqr(u);
    const scalar root = Foam::sqrt(sqr(u2) + 4*g_*h_*u2);

    return (u2 + (u < 0 ? -root : root))/(2*g_);
}


Foam::scalar Foam::regularWave::seriesEta(const scalar tau) const
{
    if (paddle_)
    {
        return etaFromVelocity(seriesVelocity(tau));
    }

    return scale_*series_->value(tau);
}


Foam::scalar Foam::regularWave::seriesVelocity(const scalar tau) const
{
    if (paddle_)
    {
        // Paddle velocity: central difference of the displacement
        return
            scale_
           *(series_->value(tau + dtDiff_) - series_->value(tau - dtDiff_))
           /(2*dtDiff_);
    }

    return velocityFromEta(seriesEta(tau));
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

    return 0.5*(1 - Foam::cos(pi*t/rampTime_));
}


Foam::scalar Foam::regularWave::phaseSpeed() const
{
    if (timeSeries_ || k_ <= 0)
    {
        return celerity_ > 0 ? celerity_ : Foam::sqrt(g_*h_);
    }

    return omega_/k_;
}


Foam::scalar Foam::regularWave::eta(const point& x, const scalar t) const
{
    if (timeSeries_)
    {
        return ramp(t)*seriesEta(recordTime(x, t));
    }

    if (streamFunction_)
    {
        const scalar th = theta(x, t);

        scalar e = 0;
        forAll(E_, j)
        {
            e += E_[j]*Foam::cos(j*th);
        }
        return ramp(t)*e;
    }

    const vector dx(x - origin_);

    scalar e = 0;
    forAll(amp_, i)
    {
        e += amp_[i]*Foam::cos(ks_[i]*(dx & dirs_[i]) - omegas_[i]*t + phases_[i]);
    }

    return ramp(t)*e;
}


Foam::vector Foam::regularWave::velocity(const point& x, const scalar t) const
{
    // Height above the bottom, limited to the water column
    const scalar zb = max(x.z() - swl_ + h_, 0);
    const scalar r = ramp(t);

    if (timeSeries_)
    {
        const scalar tau = recordTime(x, t);
        const scalar u = seriesVelocity(tau);

        // Continuity for a wave of permanent form travelling along the
        // direction at the celerity c: w = -z_b du/dx = (z_b/c) du/dt
        scalar w = 0;
        if (verticalVelocity_)
        {
            const scalar dudt =
                (seriesVelocity(tau + dtDiff_) - seriesVelocity(tau - dtDiff_))
               /(2*dtDiff_);
            const scalar D = max(h_ + seriesEta(tau), 0.01*h_);
            const scalar c = celerity_ > 0 ? celerity_ : Foam::sqrt(g_*D);

            w = min(zb, D)*dudt/c;
        }

        return r*(u*direction_ + vector(0, 0, w));
    }

    if (streamFunction_)
    {
        const scalar th = theta(x, t);

        scalar u = uMean_, w = 0;

        forAll(B_, i)
        {
            const label j = i + 1;
            const scalar jk = j*k_;
            const scalar ch = Foam::cosh(jk*h_);

            u += jk*B_[i]*Foam::cosh(jk*zb)/ch*Foam::cos(j*th);
            w += jk*B_[i]*Foam::sinh(jk*zb)/ch*Foam::sin(j*th);
        }

        return r*(u*direction_ + vector(0, 0, w));
    }

    const vector dx(x - origin_);

    vector uh(Zero);
    scalar w = 0;

    forAll(amp_, i)
    {
        const scalar k = ks_[i];
        const scalar th = k*(dx & dirs_[i]) - omegas_[i]*t + phases_[i];

        // cosh(k zb)/sinh(k h) and sinh(k zb)/sinh(k h) without overflow
        const scalar ep = Foam::exp(k*(zb - h_));
        const scalar em = Foam::exp(-k*(zb + h_));
        const scalar aw = amp_[i]*omegas_[i]/(1 - Foam::exp(-2*k*h_));

        uh += aw*(ep + em)*Foam::cos(th)*dirs_[i];
        w += aw*(ep - em)*Foam::sin(th);
    }

    return r*(uh + vector(0, 0, w)) + sqr(r)*uReturn_;
}


// ************************************************************************* //
