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
    Foam::bottomMotion

Description
    Prescribed vertical displacement of the bottom, zb = zb0 + dz(x, y, t),
    for tsunami generation by seabed deformation or solid slides.

    Optional sub-dictionary bottomMotion of constant/NHFoamProperties.

    type uplift: a horizontal box raised (or lowered) by 'amplitude'
    \verbatim
    bottomMotion
    {
        type        uplift;
        min         (-0.305 -1 0);      // box corners (z ignored)
        max         ( 0.305  1 0);
        edgeWidth   0;                  // optional smooth (cosine) edge [m]
        amplitude   0.005;              // final displacement, + up [m]
        history     exponential;        // exponential | halfCosine | linear
        riseTime    0.2;                // exponential: dz = a(1 - exp(-1.11 t/riseTime))
                                        //   (Hammack 1973: 2/3 of the motion at riseTime)
                                        // halfCosine, linear: duration of the motion
    }
    \endverbatim

    type slide: rigid quasi-Gaussian solid slide moving down a plane slope
    (Grilli & Watts 2005; Enet & Grilli 2007), added on top of the initial
    bathymetry:
        dz = T/(1 - eps) [sech(kb xi) sech(kw chi) - eps]+ / cos(theta)
        kb = 2C/b, kw = 2C/w, C = acosh(1/eps)
    with xi along the slope from the slide centre and chi across it, and the
    law of motion s(t) = s0 ln(cosh(t/t0)), s0 = ut^2/a0, t0 = ut/a0.
    \verbatim
    bottomMotion
    {
        type        slide;
        centre      (0.551 0 0);        // initial centre (horizontal)
        direction   (1 0 0);            // horizontal downslope direction
        slopeAngle  15;                 // [deg]
        thickness   0.082;              // T [m]
        length      0.395;              // b [m]
        width       0.680;              // w [m]; <= 0 for a 2D slide
        epsilon     0.717;              // optional
        terminalVelocity     1.70;      // ut [m/s]
        initialAcceleration  1.21;      // a0 [m/s2]
    }
    \endverbatim

SourceFiles
    bottomMotion.C

\*---------------------------------------------------------------------------*/

#ifndef bottomMotion_H
#define bottomMotion_H

#include "dictionary.H"
#include "vectorField.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class bottomMotion Declaration
\*---------------------------------------------------------------------------*/

class bottomMotion
{
public:

    enum motionType { NONE, UPLIFT, SLIDE };


private:

    // Private data

        motionType type_;

        // Uplift
        point min_;
        point max_;
        scalar edgeWidth_;
        scalar amplitude_;
        word history_;
        scalar riseTime_;

        // Slide
        point centre0_;
        vector dir_;
        scalar theta_;
        scalar T_;
        scalar b_;
        scalar w_;
        scalar eps_;
        scalar ut_;
        scalar a0_;


    // Private Member Functions

        //- Time factor of the uplift (0 to 1)
        scalar upliftHistory(const scalar t) const;

        //- Horizontal shape factor of the uplift box (0 to 1)
        scalar upliftShape(const point& x) const;

        //- Displacement of the slide centre along the slope
        scalar slideDistance(const scalar t) const;


public:

    // Constructors

        //- Construct from the NHFoam dictionary (optional sub-dictionary
        //  bottomMotion)
        explicit bottomMotion(const dictionary& nhDict);


    // Member Functions

        bool active() const
        {
            return type_ != NONE;
        }

        //- Vertical bottom displacement at a horizontal position and time
        scalar displacement(const point& x, const scalar t) const;

        //- Displacement at a list of horizontal positions
        tmp<scalarField> displacement
        (
            const vectorField& x,
            const scalar t
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
