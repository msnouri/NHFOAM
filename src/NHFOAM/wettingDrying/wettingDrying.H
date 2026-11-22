/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | NHFOAM: non-hydrostatic sigma-layer wave model
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of NHFOAM, built on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::wettingDrying

Description
    Wetting and drying of the water columns, after NHWAVE (Ma, Shi & Kirby
    2012, section 3.4), adapted to the implicit free surface of NHFOAM.

    Mask.  Dry columns keep a thin film of depth minDepth
    (eta = z_b + minDepth) and do not take part in the flow: all their side
    faces are closed and their velocity is held at zero.  The mask is
    updated at the start of every time step from the free surface of
    level n, as NHWAVE's update_mask:
      - flooding: a dry column becomes wet if a wet neighbour has a higher
        free surface and more than the film;
      - drying: a wet column whose depth is down to the film becomes dry.
    A side face is open when the columns on both sides are wet.  A closed
    face carries no flux, and the pressure gradient treats it as a wall,
    so still water stays at rest against a dry beach.  A wet column whose
    side faces are all closed (a puddle) is also held at rest.

    Positivity.  Within the time step the implicit pressure solve can draw
    more water from a shallow column than it holds; the side outflows of
    such a column are then scaled down (positivity limiter,
    interfacePressure::fluxes), which keeps eta^{n+1} >= z_b + minDepth and
    conserves the volume.

    Water depth at the side faces.  The mesh face between two columns has
    the depth of the point lines, about the mean of the two columns.  With
    that centred depth the front of a flow onto a dry bed advances too
    slowly and bores in thin water oscillate.  As in the semi-implicit
    schemes of Stelling and Duinmeijer (2003) and SWASH, the side fluxes use
    the depth on the upwind side instead,
        D_f = eta_f - z_b,f,
    with eta_f the upwind surface (faceDepth upwind) or its TVD
    reconstruction eta_C + psi(r)/2 (eta_D - eta_C) (minmod, vanLeer, MC;
    second order where the surface is smooth; upwind next to a column with
    less than two films, where the surface ends), and z_b,f the bottom of
    the face.  The fluxes of all layers are scaled by f = D_f/D_mesh,
    consistently in the pressure equation, the convecting flux and the face
    velocity (which stays the flux over the face depth).  NHWAVE gets the
    same effect from its MUSCL-TVD reconstruction and HLL fluxes.

    Supercritical faces.  Where the upwind column is supercritical
    (depth-averaged normal velocity above upwindFroude sqrt(g D)) the
    predicted flux of a side face is that of the upwind cell, as an HLL
    flux whose wave speeds are both positive.

    Faces opened by flooding.  In the step after a column floods, the
    predicted flux of the face may only run into it, and where it is lower
    than the wet column the face velocity is at least that of a dam break
    from rest, 2/3 sqrt(g (D - minDepth)).

    Dictionary (in NHFOAMProperties):

        wettingDrying
        {
            minDepth          0.005;  // film of the dry columns [m]
            tolerance         1e-6;   // drying when D - minDepth <= this
            faceDepth         vanLeer;// depth of the side fluxes: mesh |
                                      // upwind | minmod | vanLeer | MC
            faceDepthBelow    1e15;   // ... only where a column is
                                      // shallower than this [m]
            upwindFroude      1;      // upwind predicted flux above this
                                      // Froude number (0: off)
            hydrostaticDepth  0;      // q = 0 in columns shallower than
                                      // this (0: never)
        }

    Wetting and drying needs pressureLocation interface.

SourceFiles
    wettingDrying.C

\*---------------------------------------------------------------------------*/

#ifndef wettingDrying_H
#define wettingDrying_H

#include "fvMesh.H"
#include "sigmaColumns.H"
#include "boolList.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class wettingDrying Declaration
\*---------------------------------------------------------------------------*/

class wettingDrying
{
    // Private data

        const fvMesh& mesh_;

        const sigmaColumns& columns_;

        //- Wetting and drying switched on
        bool active_;

        //- Film depth of the dry columns
        scalar minDepth_;

        //- Tolerance of the drying test, D - minDepth <= tol
        scalar tol_;

        //- Columns shallower than this are hydrostatic (0: never)
        scalar hydrostaticDepth_;

        //- Water depth of the side-face fluxes
        enum class faceDepthType { mesh, upwind, minmod, vanLeer, MC };
        faceDepthType faceDepth_;

        //- The face depth is changed only where one of the columns is
        //  shallower than this
        scalar faceDepthBelow_;

        //- Side faces whose upwind column is supercritical (depth-averaged
        //  Froude number above this) take the predicted flux from the
        //  upwind cell, as the HLL flux of a supercritical face (0: off)
        scalar upwindFroude_;

        //- Wet columns
        boolList colWet_;

        //- Side faces between two columns (internal faces)
        boolList isSideFace_;

        //- Open internal faces (layer faces are always open)
        boolList faceOpen_;

        //- Open faces of the coupled (processor) side patches
        List<boolList> patchFaceOpen_;

        //- Number of open side faces of each column
        labelList nOpenSides_;

        //- Number of columns that flooded / dried in the last update
        label nFlooded_;
        label nDried_;


    // Private Member Functions

        //- Column values of the neighbours across the coupled patches
        //  (per patch face), from the column values of this processor
        List<scalarField> patchNeighbourColumnValues
        (
            const scalarField& colValues
        ) const;

        //- Open/closed state of the side faces from the column mask
        void setFaceMask();

        //- No copy construct
        wettingDrying(const wettingDrying&) = delete;

        //- No copy assignment
        void operator=(const wettingDrying&) = delete;


public:

    // Constructors

        //- Construct from the mesh, the columns and NHFOAMProperties
        wettingDrying
        (
            const fvMesh& mesh,
            const sigmaColumns& columns,
            const dictionary& nhProperties
        );


    // Member Functions

        bool active() const
        {
            return active_;
        }

        scalar minDepth() const
        {
            return minDepth_;
        }

        scalar hydrostaticDepth() const
        {
            return hydrostaticDepth_;
        }

        //- Side-face fluxes scaled to the upwind (TVD) water depth
        bool scaleFaceDepth() const
        {
            return active_ && faceDepth_ != faceDepthType::mesh;
        }

        const boolList& colWet() const
        {
            return colWet_;
        }

        const boolList& faceOpen() const
        {
            return faceOpen_;
        }

        const List<boolList>& patchFaceOpen() const
        {
            return patchFaceOpen_;
        }

        label nFlooded() const
        {
            return nFlooded_;
        }

        label nDried() const
        {
            return nDried_;
        }

        //- Initial state: columns with no more than the film are dry and
        //  get eta = z_b + minDepth
        void initialise(scalarField& etaCol);

        //- Update the mask from eta of level n (flooding and drying, after
        //  NHWAVE's update_mask); the columns that dry are set to the film.
        //  Returns the columns that changed state.
        boolList update(scalarField& etaCol);

        //- Lowest admissible free surface per column, z_b + minDepth
        tmp<scalarField> etaMin() const;

        //- Mark (1) the cells of the columns without horizontal flow: dry
        //  columns, and wet columns whose side faces are all closed (a
        //  puddle surrounded by dry columns, whose velocity could only
        //  carry it away from the neighbour it floods next)
        void dryCells(scalarField& cellFlag) const;

        //- Side faces whose upwind column is supercritical (depth-averaged
        //  velocity normal to the face above upwindFroude sqrt(g D)) take
        //  the predicted flux HbyA.Sf of the upwind cell instead of the
        //  interpolated one, as an HLL flux whose wave speeds are both
        //  positive.  Without it the flux into the thin front of a run-up
        //  carries the mean of the fast and the slow (or resting) column.
        void upwindSupercritical
        (
            surfaceScalarField& phiHbyA,
            const volVectorField& HbyA,
            const volVectorField& U,
            const scalarField& etaCol,
            const scalar g
        ) const;

        //- On the side faces that opened in the last update (one column
        //  newly wet), the predicted flux may only run from the column that
        //  was wet into the newly wet one: the film of the newly wet column
        //  cannot feed the flow, and a velocity of the wet column pointing
        //  away from the face must not keep the face dry (the positivity
        //  limiter would cut the flux to zero and the column would dry
        //  again).  Where the newly wet column is lower, the predicted face
        //  velocity is at least that of a dam break from rest at the dam,
        //  2/3 sqrt(g (D - minDepth)) of the wet column (NHWAVE's HLL solver
        //  gives a wet-dry face the corresponding flux): a puddle left on a
        //  slope drains downhill instead of sitting there.
        //  changed: columns that changed state in the update; etaCol: the
        //  surface at level n.
        void limitOpenedFaces
        (
            surfaceScalarField& phiHbyA,
            const boolList& changed,
            const scalarField& etaCol,
            const scalar g
        ) const;

        //- Scale factors f = D_f/D_mesh of the side-face fluxes (internal
        //  faces, and per patch; 1 on the layer faces and on patches that
        //  are not coupled, 0 on closed faces).  etaCol: the surface of the
        //  current mesh; phiRef: flux that gives the upwind side.
        void faceFactors
        (
            const scalarField& etaCol,
            const surfaceScalarField& phiRef,
            scalarField& f,
            List<scalarField>& fp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
