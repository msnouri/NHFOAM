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
    Foam::sigmaColumns

Description
    Vertical column structure of an extruded ("sigma") mesh.

    Every face of the free-surface patch defines one water column.  The cells
    below it, down to the bottom patch, are collected top to bottom.  Every
    point of the free-surface patch defines a vertical line of mesh points.

    The relative height of each point on its line (sigma: 0 at the bottom,
    1 at the free surface) is taken from the initial mesh, so any vertical
    grading of the mesh is kept when the free surface moves:

        z = zb + sigma*(eta - zb)

    The free-surface elevation eta is stored per column (one value per
    surface-patch face).  It is interpolated to the surface points with
    horizontal inverse-distance weights, synchronised across processors.

    Requirements on the mesh
      - extruded in the vertical (z) direction: hexahedra, prisms or
        extruded polyhedra stacked in columns between the surface and the
        bottom patch, with mesh points on vertical lines;
      - in parallel, whole columns on one processor (use the 'structured'
        decomposition method with the bottom or surface patch).

SourceFiles
    sigmaColumns.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaColumns_H
#define sigmaColumns_H

#include "fvMesh.H"
#include "dictionary.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sigmaColumns Declaration
\*---------------------------------------------------------------------------*/

class sigmaColumns
{
    // Private data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Name of the free-surface patch
        word surfacePatchName_;

        //- Name of the bottom patch
        word bottomPatchName_;

        //- Index of the free-surface patch
        label surfacePatchID_;

        //- Index of the bottom patch
        label bottomPatchID_;

        //- Cells of each column, from the surface downwards
        labelListList columnCells_;

        //- Column of each cell
        labelList cellColumn_;

        //- Layer of each cell (0 = surface layer)
        labelList cellLayer_;

        //- Bottom-patch mesh face of each column
        labelList columnBottomFace_;

        //- Mesh points of each vertical line (one line per surface-patch
        //  point, same order as the patch's local points), surface first
        labelListList lineMeshPoints_;

        //- Line of each mesh point
        labelList pointLine_;

        //- Sigma coordinate of each mesh point (0 bottom, 1 surface)
        scalarField pointSigma_;

        //- Bottom elevation of each line
        scalarField lineBottomZ_;

        //- Horizontal (projected) area of each column
        scalarField columnArea_;

        //- Horizontal position of each column (surface-face centre, z = 0)
        vectorField columnCentre_;

        //- Bottom elevation of each column (bottom-face centre)
        scalarField columnBottomZ_;

        //- Inverse-distance weights for column-to-point interpolation,
        //  per surface-patch face and face point
        List<scalarList> faceToPointWeights_;


    // Private Member Functions

        //- Collect the cells of every column
        void findColumns();

        //- Collect the vertical lines of points and their sigma values
        void findLines();

        //- Column geometry and interpolation weights
        void calcGeometry();

        //- No copy construct
        sigmaColumns(const sigmaColumns&) = delete;

        //- No copy assignment
        void operator=(const sigmaColumns&) = delete;


public:

    // Constructors

        //- Construct from mesh and dictionary
        //  (keywords surfacePatch, bottomPatch)
        sigmaColumns(const fvMesh& mesh, const dictionary& dict);


    // Member Functions

        // Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            label surfacePatchID() const
            {
                return surfacePatchID_;
            }

            label bottomPatchID() const
            {
                return bottomPatchID_;
            }

            const polyPatch& surfacePatch() const
            {
                return mesh_.boundaryMesh()[surfacePatchID_];
            }

            //- Number of local columns
            label nColumns() const
            {
                return columnCells_.size();
            }

            const labelListList& columnCells() const
            {
                return columnCells_;
            }

            const labelList& cellColumn() const
            {
                return cellColumn_;
            }

            const labelList& cellLayer() const
            {
                return cellLayer_;
            }

            const scalarField& pointSigma() const
            {
                return pointSigma_;
            }

            const scalarField& columnArea() const
            {
                return columnArea_;
            }

            const vectorField& columnCentre() const
            {
                return columnCentre_;
            }

            const scalarField& columnBottomZ() const
            {
                return columnBottomZ_;
            }


        // Mapping

            //- Column values from a cell field (value of the surface cell)
            tmp<scalarField> columnValues(const scalarField& cellValues) const;

            //- Copy column values to every cell of each column
            void setCellValues
            (
                const scalarField& colValues,
                scalarField& cellValues
            ) const;

            //- Interpolate column values to the surface-patch points
            //  (parallel consistent)
            tmp<scalarField> surfacePointValues
            (
                const scalarField& colValues
            ) const;

            //- Mesh points for the free-surface elevation eta (per column).
            //  The water depth used for the point distribution is limited
            //  from below by minDepth.
            tmp<pointField> newPoints
            (
                const scalarField& eta,
                const scalar minDepth
            ) const;

            //- Local column nearest (horizontally) to a point,
            //  -1 if there are no local columns.  Returns the distance.
            label findNearestColumn(const point& pt, scalar& dist) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
