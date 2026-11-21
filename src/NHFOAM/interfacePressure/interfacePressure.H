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
    Foam::interfacePressure

Description
    Pressure on the layer interfaces of a sigma mesh, as in NHWAVE
    (Ma, Shi & Kirby 2012): the kinematic pressure p_rgh = g eta + q is
    defined on the sigma levels, with q = 0 exactly on the free surface.

    Layout (layers k = 0 at the surface to N-1 at the bottom)
      - unknowns: pI of cell k, p_rgh on the bottom interface of layer k,
        and the surface pressure p_s = g eta^{n+1} of every column;
      - velocities at the cell centres; the vertical pressure gradient of a
        cell is the compact difference of its two interfaces, the horizontal
        one comes from the pressure at the centre level, the mean of the two
        interfaces (Gauss gradient with these face values);
      - continuity on control volumes centred on the interfaces, from the
        centre of layer k to the centre of layer k+1 (half a layer at the
        bottom), with the vertical fluxes U.S at the cell centres and half of
        the side fluxes of the two layers;
      - the free surface from the top half layer (surface to the centre of
        the top layer):  A d(eta)/dt = U_0.S_0 - (side outflow of its upper
        half), solved implicitly with the rest (no gravity-wave time-step
        limit).

    The side flux of an interface control volume uses the pressure at the
    centre levels of its two half layers, so it couples interface k of a
    column with interfaces k-1, k and k+1 of the neighbouring column.  The
    matrix is therefore built on an extended addressing (the cells plus one
    surface unknown per column; mesh faces plus these "diagonal"
    couplings).  It is the Hessian of a sum of squares of fluxes plus the
    free-surface term, i.e. symmetric positive definite, but not diagonally
    dominant where layers are thicker than the horizontal spacing: it is
    solved with PCG and the sigmaColumn preconditioner (exact solves of the
    vertical lines; incomplete factorisations break down).

    Lagged (deferred) and converged over the pressure correctors: at
    processor boundaries the side faces are written in a compact form
    (couplings within one level only, which bounds the exact form so the
    iteration converges) and the rest is lagged; also lagged are the
    non-orthogonal correction and the slope (sigma) terms of the vertical
    fluxes.  The vertical flux through the sloping mid-plane of a cell is
    split as in the over-relaxed corrected Laplacian: the implicit
    coefficient carries |S_mid|^2, so the lagged part holds only the
    gradient along the layer and the iteration converges for steep layers
    (bore fronts).  The face fluxes returned are conservative in every
    cell: side fluxes from the pressure, layer fluxes integrated from the
    bottom.

    Parallel: columns must not be split between processors (structured
    decomposition); the processor patches are the interfaces of the solve.

SourceFiles
    interfacePressure.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePressure_H
#define interfacePressure_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "lduPrimitiveMesh.H"
#include "sigmaColumns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class interfacePressure Declaration
\*---------------------------------------------------------------------------*/

class interfacePressure
{
    // Private data

        const fvMesh& mesh_;

        const sigmaColumns& columns_;

        //- Unknown above / cell below each cell: the cell above, or the
        //  surface unknown nCells + column; -1 below the bottom
        labelList above_;
        labelList cellBelow_;

        //- Mesh face above / below each cell (layer face, or the surface /
        //  bottom patch face)
        labelList cellTopFace_;
        labelList cellBottomFace_;

        //- Internal faces between two cells of one column
        boolList isLayerFace_;

        //- Internal side faces: the side face directly above (-1 in the
        //  top layer)
        labelList faceAbove_;

        //- Coupled side patches: the patch face directly above (-1 in the
        //  top layer)
        List<labelList> patchFaceAbove_;


        // Extended addressing: unknowns are the cells and the surfaces

            //- Coupling of each internal mesh face
            labelList faceToExt_;

            //- Coupling of each cell with the unknown above
            labelList aboveExt_;

            //- Internal side faces: coupling of the two unknowns above,
            //  of the owner with the unknown above the neighbour, and of
            //  the unknown above the owner with the neighbour
            labelList pairAboveExt_;
            labelList diagOwnNbrAbove_;
            labelList diagOwnAboveNbr_;

            //- Lower and upper unknown of each coupling
            labelList extLower_;
            labelList extUpper_;

            //- Addressing of the extended matrix
            autoPtr<lduPrimitiveMesh> extMeshPtr_;


        // Coefficients of the last solve, kept for the reconstruction

            //- Implicit side coefficient rAUf |Sf| deltaCoeff
            scalarField cf_;
            FieldField<Field, scalar> cfB_;

            //- Lagged non-orthogonal correction of the side fluxes
            scalarField corr_;
            FieldField<Field, scalar> corrB_;


    // Private Member Functions

        //- Upward area vector of a horizontal (layer, surface, bottom) face
        vector upwardSf(const label facei) const;

        //- No copy construct
        interfacePressure(const interfacePressure&) = delete;

        //- No copy assignment
        void operator=(const interfacePressure&) = delete;


public:

    // Constructors

        //- Construct from mesh and column structure
        interfacePressure(const fvMesh& mesh, const sigmaColumns& columns);


    // Member Functions

        //- Assemble and solve the pressure equation for pI and the surface
        //  pressure ps (per column).  Lagged terms come from the latest
        //  state (pI, ps, pCell, gradP), which must be consistent.
        //      etaHatEff: etaHat + aEff*sigma*etaTarget  (per column)
        //      colS     : 1 + aEff*sigma                 (per column)
        //      factor   : vertical coefficient factor per column (HFA)
        //  The linear solver is taken from controls: PCG or PBiCGStab with
        //  the sigmaColumn (default) or diagonal preconditioner; other
        //  solvers and preconditioners are replaced, keeping the
        //  tolerances.
        solverPerformance solve
        (
            volScalarField& pI,
            scalarField& ps,
            const volScalarField& pCell,
            const volVectorField& gradP,
            const volScalarField& rAU,
            const surfaceScalarField& rAUf,
            const volVectorField& HbyA,
            const surfaceScalarField& phiHbyA,
            const volVectorField& U,
            const scalarField& etaHatEff,
            const scalarField& colS,
            const scalar aEff,
            const scalar g,
            const scalarField& factor,
            const dictionary& controls
        );

        //- Cell-centre pressure (mean of the two interfaces); the surface
        //  patch of pCell (mixed) is set to ps
        void cellPressure
        (
            const volScalarField& pI,
            const scalarField& ps,
            volScalarField& pCell
        ) const;

        //- Gauss gradient with the interface values on the layer faces,
        //  the surface and the bottom, and the cell-centre pressure
        //  (interpolated) on the side faces
        tmp<volVectorField> grad
        (
            const volScalarField& pI,
            const scalarField& ps,
            const volScalarField& pCell
        ) const;

        //- Absolute face fluxes after the solve (pCell: new cell pressure).
        //  Every cell is conservative; returns the surface flux per column.
        tmp<scalarField> fluxes
        (
            const volScalarField& pCell,
            const surfaceScalarField& phiHbyA,
            const volVectorField& U,
            surfaceScalarField& phi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
