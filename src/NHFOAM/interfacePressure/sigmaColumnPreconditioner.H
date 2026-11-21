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
    Foam::sigmaColumnPreconditioner

Description
    Block-Jacobi preconditioner with one block per water column (line
    relaxation in the vertical) for the interface-pressure matrix: every
    column's tridiagonal block (surface unknown and the layer interfaces) is
    solved exactly.  The blocks of a symmetric positive definite matrix are
    symmetric positive definite, so the preconditioner can be used with PCG
    whether or not the matrix is diagonally dominant.  It captures the
    strong vertical coupling of thin layers and of breaking (HFA) columns.

    On any other matrix (e.g. the cell-centred pressure of
    pressureLocation cell) it falls back to DIC (symmetric matrices) or
    DILU (asymmetric), so a case can name it in both layouts.

SourceFiles
    sigmaColumnPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaColumnPreconditioner_H
#define sigmaColumnPreconditioner_H

#include "lduMatrix.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class sigmaColumnPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class sigmaColumnPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- Unknowns of each column, surface first (null on other meshes)
        const labelListList* columnUnknownsPtr_;

        //- Thomas algorithm: reciprocal of the modified diagonal and the
        //  off-diagonal coupling to the next unknown down the column
        solveScalarField rDiag_;
        solveScalarField offDiag_;

        //- DIC or DILU on matrices not built by interfacePressure
        autoPtr<lduMatrix::preconditioner> fallback_;


    // Private Member Functions

        //- No copy construct
        sigmaColumnPreconditioner(const sigmaColumnPreconditioner&) = delete;

        //- No copy assignment
        void operator=(const sigmaColumnPreconditioner&) = delete;


public:

    //- Runtime type information
    TypeName("sigmaColumn");


    // Constructors

        sigmaColumnPreconditioner
        (
            const lduMatrix::solver& sol,
            const dictionary& solverControlsUnused
        );


    //- Destructor
    virtual ~sigmaColumnPreconditioner() = default;


    // Member Functions

        //- Return wA, the preconditioned form of residual rA
        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt = 0
        ) const;

        //- Transposed: the blocks are symmetric
        virtual void preconditionT
        (
            solveScalarField& wT,
            const solveScalarField& rT,
            const direction cmpt = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
