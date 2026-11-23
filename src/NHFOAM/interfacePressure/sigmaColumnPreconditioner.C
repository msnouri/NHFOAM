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

#include "sigmaColumnPreconditioner.H"
#include "interfacePressureLduMesh.H"
#include "DICPreconditioner.H"
#include "DILUPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sigmaColumnPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<sigmaColumnPreconditioner>
        addsigmaColumnPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<sigmaColumnPreconditioner>
        addsigmaColumnPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace
{
    const Foam::interfacePressureLduMesh* columnMesh
    (
        const Foam::lduMatrix& matrix
    )
    {
        return dynamic_cast<const Foam::interfacePressureLduMesh*>
        (
            &matrix.mesh()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaColumnPreconditioner::sigmaColumnPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& solverControls
)
:
    lduMatrix::preconditioner(sol),
    columnUnknownsPtr_(nullptr),
    rDiag_(),
    offDiag_(),
    fallback_()
{
    const lduMatrix& matrix = solver_.matrix();
    const interfacePressureLduMesh* cmPtr = columnMesh(matrix);

    if (!cmPtr)
    {
        if (matrix.symmetric())
        {
            fallback_.reset(new DICPreconditioner(sol, solverControls));
        }
        else
        {
            fallback_.reset(new DILUPreconditioner(sol, solverControls));
        }
        return;
    }

    columnUnknownsPtr_ = &cmPtr->columnUnknowns();
    rDiag_.setSize(matrix.diag().size(), 0);
    offDiag_.setSize(matrix.diag().size(), 0);

    const labelListList& columnUnknowns = *columnUnknownsPtr_;
    const scalarField& diag = matrix.diag();
    const scalarField& upper = matrix.upper();
    const labelListList& couplings = cmPtr->columnCouplings();

    // LDL^T factorisation of each (symmetric) tridiagonal block, from the
    // surface down: d'_0 = d_0, d'_i = d_i - e_{i-1}^2/d'_{i-1}, with e_i
    // the coefficient between unknowns i and i+1.  Stored: 1/d'_i, e_i.
    forAll(columnUnknowns, coli)
    {
        const labelList& unk = columnUnknowns[coli];
        const labelList& cpl = couplings[coli];

        solveScalar dPrev = 1;
        solveScalar ePrev = 0;

        forAll(unk, i)
        {
            const solveScalar d = diag[unk[i]] - ePrev*ePrev/dPrev;
            rDiag_[unk[i]] = 1/d;
            ePrev = (i < unk.size() - 1 ? upper[cpl[i]] : 0);
            offDiag_[unk[i]] = ePrev;
            dPrev = d;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaColumnPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction cmpt
) const
{
    if (fallback_)
    {
        fallback_->precondition(wA, rA, cmpt);
        return;
    }

    const labelListList& columnUnknowns = *columnUnknownsPtr_;

    forAll(columnUnknowns, coli)
    {
        const labelList& unk = columnUnknowns[coli];
        const label n = unk.size();

        // L y = r
        wA[unk[0]] = rA[unk[0]];

        for (label i = 1; i < n; ++i)
        {
            const label p = unk[i-1];
            wA[unk[i]] = rA[unk[i]] - offDiag_[p]*rDiag_[p]*wA[p];
        }

        // D L^T w = y
        wA[unk[n-1]] *= rDiag_[unk[n-1]];

        for (label i = n - 2; i >= 0; --i)
        {
            const label u = unk[i];
            wA[u] = (wA[u] - offDiag_[u]*wA[unk[i+1]])*rDiag_[u];
        }
    }
}


void Foam::sigmaColumnPreconditioner::preconditionT
(
    solveScalarField& wT,
    const solveScalarField& rT,
    const direction cmpt
) const
{
    if (fallback_)
    {
        fallback_->preconditionT(wT, rT, cmpt);
        return;
    }

    precondition(wT, rT, cmpt);
}


// ************************************************************************* //
