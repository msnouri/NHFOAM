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
    Foam::interfacePressureLduMesh

Description
    Addressing of the interface-pressure matrix (see interfacePressure):
    the unknowns are the cells and one surface unknown per column.  The
    patches (processor interfaces) are those of the fvMesh, forwarded
    without taking ownership.  The vertical lines of unknowns (surface
    first) are kept for the sigmaColumn preconditioner.

\*---------------------------------------------------------------------------*/

#ifndef interfacePressureLduMesh_H
#define interfacePressureLduMesh_H

#include "lduPrimitiveMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class interfacePressureLduMesh Declaration
\*---------------------------------------------------------------------------*/

class interfacePressureLduMesh
:
    public lduPrimitiveMesh
{
    // Private data

        const fvMesh& mesh_;

        //- Unknowns of each column, surface first
        labelListList columnUnknowns_;

        //- Coupling between consecutive unknowns of each column
        labelListList columnCouplings_;


public:

    // Constructors

        interfacePressureLduMesh
        (
            const fvMesh& mesh,
            const label nUnknowns,
            labelList& l,
            labelList& u,
            labelListList&& columnUnknowns,
            labelListList&& columnCouplings
        )
        :
            lduPrimitiveMesh(nUnknowns, l, u, mesh.comm(), true),
            mesh_(mesh),
            columnUnknowns_(std::move(columnUnknowns)),
            columnCouplings_(std::move(columnCouplings))
        {}


    // Member Functions

        const labelListList& columnUnknowns() const
        {
            return columnUnknowns_;
        }

        const labelListList& columnCouplings() const
        {
            return columnCouplings_;
        }

        virtual lduInterfacePtrsList interfaces() const
        {
            return mesh_.interfaces();
        }

        virtual const labelUList& patchAddr(const label i) const
        {
            return mesh_.lduAddr().patchAddr(i);
        }

        virtual const lduSchedule& patchSchedule() const
        {
            return mesh_.lduAddr().patchSchedule();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
