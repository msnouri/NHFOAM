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

Application
    nhSetWaves

Description
    Initial free surface (eta) and velocity (U) for NHFoam.

    Reads system/nhSetWavesDict, overwrites the values of eta and U in the
    start-time directory (boundary conditions are kept):

    \verbatim
        type            standing;   // standing|solitary|linear|circularMode
        stillWaterLevel 0;          // optional

        // standing: eta = a cos(k (x - x0).d), U = 0
        amplitude       0.01;
        wavelength      2;
        origin          (0 0 0);
        direction       (1 0 0);

        // solitary: eta = H sech^2(K (x - x0).d), K = sqrt(3H/(4h^3)),
        //   u = c eta/(h + eta), w = -(z + h) du/dx, c = sqrt(g (h + H))
        waveHeight      0.1;
        depth           1;
        origin          (10 0 0);
        direction       (1 0 0);

        // linear: Airy wave at t = 0 (entries of linearWave)

        // circularMode: eta = a J_m(k r) cos(m theta) about 'centre'
        amplitude       0.005;
        waveNumber      1.841183781;    // k R = root of J_m'
        mode            1;
        centre          (0 0 0);
    \endverbatim

    Velocities are evaluated at the cell heights the sigma mesh will have
    once NHFoam places it on the new free surface.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sigmaColumns.H"
#include "linearWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set the initial free surface and velocity for NHFoam"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "nhSetWavesDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary nhProperties
    (
        IOobject
        (
            "NHFoamProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    #include "readGravitationalAcceleration.H"
    const scalar gMag = mag(g.value());

    sigmaColumns columns(mesh, nhProperties);

    volScalarField eta
    (
        IOobject
        (
            "eta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const word type(dict.get<word>("type"));
    const scalar swl
    (
        dict.getOrDefault<scalar>
        (
            "stillWaterLevel",
            nhProperties.getOrDefault<scalar>("stillWaterLevel", 0)
        )
    );

    vector dir(dict.getOrDefault<vector>("direction", vector(1, 0, 0)));
    dir.z() = 0;
    dir /= mag(dir);
    const point x0(dict.getOrDefault<point>("origin", point::zero));

    const vectorField& colC = columns.columnCentre();
    const scalarField& colZb = columns.columnBottomZ();
    const labelList& cellColumn = columns.cellColumn();
    const vectorField& C = mesh.C().primitiveField();

    // Current surface height of each column (surface-face centre)
    const scalarField colZs
    (
        mesh.boundaryMesh()[columns.surfacePatchID()].faceCentres()
        .component(vector::Z)
    );

    scalarField etaCol(columns.nColumns(), swl);
    vectorField& Uc = U.primitiveFieldRef();
    Uc = Zero;

    // Height of a cell centre once the mesh follows the new surface
    auto newCellZ = [&](const label celli, const scalar etaNew) -> scalar
    {
        const label coli = cellColumn[celli];
        const scalar sigma =
            (C[celli].z() - colZb[coli])/(colZs[coli] - colZb[coli]);

        return colZb[coli] + sigma*(etaNew - colZb[coli]);
    };

    if (type == "standing")
    {
        const scalar a(dict.get<scalar>("amplitude"));
        const scalar L(dict.get<scalar>("wavelength"));
        const scalar k = constant::mathematical::twoPi/L;

        forAll(etaCol, coli)
        {
            etaCol[coli] = swl + a*Foam::cos(k*((colC[coli] - x0) & dir));
        }

        Info<< "Standing wave: amplitude " << a << ", wavelength " << L
            << endl;
    }
    else if (type == "solitary")
    {
        const scalar H(dict.get<scalar>("waveHeight"));
        const scalar h(dict.get<scalar>("depth"));
        const scalar K = Foam::sqrt(3*H/(4*pow3(h)));
        const scalar c = Foam::sqrt(gMag*(h + H));

        forAll(etaCol, coli)
        {
            const scalar xi = (colC[coli] - x0) & dir;
            etaCol[coli] = swl + H/sqr(Foam::cosh(K*xi));
        }

        forAll(C, celli)
        {
            const scalar xi = (C[celli] - x0) & dir;
            const scalar sech2 = 1/sqr(Foam::cosh(K*xi));
            const scalar e = H*sech2;
            const scalar dedx = -2*H*K*sech2*Foam::tanh(K*xi);
            const scalar u = c*e/(h + e);
            const scalar dudx = c*h*dedx/sqr(h + e);
            const scalar zNew = newCellZ(celli, swl + e);

            Uc[celli] = u*dir + vector(0, 0, -(zNew - swl + h)*dudx);
        }

        Info<< "Solitary wave: height " << H << ", depth " << h
            << ", celerity " << c << endl;
    }
    else if (type == "linear")
    {
        dictionary waveDict(dict);
        waveDict.set("rampTime", 0);

        linearWave wave(waveDict, gMag, swl);

        forAll(etaCol, coli)
        {
            etaCol[coli] = swl + wave.eta(colC[coli], 0);
        }

        forAll(C, celli)
        {
            const label coli = cellColumn[celli];
            point pt(C[celli]);
            pt.z() = newCellZ(celli, etaCol[coli]);
            Uc[celli] = wave.velocity(pt, 0);
        }

        Info<< "Linear wave: k = " << wave.k() << ", omega = "
            << wave.omega() << endl;
    }
    else if (type == "circularMode")
    {
        // Standing mode of a circular basin: eta = a J_m(k r) cos(m theta)
        const scalar a(dict.get<scalar>("amplitude"));
        const scalar k(dict.get<scalar>("waveNumber"));
        const label m(dict.getOrDefault<label>("mode", 1));
        const point centre(dict.getOrDefault<point>("centre", point::zero));

        forAll(etaCol, coli)
        {
            const vector d(colC[coli] - centre);
            const scalar r = Foam::sqrt(sqr(d.x()) + sqr(d.y()));
            const scalar theta = Foam::atan2(d.y(), d.x());

            etaCol[coli] = swl + a*Foam::jn(m, k*r)*Foam::cos(m*theta);
        }

        Info<< "Circular basin mode m = " << m << ", k = " << k << endl;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown type " << type
            << ". Valid types: standing solitary linear circularMode"
            << exit(FatalIOError);
    }

    columns.setCellValues(etaCol, eta.primitiveFieldRef());
    eta.correctBoundaryConditions();
    U.correctBoundaryConditions();

    Info<< "Writing eta and U to " << runTime.timeName() << nl << endl;
    eta.write();
    U.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
