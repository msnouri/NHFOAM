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

\*---------------------------------------------------------------------------*/

#include "wettingDrying.H"
#include "volFields.H"
#include "calculatedFvPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::scalarField>
Foam::wettingDrying::patchNeighbourColumnValues
(
    const scalarField& colValues
) const
{
    volScalarField cellValues
    (
        IOobject
        (
            "wdColumnValues",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        calculatedFvPatchScalarField::typeName
    );

    columns_.setCellValues(colValues, cellValues.primitiveFieldRef());
    cellValues.correctBoundaryConditions();

    List<scalarField> nb(mesh_.boundary().size());

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (fvp.coupled())
        {
            nb[patchi] = cellValues.boundaryField()[patchi].patchNeighbourField();
        }
    }

    return nb;
}


void Foam::wettingDrying::setFaceMask()
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();

    nOpenSides_ = 0;

    forAll(nei, facei)
    {
        faceOpen_[facei] =
            !isSideFace_[facei]
         || (colWet_[cellCol[own[facei]]] && colWet_[cellCol[nei[facei]]]);

        if (isSideFace_[facei] && faceOpen_[facei])
        {
            ++nOpenSides_[cellCol[own[facei]]];
            ++nOpenSides_[cellCol[nei[facei]]];
        }
    }

    scalarField wetF(colWet_.size(), 0);
    forAll(colWet_, coli)
    {
        wetF[coli] = colWet_[coli] ? 1 : 0;
    }

    const List<scalarField> nbWet(patchNeighbourColumnValues(wetF));

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (fvp.coupled())
        {
            const labelUList& fc = fvp.faceCells();
            boolList& open = patchFaceOpen_[patchi];

            forAll(fc, i)
            {
                open[i] = colWet_[cellCol[fc[i]]] && nbWet[patchi][i] > 0.5;

                if (open[i])
                {
                    ++nOpenSides_[cellCol[fc[i]]];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wettingDrying::wettingDrying
(
    const fvMesh& mesh,
    const sigmaColumns& columns,
    const dictionary& nhProperties
)
:
    mesh_(mesh),
    columns_(columns),
    active_(nhProperties.found("wettingDrying")),
    minDepth_(1e-4),
    tol_(1e-6),
    hydrostaticDepth_(0),
    faceDepth_(faceDepthType::mesh),
    faceDepthBelow_(GREAT),
    upwindFroude_(0),
    colWet_(columns.nColumns(), true),
    isSideFace_(mesh.nInternalFaces(), false),
    faceOpen_(mesh.nInternalFaces(), true),
    patchFaceOpen_(mesh.boundary().size()),
    nOpenSides_(columns.nColumns(), 0),
    nFlooded_(0),
    nDried_(0)
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();

    forAll(nei, facei)
    {
        isSideFace_[facei] = (cellCol[own[facei]] != cellCol[nei[facei]]);
    }

    forAll(mesh_.boundary(), patchi)
    {
        patchFaceOpen_[patchi].setSize(mesh_.boundary()[patchi].size(), true);
    }

    if (active_)
    {
        const dictionary& dict = nhProperties.subDict("wettingDrying");

        minDepth_ = dict.get<scalar>("minDepth");
        tol_ = dict.getOrDefault<scalar>("tolerance", 1e-6);
        hydrostaticDepth_ =
            dict.getOrDefault<scalar>("hydrostaticDepth", 0);
        const word fd(dict.getOrDefault<word>("faceDepth", "vanLeer"));
        if (fd == "mesh") faceDepth_ = faceDepthType::mesh;
        else if (fd == "upwind") faceDepth_ = faceDepthType::upwind;
        else if (fd == "minmod") faceDepth_ = faceDepthType::minmod;
        else if (fd == "vanLeer") faceDepth_ = faceDepthType::vanLeer;
        else if (fd == "MC") faceDepth_ = faceDepthType::MC;
        else
        {
            FatalIOErrorInFunction(dict)
                << "faceDepth " << fd << " unknown; use mesh, upwind,"
                << " minmod, vanLeer or MC" << exit(FatalIOError);
        }
        faceDepthBelow_ = dict.getOrDefault<scalar>("faceDepthBelow", GREAT);
        upwindFroude_ = dict.getOrDefault<scalar>("upwindFroude", 1);

        if (minDepth_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "minDepth must be positive" << exit(FatalIOError);
        }

        Info<< "Wetting and drying: film depth " << minDepth_
            << " m, hydrostatic below " << hydrostaticDepth_
            << " m, face depth " << fd;
        if (faceDepthBelow_ < 0.1*GREAT)
        {
            Info<< " where a column is shallower than " << faceDepthBelow_
                << " m";
        }
        if (upwindFroude_ > 0)
        {
            Info<< ", upwind predicted flux above Froude " << upwindFroude_;
        }
        Info<< nl << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wettingDrying::initialise(scalarField& etaCol)
{
    if (!active_)
    {
        return;
    }

    const scalarField& zb = columns_.columnBottomZ();
    label nDry = 0;

    forAll(etaCol, coli)
    {
        if (etaCol[coli] - zb[coli] - minDepth_ <= tol_)
        {
            colWet_[coli] = false;
            etaCol[coli] = zb[coli] + minDepth_;
            ++nDry;
        }
        else
        {
            colWet_[coli] = true;
        }
    }

    setFaceMask();

    reduce(nDry, sumOp<label>());
    Info<< "Wetting and drying: " << nDry << " dry columns at the start"
        << nl << endl;
}


Foam::boolList Foam::wettingDrying::update(scalarField& etaCol)
{
    boolList changed(etaCol.size(), false);
    nFlooded_ = 0;
    nDried_ = 0;

    if (!active_)
    {
        return changed;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();
    const scalarField& zb = columns_.columnBottomZ();

    const boolList wetOld(colWet_);
    boolList wetNew(colWet_);
    auto flood = [&](const label coli)
    {
        wetNew[coli] = true;
    };

    // A wet column that can flood a neighbour: higher surface, more than
    // the film
    auto canFlood = [&](const scalar etaW, const scalar zbW, const scalar etaD)
    {
        return etaW > etaD && etaW - zbW - minDepth_ > tol_;
    };

    // Flooding (dry -> wet) across internal side faces
    forAll(nei, facei)
    {
        if (!isSideFace_[facei])
        {
            continue;
        }

        const label co = cellCol[own[facei]];
        const label cn = cellCol[nei[facei]];

        if (!wetOld[co] && wetOld[cn] && canFlood(etaCol[cn], zb[cn], etaCol[co]))
        {
            flood(co);
        }
        if (!wetOld[cn] && wetOld[co] && canFlood(etaCol[co], zb[co], etaCol[cn]))
        {
            flood(cn);
        }
    }

    // Flooding across processor boundaries
    {
        scalarField wetF(etaCol.size(), 0);
        forAll(wetOld, coli)
        {
            wetF[coli] = wetOld[coli] ? 1 : 0;
        }

        const List<scalarField> nbWet(patchNeighbourColumnValues(wetF));
        const List<scalarField> nbEta(patchNeighbourColumnValues(etaCol));
        const List<scalarField> nbZb(patchNeighbourColumnValues(zb));

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!fvp.coupled())
            {
                continue;
            }

            const labelUList& fc = fvp.faceCells();

            forAll(fc, i)
            {
                const label co = cellCol[fc[i]];

                if
                (
                    !wetOld[co]
                 && nbWet[patchi][i] > 0.5
                 && canFlood(nbEta[patchi][i], nbZb[patchi][i], etaCol[co])
                )
                {
                    wetNew[co] = true;
                }
            }
        }
    }

    // Drying (wet -> dry): down to the film
    forAll(wetOld, coli)
    {
        if (wetOld[coli] && etaCol[coli] - zb[coli] - minDepth_ <= tol_)
        {
            wetNew[coli] = false;
            etaCol[coli] = zb[coli] + minDepth_;
        }
    }

    forAll(wetOld, coli)
    {
        if (wetNew[coli] != wetOld[coli])
        {
            changed[coli] = true;
            wetNew[coli] ? ++nFlooded_ : ++nDried_;
        }
    }

    colWet_ = wetNew;
    setFaceMask();

    reduce(nFlooded_, sumOp<label>());
    reduce(nDried_, sumOp<label>());

    return changed;
}


Foam::tmp<Foam::scalarField> Foam::wettingDrying::etaMin() const
{
    return columns_.columnBottomZ() + minDepth_;
}


void Foam::wettingDrying::dryCells(scalarField& cellFlag) const
{
    const labelList& cellCol = columns_.cellColumn();

    forAll(cellFlag, celli)
    {
        const label coli = cellCol[celli];
        cellFlag[celli] = (colWet_[coli] && nOpenSides_[coli] > 0) ? 0 : 1;
    }
}


void Foam::wettingDrying::faceFactors
(
    const scalarField& etaCol,
    const surfaceScalarField& phiRef,
    scalarField& f,
    List<scalarField>& fp
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();
    const faceList& faces = mesh_.faces();
    const pointField& pts = mesh_.points();
    const labelList& pointLine = columns_.pointLine();
    const labelListList& linePts = columns_.lineMeshPoints();
    const scalarField& lineZb = columns_.lineBottomZ();
    const scalarField& zb = columns_.columnBottomZ();
    const vectorField& colC = columns_.columnCentre();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& V = mesh_.V();

    f.setSize(mesh_.nInternalFaces());
    f = 1;
    fp.setSize(mesh_.boundary().size());
    forAll(fp, patchi)
    {
        fp[patchi].setSize(mesh_.boundary()[patchi].size());
        fp[patchi] = 1;
    }

    // Depth and bottom of a mesh face: mean over the vertical point lines
    // of its points
    auto faceGeometry = [&](const face& fc, scalar& Dm, scalar& zbF)
    {
        scalar sD = 0, sZ = 0;
        forAll(fc, i)
        {
            const label l = pointLine[fc[i]];
            sD += pts[linePts[l][0]].z() - lineZb[l];
            sZ += lineZb[l];
        }
        Dm = sD/fc.size();
        zbF = sZ/fc.size();
    };

    const bool tvd =
        faceDepth_ != faceDepthType::mesh && faceDepth_ != faceDepthType::upwind;

    auto limiter = [&](const scalar r)
    {
        switch (faceDepth_)
        {
            case faceDepthType::minmod:
                return max(scalar(0), min(r, scalar(1)));
            case faceDepthType::vanLeer:
                return (r + mag(r))/(1 + mag(r));
            case faceDepthType::MC:
                return max(scalar(0), min(min(2*r, 0.5*(1 + r)), scalar(2)));
            default:
                return scalar(0);
        }
    };

    // Cell gradient of the column surface (Gauss, central, over the open
    // side faces; closed faces and the other boundaries are zero-gradient)
    const List<scalarField> nbEta(patchNeighbourColumnValues(etaCol));
    vectorField gradEta(mesh_.nCells(), Zero);

    if (tvd)
    {
        forAll(nei, facei)
        {
            if (isSideFace_[facei] && faceOpen_[facei])
            {
                const label o = own[facei];
                const label n = nei[facei];
                const scalar dEta =
                    0.5*(etaCol[cellCol[n]] - etaCol[cellCol[o]]);

                gradEta[o] += Sf[facei]*dEta/V[o];
                gradEta[n] += Sf[facei]*dEta/V[n];
            }
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!fvp.coupled())
            {
                continue;
            }

            const labelUList& fc = fvp.faceCells();
            const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];

            forAll(fc, i)
            {
                if (patchFaceOpen_[patchi][i])
                {
                    const label o = fc[i];
                    gradEta[o] +=
                        Sfp[i]*0.5*(nbEta[patchi][i] - etaCol[cellCol[o]])
                       /V[o];
                }
            }
        }
    }

    // Face surface from the upwind side C towards D.  Next to a column
    // within a film of drying the surface ends, and its film level is no
    // water surface to reconstruct towards: upwind there.
    auto etaFace = [&]
    (
        const scalar etaC,
        const scalar etaD,
        const vector& gradC,
        const vector& d,
        const scalar Dmin
    )
    {
        if (!tvd || Dmin < 2*minDepth_)
        {
            return etaC;
        }

        const scalar dEta = etaD - etaC;
        const scalar gd = 2*(d & gradC);

        // r = 2 d.grad/(etaD - etaC) - 1 (OpenFOAM's limitedScheme form)
        scalar r;
        if (mag(dEta) > SMALL*max(mag(gd), SMALL))
        {
            r = gd/dEta - 1;
        }
        else
        {
            r = (gd*dEta >= 0 ? 2 : -1);
        }

        return etaC + 0.5*limiter(r)*dEta;
    };

    auto factor = [&]
    (
        const face& fc,
        const scalar phiF,
        const label colO,
        const scalar etaO,
        const scalar zbO,
        const vector& gradO,
        const vector& cO,
        const scalar etaN,
        const scalar zbN,
        const vector& gradN,
        const vector& cN
    )
    {
        if (min(etaO - zbO, etaN - zbN) >= faceDepthBelow_)
        {
            return scalar(1);
        }

        scalar Dm, zbF;
        faceGeometry(fc, Dm, zbF);

        vector d(cN - cO);
        d.z() = 0;

        const scalar Dmin = min(etaO - zbO, etaN - zbN);

        scalar etaF;
        if (phiF > 0)
        {
            etaF = etaFace(etaO, etaN, gradO, d, Dmin);
        }
        else if (phiF < 0)
        {
            etaF = etaFace(etaN, etaO, gradN, -d, Dmin);
        }
        else
        {
            etaF = 0.5*(etaO + etaN);
        }

        return max(etaF - zbF, scalar(0))/max(Dm, VSMALL);
    };

    const scalarField& phiI = phiRef.primitiveField();

    forAll(nei, facei)
    {
        if (!isSideFace_[facei])
        {
            continue;
        }

        if (!faceOpen_[facei])
        {
            f[facei] = 0;
            continue;
        }

        const label o = own[facei];
        const label n = nei[facei];
        const label co = cellCol[o];
        const label cn = cellCol[n];

        f[facei] = factor
        (
            faces[facei], phiI[facei], co,
            etaCol[co], zb[co], gradEta[o], colC[co],
            etaCol[cn], zb[cn], gradEta[n], colC[cn]
        );
    }

    // Coupled patches: the neighbour column's surface, bottom, centre and
    // surface gradient
    {
        volVectorField gradV
        (
            IOobject
            (
                "wdGradEta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedVector(dimless, Zero),
            calculatedFvPatchVectorField::typeName
        );
        gradV.primitiveFieldRef() = gradEta;
        gradV.correctBoundaryConditions();

        const List<scalarField> nbZb(patchNeighbourColumnValues(zb));

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!fvp.coupled())
            {
                continue;
            }

            const labelUList& fc = fvp.faceCells();
            const label start = fvp.start();
            const scalarField& phip = phiRef.boundaryField()[patchi];
            const vectorField gradNb
            (
                gradV.boundaryField()[patchi].patchNeighbourField()
            );

            // Owner-to-neighbour cell-centre vector, transformed across
            // cyclic patches
            const vectorField delta(fvp.delta());

            forAll(fc, i)
            {
                if (!patchFaceOpen_[patchi][i])
                {
                    fp[patchi][i] = 0;
                    continue;
                }

                const label o = fc[i];
                const label co = cellCol[o];

                fp[patchi][i] = factor
                (
                    faces[start + i], phip[i], co,
                    etaCol[co], zb[co], gradEta[o], colC[co],
                    nbEta[patchi][i], nbZb[patchi][i], gradNb[i],
                    colC[co] + delta[i]
                );
            }
        }
    }
}


void Foam::wettingDrying::limitOpenedFaces
(
    surfaceScalarField& phiHbyA,
    const boolList& changed,
    const scalarField& etaCol,
    const scalar g
) const
{
    if (!active_)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();
    const scalarField& zb = columns_.columnBottomZ();
    const scalarField& magSf = mesh_.magSf();
    scalarField& phiI = phiHbyA.primitiveFieldRef();

    // Dam-break velocity at the dam of a wet column, for a newly wet
    // column lower than it (a puddle left on a slope drains downhill; up a
    // slope the flow has its own momentum, and a kick there would make a
    // film creep uphill)
    auto uDam = [&](const scalar D, const scalar zbW, const scalar zbD)
    {
        if (zbD >= zbW - SMALL)
        {
            return scalar(0);
        }
        return (2.0/3.0)*Foam::sqrt(g*max(D - minDepth_, scalar(0)));
    };

    forAll(nei, facei)
    {
        if (!isSideFace_[facei] || !faceOpen_[facei])
        {
            continue;
        }

        const label co = cellCol[own[facei]];
        const label cn = cellCol[nei[facei]];
        const bool chO = changed[co];
        const bool chN = changed[cn];

        if (chO && !chN)
        {
            // Owner newly wet: flux from the neighbour only (phi < 0)
            phiI[facei] = min
            (
                min(phiI[facei], scalar(0)),
               -uDam(etaCol[cn] - zb[cn], zb[cn], zb[co])*magSf[facei]
            );
        }
        else if (chN && !chO)
        {
            phiI[facei] = max
            (
                max(phiI[facei], scalar(0)),
                uDam(etaCol[co] - zb[co], zb[co], zb[cn])*magSf[facei]
            );
        }
    }

    scalarField chF(changed.size(), 0);
    forAll(changed, coli)
    {
        chF[coli] = changed[coli] ? 1 : 0;
    }
    const List<scalarField> nbCh(patchNeighbourColumnValues(chF));
    const List<scalarField> nbD(patchNeighbourColumnValues(etaCol - zb));
    const List<scalarField> nbZb(patchNeighbourColumnValues(zb));

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (!fvp.coupled())
        {
            continue;
        }

        const labelUList& fc = fvp.faceCells();
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        scalarField& phip = phiHbyA.boundaryFieldRef()[patchi];

        forAll(fc, i)
        {
            if (!patchFaceOpen_[patchi][i])
            {
                continue;
            }

            const label co = cellCol[fc[i]];
            const bool chO = changed[co];
            const bool chN = nbCh[patchi][i] > 0.5;

            if (chO && !chN)
            {
                phip[i] = min
                (
                    min(phip[i], scalar(0)),
                   -uDam(nbD[patchi][i], nbZb[patchi][i], zb[co])*magSfp[i]
                );
            }
            else if (chN && !chO)
            {
                phip[i] = max
                (
                    max(phip[i], scalar(0)),
                    uDam(etaCol[co] - zb[co], zb[co], nbZb[patchi][i])*magSfp[i]
                );
            }
        }
    }
}


void Foam::wettingDrying::upwindSupercritical
(
    surfaceScalarField& phiHbyA,
    const volVectorField& HbyA,
    const volVectorField& U,
    const scalarField& etaCol,
    const scalar g
) const
{
    if (!active_ || upwindFroude_ <= 0)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellCol = columns_.cellColumn();
    const labelListList& colCells = columns_.columnCells();
    const scalarField& zb = columns_.columnBottomZ();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();
    const scalarField& V = mesh_.V();
    scalarField& phiI = phiHbyA.primitiveFieldRef();

    // Depth-averaged velocity and wave speed of each column
    vectorField colU(columns_.nColumns(), Zero);
    scalarField colC(columns_.nColumns(), 0);

    forAll(colCells, coli)
    {
        scalar sV = 0;
        for (const label c : colCells[coli])
        {
            colU[coli] += V[c]*U[c];
            sV += V[c];
        }
        colU[coli] /= max(sV, VSMALL);
        colC[coli] = Foam::sqrt(g*max(etaCol[coli] - zb[coli], scalar(0)));
    }

    // Weight of the upwind flux: 1 above the Froude number upwindFroude
    // (HLL takes the upwind flux once both wave speeds have the sign of the
    // flow).  Blending it in below that, partly upwind fluxes in subcritical
    // flow, lets grid-scale oscillations grow (the interpolated flux is what
    // couples the pressure and the velocity on the collocated mesh); a
    // blend above it lowered the run-up of the solitary wave by 10%.
    auto weight = [&](const scalar un, const scalar c)
    {
        return (un > upwindFroude_*c ? scalar(1) : scalar(0));
    };

    forAll(nei, facei)
    {
        if (!isSideFace_[facei] || !faceOpen_[facei])
        {
            continue;
        }

        const bool fromOwner = (phiI[facei] >= 0);
        const label donor = (fromOwner ? own[facei] : nei[facei]);
        const label cd = cellCol[donor];
        const scalar un = (colU[cd] & Sf[facei])/magSf[facei];
        const scalar w = weight(fromOwner ? un : -un, colC[cd]);

        if (w > 0)
        {
            phiI[facei] += w*((HbyA[donor] & Sf[facei]) - phiI[facei]);
        }
    }

    // Coupled patches: the donor may be the cell across the patch
    {
        const List<scalarField> nbUx(patchNeighbourColumnValues(colU.component(0)));
        const List<scalarField> nbUy(patchNeighbourColumnValues(colU.component(1)));
        const List<scalarField> nbC(patchNeighbourColumnValues(colC));

        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];

            if (!fvp.coupled())
            {
                continue;
            }

            const labelUList& fc = fvp.faceCells();
            const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
            const vectorField HbyAnb
            (
                HbyA.boundaryField()[patchi].patchNeighbourField()
            );
            scalarField& phip = phiHbyA.boundaryFieldRef()[patchi];

            forAll(fc, i)
            {
                if (!patchFaceOpen_[patchi][i])
                {
                    continue;
                }

                if (phip[i] >= 0)
                {
                    const label cd = cellCol[fc[i]];
                    const scalar un = (colU[cd] & Sfp[i])/magSfp[i];
                    const scalar w = weight(un, colC[cd]);

                    phip[i] += w*((HbyA[fc[i]] & Sfp[i]) - phip[i]);
                }
                else
                {
                    const vector Unb(nbUx[patchi][i], nbUy[patchi][i], 0);
                    const scalar un = (Unb & Sfp[i])/magSfp[i];
                    const scalar w = weight(-un, nbC[patchi][i]);

                    phip[i] += w*((HbyAnb[i] & Sfp[i]) - phip[i]);
                }
            }
        }
    }
}


// ************************************************************************* //
