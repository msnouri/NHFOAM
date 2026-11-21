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

Notation (one column, layers k = 0 at the surface to N-1 at the bottom)
    p_k     pI of cell k: pressure on the bottom interface of layer k
    p_s     surface pressure g eta^{n+1} (p_{-1}), unknown nCells + column
    pc_k    pressure at the centre level of layer k, (p_{k-1} + p_k)/2
    F_k     vertical flux U_k.S_k at the centre of layer k, with S_k the
            upward mid-layer area vector:
                F_k = H_k - B_k (p_{k-1} - p_k) + L_k
            H_k = HbyA_k.S_k, B_k = rAU_k A^2/V_k (compact vertical
            gradient), L_k = -rAU_k gradP_h.S_k,h (slope terms, lagged)
    phi_f   flux of a side face f of layer k (P -> N):
                phi_f = phiHbyA_f - c_f (pc_k(N) - pc_k(P)) - corr_f
    D_k     control volume of interface k (centre of layer k to centre of
            layer k+1, or to the bottom):
                F_k - F_{k+1} + phi(lower half of layer k)
                              + phi(upper half of layer k+1) = 0
    S       top half layer (surface to the centre of layer 0):
                phi_s - F_0 + phi(upper half of layer 0) = 0,
                phi_s = A/aEff ((1 + s) p_s/g - etaHatEff)
    A half face carries half of phi_f, so a side face of layer k adds
    c/4 (e e^T), e = (+1, +1, -1, -1) on (p_k(P), p_{k-1}(P), p_k(N),
    p_{k-1}(N)), to the rows of its two lower and two upper unknowns; F_k
    adds B_k (d d^T), d = (+1, -1) on (p_k, p_{k-1}).

    Rows are assembled in outflow form, sum_j a_ij p_j + b_i = 0, and stored
    with the OpenFOAM sign convention (negative diagonal, like
    fvm::laplacian):  diag = -a_ii, offdiag = -a_ij, source = b_i,
    interface coefficients (coupled patches) = a_i,remote.

\*---------------------------------------------------------------------------*/

#include "interfacePressure.H"
#include "interfacePressureLduMesh.H"
#include "fvcSurfaceIntegrate.H"
#include "surfaceInterpolate.H"
#include "mixedFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "calculatedFvPatchFields.H"
#include "ListOps.H"
#include "SortableList.H"
#include "EdgeMap.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::interfacePressure::upwardSf(const label facei) const
{
    vector S;

    if (mesh_.isInternalFace(facei))
    {
        S = mesh_.Sf()[facei];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.whichPatch(facei);
        S = mesh_.Sf().boundaryField()[patchi][facei - pbm[patchi].start()];
    }

    return (S.z() < 0 ? -S : S);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePressure::interfacePressure
(
    const fvMesh& mesh,
    const sigmaColumns& columns
)
:
    mesh_(mesh),
    columns_(columns),
    above_(mesh.nCells(), -1),
    cellBelow_(mesh.nCells(), -1),
    cellTopFace_(mesh.nCells(), -1),
    cellBottomFace_(mesh.nCells(), -1),
    isLayerFace_(mesh.nInternalFaces(), false),
    faceAbove_(mesh.nInternalFaces(), -1),
    patchFaceAbove_(mesh.boundary().size()),
    faceToExt_(mesh.nInternalFaces(), -1),
    aboveExt_(mesh.nCells(), -1),
    pairAboveExt_(mesh.nInternalFaces(), -1),
    diagOwnNbrAbove_(mesh.nInternalFaces(), -1),
    diagOwnAboveNbr_(mesh.nInternalFaces(), -1),
    cf_(mesh.nInternalFaces(), 0),
    cfB_(mesh.boundary().size()),
    corr_(mesh.nInternalFaces(), 0),
    corrB_(mesh.boundary().size())
{
    const label nCells = mesh.nCells();
    const labelListList& colCells = columns.columnCells();
    const labelListList& colLayerFaces = columns.columnLayerFaces();
    const labelList& colBottomFace = columns.columnBottomFace();
    const label surfStart = columns.surfacePatch().start();

    forAll(colCells, coli)
    {
        const labelList& cc = colCells[coli];
        const labelList& lf = colLayerFaces[coli];

        forAll(cc, k)
        {
            const label celli = cc[k];

            if (k > 0)
            {
                above_[celli] = cc[k-1];
                cellTopFace_[celli] = lf[k-1];
            }
            else
            {
                above_[celli] = nCells + coli;
                cellTopFace_[celli] = surfStart + coli;
            }

            if (k < cc.size() - 1)
            {
                cellBelow_[celli] = cc[k+1];
                cellBottomFace_[celli] = lf[k];
            }
            else
            {
                cellBottomFace_[celli] = colBottomFace[coli];
            }
        }

        forAll(lf, i)
        {
            isLayerFace_[lf[i]] = true;
        }
    }

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const cellList& cells = mesh.cells();

    // Internal side faces: the side face between the two cells above
    forAll(nei, facei)
    {
        if (isLayerFace_[facei])
        {
            continue;
        }

        const label oa = above_[own[facei]];
        const label na = above_[nei[facei]];

        if ((oa < nCells) != (na < nCells))
        {
            FatalErrorInFunction
                << "Side face " << facei << " at " << mesh.Cf()[facei]
                << " joins cells of different layers." << nl
                << "All columns must have the same number of layers."
                << exit(FatalError);
        }

        if (oa >= nCells)
        {
            continue;
        }

        const cell& c = cells[oa];

        forAll(c, i)
        {
            const label fj = c[i];

            if
            (
                mesh.isInternalFace(fj)
             && (
                    (own[fj] == oa && nei[fj] == na)
                 || (own[fj] == na && nei[fj] == oa)
                )
            )
            {
                faceAbove_[facei] = fj;
                break;
            }
        }

        if (faceAbove_[facei] < 0)
        {
            FatalErrorInFunction
                << "No side face above face " << facei << " at "
                << mesh.Cf()[facei] << exit(FatalError);
        }
    }

    // Patches: coefficient storage; coupled side patches: face above
    const faceList& faces = mesh.faces();
    const label surfID = columns.surfacePatchID();
    const label botID = columns.bottomPatchID();

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];

        cfB_.set(patchi, new scalarField(fvp.size(), 0));
        corrB_.set(patchi, new scalarField(fvp.size(), 0));
        patchFaceAbove_[patchi].setSize(fvp.size(), -1);

        if (patchi == surfID || patchi == botID || !fvp.coupled())
        {
            continue;
        }

        const label start = fvp.start();
        const labelUList& fc = fvp.faceCells();
        const labelListList cellPatchFaces(invertOneToMany(nCells, fc));
        labelList& aboveFace = patchFaceAbove_[patchi];

        forAll(fc, i)
        {
            const label ca = above_[fc[i]];

            if (ca >= nCells)
            {
                continue;
            }

            const face& fi = faces[start + i];
            const labelList& cand = cellPatchFaces[ca];

            forAll(cand, j)
            {
                const face& fj = faces[start + cand[j]];

                label nShared = 0;
                forAll(fj, fpi)
                {
                    if (fi.found(fj[fpi]))
                    {
                        ++nShared;
                    }
                }

                if (nShared >= 2)
                {
                    aboveFace[i] = cand[j];
                    break;
                }
            }

            if (aboveFace[i] < 0)
            {
                FatalErrorInFunction
                    << "No face above face " << i << " of coupled patch "
                    << fvp.name() << " at " << fvp.Cf()[i]
                    << exit(FatalError);
            }
        }
    }


    // Extended addressing: unknowns are the cells and the surfaces of the
    // columns; couplings are the mesh faces, every unknown with the one
    // above, and for every side face the unknowns above (the two diagonals
    // and the pair above)
    {
        EdgeMap<label> pairIndex(6*nei.size());
        DynamicList<edge> pairs(6*nei.size());

        auto addPair = [&](const label a, const label b) -> label
        {
            const edge e(a, b);
            auto iter = pairIndex.cfind(e);

            if (iter.found())
            {
                return *iter;
            }

            pairIndex.insert(e, pairs.size());
            pairs.append(e);
            return pairs.size() - 1;
        };

        labelList facePair(nei.size(), -1);
        labelList abovePair(nCells, -1);
        labelList pairAbove(nei.size(), -1);
        labelList diagA(nei.size(), -1);
        labelList diagB(nei.size(), -1);

        forAll(nei, facei)
        {
            facePair[facei] = addPair(own[facei], nei[facei]);
        }

        forAll(above_, celli)
        {
            abovePair[celli] = addPair(celli, above_[celli]);
        }

        forAll(nei, facei)
        {
            if (!isLayerFace_[facei])
            {
                const label o = own[facei];
                const label n = nei[facei];

                pairAbove[facei] = addPair(above_[o], above_[n]);
                diagA[facei] = addPair(o, above_[n]);
                diagB[facei] = addPair(above_[o], n);
            }
        }

        // Sort by lower then upper unknown (upper-triangular order)
        const label nUnknowns = nCells + columns.nColumns();
        SortableList<scalar> order(pairs.size());

        forAll(pairs, i)
        {
            order[i] =
                scalar(pairs[i].minVertex())*scalar(nUnknowns)
              + scalar(pairs[i].maxVertex());
        }
        order.sort();

        const labelList& idx = order.indices();
        labelList newIndex(pairs.size());

        extLower_.setSize(pairs.size());
        extUpper_.setSize(pairs.size());

        forAll(idx, i)
        {
            newIndex[idx[i]] = i;
            extLower_[i] = pairs[idx[i]].minVertex();
            extUpper_[i] = pairs[idx[i]].maxVertex();
        }

        forAll(nei, facei)
        {
            faceToExt_[facei] = newIndex[facePair[facei]];

            if (pairAbove[facei] >= 0)
            {
                pairAboveExt_[facei] = newIndex[pairAbove[facei]];
                diagOwnNbrAbove_[facei] = newIndex[diagA[facei]];
                diagOwnAboveNbr_[facei] = newIndex[diagB[facei]];
            }
        }

        forAll(above_, celli)
        {
            aboveExt_[celli] = newIndex[abovePair[celli]];
        }

        labelList lAddr(extLower_);
        labelList uAddr(extUpper_);

        // Vertical lines of unknowns (surface first) and their couplings,
        // for the sigmaColumn preconditioner
        labelListList colUnknowns(colCells.size());
        labelListList colCouplings(colCells.size());

        forAll(colCells, coli)
        {
            const labelList& cc = colCells[coli];
            labelList& unk = colUnknowns[coli];
            labelList& cpl = colCouplings[coli];

            unk.setSize(cc.size() + 1);
            cpl.setSize(cc.size());

            unk[0] = nCells + coli;

            forAll(cc, k)
            {
                unk[k+1] = cc[k];
                cpl[k] = newIndex[abovePair[cc[k]]];
            }
        }

        extMeshPtr_.reset
        (
            new interfacePressureLduMesh
            (
                mesh,
                nUnknowns,
                lAddr,
                uAddr,
                std::move(colUnknowns),
                std::move(colCouplings)
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::interfacePressure::solve
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
)
{
    const label nCells = mesh_.nCells();
    const label nCols = columns_.nColumns();
    const label nUnknowns = nCells + nCols;
    const label nPatches = mesh_.boundary().size();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& V = mesh_.V();
    const scalarField& magSf = mesh_.magSf();
    const scalarField& dc = mesh_.nonOrthDeltaCoeffs();
    const vectorField& cv = mesh_.nonOrthCorrectionVectors();
    const scalarField& colArea = columns_.columnArea();
    const labelList& cellCol = columns_.cellColumn();
    const label surfID = columns_.surfacePatchID();
    const label botID = columns_.bottomPatchID();

    // Current values of all unknowns (initial guess and lagged state)
    scalarField x(nUnknowns);
    SubList<scalar>(x, nCells) = pI.primitiveField();
    SubList<scalar>(x, nCols, nCells) = ps;

    // Gradient at the faces for the non-orthogonal correction
    const surfaceVectorField gradPf(fvc::interpolate(gradP));

    // Value of the unknown above each cell, for the couplings across
    // processor boundaries that are lagged
    volScalarField pAbove
    (
        IOobject
        (
            "pAboveI",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(pI.dimensions(), Zero),
        calculatedFvPatchScalarField::typeName
    );

    forAll(pAbove, celli)
    {
        pAbove[celli] = x[above_[celli]];
    }
    pAbove.correctBoundaryConditions();

    // Outflow form: sum_j a_ij p_j + b_i = 0
    const label nExt = extLower_.size();
    scalarField aDiag(nUnknowns, 0);
    scalarField aUp(nExt, 0);           // row extLower, column extUpper
    scalarField aLo(nExt, 0);           // row extUpper, column extLower
    scalarField b(nUnknowns, 0);
    FieldField<Field, scalar> aBou(nPatches);

    forAll(mesh_.boundary(), patchi)
    {
        aBou.set
        (
            patchi,
            new scalarField(mesh_.boundary()[patchi].size(), 0)
        );
    }

    auto addA = [&](const label rowi, const label exti, const scalar a)
    {
        if (rowi == extLower_[exti])
        {
            aUp[exti] += a;
        }
        else
        {
            aLo[exti] += a;
        }
    };


    // --- Side faces, internal: c/4 (e e^T) on (o, oa, n, na)

    forAll(nei, facei)
    {
        if (isLayerFace_[facei])
        {
            continue;
        }

        const label o = own[facei];
        const label n = nei[facei];
        const label oa = above_[o];
        const label na = above_[n];

        const scalar c = rAUf[facei]*magSf[facei]*dc[facei];
        const scalar corr =
            rAUf[facei]*magSf[facei]*(cv[facei] & gradPf[facei]);
        const scalar c4 = 0.25*c;
        const scalar Xh = 0.5*(phiHbyA[facei] - corr);

        cf_[facei] = c;
        corr_[facei] = corr;

        const label eOn = faceToExt_[facei];
        const label eOA = aboveExt_[o];
        const label eNA = aboveExt_[n];
        const label eAA = pairAboveExt_[facei];
        const label eOna = diagOwnNbrAbove_[facei];
        const label eOAn = diagOwnAboveNbr_[facei];

        // Lower half: rows o and n
        aDiag[o] += c4;
        addA(o, eOA, c4);
        addA(o, eOn, -c4);
        addA(o, eOna, -c4);
        b[o] += Xh;

        aDiag[n] += c4;
        addA(n, eNA, c4);
        addA(n, eOn, -c4);
        addA(n, eOAn, -c4);
        b[n] -= Xh;

        // Upper half: rows oa and na (cells above, or the surfaces)
        aDiag[oa] += c4;
        addA(oa, eOA, c4);
        addA(oa, eAA, -c4);
        addA(oa, eOAn, -c4);
        b[oa] += Xh;

        aDiag[na] += c4;
        addA(na, eNA, c4);
        addA(na, eAA, -c4);
        addA(na, eOna, -c4);
        b[na] -= Xh;
    }


    // --- Side faces on patches

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (patchi == surfID || patchi == botID)
        {
            continue;
        }

        const labelUList& fc = fvp.faceCells();

        if (fvp.coupled())
        {
            const scalarField& rAUfp = rAUf.boundaryField()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
            const scalarField& dcp =
                mesh_.nonOrthDeltaCoeffs().boundaryField()[patchi];
            const vectorField& cvp =
                mesh_.nonOrthCorrectionVectors().boundaryField()[patchi];
            const vectorField& gPf = gradPf.boundaryField()[patchi];
            const scalarField& phiHp = phiHbyA.boundaryField()[patchi];
            const scalarField pInb
            (
                pI.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField pAnb
            (
                pAbove.boundaryField()[patchi].patchNeighbourField()
            );
            const labelList& aboveFace = patchFaceAbove_[patchi];

            scalarField& cfp = cfB_[patchi];
            scalarField& corrp = corrB_[patchi];
            scalarField& aBp = aBou[patchi];

            forAll(fc, i)
            {
                const label o = fc[i];
                const label oa = above_[o];

                const scalar c = rAUfp[i]*magSfp[i]*dcp[i];
                const scalar corr = rAUfp[i]*magSfp[i]*(cvp[i] & gPf[i]);
                const scalar c4 = 0.25*c;
                const scalar Xh = 0.5*(phiHp[i] - corr);

                cfp[i] = c;
                corrp[i] = corr;

                // Processor faces: the remote unknowns above the cells
                // cannot be coupled implicitly.  Dropping only those
                // couplings would make the matrix indefinite, so the face
                // is written in compact form, c/2 ((o - n)^2 + (oa - na)^2),
                // with the difference to c/4 (o + oa - n - na)^2 lagged.
                // The compact form bounds the exact one, so the lagged
                // iteration converges.  At the surface (oa, na not cells)
                // c/2 oa^2 stays on the diagonal and na is lagged.
                const scalar c2 = 0.5*c;
                const scalar xo = x[o];
                const scalar xoa = x[oa];
                const scalar xn = pInb[i];
                const scalar xna = pAnb[i];

                // Lower half: row o
                aDiag[o] += c2;
                aBp[i] -= c2;
                b[o] += Xh + c4*((xoa - xo) - (xna - xn));

                // Upper half: row oa
                aDiag[oa] += c2;

                if (aboveFace[i] >= 0)
                {
                    aBp[aboveFace[i]] -= c2;
                    b[oa] += Xh + c4*((xo - xoa) - (xn - xna));
                }
                else
                {
                    b[oa] += Xh + c4*(xo - xoa - xna - xn);
                }
            }
        }
        else
        {
            // Walls and other flux boundaries: the boundary velocity flux
            // (fixedFluxPressure semantics)
            const vectorField& Up = U.boundaryField()[patchi];
            const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];

            forAll(fc, i)
            {
                const label o = fc[i];
                const scalar halfF = 0.5*(Up[i] & Sfp[i]);

                b[o] += halfF;
                b[above_[o]] += halfF;
            }
        }
    }


    // --- Vertical fluxes at the cell centres: B_k (d d^T) on (p_k, above)

    const vectorField& Ubot = U.boundaryField()[botID];
    const vectorField& Sbot = mesh_.Sf().boundaryField()[botID];
    const label botStart = mesh_.boundaryMesh()[botID].start();

    forAll(above_, celli)
    {
        const label col = cellCol[celli];
        const label ca = above_[celli];
        const label e = aboveExt_[celli];

        const vector Smid =
            0.5
           *(
                upwardSf(cellTopFace_[celli])
              + upwardSf(cellBottomFace_[celli])
            );

        // Over-relaxed split of the flux through the sloping mid-plane,
        // as OpenFOAM's corrected Laplacian: the implicit vertical
        // coefficient carries |S_mid|^2 = A^2 + |S_h|^2, and the lagged part
        // keeps only the gradient along the layer.  With A^2 alone, the
        // lagged part would contain -|S_h|^2/A dp/dz, and the corrector
        // iteration can diverge where the layers are steep (it did under
        // bore fronts).  The converged flux is the same.
        const scalar A = colArea[col];
        const scalar Sh2 = sqr(Smid.x()) + sqr(Smid.y());

        const scalar B = rAU[celli]*(factor[col]*sqr(A) + Sh2)/V[celli];

        const scalar F0 =
            (HbyA[celli] & Smid)
          - rAU[celli]
           *(
                gradP[celli].x()*Smid.x() + gradP[celli].y()*Smid.y()
              - gradP[celli].z()*Sh2/A
            );

        // + F_k in the row of this cell, - F_k in the row above
        aDiag[celli] += B;
        aDiag[ca] += B;
        addA(celli, e, -B);
        addA(ca, e, -B);
        b[celli] += F0;
        b[ca] -= F0;

        // Flux through the bottom
        if (cellBelow_[celli] < 0)
        {
            const label i = cellBottomFace_[celli] - botStart;
            b[celli] += Ubot[i] & Sbot[i];
        }
    }


    // --- Free surface: phi_s = A/aEff ((1 + s) p_s/g - etaHatEff)

    forAll(colArea, col)
    {
        const label s = nCells + col;

        aDiag[s] += colArea[col]*colS[col]/(g*aEff);
        b[s] -= colArea[col]*etaHatEff[col]/aEff;
    }


    // --- Solve (symmetric: only the upper triangle is stored)

    lduMatrix M(extMeshPtr_());
    M.upper() = -0.5*(aUp + aLo);
    M.diag() = -aDiag;

    FieldField<Field, scalar> intCoeffs(nPatches);
    forAll(mesh_.boundary(), patchi)
    {
        intCoeffs.set
        (
            patchi,
            new scalarField(mesh_.boundary()[patchi].size(), 0)
        );
    }

    // PCG (or PBiCGStab) with the sigmaColumn or diagonal preconditioner:
    // incomplete factorisations (DIC, DILU) and GAMG are replaced, since
    // the rows are not diagonally dominant where layers are thicker than
    // the horizontal spacing
    dictionary solverControls(controls);
    {
        const word solverName(solverControls.get<word>("solver"));

        if (solverName != "PCG" && solverName != "PBiCGStab")
        {
            solverControls.set("solver", word("PCG"));
        }

        word precon("sigmaColumn");

        if (solverControls.found("preconditioner"))
        {
            const entry& e = solverControls.lookupEntry
            (
                "preconditioner",
                keyType::LITERAL
            );

            if (!e.isDict())
            {
                const word name(e.stream());

                if (name == "diagonal" || name == "none")
                {
                    precon = name;
                }
            }
        }

        solverControls.set("preconditioner", precon);
    }

    solverPerformance solverPerf =
        lduMatrix::solver::New
        (
            pI.name(),
            M,
            aBou,
            intCoeffs,
            pI.boundaryField().scalarInterfaces(),
            solverControls
        )->solve(x, b);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(mesh_.comm()));
    }

    mesh_.setSolverPerformance(pI.name(), solverPerf);

    pI.primitiveFieldRef() = SubList<scalar>(x, nCells);
    ps = SubList<scalar>(x, nCols, nCells);
    pI.correctBoundaryConditions();

    return solverPerf;
}


void Foam::interfacePressure::cellPressure
(
    const volScalarField& pI,
    const scalarField& ps,
    volScalarField& pCell
) const
{
    const label nCells = mesh_.nCells();
    scalarField& pc = pCell.primitiveFieldRef();

    forAll(pc, celli)
    {
        const label ca = above_[celli];
        const scalar top = (ca < nCells ? pI[ca] : ps[ca - nCells]);

        pc[celli] = 0.5*(top + pI[celli]);
    }

    // Free surface: p_rgh = p_s
    fvPatchScalarField& pSurf =
        pCell.boundaryFieldRef()[columns_.surfacePatchID()];

    if (isA<mixedFvPatchScalarField>(pSurf))
    {
        mixedFvPatchScalarField& m =
            refCast<mixedFvPatchScalarField>(pSurf);

        m.refValue() = ps;
        m.refGrad() = 0;
        m.valueFraction() = 1;
    }

    // Flux boundaries keep their gradient (set by constrainPressure); mark
    // them updated so they can be evaluated at any point of the time step
    forAll(pCell.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pCell.boundaryFieldRef()[patchi];

        if (isA<fixedFluxPressureFvPatchScalarField>(pp))
        {
            fixedFluxPressureFvPatchScalarField& ffp =
                refCast<fixedFluxPressureFvPatchScalarField>(pp);

            ffp.updateSnGrad(ffp.gradient());
        }
    }

    pCell.correctBoundaryConditions();
}


Foam::tmp<Foam::volVectorField> Foam::interfacePressure::grad
(
    const volScalarField& pI,
    const scalarField& ps,
    const volScalarField& pCell
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelList& cellLay = columns_.cellLayer();
    const label surfID = columns_.surfacePatchID();
    const label botID = columns_.bottomPatchID();

    // Side faces: the pressure at the centre level, interpolated
    surfaceScalarField pf(fvc::interpolate(pCell));

    // Layer faces: the interface value (bottom interface of the upper cell)
    scalarField& pfi = pf.primitiveFieldRef();

    forAll(nei, facei)
    {
        if (isLayerFace_[facei])
        {
            const label o = own[facei];
            const label n = nei[facei];
            pfi[facei] = (cellLay[o] < cellLay[n] ? pI[o] : pI[n]);
        }
    }

    // Surface and bottom
    surfaceScalarField::Boundary& pfB = pf.boundaryFieldRef();
    pfB[surfID] = ps;

    {
        scalarField& pb = pfB[botID];
        const labelUList& fc = mesh_.boundary()[botID].faceCells();

        forAll(pb, i)
        {
            pb[i] = pI[fc[i]];
        }
    }

    return fvc::surfaceIntegrate(mesh_.Sf()*pf);
}


Foam::tmp<Foam::scalarField> Foam::interfacePressure::fluxes
(
    const volScalarField& pCell,
    const surfaceScalarField& phiHbyA,
    const volVectorField& U,
    surfaceScalarField& phi
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelListList& colCells = columns_.columnCells();
    const label surfID = columns_.surfacePatchID();
    const label botID = columns_.bottomPatchID();
    const label botStart = mesh_.boundaryMesh()[botID].start();

    // Outward side flux of every cell
    scalarField sideOut(mesh_.nCells(), 0);

    scalarField& phii = phi.primitiveFieldRef();

    forAll(nei, facei)
    {
        if (isLayerFace_[facei])
        {
            continue;
        }

        const label o = own[facei];
        const label n = nei[facei];

        const scalar F =
            phiHbyA[facei]
          - cf_[facei]*(pCell[n] - pCell[o])
          - corr_[facei];

        phii[facei] = F;
        sideOut[o] += F;
        sideOut[n] -= F;
    }

    surfaceScalarField::Boundary& phiB = phi.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];

        if (patchi == surfID)
        {
            continue;
        }

        const labelUList& fc = fvp.faceCells();
        scalarField& phip = phiB[patchi];

        if (patchi != botID && fvp.coupled())
        {
            const scalarField& phiHp = phiHbyA.boundaryField()[patchi];
            const scalarField pCnb
            (
                pCell.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& cfp = cfB_[patchi];
            const scalarField& corrp = corrB_[patchi];

            forAll(fc, i)
            {
                const label o = fc[i];

                phip[i] = phiHp[i] - cfp[i]*(pCnb[i] - pCell[o]) - corrp[i];
                sideOut[o] += phip[i];
            }
        }
        else
        {
            const vectorField& Up = U.boundaryField()[patchi];
            const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];

            forAll(fc, i)
            {
                phip[i] = Up[i] & Sfp[i];

                if (patchi != botID)
                {
                    sideOut[fc[i]] += phip[i];
                }
            }
        }
    }

    // Layer and surface fluxes: continuity integrated from the bottom
    tmp<scalarField> tphis(new scalarField(columns_.nColumns(), 0));
    scalarField& phis = tphis.ref();
    const scalarField& phiBot = phiB[botID];

    forAll(colCells, col)
    {
        const labelList& cc = colCells[col];
        const label nl = cc.size();

        // Upward flux into the bottom cell
        scalar up = -phiBot[cellBottomFace_[cc[nl-1]] - botStart];

        for (label k = nl - 1; k >= 0; --k)
        {
            const label c = cc[k];

            up -= sideOut[c];

            if (k > 0)
            {
                const label f = cellTopFace_[c];
                phii[f] = (own[f] == c ? up : -up);
            }
            else
            {
                phis[col] = up;
            }
        }
    }

    phiB[surfID] = phis;

    return tphis;
}


// ************************************************************************* //
