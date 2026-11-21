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

#include "sigmaColumns.H"
#include "syncTools.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
    inline Foam::scalar horizontalDistance
    (
        const Foam::vector& a,
        const Foam::vector& b
    )
    {
        return Foam::sqrt(Foam::sqr(a.x() - b.x()) + Foam::sqr(a.y() - b.y()));
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigmaColumns::findColumns()
{
    const polyPatch& sp = surfacePatch();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    cellColumn_.setSize(mesh_.nCells(), -1);
    cellLayer_.setSize(mesh_.nCells(), -1);
    columnCells_.setSize(sp.size());
    columnBottomFace_.setSize(sp.size(), -1);
    columnLayerFaces_.setSize(sp.size());

    forAll(sp, coli)
    {
        label facei = sp.start() + coli;
        label celli = own[facei];
        label layer = 0;

        DynamicList<label> colCells(16);
        DynamicList<label> layerFaces(16);

        while (true)
        {
            if (cellColumn_[celli] != -1)
            {
                FatalErrorInFunction
                    << "Cell " << celli << " belongs to columns "
                    << cellColumn_[celli] << " and " << coli << nl
                    << "The mesh is not a vertical extrusion."
                    << exit(FatalError);
            }

            cellColumn_[celli] = coli;
            cellLayer_[celli] = layer++;
            colCells.append(celli);

            // Face of the cell opposite to the face we came through:
            // the one that shares no point with it
            const face& enterFace = faces[facei];
            const cell& c = cells[celli];

            label oppFacei = -1;

            forAll(c, cfi)
            {
                const label fj = c[cfi];

                if (fj == facei)
                {
                    continue;
                }

                const face& fc = faces[fj];

                bool shared = false;
                forAll(fc, fpi)
                {
                    if (enterFace.found(fc[fpi]))
                    {
                        shared = true;
                        break;
                    }
                }

                if (!shared)
                {
                    oppFacei = fj;
                    break;
                }
            }

            if (oppFacei < 0)
            {
                FatalErrorInFunction
                    << "Cell " << celli << " at " << mesh_.cellCentres()[celli]
                    << " has no face opposite to face " << facei << nl
                    << "The mesh is not a vertical extrusion of a 2D mesh."
                    << exit(FatalError);
            }

            if (mesh_.isInternalFace(oppFacei))
            {
                layerFaces.append(oppFacei);
                celli =
                (
                    own[oppFacei] == celli
                  ? nei[oppFacei]
                  : own[oppFacei]
                );
                facei = oppFacei;
            }
            else
            {
                const label patchi = pbm.whichPatch(oppFacei);

                if (patchi != bottomPatchID_)
                {
                    FatalErrorInFunction
                        << "The column below surface face " << coli
                        << " (at " << sp.faceCentres()[coli] << ")"
                        << " ends on patch " << pbm[patchi].name()
                        << " instead of the bottom patch "
                        << bottomPatchName_ << nl
                        << "In parallel, keep columns on one processor with"
                        << " the 'structured' decomposition method."
                        << exit(FatalError);
                }

                columnBottomFace_[coli] = oppFacei;
                break;
            }
        }

        columnCells_[coli].transfer(colCells);
        columnLayerFaces_[coli].transfer(layerFaces);
    }

    label nMissing = 0;
    forAll(cellColumn_, celli)
    {
        if (cellColumn_[celli] < 0)
        {
            ++nMissing;
        }
    }
    reduce(nMissing, sumOp<label>());

    if (nMissing)
    {
        FatalErrorInFunction
            << nMissing << " cells do not belong to any column below patch "
            << surfacePatchName_ << nl
            << "The mesh must be a vertical extrusion between the surface"
            << " and the bottom patch." << exit(FatalError);
    }
}


void Foam::sigmaColumns::findLines()
{
    const polyPatch& sp = surfacePatch();
    const polyPatch& bp = mesh_.boundaryMesh()[bottomPatchID_];
    const labelList& spMeshPoints = sp.meshPoints();
    const pointField& pts = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const labelListList& pEdges = mesh_.pointEdges();

    boolList isBottomPoint(mesh_.nPoints(), false);
    {
        const labelList& bpMeshPoints = bp.meshPoints();
        forAll(bpMeshPoints, i)
        {
            isBottomPoint[bpMeshPoints[i]] = true;
        }
    }

    lineMeshPoints_.setSize(spMeshPoints.size());
    lineBottomZ_.setSize(spMeshPoints.size(), 0);
    lineCentre_.setSize(spMeshPoints.size(), Zero);
    pointLine_.setSize(mesh_.nPoints(), -1);
    pointSigma_.setSize(mesh_.nPoints(), 0);

    forAll(spMeshPoints, linei)
    {
        DynamicList<label> line(16);

        label pointi = spMeshPoints[linei];
        line.append(pointi);

        label guard = 0;

        while (!isBottomPoint[pointi])
        {
            // Neighbouring point straight below
            label below = -1;
            scalar bestRatio = GREAT;

            const labelList& pe = pEdges[pointi];

            forAll(pe, pei)
            {
                const label pointj = edges[pe[pei]].otherVertex(pointi);
                const vector d = pts[pointj] - pts[pointi];

                if (d.z() < 0)
                {
                    const scalar ratio =
                        Foam::sqrt(sqr(d.x()) + sqr(d.y()))/(-d.z());

                    if (ratio < bestRatio)
                    {
                        bestRatio = ratio;
                        below = pointj;
                    }
                }
            }

            if (below < 0 || bestRatio > 1e-6)
            {
                FatalErrorInFunction
                    << "Point " << pointi << " at " << pts[pointi]
                    << " has no point vertically below it." << nl
                    << "Mesh points must lie on vertical lines between the"
                    << " surface and the bottom patch." << exit(FatalError);
            }

            pointi = below;
            line.append(pointi);

            if (++guard > mesh_.nPoints())
            {
                FatalErrorInFunction
                    << "Endless vertical line from surface point "
                    << spMeshPoints[linei] << exit(FatalError);
            }
        }

        lineMeshPoints_[linei].transfer(line);
    }

    forAll(lineMeshPoints_, linei)
    {
        const labelList& line = lineMeshPoints_[linei];
        const scalar zt = pts[line.first()].z();
        const scalar zb = pts[line.last()].z();

        if (zt - zb < VSMALL)
        {
            FatalErrorInFunction
                << "Zero water depth on the line below point "
                << pts[line.first()] << exit(FatalError);
        }

        lineBottomZ_[linei] = zb;
        lineCentre_[linei] =
            vector(pts[line.first()].x(), pts[line.first()].y(), 0);

        forAll(line, i)
        {
            const label pointi = line[i];

            if (pointLine_[pointi] != -1 && pointLine_[pointi] != linei)
            {
                FatalErrorInFunction
                    << "Point " << pointi << " lies on two vertical lines"
                    << exit(FatalError);
            }

            pointLine_[pointi] = linei;
            pointSigma_[pointi] = (pts[pointi].z() - zb)/(zt - zb);
        }
    }

    label nMissing = 0;
    forAll(pointLine_, pointi)
    {
        if (pointLine_[pointi] < 0)
        {
            ++nMissing;
        }
    }
    reduce(nMissing, sumOp<label>());

    if (nMissing)
    {
        FatalErrorInFunction
            << nMissing << " mesh points are not on a vertical line below the"
            << " surface patch " << surfacePatchName_ << exit(FatalError);
    }
}


void Foam::sigmaColumns::calcGeometry()
{
    const polyPatch& sp = surfacePatch();
    // Copies: polyPatch returns sub-fields of the mesh geometry
    const vectorField Sf(sp.faceAreas());
    const vectorField Cf(sp.faceCentres());
    const faceList& lf = sp.localFaces();
    const pointField& lp = sp.localPoints();

    columnArea_.setSize(sp.size());
    columnCentre_.setSize(sp.size());
    columnBottomZ_.setSize(sp.size());
    faceToPointWeights_.setSize(sp.size());

    forAll(sp, coli)
    {
        // Projected (horizontal) area: invariant under vertical point motion
        columnArea_[coli] = Sf[coli].z();

        if (columnArea_[coli] <= VSMALL)
        {
            FatalErrorInFunction
                << "Surface face " << coli << " at " << Cf[coli]
                << " does not point upwards (area vector " << Sf[coli]
                << ")" << exit(FatalError);
        }

        columnCentre_[coli] = vector(Cf[coli].x(), Cf[coli].y(), 0);
        columnBottomZ_[coli] =
            mesh_.faceCentres()[columnBottomFace_[coli]].z();

        const face& f = lf[coli];
        scalarList& w = faceToPointWeights_[coli];
        w.setSize(f.size());

        forAll(f, fpi)
        {
            w[fpi] = 1.0/max(horizontalDistance(Cf[coli], lp[f[fpi]]), VSMALL);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigmaColumns::sigmaColumns(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    surfacePatchName_(dict.getOrDefault<word>("surfacePatch", "surface")),
    bottomPatchName_(dict.getOrDefault<word>("bottomPatch", "bottom")),
    surfacePatchID_(mesh.boundaryMesh().findPatchID(surfacePatchName_)),
    bottomPatchID_(mesh.boundaryMesh().findPatchID(bottomPatchName_))
{
    if (surfacePatchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Free-surface patch " << surfacePatchName_ << " not found."
            << " Available patches: " << mesh.boundaryMesh().names()
            << exit(FatalIOError);
    }

    if (bottomPatchID_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Bottom patch " << bottomPatchName_ << " not found."
            << " Available patches: " << mesh.boundaryMesh().names()
            << exit(FatalIOError);
    }

    findColumns();
    findLines();
    calcGeometry();

    lineBottomZ0_ = lineBottomZ_;

    // Neighbouring columns across the surface-patch edges
    {
        const labelListList& ff = surfacePatch().faceFaces();
        columnNeighbours_.setSize(ff.size());
        forAll(ff, coli)
        {
            columnNeighbours_[coli] = ff[coli];
        }
    }

    label minLayers = labelMax;
    label maxLayers = 0;
    forAll(columnCells_, coli)
    {
        minLayers = min(minLayers, columnCells_[coli].size());
        maxLayers = max(maxLayers, columnCells_[coli].size());
    }
    reduce(minLayers, minOp<label>());
    reduce(maxLayers, maxOp<label>());

    Info<< "Sigma columns: " << returnReduce(nColumns(), sumOp<label>())
        << " columns below patch " << surfacePatchName_
        << ", layers per column " << minLayers;
    if (maxLayers != minLayers)
    {
        Info<< " to " << maxLayers;
    }
    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigmaColumns::setLineBottomZ(const scalarField& zb)
{
    if (zb.size() != lineBottomZ_.size())
    {
        FatalErrorInFunction
            << "Size " << zb.size() << " differs from the number of lines "
            << lineBottomZ_.size() << exit(FatalError);
    }

    lineBottomZ_ = zb;
}


void Foam::sigmaColumns::updateGeometry()
{
    const vectorField& Cf = mesh_.faceCentres();

    forAll(columnBottomFace_, coli)
    {
        columnBottomZ_[coli] = Cf[columnBottomFace_[coli]].z();
    }
}


Foam::tmp<Foam::scalarField> Foam::sigmaColumns::columnValues
(
    const scalarField& cellValues
) const
{
    tmp<scalarField> tvals(new scalarField(nColumns()));
    scalarField& vals = tvals.ref();

    forAll(columnCells_, coli)
    {
        vals[coli] = cellValues[columnCells_[coli].first()];
    }

    return tvals;
}


void Foam::sigmaColumns::setCellValues
(
    const scalarField& colValues,
    scalarField& cellValues
) const
{
    forAll(columnCells_, coli)
    {
        const labelList& cc = columnCells_[coli];

        forAll(cc, i)
        {
            cellValues[cc[i]] = colValues[coli];
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::sigmaColumns::surfacePointValues
(
    const scalarField& colValues
) const
{
    const polyPatch& sp = surfacePatch();
    const faceList& lf = sp.localFaces();
    const label nPts = sp.nPoints();

    scalarField sumW(nPts, 0);
    scalarField sumWV(nPts, 0);

    forAll(lf, coli)
    {
        const face& f = lf[coli];
        const scalarList& w = faceToPointWeights_[coli];

        forAll(f, fpi)
        {
            sumW[f[fpi]] += w[fpi];
            sumWV[f[fpi]] += w[fpi]*colValues[coli];
        }
    }

    // Points shared with other processors (or cyclics)
    syncTools::syncPointList
    (
        mesh_,
        sp.meshPoints(),
        sumW,
        plusEqOp<scalar>(),
        scalar(0)
    );
    syncTools::syncPointList
    (
        mesh_,
        sp.meshPoints(),
        sumWV,
        plusEqOp<scalar>(),
        scalar(0)
    );

    tmp<scalarField> tvals(new scalarField(nPts, 0));
    scalarField& vals = tvals.ref();

    forAll(vals, i)
    {
        vals[i] = sumWV[i]/max(sumW[i], VSMALL);
    }

    return tvals;
}


Foam::tmp<Foam::pointField> Foam::sigmaColumns::newPoints
(
    const scalarField& eta,
    const scalar minDepth
) const
{
    tmp<scalarField> tetaPoint = surfacePointValues(eta);
    const scalarField& etaPoint = tetaPoint();

    tmp<pointField> tpoints(new pointField(mesh_.points()));
    pointField& pts = tpoints.ref();

    forAll(lineMeshPoints_, linei)
    {
        const scalar zb = lineBottomZ_[linei];
        const scalar zt = max(etaPoint[linei], zb + minDepth);
        const labelList& line = lineMeshPoints_[linei];

        forAll(line, i)
        {
            const label pointi = line[i];
            pts[pointi].z() = zb + pointSigma_[pointi]*(zt - zb);
        }
    }

    return tpoints;
}


Foam::label Foam::sigmaColumns::findNearestColumn
(
    const point& pt,
    scalar& dist
) const
{
    label nearest = -1;
    dist = GREAT;

    forAll(columnCentre_, coli)
    {
        const scalar d = horizontalDistance(columnCentre_[coli], pt);

        if (d < dist)
        {
            dist = d;
            nearest = coli;
        }
    }

    return nearest;
}


// ************************************************************************* //
