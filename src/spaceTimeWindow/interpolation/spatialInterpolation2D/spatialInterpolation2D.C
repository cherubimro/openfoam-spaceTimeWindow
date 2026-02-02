/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "spatialInterpolation2D.H"
#include "DynamicList.H"
#include "SortableList.H"
#include "triFace.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    // Tolerance for point comparisons
    const scalar pointTol = 1e-10;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::spatialInterpolation2D::identifyBoxFace
(
    const point& p,
    const boundBox& bb
)
{
    // Find which face of the bounding box this point is closest to
    // Returns: 0=minX, 1=maxX, 2=minY, 3=maxY, 4=minZ, 5=maxZ

    const point& mn = bb.min();
    const point& mx = bb.max();

    scalar minDist = GREAT;
    label faceIdx = -1;

    // Check distance to each face
    scalar distMinX = mag(p.x() - mn.x());
    scalar distMaxX = mag(p.x() - mx.x());
    scalar distMinY = mag(p.y() - mn.y());
    scalar distMaxY = mag(p.y() - mx.y());
    scalar distMinZ = mag(p.z() - mn.z());
    scalar distMaxZ = mag(p.z() - mx.z());

    if (distMinX < minDist) { minDist = distMinX; faceIdx = 0; }
    if (distMaxX < minDist) { minDist = distMaxX; faceIdx = 1; }
    if (distMinY < minDist) { minDist = distMinY; faceIdx = 2; }
    if (distMaxY < minDist) { minDist = distMaxY; faceIdx = 3; }
    if (distMinZ < minDist) { minDist = distMinZ; faceIdx = 4; }
    if (distMaxZ < minDist) { minDist = distMaxZ; faceIdx = 5; }

    return faceIdx;
}


Foam::vector2D Foam::spatialInterpolation2D::projectTo2D
(
    const point& p,
    label faceIndex
)
{
    // Project 3D point to 2D by dropping the constant coordinate
    // faceIndex: 0,1 = X faces (use Y,Z), 2,3 = Y faces (use X,Z), 4,5 = Z faces (use X,Y)

    switch (faceIndex)
    {
        case 0:
        case 1:
            return vector2D(p.y(), p.z());
        case 2:
        case 3:
            return vector2D(p.x(), p.z());
        case 4:
        case 5:
        default:
            return vector2D(p.x(), p.y());
    }
}


Foam::labelListList Foam::spatialInterpolation2D::triangulate2D
(
    const vectorField& points,
    const labelList& indices,
    label faceIndex
)
{
    // Simple 2D Delaunay triangulation using Bowyer-Watson algorithm
    // Returns list of triangles, each containing 3 indices into original indices list

    const label n = indices.size();

    if (n < 3)
    {
        return labelListList();
    }

    // Project points to 2D
    List<vector2D> pts2D(n);
    forAll(indices, i)
    {
        pts2D[i] = projectTo2D(points[indices[i]], faceIndex);
    }

    // Find bounding box for super-triangle
    vector2D minPt(GREAT, GREAT);
    vector2D maxPt(-GREAT, -GREAT);

    forAll(pts2D, i)
    {
        minPt.x() = min(minPt.x(), pts2D[i].x());
        minPt.y() = min(minPt.y(), pts2D[i].y());
        maxPt.x() = max(maxPt.x(), pts2D[i].x());
        maxPt.y() = max(maxPt.y(), pts2D[i].y());
    }

    // Create super-triangle that contains all points
    scalar dx = maxPt.x() - minPt.x();
    scalar dy = maxPt.y() - minPt.y();
    scalar dmax = max(dx, dy) * 3.0;

    vector2D mid((minPt.x() + maxPt.x()) / 2, (minPt.y() + maxPt.y()) / 2);

    // Super-triangle vertices (indices n, n+1, n+2)
    List<vector2D> allPts(n + 3);
    forAll(pts2D, i)
    {
        allPts[i] = pts2D[i];
    }
    allPts[n] = vector2D(mid.x() - 2*dmax, mid.y() - dmax);
    allPts[n+1] = vector2D(mid.x() + 2*dmax, mid.y() - dmax);
    allPts[n+2] = vector2D(mid.x(), mid.y() + 2*dmax);

    // Initial triangulation with super-triangle
    DynamicList<triFace> triangles;
    triangles.append(triFace(n, n+1, n+2));

    // Lambda to check if point is inside circumcircle
    auto inCircumcircle = [&](const triFace& tri, const vector2D& p) -> bool
    {
        const vector2D& a = allPts[tri[0]];
        const vector2D& b = allPts[tri[1]];
        const vector2D& c = allPts[tri[2]];

        // Circumcircle test using determinant
        scalar ax = a.x() - p.x();
        scalar ay = a.y() - p.y();
        scalar bx = b.x() - p.x();
        scalar by = b.y() - p.y();
        scalar cx = c.x() - p.x();
        scalar cy = c.y() - p.y();

        scalar det =
            (ax*ax + ay*ay) * (bx*cy - cx*by)
          - (bx*bx + by*by) * (ax*cy - cx*ay)
          + (cx*cx + cy*cy) * (ax*by - bx*ay);

        // Check triangle orientation
        scalar orient = (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x());

        return (orient > 0) ? (det > 0) : (det < 0);
    };

    // Insert each point
    for (label i = 0; i < n; ++i)
    {
        const vector2D& p = allPts[i];

        // Find all triangles whose circumcircle contains point
        DynamicList<triFace> badTriangles;
        DynamicList<label> badIndices;

        forAll(triangles, ti)
        {
            if (inCircumcircle(triangles[ti], p))
            {
                badTriangles.append(triangles[ti]);
                badIndices.append(ti);
            }
        }

        // Find polygon hole boundary
        DynamicList<edge> polygon;

        forAll(badTriangles, bi)
        {
            const triFace& tri = badTriangles[bi];

            for (label ei = 0; ei < 3; ++ei)
            {
                edge e(tri[ei], tri[(ei+1)%3]);

                // Check if edge is shared with another bad triangle
                bool shared = false;
                forAll(badTriangles, bj)
                {
                    if (bi != bj)
                    {
                        const triFace& other = badTriangles[bj];
                        for (label ej = 0; ej < 3; ++ej)
                        {
                            edge e2(other[ej], other[(ej+1)%3]);
                            if (e == e2)
                            {
                                shared = true;
                                break;
                            }
                        }
                    }
                    if (shared) break;
                }

                if (!shared)
                {
                    polygon.append(e);
                }
            }
        }

        // Remove bad triangles (in reverse order to preserve indices)
        SortableList<label> sortedBad(badIndices);
        for (label j = sortedBad.size() - 1; j >= 0; --j)
        {
            triangles.remove(sortedBad[j]);
        }

        // Create new triangles from polygon edges to new point
        forAll(polygon, pi)
        {
            triangles.append(triFace(polygon[pi][0], polygon[pi][1], i));
        }
    }

    // Remove triangles that use super-triangle vertices
    DynamicList<triFace> finalTriangles;
    forAll(triangles, ti)
    {
        const triFace& tri = triangles[ti];
        if (tri[0] < n && tri[1] < n && tri[2] < n)
        {
            finalTriangles.append(tri);
        }
    }

    // Convert to labelListList with original indices
    labelListList result(finalTriangles.size());
    forAll(finalTriangles, ti)
    {
        result[ti].setSize(3);
        result[ti][0] = indices[finalTriangles[ti][0]];
        result[ti][1] = indices[finalTriangles[ti][1]];
        result[ti][2] = indices[finalTriangles[ti][2]];
    }

    return result;
}


Foam::vector Foam::spatialInterpolation2D::barycentricCoords2D
(
    const vector2D& p,
    const vector2D& a,
    const vector2D& b,
    const vector2D& c
)
{
    // Compute barycentric coordinates for point p in triangle abc
    scalar denom = (b.y() - c.y())*(a.x() - c.x()) + (c.x() - b.x())*(a.y() - c.y());

    if (mag(denom) < SMALL)
    {
        // Degenerate triangle
        return vector(1.0/3.0, 1.0/3.0, 1.0/3.0);
    }

    scalar w1 = ((b.y() - c.y())*(p.x() - c.x()) + (c.x() - b.x())*(p.y() - c.y())) / denom;
    scalar w2 = ((c.y() - a.y())*(p.x() - c.x()) + (a.x() - c.x())*(p.y() - c.y())) / denom;
    scalar w3 = 1.0 - w1 - w2;

    return vector(w1, w2, w3);
}


bool Foam::spatialInterpolation2D::findTriangleAndBaryCoords
(
    const vectorField& srcPoints,
    const labelListList& triangles,
    const point& tgtPoint,
    label faceIndex,
    labelList& triIndices,
    scalarList& baryCoords
)
{
    const vector2D p2D = projectTo2D(tgtPoint, faceIndex);

    // Search through triangles
    forAll(triangles, ti)
    {
        const labelList& tri = triangles[ti];

        vector2D a = projectTo2D(srcPoints[tri[0]], faceIndex);
        vector2D b = projectTo2D(srcPoints[tri[1]], faceIndex);
        vector2D c = projectTo2D(srcPoints[tri[2]], faceIndex);

        vector bary = barycentricCoords2D(p2D, a, b, c);

        // Check if point is inside triangle (all barycentric coords >= 0)
        // Use small negative tolerance to handle edge cases
        if (bary.x() >= -pointTol && bary.y() >= -pointTol && bary.z() >= -pointTol)
        {
            triIndices.setSize(3);
            triIndices[0] = tri[0];
            triIndices[1] = tri[1];
            triIndices[2] = tri[2];

            baryCoords.setSize(3);
            baryCoords[0] = max(0.0, bary.x());
            baryCoords[1] = max(0.0, bary.y());
            baryCoords[2] = max(0.0, bary.z());

            // Normalize to ensure sum = 1
            scalar sum = baryCoords[0] + baryCoords[1] + baryCoords[2];
            if (sum > SMALL)
            {
                baryCoords[0] /= sum;
                baryCoords[1] /= sum;
                baryCoords[2] /= sum;
            }

            return true;
        }
    }

    // Point not found in any triangle - find nearest triangle
    scalar minDist = GREAT;
    label nearestTri = 0;

    forAll(triangles, ti)
    {
        const labelList& tri = triangles[ti];

        // Compute centroid
        vector2D centroid(0, 0);
        for (label i = 0; i < 3; ++i)
        {
            centroid += projectTo2D(srcPoints[tri[i]], faceIndex);
        }
        centroid /= 3.0;

        scalar dist = mag(p2D - centroid);
        if (dist < minDist)
        {
            minDist = dist;
            nearestTri = ti;
        }
    }

    // Use nearest triangle with clamped barycentric coords
    const labelList& tri = triangles[nearestTri];
    vector2D a = projectTo2D(srcPoints[tri[0]], faceIndex);
    vector2D b = projectTo2D(srcPoints[tri[1]], faceIndex);
    vector2D c = projectTo2D(srcPoints[tri[2]], faceIndex);

    vector bary = barycentricCoords2D(p2D, a, b, c);

    triIndices.setSize(3);
    triIndices[0] = tri[0];
    triIndices[1] = tri[1];
    triIndices[2] = tri[2];

    baryCoords.setSize(3);
    baryCoords[0] = max(0.0, bary.x());
    baryCoords[1] = max(0.0, bary.y());
    baryCoords[2] = max(0.0, bary.z());

    scalar sum = baryCoords[0] + baryCoords[1] + baryCoords[2];
    if (sum > SMALL)
    {
        baryCoords[0] /= sum;
        baryCoords[1] /= sum;
        baryCoords[2] /= sum;
    }

    return false;
}


void Foam::spatialInterpolation2D::buildBarycentricInterpolation
(
    const vectorField& srcPoints,
    const vectorField& tgtPoints,
    const boundBox& bb
)
{
    // Group source points by box face
    labelListList srcByFace(6);
    {
        DynamicList<label> facePoints[6];
        forAll(srcPoints, i)
        {
            label fi = identifyBoxFace(srcPoints[i], bb);
            facePoints[fi].append(i);
        }
        for (label fi = 0; fi < 6; ++fi)
        {
            srcByFace[fi] = facePoints[fi];
        }
    }

    // Group target points by box face
    labelListList tgtByFace(6);
    {
        DynamicList<label> facePoints[6];
        forAll(tgtPoints, i)
        {
            label fi = identifyBoxFace(tgtPoints[i], bb);
            facePoints[fi].append(i);
        }
        for (label fi = 0; fi < 6; ++fi)
        {
            tgtByFace[fi] = facePoints[fi];
        }
    }

    // Triangulate each face and compute interpolation weights
    sourceIndices_.setSize(tgtPoints.size());
    weights_.setSize(tgtPoints.size());

    for (label fi = 0; fi < 6; ++fi)
    {
        const labelList& srcIdx = srcByFace[fi];
        const labelList& tgtIdx = tgtByFace[fi];

        if (srcIdx.size() < 3 || tgtIdx.empty())
        {
            continue;
        }

        // Triangulate source points on this face
        labelListList triangles = triangulate2D(srcPoints, srcIdx, fi);

        if (triangles.empty())
        {
            WarningInFunction
                << "Failed to triangulate face " << fi
                << " with " << srcIdx.size() << " source points"
                << endl;
            continue;
        }

        // For each target point on this face, find containing triangle
        forAll(tgtIdx, ti)
        {
            label tgtI = tgtIdx[ti];

            labelList triIndices;
            scalarList baryCoords;

            findTriangleAndBaryCoords
            (
                srcPoints,
                triangles,
                tgtPoints[tgtI],
                fi,
                triIndices,
                baryCoords
            );

            sourceIndices_[tgtI] = triIndices;
            weights_[tgtI] = baryCoords;
        }
    }
}


void Foam::spatialInterpolation2D::buildAreaWeightedInterpolation
(
    const vectorField& srcPoints,
    const vectorField& tgtPoints
)
{
    // Build octree for source points
    treeBoundBox overallBb(srcPoints);
    overallBb.grow(1e-4 * overallBb.mag());

    octreePtr_.reset
    (
        new indexedOctree<treeDataPoint>
        (
            treeDataPoint(srcPoints),
            overallBb,
            8,      // maxLevel
            10,     // leafSize
            3.0     // duplicity
        )
    );

    // Estimate typical spacing from source points
    scalar avgSpacing = Foam::pow(overallBb.volume() / srcPoints.size(), 1.0/3.0);
    scalar searchRadius = 2.0 * avgSpacing;

    sourceIndices_.setSize(tgtPoints.size());
    weights_.setSize(tgtPoints.size());

    forAll(tgtPoints, ti)
    {
        const point& tgt = tgtPoints[ti];

        // Find all source points within search radius
        DynamicList<label> foundIndices;
        DynamicList<scalar> foundDists;

        // Use sphere search
        labelList candidates = octreePtr_->findSphere(tgt, searchRadius);

        if (candidates.empty())
        {
            // Expand search radius
            scalar expandedRadius = searchRadius;
            while (candidates.empty() && expandedRadius < 10*searchRadius)
            {
                expandedRadius *= 2.0;
                candidates = octreePtr_->findSphere(tgt, expandedRadius);
            }
        }

        if (candidates.empty())
        {
            // Fall back to nearest point
            pointIndexHit nearest = octreePtr_->findNearest(tgt, GREAT);
            if (nearest.hit())
            {
                sourceIndices_[ti].setSize(1);
                sourceIndices_[ti][0] = nearest.index();
                weights_[ti].setSize(1);
                weights_[ti][0] = 1.0;
            }
            continue;
        }

        // Area-weighted average: equal weights for all source points in range
        // This is appropriate for coarsening where we want to average fine data
        sourceIndices_[ti] = candidates;
        weights_[ti].setSize(candidates.size());

        scalar weight = 1.0 / candidates.size();
        forAll(candidates, ci)
        {
            weights_[ti][ci] = weight;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spatialInterpolation2D::spatialInterpolation2D()
:
    mode_(interpMode::NONE),
    nSrcPoints_(0),
    nTgtPoints_(0),
    sourceIndices_(),
    weights_(),
    octreePtr_(nullptr)
{}


Foam::spatialInterpolation2D::spatialInterpolation2D
(
    const vectorField& srcPoints,
    const vectorField& tgtPoints,
    const boundBox& bb
)
:
    mode_(interpMode::NONE),
    nSrcPoints_(srcPoints.size()),
    nTgtPoints_(tgtPoints.size()),
    sourceIndices_(),
    weights_(),
    octreePtr_(nullptr)
{
    if (nSrcPoints_ == nTgtPoints_)
    {
        // Check if points are identical
        bool identical = true;
        forAll(srcPoints, i)
        {
            if (mag(srcPoints[i] - tgtPoints[i]) > pointTol)
            {
                identical = false;
                break;
            }
        }

        if (identical)
        {
            mode_ = interpMode::NONE;
            return;
        }
    }

    if (nTgtPoints_ > nSrcPoints_)
    {
        // Refinement: use barycentric interpolation
        mode_ = interpMode::BARYCENTRIC;
        buildBarycentricInterpolation(srcPoints, tgtPoints, bb);
    }
    else
    {
        // Coarsening: use area-weighted averaging
        mode_ = interpMode::AREA_WEIGHTED;
        buildAreaWeightedInterpolation(srcPoints, tgtPoints);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::spatialInterpolation2D::interpolate
(
    const scalarField& srcField
) const
{
    if (mode_ == interpMode::NONE)
    {
        return tmp<scalarField>(new scalarField(srcField));
    }

    tmp<scalarField> tresult(new scalarField(nTgtPoints_, 0.0));
    scalarField& result = tresult.ref();

    forAll(result, ti)
    {
        const labelList& indices = sourceIndices_[ti];
        const scalarList& wts = weights_[ti];

        forAll(indices, i)
        {
            result[ti] += wts[i] * srcField[indices[i]];
        }
    }

    return tresult;
}


Foam::tmp<Foam::vectorField> Foam::spatialInterpolation2D::interpolate
(
    const vectorField& srcField
) const
{
    if (mode_ == interpMode::NONE)
    {
        return tmp<vectorField>(new vectorField(srcField));
    }

    tmp<vectorField> tresult(new vectorField(nTgtPoints_, Zero));
    vectorField& result = tresult.ref();

    forAll(result, ti)
    {
        const labelList& indices = sourceIndices_[ti];
        const scalarList& wts = weights_[ti];

        forAll(indices, i)
        {
            result[ti] += wts[i] * srcField[indices[i]];
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Template specialization for generic Field types
template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::spatialInterpolation2D::interpolate
(
    const Field<Type>& srcField
) const
{
    if (mode_ == interpMode::NONE)
    {
        return tmp<Field<Type>>(new Field<Type>(srcField));
    }

    tmp<Field<Type>> tresult(new Field<Type>(nTgtPoints_, Zero));
    Field<Type>& result = tresult.ref();

    forAll(result, ti)
    {
        const labelList& indices = sourceIndices_[ti];
        const scalarList& wts = weights_[ti];

        forAll(indices, i)
        {
            result[ti] += wts[i] * srcField[indices[i]];
        }
    }

    return tresult;
}


// Explicit instantiations
template Foam::tmp<Foam::Field<Foam::scalar>>
Foam::spatialInterpolation2D::interpolate(const Field<scalar>&) const;

template Foam::tmp<Foam::Field<Foam::vector>>
Foam::spatialInterpolation2D::interpolate(const Field<vector>&) const;

template Foam::tmp<Foam::Field<Foam::tensor>>
Foam::spatialInterpolation2D::interpolate(const Field<tensor>&) const;

template Foam::tmp<Foam::Field<Foam::symmTensor>>
Foam::spatialInterpolation2D::interpolate(const Field<symmTensor>&) const;

template Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::spatialInterpolation2D::interpolate(const Field<sphericalTensor>&) const;


// ************************************************************************* //
