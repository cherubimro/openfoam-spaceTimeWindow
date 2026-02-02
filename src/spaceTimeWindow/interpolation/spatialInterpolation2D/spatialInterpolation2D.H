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

Class
    Foam::spatialInterpolation2D

Description
    Spatial interpolation for 2D point clouds on planar surfaces.

    Supports two modes:
    - Refinement (more target points than source): Barycentric interpolation
      using 2D Delaunay triangulation of source points
    - Coarsening (fewer target points than source): Area-weighted averaging
      of all source points within a search radius

    Designed for interpolating boundary data on box faces (rectangular planes).
    Each box face is processed independently as a 2D problem.

SourceFiles
    spatialInterpolation2D.C

\*---------------------------------------------------------------------------*/

#ifndef spatialInterpolation2D_H
#define spatialInterpolation2D_H

#include "vectorField.H"
#include "scalarField.H"
#include "labelList.H"
#include "boundBox.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class spatialInterpolation2D Declaration
\*---------------------------------------------------------------------------*/

class spatialInterpolation2D
{
public:

    //- Interpolation mode enumeration
    enum class interpMode
    {
        NONE,           //!< No interpolation (same resolution)
        BARYCENTRIC,    //!< Barycentric interpolation (refinement)
        AREA_WEIGHTED   //!< Area-weighted averaging (coarsening)
    };

private:

    // Private Data

        //- Interpolation mode
        interpMode mode_;

        //- Number of source points
        label nSrcPoints_;

        //- Number of target points
        label nTgtPoints_;

        //- Source point indices for each target (3 for barycentric, N for area-weighted)
        labelListList sourceIndices_;

        //- Interpolation weights for each target
        List<scalarList> weights_;

        //- Octree for source point lookup (coarsening mode)
        autoPtr<indexedOctree<treeDataPoint>> octreePtr_;


    // Private Member Functions

        //- Identify which box face a point belongs to (0-5 for +-X, +-Y, +-Z)
        static label identifyBoxFace(const point& p, const boundBox& bb);

        //- Project 3D point to 2D for triangulation (drops the constant coordinate)
        static vector2D projectTo2D(const point& p, label faceIndex);

        //- Compute 2D Delaunay triangulation of points
        static labelListList triangulate2D
        (
            const vectorField& points,
            const labelList& indices,
            label faceIndex
        );

        //- Find containing triangle and compute barycentric coordinates
        static bool findTriangleAndBaryCoords
        (
            const vectorField& srcPoints,
            const labelListList& triangles,
            const point& tgtPoint,
            label faceIndex,
            labelList& triIndices,
            scalarList& baryCoords
        );

        //- Compute 2D barycentric coordinates
        static vector barycentricCoords2D
        (
            const vector2D& p,
            const vector2D& a,
            const vector2D& b,
            const vector2D& c
        );

        //- Build interpolation for refinement (barycentric)
        void buildBarycentricInterpolation
        (
            const vectorField& srcPoints,
            const vectorField& tgtPoints,
            const boundBox& bb
        );

        //- Build interpolation for coarsening (area-weighted)
        void buildAreaWeightedInterpolation
        (
            const vectorField& srcPoints,
            const vectorField& tgtPoints
        );


public:

    // Constructors

        //- Default constructor (no interpolation)
        spatialInterpolation2D();

        //- Construct from source and target points
        //  Automatically determines mode (none/barycentric/area-weighted)
        spatialInterpolation2D
        (
            const vectorField& srcPoints,
            const vectorField& tgtPoints,
            const boundBox& bb
        );


    // Member Functions

        //- Return interpolation mode
        interpMode mode() const { return mode_; }

        //- Return true if interpolation is needed
        bool needsInterpolation() const { return mode_ != interpMode::NONE; }

        //- Return true if this is refinement (barycentric)
        bool isRefinement() const { return mode_ == interpMode::BARYCENTRIC; }

        //- Return true if this is coarsening (area-weighted)
        bool isCoarsening() const { return mode_ == interpMode::AREA_WEIGHTED; }

        //- Return number of source points
        label nSrcPoints() const { return nSrcPoints_; }

        //- Return number of target points
        label nTgtPoints() const { return nTgtPoints_; }

        //- Interpolate scalar field
        tmp<scalarField> interpolate(const scalarField& srcField) const;

        //- Interpolate vector field
        tmp<vectorField> interpolate(const vectorField& srcField) const;

        //- Templated interpolation
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& srcField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
