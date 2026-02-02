/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "spaceTimeWindowInletOutletFvPatchField.H"
#include "Time.H"
#include "IFstream.H"
#include "dictionary.H"
#include "foamVersion.H"
#include "token.H"
#include "IOstreamOption.H"
#include "deltaVarintCodec.H"
#include "pTraits.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "boundBox.H"
#include <type_traits>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::findSampleTimes() const
{
    const fileName dataPath
    (
        this->db().time().globalPath()
      / dataDir_
      / this->patch().name()
    );

    // Get list of time directories
    fileNameList timeDirs = Foam::readDir(dataPath, fileName::DIRECTORY);

    // Convert to instants and sort
    DynamicList<instant> times;

    for (const fileName& dir : timeDirs)
    {
        scalar timeValue = 0;
        if (readScalar(dir, timeValue))
        {
            times.append(instant(timeValue, dir));
        }
    }

    // Sort by time value
    Foam::stableSort(times);

    sampleTimes_.transfer(times);

    if (sampleTimes_.empty())
    {
        FatalErrorInFunction
            << "No time directories found in " << dataPath << nl
            << "    Patch: " << this->patch().name() << nl
            << exit(FatalError);
    }

    DebugInfo
        << "Found " << sampleTimes_.size() << " sample times in "
        << dataPath << endl;
}


template<class Type>
Foam::Field<Type> Foam::spaceTimeWindowInletOutletFvPatchField<Type>::readFieldData
(
    const instant& timeDir
) const
{
    // Initialize spatial coarsening map if needed (reads points file once)
    initSpatialCoarseningMap();

    const fileName baseDir
    (
        this->db().time().globalPath()
      / dataDir_
      / this->patch().name()
      / timeDir.name()
    );

    // Check for delta-varint compressed file first (.dvz extension)
    const fileName dvzPath = baseDir / (fieldTableName_ + "." + deltaVarintCodec::fileExtension());

    Field<Type> fineData;

    if (isFile(dvzPath))
    {
        // Read using delta-varint codec
        if constexpr (std::is_same_v<Type, scalar>)
        {
            fineData = deltaVarintCodec::readScalar(dvzPath);
        }
        else if constexpr (std::is_same_v<Type, vector>)
        {
            fineData = deltaVarintCodec::readVector(dvzPath);
        }
        else
        {
            FatalErrorInFunction
                << "Delta-varint codec only supports scalar and vector fields" << nl
                << "    File: " << dvzPath << nl
                << "    Type: " << pTraits<Type>::typeName << nl
                << exit(FatalError);
        }
    }
    else
    {
        // Fall back to standard OpenFOAM format
        const fileName dataPath = baseDir / fieldTableName_;

        IFstream is(dataPath);

        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open file " << dataPath << nl
                << "    (also checked for " << dvzPath << ")" << nl
                << exit(FatalError);
        }

        // Check if file has FoamFile header
        token firstToken(is);

        if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
        {
            // File has FoamFile header - read header dictionary content
            dictionary headerDict;
            headerDict.read(is, false);

            // Parse format from header and set stream format
            if (headerDict.found("format"))
            {
                const word formatStr = headerDict.get<word>("format");
                if (formatStr == "binary")
                {
                    is.format(IOstreamOption::BINARY);
                }
                else
                {
                    is.format(IOstreamOption::ASCII);
                }
            }

            // Read field data
            fineData = Field<Type>(is);
        }
        else
        {
            // No header - put token back and read raw field data
            is.putBack(firstToken);
            fineData = Field<Type>(is);
        }
    }

    // Check if data size matches patch (no coarsening needed)
    if (fineData.size() == this->patch().size())
    {
        return fineData;
    }

    // Check for refinement mode (barycentric interpolation from coarse to fine)
    if (refinementMode_ && fineToTriangle_.size() == this->patch().size())
    {
        // Interpolate from coarse boundaryData to fine mesh using barycentric weights
        Field<Type> fineResult(this->patch().size(), Zero);

        forAll(fineResult, fineI)
        {
            const label triI = fineToTriangle_[fineI];
            const FixedList<label, 3>& tri = triangles_[triI];
            const barycentric2D& bary = fineBaryWeights_[fineI];

            // Barycentric interpolation: value = u*v0 + v*v1 + w*v2
            fineResult[fineI] =
                bary[0] * fineData[tri[0]]
              + bary[1] * fineData[tri[1]]
              + bary[2] * fineData[tri[2]];
        }

        return fineResult;
    }

    // Check if spatial coarsening map is available (coarsening mode)
    if (coarseToFineMap_.empty())
    {
        FatalErrorInFunction
            << "Field size " << fineData.size()
            << " does not match patch size " << this->patch().size() << nl
            << "    and no spatial mapping available" << nl
            << "    Patch: " << this->patch().name() << nl
            << "    Ensure boundaryData/points file exists for mesh resolution support" << nl
            << exit(FatalError);
    }

    // Apply area-weighted averaging using the coarsening map
    Field<Type> coarseData(this->patch().size(), Zero);

    forAll(coarseToFineMap_, coarseI)
    {
        const labelList& fineFaces = coarseToFineMap_[coarseI];

        if (fineFaces.empty())
        {
            continue;
        }

        // Simple average of fine face values
        // (area-weighted would require face areas, using equal weights for now)
        Type sum = Zero;
        for (const label fineI : fineFaces)
        {
            sum += fineData[fineI];
        }
        coarseData[coarseI] = sum / scalar(fineFaces.size());
    }

    return coarseData;
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::checkTable()
{
    if (sampleTimes_.empty())
    {
        findSampleTimes();
        sampleTimeIndex0_ = -1;
        sampleTimeIndex1_ = -1;
        sampleTimeIndex2_ = -1;
        sampleTimeIndex3_ = -1;
    }

    const scalar currentTime = this->db().time().value();

    // Find bracketing times (i1 <= currentTime <= i2)
    label i1 = -1;
    label i2 = -1;

    for (label i = 0; i < sampleTimes_.size(); ++i)
    {
        if (sampleTimes_[i].value() <= currentTime)
        {
            i1 = i;
        }
        if (sampleTimes_[i].value() >= currentTime && i2 < 0)
        {
            i2 = i;
        }
    }

    // Strict time bounds - NO extrapolation allowed
    if (i1 < 0)
    {
        FatalErrorInFunction
            << "Current time " << currentTime
            << " is before first available sample time "
            << sampleTimes_[0].value() << nl
            << "    Patch: " << this->patch().name() << nl
            << "    This BC requires boundary data to be pre-computed" << nl
            << "    (run spaceTimeWindowExtract first)" << nl
            << exit(FatalError);
    }
    if (i2 < 0)
    {
        FatalErrorInFunction
            << "Current time " << currentTime
            << " is after last available sample time "
            << sampleTimes_.last().value() << nl
            << "    Patch: " << this->patch().name() << nl
            << "    This BC requires boundary data to be pre-computed" << nl
            << "    (ensure extraction covers full simulation time)" << nl
            << exit(FatalError);
    }

    // Determine indices for interpolation
    // For linear: need i1, i2
    // For cubic: need i0, i1, i2, i3 (clamp at boundaries)
    label i0 = max(i1 - 1, label(0));
    label i3 = min(i2 + 1, sampleTimes_.size() - 1);

    // Load data if indices changed
    if (i0 != sampleTimeIndex0_)
    {
        sampleTimeIndex0_ = i0;
        sampledValues0_ = readFieldData(sampleTimes_[i0]);
        DebugInfo << "Loaded time[0] data: " << sampleTimes_[i0].name() << endl;
    }

    if (i1 != sampleTimeIndex1_)
    {
        sampleTimeIndex1_ = i1;
        sampledValues1_ = readFieldData(sampleTimes_[i1]);
        DebugInfo << "Loaded time[1] data: " << sampleTimes_[i1].name() << endl;
    }

    if (i2 != sampleTimeIndex2_)
    {
        sampleTimeIndex2_ = i2;
        sampledValues2_ = readFieldData(sampleTimes_[i2]);
        DebugInfo << "Loaded time[2] data: " << sampleTimes_[i2].name() << endl;
    }

    if (i3 != sampleTimeIndex3_)
    {
        sampleTimeIndex3_ = i3;
        sampledValues3_ = readFieldData(sampleTimes_[i3]);
        DebugInfo << "Loaded time[3] data: " << sampleTimes_[i3].name() << endl;
    }
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::validateMetadata() const
{
    if (metadataValidated_)
    {
        return;
    }

    const fileName metadataPath
    (
        this->db().time().globalPath()
      / dataDir_
      / this->patch().name()
      / "extractionMetadata"
    );

    IFstream is(metadataPath);

    if (!is.good())
    {
        WarningInFunction
            << "No extraction metadata found at " << metadataPath << nl
            << "    Cannot validate time settings." << nl
            << "    Ensure extraction was done with compatible time settings." << endl;

        metadataValidated_ = true;
        return;
    }

    // Read metadata dictionary
    dictionary metadata(is);

    const Time& runTime = this->db().time();

    // Validate OpenFOAM version (warning only)
    if (metadata.found("openfoamApi"))
    {
        const label extractionApi = metadata.get<label>("openfoamApi");
        const label currentApi = foamVersion::api;

        if (extractionApi != currentApi)
        {
            WarningInFunction
                << "OpenFOAM API version mismatch!" << nl
                << "    Extraction openfoamApi: " << extractionApi << nl
                << "    Current openfoamApi:    " << currentApi << nl
                << "    Patch: " << this->patch().name() << endl;
        }
    }

    // Validate solver name (warning only)
    if (metadata.found("solver"))
    {
        const word extractionSolver = metadata.get<word>("solver");
        const word currentSolver = runTime.controlDict().getOrDefault<word>
        (
            "application",
            word("unknown")
        );

        if (extractionSolver != currentSolver && extractionSolver != "unknown")
        {
            WarningInFunction
                << "Solver mismatch!" << nl
                << "    Extraction solver: " << extractionSolver << nl
                << "    Current solver:    " << currentSolver << nl
                << "    Patch: " << this->patch().name() << endl;
        }
    }

    // Validate deltaT if time interpolation is disabled
    if (!allowTimeInterpolation_ && metadata.found("deltaT"))
    {
        const scalar extractionDeltaT = metadata.get<scalar>("deltaT");
        const scalar currentDeltaT = runTime.deltaTValue();

        if (mag(extractionDeltaT - currentDeltaT) > SMALL * extractionDeltaT)
        {
            FatalErrorInFunction
                << "Time step mismatch!" << nl
                << "    Extraction deltaT: " << extractionDeltaT << nl
                << "    Current deltaT:    " << currentDeltaT << nl
                << "    Patch: " << this->patch().name() << nl << nl
                << "    Options:" << nl
                << "      1. Set deltaT = " << extractionDeltaT << nl
                << "      2. Enable time interpolation: allowTimeInterpolation true" << nl
                << exit(FatalError);
        }
    }

    Info<< "spaceTimeWindowInletOutlet BC on patch " << this->patch().name()
        << ": validated extraction metadata" << nl;

    if (metadata.found("openfoamVersion"))
    {
        Info<< "    openfoamVersion = " << metadata.get<word>("openfoamVersion") << nl;
    }
    Info<< "    deltaT = " << runTime.deltaTValue() << endl;

    metadataValidated_ = true;
}


template<class Type>
Type Foam::spaceTimeWindowInletOutletFvPatchField<Type>::catmullRomCentripetal
(
    const Type& p0,
    const Type& p1,
    const Type& p2,
    const Type& p3,
    scalar t0,
    scalar t1,
    scalar t2,
    scalar t3,
    scalar t
) const
{
    // Centripetal Catmull-Rom spline interpolation
    // Handles non-uniform time spacing correctly

    // Compute knot intervals using centripetal parameterization
    const scalar dt01 = Foam::sqrt(Foam::mag(t1 - t0));
    const scalar dt12 = Foam::sqrt(Foam::mag(t2 - t1));
    const scalar dt23 = Foam::sqrt(Foam::mag(t3 - t2));

    const scalar eps = SMALL;

    // Knot values (cumulative)
    const scalar k0 = 0.0;
    const scalar k1 = k0 + max(dt01, eps);
    const scalar k2 = k1 + max(dt12, eps);
    const scalar k3 = k2 + max(dt23, eps);

    // Map input time t to knot parameter u
    const scalar alpha = (t - t1) / (t2 - t1 + eps);
    const scalar u = k1 + alpha * (k2 - k1);

    // Barry-Goldman pyramid algorithm
    const Type A1 = (k1 - u)/(k1 - k0 + eps) * p0 + (u - k0)/(k1 - k0 + eps) * p1;
    const Type A2 = (k2 - u)/(k2 - k1 + eps) * p1 + (u - k1)/(k2 - k1 + eps) * p2;
    const Type A3 = (k3 - u)/(k3 - k2 + eps) * p2 + (u - k2)/(k3 - k2 + eps) * p3;

    const Type B1 = (k2 - u)/(k2 - k0 + eps) * A1 + (u - k0)/(k2 - k0 + eps) * A2;
    const Type B2 = (k3 - u)/(k3 - k1 + eps) * A2 + (u - k1)/(k3 - k1 + eps) * A3;

    const Type C = (k2 - u)/(k2 - k1 + eps) * B1 + (u - k1)/(k2 - k1 + eps) * B2;

    return C;
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::buildRefinementTriangulation
(
    const vectorField& coarsePoints,
    const vectorField& finePoints
) const
{
    // Simple approach: for each fine face, find the 3 nearest coarse points
    // that form a valid (non-degenerate) triangle containing or near the point.
    // Use barycentric interpolation from those 3 points.
    //
    // This works for arbitrary point distributions without assuming any
    // grid structure or planar organization.

    const label nCoarse = coarsePoints.size();
    const label nFine = finePoints.size();

    Info<< "    Building nearest-3-point triangulation" << nl
        << "    Coarse points: " << nCoarse << ", Fine points: " << nFine << endl;

    // For each fine point, we store the 3 coarse point indices and weights
    // We don't need a global triangle list - just per-fine-point data
    triangles_.setSize(nFine);
    fineToTriangle_.setSize(nFine);
    fineBaryWeights_.setSize(nFine);

    label perfectMatches = 0;
    label goodMatches = 0;
    label fallbackMatches = 0;

    forAll(finePoints, fineI)
    {
        const point& fp = finePoints[fineI];

        // Find the k nearest coarse points (k=6 to have options for triangles)
        const label kNearest = 6;
        labelList nearest(kNearest, -1);
        scalarList nearDist(kNearest, GREAT);

        forAll(coarsePoints, coarseI)
        {
            scalar d = magSqr(fp - coarsePoints[coarseI]);

            for (label k = 0; k < kNearest; k++)
            {
                if (d < nearDist[k])
                {
                    // Shift others down
                    for (label m = kNearest - 1; m > k; m--)
                    {
                        nearest[m] = nearest[m-1];
                        nearDist[m] = nearDist[m-1];
                    }
                    nearest[k] = coarseI;
                    nearDist[k] = d;
                    break;
                }
            }
        }

        // Try to form a triangle from the nearest points that contains fp
        // or has fp close to it with valid barycentric coords
        bool found = false;
        label bestI0 = nearest[0];
        label bestI1 = nearest[1];
        label bestI2 = nearest[2];
        scalar bestU = 1.0/3.0;
        scalar bestV = 1.0/3.0;
        scalar bestW = 1.0/3.0;
        scalar bestError = GREAT;

        // Try combinations of the k nearest points
        for (label i = 0; i < kNearest && !found; i++)
        {
            for (label j = i + 1; j < kNearest && !found; j++)
            {
                for (label k = j + 1; k < kNearest && !found; k++)
                {
                    if (nearest[i] < 0 || nearest[j] < 0 || nearest[k] < 0)
                        continue;

                    const point& p0 = coarsePoints[nearest[i]];
                    const point& p1 = coarsePoints[nearest[j]];
                    const point& p2 = coarsePoints[nearest[k]];

                    // Compute barycentric coordinates in 3D
                    // Project onto the plane of the triangle
                    vector v0 = p1 - p0;
                    vector v1 = p2 - p0;
                    vector v2 = fp - p0;

                    scalar d00 = v0 & v0;
                    scalar d01 = v0 & v1;
                    scalar d11 = v1 & v1;
                    scalar d20 = v2 & v0;
                    scalar d21 = v2 & v1;

                    scalar denom = d00 * d11 - d01 * d01;

                    // Skip degenerate triangles
                    if (mag(denom) < VSMALL) continue;

                    scalar v = (d11 * d20 - d01 * d21) / denom;
                    scalar w = (d00 * d21 - d01 * d20) / denom;
                    scalar u = 1.0 - v - w;

                    // Check if point is inside triangle
                    const scalar tol = 0.001;
                    if (u >= -tol && v >= -tol && w >= -tol &&
                        u <= 1.0 + tol && v <= 1.0 + tol && w <= 1.0 + tol)
                    {
                        // Found a containing triangle
                        bestI0 = nearest[i];
                        bestI1 = nearest[j];
                        bestI2 = nearest[k];
                        bestU = u;
                        bestV = v;
                        bestW = w;
                        found = true;
                        perfectMatches++;
                        break;
                    }

                    // Track the best (closest to valid) triangle
                    scalar error = 0;
                    if (u < 0) error += u * u;
                    if (v < 0) error += v * v;
                    if (w < 0) error += w * w;
                    if (u > 1) error += (u - 1) * (u - 1);
                    if (v > 1) error += (v - 1) * (v - 1);
                    if (w > 1) error += (w - 1) * (w - 1);

                    if (error < bestError)
                    {
                        bestError = error;
                        bestI0 = nearest[i];
                        bestI1 = nearest[j];
                        bestI2 = nearest[k];
                        bestU = u;
                        bestV = v;
                        bestW = w;
                    }
                }
            }
        }

        if (!found)
        {
            if (bestError < 0.01)
            {
                goodMatches++;
            }
            else
            {
                fallbackMatches++;
            }

            // Clamp barycentric coordinates
            bestU = max(0.0, min(1.0, bestU));
            bestV = max(0.0, min(1.0, bestV));
            bestW = max(0.0, min(1.0, bestW));

            // Renormalize
            scalar sum = bestU + bestV + bestW;
            if (sum > SMALL)
            {
                bestU /= sum;
                bestV /= sum;
                bestW /= sum;
            }
            else
            {
                bestU = bestV = bestW = 1.0/3.0;
            }
        }

        // Store the triangle and weights
        triangles_[fineI][0] = bestI0;
        triangles_[fineI][1] = bestI1;
        triangles_[fineI][2] = bestI2;
        fineToTriangle_[fineI] = fineI;  // Each fine face has its own "triangle"
        fineBaryWeights_[fineI] = barycentric2D(bestU, bestV, bestW);
    }

    Info<< "    Interpolation quality:" << nl
        << "      Perfect (inside triangle): " << perfectMatches << nl
        << "      Good (small extrapolation): " << goodMatches << nl
        << "      Fallback (clamped): " << fallbackMatches << endl;
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::initSpatialCoarseningMap() const
{
    if (spatialMapInitialized_)
    {
        return;
    }

    const label patchSize = this->patch().size();

    // Read points file to get fine face centers
    const fileName pointsPath
    (
        this->db().time().globalPath()
      / dataDir_
      / this->patch().name()
      / "points"
    );

    if (!isFile(pointsPath))
    {
        // No points file - can't build spatial mapping
        spatialCoarsenFactor_ = 1;
        spatialMapInitialized_ = true;
        return;
    }

    IFstream pointsIs(pointsPath);
    if (!pointsIs.good())
    {
        spatialCoarsenFactor_ = 1;
        spatialMapInitialized_ = true;
        return;
    }

    // Skip FoamFile header if present
    token firstToken(pointsIs);
    if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
    {
        dictionary headerDict;
        headerDict.read(pointsIs, false);
    }
    else
    {
        pointsIs.putBack(firstToken);
    }

    // Read face centers from points file
    vectorField finePoints(pointsIs);
    const label nFineFaces = finePoints.size();

    if (nFineFaces == patchSize)
    {
        // No coarsening needed
        spatialCoarsenFactor_ = 1;
        spatialMapInitialized_ = true;
        return;
    }

    if (nFineFaces < patchSize)
    {
        // Mesh refinement case: need to interpolate from coarse boundaryData
        // to fine mesh using barycentric interpolation on triangulated surface
        refinementMode_ = true;
        spatialCoarsenFactor_ = -1;  // Negative indicates refinement mode

        Info<< "spaceTimeWindowInletOutlet BC on patch " << this->patch().name()
            << ": detected mesh refinement" << nl
            << "    BoundaryData faces (coarse): " << nFineFaces << nl
            << "    Patch faces (fine): " << patchSize << nl
            << "    Using barycentric triangulation interpolation" << endl;

        // Store coarse points and build triangulation
        coarsePoints_ = finePoints;  // finePoints from file = coarse data

        // Get fine face centers from mesh patch
        const vectorField& fineMeshPoints = this->patch().Cf();

        // Build triangulation and compute barycentric weights for all fine faces
        buildRefinementTriangulation(coarsePoints_, fineMeshPoints);

        spatialMapInitialized_ = true;
        return;
    }

    // Calculate coarsening ratio
    const scalar ratio = scalar(nFineFaces) / scalar(patchSize);

    // For 2:1 coarsening in each dimension of a 2D face grid,
    // we expect ratio ~ 4. For 3:1, ratio ~ 9, etc.
    spatialCoarsenFactor_ = label(Foam::sqrt(ratio) + 0.5);

    if (spatialCoarsenFactor_ < 2)
    {
        spatialCoarsenFactor_ = 2;
    }

    Info<< "spaceTimeWindowInletOutlet BC on patch " << this->patch().name()
        << ": detected spatial coarsening" << nl
        << "    BoundaryData faces: " << nFineFaces << nl
        << "    Patch faces: " << patchSize << nl
        << "    Coarsening factor: " << spatialCoarsenFactor_ << "x" << spatialCoarsenFactor_
        << " (ratio " << ratio << ")" << endl;

    // Get coarse face centers from the mesh patch
    const vectorField& coarsePoints = this->patch().Cf();

    // Build mapping: for each coarse face, find the closest fine faces
    // Use a cell-based approach: group fine faces whose centers fall within
    // the vicinity of each coarse face

    coarseToFineMap_.setSize(patchSize);

    // Calculate typical coarse face spacing
    boundBox coarseBb(coarsePoints);
    const scalar coarseSpacing = Foam::cbrt(coarseBb.volume() / max(patchSize, 1));
    const scalar searchRadius = coarseSpacing * 0.6 * spatialCoarsenFactor_;

    // For each coarse face, find nearby fine faces
    labelList fineFaceUsed(nFineFaces, 0);

    forAll(coarsePoints, coarseI)
    {
        const point& coarseCf = coarsePoints[coarseI];
        DynamicList<label> nearbyFine;

        forAll(finePoints, fineI)
        {
            if (mag(finePoints[fineI] - coarseCf) < searchRadius)
            {
                nearbyFine.append(fineI);
                fineFaceUsed[fineI]++;
            }
        }

        coarseToFineMap_[coarseI].transfer(nearbyFine);
    }

    // Check coverage and warn if some fine faces weren't mapped
    label unmapped = 0;
    forAll(fineFaceUsed, fineI)
    {
        if (fineFaceUsed[fineI] == 0)
        {
            unmapped++;
        }
    }

    if (unmapped > 0)
    {
        WarningInFunction
            << unmapped << " fine faces were not mapped to any coarse face" << nl
            << "    This may indicate an unsupported coarsening ratio" << nl
            << "    Patch: " << this->patch().name() << endl;
    }

    // Check that each coarse face has at least one fine face
    label emptyCoarse = 0;
    forAll(coarseToFineMap_, coarseI)
    {
        if (coarseToFineMap_[coarseI].empty())
        {
            emptyCoarse++;
        }
    }

    if (emptyCoarse > 0)
    {
        FatalErrorInFunction
            << emptyCoarse << " coarse faces have no corresponding fine faces" << nl
            << "    Check that boundaryData points file matches extraction mesh" << nl
            << "    Patch: " << this->patch().name() << nl
            << exit(FatalError);
    }

    spatialMapInitialized_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::spaceTimeWindowInletOutletFvPatchField<Type>::spaceTimeWindowInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    dataDir_("constant/boundaryData"),
    fieldTableName_(iF.name()),
    phiName_("phi"),
    fixesValue_(false),
    allowTimeInterpolation_(true),
    timeInterpolationScheme_(timeInterpScheme::LINEAR),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(false),
    spatialCoarsenFactor_(1),
    coarseToFineMap_(),
    spatialMapInitialized_(false),
    refinementMode_(false),
    coarsePoints_(),
    triangles_(),
    fineToTriangle_(),
    fineBaryWeights_()
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::spaceTimeWindowInletOutletFvPatchField<Type>::spaceTimeWindowInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    dataDir_(dict.getOrDefault<fileName>("dataDir", "constant/boundaryData")),
    fieldTableName_(dict.getOrDefault<word>("fieldTableName", iF.name())),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    fixesValue_(dict.getOrDefault("fixesValue", false)),
    allowTimeInterpolation_(dict.getOrDefault("allowTimeInterpolation", true)),
    timeInterpolationScheme_(timeInterpScheme::LINEAR),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(false),
    spatialCoarsenFactor_(1),
    coarseToFineMap_(),
    spatialMapInitialized_(false),
    refinementMode_(false),
    coarsePoints_(),
    triangles_(),
    fineToTriangle_(),
    fineBaryWeights_()
{
    // Parse time interpolation scheme
    const word schemeStr = dict.getOrDefault<word>("timeInterpolationScheme", "linear");
    if (schemeStr == "cubic")
    {
        timeInterpolationScheme_ = timeInterpScheme::CUBIC;
    }
    else if (schemeStr == "none")
    {
        timeInterpolationScheme_ = timeInterpScheme::NONE;
        allowTimeInterpolation_ = false;
    }
    else if (schemeStr != "linear")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown timeInterpolationScheme: " << schemeStr << nl
            << "    Valid options: none, linear, cubic" << nl
            << exit(FatalIOError);
    }

    // Read initial value
    this->patchType() = dict.getOrDefault<word>("patchType", word::null);

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    // Initialize mixed BC fields
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;  // Start with zeroGradient (outflow assumption)
}


template<class Type>
Foam::spaceTimeWindowInletOutletFvPatchField<Type>::spaceTimeWindowInletOutletFvPatchField
(
    const spaceTimeWindowInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    phiName_(ptf.phiName_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(ptf.metadataValidated_),
    spatialCoarsenFactor_(1),
    coarseToFineMap_(),
    spatialMapInitialized_(false),
    refinementMode_(false),
    coarsePoints_(),
    triangles_(),
    fineToTriangle_(),
    fineBaryWeights_()
{}


template<class Type>
Foam::spaceTimeWindowInletOutletFvPatchField<Type>::spaceTimeWindowInletOutletFvPatchField
(
    const spaceTimeWindowInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    phiName_(ptf.phiName_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledValues0_(ptf.sampledValues0_),
    sampledValues1_(ptf.sampledValues1_),
    sampledValues2_(ptf.sampledValues2_),
    sampledValues3_(ptf.sampledValues3_),
    metadataValidated_(ptf.metadataValidated_),
    spatialCoarsenFactor_(ptf.spatialCoarsenFactor_),
    coarseToFineMap_(ptf.coarseToFineMap_),
    spatialMapInitialized_(ptf.spatialMapInitialized_),
    refinementMode_(ptf.refinementMode_),
    coarsePoints_(ptf.coarsePoints_),
    triangles_(ptf.triangles_),
    fineToTriangle_(ptf.fineToTriangle_),
    fineBaryWeights_(ptf.fineBaryWeights_)
{}


template<class Type>
Foam::spaceTimeWindowInletOutletFvPatchField<Type>::spaceTimeWindowInletOutletFvPatchField
(
    const spaceTimeWindowInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    phiName_(ptf.phiName_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledValues0_(ptf.sampledValues0_),
    sampledValues1_(ptf.sampledValues1_),
    sampledValues2_(ptf.sampledValues2_),
    sampledValues3_(ptf.sampledValues3_),
    metadataValidated_(ptf.metadataValidated_),
    spatialCoarsenFactor_(ptf.spatialCoarsenFactor_),
    coarseToFineMap_(ptf.coarseToFineMap_),
    spatialMapInitialized_(ptf.spatialMapInitialized_),
    refinementMode_(ptf.refinementMode_),
    coarsePoints_(ptf.coarsePoints_),
    triangles_(ptf.triangles_),
    fineToTriangle_(ptf.fineToTriangle_),
    fineBaryWeights_(ptf.fineBaryWeights_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);

    // Clear cached data
    sampledValues0_.clear();
    sampledValues1_.clear();
    sampledValues2_.clear();
    sampledValues3_.clear();
    sampleTimeIndex0_ = -1;
    sampleTimeIndex1_ = -1;
    sampleTimeIndex2_ = -1;
    sampleTimeIndex3_ = -1;
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);

    // Clear cached data
    sampledValues0_.clear();
    sampledValues1_.clear();
    sampledValues2_.clear();
    sampledValues3_.clear();
    sampleTimeIndex0_ = -1;
    sampleTimeIndex1_ = -1;
    sampleTimeIndex2_ = -1;
    sampleTimeIndex3_ = -1;
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Validate extraction metadata (once)
    validateMetadata();

    // Load/check sample times and data
    checkTable();

    const scalar currentTime = this->db().time().value();

    // Compute interpolated boundary values for inflow faces
    Field<Type> interpolatedValues(this->patch().size(), Zero);

    // Get bracketing times t1 and t2
    const scalar t1 = sampleTimes_[sampleTimeIndex1_].value();
    const scalar t2 = sampleTimes_[sampleTimeIndex2_].value();

    if (sampleTimeIndex1_ == sampleTimeIndex2_)
    {
        // Exact time match
        interpolatedValues = sampledValues1_;
    }
    else if (!allowTimeInterpolation_)
    {
        // Time interpolation disabled - require exact match
        const scalar timeTol = 0.01 * (t2 - t1);
        const bool closeToT1 = (mag(currentTime - t1) <= timeTol);
        const bool closeToT2 = (mag(currentTime - t2) <= timeTol);

        if (!closeToT1 && !closeToT2)
        {
            FatalErrorInFunction
                << "Current time " << currentTime
                << " does not match any available sample time." << nl
                << "    Nearest lower sample: " << t1 << nl
                << "    Nearest upper sample: " << t2 << nl
                << "    Patch: " << this->patch().name() << nl << nl
                << "    Time interpolation is disabled." << nl
                << "    Options:" << nl
                << "      1. Use identical deltaT as extraction" << nl
                << "      2. Set allowTimeInterpolation true" << nl
                << exit(FatalError);
        }

        interpolatedValues = closeToT2 ? sampledValues2_ : sampledValues1_;
    }
    else
    {
        // Temporal interpolation enabled
        const scalar alpha = (currentTime - t1) / (t2 - t1 + SMALL);

        const bool useCubic =
            (timeInterpolationScheme_ == timeInterpScheme::CUBIC)
         && (sampleTimeIndex0_ != sampleTimeIndex1_)
         && (sampleTimeIndex2_ != sampleTimeIndex3_);

        if (useCubic)
        {
            // Centripetal Catmull-Rom cubic spline
            const scalar t0 = sampleTimes_[sampleTimeIndex0_].value();
            const scalar t3 = sampleTimes_[sampleTimeIndex3_].value();

            forAll(interpolatedValues, facei)
            {
                interpolatedValues[facei] = catmullRomCentripetal
                (
                    sampledValues0_[facei],
                    sampledValues1_[facei],
                    sampledValues2_[facei],
                    sampledValues3_[facei],
                    t0, t1, t2, t3,
                    currentTime
                );
            }

            DebugInfo
                << "Cubic interpolation: t=" << currentTime
                << " using t0=" << t0 << " t1=" << t1
                << " t2=" << t2 << " t3=" << t3 << endl;
        }
        else
        {
            // Linear interpolation
            interpolatedValues =
                (1.0 - alpha) * sampledValues1_
              + alpha * sampledValues2_;

            DebugInfo
                << "Linear interpolation: t=" << currentTime
                << " between " << t1 << " and " << t2
                << " alpha=" << alpha << endl;
        }
    }

    // Set refValue to interpolated boundary data
    this->refValue() = interpolatedValues;

    // Set refGrad to zero (zeroGradient for outflow)
    this->refGrad() = Zero;

    // Determine valueFraction based on flux direction
    // valueFraction = 1 -> Dirichlet (inflow: use refValue)
    // valueFraction = 0 -> Neumann (outflow: use refGrad = 0 = zeroGradient)
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        const scalarField& phip = phi.boundaryField()[this->patch().index()];

        // Inflow: phi < 0 (flux into domain) -> valueFraction = 1
        // Outflow: phi >= 0 (flux out of domain) -> valueFraction = 0
        this->valueFraction() = 1.0 - pos(phip);
    }
    else
    {
        // No flux field available - assume all inflow (Dirichlet everywhere)
        WarningInFunction
            << "Flux field " << phiName_ << " not found." << nl
            << "    Assuming all faces are inflow (Dirichlet)." << nl
            << "    Patch: " << this->patch().name() << endl;

        this->valueFraction() = 1.0;
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::spaceTimeWindowInletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntry("dataDir", dataDir_);

    if (fieldTableName_ != this->internalField().name())
    {
        os.writeEntry("fieldTableName", fieldTableName_);
    }

    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);

    if (fixesValue_)
    {
        os.writeEntry("fixesValue", fixesValue_);
    }

    if (!allowTimeInterpolation_)
    {
        os.writeEntry("allowTimeInterpolation", allowTimeInterpolation_);
    }

    // Write interpolation scheme
    if (timeInterpolationScheme_ == timeInterpScheme::CUBIC)
    {
        os.writeEntry("timeInterpolationScheme", "cubic");
    }
    else if (timeInterpolationScheme_ == timeInterpScheme::NONE)
    {
        os.writeEntry("timeInterpolationScheme", "none");
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
