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

#include "spaceTimeWindowCoupledPressureFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IFstream.H"
#include "Time.H"
#include "deltaVarintCodec.H"
#include "deltaVarintTemporalCodec.H"
#include "zstdWrapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        spaceTimeWindowCoupledPressureFvPatchScalarField
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::findSampleTimes() const
{
    if (!sampleTimes_.empty())
    {
        return;
    }

    // Construct the boundaryData path
    fileName boundaryPath = dataDir_ / patch().name();

    // Find all time directories
    fileNameList timeDirs = readDir(boundaryPath, fileName::DIRECTORY);

    DynamicList<instant> times;
    for (const fileName& dir : timeDirs)
    {
        scalar timeValue;
        if (readScalar(dir, timeValue))
        {
            times.append(instant(timeValue, dir));
        }
    }

    // Sort by time value
    sampleTimes_.transfer(times);
    std::sort
    (
        sampleTimes_.begin(),
        sampleTimes_.end(),
        [](const instant& a, const instant& b)
        {
            return a.value() < b.value();
        }
    );

    if (sampleTimes_.empty())
    {
        FatalErrorInFunction
            << "No time directories found in " << boundaryPath << nl
            << exit(FatalError);
    }

    // Initialize time indices
    sampleTimeIndex0_ = -1;
    sampleTimeIndex1_ = -1;
    sampleTimeIndex2_ = -1;
    sampleTimeIndex3_ = -1;

    Info<< "spaceTimeWindowCoupledPressure: Found "
        << sampleTimes_.size() << " time directories" << nl
        << "    Time range: " << sampleTimes_.first().value()
        << " to " << sampleTimes_.last().value() << endl;
}


Foam::scalarField
Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::readPressureData
(
    const instant& timeDir
) const
{
    fileName dataPath = dataDir_ / patch().name() / timeDir.name();

    // Try different file formats
    // 1. Standard OpenFOAM format
    fileName filePath = dataPath / pressureFieldName_;
    if (isFile(filePath))
    {
        IFstream is(filePath);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open file " << filePath << nl
                << exit(FatalError);
        }

        // Parse FoamFile header if present (handles binary format)
        token t(is);
        if (t.isWord() && t.wordToken() == "FoamFile")
        {
            dictionary headerDict;
            headerDict.read(is, false);
            if (headerDict.found("format"))
            {
                const word formatStr = headerDict.get<word>("format");
                if (formatStr == "binary")
                {
                    is.format(IOstreamOption::BINARY);
                }
            }
        }
        else
        {
            is.putBack(t);
        }

        scalarField data(is);
        return data;
    }

#ifdef FOAM_USE_ZSTD
    // 2a. Delta-varint + zstd format (.dvz.zstd)
    filePath = dataPath / (pressureFieldName_ + "."
        + deltaVarintCodec::fileExtension() + "." + zstdWrapper::fileExtension);
    if (isFile(filePath))
    {
        std::vector<uint8_t> buf = zstdWrapper::decompressFromFile(filePath);
        return deltaVarintCodec::decodeScalar(buf);
    }

    // 2b. Delta-varint-temporal + zstd format (.dvzt.zstd)
    filePath = dataPath / (pressureFieldName_ + "."
        + deltaVarintTemporalCodec::fileExtension() + "." + zstdWrapper::fileExtension);
    if (isFile(filePath))
    {
        std::vector<uint8_t> buf = zstdWrapper::decompressFromFile(filePath);
        return deltaVarintTemporalCodec::decodeScalar(buf);
    }
#endif

    // 3. Delta-varint format (.dvz)
    filePath = dataPath / (pressureFieldName_ + "." + deltaVarintCodec::fileExtension());
    if (isFile(filePath))
    {
        return deltaVarintCodec::readScalar(filePath);
    }

    // 4. Delta-varint-temporal format (.dvzt)
    filePath = dataPath / (pressureFieldName_ + "." + deltaVarintTemporalCodec::fileExtension());
    if (isFile(filePath))
    {
        return deltaVarintTemporalCodec::readScalar(filePath);
    }

    FatalErrorInFunction
        << "Cannot find pressure data file for time " << timeDir.name() << nl
        << "Searched for: " << pressureFieldName_ << ", "
        << pressureFieldName_ + "." + deltaVarintCodec::fileExtension() + "." + zstdWrapper::fileExtension << ", "
        << pressureFieldName_ + "." + deltaVarintTemporalCodec::fileExtension() + "." + zstdWrapper::fileExtension << ", "
        << pressureFieldName_ + "." + deltaVarintCodec::fileExtension() << ", "
        << pressureFieldName_ + "." + deltaVarintTemporalCodec::fileExtension() << nl
        << "in directory: " << dataPath << nl
        << exit(FatalError);

    return scalarField();
}


Foam::scalarField
Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::readGradientData
(
    const instant& timeDir
) const
{
    fileName dataPath = dataDir_ / patch().name() / timeDir.name();

    // Try different file formats
    // 1. Standard OpenFOAM format
    fileName filePath = dataPath / gradientFieldName_;
    if (isFile(filePath))
    {
        IFstream is(filePath);
        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot open file " << filePath << nl
                << exit(FatalError);
        }

        // Parse FoamFile header if present (handles binary format)
        token t(is);
        if (t.isWord() && t.wordToken() == "FoamFile")
        {
            dictionary headerDict;
            headerDict.read(is, false);
            if (headerDict.found("format"))
            {
                const word formatStr = headerDict.get<word>("format");
                if (formatStr == "binary")
                {
                    is.format(IOstreamOption::BINARY);
                }
            }
        }
        else
        {
            is.putBack(t);
        }

        scalarField data(is);
        return data;
    }

#ifdef FOAM_USE_ZSTD
    // 2a. Delta-varint + zstd format (.dvz.zstd)
    filePath = dataPath / (gradientFieldName_ + "."
        + deltaVarintCodec::fileExtension() + "." + zstdWrapper::fileExtension);
    if (isFile(filePath))
    {
        std::vector<uint8_t> buf = zstdWrapper::decompressFromFile(filePath);
        return deltaVarintCodec::decodeScalar(buf);
    }

    // 2b. Delta-varint-temporal + zstd format (.dvzt.zstd)
    filePath = dataPath / (gradientFieldName_ + "."
        + deltaVarintTemporalCodec::fileExtension() + "." + zstdWrapper::fileExtension);
    if (isFile(filePath))
    {
        std::vector<uint8_t> buf = zstdWrapper::decompressFromFile(filePath);
        return deltaVarintTemporalCodec::decodeScalar(buf);
    }
#endif

    // 3. Delta-varint format (.dvz)
    filePath = dataPath / (gradientFieldName_ + "." + deltaVarintCodec::fileExtension());
    if (isFile(filePath))
    {
        return deltaVarintCodec::readScalar(filePath);
    }

    // 4. Delta-varint-temporal format (.dvzt)
    filePath = dataPath / (gradientFieldName_ + "." + deltaVarintTemporalCodec::fileExtension());
    if (isFile(filePath))
    {
        return deltaVarintTemporalCodec::readScalar(filePath);
    }

    FatalErrorInFunction
        << "Cannot find gradient data file for time " << timeDir.name() << nl
        << "Searched for: " << gradientFieldName_ << ", "
        << gradientFieldName_ + "." + deltaVarintCodec::fileExtension() + "." + zstdWrapper::fileExtension << ", "
        << gradientFieldName_ + "." + deltaVarintTemporalCodec::fileExtension() + "." + zstdWrapper::fileExtension << ", "
        << gradientFieldName_ + "." + deltaVarintCodec::fileExtension() << ", "
        << gradientFieldName_ + "." + deltaVarintTemporalCodec::fileExtension() << nl
        << "in directory: " << dataPath << nl
        << exit(FatalError);

    return scalarField();
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::checkTable()
{
    const scalar t = this->db().time().value();

    // Find sample times if not already done
    findSampleTimes();

    // Check if we need to update cached data
    // Find bracketing times for current t
    label newIndex1 = -1;
    label newIndex2 = -1;

    for (label i = 0; i < sampleTimes_.size() - 1; ++i)
    {
        if (t >= sampleTimes_[i].value() && t <= sampleTimes_[i+1].value())
        {
            newIndex1 = i;
            newIndex2 = i + 1;
            break;
        }
    }

    // Handle exact matches at boundaries
    if (newIndex1 < 0)
    {
        if (mag(t - sampleTimes_.first().value()) < SMALL)
        {
            newIndex1 = 0;
            newIndex2 = 0;
        }
        else if (mag(t - sampleTimes_.last().value()) < SMALL)
        {
            newIndex1 = sampleTimes_.size() - 1;
            newIndex2 = sampleTimes_.size() - 1;
        }
        else
        {
            FatalErrorInFunction
                << "Time " << t << " is outside the extraction window ["
                << sampleTimes_.first().value() << ", "
                << sampleTimes_.last().value() << "]" << nl
                << "Extrapolation is not allowed." << nl
                << exit(FatalError);
        }
    }

    // For cubic interpolation, we need 4 points
    label newIndex0 = max(0, newIndex1 - 1);
    label newIndex3 = min(sampleTimes_.size() - 1, newIndex2 + 1);

    // Check if we need to reload data
    bool needReload = (newIndex1 != sampleTimeIndex1_ || newIndex2 != sampleTimeIndex2_);

    if (needReload)
    {
        // Read pressure data
        if (timeInterpolationScheme_ == timeInterpScheme::CUBIC && sampleTimes_.size() >= 4)
        {
            sampledPressure0_ = readPressureData(sampleTimes_[newIndex0]);
            sampledGradient0_ = readGradientData(sampleTimes_[newIndex0]);
        }
        sampledPressure1_ = readPressureData(sampleTimes_[newIndex1]);
        sampledGradient1_ = readGradientData(sampleTimes_[newIndex1]);

        if (newIndex2 != newIndex1)
        {
            sampledPressure2_ = readPressureData(sampleTimes_[newIndex2]);
            sampledGradient2_ = readGradientData(sampleTimes_[newIndex2]);
        }
        else
        {
            sampledPressure2_ = sampledPressure1_;
            sampledGradient2_ = sampledGradient1_;
        }

        if (timeInterpolationScheme_ == timeInterpScheme::CUBIC && sampleTimes_.size() >= 4)
        {
            sampledPressure3_ = readPressureData(sampleTimes_[newIndex3]);
            sampledGradient3_ = readGradientData(sampleTimes_[newIndex3]);
        }

        sampleTimeIndex0_ = newIndex0;
        sampleTimeIndex1_ = newIndex1;
        sampleTimeIndex2_ = newIndex2;
        sampleTimeIndex3_ = newIndex3;
    }
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::validateMetadata() const
{
    if (metadataValidated_)
    {
        return;
    }

    // Try to read extraction metadata
    fileName metadataPath = dataDir_ / patch().name() / "extractionMetadata";

    if (isFile(metadataPath))
    {
        IFstream is(metadataPath);
        if (is.good())
        {
            dictionary metaDict(is);

            // Validate deltaT match
            scalar extractedDeltaT = metaDict.getOrDefault<scalar>("deltaT", 0);
            scalar currentDeltaT = this->db().time().deltaTValue();

            if (extractedDeltaT > 0 && mag(extractedDeltaT - currentDeltaT) > 1e-10)
            {
                WarningInFunction
                    << "Time step mismatch: extraction used deltaT = " << extractedDeltaT
                    << ", current simulation uses deltaT = " << currentDeltaT << nl
                    << "This may cause interpolation artifacts" << endl;
            }
        }
    }

    metadataValidated_ = true;
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::initSpatialInterpolation() const
{
    if (spatialInterpInitialized_)
    {
        return;
    }

    // Read source points from boundaryData
    fileName pointsPath = dataDir_ / patch().name() / "points";

    if (isFile(pointsPath))
    {
        IFstream is(pointsPath);
        if (is.good())
        {
            // Skip FoamFile header if present
            token t(is);
            if (t.isWord() && t.wordToken() == "FoamFile")
            {
                dictionary headerDict(is);
            }
            else
            {
                is.putBack(t);
            }

            vectorField sourcePoints(is);

            // Get target face centers
            const vectorField& targetPoints = patch().Cf();

            // Check if spatial interpolation is needed
            if (sourcePoints.size() != targetPoints.size())
            {
                // Compute bounding box from source points
                boundBox bb(sourcePoints);
                boundingBox_ = bb;

                // Initialize spatial interpolation
                spatialInterpPtr_.reset
                (
                    new spatialInterpolation2D(sourcePoints, targetPoints, bb)
                );

                Info<< "spaceTimeWindowCoupledPressure: Spatial interpolation enabled" << nl
                    << "    Source points: " << sourcePoints.size() << nl
                    << "    Target points: " << targetPoints.size() << endl;
            }
        }
    }

    spatialInterpInitialized_ = true;
}


Foam::tmp<Foam::scalarField>
Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::computeBlendingWeights() const
{
    // Get flux field
    const surfaceScalarField& phi =
        this->db().lookupObject<surfaceScalarField>(phiName_);

    // Get patch flux
    const scalarField& phip = phi.boundaryField()[patch().index()];

    // Compute blending weight for each face
    auto tweights = tmp<scalarField>::New(patch().size());
    scalarField& weights = tweights.ref();

    // Get reference flux magnitude for normalization
    scalar phiMax = max(mag(phip));
    reduce(phiMax, maxOp<scalar>());
    phiMax = max(phiMax, SMALL);

    if (blendingMode_ == blendingMode::FLUX_MAGNITUDE)
    {
        // Flux magnitude mode (recommended):
        // Dirichlet at stagnant faces (zero flux), Neumann at active flow faces
        forAll(weights, facei)
        {
            // Normalized flux magnitude (0 = no flow, 1 = max flow)
            scalar phiMagNorm = mag(phip[facei]) / phiMax;

            // Smooth transition based on flux magnitude using tanh
            scalar transition = Foam::tanh(phiMagNorm / transitionWidth_);

            // At zero flux: weight = dirichletWeight_ (high Dirichlet)
            // At high flux: weight = 1 - neumannWeight_ (low Dirichlet / high Neumann)
            scalar wDir = dirichletWeight_;
            scalar wNeu = 1.0 - neumannWeight_;

            weights[facei] = wDir + (wNeu - wDir) * transition;
            weights[facei] = max(0.0, min(1.0, weights[facei]));
        }
    }
    else // blendingMode::FLOW_DIRECTION
    {
        // Flow direction mode (original):
        // Dirichlet at inflow, Neumann at outflow
        forAll(weights, facei)
        {
            // Normalized flux: positive = outflow, negative = inflow
            scalar phiNorm = phip[facei] / phiMax;

            // Smooth transition using tanh
            // At strong inflow (phi < -transitionWidth_): weight -> dirichletWeight_
            // At strong outflow (phi > transitionWidth_): weight -> 1 - neumannWeight_
            scalar transition = Foam::tanh(phiNorm / transitionWidth_);

            // Map from [-1, 1] to [dirichletWeight_, 1-neumannWeight_]
            scalar wDir = dirichletWeight_;        // Weight at inflow
            scalar wNeu = 1.0 - neumannWeight_;    // Weight at outflow

            weights[facei] = 0.5 * ((wDir + wNeu) - (wDir - wNeu) * transition);
            weights[facei] = max(0.0, min(1.0, weights[facei]));
        }
    }

    return tweights;
}


Foam::scalar
Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::catmullRomCentripetal
(
    scalar p0,
    scalar p1,
    scalar p2,
    scalar p3,
    scalar t0,
    scalar t1,
    scalar t2,
    scalar t3,
    scalar t
) const
{
    // Centripetal Catmull-Rom spline interpolation
    // Handles non-uniform time spacing correctly

    // Compute segment lengths using centripetal parameterization
    // For scalar time values, this reduces to simple differences
    scalar dt01 = Foam::sqrt(mag(t1 - t0));
    scalar dt12 = Foam::sqrt(mag(t2 - t1));
    scalar dt23 = Foam::sqrt(mag(t3 - t2));

    // Avoid division by zero
    dt01 = max(dt01, SMALL);
    dt12 = max(dt12, SMALL);
    dt23 = max(dt23, SMALL);

    // Compute tangents
    scalar m1 = (p2 - p0) / (dt01 + dt12) * dt12;
    scalar m2 = (p3 - p1) / (dt12 + dt23) * dt12;

    // Compute interpolation parameter
    scalar s = (t - t1) / (t2 - t1);

    // Hermite basis functions
    scalar s2 = s * s;
    scalar s3 = s2 * s;

    scalar h00 = 2*s3 - 3*s2 + 1;
    scalar h10 = s3 - 2*s2 + s;
    scalar h01 = -2*s3 + 3*s2;
    scalar h11 = s3 - s2;

    // Interpolate
    return h00*p1 + h10*m1 + h01*p2 + h11*m2;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::
spaceTimeWindowCoupledPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    dataDir_(this->db().time().path()/"constant"/"boundaryData"),
    pressureFieldName_("p"),
    gradientFieldName_("gradp"),
    phiName_("phi"),
    dirichletWeight_(0.8),
    neumannWeight_(0.8),
    transitionWidth_(0.1),
    blendingMode_(blendingMode::FLUX_MAGNITUDE),
    allowTimeInterpolation_(true),
    timeInterpolationScheme_(timeInterpScheme::CUBIC),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledPressure0_(),
    sampledPressure1_(),
    sampledPressure2_(),
    sampledPressure3_(),
    sampledGradient0_(),
    sampledGradient1_(),
    sampledGradient2_(),
    sampledGradient3_(),
    metadataValidated_(false),
    spatialInterpPtr_(),
    boundingBox_(),
    spatialInterpInitialized_(false)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0.5;
}


Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::
spaceTimeWindowCoupledPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    dataDir_
    (
        dict.getOrDefault<fileName>
        (
            "dataDir",
            this->db().time().path()/"constant"/"boundaryData"
        )
    ),
    pressureFieldName_(dict.getOrDefault<word>("pressureField", "p")),
    gradientFieldName_(dict.getOrDefault<word>("gradientField", "gradp")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    dirichletWeight_(dict.getOrDefault<scalar>("dirichletWeight", 0.8)),
    neumannWeight_(dict.getOrDefault<scalar>("neumannWeight", 0.8)),
    transitionWidth_(dict.getOrDefault<scalar>("transitionWidth", 0.1)),
    blendingMode_(blendingMode::FLUX_MAGNITUDE),
    allowTimeInterpolation_(dict.getOrDefault<bool>("allowTimeInterpolation", true)),
    timeInterpolationScheme_(timeInterpScheme::CUBIC),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledPressure0_(),
    sampledPressure1_(),
    sampledPressure2_(),
    sampledPressure3_(),
    sampledGradient0_(),
    sampledGradient1_(),
    sampledGradient2_(),
    sampledGradient3_(),
    metadataValidated_(false),
    spatialInterpPtr_(),
    boundingBox_(),
    spatialInterpInitialized_(false)
{
    // Read blending mode
    const word modeStr = dict.getOrDefault<word>("blendingMode", "fluxMagnitude");
    if (modeStr == "flowDirection")
    {
        blendingMode_ = blendingMode::FLOW_DIRECTION;
    }
    else if (modeStr == "fluxMagnitude")
    {
        blendingMode_ = blendingMode::FLUX_MAGNITUDE;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown blending mode: " << modeStr << nl
            << "Valid modes are: flowDirection, fluxMagnitude" << nl
            << exit(FatalIOError);
    }

    // Read time interpolation scheme
    const word schemeStr = dict.getOrDefault<word>("timeInterpolationScheme", "cubic");
    if (schemeStr == "none")
    {
        timeInterpolationScheme_ = timeInterpScheme::NONE;
    }
    else if (schemeStr == "linear")
    {
        timeInterpolationScheme_ = timeInterpScheme::LINEAR;
    }
    else if (schemeStr == "cubic")
    {
        timeInterpolationScheme_ = timeInterpScheme::CUBIC;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Unknown time interpolation scheme: " << schemeStr << nl
            << "Valid schemes are: none, linear, cubic" << nl
            << exit(FatalIOError);
    }

    // Initialize mixed BC values
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = *this;
    this->refGrad() = 0;
    this->valueFraction() = 0.5;

    if (dict.found("refValue"))
    {
        this->refValue() = scalarField("refValue", dict, p.size());
    }

    if (dict.found("refGradient"))
    {
        this->refGrad() = scalarField("refGradient", dict, p.size());
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() = scalarField("valueFraction", dict, p.size());
    }
}


Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::
spaceTimeWindowCoupledPressureFvPatchScalarField
(
    const spaceTimeWindowCoupledPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    dataDir_(ptf.dataDir_),
    pressureFieldName_(ptf.pressureFieldName_),
    gradientFieldName_(ptf.gradientFieldName_),
    phiName_(ptf.phiName_),
    dirichletWeight_(ptf.dirichletWeight_),
    neumannWeight_(ptf.neumannWeight_),
    transitionWidth_(ptf.transitionWidth_),
    blendingMode_(ptf.blendingMode_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledPressure0_(),
    sampledPressure1_(),
    sampledPressure2_(),
    sampledPressure3_(),
    sampledGradient0_(),
    sampledGradient1_(),
    sampledGradient2_(),
    sampledGradient3_(),
    metadataValidated_(ptf.metadataValidated_),
    spatialInterpPtr_(),
    boundingBox_(ptf.boundingBox_),
    spatialInterpInitialized_(false)
{}


Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::
spaceTimeWindowCoupledPressureFvPatchScalarField
(
    const spaceTimeWindowCoupledPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    dataDir_(ptf.dataDir_),
    pressureFieldName_(ptf.pressureFieldName_),
    gradientFieldName_(ptf.gradientFieldName_),
    phiName_(ptf.phiName_),
    dirichletWeight_(ptf.dirichletWeight_),
    neumannWeight_(ptf.neumannWeight_),
    transitionWidth_(ptf.transitionWidth_),
    blendingMode_(ptf.blendingMode_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledPressure0_(ptf.sampledPressure0_),
    sampledPressure1_(ptf.sampledPressure1_),
    sampledPressure2_(ptf.sampledPressure2_),
    sampledPressure3_(ptf.sampledPressure3_),
    sampledGradient0_(ptf.sampledGradient0_),
    sampledGradient1_(ptf.sampledGradient1_),
    sampledGradient2_(ptf.sampledGradient2_),
    sampledGradient3_(ptf.sampledGradient3_),
    metadataValidated_(ptf.metadataValidated_),
    spatialInterpPtr_(),
    boundingBox_(ptf.boundingBox_),
    spatialInterpInitialized_(ptf.spatialInterpInitialized_)
{}


Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::
spaceTimeWindowCoupledPressureFvPatchScalarField
(
    const spaceTimeWindowCoupledPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    dataDir_(ptf.dataDir_),
    pressureFieldName_(ptf.pressureFieldName_),
    gradientFieldName_(ptf.gradientFieldName_),
    phiName_(ptf.phiName_),
    dirichletWeight_(ptf.dirichletWeight_),
    neumannWeight_(ptf.neumannWeight_),
    transitionWidth_(ptf.transitionWidth_),
    blendingMode_(ptf.blendingMode_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledPressure0_(ptf.sampledPressure0_),
    sampledPressure1_(ptf.sampledPressure1_),
    sampledPressure2_(ptf.sampledPressure2_),
    sampledPressure3_(ptf.sampledPressure3_),
    sampledGradient0_(ptf.sampledGradient0_),
    sampledGradient1_(ptf.sampledGradient1_),
    sampledGradient2_(ptf.sampledGradient2_),
    sampledGradient3_(ptf.sampledGradient3_),
    metadataValidated_(ptf.metadataValidated_),
    spatialInterpPtr_(),
    boundingBox_(ptf.boundingBox_),
    spatialInterpInitialized_(ptf.spatialInterpInitialized_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    // Clear cached data - will be reloaded
    sampledPressure0_.clear();
    sampledPressure1_.clear();
    sampledPressure2_.clear();
    sampledPressure3_.clear();
    sampledGradient0_.clear();
    sampledGradient1_.clear();
    sampledGradient2_.clear();
    sampledGradient3_.clear();

    spatialInterpInitialized_ = false;
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    // Unused but kept for potential future use
    // const spaceTimeWindowCoupledPressureFvPatchScalarField& tiptf =
    //     refCast<const spaceTimeWindowCoupledPressureFvPatchScalarField>(ptf);

    // Clear cached data
    sampledPressure0_.clear();
    sampledPressure1_.clear();
    sampledPressure2_.clear();
    sampledPressure3_.clear();
    sampledGradient0_.clear();
    sampledGradient1_.clear();
    sampledGradient2_.clear();
    sampledGradient3_.clear();

    spatialInterpInitialized_ = false;
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Validate metadata on first call
    validateMetadata();

    // Initialize spatial interpolation if needed
    initSpatialInterpolation();

    // Update sample data
    checkTable();

    const scalar t = this->db().time().value();

    // Get interpolation times
    scalar t1 = sampleTimes_[sampleTimeIndex1_].value();
    scalar t2 = sampleTimes_[sampleTimeIndex2_].value();

    // Interpolate pressure values
    scalarField interpolatedPressure(patch().size());
    scalarField interpolatedGradient(patch().size());

    if (timeInterpolationScheme_ == timeInterpScheme::CUBIC
        && sampleTimes_.size() >= 4
        && sampleTimeIndex0_ >= 0
        && sampleTimeIndex3_ < sampleTimes_.size())
    {
        scalar t0 = sampleTimes_[sampleTimeIndex0_].value();
        scalar t3 = sampleTimes_[sampleTimeIndex3_].value();

        forAll(interpolatedPressure, facei)
        {
            interpolatedPressure[facei] = catmullRomCentripetal
            (
                sampledPressure0_[facei],
                sampledPressure1_[facei],
                sampledPressure2_[facei],
                sampledPressure3_[facei],
                t0, t1, t2, t3, t
            );

            interpolatedGradient[facei] = catmullRomCentripetal
            (
                sampledGradient0_[facei],
                sampledGradient1_[facei],
                sampledGradient2_[facei],
                sampledGradient3_[facei],
                t0, t1, t2, t3, t
            );
        }
    }
    else
    {
        // Linear interpolation
        scalar alpha = 0;
        if (mag(t2 - t1) > SMALL)
        {
            alpha = (t - t1) / (t2 - t1);
        }

        interpolatedPressure = (1 - alpha) * sampledPressure1_ + alpha * sampledPressure2_;
        interpolatedGradient = (1 - alpha) * sampledGradient1_ + alpha * sampledGradient2_;
    }

    // Apply spatial interpolation if needed
    if (spatialInterpPtr_.valid())
    {
        interpolatedPressure = spatialInterpPtr_->interpolate(interpolatedPressure);
        interpolatedGradient = spatialInterpPtr_->interpolate(interpolatedGradient);
    }

    // Compute blending weights based on flux direction
    tmp<scalarField> tweights = computeBlendingWeights();
    const scalarField& weights = tweights();

    // Set mixed BC coefficients
    // refValue = Dirichlet value (extracted pressure)
    // refGrad = Neumann value (extracted gradient)
    // valueFraction = blending weight (1 = pure Dirichlet, 0 = pure Neumann)

    this->refValue() = interpolatedPressure;
    this->refGrad() = interpolatedGradient;
    this->valueFraction() = weights;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::spaceTimeWindowCoupledPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);

    os.writeEntry("dataDir", dataDir_);
    os.writeEntry("pressureField", pressureFieldName_);
    os.writeEntry("gradientField", gradientFieldName_);
    os.writeEntry("phi", phiName_);

    // Write blending mode
    word modeStr;
    switch (blendingMode_)
    {
        case blendingMode::FLOW_DIRECTION:
            modeStr = "flowDirection";
            break;
        case blendingMode::FLUX_MAGNITUDE:
            modeStr = "fluxMagnitude";
            break;
    }
    os.writeEntry("blendingMode", modeStr);

    os.writeEntry("dirichletWeight", dirichletWeight_);
    os.writeEntry("neumannWeight", neumannWeight_);
    os.writeEntry("transitionWidth", transitionWidth_);
    os.writeEntry("allowTimeInterpolation", allowTimeInterpolation_);

    word schemeStr;
    switch (timeInterpolationScheme_)
    {
        case timeInterpScheme::NONE:
            schemeStr = "none";
            break;
        case timeInterpScheme::LINEAR:
            schemeStr = "linear";
            break;
        case timeInterpScheme::CUBIC:
            schemeStr = "cubic";
            break;
    }
    os.writeEntry("timeInterpolationScheme", schemeStr);
}


// ************************************************************************* //
