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
\*---------------------------------------------------------------------------*/

#include "spaceTimeWindowFvPatchField.H"
#include "Time.H"
#include "IFstream.H"
#include "dictionary.H"
#include "foamVersion.H"
#include "token.H"
#include "IOstreamOption.H"
#include "deltaVarintCodec.H"
#include "pTraits.H"
#include "volFields.H"
#include <type_traits>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::spaceTimeWindowFvPatchField<Type>::findSampleTimes() const
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
Foam::Field<Type> Foam::spaceTimeWindowFvPatchField<Type>::readFieldData
(
    const instant& timeDir
) const
{
    const fileName baseDir
    (
        this->db().time().globalPath()
      / dataDir_
      / this->patch().name()
      / timeDir.name()
    );

    // Check for delta-varint compressed file first (.dvz extension)
    const fileName dvzPath = baseDir / (fieldTableName_ + "." + deltaVarintCodec::fileExtension());

    if (isFile(dvzPath))
    {
        // Read using delta-varint codec
        // Use type traits to determine scalar vs vector
        Field<Type> data;

        if constexpr (std::is_same_v<Type, scalar>)
        {
            data = deltaVarintCodec::readScalar(dvzPath);
        }
        else if constexpr (std::is_same_v<Type, vector>)
        {
            data = deltaVarintCodec::readVector(dvzPath);
        }
        else
        {
            FatalErrorInFunction
                << "Delta-varint codec only supports scalar and vector fields" << nl
                << "    File: " << dvzPath << nl
                << "    Type: " << pTraits<Type>::typeName << nl
                << exit(FatalError);
        }

        // Verify size matches patch
        if (data.size() != this->patch().size())
        {
            FatalErrorInFunction
                << "Field size " << data.size()
                << " does not match patch size " << this->patch().size() << nl
                << "    File: " << dvzPath << nl
                << "    Patch: " << this->patch().name() << nl
                << "    This BC requires pre-computed face values with NO spatial mapping" << nl
                << exit(FatalError);
        }

        return data;
    }

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

    // Check if file has FoamFile header (supports both ASCII and binary formats)
    // Peek at first token to detect header
    token firstToken(is);

    if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
    {
        // File has FoamFile header - read header dictionary content
        dictionary headerDict;
        headerDict.read(is, false);  // Read sub-dictionary content (without leading '{')

        // Parse format from header and set stream format (critical for binary)
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

        // Read field data (format now determined by header)
        Field<Type> data(is);

        // Verify size matches patch
        if (data.size() != this->patch().size())
        {
            FatalErrorInFunction
                << "Field size " << data.size()
                << " does not match patch size " << this->patch().size() << nl
                << "    File: " << dataPath << nl
                << "    Patch: " << this->patch().name() << nl
                << "    This BC requires pre-computed face values with NO spatial mapping" << nl
                << exit(FatalError);
        }

        return data;
    }

    // No header - put token back and read raw field data (backwards compatibility)
    is.putBack(firstToken);
    Field<Type> data(is);

    // Verify size matches patch
    if (data.size() != this->patch().size())
    {
        FatalErrorInFunction
            << "Field size " << data.size()
            << " does not match patch size " << this->patch().size() << nl
            << "    File: " << dataPath << nl
            << "    Patch: " << this->patch().name() << nl
            << "    This BC requires pre-computed face values with NO spatial mapping" << nl
            << exit(FatalError);
    }

    return data;
}


template<class Type>
void Foam::spaceTimeWindowFvPatchField<Type>::checkTable()
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
void Foam::spaceTimeWindowFvPatchField<Type>::validateMetadata() const
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

    // Validate OpenFOAM version (warning only, not fatal)
    if (metadata.found("openfoamApi"))
    {
        const label extractionApi = metadata.get<label>("openfoamApi");
        const label currentApi = foamVersion::api;

        if (extractionApi != currentApi)
        {
            WarningInFunction
                << "OpenFOAM API version mismatch between extraction and reconstruction!" << nl
                << "    Extraction openfoamApi: " << extractionApi << nl
                << "    Current openfoamApi:    " << currentApi << nl
                << "    Patch: " << this->patch().name() << nl
                << "    This may cause compatibility issues." << endl;
        }
    }

    // Validate solver name (warning only, not fatal)
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
                << "Solver mismatch between extraction and reconstruction!" << nl
                << "    Extraction solver: " << extractionSolver << nl
                << "    Current solver:    " << currentSolver << nl
                << "    Patch: " << this->patch().name() << nl
                << "    Ensure the same solver is used for accuracy." << endl;
        }
    }

    // Validate deltaT
    if (metadata.found("deltaT"))
    {
        const scalar extractionDeltaT = metadata.get<scalar>("deltaT");
        const scalar currentDeltaT = runTime.deltaTValue();

        if (mag(extractionDeltaT - currentDeltaT) > SMALL * extractionDeltaT)
        {
            FatalErrorInFunction
                << "Time step mismatch between extraction and reconstruction!" << nl
                << "    Extraction deltaT: " << extractionDeltaT << nl
                << "    Current deltaT:    " << currentDeltaT << nl
                << "    Patch: " << this->patch().name() << nl << nl
                << "    The spaceTimeWindow BC requires identical time settings" << nl
                << "    to avoid temporal interpolation/extrapolation errors." << nl
                << "    Please set deltaT = " << extractionDeltaT
                << " in system/controlDict" << nl
                << exit(FatalError);
        }
    }

    // Validate adjustTimeStep
    if (metadata.found("adjustTimeStep"))
    {
        const bool extractionAdjustTimeStep = metadata.get<bool>("adjustTimeStep");
        const bool currentAdjustTimeStep =
            runTime.controlDict().getOrDefault<Switch>("adjustTimeStep", false);

        if (extractionAdjustTimeStep != currentAdjustTimeStep)
        {
            FatalErrorInFunction
                << "adjustTimeStep mismatch between extraction and reconstruction!" << nl
                << "    Extraction adjustTimeStep: " << extractionAdjustTimeStep << nl
                << "    Current adjustTimeStep:    " << currentAdjustTimeStep << nl
                << "    Patch: " << this->patch().name() << nl << nl
                << "    The spaceTimeWindow BC requires identical time settings." << nl
                << "    Please set adjustTimeStep = "
                << (extractionAdjustTimeStep ? "yes" : "no")
                << " in system/controlDict" << nl
                << exit(FatalError);
        }

        // If adjustTimeStep is enabled, warn that this may cause issues
        if (extractionAdjustTimeStep)
        {
            WarningInFunction
                << "Extraction was done with adjustTimeStep = yes" << nl
                << "    This may cause temporal interpolation if time steps don't match exactly." << nl
                << "    For best accuracy, use fixed time step (adjustTimeStep = no)." << endl;
        }
    }

    // Report successful validation
    Info<< "spaceTimeWindow BC on patch " << this->patch().name()
        << ": validated extraction metadata" << nl;

    if (metadata.found("openfoamVersion"))
    {
        Info<< "    openfoamVersion = " << metadata.get<word>("openfoamVersion") << nl;
    }
    if (metadata.found("solver"))
    {
        Info<< "    solver = " << metadata.get<word>("solver") << nl;
    }
    Info<< "    deltaT = " << runTime.deltaTValue() << endl;

    metadataValidated_ = true;
}


template<class Type>
Type Foam::spaceTimeWindowFvPatchField<Type>::catmullRomCentripetal
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
    // Handles non-uniform time spacing correctly (important for adaptive dt)
    //
    // Reference: Barry & Goldman (1988) "A recursive evaluation algorithm
    //            for a class of Catmull-Rom splines"
    //
    // The centripetal parameterization uses alpha=0.5 for the knot sequence,
    // which prevents cusps and self-intersections with non-uniform data.
    //
    // Knot parameterization with alpha=0.5 (centripetal):
    //   For 1D time data, |t_{i+1} - t_i|^0.5 simplifies since time is monotonic

    // Compute knot intervals using centripetal parameterization
    // For time values (1D), the "distance" is just the time difference
    const scalar dt01 = Foam::sqrt(Foam::mag(t1 - t0));
    const scalar dt12 = Foam::sqrt(Foam::mag(t2 - t1));
    const scalar dt23 = Foam::sqrt(Foam::mag(t3 - t2));

    // Avoid division by zero for duplicate times
    const scalar eps = SMALL;

    // Knot values (cumulative)
    const scalar k0 = 0.0;
    const scalar k1 = k0 + max(dt01, eps);
    const scalar k2 = k1 + max(dt12, eps);
    const scalar k3 = k2 + max(dt23, eps);

    // Map input time t to knot parameter u
    // t is in [t1, t2], map to u in [k1, k2]
    const scalar alpha = (t - t1) / (t2 - t1 + eps);
    const scalar u = k1 + alpha * (k2 - k1);

    // Barry-Goldman pyramid algorithm for Catmull-Rom evaluation
    // First level interpolations
    const Type A1 = (k1 - u)/(k1 - k0 + eps) * p0 + (u - k0)/(k1 - k0 + eps) * p1;
    const Type A2 = (k2 - u)/(k2 - k1 + eps) * p1 + (u - k1)/(k2 - k1 + eps) * p2;
    const Type A3 = (k3 - u)/(k3 - k2 + eps) * p2 + (u - k2)/(k3 - k2 + eps) * p3;

    // Second level interpolations
    const Type B1 = (k2 - u)/(k2 - k0 + eps) * A1 + (u - k0)/(k2 - k0 + eps) * A2;
    const Type B2 = (k3 - u)/(k3 - k1 + eps) * A2 + (u - k1)/(k3 - k1 + eps) * A3;

    // Final interpolation
    const Type C = (k2 - u)/(k2 - k1 + eps) * B1 + (u - k1)/(k2 - k1 + eps) * B2;

    return C;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::spaceTimeWindowFvPatchField<Type>::spaceTimeWindowFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    dataDir_("constant/boundaryData"),
    fieldTableName_(iF.name()),
    setAverage_(false),
    offset_(Zero),
    fixesValue_(true),
    allowTimeInterpolation_(false),
    timeInterpolationScheme_(timeInterpScheme::LINEAR),
    reportFlux_(false),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(false)
{}


template<class Type>
Foam::spaceTimeWindowFvPatchField<Type>::spaceTimeWindowFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, IOobjectOption::MUST_READ),
    dataDir_(dict.getOrDefault<fileName>("dataDir", "constant/boundaryData")),
    fieldTableName_(dict.getOrDefault<word>("fieldTableName", iF.name())),
    setAverage_(dict.getOrDefault("setAverage", false)),
    offset_(dict.getOrDefault<Type>("offset", Zero)),
    fixesValue_(dict.getOrDefault("fixesValue", true)),
    allowTimeInterpolation_(dict.getOrDefault("allowTimeInterpolation", false)),
    timeInterpolationScheme_(timeInterpScheme::LINEAR),
    reportFlux_(dict.getOrDefault("reportFlux", false)),
    sampleTimes_(),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(false)
{
    // Parse time interpolation scheme
    const word schemeStr = dict.getOrDefault<word>("timeInterpolationScheme", "linear");
    if (schemeStr == "cubic")
    {
        timeInterpolationScheme_ = timeInterpScheme::CUBIC;
    }
    else if (schemeStr != "linear")
    {
        FatalIOErrorInFunction(dict)
            << "Unknown timeInterpolationScheme: " << schemeStr << nl
            << "    Valid options: linear, cubic" << nl
            << exit(FatalIOError);
    }

    // Note: value is already read by parent fixedValueFvPatchField constructor
    // with IOobjectOption::MUST_READ - do NOT read it again here as compound
    // tokens can only be transferred once
}


template<class Type>
Foam::spaceTimeWindowFvPatchField<Type>::spaceTimeWindowFvPatchField
(
    const spaceTimeWindowFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    offset_(ptf.offset_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    reportFlux_(ptf.reportFlux_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(-1),
    sampleTimeIndex1_(-1),
    sampleTimeIndex2_(-1),
    sampleTimeIndex3_(-1),
    sampledValues0_(),
    sampledValues1_(),
    sampledValues2_(),
    sampledValues3_(),
    metadataValidated_(ptf.metadataValidated_)
{}


template<class Type>
Foam::spaceTimeWindowFvPatchField<Type>::spaceTimeWindowFvPatchField
(
    const spaceTimeWindowFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    offset_(ptf.offset_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    reportFlux_(ptf.reportFlux_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledValues0_(ptf.sampledValues0_),
    sampledValues1_(ptf.sampledValues1_),
    sampledValues2_(ptf.sampledValues2_),
    sampledValues3_(ptf.sampledValues3_),
    metadataValidated_(ptf.metadataValidated_)
{}


template<class Type>
Foam::spaceTimeWindowFvPatchField<Type>::spaceTimeWindowFvPatchField
(
    const spaceTimeWindowFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    dataDir_(ptf.dataDir_),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    offset_(ptf.offset_),
    fixesValue_(ptf.fixesValue_),
    allowTimeInterpolation_(ptf.allowTimeInterpolation_),
    timeInterpolationScheme_(ptf.timeInterpolationScheme_),
    reportFlux_(ptf.reportFlux_),
    sampleTimes_(ptf.sampleTimes_),
    sampleTimeIndex0_(ptf.sampleTimeIndex0_),
    sampleTimeIndex1_(ptf.sampleTimeIndex1_),
    sampleTimeIndex2_(ptf.sampleTimeIndex2_),
    sampleTimeIndex3_(ptf.sampleTimeIndex3_),
    sampledValues0_(ptf.sampledValues0_),
    sampledValues1_(ptf.sampledValues1_),
    sampledValues2_(ptf.sampledValues2_),
    sampledValues3_(ptf.sampledValues3_),
    metadataValidated_(ptf.metadataValidated_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::spaceTimeWindowFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);

    // Clear cached data - will reload on next updateCoeffs
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
void Foam::spaceTimeWindowFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

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
void Foam::spaceTimeWindowFvPatchField<Type>::updateCoeffs()
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

    // Compute interpolated values
    Field<Type>& patchValues = *this;

    // Get bracketing times t1 and t2 (currentTime is between them)
    const scalar t1 = sampleTimes_[sampleTimeIndex1_].value();
    const scalar t2 = sampleTimes_[sampleTimeIndex2_].value();

    if (sampleTimeIndex1_ == sampleTimeIndex2_)
    {
        // No interpolation needed - exact time match or at boundary
        patchValues = sampledValues1_ + offset_;
    }
    else
    {
        // Time does not exactly match a sample time
        // Check if interpolation is allowed
        if (!allowTimeInterpolation_)
        {
            // Use tolerance based on deltaT (1% of timestep interval)
            const scalar timeTol = 0.01 * (t2 - t1);

            // Check if we're close enough to t1 or t2 to consider it a match
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
                    << "    Time interpolation is disabled (allowTimeInterpolation = false)." << nl
                    << "    This BC requires exact timestep matching with the extraction." << nl << nl
                    << "    Options:" << nl
                    << "      1. Ensure reconstruction uses identical deltaT as extraction" << nl
                    << "      2. Set allowTimeInterpolation = true in the BC configuration" << nl
                    << "         to permit linear interpolation between available timesteps" << nl
                    << exit(FatalError);
            }

            // Close enough to t1 or t2 - use the nearest without interpolation
            if (closeToT2)
            {
                patchValues = sampledValues2_ + offset_;
            }
            else
            {
                patchValues = sampledValues1_ + offset_;
            }
        }
        else
        {
            // Temporal interpolation is allowed
            const scalar alpha = (currentTime - t1) / (t2 - t1 + SMALL);

            // Check if cubic interpolation is possible and requested
            // Cubic requires 4 distinct points; if i0==i1 or i2==i3, fall back to linear
            const bool useCubic =
                (timeInterpolationScheme_ == timeInterpScheme::CUBIC)
             && (sampleTimeIndex0_ != sampleTimeIndex1_)
             && (sampleTimeIndex2_ != sampleTimeIndex3_);

            if (useCubic)
            {
                // Centripetal Catmull-Rom cubic spline interpolation (per-face)
                // Get all 4 time values for proper non-uniform spacing handling
                const scalar t0 = sampleTimes_[sampleTimeIndex0_].value();
                const scalar t3 = sampleTimes_[sampleTimeIndex3_].value();

                patchValues.setSize(sampledValues1_.size());

                forAll(patchValues, facei)
                {
                    patchValues[facei] = catmullRomCentripetal
                    (
                        sampledValues0_[facei],
                        sampledValues1_[facei],
                        sampledValues2_[facei],
                        sampledValues3_[facei],
                        t0, t1, t2, t3,
                        currentTime
                    ) + offset_;
                }

                DebugInfo
                    << "Cubic time interpolation (centripetal): t=" << currentTime
                    << " using t0=" << t0 << " t1=" << t1
                    << " t2=" << t2 << " t3=" << t3 << endl;
            }
            else
            {
                // Linear temporal interpolation
                patchValues =
                    (1.0 - alpha) * sampledValues1_
                  + alpha * sampledValues2_
                  + offset_;

                DebugInfo
                    << "Linear time interpolation: t=" << currentTime
                    << " between " << t1 << " and " << t2
                    << " alpha=" << alpha << endl;
            }
        }
    }

    // Report flux through patch (vector fields only, i.e., velocity)
    // This reports BEFORE adjustPhi, showing the BC's prescribed flux contribution
    // and the total mesh flux imbalance that adjustPhi will need to correct
    if (reportFlux_)
    {
        if constexpr (std::is_same_v<Type, vector>)
        {
            // Compute flux through this patch: phi = U & Sf (face area vector)
            const vectorField& Sf = this->patch().Sf();
            const scalarField flux(patchValues & Sf);
            const scalar thisNetFlux = gSum(flux);
            const scalar thisInFlux = gSum(neg(flux) * flux);   // negative = into domain
            const scalar thisOutFlux = gSum(pos(flux) * flux);  // positive = out of domain

            // Compute total mesh flux from ALL boundary patches
            // This shows the imbalance that adjustPhi will correct
            const fvMesh& mesh = this->patch().boundaryMesh().mesh();
            const volVectorField& U =
                mesh.template lookupObject<volVectorField>(this->internalField().name());

            scalar totalMeshFlux = 0;
            scalar fixedFlux = 0;      // Flux from patches where fixesValue=true
            scalar adjustableFlux = 0; // Flux from patches where fixesValue=false

            forAll(U.boundaryField(), patchi)
            {
                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                const vectorField& pSf = Up.patch().Sf();
                const scalar patchFlux = gSum(Up & pSf);
                totalMeshFlux += patchFlux;

                if (Up.fixesValue())
                {
                    fixedFlux += patchFlux;
                }
                else
                {
                    adjustableFlux += patchFlux;
                }
            }

            Info<< "spaceTimeWindow flux [" << this->patch().name() << "]"
                << " t=" << this->db().time().value()
                << " thisPatch=" << thisNetFlux
                << " (in=" << thisInFlux << " out=" << thisOutFlux << ")"
                << " | MESH TOTAL=" << totalMeshFlux
                << " (fixed=" << fixedFlux
                << " adjustable=" << adjustableFlux << ")"
                << " | adjustPhi will correct: " << -totalMeshFlux
                << endl;
        }
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::spaceTimeWindowFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);

    os.writeEntry("dataDir", dataDir_);

    if (fieldTableName_ != this->internalField().name())
    {
        os.writeEntry("fieldTableName", fieldTableName_);
    }

    if (setAverage_)
    {
        os.writeEntry("setAverage", setAverage_);
    }

    if (mag(offset_) > SMALL)
    {
        os.writeEntry("offset", offset_);
    }

    if (!fixesValue_)
    {
        os.writeEntry("fixesValue", fixesValue_);
    }

    if (allowTimeInterpolation_)
    {
        os.writeEntry("allowTimeInterpolation", allowTimeInterpolation_);

        // Only write scheme if interpolation is enabled (otherwise irrelevant)
        if (timeInterpolationScheme_ == timeInterpScheme::CUBIC)
        {
            os.writeEntry("timeInterpolationScheme", "cubic");
        }
        // Don't write "linear" - it's the default
    }

    if (reportFlux_)
    {
        os.writeEntry("reportFlux", reportFlux_);
    }
}


// ************************************************************************* //
