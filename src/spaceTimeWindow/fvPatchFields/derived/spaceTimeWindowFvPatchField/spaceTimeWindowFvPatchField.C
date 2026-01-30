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
        startSampleTime_ = -1;
        endSampleTime_ = -1;
    }

    const scalar currentTime = this->db().time().value();

    // Find bracketing times
    label lo = -1;
    label hi = -1;

    for (label i = 0; i < sampleTimes_.size(); ++i)
    {
        if (sampleTimes_[i].value() <= currentTime)
        {
            lo = i;
        }
        if (sampleTimes_[i].value() >= currentTime && hi < 0)
        {
            hi = i;
        }
    }

    // Strict time bounds - NO extrapolation allowed
    if (lo < 0)
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
    if (hi < 0)
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

    // Check if we need to load new data
    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;
        startSampledValues_ = readFieldData(sampleTimes_[lo]);

        DebugInfo
            << "Loaded start time data: " << sampleTimes_[lo].name() << endl;
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;
        endSampledValues_ = readFieldData(sampleTimes_[hi]);

        DebugInfo
            << "Loaded end time data: " << sampleTimes_[hi].name() << endl;
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
    sampleTimes_(),
    startSampleTime_(-1),
    endSampleTime_(-1),
    startSampledValues_(),
    endSampledValues_(),
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
    sampleTimes_(),
    startSampleTime_(-1),
    endSampleTime_(-1),
    startSampledValues_(),
    endSampledValues_(),
    metadataValidated_(false)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
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
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(-1),
    endSampleTime_(-1),
    startSampledValues_(),
    endSampledValues_(),
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
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    endSampleTime_(ptf.endSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampledValues_(ptf.endSampledValues_),
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
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    endSampleTime_(ptf.endSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    endSampledValues_(ptf.endSampledValues_),
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
    startSampledValues_.clear();
    endSampledValues_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
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
    startSampledValues_.clear();
    endSampledValues_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
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

    if (startSampleTime_ == endSampleTime_)
    {
        // No interpolation needed - exact time match or at boundary
        patchValues = startSampledValues_ + offset_;
    }
    else
    {
        // Linear temporal interpolation
        const scalar t0 = sampleTimes_[startSampleTime_].value();
        const scalar t1 = sampleTimes_[endSampleTime_].value();
        const scalar alpha = (currentTime - t0) / (t1 - t0 + SMALL);

        patchValues =
            (1.0 - alpha) * startSampledValues_
          + alpha * endSampledValues_
          + offset_;

        DebugInfo
            << "Time interpolation: t=" << currentTime
            << " between " << t0 << " and " << t1
            << " alpha=" << alpha << endl;
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
}


// ************************************************************************* //
