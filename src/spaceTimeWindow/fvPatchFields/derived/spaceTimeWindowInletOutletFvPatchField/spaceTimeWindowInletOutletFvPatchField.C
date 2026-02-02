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

    // No header - put token back and read raw field data
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
    metadataValidated_(false)
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
    metadataValidated_(false)
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
    metadataValidated_(ptf.metadataValidated_)
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
    metadataValidated_(ptf.metadataValidated_)
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
    metadataValidated_(ptf.metadataValidated_)
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
