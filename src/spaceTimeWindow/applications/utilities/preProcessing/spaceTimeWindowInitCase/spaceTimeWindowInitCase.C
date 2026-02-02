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

Application
    spaceTimeWindowInitCase

Description
    Initialize a reconstruction case from spaceTimeWindowExtract output.

    Reads the extraction metadata and source case configuration to create
    a fully configured reconstruction case that can be run directly with
    the solver.

    Creates:
    - system/controlDict with matching time settings and solver
    - system/fvSchemes, fvSolution copied from source
    - constant/ files (turbulenceProperties, transportProperties, etc.)
    - Initial fields with spaceTimeWindow BC on oldInternalFaces

Usage
    spaceTimeWindowInitCase [OPTIONS]

    Options:
        -sourceCase <dir>   Source case directory (where extraction ran)
        -extractDir <dir>   Directory containing extracted data (default: current dir)
        -overwrite          Overwrite existing files
        -outletDirection <vector>  Normal direction for outlet patch (e.g., "(1 0 0)")
        -outletFraction <scalar>   Fraction of box face to use as outlet (default: 0.1)

    The outlet patch provides a pressure relief for incompressible solvers.
    Without an outlet, the all-Dirichlet velocity BC can cause pressure buildup.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OSspecific.H"
#include "foamVersion.H"
#include "Tuple2.H"
#include "DynamicList.H"
#include "SortableList.H"
#include "boundBox.H"
#include "deltaVarintCodec.H"
#include "deltaVarintTemporalCodec.H"
#include "sodiumCrypto.H"
#include <fstream>
#include <sstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void copyFile(const fileName& src, const fileName& dst)
{
    if (isFile(src))
    {
        cp(src, dst);
        Info<< "    Copied: " << src.name() << endl;
    }
    else
    {
        WarningInFunction
            << "Source file not found: " << src << endl;
    }
}

void copyDirectory(const fileName& src, const fileName& dst)
{
    if (isDir(src))
    {
        cp(src, dst);
        Info<< "    Copied directory: " << src.name() << endl;
    }
}


// Decrypt and decompress boundary data files
// Returns the base field name (without .enc or .dvz.enc extensions)
word decryptBoundaryFile
(
    const fileName& encryptedPath,
    const fileName& outputDir,
    const std::vector<uint8_t>& publicKey,
    const std::vector<uint8_t>& privateKey,
    bool verbose = true
)
{
    word fieldName = encryptedPath.name();

    // Determine file type from extension
    // Possible patterns: fieldName.enc, fieldName.dvz.enc
    const word encExt = "." + sodiumCrypto::fileExtension;
    const word dvzExt = "." + deltaVarintCodec::fileExtension();
    const word dvzEncExt = dvzExt + encExt;

    bool isDvzEnc = fieldName.ends_with(dvzEncExt);
    bool isEnc = fieldName.ends_with(encExt);

    if (!isEnc)
    {
        // Not an encrypted file
        return word::null;
    }

    // Strip extensions to get base field name
    word baseFieldName;
    if (isDvzEnc)
    {
        baseFieldName = fieldName.substr(0, fieldName.size() - dvzEncExt.size());
    }
    else
    {
        baseFieldName = fieldName.substr(0, fieldName.size() - encExt.size());
    }

#ifdef FOAM_USE_SODIUM
    // Decrypt the file
    std::vector<uint8_t> decrypted = sodiumCrypto::decryptFromFile
    (
        encryptedPath,
        publicKey,
        privateKey
    );

    if (decrypted.empty())
    {
        FatalErrorInFunction
            << "Failed to decrypt file: " << encryptedPath << nl
            << exit(FatalError);
    }

    fileName outputPath = outputDir / baseFieldName;

    if (isDvzEnc)
    {
        // Decrypted data is dvz-compressed - decompress and write as OpenFOAM format
        // Detect if it's scalar or vector from the dvz header
        if (deltaVarintCodec::isDeltaVarintBuffer(decrypted))
        {
            // Read header to determine type
            // DVZ format: magic(4) + nFaces(4) + nComponents(4) + precision(4) + data
            uint32_t nComponents = 0;
            if (decrypted.size() >= 12)
            {
                std::memcpy(&nComponents, decrypted.data() + 8, sizeof(uint32_t));
            }

            if (nComponents == 3)
            {
                // Vector field
                vectorField field = deltaVarintCodec::decodeVector(decrypted);

                OFstream os(outputPath);
                os  << "FoamFile" << nl
                    << "{" << nl
                    << "    version     2.0;" << nl
                    << "    format      ascii;" << nl
                    << "    class       vectorField;" << nl
                    << "    object      " << baseFieldName << ";" << nl
                    << "}" << nl << nl;
                os << field;

                if (verbose)
                {
                    Info<< "    Decrypted+decompressed: " << fieldName
                        << " -> " << baseFieldName << " (vector, "
                        << field.size() << " values)" << endl;
                }
            }
            else
            {
                // Scalar field
                scalarField field = deltaVarintCodec::decodeScalar(decrypted);

                OFstream os(outputPath);
                os  << "FoamFile" << nl
                    << "{" << nl
                    << "    version     2.0;" << nl
                    << "    format      ascii;" << nl
                    << "    class       scalarField;" << nl
                    << "    object      " << baseFieldName << ";" << nl
                    << "}" << nl << nl;
                os << field;

                if (verbose)
                {
                    Info<< "    Decrypted+decompressed: " << fieldName
                        << " -> " << baseFieldName << " (scalar, "
                        << field.size() << " values)" << endl;
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Decrypted data is not valid DVZ format: " << encryptedPath << nl
                << exit(FatalError);
        }
    }
    else
    {
        // Decrypted data is raw binary - write as-is
        // Format: label(size) + raw Type data
        // Determine type from size: scalar=8 bytes per value, vector=24 bytes per value
        label fieldSize = 0;
        std::memcpy(&fieldSize, decrypted.data(), sizeof(label));

        size_t dataSize = decrypted.size() - sizeof(label);
        size_t expectedScalarSize = fieldSize * sizeof(scalar);
        size_t expectedVectorSize = fieldSize * sizeof(vector);

        OFstream os(outputPath);
        os  << "FoamFile" << nl
            << "{" << nl
            << "    version     2.0;" << nl
            << "    format      ascii;" << nl;

        if (dataSize == expectedVectorSize)
        {
            // Vector field
            vectorField field(fieldSize);
            std::memcpy(field.data(), decrypted.data() + sizeof(label), dataSize);

            os  << "    class       vectorField;" << nl
                << "    object      " << baseFieldName << ";" << nl
                << "}" << nl << nl;
            os << field;

            if (verbose)
            {
                Info<< "    Decrypted: " << fieldName
                    << " -> " << baseFieldName << " (vector, "
                    << fieldSize << " values)" << endl;
            }
        }
        else if (dataSize == expectedScalarSize)
        {
            // Scalar field
            scalarField field(fieldSize);
            std::memcpy(field.data(), decrypted.data() + sizeof(label), dataSize);

            os  << "    class       scalarField;" << nl
                << "    object      " << baseFieldName << ";" << nl
                << "}" << nl << nl;
            os << field;

            if (verbose)
            {
                Info<< "    Decrypted: " << fieldName
                    << " -> " << baseFieldName << " (scalar, "
                    << fieldSize << " values)" << endl;
            }
        }
        else
        {
            FatalErrorInFunction
                << "Cannot determine field type from decrypted data size" << nl
                << "    File: " << encryptedPath << nl
                << "    Size: " << fieldSize << " values, " << dataSize << " bytes" << nl
                << exit(FatalError);
        }
    }

    // Remove the encrypted file after successful decryption
    rm(encryptedPath);

    return baseFieldName;
#else
    FatalErrorInFunction
        << "Cannot decrypt file - library built without libsodium support" << nl
        << "    File: " << encryptedPath << nl
        << "    Rebuild with FOAM_USE_SODIUM=1 to enable decryption" << nl
        << exit(FatalError);

    return word::null;
#endif
}


// Convert all DVZT (delta-varint-temporal) files to DVZ (delta-varint) in boundaryData
// DVZT requires sequential processing because delta frames need previous timestep data
// Output is DVZ format which the spaceTimeWindow BC can handle at runtime
// Returns the number of files converted
label convertDvztToDvzBoundaryData
(
    const fileName& patchBoundaryDir,
    const DynamicList<Tuple2<scalar, word>>& timeList,  // Must be sorted by time!
    label dvzPrecision = 6,
    bool verbose = true
)
{
    const word dvztExt = "." + deltaVarintTemporalCodec::fileExtension();
    const word dvzExt = "." + deltaVarintCodec::fileExtension();
    label totalConverted = 0;

    // Track previous timestep fields for each field name
    // Needed for delta frame decompression
    HashTable<scalarField> prevScalarFields;
    HashTable<vectorField> prevVectorFields;

    Info<< nl << "Converting DVZT to DVZ boundary data files..." << endl;
    Info<< "    (Sequential processing required for temporal delta frames)" << endl;
    Info<< "    Output: DVZ format (spaceTimeWindow BC handles at runtime)" << endl;

    // Process timesteps in chronological order (essential for DVZT!)
    for (const auto& timePair : timeList)
    {
        const word& timeName = timePair.second();
        fileName timeDir = patchBoundaryDir / timeName;

        // Get all files in this time directory
        fileNameList files = readDir(timeDir, fileName::FILE);

        for (const fileName& f : files)
        {
            // Check if it's a DVZT file
            if (!f.ends_with(dvztExt))
            {
                continue;
            }

            // Get base field name (strip .dvzt extension)
            word baseFieldName = f.substr(0, f.size() - dvztExt.size());
            fileName dvztPath = timeDir / f;
            fileName dvzOutputPath = timeDir / (baseFieldName + dvzExt);

            // Check if file is a keyframe or delta frame
            bool isKeyframe = deltaVarintTemporalCodec::isKeyframe(dvztPath);

            // Read the DVZT file to determine type (scalar or vector)
            // Read header to get nComponents
            std::ifstream ifs(dvztPath.c_str(), std::ios::binary);
            if (!ifs)
            {
                WarningInFunction
                    << "Cannot open DVZT file: " << dvztPath << endl;
                continue;
            }

            // Read header: magic(4) + nFaces(4) + nComponents(4) + precision(4) + flags(4)
            uint32_t magic, nFaces, nComponents, precision, flags;
            ifs.read(reinterpret_cast<char*>(&magic), sizeof(magic));
            ifs.read(reinterpret_cast<char*>(&nFaces), sizeof(nFaces));
            ifs.read(reinterpret_cast<char*>(&nComponents), sizeof(nComponents));
            ifs.read(reinterpret_cast<char*>(&precision), sizeof(precision));
            ifs.close();

            // Verify magic number
            if (magic != deltaVarintTemporalCodec::magic())
            {
                WarningInFunction
                    << "Invalid DVZT magic number in: " << dvztPath << endl;
                continue;
            }

            // Use precision from DVZT file for DVZ output
            uint32_t outputPrecision = (precision > 0) ? precision : dvzPrecision;

            if (nComponents == 3)
            {
                // Vector field
                const vectorField* prevPtr = nullptr;
                if (!isKeyframe && prevVectorFields.found(baseFieldName))
                {
                    prevPtr = &prevVectorFields[baseFieldName];
                }
                else if (!isKeyframe)
                {
                    FatalErrorInFunction
                        << "DVZT delta frame requires previous timestep data" << nl
                        << "    File: " << dvztPath << nl
                        << "    Field: " << baseFieldName << nl
                        << "    This usually means timesteps are out of order" << nl
                        << exit(FatalError);
                }

                // Decompress DVZT
                vectorField field = deltaVarintTemporalCodec::readVector(dvztPath, prevPtr);

                // Store for next timestep
                prevVectorFields.set(baseFieldName, field);

                // Write as DVZ format
                deltaVarintCodec::write(dvzOutputPath, field, outputPrecision);

                if (verbose)
                {
                    Info<< "    Converted: " << f << " -> " << baseFieldName << dvzExt
                        << " (vector, " << field.size() << " values, "
                        << (isKeyframe ? "keyframe" : "delta") << ")" << endl;
                }
            }
            else
            {
                // Scalar field
                const scalarField* prevPtr = nullptr;
                if (!isKeyframe && prevScalarFields.found(baseFieldName))
                {
                    prevPtr = &prevScalarFields[baseFieldName];
                }
                else if (!isKeyframe)
                {
                    FatalErrorInFunction
                        << "DVZT delta frame requires previous timestep data" << nl
                        << "    File: " << dvztPath << nl
                        << "    Field: " << baseFieldName << nl
                        << "    This usually means timesteps are out of order" << nl
                        << exit(FatalError);
                }

                // Decompress DVZT
                scalarField field = deltaVarintTemporalCodec::readScalar(dvztPath, prevPtr);

                // Store for next timestep
                prevScalarFields.set(baseFieldName, field);

                // Write as DVZ format
                deltaVarintCodec::write(dvzOutputPath, field, outputPrecision);

                if (verbose)
                {
                    Info<< "    Converted: " << f << " -> " << baseFieldName << dvzExt
                        << " (scalar, " << field.size() << " values, "
                        << (isKeyframe ? "keyframe" : "delta") << ")" << endl;
                }
            }

            // Remove the DVZT file after successful conversion
            rm(dvztPath);
            totalConverted++;
        }
    }

    Info<< "    Converted " << totalConverted << " DVZT files to DVZ across "
        << timeList.size() << " timesteps" << nl << endl;

    return totalConverted;
}


void writeControlDict
(
    const fileName& targetDir,
    const dictionary& metadata,
    const dictionary& sourceControlDict,
    const word& startTimeName,
    const word& endTimeName
)
{
    fileName controlDictPath = targetDir / "system" / "controlDict";
    mkDir(targetDir / "system");

    OFstream os(controlDictPath);

    // Write header
    os  << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       dictionary;" << nl
        << "    object      controlDict;" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << "// Generated by spaceTimeWindowInitCase" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << nl;

    // Application - from metadata
    word solver = metadata.getOrDefault<word>("solver", "pimpleFoam");
    os.writeEntry("application", solver);
    os << nl;

    // Time settings from metadata
    scalar deltaT = metadata.get<scalar>("deltaT");
    bool adjustTimeStep = metadata.get<bool>("adjustTimeStep");
    label timePrecision = metadata.getOrDefault<label>("timePrecision", 8);

    // Write time values as strings to match directory names exactly
    os.writeEntry("startFrom", word("startTime"));
    os  << "startTime       " << startTimeName << ";" << nl;
    os << nl;

    // End time from boundaryData time range
    os.writeEntry("stopAt", word("endTime"));
    os  << "endTime         " << endTimeName << ";" << nl;
    os << nl;

    os.writeEntry("deltaT", deltaT);
    os << nl;

    // Write control
    os.writeEntry("writeControl", word("timeStep"));
    os.writeEntry("writeInterval", 100);
    os << nl;

    os.writeEntry("purgeWrite", 0);
    os << nl;

    os.writeEntry("writeFormat", word("ascii"));
    os.writeEntry("writePrecision", 8);
    os.writeEntry("writeCompression", word("off"));
    os << nl;

    os.writeEntry("timeFormat", word("fixed"));
    os.writeEntry("timePrecision", timePrecision);
    os << nl;

    os.writeEntry("runTimeModifiable", true);
    os << nl;

    os.writeEntry("adjustTimeStep", adjustTimeStep ? "yes" : "no");
    os << nl;

    // Copy maxCo and related settings if they exist
    if (sourceControlDict.found("maxCo"))
    {
        os.writeEntry("maxCo", sourceControlDict.get<scalar>("maxCo"));
    }
    if (sourceControlDict.found("maxDeltaT"))
    {
        os.writeEntry("maxDeltaT", sourceControlDict.get<scalar>("maxDeltaT"));
    }

    // Add libs entry for spaceTimeWindow
    os  << nl
        << "libs" << nl
        << "(" << nl
        << "    spaceTimeWindow" << nl
        << ");" << nl;

    Info<< "    Created: system/controlDict" << nl
        << "        solver = " << solver << nl
        << "        startTime = " << startTimeName << nl
        << "        endTime = " << endTimeName << nl
        << "        deltaT = " << deltaT << nl
        << "        timePrecision = " << timePrecision << nl
        << "        adjustTimeStep = " << (adjustTimeStep ? "yes" : "no") << endl;
}


void writeInitialField
(
    const fileName& targetDir,
    const word& fieldName,
    const word& fieldType,
    const dictionary& fieldDict,
    const word& timeName,
    const polyBoundaryMesh& bMesh
)
{
    fileName fieldPath = targetDir / timeName / fieldName;
    mkDir(targetDir / timeName);

    OFstream os(fieldPath);

    // Determine class name
    word className;
    if (fieldType == "scalar")
    {
        className = "volScalarField";
    }
    else if (fieldType == "vector")
    {
        className = "volVectorField";
    }
    else
    {
        className = "volScalarField";  // Default
    }

    // Write header
    os  << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       " << className << ";" << nl
        << "    object      " << fieldName << ";" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << nl;

    // Get dimensions from source field if available
    if (fieldDict.found("dimensions"))
    {
        os  << "dimensions      " << fieldDict.lookup("dimensions") << ";" << nl << nl;
    }
    else
    {
        // Default dimensions
        if (fieldName == "U")
        {
            os << "dimensions      [0 1 -1 0 0 0 0];" << nl << nl;
        }
        else if (fieldName == "p")
        {
            os << "dimensions      [0 2 -2 0 0 0 0];" << nl << nl;
        }
        else
        {
            os << "dimensions      [0 0 0 0 0 0 0];" << nl << nl;
        }
    }

    // Internal field - use source if available
    if (fieldDict.found("internalField"))
    {
        os  << "internalField   " << fieldDict.lookup("internalField") << ";" << nl << nl;
    }
    else
    {
        if (fieldType == "vector")
        {
            os << "internalField   uniform (0 0 0);" << nl << nl;
        }
        else
        {
            os << "internalField   uniform 0;" << nl << nl;
        }
    }

    // Boundary field
    os << "boundaryField" << nl << "{" << nl;

    // Get source boundary field if available
    const dictionary* bfDictPtr = fieldDict.findDict("boundaryField");

    forAll(bMesh, patchi)
    {
        const word& patchName = bMesh[patchi].name();

        os << "    " << patchName << nl
           << "    {" << nl;

        if (patchName == "oldInternalFaces")
        {
            // Use spaceTimeWindow BC
            os  << "        type            spaceTimeWindow;" << nl
                << "        dataDir         \"constant/boundaryData\";" << nl
                << "        fixesValue      true;" << nl;

            if (fieldType == "vector")
            {
                os << "        value           uniform (0 0 0);" << nl;
            }
            else
            {
                os << "        value           uniform 0;" << nl;
            }
        }
        else if (bfDictPtr && bfDictPtr->found(patchName))
        {
            // Copy from source
            const dictionary& patchDict = bfDictPtr->subDict(patchName);
            for (const entry& e : patchDict)
            {
                if (e.isStream())
                {
                    os << "        " << e.keyword() << "    " << e.stream() << ";" << nl;
                }
                else if (e.isDict())
                {
                    os << "        " << e.keyword() << nl;
                    e.dict().write(os, false);
                }
            }
        }
        else
        {
            // Default BC based on patch type
            const word& patchType = bMesh[patchi].type();
            if (patchType == "wall")
            {
                if (fieldName == "U")
                {
                    os  << "        type            noSlip;" << nl;
                }
                else if (fieldName == "p")
                {
                    os  << "        type            zeroGradient;" << nl;
                }
                else if (fieldName == "nut" || fieldName == "nuTilda")
                {
                    os  << "        type            nutWallFunction;" << nl
                        << "        value           uniform 0;" << nl;
                }
                else if (fieldName == "k")
                {
                    os  << "        type            kqRWallFunction;" << nl
                        << "        value           uniform 0;" << nl;
                }
                else if (fieldName == "epsilon")
                {
                    os  << "        type            epsilonWallFunction;" << nl
                        << "        value           uniform 0;" << nl;
                }
                else if (fieldName == "omega")
                {
                    os  << "        type            omegaWallFunction;" << nl
                        << "        value           uniform 0;" << nl;
                }
                else
                {
                    os  << "        type            zeroGradient;" << nl;
                }
            }
            else if (patchType == "empty")
            {
                os  << "        type            empty;" << nl;
            }
            else if (patchType == "symmetry" || patchType == "symmetryPlane")
            {
                os  << "        type            symmetry;" << nl;
            }
            else
            {
                os  << "        type            zeroGradient;" << nl;
            }
        }

        os << "    }" << nl;
    }

    os << "}" << nl;
    os << nl << "// ************************************************************************* //" << nl;

    Info<< "    Created: " << timeName << "/" << fieldName << endl;
}


// Correct velocity field for exact mass conservation using least-squares approach
// The correction is: U_corrected = U - (imbalance / totalSfMag) * n
// where n is the face normal (Sf / |Sf|)
// This minimizes sum(|U_corrected - U|^2) subject to sum(U_corrected . Sf) = 0
void correctMassFluxInBoundaryData
(
    const fileName& patchBoundaryDir,
    const DynamicList<Tuple2<scalar, word>>& timeList,
    const vectorField& faceSf  // Face area vectors (Sf)
)
{
    Info<< nl << "Applying mass flux correction to boundaryData..." << endl;

    // Compute total face area magnitude for normalization
    scalar totalSfMag = 0;
    forAll(faceSf, facei)
    {
        totalSfMag += mag(faceSf[facei]);
    }

    if (totalSfMag < SMALL)
    {
        WarningInFunction
            << "Total face area is zero, skipping mass flux correction" << endl;
        return;
    }

    // Statistics tracking
    scalar maxImbalanceBefore = 0;
    scalar maxImbalanceAfter = 0;
    scalar sumImbalanceBefore = 0;
    label nTimesteps = 0;

    // Process each time directory
    for (const auto& timePair : timeList)
    {
        const word& timeName = timePair.second();
        fileName timeDir = patchBoundaryDir / timeName;

        // Look for velocity field (U or U.dvz)
        fileName velocityPath;
        word velocityFileName;

        if (isFile(timeDir / "U"))
        {
            velocityPath = timeDir / "U";
            velocityFileName = "U";
        }
        else
        {
            // No velocity field in this timestep (unusual but possible)
            continue;
        }

        // Read the velocity field
        std::ifstream ifs(velocityPath.c_str());
        std::string content((std::istreambuf_iterator<char>(ifs)),
                           std::istreambuf_iterator<char>());
        ifs.close();

        // Find the data section after the FoamFile header
        std::size_t headerEnd = content.find('}');
        if (headerEnd == std::string::npos)
        {
            WarningInFunction
                << "Cannot find FoamFile header end in " << velocityPath << nl
                << "    Skipping..." << endl;
            continue;
        }

        // Parse the velocity field
        std::string dataSection = content.substr(headerEnd + 1);
        IStringStream fieldIs(dataSection);
        vectorField U(fieldIs);

        if (U.size() != faceSf.size())
        {
            WarningInFunction
                << "Velocity field size " << U.size()
                << " doesn't match face count " << faceSf.size() << nl
                << "    Skipping timestep " << timeName << endl;
            continue;
        }

        // Compute current mass flux imbalance: sum(U . Sf)
        scalar imbalance = 0;
        forAll(U, facei)
        {
            imbalance += (U[facei] & faceSf[facei]);
        }

        // Track statistics
        maxImbalanceBefore = max(maxImbalanceBefore, mag(imbalance));
        sumImbalanceBefore += mag(imbalance);
        nTimesteps++;

        // Apply correction: U_corrected = U - (imbalance / totalSfMag) * n
        // where n = Sf / |Sf| is the unit face normal
        // This distributes the correction uniformly by face area
        vectorField U_corrected(U.size());
        forAll(U, facei)
        {
            scalar sfMag = mag(faceSf[facei]);
            if (sfMag > SMALL)
            {
                vector n = faceSf[facei] / sfMag;
                U_corrected[facei] = U[facei] - (imbalance / totalSfMag) * n;
            }
            else
            {
                U_corrected[facei] = U[facei];
            }
        }

        // Verify correction (should be ~0)
        scalar imbalanceAfter = 0;
        forAll(U_corrected, facei)
        {
            imbalanceAfter += (U_corrected[facei] & faceSf[facei]);
        }
        maxImbalanceAfter = max(maxImbalanceAfter, mag(imbalanceAfter));

        // Write corrected velocity field
        OFstream os(velocityPath);
        os  << "FoamFile" << nl
            << "{" << nl
            << "    version     2.0;" << nl
            << "    format      ascii;" << nl
            << "    class       vectorField;" << nl
            << "    object      U;" << nl
            << "}" << nl << nl;
        os << U_corrected;
    }

    // Report statistics
    if (nTimesteps > 0)
    {
        Info<< "    Processed " << nTimesteps << " timesteps" << nl
            << "    Max imbalance before: " << maxImbalanceBefore << nl
            << "    Max imbalance after:  " << maxImbalanceAfter << nl
            << "    Avg imbalance before: " << sumImbalanceBefore / nTimesteps << nl
            << "    Correction formula: U_corrected = U - (imbalance/totalSfMag) * n" << endl;
    }
    else
    {
        WarningInFunction
            << "No velocity fields found in boundaryData" << endl;
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Initialize a reconstruction case from spaceTimeWindowExtract output"
    );

    argList::addOption
    (
        "sourceCase",
        "dir",
        "Source case directory (where extraction ran)"
    );

    argList::addOption
    (
        "extractDir",
        "dir",
        "Directory containing extracted data (default: current dir)"
    );

    argList::addBoolOption
    (
        "overwrite",
        "Overwrite existing files"
    );

    argList::addOption
    (
        "outletDirection",
        "vector",
        "Normal direction for outlet patch (e.g., \"(1 0 0)\" for +x)"
    );

    argList::addOption
    (
        "outletFraction",
        "scalar",
        "Fraction of box extent for outlet region (default: 0.1)"
    );

    argList::addBoolOption
    (
        "correctMassFlux",
        "Apply least-squares mass flux correction to boundaryData velocity"
    );

    argList::addBoolOption
    (
        "inletOutletBC",
        "Use spaceTimeWindowInletOutlet BC for U (flux-based Dirichlet/zeroGradient) "
        "and zeroGradient for all scalar fields (p, turbulence). Recommended for vortex shedding flows."
    );

    argList::addOption
    (
        "initialFields",
        "wordList",
        "Fields to extract for initial conditions (e.g., \"(U p nut)\"). "
        "Default: same as fields in boundaryData"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"

    // Get arguments
    fileName sourceCase;
    if (!args.readIfPresent("sourceCase", sourceCase))
    {
        FatalErrorInFunction
            << "Must specify -sourceCase <dir>" << nl
            << exit(FatalError);
    }

    fileName extractDir = args.getOrDefault<fileName>("extractDir", args.path());
    bool overwrite = args.found("overwrite");

    // Outlet configuration
    vector outletDirection = Zero;
    scalar outletFraction = 0.1;
    bool createOutlet = false;

    if (args.found("outletDirection"))
    {
        outletDirection = args.get<vector>("outletDirection");
        outletDirection /= mag(outletDirection) + SMALL;  // Normalize
        createOutlet = true;
        outletFraction = args.getOrDefault<scalar>("outletFraction", 0.1);
    }

    // Mass flux correction option
    bool correctMassFlux = args.found("correctMassFlux");

    // Inlet-outlet BC option: U uses spaceTimeWindowInletOutlet, p/nut use zeroGradient
    bool useInletOutletBC = args.found("inletOutletBC");

    // Validate incompatible options: -inletOutletBC is incompatible with -outletDirection/-outletFraction
    // The inlet-outlet BC approach uses flux-based switching on oldInternalFaces
    // (Dirichlet at inflow, zeroGradient at outflow), so a separate outlet patch is unnecessary
    // and would create conflicting boundary conditions
    if (useInletOutletBC && createOutlet)
    {
        FatalErrorInFunction
            << "Options -inletOutletBC and -outletDirection are incompatible." << nl
            << nl
            << "  -inletOutletBC: Uses spaceTimeWindowInletOutlet BC which automatically" << nl
            << "                  applies Dirichlet at inflow faces and zeroGradient at" << nl
            << "                  outflow faces based on the flux direction." << nl
            << nl
            << "  -outletDirection: Creates a separate outlet patch for pressure relief" << nl
            << "                    and applies Dirichlet on all remaining faces." << nl
            << nl
            << "Choose one approach:" << nl
            << "  1. For vortex shedding/unsteady flows: use -inletOutletBC (recommended)" << nl
            << "     This allows natural outflow without a fixed outlet region." << nl
            << nl
            << "  2. For flows with known outlet direction: use -outletDirection" << nl
            << "     This creates a fixed outlet patch for pressure relief." << nl
            << exit(FatalError);
    }

    // User-specified initial fields (if not specified, default to fields in boundaryData)
    wordList userInitialFields;
    bool hasUserInitialFields = false;
    if (args.found("initialFields"))
    {
        IStringStream is(args["initialFields"]);
        is >> userInitialFields;
        hasUserInitialFields = true;
    }

    Info<< "spaceTimeWindowInitCase" << nl
        << "    Source case:    " << sourceCase << nl
        << "    Extract dir:    " << extractDir << nl;

    if (createOutlet)
    {
        Info<< "    Outlet direction: " << outletDirection << nl
            << "    Outlet fraction:  " << outletFraction << nl;
    }
    if (correctMassFlux)
    {
        Info<< "    Mass flux correction: enabled" << nl;
    }
    if (useInletOutletBC)
    {
        Info<< "    Inlet-outlet BC mode: enabled" << nl
            << "        U: spaceTimeWindowInletOutlet (flux-based Dirichlet/zeroGradient)" << nl
            << "        All scalar fields (p, nut, k, epsilon, omega, etc.): zeroGradient" << nl;
    }
    if (hasUserInitialFields)
    {
        Info<< "    Initial fields: " << userInitialFields << nl;
    }
    Info<< endl;

    // Verify source case exists
    if (!isDir(sourceCase))
    {
        FatalErrorInFunction
            << "Source case directory not found: " << sourceCase << nl
            << exit(FatalError);
    }

    // Find and read extraction metadata
    // Look for oldInternalFaces in boundaryData
    fileName boundaryDataDir = extractDir / "constant" / "boundaryData";
    fileNameList patches = readDir(boundaryDataDir, fileName::DIRECTORY);

    fileName metadataPath;
    word oldInternalPatchName;

    for (const fileName& patchDir : patches)
    {
        fileName testPath = boundaryDataDir / patchDir / "extractionMetadata";
        if (isFile(testPath))
        {
            metadataPath = testPath;
            oldInternalPatchName = patchDir;
            break;
        }
    }

    if (metadataPath.empty())
    {
        FatalErrorInFunction
            << "No extractionMetadata found in " << boundaryDataDir << nl
            << "    Run spaceTimeWindowExtract first to create extraction data" << nl
            << exit(FatalError);
    }

    Info<< "Found extraction metadata: " << metadataPath << nl << endl;

    // Read metadata
    IFstream metaIs(metadataPath);
    dictionary metadata(metaIs);

    Info<< "Extraction metadata:" << nl
        << "    openfoamVersion = " << metadata.getOrDefault<word>("openfoamVersion", "unknown") << nl
        << "    solver = " << metadata.getOrDefault<word>("solver", "unknown") << nl
        << "    deltaT = " << metadata.get<scalar>("deltaT") << nl
        << "    adjustTimeStep = " << metadata.get<bool>("adjustTimeStep") << nl
        << "    extractionStartTime = " << metadata.get<scalar>("extractionStartTime") << nl
        << endl;

    // Check if boundary data is encrypted
    bool isEncrypted = metadata.getOrDefault<bool>("encrypted", false);
    std::vector<uint8_t> publicKey;
    std::vector<uint8_t> privateKey;

    if (isEncrypted)
    {
#ifdef FOAM_USE_SODIUM
        Info<< nl << "*** Boundary data is ENCRYPTED ***" << nl << endl;

        if (!sodiumCrypto::available())
        {
            FatalErrorInFunction
                << "Encrypted boundary data detected but libsodium failed to initialize" << nl
                << exit(FatalError);
        }

        // Prompt user for private key (no echo)
        std::string privateKeyB64 = sodiumCrypto::readPrivateKeyFromStdin
        (
            "Enter private key (base64): "
        );

        if (privateKeyB64.empty())
        {
            FatalErrorInFunction
                << "No private key provided" << nl
                << exit(FatalError);
        }

        privateKey = sodiumCrypto::fromBase64(privateKeyB64);
        if (privateKey.size() != sodiumCrypto::PRIVATE_KEY_SIZE)
        {
            FatalErrorInFunction
                << "Invalid private key size. Expected "
                << sodiumCrypto::PRIVATE_KEY_SIZE << " bytes, got "
                << privateKey.size() << " bytes" << nl
                << exit(FatalError);
        }

        // Derive public key from private key
        // For X25519 (used by sealed boxes), public = scalarmult_base(private)
        publicKey = sodiumCrypto::derivePublicKey(privateKey);

        Info<< "Public key derived. Decrypting boundary data..." << nl << endl;
#else
        FatalErrorInFunction
            << "Boundary data is encrypted but library was built without libsodium" << nl
            << "Rebuild with FOAM_USE_SODIUM=1 to enable decryption" << nl
            << exit(FatalError);
#endif
    }

    // Read source controlDict
    fileName sourceControlDictPath = sourceCase / "system" / "controlDict";
    IFstream sourceControlDictIs(sourceControlDictPath);
    dictionary sourceControlDict(sourceControlDictIs);

    // Find available time directories in boundaryData (for endTime)
    fileName patchBoundaryDir = boundaryDataDir / oldInternalPatchName;
    fileNameList timeDirs = readDir(patchBoundaryDir, fileName::DIRECTORY);

    // Collect all valid time directories with their scalar values
    DynamicList<Tuple2<scalar, word>> timeList;

    for (const fileName& dir : timeDirs)
    {
        scalar t;
        if (readScalar(dir, t))
        {
            timeList.append(Tuple2<scalar, word>(t, dir));
        }
    }

    if (timeList.size() < 5)
    {
        FatalErrorInFunction
            << "Need at least 5 time directories in boundaryData, found "
            << timeList.size() << nl
            << "    The spaceTimeWindow BC with cubic interpolation requires" << nl
            << "    4 timesteps for interpolation: t_{i-1}, t_i, t_{i+1}, t_{i+2}" << nl
            << "    Plus buffer at start (t_0, t_1) and end (t_{n-1}, t_n)" << nl
            << exit(FatalError);
    }

    // Sort by time value
    Foam::sort(timeList, [](const Tuple2<scalar, word>& a, const Tuple2<scalar, word>& b)
    {
        return a.first() < b.first();
    });

    // First timestep (minimum) - t_0, used only for cubic interpolation lookback
    word minTimeName = timeList.first().second();
    scalar minTime = timeList.first().first();

    // Last timestep (maximum) - for information only
    word maxTimeName = timeList.last().second();

    // Third timestep (t_2) - this is the usable startTime
    // Cubic interpolation in interval [t_2, t_3] needs t_1, t_2, t_3, t_4
    // so t_0 and t_1 are buffer for lookback
    word startTimeName = timeList[2].second();
    scalar startTime = timeList[2].first();

    // Third-to-last timestep (t_{n-2}) - this is the usable endTime
    // Cubic interpolation in interval [t_{n-3}, t_{n-2}] needs t_{n-4}, t_{n-3}, t_{n-2}, t_{n-1}
    // so t_{n-1} and t_n are buffer for lookahead
    word endTimeName = timeList[timeList.size() - 3].second();

    Info<< "Boundary data time range: " << minTimeName << " to " << maxTimeName << nl
        << "    Usable startTime: " << startTimeName << " (cubic needs t_0, t_1 for lookback)" << nl
        << "    Usable endTime: " << endTimeName << " (cubic needs t_{n-1}, t_n for lookahead)" << nl
        << endl;

    // If encrypted, decrypt all boundary data files now
    if (isEncrypted)
    {
#ifdef FOAM_USE_SODIUM
        Info<< nl << "Decrypting boundary data files..." << endl;

        label totalDecrypted = 0;

        // Process each time directory
        for (const auto& timePair : timeList)
        {
            const word& timeName = timePair.second();
            fileName timeDir = patchBoundaryDir / timeName;

            // Get all files in this time directory
            fileNameList files = readDir(timeDir, fileName::FILE);

            for (const fileName& f : files)
            {
                // Check if it's an encrypted file
                if (f.ends_with("." + sodiumCrypto::fileExtension))
                {
                    decryptBoundaryFile
                    (
                        timeDir / f,
                        timeDir,
                        publicKey,
                        privateKey,
                        false  // Not verbose for each file
                    );
                    totalDecrypted++;
                }
            }
        }

        Info<< "    Decrypted " << totalDecrypted << " files across "
            << timeList.size() << " timesteps" << nl << endl;
#endif
    }

    // Check for DVZT (delta-varint-temporal) compressed files and decompress them
    // DVZT files require sequential decompression because delta frames depend on previous timesteps
    // Check by looking for .dvzt files in the first time directory
    {
        fileName firstTimeDir = patchBoundaryDir / timeList.first().second();
        fileNameList files = readDir(firstTimeDir, fileName::FILE);
        bool hasDvzt = false;
        const word dvztExt = "." + deltaVarintTemporalCodec::fileExtension();

        for (const fileName& f : files)
        {
            if (f.ends_with(dvztExt))
            {
                hasDvzt = true;
                break;
            }
        }

        if (hasDvzt)
        {
            // Convert all DVZT files to DVZ (must be done sequentially in time order)
            // The resulting DVZ files are identical to what direct DVZ extraction produces
            convertDvztToDvzBoundaryData(patchBoundaryDir, timeList, 6, false);
        }
    }

    // Create target directories
    mkDir(extractDir / "system");
    mkDir(extractDir / "constant");

    // 1. Create controlDict
    Info<< "Creating system files..." << endl;
    writeControlDict(extractDir, metadata, sourceControlDict, startTimeName, endTimeName);

    // 2. Copy system files
    copyFile(sourceCase / "system" / "fvSchemes", extractDir / "system" / "fvSchemes");
    copyFile(sourceCase / "system" / "decomposeParDict", extractDir / "system" / "decomposeParDict");

    // Copy and modify fvSolution to add pRefPoint at center of extraction box
    {
        fileName srcFvSolution = sourceCase / "system" / "fvSolution";
        fileName dstFvSolution = extractDir / "system" / "fvSolution";

        if (isFile(srcFvSolution))
        {
            // Read source fvSolution
            IFstream srcIs(srcFvSolution);
            dictionary fvSolutionDict(srcIs);

            // Get pRefPoint from metadata (center of extraction box)
            vector pRefPoint = metadata.getOrDefault<vector>
            (
                "pRefPoint",
                0.5 * (metadata.get<vector>("boxMin") + metadata.get<vector>("boxMax"))
            );

            // Add pRefPoint to PIMPLE subdictionary (most common for incompressible)
            if (fvSolutionDict.found("PIMPLE"))
            {
                dictionary& pimpleDict = fvSolutionDict.subDict("PIMPLE");
                pimpleDict.set("pRefPoint", pRefPoint);
                pimpleDict.set("pRefValue", 0);
            }
            // Also check for SIMPLE
            if (fvSolutionDict.found("SIMPLE"))
            {
                dictionary& simpleDict = fvSolutionDict.subDict("SIMPLE");
                simpleDict.set("pRefPoint", pRefPoint);
                simpleDict.set("pRefValue", 0);
            }

            // Write modified fvSolution
            OFstream dstOs(dstFvSolution);
            dstOs << "FoamFile" << nl
                  << "{" << nl
                  << "    version     2.0;" << nl
                  << "    format      ascii;" << nl
                  << "    class       dictionary;" << nl
                  << "    object      fvSolution;" << nl
                  << "}" << nl
                  << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
                  << "// Modified by spaceTimeWindowInitCase - added pRefPoint" << nl
                  << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
                  << nl;
            fvSolutionDict.write(dstOs, false);

            Info<< "    Copied and modified: fvSolution" << nl
                << "        pRefPoint = " << pRefPoint << endl;
        }
        else
        {
            WarningInFunction
                << "Source file not found: " << srcFvSolution << endl;
        }
    }

    // 3. Copy constant files (mandatory for physical fidelity)
    Info<< nl << "Copying constant files (mandatory for fidelity)..." << endl;
    copyFile(sourceCase / "constant" / "turbulenceProperties", extractDir / "constant" / "turbulenceProperties");
    copyFile(sourceCase / "constant" / "transportProperties", extractDir / "constant" / "transportProperties");
    copyFile(sourceCase / "constant" / "momentumTransport", extractDir / "constant" / "momentumTransport");
    copyFile(sourceCase / "constant" / "thermophysicalProperties", extractDir / "constant" / "thermophysicalProperties");
    copyFile(sourceCase / "constant" / "g", extractDir / "constant" / "g");
    copyFile(sourceCase / "constant" / "RASProperties", extractDir / "constant" / "RASProperties");
    copyFile(sourceCase / "constant" / "LESProperties", extractDir / "constant" / "LESProperties");

    // 4. Read or create subset mesh
    // Check if mesh already exists (serial extraction) or needs to be created (parallel extraction)
    fileName meshDir = extractDir / "constant" / "polyMesh";
    fileName extractionBoxPath = meshDir / "extractionBox";

    bool needToCreateMesh = false;

    if (isFile(extractionBoxPath))
    {
        // Parallel extraction: mesh needs to be created from source case
        Info<< nl << "Found extractionBox - creating subset mesh from source case..." << endl;
        needToCreateMesh = true;
    }
    else if (!isFile(meshDir / "points") || !isFile(meshDir / "faces"))
    {
        FatalErrorInFunction
            << "No mesh found in " << meshDir << nl
            << "    Expected either a complete mesh (serial extraction) or" << nl
            << "    extractionBox file (parallel extraction)" << nl
            << exit(FatalError);
    }

    if (needToCreateMesh)
    {
        // Read extraction box parameters
        IFstream boxIs(extractionBoxPath);
        dictionary extractionBoxDict(boxIs);

        vector boxMin = extractionBoxDict.get<vector>("boxMin");
        vector boxMax = extractionBoxDict.get<vector>("boxMax");
        boundBox bb(boxMin, boxMax);

        Info<< "    Extraction box: " << boxMin << " to " << boxMax << endl;

        // Create Time object for source case
        Time sourceTime
        (
            Time::controlDictName,
            sourceCase,
            ".",
            false,  // enableFunctionObjects
            false   // enableLibs
        );

        // Read source mesh
        fvMesh sourceMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                sourceTime.constant(),
                sourceTime,
                IOobject::MUST_READ
            )
        );

        Info<< "    Source mesh has " << sourceMesh.nCells() << " cells" << endl;

        // Find cells where ALL vertices are inside the box
        // This ensures fvMeshSubset produces topologically valid cells
        // (no partial cells at boundaries that cause checkMesh failures)
        const pointField& points = sourceMesh.points();
        const cellList& cells = sourceMesh.cells();
        const faceList& faces = sourceMesh.faces();
        labelList cellsInBox;

        forAll(cells, celli)
        {
            const cell& c = cells[celli];
            bool allVerticesInside = true;

            // Check all vertices of all faces of this cell
            for (const label facei : c)
            {
                const face& f = faces[facei];
                for (const label pointi : f)
                {
                    if (!bb.contains(points[pointi]))
                    {
                        allVerticesInside = false;
                        break;
                    }
                }
                if (!allVerticesInside) break;
            }

            if (allVerticesInside)
            {
                cellsInBox.append(celli);
            }
        }

        Info<< "    Found " << cellsInBox.size() << " cells inside box" << endl;

        if (cellsInBox.empty())
        {
            FatalErrorInFunction
                << "No cells found inside extraction box" << nl
                << "    Box: " << boxMin << " to " << boxMax << nl
                << exit(FatalError);
        }

        // Create mesh subset
        fvMeshSubset meshSubset(sourceMesh);
        meshSubset.setCellSubset(cellsInBox);

        const fvMesh& subMesh = meshSubset.subMesh();
        const polyMesh& pm = subMesh;

        Info<< "    Created subset mesh with " << subMesh.nCells() << " cells" << endl;

        // Write the subset mesh to extractDir
        mkDir(meshDir);

        // Helper lambda to write FoamFile header
        auto writeHeader = [](Ostream& os, const word& className, const word& objectName)
        {
            os  << "FoamFile" << nl
                << "{" << nl
                << "    version     2.0;" << nl
                << "    format      ascii;" << nl
                << "    class       " << className << ";" << nl
                << "    object      " << objectName << ";" << nl
                << "}" << nl << nl;
        };

        // Points
        {
            OFstream os(meshDir / "points");
            writeHeader(os, "vectorField", "points");
            os << pm.points();
        }

        // Faces
        {
            OFstream os(meshDir / "faces");
            writeHeader(os, "faceList", "faces");
            os << pm.faces();
        }

        // Owner
        {
            OFstream os(meshDir / "owner");
            writeHeader(os, "labelList", "owner");
            os << pm.faceOwner();
        }

        // Neighbour
        {
            OFstream os(meshDir / "neighbour");
            writeHeader(os, "labelList", "neighbour");
            os << pm.faceNeighbour();
        }

        // Boundary - write oldInternalFaces and optionally outlet patch
        // If createOutlet is true, identify outlet faces and split the patch
        // Declare face lists outside the block so they can be used for field extraction
        labelList outletFaceIndices;
        labelList inletFaceIndices;
        {
            const polyBoundaryMesh& bm = pm.boundaryMesh();
            label oldInternalPatchI = -1;

            forAll(bm, patchi)
            {
                if (bm[patchi].name() == "oldInternalFaces")
                {
                    oldInternalPatchI = patchi;
                    break;
                }
            }

            if (oldInternalPatchI < 0)
            {
                FatalErrorInFunction
                    << "oldInternalFaces patch not found in subset mesh" << nl
                    << exit(FatalError);
            }

            const polyPatch& oifPatch = bm[oldInternalPatchI];
            label nTotalBndFaces = oifPatch.size();
            label startFace = pm.nInternalFaces();

            if (createOutlet)
            {
                // Calculate threshold position for outlet
                // Outlet faces are those with face centre in the direction of outletDirection
                // within outletFraction of the box extent

                const vectorField& faceCentres = pm.faceCentres();

                // Find extent of patch faces in outlet direction
                scalar minProj = GREAT;
                scalar maxProj = -GREAT;

                for (label i = 0; i < nTotalBndFaces; ++i)
                {
                    label facei = startFace + i;
                    scalar proj = faceCentres[facei] & outletDirection;
                    minProj = min(minProj, proj);
                    maxProj = max(maxProj, proj);
                }

                // Outlet threshold: faces in top outletFraction of extent
                scalar extent = maxProj - minProj;
                scalar outletThreshold = maxProj - outletFraction * extent;

                Info<< "    Outlet identification:" << nl
                    << "        Direction: " << outletDirection << nl
                    << "        Extent range: " << minProj << " to " << maxProj << nl
                    << "        Outlet threshold: > " << outletThreshold << endl;

                // Classify faces
                for (label i = 0; i < nTotalBndFaces; ++i)
                {
                    label facei = startFace + i;
                    scalar proj = faceCentres[facei] & outletDirection;

                    if (proj >= outletThreshold)
                    {
                        outletFaceIndices.append(i);
                    }
                    else
                    {
                        inletFaceIndices.append(i);
                    }
                }

                Info<< "        oldInternalFaces: " << inletFaceIndices.size() << " faces" << nl
                    << "        outlet: " << outletFaceIndices.size() << " faces" << endl;
            }
            else
            {
                // No outlet - all faces are oldInternalFaces
                for (label i = 0; i < nTotalBndFaces; ++i)
                {
                    inletFaceIndices.append(i);
                }
            }

            // If we have outlet faces, we need to reorder the mesh faces
            // so that oldInternalFaces come first, then outlet faces
            if (outletFaceIndices.size() > 0)
            {
                // Build reordering map: new face index -> old face index
                labelList faceOrder(nTotalBndFaces);
                label newIdx = 0;

                // First: oldInternalFaces (inlet)
                for (label oldIdx : inletFaceIndices)
                {
                    faceOrder[newIdx++] = oldIdx;
                }

                // Then: outlet faces
                for (label oldIdx : outletFaceIndices)
                {
                    faceOrder[newIdx++] = oldIdx;
                }

                // Build inverse map for field reordering
                labelList invFaceOrder(nTotalBndFaces);
                forAll(faceOrder, newi)
                {
                    invFaceOrder[faceOrder[newi]] = newi;
                }

                // Reorder faces - need to rewrite the faces file with reordered boundary faces
                faceList allFaces = pm.faces();
                labelList allOwner = pm.faceOwner();

                faceList reorderedFaces(allFaces.size());
                labelList reorderedOwner(allOwner.size());

                // Copy internal faces unchanged
                for (label i = 0; i < startFace; ++i)
                {
                    reorderedFaces[i] = allFaces[i];
                    reorderedOwner[i] = allOwner[i];
                }

                // Copy boundary faces in new order
                for (label newi = 0; newi < nTotalBndFaces; ++newi)
                {
                    label oldi = faceOrder[newi];
                    reorderedFaces[startFace + newi] = allFaces[startFace + oldi];
                    reorderedOwner[startFace + newi] = allOwner[startFace + oldi];
                }

                // Write reordered faces
                {
                    OFstream os(meshDir / "faces");
                    writeHeader(os, "faceList", "faces");
                    os << reorderedFaces;
                }

                // Write reordered owner
                {
                    OFstream os(meshDir / "owner");
                    writeHeader(os, "labelList", "owner");
                    os << reorderedOwner;
                }

                // Write boundary with two patches
                {
                    OFstream os(meshDir / "boundary");
                    writeHeader(os, "polyBoundaryMesh", "boundary");

                    os << "2" << nl << "(" << nl;

                    // oldInternalFaces patch
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl
                       << "        type            patch;" << nl
                       << "        nFaces          " << inletFaceIndices.size() << ";" << nl
                       << "        startFace       " << startFace << ";" << nl
                       << "    }" << nl;

                    // outlet patch
                    os << "    outlet" << nl
                       << "    {" << nl
                       << "        type            patch;" << nl
                       << "        nFaces          " << outletFaceIndices.size() << ";" << nl
                       << "        startFace       " << startFace + inletFaceIndices.size() << ";" << nl
                       << "    }" << nl;

                    os << ")" << nl;
                }

                // Store reordering info for field update
                // Write to a file so we can use it later when updating fields
                {
                    OFstream os(meshDir / "outletFaceReorder");
                    writeHeader(os, "labelList", "outletFaceReorder");
                    os << "// Maps new boundary face index to old boundary face index" << nl;
                    os << "// Use this to reorder boundary field values" << nl;
                    os << faceOrder;
                }

                // Update boundaryData to match new patch size (remove outlet faces)
                // The original boundaryData has data for all faces, but now
                // oldInternalFaces only has inletFaceIndices.size() faces
                Info<< nl << "    Updating boundaryData to match split patch..." << endl;

                // Process each time directory in boundaryData
                for (const auto& timePair : timeList)
                {
                    const word& timeName = timePair.second();
                    fileName timeDir = patchBoundaryDir / timeName;

                    // Get all field files in this time directory
                    fileNameList fieldFiles = readDir(timeDir, fileName::FILE);

                    for (const fileName& fieldFile : fieldFiles)
                    {
                        // Skip metadata
                        if (fieldFile == "extractionMetadata")
                        {
                            continue;
                        }

                        fileName fieldPath = timeDir / fieldFile;

                        // Read the file content to determine field type and read data
                        std::ifstream ifs(fieldPath.c_str());
                        std::string content((std::istreambuf_iterator<char>(ifs)),
                                           std::istreambuf_iterator<char>());
                        ifs.close();

                        // Find the data section after the FoamFile header
                        // Look for the closing brace of FoamFile, then find the size
                        std::size_t headerEnd = content.find('}');
                        if (headerEnd == std::string::npos)
                        {
                            WarningInFunction
                                << "Cannot find FoamFile header end in " << fieldFile << nl
                                << "    Skipping..." << endl;
                            continue;
                        }

                        // Check if it's a vector field by looking at data format
                        // Vectors have ((x y z) format after the size
                        // Or look for "vectorField" class
                        bool isVector = (content.find("vectorField") != std::string::npos
                                      || content.find("Field<vector>") != std::string::npos);

                        // If class is just "Field", check the data pattern
                        if (!isVector && content.find("class       Field;") != std::string::npos)
                        {
                            // Look for "((" pattern which indicates vector data
                            std::size_t dataStart = content.find('(', headerEnd);
                            if (dataStart != std::string::npos)
                            {
                                // Skip whitespace after first '('
                                std::size_t pos = dataStart + 1;
                                while (pos < content.size() && std::isspace(content[pos]))
                                {
                                    ++pos;
                                }
                                // If next char is '(' then it's a vector
                                if (pos < content.size() && content[pos] == '(')
                                {
                                    isVector = true;
                                }
                            }
                        }

                        // Create a string stream from the data portion
                        std::string dataSection = content.substr(headerEnd + 1);
                        IStringStream fieldIs(dataSection);

                        // Read the field data
                        if (isVector)
                        {
                            vectorField fullField(fieldIs);

                            if (fullField.size() != nTotalBndFaces)
                            {
                                WarningInFunction
                                    << "Field " << fieldFile << " has " << fullField.size()
                                    << " values but expected " << nTotalBndFaces << nl
                                    << "    Skipping..." << endl;
                                continue;
                            }

                            vectorField reducedField(inletFaceIndices.size());

                            // Extract only the faces that remain in oldInternalFaces
                            forAll(inletFaceIndices, i)
                            {
                                reducedField[i] = fullField[inletFaceIndices[i]];
                            }

                            // Write the reduced field
                            OFstream os(fieldPath);
                            os  << "FoamFile" << nl
                                << "{" << nl
                                << "    version     2.0;" << nl
                                << "    format      ascii;" << nl
                                << "    class       vectorField;" << nl
                                << "    object      " << fieldFile << ";" << nl
                                << "}" << nl << nl;
                            os << reducedField;
                        }
                        else
                        {
                            scalarField fullField(fieldIs);

                            if (fullField.size() != nTotalBndFaces)
                            {
                                WarningInFunction
                                    << "Field " << fieldFile << " has " << fullField.size()
                                    << " values but expected " << nTotalBndFaces << nl
                                    << "    Skipping..." << endl;
                                continue;
                            }

                            scalarField reducedField(inletFaceIndices.size());

                            // Extract only the faces that remain in oldInternalFaces
                            forAll(inletFaceIndices, i)
                            {
                                reducedField[i] = fullField[inletFaceIndices[i]];
                            }

                            // Write the reduced field
                            OFstream os(fieldPath);
                            os  << "FoamFile" << nl
                                << "{" << nl
                                << "    version     2.0;" << nl
                                << "    format      ascii;" << nl
                                << "    class       scalarField;" << nl
                                << "    object      " << fieldFile << ";" << nl
                                << "}" << nl << nl;
                            os << reducedField;
                        }
                    }
                }

                // Update the extractionMetadata to reflect new face count
                {
                    fileName metaPath = patchBoundaryDir / "extractionMetadata";
                    IFstream metaIs(metaPath);
                    dictionary metaDict(metaIs);
                    metaDict.set("nGlobalFaces", inletFaceIndices.size());
                    metaDict.set("originalNGlobalFaces", nTotalBndFaces);
                    metaDict.set("outletFacesRemoved", outletFaceIndices.size());

                    OFstream os(metaPath);
                    os  << "FoamFile" << nl
                        << "{" << nl
                        << "    version     2.0;" << nl
                        << "    format      ascii;" << nl
                        << "    class       dictionary;" << nl
                        << "    object      extractionMetadata;" << nl
                        << "}" << nl << nl;
                    metaDict.write(os, false);
                }

                Info<< "        Updated " << timeList.size() << " timesteps" << nl
                    << "        New oldInternalFaces size: " << inletFaceIndices.size() << endl;
            }
            else
            {
                // No reordering needed - write single patch boundary
                OFstream os(meshDir / "boundary");
                writeHeader(os, "polyBoundaryMesh", "boundary");

                os << "1" << nl << "(" << nl;
                os << "    oldInternalFaces" << nl
                   << "    {" << nl
                   << "        type            patch;" << nl
                   << "        nFaces          " << nTotalBndFaces << ";" << nl
                   << "        startFace       " << startFace << ";" << nl
                   << "    }" << nl;
                os << ")" << nl;
            }
        }

        Info<< "    Written subset mesh to " << meshDir << endl;

        // Also extract initial fields from source case at the reconstruction start time (t_1)
        // This is needed because parallel extraction doesn't write initial fields
        // (cell ordering would not match the mesh created here)
        // Note: We use startTime (t_1), not minTime (t_0), because reconstruction
        // starts at t_1 to ensure cubic interpolation has lookback data at t_0

        Info<< nl << "Extracting initial fields from source case at t=" << startTimeName << "..." << endl;

        // Set source time to reconstruction start time (t_1)
        sourceTime.setTime(startTime, 0);

        // Create initial field time directory
        fileName initFieldDir = extractDir / startTimeName;
        mkDir(initFieldDir);

        // Helper to write FoamFile header for fields
        auto writeFieldHeader = [](Ostream& os, const word& className, const word& objectName)
        {
            os  << "FoamFile" << nl
                << "{" << nl
                << "    version     2.0;" << nl
                << "    format      ascii;" << nl
                << "    class       " << className << ";" << nl
                << "    object      " << objectName << ";" << nl
                << "}" << nl << nl;
        };

        // Get field names for initial condition extraction
        // Use user-specified fields if provided, otherwise fall back to boundaryData fields
        wordList fieldsToExtract;
        wordList boundaryDataFields;

        // First, get the fields that are actually in boundaryData
        fileNameList firstTimeFiles;
        for (const fileName& dir : timeDirs)
        {
            scalar t;
            if (readScalar(dir, t) && mag(t - minTime) < SMALL)
            {
                firstTimeFiles = readDir(patchBoundaryDir / dir, fileName::FILE);
                break;
            }
        }

        for (const fileName& f : firstTimeFiles)
        {
            if (f != "points" && f != "extractionMetadata")
            {
                // Strip compression extensions if present (.dvz or .dvzt)
                word fieldName = f;
                const word dvztExt = "." + deltaVarintTemporalCodec::fileExtension();
                const word dvzExt = "." + deltaVarintCodec::fileExtension();
                if (fieldName.ends_with(dvztExt))
                {
                    fieldName = fieldName.substr(0, fieldName.size() - dvztExt.size());
                }
                else if (fieldName.ends_with(dvzExt))
                {
                    fieldName = fieldName.substr(0, fieldName.size() - dvzExt.size());
                }
                boundaryDataFields.append(fieldName);
            }
        }

        // Decide which fields to extract for initial conditions
        if (hasUserInitialFields)
        {
            fieldsToExtract = userInitialFields;
            Info<< "    Using user-specified initial fields: " << fieldsToExtract << nl
                << "    BoundaryData fields (for BC): " << boundaryDataFields << endl;
        }
        else
        {
            fieldsToExtract = boundaryDataFields;
            Info<< "    Using boundaryData fields for initial conditions: " << fieldsToExtract << endl;
        }

        // Read face reordering info if it exists (created when outlet patch is made)
        // Use the face order list we already have (still in scope)
        labelList faceReorder;
        label nOldInternalFaces = 0;
        label nOutletFaces = 0;
        bool hasOutlet = createOutlet && (outletFaceIndices.size() > 0);

        if (hasOutlet)
        {
            // Build faceReorder from inletFaceIndices and outletFaceIndices
            // (still in scope from boundary creation above)
            faceReorder.setSize(inletFaceIndices.size() + outletFaceIndices.size());
            label idx = 0;
            for (label oldIdx : inletFaceIndices)
            {
                faceReorder[idx++] = oldIdx;
            }
            for (label oldIdx : outletFaceIndices)
            {
                faceReorder[idx++] = oldIdx;
            }

            nOldInternalFaces = inletFaceIndices.size();
            nOutletFaces = outletFaceIndices.size();

            Info<< "    Using face reordering for outlet patch" << nl
                << "        oldInternalFaces: " << nOldInternalFaces << nl
                << "        outlet: " << nOutletFaces << endl;
        }

        Info<< "    Fields to extract: " << fieldsToExtract << endl;

        if (hasOutlet)
        {
            // With outlet patch, use cellMap() to extract real cell values
            // This avoids any interpolation while still providing accurate initial values
            Info<< "    Extracting fields with cellMap (outlet patch mode)..." << endl;

            const labelList& cellMap = meshSubset.cellMap();

            for (const word& fieldName : fieldsToExtract)
            {
                Info<< "    Processing field: " << fieldName << endl;

                bool isVector = (fieldName == "U" || fieldName == "U_0");

                if (isVector)
                {
                    IOobject vfIO
                    (
                        fieldName,
                        sourceTime.timeName(),
                        sourceMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (!vfIO.typeHeaderOk<volVectorField>(false))
                    {
                        WarningInFunction
                            << "Cannot find vector field: " << fieldName << nl
                            << "    Skipping..." << endl;
                        continue;
                    }

                    volVectorField vf(vfIO, sourceMesh);

                    // Extract internal field using cellMap (direct mapping, no interpolation)
                    vectorField subInternal(cellMap.size());
                    forAll(cellMap, subCelli)
                    {
                        subInternal[subCelli] = vf[cellMap[subCelli]];
                    }

                    OFstream os(initFieldDir / fieldName);
                    writeFieldHeader(os, "volVectorField", fieldName);

                    os << "dimensions      " << vf.dimensions() << ";" << nl << nl;
                    os << "internalField   nonuniform List<vector>" << nl
                       << subInternal.size() << nl << '(' << nl;
                    forAll(subInternal, i)
                    {
                        os << subInternal[i] << nl;
                    }
                    os << ')' << ';' << nl << nl;

                    // Check if field is in boundaryData (has time-varying BC data)
                    bool fieldInBoundaryData = boundaryDataFields.found(fieldName);

                    os << "boundaryField" << nl << '{' << nl;
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl;

                    if (fieldInBoundaryData)
                    {
                        // Field is in boundaryData - use spaceTimeWindow BC
                        if (useInletOutletBC)
                        {
                            os << "        type            spaceTimeWindowInletOutlet;" << nl
                               << "        dataDir         \"constant/boundaryData\";" << nl
                               << "        phi             phi;" << nl
                               << "        allowTimeInterpolation  true;" << nl
                               << "        timeInterpolationScheme cubic;" << nl
                               << "        value           uniform (0 0 0);" << nl;
                        }
                        else
                        {
                            os << "        type            spaceTimeWindow;" << nl
                               << "        dataDir         \"constant/boundaryData\";" << nl
                               << "        fixesValue      true;" << nl
                               << "        value           uniform (0 0 0);" << nl;
                        }
                    }
                    else
                    {
                        // Field NOT in boundaryData - use zeroGradient (no time-varying data)
                        os << "        type            zeroGradient;" << nl;
                        Info<< "        Note: " << fieldName << " not in boundaryData, using zeroGradient BC" << endl;
                    }

                    os << "    }" << nl;
                    os << "    outlet" << nl
                       << "    {" << nl
                       << "        type            inletOutlet;" << nl
                       << "        inletValue      uniform (0 0 0);" << nl
                       << "        value           uniform (0 0 0);" << nl
                       << "    }" << nl;
                    os << '}' << nl;

                    Info<< "    Extracted: " << fieldName << " (vector, "
                        << subInternal.size() << " cells)" << endl;
                }
                else
                {
                    IOobject sfIO
                    (
                        fieldName,
                        sourceTime.timeName(),
                        sourceMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (!sfIO.typeHeaderOk<volScalarField>(false))
                    {
                        WarningInFunction
                            << "Cannot find scalar field: " << fieldName << nl
                            << "    Skipping..." << endl;
                        continue;
                    }

                    volScalarField sf(sfIO, sourceMesh);

                    // Extract internal field using cellMap (direct mapping, no interpolation)
                    scalarField subInternal(cellMap.size());
                    forAll(cellMap, subCelli)
                    {
                        subInternal[subCelli] = sf[cellMap[subCelli]];
                    }

                    OFstream os(initFieldDir / fieldName);
                    writeFieldHeader(os, "volScalarField", fieldName);

                    os << "dimensions      " << sf.dimensions() << ";" << nl << nl;
                    os << "internalField   nonuniform List<scalar>" << nl
                       << subInternal.size() << nl << '(' << nl;
                    forAll(subInternal, i)
                    {
                        os << subInternal[i] << nl;
                    }
                    os << ')' << ';' << nl << nl;

                    // Check if field is in boundaryData (has time-varying BC data)
                    bool fieldInBoundaryData = boundaryDataFields.found(fieldName);

                    os << "boundaryField" << nl << '{' << nl;
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl;

                    if (fieldInBoundaryData && !useInletOutletBC)
                    {
                        // Field is in boundaryData and not using inletOutletBC - use spaceTimeWindow
                        os << "        type            spaceTimeWindow;" << nl
                           << "        dataDir         \"constant/boundaryData\";" << nl
                           << "        fixesValue      true;" << nl
                           << "        value           uniform 0;" << nl;
                    }
                    else
                    {
                        // Field NOT in boundaryData, or using inletOutletBC - use zeroGradient
                        os << "        type            zeroGradient;" << nl;
                        if (!fieldInBoundaryData)
                        {
                            Info<< "        Note: " << fieldName << " not in boundaryData, using zeroGradient BC" << endl;
                        }
                    }

                    os << "    }" << nl;
                    os << "    outlet" << nl
                       << "    {" << nl
                       << "        type            zeroGradient;" << nl
                       << "    }" << nl;
                    os << '}' << nl;

                    Info<< "    Extracted: " << fieldName << " (scalar, "
                        << subInternal.size() << " cells)" << endl;
                }
            }
        }
        else
        {
            // No outlet - use cellMap() for direct cell value extraction (no interpolation)
            // The subset mesh has only oldInternalFaces patch
            Info<< "    Extracting fields with cellMap (no outlet)..." << endl;

            const labelList& cellMap = meshSubset.cellMap();

            for (const word& fieldName : fieldsToExtract)
            {
                Info<< "    Processing field: " << fieldName << endl;

                bool isVector = (fieldName == "U" || fieldName == "U_0");

                // Check if field is in boundaryData (has time-varying BC data)
                bool fieldInBoundaryData = boundaryDataFields.found(fieldName);

                if (!isVector)
                {
                    IOobject sfIO
                    (
                        fieldName,
                        sourceTime.timeName(),
                        sourceMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (!sfIO.typeHeaderOk<volScalarField>(false))
                    {
                        WarningInFunction
                            << "Cannot find scalar field: " << fieldName << nl
                            << "    Skipping..." << endl;
                        continue;
                    }

                    volScalarField sf(sfIO, sourceMesh);

                    // Extract internal field using cellMap (direct mapping, no interpolation)
                    scalarField subInternal(cellMap.size());
                    forAll(cellMap, subCelli)
                    {
                        subInternal[subCelli] = sf[cellMap[subCelli]];
                    }

                    OFstream os(initFieldDir / fieldName);
                    writeFieldHeader(os, "volScalarField", fieldName);

                    os << "dimensions      " << sf.dimensions() << ";" << nl << nl;
                    os << "internalField   nonuniform List<scalar>" << nl
                       << subInternal.size() << nl << '(' << nl;
                    forAll(subInternal, i)
                    {
                        os << subInternal[i] << nl;
                    }
                    os << ')' << ';' << nl << nl;

                    // Write boundary field - only oldInternalFaces exists in subset mesh
                    os << "boundaryField" << nl << '{' << nl;
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl;

                    if (fieldInBoundaryData && !useInletOutletBC)
                    {
                        // Field is in boundaryData and not using inletOutletBC - use spaceTimeWindow
                        os << "        type            spaceTimeWindow;" << nl
                           << "        dataDir         \"constant/boundaryData\";" << nl
                           << "        fixesValue      true;" << nl
                           << "        value           uniform 0;" << nl;
                    }
                    else
                    {
                        // Field NOT in boundaryData, or using inletOutletBC - use zeroGradient
                        os << "        type            zeroGradient;" << nl;
                    }
                    os << "    }" << nl;
                    os << '}' << nl;

                    if (!fieldInBoundaryData)
                    {
                        Info<< "        Note: " << fieldName << " not in boundaryData, using zeroGradient BC" << endl;
                    }

                    Info<< "    Extracted: " << fieldName << " (scalar, "
                        << subInternal.size() << " cells)" << endl;
                }
                else
                {
                    IOobject vfIO
                    (
                        fieldName,
                        sourceTime.timeName(),
                        sourceMesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    );

                    if (!vfIO.typeHeaderOk<volVectorField>(false))
                    {
                        WarningInFunction
                            << "Cannot find vector field: " << fieldName << nl
                            << "    Skipping..." << endl;
                        continue;
                    }

                    volVectorField vf(vfIO, sourceMesh);

                    // Extract internal field using cellMap (direct mapping, no interpolation)
                    vectorField subInternal(cellMap.size());
                    forAll(cellMap, subCelli)
                    {
                        subInternal[subCelli] = vf[cellMap[subCelli]];
                    }

                    OFstream os(initFieldDir / fieldName);
                    writeFieldHeader(os, "volVectorField", fieldName);

                    os << "dimensions      " << vf.dimensions() << ";" << nl << nl;
                    os << "internalField   nonuniform List<vector>" << nl
                       << subInternal.size() << nl << '(' << nl;
                    forAll(subInternal, i)
                    {
                        os << subInternal[i] << nl;
                    }
                    os << ')' << ';' << nl << nl;

                    // Write boundary field - only oldInternalFaces exists in subset mesh
                    os << "boundaryField" << nl << '{' << nl;
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl;

                    if (fieldInBoundaryData)
                    {
                        // Field is in boundaryData - use spaceTimeWindow BC
                        if (useInletOutletBC)
                        {
                            // With -inletOutletBC: U uses spaceTimeWindowInletOutlet
                            os << "        type            spaceTimeWindowInletOutlet;" << nl
                               << "        dataDir         \"constant/boundaryData\";" << nl
                               << "        phi             phi;" << nl
                               << "        allowTimeInterpolation  true;" << nl
                               << "        timeInterpolationScheme cubic;" << nl
                               << "        value           uniform (0 0 0);" << nl;
                        }
                        else
                        {
                            // Default: all fields use spaceTimeWindow
                            os << "        type            spaceTimeWindow;" << nl
                               << "        dataDir         \"constant/boundaryData\";" << nl
                               << "        fixesValue      true;" << nl
                               << "        value           uniform (0 0 0);" << nl;
                        }
                    }
                    else
                    {
                        // Field NOT in boundaryData - use zeroGradient (no time-varying data)
                        os << "        type            zeroGradient;" << nl;
                        Info<< "        Note: " << fieldName << " not in boundaryData, using zeroGradient BC" << endl;
                    }
                    os << "    }" << nl;
                    os << '}' << nl;

                    Info<< "    Extracted: " << fieldName << " (vector, "
                        << subInternal.size() << " cells)" << endl;
                }
            }
        }

        // Keep extractionBox file for reference (documents extraction parameters)
    }

    Info<< nl << "Reading subset mesh..." << endl;

    // Create a minimal Time object for the extract case
    Time runTime
    (
        Time::controlDictName,
        extractDir,
        ".",
        false,  // enableFunctionObjects
        false   // enableLibs
    );

    // Read the mesh
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::MUST_READ
        )
    );

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    Info<< "Subset mesh patches:" << endl;
    forAll(bMesh, patchi)
    {
        Info<< "    " << bMesh[patchi].name()
            << " (" << bMesh[patchi].type() << "): "
            << bMesh[patchi].size() << " faces" << endl;
    }

    // Apply mass flux correction if requested
    // Must be done BEFORE outlet faces are removed from boundaryData
    // so that the correction is computed over all original boundary faces
    if (correctMassFlux)
    {
        // Find the oldInternalFaces patch to get face area vectors
        label oifPatchI = bMesh.findPatchID("oldInternalFaces");
        if (oifPatchI >= 0)
        {
            // Get face area vectors for oldInternalFaces patch
            const polyPatch& oifPatch = bMesh[oifPatchI];
            vectorField faceSf(oifPatch.size());

            // Compute Sf from face areas and normals
            const vectorField& faceAreas = mesh.faceAreas();
            label startFace = oifPatch.start();

            forAll(oifPatch, i)
            {
                faceSf[i] = faceAreas[startFace + i];
            }

            // Apply mass flux correction to all timesteps in boundaryData
            correctMassFluxInBoundaryData(patchBoundaryDir, timeList, faceSf);
        }
        else
        {
            WarningInFunction
                << "oldInternalFaces patch not found, skipping mass flux correction" << endl;
        }
    }

    // 5. Find fields in boundaryData (for BC assignment)
    fileNameList firstTimeFiles;
    for (const fileName& dir : timeDirs)
    {
        scalar t;
        if (readScalar(dir, t) && mag(t - minTime) < SMALL)
        {
            firstTimeFiles = readDir(patchBoundaryDir / dir, fileName::FILE);
            break;
        }
    }

    wordList boundaryDataFields;
    for (const fileName& f : firstTimeFiles)
    {
        if (f != "points" && f != "extractionMetadata")
        {
            // Strip compression extensions if present (.dvz or .dvzt)
            word fieldName = f;
            const word dvztExt = "." + deltaVarintTemporalCodec::fileExtension();
            const word dvzExt = "." + deltaVarintCodec::fileExtension();
            if (fieldName.ends_with(dvztExt))
            {
                fieldName = fieldName.substr(0, fieldName.size() - dvztExt.size());
            }
            else if (fieldName.ends_with(dvzExt))
            {
                fieldName = fieldName.substr(0, fieldName.size() - dvzExt.size());
            }
            boundaryDataFields.append(fieldName);
        }
    }

    Info<< nl << "BoundaryData fields (for BC): " << boundaryDataFields << endl;

    // 6. Update initial fields with spaceTimeWindow BC
    // The spaceTimeWindowExtract already wrote subset fields to extractDir
    // We just need to update the boundary conditions on oldInternalFaces
    Info<< nl << "Updating initial fields with spaceTimeWindow BC..." << endl;

    // Find the initial time directory in extract case (written by spaceTimeWindowExtract)
    fileNameList extractTimeDirs = readDir(extractDir, fileName::DIRECTORY);
    word initialTimeDir;

    for (const fileName& dir : extractTimeDirs)
    {
        scalar t;
        if (readScalar(dir, t))
        {
            initialTimeDir = dir;
            break;  // Use first time directory found
        }
    }

    if (initialTimeDir.empty())
    {
        FatalErrorInFunction
            << "No time directory found in extract directory: " << extractDir << nl
            << "    spaceTimeWindowExtract should have created initial fields" << nl
            << exit(FatalError);
    }

    Info<< "    Using initial time: " << initialTimeDir << endl;

    // Discover initial fields from the time directory
    // (these may differ from boundaryDataFields if initialFields was configured separately)
    wordList initialFields;
    fileNameList initFieldFiles = readDir(extractDir / initialTimeDir, fileName::FILE);
    for (const fileName& f : initFieldFiles)
    {
        // Only include files that look like OpenFOAM field files
        // Skip compressed files, backup files, etc.
        if (!f.ends_with("~") && !f.ends_with(".gz") && f != "uniform")
        {
            initialFields.append(f);
        }
    }

    Info<< "    Initial fields in " << initialTimeDir << "/: " << initialFields << endl;

    // Update each extracted field - rewrite with only oldInternalFaces in boundaryField
    // Serial extraction writes fields with all source mesh patches, but the subset mesh
    // only has oldInternalFaces. We need to rewrite the entire boundaryField section.
    // - Fields in boundaryDataFields: use spaceTimeWindow or spaceTimeWindowInletOutlet BC
    // - Fields NOT in boundaryDataFields: use zeroGradient BC
    for (const word& fieldName : initialFields)
    {
        fileName fieldPath = extractDir / initialTimeDir / fieldName;

        if (!isFile(fieldPath))
        {
            WarningInFunction
                << "Field not found in extract directory: " << fieldPath << nl
                << "    Skipping..." << endl;
            continue;
        }

        // Read the entire file
        std::ifstream ifs(fieldPath.c_str());
        std::ostringstream buffer;
        buffer << ifs.rdbuf();
        std::string content = buffer.str();
        ifs.close();

        // Determine field type from class in header
        word fieldType = "scalar";
        if (content.find("volVectorField") != std::string::npos)
        {
            fieldType = "vector";
        }

        // Check if field has time-varying BC data in boundaryData
        bool fieldInBoundaryData = boundaryDataFields.found(fieldName);

        // Find boundaryField section and replace the entire section
        // to only include oldInternalFaces (remove walls, inlet, outlet, etc.)
        std::string searchBoundary = "boundaryField";
        std::size_t boundaryPos = content.find(searchBoundary);

        if (boundaryPos != std::string::npos)
        {
            // Find the opening brace of boundaryField
            std::size_t bfBraceStart = content.find('{', boundaryPos);
            if (bfBraceStart != std::string::npos)
            {
                // Find the matching closing brace of boundaryField
                int braceCount = 1;
                std::size_t bfBraceEnd = bfBraceStart + 1;
                while (bfBraceEnd < content.size() && braceCount > 0)
                {
                    if (content[bfBraceEnd] == '{') braceCount++;
                    else if (content[bfBraceEnd] == '}') braceCount--;
                    bfBraceEnd++;
                }

                // Build the new boundaryField section with only oldInternalFaces
                std::string newBoundaryField = "boundaryField\n{\n    oldInternalFaces\n    {\n";
                word bcType;

                if (fieldInBoundaryData)
                {
                    // Field is in boundaryData - use spaceTimeWindow BC
                    if (fieldType == "vector")
                    {
                        if (useInletOutletBC)
                        {
                            // With -inletOutletBC: U uses spaceTimeWindowInletOutlet
                            newBoundaryField += "        type            spaceTimeWindowInletOutlet;\n";
                            newBoundaryField += "        dataDir         \"constant/boundaryData\";\n";
                            newBoundaryField += "        phi             phi;\n";
                            newBoundaryField += "        allowTimeInterpolation  true;\n";
                            newBoundaryField += "        timeInterpolationScheme cubic;\n";
                            newBoundaryField += "        value           uniform (0 0 0);\n";
                            bcType = "spaceTimeWindowInletOutlet";
                        }
                        else
                        {
                            // Default: use spaceTimeWindow (Dirichlet)
                            newBoundaryField += "        type            spaceTimeWindow;\n";
                            newBoundaryField += "        dataDir         \"constant/boundaryData\";\n";
                            newBoundaryField += "        fixesValue      true;\n";
                            newBoundaryField += "        value           uniform (0 0 0);\n";
                            bcType = "spaceTimeWindow";
                        }
                    }
                    else
                    {
                        // Scalar field
                        if (useInletOutletBC)
                        {
                            // With -inletOutletBC: scalar fields use zeroGradient
                            newBoundaryField += "        type            zeroGradient;\n";
                            bcType = "zeroGradient";
                        }
                        else
                        {
                            // Default: use spaceTimeWindow (Dirichlet)
                            newBoundaryField += "        type            spaceTimeWindow;\n";
                            newBoundaryField += "        dataDir         \"constant/boundaryData\";\n";
                            newBoundaryField += "        fixesValue      true;\n";
                            newBoundaryField += "        value           uniform 0;\n";
                            bcType = "spaceTimeWindow";
                        }
                    }
                }
                else
                {
                    // Field NOT in boundaryData - use zeroGradient (no time-varying data available)
                    newBoundaryField += "        type            zeroGradient;\n";
                    bcType = "zeroGradient (not in boundaryData)";
                }

                newBoundaryField += "    }\n}";

                // Replace the entire boundaryField section
                content = content.substr(0, boundaryPos) + newBoundaryField + content.substr(bfBraceEnd);

                // Write back
                std::ofstream ofs(fieldPath.c_str());
                ofs << content;
                ofs.close();

                Info<< "    Updated: " << fieldName << " (oldInternalFaces -> " << bcType << ")" << endl;
            }
        }
    }

    // Note: All fields in the initial time directory are now handled by the loop above.
    // Fields with boundaryData get spaceTimeWindow BC; fields without get zeroGradient.

    Info<< nl
        << "========================================" << nl
        << "Case initialization complete!" << nl
        << "========================================" << nl
        << nl
        << "Next steps:" << nl
        << "  1. (Optional) Adjust writeInterval in system/controlDict" << nl
        << "  2. Review boundary conditions in " << initialTimeDir << "/" << nl
        << "  3. Run the solver:" << nl
        << "     " << metadata.getOrDefault<word>("solver", "pimpleFoam") << nl
        << nl;

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
