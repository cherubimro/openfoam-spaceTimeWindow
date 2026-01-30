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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "IOdictionary.H"
#include "IFstream.H"
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

    Info<< "spaceTimeWindowInitCase" << nl
        << "    Source case:    " << sourceCase << nl
        << "    Extract dir:    " << extractDir << nl
        << endl;

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

    if (timeList.size() < 2)
    {
        FatalErrorInFunction
            << "Need at least 2 time directories in boundaryData, found "
            << timeList.size() << nl
            << "    The spaceTimeWindow BC requires data beyond the current time" << nl
            << exit(FatalError);
    }

    // Sort by time value
    Foam::sort(timeList, [](const Tuple2<scalar, word>& a, const Tuple2<scalar, word>& b)
    {
        return a.first() < b.first();
    });

    // First timestep (minimum)
    word minTimeName = timeList.first().second();
    scalar minTime = timeList.first().first();

    // Last timestep (maximum) - for information only
    word maxTimeName = timeList.last().second();

    // Second-to-last timestep - this is the usable endTime
    // because the BC reads ahead to the next timestep
    word endTimeName = timeList[timeList.size() - 2].second();

    Info<< "Boundary data time range: " << minTimeName << " to " << maxTimeName << nl
        << "    Usable endTime: " << endTimeName << " (BC requires next timestep data)" << nl
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

    // Create target directories
    mkDir(extractDir / "system");
    mkDir(extractDir / "constant");

    // 1. Create controlDict
    Info<< "Creating system files..." << endl;
    writeControlDict(extractDir, metadata, sourceControlDict, minTimeName, endTimeName);

    // 2. Copy system files
    copyFile(sourceCase / "system" / "fvSchemes", extractDir / "system" / "fvSchemes");
    copyFile(sourceCase / "system" / "fvSolution", extractDir / "system" / "fvSolution");
    copyFile(sourceCase / "system" / "decomposeParDict", extractDir / "system" / "decomposeParDict");

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

        // Find cells inside the box
        const vectorField& cellCentres = sourceMesh.cellCentres();
        labelList cellsInBox;

        forAll(cellCentres, celli)
        {
            if (bb.contains(cellCentres[celli]))
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

        // Boundary - write only oldInternalFaces patch
        {
            OFstream os(meshDir / "boundary");
            writeHeader(os, "polyBoundaryMesh", "boundary");

            const polyBoundaryMesh& bm = pm.boundaryMesh();

            os << "1" << nl << "(" << nl;
            forAll(bm, patchi)
            {
                if (bm[patchi].name() == "oldInternalFaces")
                {
                    os << "    oldInternalFaces" << nl
                       << "    {" << nl
                       << "        type            patch;" << nl
                       << "        nFaces          " << bm[patchi].size() << ";" << nl
                       << "        startFace       " << bm[patchi].start() << ";" << nl
                       << "    }" << nl;
                    break;
                }
            }
            os << ")" << nl;
        }

        Info<< "    Written subset mesh to " << meshDir << endl;

        // Also extract initial fields from source case at extraction start time
        // This is needed because parallel extraction doesn't write initial fields
        // (cell ordering would not match the mesh created here)
        scalar extractionStartTime = metadata.get<scalar>("extractionStartTime");
        word startTimeName = minTimeName;

        Info<< nl << "Extracting initial fields from source case at t=" << startTimeName << "..." << endl;

        // Set source time to extraction start time
        sourceTime.setTime(extractionStartTime, 0);

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

        // Get field names from boundaryData first time directory
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

        wordList fieldsToExtract;
        for (const fileName& f : firstTimeFiles)
        {
            if (f != "points" && f != "extractionMetadata")
            {
                // Strip .dvz extension if present
                word fieldName = f;
                const word dvzExt = "." + deltaVarintCodec::fileExtension();
                if (fieldName.ends_with(dvzExt))
                {
                    fieldName = fieldName.substr(0, fieldName.size() - dvzExt.size());
                }
                fieldsToExtract.append(fieldName);
            }
        }

        // Extract each field that appears in boundaryData
        for (const word& fieldName : fieldsToExtract)
        {
            // Try to read scalar field from source
            IOobject sfHeader
            (
                fieldName,
                sourceTime.timeName(),
                sourceMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (sfHeader.typeHeaderOk<volScalarField>(true))
            {
                volScalarField sf(sfHeader, sourceMesh);

                // Interpolate to subset
                tmp<volScalarField> tsubField = meshSubset.interpolate(sf);
                const volScalarField& subField = tsubField();

                // Write to extract directory
                OFstream os(initFieldDir / fieldName);
                writeFieldHeader(os, "volScalarField", fieldName);

                os << "dimensions      " << sf.dimensions() << ";" << nl << nl;
                os << "internalField   nonuniform List<scalar>" << nl
                   << subField.primitiveField().size() << nl << '(' << nl;
                forAll(subField.primitiveField(), i)
                {
                    os << subField.primitiveField()[i] << nl;
                }
                os << ')' << ';' << nl << nl;

                // Write boundary field
                os << "boundaryField" << nl << '{' << nl;
                forAll(subMesh.boundary(), patchi)
                {
                    const fvPatch& patch = subMesh.boundary()[patchi];
                    os << "    " << patch.name() << nl
                       << "    {" << nl
                       << "        type            calculated;" << nl
                       << "        value           nonuniform List<scalar>" << nl
                       << "        " << subField.boundaryField()[patchi].size() << nl
                       << "        (" << nl;
                    forAll(subField.boundaryField()[patchi], i)
                    {
                        os << "        " << subField.boundaryField()[patchi][i] << nl;
                    }
                    os << "        );" << nl
                       << "    }" << nl;
                }
                os << '}' << nl;

                Info<< "    Extracted: " << fieldName << " (scalar)" << endl;
                continue;
            }

            // Try to read vector field from source
            IOobject vfHeader
            (
                fieldName,
                sourceTime.timeName(),
                sourceMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );

            if (vfHeader.typeHeaderOk<volVectorField>(true))
            {
                volVectorField vf(vfHeader, sourceMesh);

                // Interpolate to subset
                tmp<volVectorField> tsubField = meshSubset.interpolate(vf);
                const volVectorField& subField = tsubField();

                // Write to extract directory
                OFstream os(initFieldDir / fieldName);
                writeFieldHeader(os, "volVectorField", fieldName);

                os << "dimensions      " << vf.dimensions() << ";" << nl << nl;
                os << "internalField   nonuniform List<vector>" << nl
                   << subField.primitiveField().size() << nl << '(' << nl;
                forAll(subField.primitiveField(), i)
                {
                    os << subField.primitiveField()[i] << nl;
                }
                os << ')' << ';' << nl << nl;

                // Write boundary field
                os << "boundaryField" << nl << '{' << nl;
                forAll(subMesh.boundary(), patchi)
                {
                    const fvPatch& patch = subMesh.boundary()[patchi];
                    os << "    " << patch.name() << nl
                       << "    {" << nl
                       << "        type            calculated;" << nl
                       << "        value           nonuniform List<vector>" << nl
                       << "        " << subField.boundaryField()[patchi].size() << nl
                       << "        (" << nl;
                    forAll(subField.boundaryField()[patchi], i)
                    {
                        os << "        " << subField.boundaryField()[patchi][i] << nl;
                    }
                    os << "        );" << nl
                       << "    }" << nl;
                }
                os << '}' << nl;

                Info<< "    Extracted: " << fieldName << " (vector)" << endl;
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

    // 5. Find extracted fields from boundaryData
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

    wordList extractedFields;
    for (const fileName& f : firstTimeFiles)
    {
        if (f != "points" && f != "extractionMetadata")
        {
            // Strip .dvz extension if present
            word fieldName = f;
            const word dvzExt = "." + deltaVarintCodec::fileExtension();
            if (fieldName.ends_with(dvzExt))
            {
                fieldName = fieldName.substr(0, fieldName.size() - dvzExt.size());
            }
            extractedFields.append(fieldName);
        }
    }

    Info<< nl << "Extracted fields: " << extractedFields << endl;

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

    // Update each extracted field - modify oldInternalFaces BC to spaceTimeWindow
    // The extraction wrote fields with "calculated" BC on oldInternalFaces
    // We need to replace it with "spaceTimeWindow" BC
    for (const word& fieldName : extractedFields)
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

        // Find and replace the oldInternalFaces BC
        // Pattern to find: oldInternalFaces { type calculated; ... }
        std::string searchStart = "oldInternalFaces";
        std::size_t pos = content.find(searchStart);

        if (pos != std::string::npos)
        {
            // Find the opening brace
            std::size_t braceStart = content.find('{', pos);
            if (braceStart != std::string::npos)
            {
                // Find the matching closing brace
                int braceCount = 1;
                std::size_t braceEnd = braceStart + 1;
                while (braceEnd < content.size() && braceCount > 0)
                {
                    if (content[braceEnd] == '{') braceCount++;
                    else if (content[braceEnd] == '}') braceCount--;
                    braceEnd++;
                }

                // Build the replacement BC
                std::string newBC;
                if (fieldType == "vector")
                {
                    newBC = "oldInternalFaces\n    {\n        type            spaceTimeWindow;\n        dataDir         \"constant/boundaryData\";\n        fixesValue      true;\n        value           uniform (0 0 0);\n    }";
                }
                else
                {
                    newBC = "oldInternalFaces\n    {\n        type            spaceTimeWindow;\n        dataDir         \"constant/boundaryData\";\n        fixesValue      true;\n        value           uniform 0;\n    }";
                }

                // Replace the old BC with the new one
                content = content.substr(0, pos) + newBC + content.substr(braceEnd);

                // Write back
                std::ofstream ofs(fieldPath.c_str());
                ofs << content;
                ofs.close();

                Info<< "    Updated: " << fieldName << " (oldInternalFaces -> spaceTimeWindow)" << endl;
            }
        }
    }

    // Also handle turbulence fields that exist in extract dir but weren't in extracted list
    // (e.g., if nut was extracted separately)
    wordList commonTurbFields = {"nut", "nuTilda", "k", "epsilon", "omega"};

    for (const word& fieldName : commonTurbFields)
    {
        if (!extractedFields.found(fieldName))
        {
            // Check if the field already exists in extract dir
            fileName fieldPath = extractDir / initialTimeDir / fieldName;
            if (isFile(fieldPath))
            {
                // Apply same BC update
                std::ifstream ifs(fieldPath.c_str());
                std::ostringstream buffer;
                buffer << ifs.rdbuf();
                std::string content = buffer.str();
                ifs.close();

                std::string searchStart = "oldInternalFaces";
                std::size_t pos = content.find(searchStart);

                if (pos != std::string::npos)
                {
                    std::size_t braceStart = content.find('{', pos);
                    if (braceStart != std::string::npos)
                    {
                        int braceCount = 1;
                        std::size_t braceEnd = braceStart + 1;
                        while (braceEnd < content.size() && braceCount > 0)
                        {
                            if (content[braceEnd] == '{') braceCount++;
                            else if (content[braceEnd] == '}') braceCount--;
                            braceEnd++;
                        }

                        std::string newBC = "oldInternalFaces\n    {\n        type            spaceTimeWindow;\n        dataDir         \"constant/boundaryData\";\n        fixesValue      true;\n        value           uniform 0;\n    }";
                        content = content.substr(0, pos) + newBC + content.substr(braceEnd);

                        std::ofstream ofs(fieldPath.c_str());
                        ofs << content;
                        ofs.close();

                        Info<< "    Updated: " << fieldName << " (oldInternalFaces -> spaceTimeWindow)" << endl;
                    }
                }
            }
        }
    }

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
