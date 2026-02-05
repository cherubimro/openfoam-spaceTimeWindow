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
    compareSubsetFields

Description
    Compare internal fields between original source case and recalculated
    subset case. The subset mesh must be extracted from the source using
    the spaceTimeWindow extraction (all vertices of a cell inside box).

    Reads only internal field data (ignores boundary fields) to avoid
    issues with boundary field sizes.

Usage
    compareSubsetFields -sourceCase <path> [options]

    Options:
        -sourceCase <path>  : Path to original source case (required)
        -fields <list>      : Fields to compare (default: "U p nut")
        -time <ranges>      : Time selection (default: all common times)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "boundBox.H"
#include "IFstream.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find cells where ALL vertices are inside the bounding box
// This matches the extraction algorithm exactly
labelList findCellsInBox
(
    const fvMesh& mesh,
    const boundBox& bb
)
{
    const pointField& points = mesh.points();
    const cellList& cells = mesh.cells();
    const faceList& faces = mesh.faces();

    DynamicList<label> cellsInBox;

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

    return cellsInBox;
}


// Read internal field from file (bypasses boundary field reading)
template<class Type>
tmp<Field<Type>> readInternalField
(
    const fileName& fieldPath,
    const label expectedSize
)
{
    IFstream is(fieldPath);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open file: " << fieldPath
            << exit(FatalError);
    }

    // Read the dictionary-style field file
    dictionary fieldDict(is);

    // Get the internalField entry
    const entry& internalEntry = fieldDict.lookupEntry("internalField", keyType::LITERAL);

    // Parse the internal field
    ITstream& stream = internalEntry.stream();

    word fieldType;
    stream >> fieldType;

    auto tfield = tmp<Field<Type>>::New();
    Field<Type>& field = tfield.ref();

    if (fieldType == "uniform")
    {
        Type uniformValue;
        stream >> uniformValue;
        field.setSize(expectedSize, uniformValue);
    }
    else if (fieldType == "nonuniform")
    {
        stream >> field;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown field type: " << fieldType
            << " in file " << fieldPath
            << exit(FatalError);
    }

    return tfield;
}


// Determine field type from file header
word getFieldType(const fileName& fieldPath)
{
    IFstream is(fieldPath);
    if (!is.good())
    {
        return word::null;
    }

    // Read just enough to find the class
    dictionary header;
    header.read(is, false);  // Read without checking for FoamFile

    // Reopen and look for FoamFile header
    IFstream is2(fieldPath);
    token tok;

    while (is2.good())
    {
        is2 >> tok;
        if (tok.isWord() && tok.wordToken() == "FoamFile")
        {
            // Read the header dictionary
            dictionary foamFileDict(is2);
            return foamFileDict.getOrDefault<word>("class", word::null);
        }
        // Stop after reasonable header length
        if (is2.lineNumber() > 50)
        {
            break;
        }
    }

    return word::null;
}


// Compare scalar fields and return statistics
void compareScalarFields
(
    const scalarField& sourceValues,
    const scalarField& subsetValues,
    const word& fieldName,
    scalar& maxAbsDiff,
    scalar& meanAbsDiff,
    scalar& rmsDiff,
    scalar& maxRelErr,
    scalar& meanRelErr
)
{
    if (sourceValues.size() != subsetValues.size())
    {
        FatalErrorInFunction
            << "Field " << fieldName << " size mismatch: "
            << "source=" << sourceValues.size()
            << " subset=" << subsetValues.size()
            << exit(FatalError);
    }

    scalarField diff = subsetValues - sourceValues;
    scalarField absDiff = mag(diff);

    maxAbsDiff = max(absDiff);
    meanAbsDiff = average(absDiff);
    rmsDiff = Foam::sqrt(average(sqr(diff)));

    // Relative error (avoid division by zero)
    scalarField relErr(sourceValues.size(), 0.0);
    forAll(sourceValues, i)
    {
        if (mag(sourceValues[i]) > SMALL)
        {
            relErr[i] = absDiff[i] / mag(sourceValues[i]);
        }
    }

    maxRelErr = max(relErr);
    meanRelErr = average(relErr);
}


// Compare vector fields and return statistics
void compareVectorFields
(
    const vectorField& sourceValues,
    const vectorField& subsetValues,
    const word& fieldName,
    scalar& maxAbsDiff,
    scalar& meanAbsDiff,
    scalar& rmsDiff,
    scalar& maxRelErr,
    scalar& meanRelErr
)
{
    if (sourceValues.size() != subsetValues.size())
    {
        FatalErrorInFunction
            << "Field " << fieldName << " size mismatch: "
            << "source=" << sourceValues.size()
            << " subset=" << subsetValues.size()
            << exit(FatalError);
    }

    vectorField diff = subsetValues - sourceValues;
    scalarField diffMag = mag(diff);
    scalarField sourceMag = mag(sourceValues);

    maxAbsDiff = max(diffMag);
    meanAbsDiff = average(diffMag);
    rmsDiff = Foam::sqrt(average(magSqr(diff)));

    // Relative error
    scalarField relErr(sourceValues.size(), 0.0);
    forAll(sourceValues, i)
    {
        if (sourceMag[i] > SMALL)
        {
            relErr[i] = diffMag[i] / sourceMag[i];
        }
    }

    maxRelErr = max(relErr);
    meanRelErr = average(relErr);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Compare internal fields between original source case and "
        "recalculated subset case"
    );

    timeSelector::addOptions();

    argList::addOption
    (
        "sourceCase",
        "dir",
        "Path to original source case (required)"
    );

    argList::addOption
    (
        "fields",
        "list",
        "Fields to compare (default: \"(U p nut)\")"
    );

    #include "setRootCase.H"

    // Get source case path
    fileName sourceCasePath;
    if (!args.readIfPresent("sourceCase", sourceCasePath))
    {
        FatalErrorInFunction
            << "Missing required option -sourceCase" << nl
            << exit(FatalError);
    }

    // Make path absolute if relative
    if (!sourceCasePath.isAbsolute())
    {
        sourceCasePath = cwd() / sourceCasePath;
    }

    // Get fields to compare
    wordList fieldNames({"U", "p", "nut"});
    args.readIfPresent("fields", fieldNames);

    Info<< "Comparing fields: " << fieldNames << nl << endl;

    // Create time objects for both cases
    #include "createTime.H"

    // Create subset mesh (current case)
    #include "createMesh.H"

    Info<< "Subset case: " << runTime.path() << nl
        << "Subset mesh cells: " << mesh.nCells() << nl << endl;

    // Create source time and mesh
    Time sourceRunTime
    (
        Time::controlDictName,
        sourceCasePath.path(),
        sourceCasePath.name()
    );

    fvMesh sourceMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            sourceRunTime.timeName(),
            sourceRunTime,
            IOobject::MUST_READ
        )
    );

    Info<< "Source case: " << sourceRunTime.path() << nl
        << "Source mesh cells: " << sourceMesh.nCells() << nl << endl;

    // Read extraction box from subset case
    fileName extractionBoxPath =
        mesh.time().path() / "constant" / "polyMesh" / "extractionBox";

    point boxMin, boxMax;

    if (isFile(extractionBoxPath))
    {
        IFstream is(extractionBoxPath);
        dictionary extractionDict(is);

        boxMin = extractionDict.get<point>("boxMin");
        boxMax = extractionDict.get<point>("boxMax");

        Info<< "Extraction box from file:" << nl
            << "    min: " << boxMin << nl
            << "    max: " << boxMax << nl << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find extractionBox file at " << extractionBoxPath << nl
            << "This file is required to identify cells in source mesh"
            << exit(FatalError);
    }

    boundBox bb(boxMin, boxMax);

    // Find cells in source mesh that are inside the extraction box
    Info<< "Finding cells in source mesh inside extraction box..." << endl;
    labelList sourceCellsInBox = findCellsInBox(sourceMesh, bb);
    Info<< "Found " << sourceCellsInBox.size() << " cells in box" << nl << endl;

    if (sourceCellsInBox.size() != mesh.nCells())
    {
        WarningInFunction
            << "Cell count mismatch!" << nl
            << "    Source cells in box: " << sourceCellsInBox.size() << nl
            << "    Subset mesh cells:   " << mesh.nCells() << nl
            << "    Results may not be valid." << endl;
    }

    // Select times
    instantList subsetTimes = timeSelector::select0(runTime, args);

    Info<< "======================================================" << nl
        << "FIELD COMPARISON RESULTS" << nl
        << "======================================================" << nl
        << endl;

    // Storage for summary statistics
    HashTable<scalarList> summaryMaxAbs;
    HashTable<scalarList> summaryMeanRel;

    for (const word& fieldName : fieldNames)
    {
        summaryMaxAbs.insert(fieldName, scalarList());
        summaryMeanRel.insert(fieldName, scalarList());
    }

    // Compare at each timestep
    forAll(subsetTimes, timeI)
    {
        runTime.setTime(subsetTimes[timeI], timeI);

        // Find corresponding source time
        word timeName = runTime.timeName();

        // Check if source time exists
        fileName sourceTimeDir = sourceCasePath / timeName;
        if (!isDir(sourceTimeDir))
        {
            Info<< "Skipping time " << timeName
                << " (not found in source case)" << nl;
            continue;
        }

        sourceRunTime.setTime(subsetTimes[timeI], timeI);

        Info<< "Time: " << timeName << nl
            << "------------------------------------------------------" << endl;

        for (const word& fieldName : fieldNames)
        {
            // Check if files exist
            fileName subsetFieldPath = runTime.path() / timeName / fieldName;
            fileName sourceFieldPath = sourceCasePath / timeName / fieldName;

            if (!isFile(subsetFieldPath))
            {
                Info<< "  " << fieldName << ": not found in subset case" << nl;
                continue;
            }

            if (!isFile(sourceFieldPath))
            {
                Info<< "  " << fieldName << ": not found in source case" << nl;
                continue;
            }

            // Determine field type
            word fieldClass = getFieldType(subsetFieldPath);

            if (fieldClass == "volVectorField")
            {
                // Read internal fields only
                tmp<vectorField> tsubsetInternal = readInternalField<vector>
                (
                    subsetFieldPath,
                    mesh.nCells()
                );

                tmp<vectorField> tsourceInternal = readInternalField<vector>
                (
                    sourceFieldPath,
                    sourceMesh.nCells()
                );

                // Extract source values at cells in box
                vectorField sourceValues(sourceCellsInBox.size());
                forAll(sourceCellsInBox, i)
                {
                    sourceValues[i] = tsourceInternal()[sourceCellsInBox[i]];
                }

                // Compare
                scalar maxAbsDiff, meanAbsDiff, rmsDiff, maxRelErr, meanRelErr;
                compareVectorFields
                (
                    sourceValues,
                    tsubsetInternal(),
                    fieldName,
                    maxAbsDiff, meanAbsDiff, rmsDiff, maxRelErr, meanRelErr
                );

                Info<< "  " << fieldName << " (vector):" << nl
                    << "    Max abs diff:   " << maxAbsDiff << nl
                    << "    Mean abs diff:  " << meanAbsDiff << nl
                    << "    RMS diff:       " << rmsDiff << nl
                    << "    Max rel error:  " << maxRelErr << nl
                    << "    Mean rel error: " << meanRelErr << nl;

                summaryMaxAbs[fieldName].append(maxAbsDiff);
                summaryMeanRel[fieldName].append(meanRelErr);
            }
            else if (fieldClass == "volScalarField")
            {
                // Read internal fields only
                tmp<scalarField> tsubsetInternal = readInternalField<scalar>
                (
                    subsetFieldPath,
                    mesh.nCells()
                );

                tmp<scalarField> tsourceInternal = readInternalField<scalar>
                (
                    sourceFieldPath,
                    sourceMesh.nCells()
                );

                // Extract source values at cells in box
                scalarField sourceValues(sourceCellsInBox.size());
                forAll(sourceCellsInBox, i)
                {
                    sourceValues[i] = tsourceInternal()[sourceCellsInBox[i]];
                }

                // Compare
                scalar maxAbsDiff, meanAbsDiff, rmsDiff, maxRelErr, meanRelErr;
                compareScalarFields
                (
                    sourceValues,
                    tsubsetInternal(),
                    fieldName,
                    maxAbsDiff, meanAbsDiff, rmsDiff, maxRelErr, meanRelErr
                );

                Info<< "  " << fieldName << " (scalar):" << nl
                    << "    Max abs diff:   " << maxAbsDiff << nl
                    << "    Mean abs diff:  " << meanAbsDiff << nl
                    << "    RMS diff:       " << rmsDiff << nl
                    << "    Max rel error:  " << maxRelErr << nl
                    << "    Mean rel error: " << meanRelErr << nl;

                summaryMaxAbs[fieldName].append(maxAbsDiff);
                summaryMeanRel[fieldName].append(meanRelErr);
            }
            else
            {
                Info<< "  " << fieldName << ": unknown field type ("
                    << fieldClass << ")" << nl;
            }
        }

        Info<< endl;
    }

    // Print summary
    Info<< "======================================================" << nl
        << "SUMMARY (across all timesteps)" << nl
        << "======================================================" << nl;

    for (const word& fieldName : fieldNames)
    {
        const scalarList& maxAbsList = summaryMaxAbs[fieldName];
        const scalarList& meanRelList = summaryMeanRel[fieldName];

        if (!maxAbsList.empty())
        {
            Info<< fieldName << ":" << nl
                << "    Max abs diff (worst):     " << max(maxAbsList) << nl
                << "    Avg mean rel error:       " << average(meanRelList) << nl;
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
