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

#include "spaceTimeWindowExtract.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMeshSubset.H"
#include "OFstream.H"
#include "Time.H"
#include "surfaceFields.H"
#include "Pstream.H"
#include "foamVersion.H"
#include "argList.H"
#include "boundBox.H"
#include "polyPatch.H"
#include "Tuple2.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(spaceTimeWindowExtract, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        spaceTimeWindowExtract,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::spaceTimeWindowExtract::writeFoamFileHeader
(
    Ostream& os,
    const word& className,
    const word& objectName
) const
{
    os  << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       " << className << ";" << nl
        << "    object      " << objectName << ";" << nl
        << "}" << nl << nl;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::spaceTimeWindowExtract::initializeSubset()
{
    if (meshSubsetPtr_)
    {
        return; // Already initialized
    }

    Info<< type() << " " << name() << ": Initializing mesh subset" << endl;

    // Create bounding box from min/max points
    boundBox bb(boxMin_, boxMax_);

    Info<< "    Extraction box: " << bb.min() << " to " << bb.max() << endl;

    // Find cells whose centres are inside the box
    const vectorField& cellCentres = mesh_.cellCentres();
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
            << "No cells found inside box " << bb << nl
            << "Check that the box is inside the mesh domain" << nl
            << exit(FatalError);
    }

    // Create mesh subset
    meshSubsetPtr_.reset(new fvMeshSubset(mesh_));

    // Use setCellSubset with cells inside the box
    meshSubsetPtr_->setCellSubset(cellsInBox);

    const fvMesh& subMesh = meshSubsetPtr_->subMesh();

    Info<< "    Subset mesh has:" << nl
        << "        " << subMesh.nCells() << " cells" << nl
        << "        " << subMesh.nFaces() << " faces" << nl
        << "        " << subMesh.nPoints() << " points" << endl;

    // Find the oldInternalFaces patch
    label oldInternalPatchId = -1;
    const polyBoundaryMesh& pbm = subMesh.boundaryMesh();

    forAll(pbm, patchi)
    {
        if (pbm[patchi].name() == "oldInternalFaces")
        {
            oldInternalPatchId = patchi;
            break;
        }
    }

    if (oldInternalPatchId < 0)
    {
        FatalErrorInFunction
            << "Cannot find patch 'oldInternalFaces' in subset mesh" << nl
            << "    Available patches: " << pbm.names() << nl
            << exit(FatalError);
    }

    const polyPatch& oldInternalPatch = pbm[oldInternalPatchId];

    Info<< "    oldInternalFaces patch has "
        << oldInternalPatch.size() << " faces" << endl;

    // Get face map: subset face -> original face
    const labelList& faceMap = meshSubsetPtr_->faceMap();

    // Get cell map: subset cell -> original cell
    const labelList& cellMap = meshSubsetPtr_->cellMap();

    // Build lookup: original cell -> subset cell
    labelList origCellToSubset(mesh_.nCells(), -1);
    forAll(cellMap, subCelli)
    {
        origCellToSubset[cellMap[subCelli]] = subCelli;
    }

    // Store data for each oldInternalFace
    const label nFaces = oldInternalPatch.size();
    oldInternalFaceIndices_.setSize(nFaces);
    faceCellsInside_.setSize(nFaces);
    faceCellsOutside_.setSize(nFaces);
    faceWeights_.setSize(nFaces);
    faceCentres_.setSize(nFaces);

    const vectorField& faceCentresAll = mesh_.faceCentres();

    forAll(oldInternalPatch, facei)
    {
        // Global face index in subset mesh
        label subFaceI = oldInternalPatch.start() + facei;

        // Original face index in parent mesh
        label origFaceI = faceMap[subFaceI];
        oldInternalFaceIndices_[facei] = origFaceI;

        // Store face centre
        faceCentres_[facei] = faceCentresAll[origFaceI];

        // Get owner/neighbour cells in original mesh
        label origOwner = mesh_.faceOwner()[origFaceI];
        label origNeighbour = mesh_.faceNeighbour()[origFaceI];

        // Determine which cell is inside the subset
        label subOwner = origCellToSubset[origOwner];
        label subNeighbour = origCellToSubset[origNeighbour];

        label cellInside, cellOutside;

        if (subOwner >= 0 && subNeighbour < 0)
        {
            cellInside = origOwner;
            cellOutside = origNeighbour;
        }
        else if (subNeighbour >= 0 && subOwner < 0)
        {
            cellInside = origNeighbour;
            cellOutside = origOwner;
        }
        else
        {
            FatalErrorInFunction
                << "oldInternalFace " << facei << " (original " << origFaceI
                << ") has unexpected owner/neighbour mapping" << nl
                << "    origOwner=" << origOwner << " -> subset=" << subOwner << nl
                << "    origNeighbour=" << origNeighbour << " -> subset=" << subNeighbour
                << exit(FatalError);
        }

        faceCellsInside_[facei] = cellInside;
        faceCellsOutside_[facei] = cellOutside;

        // Compute interpolation weight
        // Face value = w * inside + (1-w) * outside
        // where w = distance(outside to face) / distance(outside to inside)
        const point& Cf = faceCentres_[facei];
        const point& Ci = cellCentres[cellInside];
        const point& Co = cellCentres[cellOutside];

        scalar d_inside = mag(Cf - Ci);
        scalar d_outside = mag(Cf - Co);
        scalar d_total = d_inside + d_outside;

        // Weight for inside cell
        faceWeights_[facei] = d_outside / (d_total + SMALL);
    }

    Info<< "    Computed interpolation weights for "
        << nFaces << " boundary faces" << endl;
}


void Foam::functionObjects::spaceTimeWindowExtract::writeSubsetMesh()
{
    if (meshWritten_)
    {
        return;
    }

    const fvMesh& subMesh = meshSubsetPtr_->subMesh();

    // Create output directory
    fileName meshDir = outputDir_ / "constant" / "polyMesh";
    mkDir(meshDir);

    Info<< type() << " " << name() << ": Writing subset mesh to "
        << meshDir << endl;

    // Write mesh in ASCII format
    const polyMesh& pm = subMesh;

    // Points
    {
        OFstream os(meshDir / "points");
        writeFoamFileHeader(os, "vectorField", "points");
        os << pm.points();
    }

    // Faces
    {
        OFstream os(meshDir / "faces");
        writeFoamFileHeader(os, "faceList", "faces");
        os << pm.faces();
    }

    // Owner
    {
        OFstream os(meshDir / "owner");
        writeFoamFileHeader(os, "labelList", "owner");
        os << pm.faceOwner();
    }

    // Neighbour
    {
        OFstream os(meshDir / "neighbour");
        writeFoamFileHeader(os, "labelList", "neighbour");
        os << pm.faceNeighbour();
    }

    // Boundary - write manually to ensure oldInternalFaces is type 'patch'
    {
        OFstream os(meshDir / "boundary");
        writeFoamFileHeader(os, "polyBoundaryMesh", "boundary");

        const polyBoundaryMesh& bm = pm.boundaryMesh();

        // Count only patches with faces (for a fully internal box, only oldInternalFaces)
        label nPatches = 0;
        forAll(bm, patchi)
        {
            if (bm[patchi].size() > 0 || bm[patchi].name() == "oldInternalFaces")
            {
                nPatches++;
            }
        }

        // For a fully internal subset, we should only have oldInternalFaces
        // Write just that patch with type 'patch'
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

        Info<< "    Written boundary with oldInternalFaces as type 'patch'" << endl;
    }

    // Write face centres for reference
    fileName boundaryDir = outputDir_ / "constant" / "boundaryData" / "oldInternalFaces";
    mkDir(boundaryDir);

    {
        OFstream os(boundaryDir / "points");
        writeFoamFileHeader(os, "vectorField", "points");
        os << faceCentres_;
    }

    // Write extraction metadata for validation during reconstruction
    {
        OFstream os(boundaryDir / "extractionMetadata");
        writeFoamFileHeader(os, "dictionary", "extractionMetadata");

        const Time& runTime = mesh_.time();

        os << "// Extraction metadata for spaceTimeWindow reconstruction" << nl
           << "// Reconstruction MUST use identical settings" << nl << nl;

        // OpenFOAM version information
        os.writeEntry("openfoamVersion", word(foamVersion::version));
        os.writeEntry("openfoamApi", foamVersion::api);

        // Solver/application name
        const word solverName = runTime.controlDict().getOrDefault<word>
        (
            "application",
            word("unknown")
        );
        os.writeEntry("solver", solverName);

        os << nl << "// Time settings" << nl;

        os.writeEntry("deltaT", runTime.deltaTValue());
        os.writeEntry("adjustTimeStep", runTime.controlDict().getOrDefault<Switch>("adjustTimeStep", false));

        // Write time precision
        const unsigned int timePrecision = IOstream::defaultPrecision();
        os.writeEntry("timePrecision", timePrecision);

        // Write startTime for reference
        os.writeEntry("extractionStartTime", runTime.value());

        // Check if using fixed time step
        if (runTime.controlDict().found("deltaT"))
        {
            os.writeEntry("fixedDeltaT", runTime.controlDict().get<scalar>("deltaT"));
        }

        Info<< "    Written extraction metadata:" << nl
            << "        openfoamVersion = " << foamVersion::version << nl
            << "        openfoamApi = " << foamVersion::api << nl
            << "        solver = " << solverName << nl
            << "        deltaT = " << runTime.deltaTValue() << nl
            << "        adjustTimeStep = "
            << runTime.controlDict().getOrDefault<Switch>("adjustTimeStep", false) << nl
            << "        timePrecision = " << timePrecision << endl;
    }

    meshWritten_ = true;
}


void Foam::functionObjects::spaceTimeWindowExtract::writeInitialFields()
{
    const fvMesh& subMesh = meshSubsetPtr_->subMesh();

    // Create time directory
    const word timeName = mesh_.time().timeName();
    fileName timeDir = outputDir_ / timeName;
    mkDir(timeDir);

    Info<< type() << " " << name() << ": Writing initial fields to "
        << timeDir << endl;

    for (const word& fieldName : fieldNames_)
    {
        // Try scalar field
        const auto* sfPtr = mesh_.findObject<volScalarField>(fieldName);
        if (sfPtr)
        {
            // Interpolate to subset mesh using fvMeshSubset
            tmp<volScalarField> tsubField = meshSubsetPtr_->interpolate(*sfPtr);
            const volScalarField& subField = tsubField();

            // Write subset field
            OFstream os(timeDir / fieldName);
            writeFoamFileHeader(os, "volScalarField", fieldName);

            os  << "dimensions      " << sfPtr->dimensions() << ";" << nl << nl
                << "internalField   ";

            // Write internal field
            os << "nonuniform List<scalar>" << nl
               << subField.primitiveField().size() << nl
               << '(' << nl;
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

            continue;
        }

        // Try vector field
        const auto* vfPtr = mesh_.findObject<volVectorField>(fieldName);
        if (vfPtr)
        {
            tmp<volVectorField> tsubField = meshSubsetPtr_->interpolate(*vfPtr);
            const volVectorField& subField = tsubField();

            OFstream os(timeDir / fieldName);
            writeFoamFileHeader(os, "volVectorField", fieldName);

            os  << "dimensions      " << vfPtr->dimensions() << ";" << nl << nl
                << "internalField   ";

            os << "nonuniform List<vector>" << nl
               << subField.primitiveField().size() << nl
               << '(' << nl;
            forAll(subField.primitiveField(), i)
            {
                os << subField.primitiveField()[i] << nl;
            }
            os << ')' << ';' << nl << nl;

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
        }
    }
}


void Foam::functionObjects::spaceTimeWindowExtract::writeBoundaryData()
{
    const word timeName = mesh_.time().timeName();
    fileName boundaryTimeDir =
        outputDir_ / "constant" / "boundaryData" / "oldInternalFaces" / timeName;
    mkDir(boundaryTimeDir);

    // Track this timestep (exact string name)
    extractedTimesteps_.append(timeName);

    DebugInfo
        << type() << " " << name() << ": Writing boundary data to "
        << boundaryTimeDir << endl;

    for (const word& fieldName : fieldNames_)
    {
        // Try scalar field
        const auto* sfPtr = mesh_.findObject<volScalarField>(fieldName);
        if (sfPtr)
        {
            tmp<scalarField> tfaceValues = interpolateToFaces(*sfPtr);
            writeField(boundaryTimeDir, fieldName, tfaceValues());
            continue;
        }

        // Try vector field
        const auto* vfPtr = mesh_.findObject<volVectorField>(fieldName);
        if (vfPtr)
        {
            tmp<vectorField> tfaceValues = interpolateToFaces(*vfPtr);
            writeField(boundaryTimeDir, fieldName, tfaceValues());
        }
    }
}


void Foam::functionObjects::spaceTimeWindowExtract::updateMetadataTimesteps()
{
    // Append timestep list to the metadata file
    fileName metadataPath = outputDir_ / "constant" / "boundaryData"
        / "oldInternalFaces" / "extractionMetadata";

    // Append to existing file
    OFstream os(metadataPath, IOstreamOption(), IOstreamOption::APPEND);

    os << nl << "// Extracted timesteps (precise directory names)" << nl;
    os << "timesteps" << nl << "(" << nl;
    for (const word& ts : extractedTimesteps_)
    {
        os << "    " << ts << nl;
    }
    os << ");" << nl;

    Info<< type() << " " << name() << ": Written " << extractedTimesteps_.size()
        << " timesteps to metadata" << endl;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::spaceTimeWindowExtract::interpolateToFaces
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    auto tfaceValues = tmp<Field<Type>>::New(faceCellsInside_.size());
    Field<Type>& faceValues = tfaceValues.ref();

    forAll(faceValues, facei)
    {
        const Type& valueInside = vf[faceCellsInside_[facei]];
        const Type& valueOutside = vf[faceCellsOutside_[facei]];
        const scalar w = faceWeights_[facei];

        faceValues[facei] = w * valueInside + (1.0 - w) * valueOutside;
    }

    return tfaceValues;
}


template<class Type>
void Foam::functionObjects::spaceTimeWindowExtract::writeField
(
    const fileName& dir,
    const word& fieldName,
    const Field<Type>& field
) const
{
    OFstream os(dir / fieldName);

    // Simple field format - just the data
    os << field;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::spaceTimeWindowExtract::~spaceTimeWindowExtract()
{
    // Write final metadata with timestep list when the function object is destroyed
    if (meshWritten_ && !extractedTimesteps_.empty())
    {
        updateMetadataTimesteps();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::spaceTimeWindowExtract::spaceTimeWindowExtract
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    boxMin_(Zero),
    boxMax_(Zero),
    outputDir_(),
    fieldNames_(),
    meshSubsetPtr_(),
    meshWritten_(false),
    oldInternalFaceIndices_(),
    faceCellsInside_(),
    faceCellsOutside_(),
    faceWeights_(),
    faceCentres_(),
    extractedTimesteps_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::spaceTimeWindowExtract::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Check for parallel - not supported yet
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "spaceTimeWindowExtract does not support parallel runs." << nl
            << "Please run in serial mode, or use reconstructPar first" << nl
            << "and then run extraction on the reconstructed case." << nl
            << exit(FatalError);
    }

    // Read box definition (two points: min and max)
    Tuple2<point, point> boxDef = dict.get<Tuple2<point, point>>("box");
    boxMin_ = boxDef.first();
    boxMax_ = boxDef.second();

    dict.readEntry("outputDir", outputDir_);
    dict.readEntry("fields", fieldNames_);

    // Make outputDir absolute if relative
    if (!outputDir_.isAbsolute())
    {
        outputDir_ = mesh_.time().globalPath() / outputDir_;
    }

    Info<< type() << " " << name() << ":" << nl
        << "    box:       " << boxMin_ << " " << boxMax_ << nl
        << "    outputDir: " << outputDir_ << nl
        << "    fields:    " << fieldNames_ << endl;

    return true;
}


bool Foam::functionObjects::spaceTimeWindowExtract::execute()
{
    // Ensure subset is initialized
    initializeSubset();

    // Write mesh once (on first execute)
    if (!meshWritten_)
    {
        writeSubsetMesh();
        writeInitialFields();
    }

    // Write boundary data at EVERY timestep (not just write intervals)
    // This is critical - the spaceTimeWindow BC needs data at every timestep
    writeBoundaryData();

    return true;
}


bool Foam::functionObjects::spaceTimeWindowExtract::write()
{
    // All work is done in execute() to ensure every timestep is captured
    // write() is only called at writeInterval, which would skip timesteps
    return true;
}


// Explicit instantiation of template functions
template Foam::tmp<Foam::Field<Foam::scalar>>
Foam::functionObjects::spaceTimeWindowExtract::interpolateToFaces<Foam::scalar>
(
    const GeometricField<scalar, fvPatchField, volMesh>&
) const;

template Foam::tmp<Foam::Field<Foam::vector>>
Foam::functionObjects::spaceTimeWindowExtract::interpolateToFaces<Foam::vector>
(
    const GeometricField<vector, fvPatchField, volMesh>&
) const;

template void Foam::functionObjects::spaceTimeWindowExtract::writeField<Foam::scalar>
(
    const fileName&,
    const word&,
    const Field<scalar>&
) const;

template void Foam::functionObjects::spaceTimeWindowExtract::writeField<Foam::vector>
(
    const fileName&,
    const word&,
    const Field<vector>&
) const;


// ************************************************************************* //
