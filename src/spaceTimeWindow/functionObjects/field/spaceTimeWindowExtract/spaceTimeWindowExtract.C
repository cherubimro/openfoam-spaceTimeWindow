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
#include "PstreamBuffers.H"
#include "foamVersion.H"
#include "argList.H"
#include "boundBox.H"
#include "polyPatch.H"
#include "processorPolyPatch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "ListOps.H"
#include "deltaVarintCodec.H"
#include "sodiumCrypto.H"
#include <fstream>
#include <type_traits>
#include <cstring>

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
    if (subsetInitialized_)
    {
        return; // Already initialized
    }

    Info<< type() << " " << name() << ": Initializing mesh subset" << endl;

    // Create bounding box from min/max points
    boundBox bb(boxMin_, boxMax_);

    Info<< "    Extraction box: " << bb.min() << " to " << bb.max() << endl;

    // Find cells whose centres are inside the box (local to this processor)
    const vectorField& cellCentres = mesh_.cellCentres();
    labelList cellsInBox;

    forAll(cellCentres, celli)
    {
        if (bb.contains(cellCentres[celli]))
        {
            cellsInBox.append(celli);
        }
    }

    // In parallel, check total cells across all processors
    label localCells = cellsInBox.size();
    label totalCells = localCells;
    reduce(totalCells, sumOp<label>());

    Info<< "    Found " << localCells << " cells inside box on this processor"
        << " (" << totalCells << " total)" << endl;


    if (totalCells == 0)
    {
        FatalErrorInFunction
            << "No cells found inside box " << bb << nl
            << "Check that the box is inside the mesh domain" << nl
            << exit(FatalError);
    }

    // Build set of cells inside box for fast lookup
    labelHashSet cellsInBoxSet(cellsInBox);

    // Identify boundary faces of the extraction region:
    // 1. Internal faces where one cell is inside, one is outside
    // 2. Processor boundary faces where local cell is inside, remote is outside
    //    (or local is outside, remote is inside - need both for proper coverage)

    DynamicList<label> boundaryFaceIndices;
    DynamicList<label> boundaryCellInside;
    DynamicList<label> boundaryCellOutside;
    DynamicList<vector> boundaryFaceCentres;
    DynamicList<scalar> boundaryWeights;

    // Also track faces where outside cell is on another processor
    DynamicList<label> procBoundaryLocalFaceIdx;  // index in our boundary list
    DynamicList<label> procBoundaryRemoteCell;    // remote cell index
    DynamicList<label> procBoundaryRemoteProc;    // which processor

    const vectorField& faceCentresAll = mesh_.faceCentres();

    // 1. Check internal faces
    for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
    {
        label owner = mesh_.faceOwner()[facei];
        label neighbour = mesh_.faceNeighbour()[facei];

        bool ownerInside = cellsInBoxSet.found(owner);
        bool neighbourInside = cellsInBoxSet.found(neighbour);

        if (ownerInside != neighbourInside)
        {
            // This is a boundary face of the extraction region
            label cellInside = ownerInside ? owner : neighbour;
            label cellOutside = ownerInside ? neighbour : owner;

            boundaryFaceIndices.append(facei);
            boundaryCellInside.append(cellInside);
            boundaryCellOutside.append(cellOutside);
            boundaryFaceCentres.append(faceCentresAll[facei]);

            // Compute weight
            const point& Cf = faceCentresAll[facei];
            const point& Ci = cellCentres[cellInside];
            const point& Co = cellCentres[cellOutside];

            scalar d_inside = mag(Cf - Ci);
            scalar d_outside = mag(Cf - Co);
            scalar d_total = d_inside + d_outside;

            boundaryWeights.append(d_outside / (d_total + SMALL));
        }
    }

    // 2. Check processor boundary faces
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            label nbrProc = procPatch.neighbProcNo();

            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label owner = mesh_.faceOwner()[facei];

                bool ownerInside = cellsInBoxSet.found(owner);

                // For processor patches, we need to know if the remote cell
                // is inside or outside. We'll exchange this info.
                // For now, mark faces where local cell is inside
                // (the remote cell data will be fetched during interpolation)

                if (ownerInside)
                {
                    // Local cell is inside, remote might be outside
                    // We'll store this and verify during comm exchange
                    boundaryFaceIndices.append(facei);
                    boundaryCellInside.append(owner);
                    boundaryCellOutside.append(-1);  // Mark as remote
                    boundaryFaceCentres.append(faceCentresAll[facei]);

                    // Weight will be computed after getting remote cell centre
                    boundaryWeights.append(-1.0);  // Placeholder

                    // Track for parallel communication
                    procBoundaryLocalFaceIdx.append(boundaryFaceIndices.size() - 1);
                    procBoundaryRemoteCell.append(i);  // Index in processor patch
                    procBoundaryRemoteProc.append(nbrProc);
                }
            }
        }
    }

    // Exchange information about which cells are inside the box
    // across processor boundaries to identify true boundary faces
    if (Pstream::parRun())
    {
        initializeParallelComm();
        identifyProcessorBoundaryFaces(bb, cellsInBox);
    }

    // Transfer to member variables
    oldInternalFaceIndices_.transfer(boundaryFaceIndices);
    faceCellsInside_.transfer(boundaryCellInside);
    faceCellsOutside_.transfer(boundaryCellOutside);
    faceCentres_.transfer(boundaryFaceCentres);
    faceWeights_.transfer(boundaryWeights);

    remoteFaceIndices_.transfer(procBoundaryLocalFaceIdx);
    remoteCellIndices_.transfer(procBoundaryRemoteCell);
    remoteCellProcs_.transfer(procBoundaryRemoteProc);

    hasProcessorBoundaryFaces_ = !remoteFaceIndices_.empty();
    reduce(hasProcessorBoundaryFaces_, orOp<bool>());

    // Create global index for face numbering
    globalFaceIndexPtr_.reset(new globalIndex(oldInternalFaceIndices_.size()));
    nGlobalFaces_ = globalFaceIndexPtr_->totalSize();

    // Now we need to actually determine which processor boundary faces
    // are true extraction boundaries (remote cell is outside box)
    if (Pstream::parRun() && hasProcessorBoundaryFaces_)
    {
        // Exchange cell-in-box status across processor boundaries
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // For each processor patch, send status of our boundary cells
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Send: for each face, is our cell inside the box?
                boolList localInside(pp.size());
                forAll(pp, i)
                {
                    label owner = mesh_.faceOwner()[pp.start() + i];
                    localInside[i] = cellsInBoxSet.found(owner);
                }

                UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
                toNbr << localInside;
            }
        }

        pBufs.finishedSends();

        // Receive and process - store results for later filtering
        Map<bool> procFaceIsBoundary;  // Store boundary status for filtering
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                boolList remoteInside(pp.size());
                UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                fromNbr >> remoteInside;

                // Store boundary status for each processor face
                forAll(pp, i)
                {
                    label facei = pp.start() + i;
                    label owner = mesh_.faceOwner()[facei];
                    bool localIn = cellsInBoxSet.found(owner);
                    bool remoteIn = remoteInside[i];

                    // True boundary if exactly one is inside
                    procFaceIsBoundary.insert(facei, localIn != remoteIn);
                }
            }
        }

        // Exchange cell centres for weight computation
        PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Send cell centres for our boundary cells
                vectorField localCentres(pp.size());
                forAll(pp, i)
                {
                    label owner = mesh_.faceOwner()[pp.start() + i];
                    localCentres[i] = cellCentres[owner];
                }

                UOPstream toNbr(procPatch.neighbProcNo(), pBufs2);
                toNbr << localCentres;
            }
        }

        pBufs2.finishedSends();

        // Receive remote cell centres and compute weights
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                vectorField remoteCentres(pp.size());
                UIPstream fromNbr(procPatch.neighbProcNo(), pBufs2);
                fromNbr >> remoteCentres;

                // Store remote centres for faces we own
                forAll(remoteFaceIndices_, idx)
                {
                    if (remoteCellProcs_[idx] == procPatch.neighbProcNo())
                    {
                        label localFaceIdx = remoteFaceIndices_[idx];
                        label patchFaceIdx = remoteCellIndices_[idx];

                        // Compute weight using remote cell centre
                        const point& Cf = faceCentres_[localFaceIdx];
                        const point& Ci = cellCentres[faceCellsInside_[localFaceIdx]];
                        const point& Co = remoteCentres[patchFaceIdx];

                        scalar d_inside = mag(Cf - Ci);
                        scalar d_outside = mag(Cf - Co);
                        scalar d_total = d_inside + d_outside;

                        faceWeights_[localFaceIdx] = d_outside / (d_total + SMALL);

                        // Store remote cell index (local index on remote proc)
                        // This will be used for field value exchange
                        faceCellsOutside_[localFaceIdx] = patchFaceIdx;
                    }
                }
            }
        }


        // Filter out faces that are not true boundaries
        // (where both cells ended up being inside the box)
        // Use procFaceIsBoundary Map computed during pBufs exchange
        DynamicList<label> validFaceIndices;
        DynamicList<label> validCellInside;
        DynamicList<label> validCellOutside;
        DynamicList<vector> validFaceCentres;
        DynamicList<scalar> validWeights;
        DynamicList<label> validRemoteFaceIdx;
        DynamicList<label> validRemoteCellIdx;
        DynamicList<label> validRemoteProc;


        // Filter face lists
        forAll(oldInternalFaceIndices_, i)
        {
            label facei = oldInternalFaceIndices_[i];

            bool keep = true;
            if (facei >= mesh_.nInternalFaces())
            {
                // Processor boundary face - check if true boundary
                auto iter = procFaceIsBoundary.find(facei);
                if (iter.good())
                {
                    keep = iter.val();
                    // Only keep if local cell is inside (we own this face)
                    if (keep)
                    {
                        label owner = mesh_.faceOwner()[facei];
                        keep = cellsInBoxSet.found(owner);
                    }
                }
            }

            if (keep)
            {
                validFaceIndices.append(oldInternalFaceIndices_[i]);
                validCellInside.append(faceCellsInside_[i]);
                validCellOutside.append(faceCellsOutside_[i]);
                validFaceCentres.append(faceCentres_[i]);
                validWeights.append(faceWeights_[i]);
            }
        }

        // Update remote face tracking
        forAll(remoteFaceIndices_, i)
        {
            label oldIdx = remoteFaceIndices_[i];
            label facei = oldInternalFaceIndices_[oldIdx];

            auto iter = procFaceIsBoundary.find(facei);
            if (iter.good() && iter.val())
            {
                label owner = mesh_.faceOwner()[facei];
                if (cellsInBoxSet.found(owner))
                {
                    // Find new index
                    label newIdx = -1;
                    forAll(validFaceIndices, j)
                    {
                        if (validFaceIndices[j] == facei)
                        {
                            newIdx = j;
                            break;
                        }
                    }

                    if (newIdx >= 0)
                    {
                        validRemoteFaceIdx.append(newIdx);
                        validRemoteCellIdx.append(remoteCellIndices_[i]);
                        validRemoteProc.append(remoteCellProcs_[i]);
                    }
                }
            }
        }

        oldInternalFaceIndices_.transfer(validFaceIndices);
        faceCellsInside_.transfer(validCellInside);
        faceCellsOutside_.transfer(validCellOutside);
        faceCentres_.transfer(validFaceCentres);
        faceWeights_.transfer(validWeights);
        remoteFaceIndices_.transfer(validRemoteFaceIdx);
        remoteCellIndices_.transfer(validRemoteCellIdx);
        remoteCellProcs_.transfer(validRemoteProc);

        // Recreate global index with filtered count
        globalFaceIndexPtr_.reset(new globalIndex(oldInternalFaceIndices_.size()));
        nGlobalFaces_ = globalFaceIndexPtr_->totalSize();
    }

    // Create mesh subset for initial field output (only needs local cells)
    // NOTE: setCellSubset is a collective operation in parallel - all procs must call it
    // So we always create the meshSubset, even if empty, to avoid deadlock
    meshSubsetPtr_.reset(new fvMeshSubset(mesh_));
    if (!cellsInBox.empty())
    {
        meshSubsetPtr_->setCellSubset(cellsInBox);
    }
    else if (Pstream::parRun())
    {
        // For processors with no cells in box, create empty subset
        labelList emptyList;
        meshSubsetPtr_->setCellSubset(emptyList);
    }

    label localFaces = oldInternalFaceIndices_.size();
    reduce(localFaces, sumOp<label>());

    Info<< "    Extraction region boundary has " << localFaces
        << " faces globally" << endl;

    if (Pstream::parRun())
    {
        Info<< "    This processor has "
            << oldInternalFaceIndices_.size() << " boundary faces" << endl;

        if (hasProcessorBoundaryFaces_)
        {
            label nProcFaces = remoteFaceIndices_.size();
            reduce(nProcFaces, sumOp<label>());
            Info<< "    " << nProcFaces
                << " faces cross processor boundaries" << endl;
        }
    }

    subsetInitialized_ = true;
}


void Foam::functionObjects::spaceTimeWindowExtract::initializeParallelComm()
{
    // Initialize parallel communication structures
    // This is called once during setup
}


void Foam::functionObjects::spaceTimeWindowExtract::identifyProcessorBoundaryFaces
(
    const boundBox& bb,
    const labelList& cellsInBox
)
{
    // Additional processor boundary handling done in initializeSubset
    // This function reserved for more complex processing if needed
}


void Foam::functionObjects::spaceTimeWindowExtract::gatherFaceDataToMaster
(
    vectorField& allFaceCentres,
    scalarField& allFaceWeights,
    labelList& allFaceCellsInside,
    labelList& allFaceCellsOutside
) const
{
    if (!Pstream::parRun())
    {
        // Serial - just copy
        allFaceCentres = faceCentres_;
        allFaceWeights = faceWeights_;
        allFaceCellsInside = faceCellsInside_;
        allFaceCellsOutside = faceCellsOutside_;
        return;
    }

    // Gather all face centres to master
    List<vectorField> gatheredCentres(Pstream::nProcs());
    gatheredCentres[Pstream::myProcNo()] = faceCentres_;
    Pstream::gatherList(gatheredCentres);

    List<scalarField> gatheredWeights(Pstream::nProcs());
    gatheredWeights[Pstream::myProcNo()] = faceWeights_;
    Pstream::gatherList(gatheredWeights);

    List<labelList> gatheredCellsInside(Pstream::nProcs());
    gatheredCellsInside[Pstream::myProcNo()] = faceCellsInside_;
    Pstream::gatherList(gatheredCellsInside);

    List<labelList> gatheredCellsOutside(Pstream::nProcs());
    gatheredCellsOutside[Pstream::myProcNo()] = faceCellsOutside_;
    Pstream::gatherList(gatheredCellsOutside);

    if (Pstream::master())
    {
        // Combine all data
        label totalSize = 0;
        forAll(gatheredCentres, proci)
        {
            totalSize += gatheredCentres[proci].size();
        }

        allFaceCentres.setSize(totalSize);
        allFaceWeights.setSize(totalSize);
        allFaceCellsInside.setSize(totalSize);
        allFaceCellsOutside.setSize(totalSize);

        label offset = 0;
        forAll(gatheredCentres, proci)
        {
            forAll(gatheredCentres[proci], i)
            {
                allFaceCentres[offset + i] = gatheredCentres[proci][i];
                allFaceWeights[offset + i] = gatheredWeights[proci][i];
                allFaceCellsInside[offset + i] = gatheredCellsInside[proci][i];
                allFaceCellsOutside[offset + i] = gatheredCellsOutside[proci][i];
            }
            offset += gatheredCentres[proci].size();
        }
    }
}


void Foam::functionObjects::spaceTimeWindowExtract::gatherAndWriteSubsetMesh()
{
    // In parallel mode, we don't write the mesh directly because combining
    // distributed meshes with processor boundaries is complex (requires
    // proper point merging at processor interfaces).
    //
    // Instead, we write the extraction box parameters and let
    // spaceTimeWindowInitCase create the subset mesh from the source case.
    //
    // The initial fields are still gathered and written by the extractor.

    if (!Pstream::parRun())
    {
        // Serial - mesh is written directly in writeSubsetMesh()
        return;
    }

    // Only master writes the extraction parameters
    if (Pstream::master())
    {
        fileName meshDir = outputDir_ / "constant" / "polyMesh";
        mkDir(meshDir);

        // Write extraction box parameters for spaceTimeWindowInitCase
        {
            OFstream os(meshDir / "extractionBox");
            writeFoamFileHeader(os, "dictionary", "extractionBox");

            os << "// Extraction bounding box for mesh creation" << nl;
            os << "// spaceTimeWindowInitCase will create subset mesh from source case" << nl << nl;

            os.writeEntry("boxMin", boxMin_);
            os.writeEntry("boxMax", boxMax_);
            os.writeEntry("nGlobalFaces", nGlobalFaces_);
        }

        Info<< "    Written extraction parameters for parallel case" << nl
            << "    spaceTimeWindowInitCase will create mesh from source case" << endl;
    }
}


void Foam::functionObjects::spaceTimeWindowExtract::gatherAndWriteInitialFields()
{
    // Gather subset field values from all processors and write combined fields
    // Uses exact cell values from fvMeshSubset (no interpolation)

    if (!Pstream::parRun())
    {
        return;  // Serial case handled elsewhere
    }

    const word timeName = mesh_.time().timeName();

    // Gather cell counts for field assembly
    label localNCells = 0;
    label localNBoundaryFaces = 0;

    if (meshSubsetPtr_)
    {
        const fvMesh& subMesh = meshSubsetPtr_->subMesh();
        localNCells = subMesh.nCells();

        const polyBoundaryMesh& bm = subMesh.boundaryMesh();
        forAll(bm, patchi)
        {
            if (bm[patchi].name() == "oldInternalFaces")
            {
                localNBoundaryFaces = bm[patchi].size();
                break;
            }
        }
    }

    labelList gatheredNCells(Pstream::nProcs());
    gatheredNCells[Pstream::myProcNo()] = localNCells;
    Pstream::gatherList(gatheredNCells);

    labelList gatheredNBoundaryFaces(Pstream::nProcs());
    gatheredNBoundaryFaces[Pstream::myProcNo()] = localNBoundaryFaces;
    Pstream::gatherList(gatheredNBoundaryFaces);

    // Calculate totals on master
    label totalCells = 0;
    label totalBoundaryFaces = 0;
    if (Pstream::master())
    {
        forAll(gatheredNCells, proci)
        {
            totalCells += gatheredNCells[proci];
            totalBoundaryFaces += gatheredNBoundaryFaces[proci];
        }
    }

    // Create output directory on master
    fileName timeDir = outputDir_ / timeName;
    if (Pstream::master())
    {
        mkDir(timeDir);
    }
    UPstream::barrier(UPstream::worldComm);

    Info<< type() << " " << name() << ": Writing initial fields to "
        << timeDir << " (gathered from " << Pstream::nProcs() << " processors)" << endl;

    for (const word& fieldName : fieldNames_)
    {
        // Try scalar field
        const auto* sfPtr = mesh_.findObject<volScalarField>(fieldName);
        if (sfPtr)
        {
            // Get local subset field values (exact, no interpolation)
            scalarField localInternal;
            scalarField localBoundary;

            if (meshSubsetPtr_)
            {
                tmp<volScalarField> tsubField = meshSubsetPtr_->interpolate(*sfPtr);
                const volScalarField& subField = tsubField();
                localInternal = subField.primitiveField();

                // Get boundary values for oldInternalFaces
                const fvMesh& subMesh = meshSubsetPtr_->subMesh();
                forAll(subMesh.boundary(), patchi)
                {
                    if (subMesh.boundary()[patchi].name() == "oldInternalFaces")
                    {
                        localBoundary = subField.boundaryField()[patchi];
                        break;
                    }
                }
            }

            // Gather to master
            List<scalarField> gatheredInternal(Pstream::nProcs());
            gatheredInternal[Pstream::myProcNo()] = localInternal;
            Pstream::gatherList(gatheredInternal);

            List<scalarField> gatheredBoundary(Pstream::nProcs());
            gatheredBoundary[Pstream::myProcNo()] = localBoundary;
            Pstream::gatherList(gatheredBoundary);

            if (Pstream::master())
            {
                // Combine fields
                scalarField allInternal(totalCells);
                scalarField allBoundary(totalBoundaryFaces);

                label cellOffset = 0;
                label bfaceOffset = 0;
                forAll(gatheredInternal, proci)
                {
                    forAll(gatheredInternal[proci], i)
                    {
                        allInternal[cellOffset + i] = gatheredInternal[proci][i];
                    }
                    cellOffset += gatheredInternal[proci].size();

                    forAll(gatheredBoundary[proci], i)
                    {
                        allBoundary[bfaceOffset + i] = gatheredBoundary[proci][i];
                    }
                    bfaceOffset += gatheredBoundary[proci].size();
                }

                // Write field
                OFstream os(timeDir / fieldName);
                writeFoamFileHeader(os, "volScalarField", fieldName);

                os << "dimensions      " << sfPtr->dimensions() << ";" << nl << nl
                   << "internalField   nonuniform List<scalar>" << nl
                   << allInternal.size() << nl << '(' << nl;
                forAll(allInternal, i)
                {
                    os << allInternal[i] << nl;
                }
                os << ')' << ';' << nl << nl;

                os << "boundaryField" << nl << '{' << nl;
                os << "    oldInternalFaces" << nl
                   << "    {" << nl
                   << "        type            calculated;" << nl
                   << "        value           nonuniform List<scalar>" << nl
                   << "        " << allBoundary.size() << nl
                   << "        (" << nl;
                forAll(allBoundary, i)
                {
                    os << "        " << allBoundary[i] << nl;
                }
                os << "        );" << nl
                   << "    }" << nl;
                os << '}' << nl;
            }
            continue;
        }

        // Try vector field
        const auto* vfPtr = mesh_.findObject<volVectorField>(fieldName);
        if (vfPtr)
        {
            vectorField localInternal;
            vectorField localBoundary;

            if (meshSubsetPtr_)
            {
                tmp<volVectorField> tsubField = meshSubsetPtr_->interpolate(*vfPtr);
                const volVectorField& subField = tsubField();
                localInternal = subField.primitiveField();

                const fvMesh& subMesh = meshSubsetPtr_->subMesh();
                forAll(subMesh.boundary(), patchi)
                {
                    if (subMesh.boundary()[patchi].name() == "oldInternalFaces")
                    {
                        localBoundary = subField.boundaryField()[patchi];
                        break;
                    }
                }
            }

            List<vectorField> gatheredInternal(Pstream::nProcs());
            gatheredInternal[Pstream::myProcNo()] = localInternal;
            Pstream::gatherList(gatheredInternal);

            List<vectorField> gatheredBoundary(Pstream::nProcs());
            gatheredBoundary[Pstream::myProcNo()] = localBoundary;
            Pstream::gatherList(gatheredBoundary);

            if (Pstream::master())
            {
                vectorField allInternal(totalCells);
                vectorField allBoundary(totalBoundaryFaces);

                label cellOffset = 0;
                label bfaceOffset = 0;
                forAll(gatheredInternal, proci)
                {
                    forAll(gatheredInternal[proci], i)
                    {
                        allInternal[cellOffset + i] = gatheredInternal[proci][i];
                    }
                    cellOffset += gatheredInternal[proci].size();

                    forAll(gatheredBoundary[proci], i)
                    {
                        allBoundary[bfaceOffset + i] = gatheredBoundary[proci][i];
                    }
                    bfaceOffset += gatheredBoundary[proci].size();
                }

                OFstream os(timeDir / fieldName);
                writeFoamFileHeader(os, "volVectorField", fieldName);

                os << "dimensions      " << vfPtr->dimensions() << ";" << nl << nl
                   << "internalField   nonuniform List<vector>" << nl
                   << allInternal.size() << nl << '(' << nl;
                forAll(allInternal, i)
                {
                    os << allInternal[i] << nl;
                }
                os << ')' << ';' << nl << nl;

                os << "boundaryField" << nl << '{' << nl;
                os << "    oldInternalFaces" << nl
                   << "    {" << nl
                   << "        type            calculated;" << nl
                   << "        value           nonuniform List<vector>" << nl
                   << "        " << allBoundary.size() << nl
                   << "        (" << nl;
                forAll(allBoundary, i)
                {
                    os << "        " << allBoundary[i] << nl;
                }
                os << "        );" << nl
                   << "    }" << nl;
                os << '}' << nl;
            }
        }
    }
}


void Foam::functionObjects::spaceTimeWindowExtract::writeSubsetMesh()
{
    if (meshWritten_)
    {
        return;
    }

    // Gather face centres to master for boundaryData output
    vectorField allFaceCentres;
    scalarField allFaceWeights;
    labelList allFaceCellsInside;
    labelList allFaceCellsOutside;

    gatherFaceDataToMaster
    (
        allFaceCentres,
        allFaceWeights,
        allFaceCellsInside,
        allFaceCellsOutside
    );

    // In parallel, gather and write combined subset mesh
    if (Pstream::parRun())
    {
        gatherAndWriteSubsetMesh();
    }

    // Only master writes (or serial)
    if (Pstream::master())
    {
        // Create output directory structure
        fileName meshDir = outputDir_ / "constant" / "polyMesh";
        mkDir(meshDir);

        Info<< type() << " " << name() << ": Writing subset mesh data to "
            << outputDir_ << endl;

        if (!Pstream::parRun() && meshSubsetPtr_)
        {
            // Serial case: write full subset mesh
            const fvMesh& subMesh = meshSubsetPtr_->subMesh();
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

            // Boundary
            {
                OFstream os(meshDir / "boundary");
                writeFoamFileHeader(os, "polyBoundaryMesh", "boundary");

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

            Info<< "    Written subset mesh with "
                << subMesh.nCells() << " cells" << endl;
        }
        // Note: In parallel mode, extractionBox is written by gatherAndWriteSubsetMesh()

        // Write face centres for reference (used by spaceTimeWindow BC)
        fileName boundaryDir = outputDir_ / "constant" / "boundaryData" / "oldInternalFaces";
        mkDir(boundaryDir);

        {
            OFstream os(boundaryDir / "points");
            writeFoamFileHeader(os, "vectorField", "points");
            os << allFaceCentres;
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

            // Write extraction box for reference
            os << nl << "// Extraction region" << nl;
            os.writeEntry("boxMin", boxMin_);
            os.writeEntry("boxMax", boxMax_);

            // Parallel info
            if (Pstream::parRun())
            {
                os << nl << "// Parallel extraction info" << nl;
                os.writeEntry("nProcs", Pstream::nProcs());
                os.writeEntry("nGlobalFaces", nGlobalFaces_);
            }

            // Write format info
            os << nl << "// Boundary data format" << nl;
            os.writeEntry("useDeltaVarint", useDeltaVarint_);
            if (useDeltaVarint_)
            {
                os.writeEntry("deltaVarintPrecision", deltaVarintPrecision_);
            }
            os.writeEntry("encrypted", useEncryption_);

            Info<< "    Written extraction metadata:" << nl
                << "        openfoamVersion = " << foamVersion::version << nl
                << "        openfoamApi = " << foamVersion::api << nl
                << "        solver = " << solverName << nl
                << "        deltaT = " << runTime.deltaTValue() << nl
                << "        adjustTimeStep = "
                << runTime.controlDict().getOrDefault<Switch>("adjustTimeStep", false) << nl
                << "        timePrecision = " << timePrecision << endl;
        }

        Info<< "    Written boundary with oldInternalFaces ("
            << allFaceCentres.size() << " faces)" << endl;
    }

    // Barrier to ensure master has finished writing before others continue
    if (Pstream::parRun())
    {
        UPstream::barrier(UPstream::worldComm);
    }

    meshWritten_ = true;
}


void Foam::functionObjects::spaceTimeWindowExtract::writeInitialFields()
{
    // In parallel mode, we skip writing initial fields here because:
    // 1. The mesh is not written by the extractor (only extractionBox parameters)
    // 2. The cell ordering in the gathered fields won't match the mesh that
    //    spaceTimeWindowInitCase creates from the source case
    //
    // Instead, spaceTimeWindowInitCase extracts initial fields when it creates
    // the mesh, ensuring consistent cell ordering.
    if (Pstream::parRun())
    {
        Info<< type() << " " << name()
            << ": Skipping initial fields in parallel mode" << nl
            << "    spaceTimeWindowInitCase will extract them from source case" << endl;
        return;
    }

    // Serial mode: write initial fields from subset mesh
    if (!meshSubsetPtr_)
    {
        WarningInFunction
            << "No mesh subset available for writing initial fields" << endl;
        return;
    }

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

    // Track this timestep (exact string name) - only on master to avoid duplicates
    if (Pstream::master())
    {
        extractedTimesteps_.append(timeName);
    }

    // Create output directory only on master
    fileName boundaryTimeDir =
        outputDir_ / "constant" / "boundaryData" / "oldInternalFaces" / timeName;

    if (Pstream::master())
    {
        mkDir(boundaryTimeDir);
    }

    // Barrier to ensure directory exists before writing
    if (Pstream::parRun())
    {
        UPstream::barrier(UPstream::worldComm);
    }

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

            // Gather to master and write
            if (Pstream::parRun())
            {
                tmp<scalarField> tgathered = gatherFieldToMaster(tfaceValues());
                if (Pstream::master())
                {
                    writeField(boundaryTimeDir, fieldName, tgathered());
                }
            }
            else
            {
                writeField(boundaryTimeDir, fieldName, tfaceValues());
            }
            continue;
        }

        // Try vector field
        const auto* vfPtr = mesh_.findObject<volVectorField>(fieldName);
        if (vfPtr)
        {
            tmp<vectorField> tfaceValues = interpolateToFaces(*vfPtr);

            // Gather to master and write
            if (Pstream::parRun())
            {
                tmp<vectorField> tgathered = gatherFieldToMaster(tfaceValues());
                if (Pstream::master())
                {
                    writeField(boundaryTimeDir, fieldName, tgathered());
                }
            }
            else
            {
                writeField(boundaryTimeDir, fieldName, tfaceValues());
            }
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

    // For parallel runs with processor boundary faces, we need to exchange
    // field values for remote cells
    Field<Type> remoteValues;

    if (Pstream::parRun() && hasProcessorBoundaryFaces_)
    {
        // Exchange field values across processor boundaries
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send our cell values to neighbours
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                // Send field values for cells adjacent to this patch
                Field<Type> localValues(pp.size());
                forAll(pp, i)
                {
                    label owner = mesh_.faceOwner()[pp.start() + i];
                    localValues[i] = vf[owner];
                }

                UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
                toNbr << localValues;
            }
        }

        pBufs.finishedSends();

        // Build map of patch face index to received value
        Map<Type> patchFaceToValue;

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                Field<Type> receivedValues(pp.size());
                UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                fromNbr >> receivedValues;

                // Store received values indexed by processor + patch face
                forAll(receivedValues, i)
                {
                    label facei = pp.start() + i;
                    patchFaceToValue.insert(facei, receivedValues[i]);
                }
            }
        }

        // Now interpolate, using remote values where needed
        forAll(faceValues, facei)
        {
            const Type& valueInside = vf[faceCellsInside_[facei]];
            const scalar w = faceWeights_[facei];

            // Check if outside cell is remote
            label origFaceIdx = oldInternalFaceIndices_[facei];
            auto iter = patchFaceToValue.find(origFaceIdx);

            if (iter.good())
            {
                // Remote cell - use exchanged value
                const Type& valueOutside = iter.val();
                faceValues[facei] = w * valueInside + (1.0 - w) * valueOutside;
            }
            else
            {
                // Local cell
                const Type& valueOutside = vf[faceCellsOutside_[facei]];
                faceValues[facei] = w * valueInside + (1.0 - w) * valueOutside;
            }
        }
    }
    else
    {
        // Serial or no processor boundary faces - simple case
        forAll(faceValues, facei)
        {
            const Type& valueInside = vf[faceCellsInside_[facei]];
            const Type& valueOutside = vf[faceCellsOutside_[facei]];
            const scalar w = faceWeights_[facei];

            faceValues[facei] = w * valueInside + (1.0 - w) * valueOutside;
        }
    }

    return tfaceValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::spaceTimeWindowExtract::gatherFieldToMaster
(
    const Field<Type>& localField
) const
{
    if (!Pstream::parRun())
    {
        return tmp<Field<Type>>::New(localField);
    }

    // Gather all local fields to master
    List<Field<Type>> gatheredFields(Pstream::nProcs());
    gatheredFields[Pstream::myProcNo()] = localField;
    Pstream::gatherList(gatheredFields);

    if (Pstream::master())
    {
        // Combine all data
        label totalSize = 0;
        forAll(gatheredFields, proci)
        {
            totalSize += gatheredFields[proci].size();
        }

        auto tresult = tmp<Field<Type>>::New(totalSize);
        Field<Type>& result = tresult.ref();

        label offset = 0;
        forAll(gatheredFields, proci)
        {
            forAll(gatheredFields[proci], i)
            {
                result[offset + i] = gatheredFields[proci][i];
            }
            offset += gatheredFields[proci].size();
        }

        return tresult;
    }

    // Non-master processors return empty field
    return tmp<Field<Type>>::New();
}


template<class Type>
void Foam::functionObjects::spaceTimeWindowExtract::writeField
(
    const fileName& dir,
    const word& fieldName,
    const Field<Type>& field
) const
{
    // Check if we should use delta-varint codec (only for scalar and vector)
    if (useDeltaVarint_)
    {
        if constexpr (std::is_same_v<Type, scalar>)
        {
#ifdef FOAM_USE_SODIUM
            if (useEncryption_)
            {
                // Compress to buffer, then encrypt
                std::vector<uint8_t> compressed = deltaVarintCodec::encode(field, deltaVarintPrecision_);
                std::vector<uint8_t> publicKey = sodiumCrypto::fromBase64(publicKeyBase64_);

                // Write encrypted file with .dvz.enc extension
                fileName outPath = dir / (fieldName + "." + deltaVarintCodec::fileExtension() + "." + sodiumCrypto::fileExtension);
                sodiumCrypto::encryptToFile(outPath, compressed, publicKey);
                return;
            }
#endif
            fileName outPath = dir / (fieldName + "." + deltaVarintCodec::fileExtension());
            deltaVarintCodec::write(outPath, field, deltaVarintPrecision_);
            return;
        }
        else if constexpr (std::is_same_v<Type, vector>)
        {
#ifdef FOAM_USE_SODIUM
            if (useEncryption_)
            {
                // Compress to buffer, then encrypt
                std::vector<uint8_t> compressed = deltaVarintCodec::encode(field, deltaVarintPrecision_);
                std::vector<uint8_t> publicKey = sodiumCrypto::fromBase64(publicKeyBase64_);

                // Write encrypted file with .dvz.enc extension
                fileName outPath = dir / (fieldName + "." + deltaVarintCodec::fileExtension() + "." + sodiumCrypto::fileExtension);
                sodiumCrypto::encryptToFile(outPath, compressed, publicKey);
                return;
            }
#endif
            fileName outPath = dir / (fieldName + "." + deltaVarintCodec::fileExtension());
            deltaVarintCodec::write(outPath, field, deltaVarintPrecision_);
            return;
        }
        // For other types, fall through to standard format
    }

#ifdef FOAM_USE_SODIUM
    // Encryption without delta-varint: encrypt raw binary data
    if (useEncryption_)
    {
        if constexpr (std::is_same_v<Type, scalar> || std::is_same_v<Type, vector>)
        {
            // Serialize field to binary buffer
            std::vector<uint8_t> buffer;

            // Write a simple binary format: size + raw data
            label fieldSize = field.size();
            buffer.resize(sizeof(label) + fieldSize * sizeof(Type));

            std::memcpy(buffer.data(), &fieldSize, sizeof(label));
            std::memcpy(buffer.data() + sizeof(label), field.cdata(), fieldSize * sizeof(Type));

            std::vector<uint8_t> publicKey = sodiumCrypto::fromBase64(publicKeyBase64_);

            // Write encrypted file with .enc extension
            fileName outPath = dir / (fieldName + "." + sodiumCrypto::fileExtension);
            sodiumCrypto::encryptToFile(outPath, buffer, publicKey);
            return;
        }
    }
#endif

    // Standard OpenFOAM format
    OFstream os
    (
        dir / fieldName,
        IOstreamOption(writeFormat_, writeCompression_, IOstreamOption::currentVersion)
    );

    const word formatStr =
        (writeFormat_ == IOstreamOption::BINARY) ? "binary" : "ascii";

    os  << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      " << formatStr << ";" << nl
        << "    class       " << Field<Type>::typeName << ";" << nl
        << "    object      " << fieldName << ";" << nl
        << "}" << nl << nl;

    // Write field data (format is determined by stream setting)
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
    writeFormat_(IOstreamOption::ASCII),
    writeCompression_(IOstreamOption::UNCOMPRESSED),
    useDeltaVarint_(false),
    deltaVarintPrecision_(6),
    useEncryption_(false),
    publicKeyBase64_(),
    meshSubsetPtr_(),
    subsetInitialized_(false),
    meshWritten_(false),
    oldInternalFaceIndices_(),
    faceCellsInside_(),
    faceCellsOutside_(),
    faceWeights_(),
    faceCentres_(),
    extractedTimesteps_(),
    globalFaceIndexPtr_(),
    nGlobalFaces_(0),
    outsideCellMapPtr_(),
    remoteFaceIndices_(),
    remoteCellIndices_(),
    remoteCellProcs_(),
    hasProcessorBoundaryFaces_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::spaceTimeWindowExtract::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    // Read box definition (two points: min and max)
    Tuple2<point, point> boxDef = dict.get<Tuple2<point, point>>("box");
    boxMin_ = boxDef.first();
    boxMax_ = boxDef.second();

    dict.readEntry("outputDir", outputDir_);
    dict.readEntry("fields", fieldNames_);

    // Read write format (ascii, binary, or deltaVarint), default to ascii
    const word writeFormatStr = dict.getOrDefault<word>("writeFormat", "ascii");
    if (writeFormatStr == "binary")
    {
        writeFormat_ = IOstreamOption::BINARY;
        useDeltaVarint_ = false;
    }
    else if (writeFormatStr == "ascii")
    {
        writeFormat_ = IOstreamOption::ASCII;
        useDeltaVarint_ = false;
    }
    else if (writeFormatStr == "deltaVarint" || writeFormatStr == "dvz")
    {
        // Delta-varint codec for high-efficiency compression
        useDeltaVarint_ = true;
        writeFormat_ = IOstreamOption::BINARY;  // Not used but set for consistency
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Invalid writeFormat '" << writeFormatStr << "'" << nl
            << "Valid options are: ascii, binary, deltaVarint" << nl
            << exit(FatalIOError);
    }

    // Read delta-varint precision (significant digits)
    deltaVarintPrecision_ = dict.getOrDefault<label>("deltaVarintPrecision", 6);

    // Read write compression (on/off/true/false/compressed/uncompressed), default to off
    const word writeCompressionStr = dict.getOrDefault<word>("writeCompression", "off");
    if (writeCompressionStr == "on" || writeCompressionStr == "true"
        || writeCompressionStr == "compressed" || writeCompressionStr == "zstd"
        || writeCompressionStr == "gzip")
    {
        writeCompression_ = IOstreamOption::COMPRESSED;
    }
    else if (writeCompressionStr == "off" || writeCompressionStr == "false"
             || writeCompressionStr == "uncompressed")
    {
        writeCompression_ = IOstreamOption::UNCOMPRESSED;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Invalid writeCompression '" << writeCompressionStr << "'" << nl
            << "Valid options are: on, off, compressed, uncompressed" << nl
            << exit(FatalIOError);
    }

    // Make outputDir absolute if relative
    if (!outputDir_.isAbsolute())
    {
        outputDir_ = mesh_.time().globalPath() / outputDir_;
    }

    Info<< type() << " " << name() << ":" << nl
        << "    box:       " << boxMin_ << " " << boxMax_ << nl
        << "    outputDir: " << outputDir_ << nl
        << "    fields:    " << fieldNames_ << nl
        << "    writeFormat: " << writeFormatStr;

    if (useDeltaVarint_)
    {
        Info<< " (precision: " << deltaVarintPrecision_ << " digits)";
    }
    Info<< nl;

    if (!useDeltaVarint_)
    {
        Info<< "    writeCompression: " << writeCompressionStr << endl;
    }
    else
    {
        Info<< endl;
    }

    // Read optional encryption settings
    useEncryption_ = dict.getOrDefault<bool>("encrypt", false);
    if (useEncryption_)
    {
#ifdef FOAM_USE_SODIUM
        if (!sodiumCrypto::available())
        {
            FatalIOErrorInFunction(dict)
                << "Encryption requested but libsodium failed to initialize" << nl
                << exit(FatalIOError);
        }

        if (!dict.readIfPresent("publicKey", publicKeyBase64_))
        {
            FatalIOErrorInFunction(dict)
                << "Encryption enabled but 'publicKey' not specified" << nl
                << "Generate a keypair with spaceTimeWindowKeygen utility" << nl
                << exit(FatalIOError);
        }

        // Validate the public key
        std::vector<uint8_t> pubKey = sodiumCrypto::fromBase64(publicKeyBase64_);
        if (pubKey.size() != sodiumCrypto::PUBLIC_KEY_SIZE)
        {
            FatalIOErrorInFunction(dict)
                << "Invalid public key size. Expected "
                << sodiumCrypto::PUBLIC_KEY_SIZE << " bytes, got "
                << pubKey.size() << " bytes" << nl
                << "Ensure you're using a valid base64-encoded 32-byte public key" << nl
                << exit(FatalIOError);
        }

        Info<< "    encryption: enabled (sealed box)" << endl;
#else
        FatalIOErrorInFunction(dict)
            << "Encryption requested but library was built without libsodium" << nl
            << "Rebuild with FOAM_USE_SODIUM=1 to enable encryption" << nl
            << exit(FatalIOError);
#endif
    }

    return true;
}


bool Foam::functionObjects::spaceTimeWindowExtract::execute()
{
    // Ensure subset is initialized
    initializeSubset();

    // Write mesh once (on first execute)
    if (!meshWritten_)
    {
        // In parallel mode, force the solver to write current fields to disk
        // at the extraction start time. This ensures spaceTimeWindowInitCase
        // can extract initial fields from the reconstructed source case.
        if (Pstream::parRun())
        {
            Info<< type() << " " << name()
                << ": Forcing field write at extraction start time "
                << mesh_.time().timeName() << endl;

            // Use writeNow() to force immediate write of all fields
            const_cast<Time&>(mesh_.time()).writeNow();

            // Barrier to ensure all processors have finished writing
            UPstream::barrier(UPstream::worldComm);

            Info<< "    Field write complete. Run reconstructPar -time "
                << mesh_.time().timeName()
                << " before spaceTimeWindowInitCase" << endl;
        }

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

template Foam::tmp<Foam::Field<Foam::scalar>>
Foam::functionObjects::spaceTimeWindowExtract::gatherFieldToMaster<Foam::scalar>
(
    const Field<scalar>&
) const;

template Foam::tmp<Foam::Field<Foam::vector>>
Foam::functionObjects::spaceTimeWindowExtract::gatherFieldToMaster<Foam::vector>
(
    const Field<vector>&
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
