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

Class
    Foam::functionObjects::spaceTimeWindowExtract

Group
    grpFieldFunctionObjects

Description
    Function object to extract space-time window data during simulation.

    Outputs:
    - Subset mesh (once, at first write)
    - Face-interpolated boundary values at each timestep

    The face values are computed using linear interpolation between
    adjacent cell centers, which is exact for OpenFOAM's cell-centered FVM.

    Output structure:
    \verbatim
    outputDir/
        constant/
            polyMesh/           # Subset mesh files
            boundaryData/
                oldInternalFaces/
                    points      # Face centres (for reference only)
                    <time>/
                        U       # Face-interpolated values
                        p
        <startTime>/            # Initial subset fields
            U
            p
    \endverbatim

Usage
    \verbatim
    spaceTimeWindowExtract1
    {
        type            spaceTimeWindowExtract;
        libs            (spaceTimeWindow);

        // Define region directly with a box (fully internal to domain)
        box             (0.05 -0.25 0.01) (0.90 0.25 0.38);

        outputDir       "../ufr2-02-subset"; // Output case directory
        fields          (U p);              // Fields to extract

        // Optional: write format for boundary data
        // Options: ascii, binary, deltaVarint
        // Default: ascii
        // deltaVarint provides ~2.7% compression ratio (vs raw binary)
        writeFormat     deltaVarint;

        // Optional: precision for deltaVarint (significant digits)
        // Default: 6 (matches OpenFOAM's default writePrecision)
        deltaVarintPrecision  6;

        // Optional: compression (gzip) - ignored for deltaVarint
        writeCompression on;

        // Optional: encryption with libsodium (requires FOAM_USE_SODIUM=1 build)
        // encrypt         true;
        // publicKey       "base64-encoded-32-byte-public-key";

        writeControl    timeStep;
        writeInterval   1;
    }
    \endverbatim

SourceFiles
    spaceTimeWindowExtract.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_spaceTimeWindowExtract_H
#define functionObjects_spaceTimeWindowExtract_H

#include "fvMeshFunctionObject.H"
#include "fvMeshSubset.H"
#include "volFields.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "IOstreamOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class spaceTimeWindowExtract Declaration
\*---------------------------------------------------------------------------*/

class spaceTimeWindowExtract
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Box min point
        point boxMin_;

        //- Box max point
        point boxMax_;

        //- Output directory
        fileName outputDir_;

        //- Fields to extract
        wordList fieldNames_;

        //- Write format for boundary data (ascii or binary)
        IOstreamOption::streamFormat writeFormat_;

        //- Write compression for boundary data (compressed or uncompressed)
        IOstreamOption::compressionType writeCompression_;

        //- Use delta-varint codec for high-efficiency compression
        bool useDeltaVarint_;

        //- Precision for delta-varint codec (number of significant digits)
        label deltaVarintPrecision_;

        //- Enable encryption of boundary data
        bool useEncryption_;

        //- Public key for encryption (base64 encoded)
        string publicKeyBase64_;

        //- Mesh subset helper
        autoPtr<fvMeshSubset> meshSubsetPtr_;

        //- Whether subset has been initialized
        bool subsetInitialized_;

        //- Whether mesh has been written
        bool meshWritten_;

        //- Face indices of oldInternalFaces (in original mesh)
        labelList oldInternalFaceIndices_;

        //- Cell indices inside subset for each oldInternalFace
        labelList faceCellsInside_;

        //- Cell indices outside subset for each oldInternalFace
        labelList faceCellsOutside_;

        //- Interpolation weights for each face
        scalarList faceWeights_;

        //- Face centres for oldInternalFaces
        vectorField faceCentres_;

        //- List of extracted timestep names (precise strings)
        wordList extractedTimesteps_;

        // Parallel support data

        //- Global index for combining face data across processors
        autoPtr<globalIndex> globalFaceIndexPtr_;

        //- Total number of boundary faces across all processors
        label nGlobalFaces_;

        //- Map for distributing "outside" cell data in parallel
        //  For processor boundary faces where the outside cell is remote
        autoPtr<mapDistribute> outsideCellMapPtr_;

        //- Local indices of faces that need remote outside cell data
        labelList remoteFaceIndices_;

        //- Remote cell indices (in source processor numbering) needed for interpolation
        labelList remoteCellIndices_;

        //- Processor IDs for remote cells
        labelList remoteCellProcs_;

        //- Whether we have processor boundary faces crossing the extraction region
        bool hasProcessorBoundaryFaces_;

        //- Counter for timesteps since extraction started
        label timestepCount_;

        //- Whether t_2 (third timestep) has been written to disk
        bool t2Written_;


    // Protected Member Functions

        //- Write FoamFile header
        void writeFoamFileHeader
        (
            Ostream& os,
            const word& className,
            const word& objectName
        ) const;

        //- Initialize the mesh subset
        void initializeSubset();

        //- Initialize parallel communication structures
        void initializeParallelComm();

        //- Identify faces on processor boundaries that cross extraction region
        void identifyProcessorBoundaryFaces
        (
            const boundBox& bb,
            const labelList& cellsInBox
        );

        //- Gather face data to master processor
        void gatherFaceDataToMaster
        (
            vectorField& allFaceCentres,
            scalarField& allFaceWeights,
            labelList& allFaceCellsInside,
            labelList& allFaceCellsOutside
        ) const;

        //- Gather subset mesh to master and write combined mesh
        void gatherAndWriteSubsetMesh();

        //- Gather initial fields to master and write combined fields
        void gatherAndWriteInitialFields();

        //- Write subset mesh
        void writeSubsetMesh();

        //- Write initial fields (subset)
        void writeInitialFields();

        //- Write boundary data for current time
        void writeBoundaryData();

        //- Update metadata with timestep list (called at end)
        void updateMetadataTimesteps();

        //- Interpolate field to face values (handles parallel)
        template<class Type>
        tmp<Field<Type>> interpolateToFaces
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Gather interpolated field values to master
        template<class Type>
        tmp<Field<Type>> gatherFieldToMaster
        (
            const Field<Type>& localField
        ) const;

        //- Write field to file in OpenFOAM format
        template<class Type>
        void writeField
        (
            const fileName& dir,
            const word& fieldName,
            const Field<Type>& field
        ) const;


private:

        //- No copy construct
        spaceTimeWindowExtract(const spaceTimeWindowExtract&) = delete;

        //- No copy assignment
        void operator=(const spaceTimeWindowExtract&) = delete;


public:

    //- Runtime type information
    TypeName("spaceTimeWindowExtract");


    // Constructors

        //- Construct from Time and dictionary
        spaceTimeWindowExtract
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor - writes final metadata
    virtual ~spaceTimeWindowExtract();


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Execute (called at each time step)
        virtual bool execute();

        //- Write (called at write intervals)
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
