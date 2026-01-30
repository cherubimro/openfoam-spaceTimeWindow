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

        //- Mesh subset helper
        autoPtr<fvMeshSubset> meshSubsetPtr_;

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

        //- Write subset mesh
        void writeSubsetMesh();

        //- Write initial fields (subset)
        void writeInitialFields();

        //- Write boundary data for current time
        void writeBoundaryData();

        //- Update metadata with timestep list (called at end)
        void updateMetadataTimesteps();

        //- Interpolate field to face values
        template<class Type>
        tmp<Field<Type>> interpolateToFaces
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
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
