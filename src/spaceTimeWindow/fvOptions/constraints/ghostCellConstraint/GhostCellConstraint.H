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
    Foam::fv::GhostCellConstraint

Group
    grpFvOptionsConstraints

Description
    Constrains fields in ghost cells to values read from extracted global
    simulation data.  Works with the spaceTimeWindow ghost cell workflow:
    spaceTimeWindowExtract writes cell-centre values for ghost cells each
    timestep; this fvOption reads those values and overrides the
    corresponding matrix rows via eqn.setValues().

    The ghost cell zone is identified by a cellZone (typically "ghostCells"
    created by spaceTimeWindowInitCase -ghostCells).

Usage
    \verbatim
    ghostU
    {
        type            vectorGhostCellConstraint;
        selectionMode   cellZone;
        cellZone        ghostCells;
        dataDir         "constant/boundaryData";
        fields          (U);
    }

    ghostP
    {
        type            scalarGhostCellConstraint;
        selectionMode   cellZone;
        cellZone        ghostCells;
        dataDir         "constant/boundaryData";
        fields          (p);
    }
    \endverbatim

SourceFiles
    GhostCellConstraint.C
    ghostCellConstraints.C

\*---------------------------------------------------------------------------*/

#ifndef GhostCellConstraint_H
#define GhostCellConstraint_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class GhostCellConstraint Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class GhostCellConstraint
:
    public fv::cellSetOption
{
    // Private Data

        //- Path to boundaryData directory (contains ghostCells/)
        fileName dataDir_;

        //- Mapping from ghost data index to cell index in the mesh
        //  Read from dataDir/ghostCells/ghostCellMap
        labelList ghostCellMap_;

        //- Cached ghost field values per field name
        mutable HashTable<Field<Type>> cachedValues_;

        //- Time name for which values are cached
        mutable word cachedTimeName_;


    // Private Member Functions

        //- Read ghost cell data for all fields at given time
        void readGhostData(const word& timeName) const;

        //- No copy construct
        GhostCellConstraint(const GhostCellConstraint&) = delete;

        //- No copy assignment
        void operator=(const GhostCellConstraint&) = delete;


public:

    //- Runtime type information
    TypeName("GhostCellConstraint");


    // Constructors

        //- Construct from components
        GhostCellConstraint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~GhostCellConstraint() = default;


    // Member Functions

        //- Read source dictionary
        virtual bool read(const dictionary& dict);

        //- Set values on ghost cells
        virtual void constrain(fvMatrix<Type>& eqn, const label fieldi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GhostCellConstraint.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
