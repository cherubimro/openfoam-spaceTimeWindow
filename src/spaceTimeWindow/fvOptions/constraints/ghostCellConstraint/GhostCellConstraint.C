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

#include "GhostCellConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::GhostCellConstraint<Type>::GhostCellConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    dataDir_(),
    ghostCellMap_(),
    cachedValues_(),
    cachedTimeName_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::GhostCellConstraint<Type>::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("dataDir", dataDir_);

        // Read field names
        coeffs_.readEntry("fields", fieldNames_);

        fv::option::resetApplied();

        // Read ghost cell map: data index -> cell index
        fileName mapPath = dataDir_ / "ghostCells" / "ghostCellMap";
        if (isFile(mapPath))
        {
            IFstream mapIs(mapPath);
            // Skip FoamFile header if present
            token firstToken(mapIs);
            if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
            {
                // Read and discard the header dictionary
                dictionary headerDict(mapIs);
            }
            else
            {
                mapIs.putBack(firstToken);
            }
            mapIs >> ghostCellMap_;

            Info<< "    GhostCellConstraint: read " << ghostCellMap_.size()
                << " ghost cell mappings from " << mapPath << endl;
        }
        else
        {
            // If no map file, use identity mapping (data[i] -> cells_[i])
            Info<< "    GhostCellConstraint: no ghostCellMap found, "
                << "using identity mapping (cells_ order)" << endl;
        }

        return true;
    }

    return false;
}


template<class Type>
void Foam::fv::GhostCellConstraint<Type>::readGhostData
(
    const word& timeName
) const
{
    fileName ghostDir = dataDir_ / "ghostCells" / timeName;

    cachedValues_.clear();

    for (const word& fieldName : fieldNames_)
    {
        fileName fieldPath = ghostDir / fieldName;

        if (!isFile(fieldPath))
        {
            FatalErrorInFunction
                << "Ghost cell data file not found: " << fieldPath << nl
                << "    Ensure spaceTimeWindowExtract was run with "
                << "extractGhostCells enabled" << nl
                << exit(FatalError);
        }

        IFstream is(fieldPath);

        // Skip FoamFile header if present
        token firstToken(is);
        if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
        {
            dictionary headerDict(is);
        }
        else
        {
            is.putBack(firstToken);
        }

        Field<Type> values(is);
        cachedValues_.set(fieldName, values);
    }

    cachedTimeName_ = timeName;
}


template<class Type>
void Foam::fv::GhostCellConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    const word timeName = mesh_.time().timeName();

    // Read ghost data if not cached for this timestep
    if (cachedTimeName_ != timeName)
    {
        readGhostData(timeName);
    }

    const word& fieldName = fieldNames_[fieldi];

    if (!cachedValues_.found(fieldName))
    {
        WarningInFunction
            << "No ghost data for field " << fieldName
            << " at time " << timeName << endl;
        return;
    }

    const Field<Type>& ghostValues = cachedValues_[fieldName];

    // Build cell values list for setValues
    if (ghostCellMap_.size() > 0)
    {
        // Use explicit mapping: ghostCellMap maps data index to mesh cell
        // cells_ from cellSetOption contains the ghost zone cells
        // We need to map: for each cell in cells_, find the corresponding
        // ghost data value

        // Build reverse map: mesh cell index -> position in cells_ list
        Map<label> cellToIdx;
        forAll(cells_, i)
        {
            cellToIdx.insert(cells_[i], i);
        }

        // ghostCellMap_[dataIdx] = mesh cell index
        // ghostValues[dataIdx] = field value
        List<Type> cellValues(cells_.size(), pTraits<Type>::zero);

        forAll(ghostCellMap_, dataIdx)
        {
            label meshCell = ghostCellMap_[dataIdx];
            auto iter = cellToIdx.find(meshCell);
            if (iter.good())
            {
                cellValues[iter.val()] = ghostValues[dataIdx];
            }
        }

        eqn.setValues(cells_, cellValues);
    }
    else
    {
        // Identity mapping: ghostValues[i] corresponds to cells_[i]
        if (ghostValues.size() != cells_.size())
        {
            FatalErrorInFunction
                << "Ghost data size (" << ghostValues.size()
                << ") does not match ghost cell count (" << cells_.size()
                << ") for field " << fieldName << nl
                << exit(FatalError);
        }

        List<Type> cellValues(ghostValues.size());
        forAll(ghostValues, i)
        {
            cellValues[i] = ghostValues[i];
        }

        eqn.setValues(cells_, cellValues);
    }
}


// ************************************************************************* //
