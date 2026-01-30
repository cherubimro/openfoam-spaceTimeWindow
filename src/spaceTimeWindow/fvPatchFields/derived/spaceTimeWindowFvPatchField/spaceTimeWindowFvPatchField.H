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

Class
    Foam::spaceTimeWindowFvPatchField

Group
    grpInletBoundaryConditions

Description
    Boundary condition for space-time window reconstruction.

    Reads pre-computed face values from boundaryData directory structure.
    Unlike timeVaryingMappedFixedValue, this BC:
    - Does NO spatial interpolation (values are pre-computed for exact faces)
    - Only does temporal interpolation between available timesteps
    - Expects face values computed via linear cell-to-face interpolation

    The boundaryData structure:
    \verbatim
        constant/boundaryData/<patchName>/
            points              (optional, not used - no spatial mapping)
            <time1>/
                <fieldName>     (face values as Field<Type>)
            <time2>/
                <fieldName>
            ...
    \endverbatim

Usage
    \table
        Property     | Description                      | Required | Default
        dataDir      | boundaryData directory path      | no       | constant/boundaryData
        setAverage   | adjust mapped field average      | no       | false
        offset       | offset value added to field      | no       | Zero
        fixesValue   | report to adjustPhi that values are fixed | no | true
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            spaceTimeWindow;
        dataDir         "constant/boundaryData";
        setAverage      false;
        offset          (0 0 0);    // For vectors
        fixesValue      true;       // Tell adjustPhi to exclude this patch
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    spaceTimeWindowFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeWindowFvPatchField_H
#define spaceTimeWindowFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class spaceTimeWindowFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class spaceTimeWindowFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Path to boundaryData directory
        fileName dataDir_;

        //- Name of the field data file (defaults to field name)
        word fieldTableName_;

        //- If true, adjust average of mapped field to match offset
        bool setAverage_;

        //- Offset value added to field
        Type offset_;

        //- If true, report to adjustPhi that this patch fixes values
        bool fixesValue_;

        //- List of available time directories
        mutable instantList sampleTimes_;

        //- Current time index in sampleTimes_
        mutable label startSampleTime_;

        //- End time index for interpolation
        mutable label endSampleTime_;

        //- Cached field data at start time
        mutable Field<Type> startSampledValues_;

        //- Cached field data at end time
        mutable Field<Type> endSampledValues_;

        //- Whether metadata has been validated
        mutable bool metadataValidated_;


    // Private Member Functions

        //- Find available time directories
        void findSampleTimes() const;

        //- Read field data from file
        Field<Type> readFieldData(const instant& timeDir) const;

        //- Update the sampled data
        void checkTable();

        //- Validate extraction metadata against current time settings
        void validateMetadata() const;


public:

    //- Runtime type information
    TypeName("spaceTimeWindow");


    // Constructors

        //- Construct from patch and internal field
        spaceTimeWindowFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaceTimeWindowFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given spaceTimeWindowFvPatchField
        //  onto a new patch
        spaceTimeWindowFvPatchField
        (
            const spaceTimeWindowFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaceTimeWindowFvPatchField
        (
            const spaceTimeWindowFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new spaceTimeWindowFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        spaceTimeWindowFvPatchField
        (
            const spaceTimeWindowFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new spaceTimeWindowFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Return true if this patch field fixes a value
            //  Used by adjustPhi to determine flux correction strategy
            virtual bool fixesValue() const
            {
                return fixesValue_;
            }

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "spaceTimeWindowFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
