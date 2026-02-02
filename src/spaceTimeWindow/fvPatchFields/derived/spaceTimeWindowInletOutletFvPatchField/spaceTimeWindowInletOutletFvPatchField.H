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
    Foam::spaceTimeWindowInletOutletFvPatchField

Group
    grpInletBoundaryConditions

Description
    Inlet-outlet boundary condition for space-time window reconstruction.

    Applies Dirichlet (time-interpolated extracted values) at inflow faces
    and zeroGradient at outflow faces, based on flux direction.

    This is the recommended BC for velocity in space-time window reconstruction:
    - Preserves upstream flow structure (vortices entering domain)
    - Allows natural advection out (vortices exiting domain)
    - Does not overconstrain the system

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
        phi          | flux field name                  | no       | phi
        fixesValue   | report to adjustPhi              | no       | false
        allowTimeInterpolation | permit interpolation   | no       | true
        timeInterpolationScheme | linear or cubic       | no       | linear
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            spaceTimeWindowInletOutlet;
        dataDir         "constant/boundaryData";
        phi             phi;
        fixesValue      false;
        allowTimeInterpolation  true;
        timeInterpolationScheme cubic;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    spaceTimeWindowInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeWindowInletOutletFvPatchField_H
#define spaceTimeWindowInletOutletFvPatchField_H

#include "mixedFvPatchFields.H"
#include "instantList.H"
#include "barycentric2D.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class spaceTimeWindowInletOutletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class spaceTimeWindowInletOutletFvPatchField
:
    public mixedFvPatchField<Type>
{
public:

    //- Time interpolation scheme enumeration
    enum class timeInterpScheme
    {
        NONE,       //!< No interpolation (exact timestep match required)
        LINEAR,     //!< Linear interpolation (2 points)
        CUBIC       //!< Catmull-Rom cubic spline (4 points)
    };

private:

    // Private Data

        //- Path to boundaryData directory
        fileName dataDir_;

        //- Name of the field data file (defaults to field name)
        word fieldTableName_;

        //- Name of flux field
        word phiName_;

        //- If true, report to adjustPhi that this patch fixes values
        //  Default is false for inlet-outlet BC
        bool fixesValue_;

        //- If true, allow temporal interpolation for missing timesteps
        bool allowTimeInterpolation_;

        //- Time interpolation scheme (none, linear, or cubic)
        timeInterpScheme timeInterpolationScheme_;

        //- List of available time directories
        mutable instantList sampleTimes_;

        //- Time indices for interpolation [i0, i1, i2, i3] for cubic
        mutable label sampleTimeIndex0_;
        mutable label sampleTimeIndex1_;
        mutable label sampleTimeIndex2_;
        mutable label sampleTimeIndex3_;

        //- Cached field data at the 4 time levels
        mutable Field<Type> sampledValues0_;
        mutable Field<Type> sampledValues1_;
        mutable Field<Type> sampledValues2_;
        mutable Field<Type> sampledValues3_;

        //- Whether metadata has been validated
        mutable bool metadataValidated_;

        //- Spatial coarsening factor (1 = no coarsening, 2 = 2:1 in each dim)
        mutable label spatialCoarsenFactor_;

        //- Mapping from fine boundary data faces to coarse mesh faces
        //  coarseToFineMap_[coarseFaceI] = list of fine face indices
        mutable List<labelList> coarseToFineMap_;

        //- Whether spatial mapping has been initialized
        mutable bool spatialMapInitialized_;

        //- Refinement mode: coarse boundaryData to fine mesh interpolation
        //  For each fine face, stores: triangle vertex indices and weights
        mutable bool refinementMode_;

        //- Coarse points from boundaryData (only stored in refinement mode)
        mutable vectorField coarsePoints_;

        //- Triangle vertex indices for refinement interpolation
        //  triangles_[triI] = FixedList<label, 3> with 3 coarse point indices
        mutable List<FixedList<label, 3>> triangles_;

        //- Triangle index for each fine face
        mutable labelList fineToTriangle_;

        //- Barycentric weights for each fine face
        mutable List<barycentric2D> fineBaryWeights_;


    // Private Member Functions

        //- Initialize spatial coarsening map if needed
        void initSpatialCoarseningMap() const;

        //- Build Delaunay triangulation for refinement interpolation
        void buildRefinementTriangulation
        (
            const vectorField& coarsePoints,
            const vectorField& finePoints
        ) const;

        //- Find available time directories
        void findSampleTimes() const;

        //- Read field data from file
        Field<Type> readFieldData(const instant& timeDir) const;

        //- Update the sampled data
        void checkTable();

        //- Validate extraction metadata
        void validateMetadata() const;

        //- Centripetal Catmull-Rom cubic spline interpolation
        Type catmullRomCentripetal
        (
            const Type& p0,
            const Type& p1,
            const Type& p2,
            const Type& p3,
            scalar t0,
            scalar t1,
            scalar t2,
            scalar t3,
            scalar t
        ) const;


public:

    //- Runtime type information
    TypeName("spaceTimeWindowInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        spaceTimeWindowInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaceTimeWindowInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        spaceTimeWindowInletOutletFvPatchField
        (
            const spaceTimeWindowInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaceTimeWindowInletOutletFvPatchField
        (
            const spaceTimeWindowInletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new spaceTimeWindowInletOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        spaceTimeWindowInletOutletFvPatchField
        (
            const spaceTimeWindowInletOutletFvPatchField<Type>&,
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
                new spaceTimeWindowInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Return true if this patch field fixes a value
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
    #include "spaceTimeWindowInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
