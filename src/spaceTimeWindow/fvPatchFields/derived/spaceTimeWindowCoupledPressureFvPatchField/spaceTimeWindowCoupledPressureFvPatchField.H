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
    Foam::spaceTimeWindowCoupledPressureFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Coupled pressure boundary condition for space-time window reconstruction.

    This BC blends Dirichlet (extracted pressure) and Neumann (extracted gradient)
    conditions based on the local flux. Two blending modes are available:

    flowDirection mode (original):
    - At inflow faces: Uses more Dirichlet weight (pressure enters domain)
    - At outflow faces: Uses more Neumann weight (pressure from interior dominates)

    fluxMagnitude mode (recommended):
    - At stagnant faces (zero flux): Uses Dirichlet (safe to prescribe pressure)
    - At active flow faces: Uses Neumann (gradient for continuity compatibility)

    The blending provides proper pressure coupling at the extraction boundary,
    which is critical because pressure is governed by an elliptic equation and
    cutting the domain artificially severs the pressure influence from outside.

    Mathematical formulation (mixed BC):
    \f[
        p_f = w \cdot p_{Dirichlet} + (1-w) \cdot (p_C + \nabla p_{Neumann} \cdot d)
    \f]

    where:
    - w = blending weight (function of flux)
    - p_Dirichlet = time-interpolated extracted pressure
    - p_C = cell center value (from solver)
    - grad p_Neumann = time-interpolated extracted normal gradient
    - d = distance from cell center to face center

Usage
    \table
        Property     | Description                      | Required | Default
        dataDir      | boundaryData directory path      | no       | constant/boundaryData
        phi          | flux field name                  | no       | phi
        blendingMode | flowDirection or fluxMagnitude   | no       | fluxMagnitude
        dirichletWeight | base Dirichlet weight         | no       | 0.8
        neumannWeight | base Neumann weight             | no       | 0.8
        transitionWidth | flux transition width         | no       | 0.1
        allowTimeInterpolation | permit interpolation   | no       | true
        timeInterpolationScheme | linear or cubic       | no       | cubic
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            spaceTimeWindowCoupledPressure;
        dataDir         "constant/boundaryData";
        phi             phi;
        blendingMode    fluxMagnitude;
        dirichletWeight 0.8;
        neumannWeight   0.8;
        transitionWidth 0.1;
        allowTimeInterpolation  true;
        timeInterpolationScheme cubic;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    spaceTimeWindowCoupledPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef spaceTimeWindowCoupledPressureFvPatchScalarField_H
#define spaceTimeWindowCoupledPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "instantList.H"
#include "spatialInterpolation2D.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class spaceTimeWindowCoupledPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class spaceTimeWindowCoupledPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
public:

    //- Time interpolation scheme enumeration
    enum class timeInterpScheme
    {
        NONE,       //!< No interpolation (exact timestep match required)
        LINEAR,     //!< Linear interpolation (2 points)
        CUBIC       //!< Catmull-Rom cubic spline (4 points)
    };

    //- Blending mode enumeration
    enum class blendingMode
    {
        FLOW_DIRECTION,  //!< Dirichlet at inflow, Neumann at outflow (original)
        FLUX_MAGNITUDE   //!< Dirichlet at stagnant faces, Neumann at flowing faces
    };

private:

    // Private Data

        //- Path to boundaryData directory
        fileName dataDir_;

        //- Name of the pressure field data file (defaults to "p")
        word pressureFieldName_;

        //- Name of the pressure gradient field data file (defaults to "gradp")
        word gradientFieldName_;

        //- Name of flux field
        word phiName_;

        //- Base Dirichlet weight at inflow (default 0.8)
        scalar dirichletWeight_;

        //- Base Neumann weight at outflow (default 0.8)
        scalar neumannWeight_;

        //- Flux transition width for blending (default 0.1)
        scalar transitionWidth_;

        //- Blending mode (flowDirection or fluxMagnitude)
        blendingMode blendingMode_;

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

        //- Cached pressure values at the 4 time levels
        mutable scalarField sampledPressure0_;
        mutable scalarField sampledPressure1_;
        mutable scalarField sampledPressure2_;
        mutable scalarField sampledPressure3_;

        //- Cached gradient values at the 4 time levels
        mutable scalarField sampledGradient0_;
        mutable scalarField sampledGradient1_;
        mutable scalarField sampledGradient2_;
        mutable scalarField sampledGradient3_;

        //- Whether metadata has been validated
        mutable bool metadataValidated_;

        //- Spatial interpolation helper (for multi-resolution support)
        mutable autoPtr<spatialInterpolation2D> spatialInterpPtr_;

        //- Bounding box for spatial interpolation
        mutable boundBox boundingBox_;

        //- Whether spatial interpolation has been initialized
        mutable bool spatialInterpInitialized_;


    // Private Member Functions

        //- Find available time directories
        void findSampleTimes() const;

        //- Read pressure field data from file
        scalarField readPressureData(const instant& timeDir) const;

        //- Read gradient field data from file
        scalarField readGradientData(const instant& timeDir) const;

        //- Update the sampled data
        void checkTable();

        //- Validate extraction metadata
        void validateMetadata() const;

        //- Initialize spatial interpolation if needed (for multi-resolution)
        void initSpatialInterpolation() const;

        //- Compute blending weights based on flux direction
        //  Returns valueFraction for mixed BC (1 = pure Dirichlet, 0 = pure Neumann)
        tmp<scalarField> computeBlendingWeights() const;

        //- Centripetal Catmull-Rom cubic spline interpolation
        scalar catmullRomCentripetal
        (
            scalar p0,
            scalar p1,
            scalar p2,
            scalar p3,
            scalar t0,
            scalar t1,
            scalar t2,
            scalar t3,
            scalar t
        ) const;


public:

    //- Runtime type information
    TypeName("spaceTimeWindowCoupledPressure");


    // Constructors

        //- Construct from patch and internal field
        spaceTimeWindowCoupledPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        spaceTimeWindowCoupledPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given field onto a new patch
        spaceTimeWindowCoupledPressureFvPatchScalarField
        (
            const spaceTimeWindowCoupledPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        spaceTimeWindowCoupledPressureFvPatchScalarField
        (
            const spaceTimeWindowCoupledPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new spaceTimeWindowCoupledPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        spaceTimeWindowCoupledPressureFvPatchScalarField
        (
            const spaceTimeWindowCoupledPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new spaceTimeWindowCoupledPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
