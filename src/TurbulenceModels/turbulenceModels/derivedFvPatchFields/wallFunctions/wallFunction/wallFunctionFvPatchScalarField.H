/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::wallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description

Usage

See also
  - Foam::fixedValueFvPatchField

SourceFiles
    wallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallFunctionFvPatchScalarField_H
#define wallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
               Class wallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Member Functions

        //- Estimate the y+ at the intersection of the two sublayers
        static scalar yPlusLam(const scalar kappa, const scalar E);

        //- Check if the type of the patch is wall
        void checkType();


protected:

    // Protected Enumerations

        //- Options for the blending treatment of viscous and inertial sublayers
        enum blendingType
        {
            STEPWISE,       //!< "Stepwise switch (discontinuous)"
            MAX,            //!< "Maximum value switch (discontinuous)"
            BINOMIAL,       //!< "Binomial blending (smooth)"
            EXPONENTIAL,    //!< "Exponential blending (smooth)"
            TANH            //!< "Tanh blending (smooth)"
        };

        //- Names for blendingType
        static const Enum<blendingType> blendingTypeNames;


    // Protected Data

        //- Empirical model coefficient
        const scalar Cmu_;

        //- von Karman constant
        const scalar kappa_;

        //- Wall roughness parameter
        const scalar E_;

        //- Estimated y+ value at the intersection
        //- of the viscous and inertial sublayers
        const scalar yPlusLam_;

        //- Blending treatment
        enum blendingType blending_;

        //- Blending exponent being used when
        //- blendingType is blendingType::BINOMIAL
        scalar n_;


    // Protected Member Functions

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("wallFunction");


    // Constructors

        //- Construct from patch and internal field
        wallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- wallFunctionFvPatchScalarField
        //- onto a new patch
        wallFunctionFvPatchScalarField
        (
            const wallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallFunctionFvPatchScalarField
        (
            const wallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        wallFunctionFvPatchScalarField
        (
            const wallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        // No clone methods - abstract class


    // Member Functions

        // Access

            //- Return the empirical model coefficient, Cmu
            scalar Cmu() const noexcept
            {
                return Cmu_;
            }

            //- Return the von Karman constant
            scalar kappa() const noexcept
            {
                return kappa_;
            }

            //- Return the wall-roughness parameter, E
            scalar E() const noexcept
            {
                return E_;
            }

            //- Return the estimated y+ at the two-sublayer intersection
            scalar yPlusLam() const noexcept
            {
                return yPlusLam_;
            }


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // IO

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
