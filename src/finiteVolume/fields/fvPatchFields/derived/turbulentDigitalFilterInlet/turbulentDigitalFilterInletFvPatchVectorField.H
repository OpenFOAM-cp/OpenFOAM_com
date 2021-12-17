/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::turbulentDigitalFilterFvPatchVectorField

Group
    grpInletBoundaryConditions

Description

Usage

Note

See also
    - turbulentDFSEMInletFvPatchVectorField.C

SourceFiles
    turbulentDigitalFilterInletFvPatchVectorField.C
    turbulentDigitalFilterInletFvPatchVectorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDigitalFilterInletFvPatchVectorField_H
#define turbulentDigitalFilterInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fieldTypes.H"
#include "PatchFunction1.H"
#include "Enum.H"
#include "randomNumberBox.H"
#include "AMIPatchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class turbulentDigitalFilterInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDigitalFilterInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Enumerations

        //- Options for the correlation function form
        enum kernelType : bool
        {
            GAUSSIAN = true,     //!< "Standard Gaussian correlation function"
            EXPONENTIAL = false  //!< "Exponential correlation function"
        };

        //- Names for kernelType
        static const Enum<kernelType> kernelTypeNames;


    // Private Data

        //- Mean inlet velocity profile field in global coordinates [m/s]
        autoPtr<PatchFunction1<vector>> UMean_;

        //- Reynolds stress tensor profile field in global coordinates [m2/s2]
        autoPtr<PatchFunction1<symmTensor>> R_;


        // Random-number box information

            //- Random-number box
            syntheticTurbulence::randomNumberBox box_;

            //- Pointer to AMI interpolator
            mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;


        // Kernel information

            //- Correlation function form of the kernel
            const enum kernelType kernelType_;

            //- Kernel coefficients corresponding to L [-]
            const List<List<scalar>> kernel_;


        //- Patch normal into the domain
        vector patchNormal_;

        //- Current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Initialise patch information
        void initialisePatch();

        //- Calculate the kernel once per simulation
        List<List<scalar>> calcKernel() const;

        //- Lund transformation
        void transform(symmTensorField& R) const;


public:

   //- Runtime type information
   TypeName("turbulentDigitalFilterInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentDigitalFilterInletFvPatchVectorField onto a new patch
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDigitalFilterInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDigitalFilterInletFvPatchVectorField
        (
            const turbulentDigitalFilterInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDigitalFilterInletFvPatchVectorField(*this, iF)
            );
        }


    //- Destructor
    virtual ~turbulentDigitalFilterInletFvPatchVectorField() = default;


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField& ptf,
                const labelList& addr
            );


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write boundary condition settings
            virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
