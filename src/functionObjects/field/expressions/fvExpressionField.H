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
    Foam::functionObjects::fvExpressionField

Group
    grpFieldFunctionObjects

Description
    Function object that generates or modifies a field based on expressions.

Usage
    A minimal example:
    \verbatim
    newField
    {
        type    exprField;
        libs    (fieldFunctionObjects);

        operation      define;  // Default
        fieldName      rhoUT;

        expression  "rho * U * T";
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type   | Reqd | Dflt
      type         | Type name: exprField               | word   | yes  |
      libs         | Libraries: fieldFunctionObjects    | wordList | yes |
      operation    | Operation type: see below          | word   | no  | define
      fieldName    | Names of the output field          | word   | yes |
      readFields   | Preload named fields (post-process mode) | wordList | no |
      expression   | Field evaluation expression        | string | yes |
      fieldMask    | Masking as logical expression      | string | no  | ""
      autowrite    | Add AUTO_WRITE to created field    | string | yes |
      dimensions   | Use specified dimensions for created field  | dim-spec | no |
      useNamePrefix | Add prefix scoping to output name | bool   | no  | false
    \endtable

    Options for the \c operation entry:
    \plaintable
      none          | no operation
      calc          | simple calculation, without storing the result
      define        | define field based on expression
      modify        | adjust field according to expression and fieldMask
    \endplaintable

Note
    The \c scopedName entry is ignored for the \c modify operation.

SourceFiles
    fvExpressionField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fvExpressionField_H
#define functionObjects_fvExpressionField_H

#include "fvMeshFunctionObject.H"
#include "volumeExprDriver.H"
#include "Enum.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class fvExpressionField Declaration
\*---------------------------------------------------------------------------*/

class fvExpressionField
:
    public fvMeshFunctionObject
{
public:

    // Public Data Types

        //- Operation type enumeration
        enum operationType : unsigned char
        {
            opNone = 0,     //!< No operation
            opCalculate,    //!< Calculate and write (no storage)
            opDefine,       //!< Define field and store (default)
            opModify        //!< Modify existing
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


protected:

    // Private Data

        //- The context dictionary
        dictionary dict_;

        //- Name of the field
        word fieldName_;

        //- Names fields to preload
        wordList preloadFields_;

        //- The field-mask expression (modify mode)
        expressions::exprString maskExpr_;

        //- Expression to evaluate
        expressions::exprString valueExpr_;

        //- Dimensions for new field
        dimensionSet dimensions_;

        //- Operation mode
        operationType operation_;

        //- Set AUTO_WRITE for new field (future?)
        bool autowrite_;

        //- True if dimensions_ should be used (creation)
        bool hasDimensions_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        autoPtr<expressions::volumeExprDriver> driver_;


    // Private Member Functions

        //- Attempt load from io, store on database if successful
        template<class FieldType>
        bool loadAndStore(const IOobject& io);

        //- Forward to loadAndStore for supported types
        template<class Type>
        bool loadField(const IOobject& io);

        //- Attempt to load specified fields
        label loadFields(const UList<word>& fieldSet_);

        template<class GeoField>
        bool setField
        (
            GeoField& output,
            const GeoField& evaluated,
            const boolField& cond
        );

        bool performAction(bool doWrite);

public:

    //- Runtime type information
    TypeName("exprField");


    // Constructors

        //- Construct from Time and dictionary
        fvExpressionField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false
        );

        //- No copy construct
        fvExpressionField(const fvExpressionField&) = delete;

        //- No copy assignment
        void operator=(const fvExpressionField&) = delete;


    //- Destructor
    virtual ~fvExpressionField();


    // Member Functions

        //- Qualified/unqualified field name (depends on operation type)
        virtual word fieldName() const;

        //- Read the data
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
