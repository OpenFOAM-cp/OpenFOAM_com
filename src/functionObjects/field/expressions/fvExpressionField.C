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

\*---------------------------------------------------------------------------*/

#include "fvExpressionField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "volumeExprDriver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fvExpressionField, 0);
    addToRunTimeSelectionTable(functionObject, fvExpressionField, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::fvExpressionField::operationType
>
Foam::functionObjects::fvExpressionField::operationTypeNames_
({
    { operationType::opNone, "none" },
    { operationType::opCalculate, "calc" },
    { operationType::opDefine, "define" },
    { operationType::opModify, "modify" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

word fieldGeoType(const expressions::FieldAssociation geoType)
{
    switch (geoType)
    {
        case expressions::FieldAssociation::POINT_DATA : return "points"; break;
        case expressions::FieldAssociation::FACE_DATA : return "faces"; break;
        case expressions::FieldAssociation::VOLUME_DATA : return "cells"; break;
        default: break;
    }
    return "unknown";
}


template<class Type>
static void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    if (correctBCs)
    {
        // Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, pointPatchField, pointMesh>& field
)
{
    if (correctBCs)
    {
        // Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FieldType>
bool Foam::functionObjects::fvExpressionField::loadAndStore(const IOobject& io)
{
    if (FieldType::typeName == io.headerClassName())
    {
        // Store field on mesh database
        Log << "    Reading " << io.name()
            << " (" << FieldType::typeName << ')' << endl;

        mesh_.objectRegistry::store(new FieldType(io, mesh_));
        return true;
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::fvExpressionField::loadField(const IOobject& io)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    // typedef typename VolFieldType::Internal IntVolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    return
    (
        loadAndStore<VolFieldType>(io)
        /// || loadAndStore<IntVolFieldType>(io)
     || loadAndStore<SurfaceFieldType>(io)
    );
}


Foam::label Foam::functionObjects::fvExpressionField::loadFields
(
    const UList<word>& fieldSet_
)
{
    label nLoaded = 0;

    for (const word& fieldName : fieldSet_)
    {
        // Already loaded?
        const auto* ptr = mesh_.cfindObject<regIOobject>(fieldName);

        if (ptr)
        {
            ++nLoaded;
            DebugInfo
                << "readFields : "
                << ptr->name() << " (" << ptr->type()
                << ") already in database" << endl;
            continue;
        }

        // Load field as necessary
        IOobject io
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        const bool ok =
        (
            io.typeHeaderOk<regIOobject>(false) // Preload header info
         && !io.headerClassName().empty()       // Extra safety
         &&
            (
                loadField<scalar>(io)
             || loadField<vector>(io)
             || loadField<sphericalTensor>(io)
             || loadField<symmTensor>(io)
             || loadField<tensor>(io)
            )
        );

        if (ok)
        {
            ++nLoaded;
        }
        else
        {
            DebugInfo
                << "readFields : failed to load " << fieldName << endl;
        }
    }

    return nLoaded;
}


template<class GeoField>
bool Foam::functionObjects::fvExpressionField::setField
(
    GeoField& output,
    const GeoField& evaluated,
    const boolField& fieldMask
)
{
    label numValuesChanged = 0;

    if (fieldMask.empty())
    {
        // No field-mask - set all
        numValuesChanged = output.size();

        const bool oldDimChecking = dimensionSet::checking(false);

        output = evaluated;

        dimensionSet::checking(oldDimChecking);
    }
    else
    {
        forAll(output, idx)
        {
            if (fieldMask[idx])
            {
                output[idx] = evaluated[idx];
                ++numValuesChanged;
            }
        }
    }

    forAll(evaluated.boundaryField(), patchi)
    {
        auto& pf = output.boundaryFieldRef()[patchi];

        if (pf.patch().coupled())
        {
            pf == evaluated.boundaryField()[patchi];
        }
    }

    doCorrectBoundaryConditions(true, output);

    const label numTotal = returnReduce(output.size(), plusOp<label>());
    reduce(numValuesChanged, plusOp<label>());

    if (numValuesChanged == numTotal)
    {
        Log<< "Set all ";
    }
    else
    {
        Log<< "Set " << numValuesChanged << " of ";
    }
    Log<< numTotal << " values" << endl;

    if (hasDimensions_)
    {
        // Info<< "Setting dimensions to " << dims << endl;
        output.dimensions().reset(dimensions_);
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fvExpressionField::fvExpressionField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fvMeshFunctionObject(name, runTime, dict),
    dict_(dict),  // Deep copy
    fieldName_(),
    preloadFields_(),
    maskExpr_(),
    valueExpr_(),
    dimensions_(),
    operation_(operationType::opDefine),
    autowrite_(false),
    hasDimensions_(false),
    loadFromFiles_(loadFromFiles)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fvExpressionField::~fvExpressionField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::functionObjects::fvExpressionField::fieldName() const
{
    switch (operation_)
    {
        case opNone:
        {
            break;  // No-op
        }
        case opCalculate:
        case opDefine:
        {
            return scopedName(fieldName_);
        }
        case opModify:
        {
            return fieldName_;
        }
    }

    return word::null;
}


bool Foam::functionObjects::fvExpressionField::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    operation_ = operationTypeNames_.getOrDefault
    (
        "operation",
        dict,
        operationType::opDefine
    );

    fieldName_ = dict.get<word>("fieldName");
    const word fldName = fieldName();

    Log << type() << ' ' << this->name() << ':' << nl
        << "    operation     = " << operationTypeNames_[operation_] << nl
        << "    field         = " << fldName << nl;

    maskExpr_.clear();
    valueExpr_.clear();

    preloadFields_.clear();
    dict.readIfPresent("readFields", preloadFields_);

    switch (operation_)
    {
        case opNone:
        {
            break;  // No-op
        }
        case opCalculate:
        case opDefine:
        case opModify:
        {
            if (operation_ == opModify)
            {
                // fieldMask is optional for modify
                maskExpr_.readEntry("fieldMask", dict, false);
            }

            valueExpr_.readEntry("expression", dict);  // mandatory
            break;
        }
    }

    autowrite_ = dict.getOrDefault<bool>("autowrite", true);

    // "dimensions" is optional
    dimensions_.clear();
    hasDimensions_ = dimensions_.readEntry("dimensions", dict, false);

    if (operation_ == opCalculate || operation_ == opDefine)
    {
        if (!hasDimensions_)
        {
            Log << "    no 'dimensions' : treat '" << fldName
                << "' as dimensionless" << endl;
        }
    }
    else
    {
        // Ignore for none/modify
        hasDimensions_  = false;
    }


    if (operation_ == opNone)
    {
        driver_.reset(nullptr);
        return true;  // Done
    }

    driver_.reset
    (
        new expressions::volumeExprDriver(mesh_, dict_)
    );

    int ctrl
    (
        expressions::exprDriver::SEARCH_REGISTRY
    );

    if (loadFromFiles_)
    {
        ctrl |= expressions::exprDriver::SEARCH_FILES;
    }

    driver_->setSearchBehaviour
    (
        expressions::exprDriver::searchControls(ctrl),
        false  // No caching
    );

    driver_->readDict(dict_);

    return true;
}


bool Foam::functionObjects::fvExpressionField::performAction(bool doWrite)
{
    using FieldAssociation = expressions::FieldAssociation;

    if (!driver_ || operation_ == opNone)
    {
        // No-op
        return true;
    }

    const word fldName = fieldName();

    if (loadFromFiles_)
    {
        loadFields(preloadFields_);
    }

    if (operation_ == opModify && loadFromFiles_)
    {
        loadFields(wordList({fldName}));
    }

    auto& driver = *driver_;


    // Current availability
    auto* regIOobjectPtr = mesh_.getObjectPtr<regIOobject>(fldName);

    if (!regIOobjectPtr && operation_ == opModify)
    {
        // Cannot continue
        FatalErrorInFunction
            << type() << ' ' << this->name() << ':' << nl
            << "    missing-field: " << fldName << nl
            << exit(FatalError);

        return false;
    }


    // Handle "field-mask" evaluation
    bool evalFieldMask =
    (
        (operation_ == opModify)
     && maskExpr_.size() && maskExpr_ != "true" && maskExpr_ != "1"
    );

    boolField fieldMask;
    FieldAssociation maskFieldAssoc(FieldAssociation::NO_DATA);

    if (evalFieldMask)
    {
        driver.parse(maskExpr_);

        if (driver.isLogical())
        {
            auto& result = driver.result();
            if (result.is_bool())
            {
                fieldMask = result.getResult<bool>();
                maskFieldAssoc = driver.fieldAssociation();
            }
        }

        // Slightly pedantic...
        driver.clearResult();
        driver.clearField();

        evalFieldMask = (maskFieldAssoc != FieldAssociation::NO_DATA);

        if (!evalFieldMask)
        {
            FatalErrorInFunction
                << "field-mask: " << maskExpr_
                << " does not evaluate to a logical expression: "
                << driver.resultType() << nl
                #ifdef FULLDEBUG
                << "contents: " << fieldMask
                #endif
                << exit(FatalError);
        }
    }


    // Start "expression" evaluation

    bool applied = false;
    autoPtr<regIOobject> toutputField;

    {
        driver.clearVariables();

        driver.parse(valueExpr_);

        if (evalFieldMask && maskFieldAssoc != driver.fieldAssociation())
        {
            FatalErrorInFunction
                << "Mismatch between field-mask geometric type ("
                << fieldGeoType(maskFieldAssoc) << ") and" << nl
                << "expression geometric type ("
                << fieldGeoType(driver.fieldAssociation()) << ')' << nl
                << nl
                << "Expression: " << valueExpr_ << nl
                << "Field-mask: " << maskExpr_ << nl
                << nl
                << exit(FatalError);
        }

        // The output field does not appear to exist
        // - create a new 'blank slate'
        if (!regIOobjectPtr)
        {
            toutputField.reset(driver.dupZeroField());

            if (toutputField)
            {
                toutputField->rename(fldName);

                if (autowrite_)
                {
                    toutputField->writeOpt(IOobject::AUTO_WRITE);
                }
            }

            if (operation_ == opCalculate)
            {
                // Local (non-registered) storage. Discard afterwards
                regIOobjectPtr = toutputField.get();
            }
            else
            {
                if (toutputField->checkIn() && toutputField->store())
                {
                    // Register and transfer ownership to registry
                    toutputField.release();
                }

                regIOobjectPtr = mesh_.getObjectPtr<regIOobject>(fldName);
            }
        }


        // Addiitional checks (TBD):

        if (!regIOobjectPtr)
        {
            // Cannot continue
            FatalErrorInFunction
                << type() << ' ' << this->name() << ':' << nl
                << "    missing-field: " << fldName << nl
                << exit(FatalError);
        }

        // const word oldFieldType = regIOobject->type();

        // if (!ctrl.createNew && driver.resultType() != oldFieldType)
        // {
        //     FatalErrorInFunction
        //         << "Inconsistent types: " << fieldName << " is  "
        //         << oldFieldType
        //         << " but the expression evaluates to "
        //         << driver.resultType()
        //         << exit(FatalError);
        //     }

        switch (driver.fieldAssociation())
        {
            #undef  doLocalCode
            #define doLocalCode(GeoField)                                    \
            {                                                                \
                /* FieldType */                                              \
                auto* outPtr = dynamic_cast<GeoField*>(regIOobjectPtr);      \
                const auto* ptr = driver.isResultType<GeoField>();           \
                                                                             \
                if (outPtr && ptr)                                           \
                {                                                            \
                    applied = setField                                       \
                    (                                                        \
                        *outPtr,                                             \
                        *ptr,                                                \
                        fieldMask                                            \
                    );                                                       \
                    if (doWrite)                                             \
                    {                                                        \
                        outPtr->write();                                     \
                    }                                                        \
                    break;                                                   \
                }                                                            \
            }

            case FieldAssociation::VOLUME_DATA:
            {
                doLocalCode(volScalarField);
                doLocalCode(volVectorField);
                doLocalCode(volTensorField);
                doLocalCode(volSymmTensorField);
                doLocalCode(volSphericalTensorField);
                break;
            }
            case FieldAssociation::FACE_DATA:
            {
                doLocalCode(surfaceScalarField);
                doLocalCode(surfaceVectorField);
                doLocalCode(surfaceTensorField);
                doLocalCode(surfaceSymmTensorField);
                doLocalCode(surfaceSphericalTensorField);
                break;
            }
            case FieldAssociation::POINT_DATA:
            {
                doLocalCode(pointScalarField);
                doLocalCode(pointVectorField);
                doLocalCode(pointTensorField);
                doLocalCode(pointSymmTensorField);
                doLocalCode(pointSphericalTensorField);
                break;
            }

            default: break;
            #undef doLocalCode
        }
    }


    // Clear out heavier data
    driver.clearResult();
    driver.clearField();

    if (!applied)
    {
        // Or error?
        WarningInFunction
            << type() << ' ' << this->name() << ": Failed to apply "
            << operationTypeNames_[operation_] << " for " << fldName
            << nl;
    }

    return true;
}


bool Foam::functionObjects::fvExpressionField::execute()
{
    return performAction(false);
}


bool Foam::functionObjects::fvExpressionField::write()
{
    return performAction(true);
}


// ************************************************************************* //
