/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "nutWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nutWallFunctionFvPatchScalarField, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::volVectorField& Foam::nutWallFunctionFvPatchScalarField::U
(
    const turbulenceModel& turb
) const
{
    if (UName_ == word::null)
    {
        return turb.U();
    }
    else
    {
        return db().lookupObject<volVectorField>(UName_);
    }
}


void Foam::nutWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    if (UName_ != word::null)
    {
        os.writeEntry("U", UName_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallFunctionFvPatchScalarField(p, iF),
    UName_(word::null)
{}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    wallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_)
{}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    wallFunctionFvPatchScalarField(p, iF, dict),
    UName_(dict.getOrDefault<word>("U", word::null))
{}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& wfpsf
)
:
    wallFunctionFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_)
{}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    wallFunctionFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::nutWallFunctionFvPatchScalarField::blend
(
    const scalar nutVis,
    const scalar nutLog,
    const scalar yPlus
) const
{
    scalar nutw = 0.0;

    switch (blending_)
    {
        case blendingType::STEPWISE:
        {
            if (yPlus > this->yPlusLam())
            {
                nutw = nutLog;
            }
            else
            {
                nutw = nutVis;
            }
            break;
        }

        case blendingType::MAX:
        {
            // (PH:Eq. 27)
            nutw = max(nutVis, nutLog);
            break;
        }

        case blendingType::BINOMIAL:
        {
            // (ME:Eqs. 15-16)
            nutw =
                pow
                (
                    pow(nutVis, n_) + pow(nutLog, n_),
                    1.0/n_
                );
            break;
        }

        case blendingType::EXPONENTIAL:
        {
            // (PH:Eq. 31)
            const scalar Gamma = 0.01*pow4(yPlus)/(1.0 + 5.0*yPlus);
            const scalar invGamma = 1.0/(Gamma + ROOTVSMALL);

            nutw = nutVis*exp(-Gamma) + nutLog*exp(-invGamma);
            break;
        }

        case blendingType::TANH:
        {
            break;
        }
    }

    return nutw;
}


void Foam::nutWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(calcNut());

    wallFunctionFvPatchScalarField::updateCoeffs();
}


void Foam::nutWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// ************************************************************************* //
