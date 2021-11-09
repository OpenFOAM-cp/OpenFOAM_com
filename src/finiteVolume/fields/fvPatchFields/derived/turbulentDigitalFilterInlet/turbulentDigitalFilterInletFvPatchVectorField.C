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

\*---------------------------------------------------------------------------*/

#include "turbulentDigitalFilterInletFvPatchVectorField.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaWeightAMI.H"
#include "turbulentDFSEMInletFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::turbulentDigitalFilterInletFvPatchVectorField::kernelType
>
Foam::turbulentDigitalFilterInletFvPatchVectorField::kernelTypeNames
({
    { kernelType::GAUSSIAN, "Gaussian" },
    { kernelType::EXPONENTIAL , "exponential" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentDigitalFilterInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    const scalar error = max(magSqr(patchNormal_ + nf));

    if (error > SMALL)
    {
        WarningInFunction
            << "Patch " << patch().name() << " is not planar"
            << endl;
    }

    patchNormal_ /= mag(patchNormal_) + ROOTVSMALL;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::calcKernel() const
{
    List<List<scalar>> kernel
    (
        pTraits<tensor>::nComponents,
        List<scalar>(1, 1.0)
    );

    const tensor L(box_.L());

    for (direction dir = 0; dir < pTraits<tensor>::nComponents; ++dir)
    {
        // The smallest filter width larger than length scale
        // (KSJ:'n' in Eq. 15)
        const label n = ceil(L[dir]);

        // Full filter-width (except mid-zero) according to the condition
        // (KSJ:Eq. 15 whereat N is minimum =2n)
        const label twiceN = 4*n;

        // Initialise filter-coeffs containers with full filter-width size
        // Extra elem is mid-zero within [-N, N]
        kernel[dir] = List<scalar>(twiceN + 1, Zero);

        // First element: -N within [-N, N]
        const scalar initElem = -2*scalar(n);

        // Container initialised with [-N, N] (KSJ:p. 658, item-b)
        std::iota
        (
            kernel[dir].begin(),
            kernel[dir].end(),
            initElem
        );

        // Compute filter-coeffs (KSJ:Eq. 14 (Gaussian))
        List<scalar> fTemp(kernel[dir]);
        scalar fSum = 0;
        const scalar nSqr = n*n;

        // Model constant shaping the autocorrelation function (KSJ:Eq. 14)
        const scalar C = -0.5*constant::mathematical::pi;

        if (kernelType_)
        {
            fTemp = sqr(exp(C*sqr(fTemp)/nSqr));
            fSum = Foam::sqrt(sum(fTemp));

            kernel[dir] =
                exp(C*sqr(kernel[dir])/nSqr)/fSum;
        }
        else
        {
            fTemp = sqr(exp(C*mag(fTemp)/n));
            fSum = Foam::sqrt(sum(fTemp));

            kernel[dir] = exp(C*mag(kernel[dir])/n)/fSum;
        }
    }

    return kernel;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::transform
(
    symmTensorField& R
) const
{
    #ifdef FULLDEBUG
    turbulentDFSEMInletFvPatchVectorField::checkStresses(R);
    #endif

    forAll(R, i)
    {
        symmTensor& r = R[i];

        // (KSJ:Eq. 5)
        r.xx() = Foam::sqrt(r.xx());
        r.xy() /= r.xx();
        r.xz() /= r.xx();
        r.yy() = Foam::sqrt(r.yy() - sqr(r.xy()));
        r.yz() = (r.yz() - r.xy()*r.xz())/r.yy();
        r.zz() = Foam::sqrt(r.zz() - sqr(r.xz()) - sqr(r.yz()));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    UMean_(nullptr),
    R_(nullptr),
    box_(),
    AMIPtr_(new faceAreaWeightAMI(true, false)),
    kernelType_(kernelType::GAUSSIAN),
    kernel_(Zero),
    curTimeIndex_(-1)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    UMean_(ptf.UMean_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    box_(ptf.box_),
    AMIPtr_(ptf.AMIPtr_->clone()),
    kernelType_(ptf.kernelType_),
    kernel_(ptf.kernel_),
    curTimeIndex_(-1)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    UMean_(PatchFunction1<vector>::New(patch().patch(), "UMean", dict)),
    R_(PatchFunction1<symmTensor>::New(patch().patch(), "R", dict)),
    box_(p, dict),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            true // flipNormals
        )
    ),
    kernelType_
    (
        kernelTypeNames.getOrDefault
        (
            "kernelType",
            dict,
            kernelType::GAUSSIAN
        )
    ),
    kernel_(calcKernel()),
    curTimeIndex_(-1)
{
    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t));

    turbulentDFSEMInletFvPatchVectorField::checkStresses(R);

    // Check if varying or fixed time-step computation
    if (db().time().isAdjustTimeStep())
    {
        WarningInFunction
            << "  # Varying time-step computations are not fully supported #"
            << endl;
    }

    AMIPtr_->calculate(p.patch(), box_.patch());
}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    UMean_(ptf.UMean_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    box_(ptf.box_),
    AMIPtr_(ptf.AMIPtr_->clone()),
    kernelType_(ptf.kernelType_),
    kernel_(ptf.kernel_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    UMean_(ptf.UMean_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    box_(ptf.box_),
    AMIPtr_(ptf.AMIPtr_->clone()),
    kernelType_(ptf.kernelType_),
    kernel_(ptf.kernel_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDigitalFilterInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);

    if (UMean_)
    {
        UMean_->autoMap(m);
    }
    if (R_)
    {
        R_->autoMap(m);
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const auto& dfsemptf =
        refCast<const turbulentDigitalFilterInletFvPatchVectorField>(ptf);

    if (UMean_)
    {
        UMean_->rmap(dfsemptf.UMean_(), addr);
    }
    if (R_)
    {
        R_->rmap(dfsemptf.R_(), addr);
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        vectorField& U = *this;
        /*if (U.size())
        {
            U = Zero;
        }*/
        U = Zero;
        vectorField virtualField;

        if (Pstream::master())
        {
            // Embed two-point correlations
            virtualField = box_.convolve(kernel_);

            box_.shift();

            box_.refill();
        }

        plusEqOp<vector> cop;

        AMIPtr_->interpolateToSource
        (
            virtualField,
            multiplyWeightedOp<vector, plusEqOp<vector>>(cop),
            U,
            UList<vector>::null()
        );

        const scalar t = db().time().timeOutputValue();
        tmp<symmTensorField> tR = R_->value(t);

        // Apply Lund transformation
        transform(tR.ref());

        // Embed one-point correlations (KSJ:p. 658, item-e)
        U = U & tR;

        tmp<vectorField> tUMean = UMean_->value(t);

        // (PCR:p. 522)
        const vector UBulk
        (
            gSum(tUMean.cref()*patch().magSf())
           /(gSum(patch().magSf()) + ROOTVSMALL)
        );

        U += tUMean;

        // (KCX:Eq. 8)
        // Re-scale to ensure correct flow rate
        const scalar fCorr =
            gSum((UBulk & patchNormal_)*patch().magSf())
           /gSum(U & -patch().Sf());

        U *= fCorr;

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);

    if (UMean_)
    {
        UMean_->writeData(os);
    }
    if (R_)
    {
        R_->writeData(os);
    }
    os.writeEntry("kernelType", kernelTypeNames[kernelType_]);
    os.writeEntry("AMIMethod", AMIPtr_->type());
    box_.write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentDigitalFilterInletFvPatchVectorField
   );
}


// ************************************************************************* //
