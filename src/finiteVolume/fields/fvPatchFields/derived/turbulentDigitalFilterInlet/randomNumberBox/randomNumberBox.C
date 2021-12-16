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

#include "randomNumberBox.H"
#include "Tuple2.H"
#include "vector.H"
#include "cartesianCS.H"
#include "SubList.H"

#include "OBJstream.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::syntheticTurbulence::randomNumberBox::debug = 0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::syntheticTurbulence::randomNumberBox::initialiseBoundBox()
{
    // Get patch normal direction into the domain
    const vector nf((-gAverage(p_.nf())).normalise());

    // Find the second local coordinate direction
    direction minCmpt = 0;
    scalar minMag = mag(nf[minCmpt]);
    for (direction cmpt = 1; cmpt < pTraits<vector>::nComponents; ++cmpt)
    {
        const scalar s = mag(nf[cmpt]);
        if (s < minMag)
        {
            minMag = s;
            minCmpt = cmpt;
        }
    }

    // Create the second local coordinate direction
    vector e2(Zero);
    e2[minCmpt] = 1;

    // Remove normal component
    e2 -= (nf&e2)*nf;

    // Create the local coordinate system
    coordSystem::cartesian cs
    (
        Zero,   // origin
        nf,     // normal
        e2      // 0-axis
    );

    // Convert patch points into local coordinate system
    const pointField localPos
    (
        cs.localPosition
        (
            pointField
            (
                p_.patch().points(),
                p_.patch().meshPoints()
            )
        )
    );

    const bool globalReduction = true;
    const boundBox bb(localPos, globalReduction);

    bbSpan_ = bb.span();

    return bb.min();
}


Foam::Pair<Foam::scalar>
Foam::syntheticTurbulence::randomNumberBox::initialiseDelta() const
{
    return Pair<scalar>(bbSpan_[0]/n_.first(), bbSpan_[1]/n_.second());
}


Foam::List<Foam::label>
Foam::syntheticTurbulence::randomNumberBox::initialiseRanges() const
{
    if (!Pstream::master())
    {
        return labelList();
    }

    List<label> ranges(pTraits<tensor>::nComponents, 0);
    const Vector<label> slice(1, n_.first(), n_.second());

    const tensor L(meterToCell(L_));

    for (label i = 0; i < pTraits<tensor>::nComponents; ++i)
    {
        const label slicei = label(i/pTraits<vector>::nComponents);

        const label n = ceil(L[i]);
        const label twiceN = 4*n;

        ranges[i] = slice[slicei] + twiceN;
    }

    return ranges;
}


Foam::List<Foam::label>
Foam::syntheticTurbulence::randomNumberBox::initialiseSliceRanges() const
{
    if (!Pstream::master())
    {
        return labelList();
    }

    List<label> sliceRanges(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        sliceRanges[dir] =
            ranges_[pTraits<vector>::nComponents + dir]
           *ranges_[pTraits<symmTensor>::nComponents + dir];
    }

    return sliceRanges;
}


Foam::List<Foam::label>
Foam::syntheticTurbulence::randomNumberBox::initialiseComponentRanges() const
{
    if (!Pstream::master())
    {
        return labelList();
    }

    List<label> componentRanges(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        componentRanges[dir] = sliceRanges_[dir]*ranges_[dir];
    }

    return componentRanges;
}


Foam::List<Foam::label>
Foam::syntheticTurbulence::randomNumberBox::initialiseLastSlice() const
{
    if (!Pstream::master())
    {
        return labelList();
    }

    List<label> lastSlice(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        lastSlice[dir] = sliceRanges_[dir]*(ranges_[dir] - 1);
    }

    return lastSlice;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::syntheticTurbulence::randomNumberBox::initialiseBox()
{
    if (!Pstream::master())
    {
        return scalarListList();
    }

    List<List<scalar>> box(pTraits<vector>::nComponents, List<scalar>());

        // Initialise: Remaining convenience factors for (e1 e2 e3)
        for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
        {
            List<scalar>& randomSet = box[dir];

            randomSet = List<scalar>(componentRanges_[dir]);

            // Initialise: Random-box content with
            // random-number sets, generate random-number sets
            // obeying the standard normal distribution
            std::generate
            (
                randomSet.begin(),
                randomSet.end(),
                [&]{ return rndGen_.GaussNormal<scalar>(); }
            );
        }

    return box;
}


Foam::Field<Foam::point>
Foam::syntheticTurbulence::randomNumberBox::initialisePatchPoints() const
{
    if (!Pstream::master())
    {
        return pointField();
    }

    // List of vertex points of the virtual patch in the local coordinate system
    const label nx = n_.first();
    const label ny = n_.second();
    const label nPoints = (nx + 1)*(ny + 1);
    pointField points(nPoints, Zero);

    label pointi = 0;
    for (label i = 0; i <= nx; ++i) // height
    {
        for (label j = 0; j <= ny; ++j) // width
        {
            const point p
            (
                Zero, // TODO: third direction
                bbMin_[0] + i*delta_.first(),
                bbMin_[1] + j*delta_.second()
            );
            points[pointi] = p;
            ++pointi;
        }
    }

    // Projection of virtual patch points onto patch plane
    /* const point& pointOnPatch = p_.patch().points()[0];
    for (point& p : points)
    {
        p -= (nf & (p - pointOnPatch))*nf;
    } // TODO: WRONG - look for inverse transformation*/

    return points;
}


Foam::List<Foam::face>
Foam::syntheticTurbulence::randomNumberBox::initialisePatchFaces() const
{
    if (!Pstream::master())
    {
        return List<face>();
    }

     // List of faces of the virtual patch
    const label nx = n_.first();
    const label ny = n_.second();
    const label nFaces = nx*ny;
    List<face> faces(nFaces);

    label m = 0;
    for (label i = 0; i < nx; ++i)
    {
        for (label j = 0; j < ny; ++j)
        {
            const label k = i*(ny+1) + j;
            faces[m] = face({k, k+1, k+(ny+2), k+(ny+1)});
            ++m;
        }
    }

    return faces;
}


Foam::primitivePatch
Foam::syntheticTurbulence::randomNumberBox::initialisePatch() const
{
    //if (debug)
    {
        const auto& tm = p_.patch().boundaryMesh().mesh().time();
        OBJstream os(tm.path()/"patch.obj");
        os.write(patchFaces_, patchPoints_, false);
    }

    return
        primitivePatch
        (
            SubList<face>(patchFaces_, patchFaces_.size()),
            patchPoints_
        );
}


Foam::tensor Foam::syntheticTurbulence::randomNumberBox::meterToCell
(
    const tensor& L
) const
{
    tensor Ls(L);

    const scalar deltaT =
        p_.patch().boundaryMesh().mesh().time().deltaTValue();

    //  (KSJ:Eq. 13)
    // Integral time scales
    Ls.row(0, Ls.x()/deltaT);
    // Integral length scales
    Ls.row(1, Ls.y()/delta_.first());
    Ls.row(2, Ls.z()/delta_.second());

    return Ls;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::syntheticTurbulence::randomNumberBox::randomNumberBox
(
    const fvPatch& p,
    const dictionary& dict
)
:
    p_(p),
    bbSpan_(Zero),
    bbMin_(initialiseBoundBox()),
    fixSeed_(dict.getOrDefault("fixSeed", true)),
    continuous_(dict.getOrDefault("continuous", false)),
    L_
    (
        dict.getCheck<tensor>
        (
            "L",
            [=](const tensor x)
            {
                return cmptMin(x) > ROOTSMALL ? true : false;
            }
        )
    ),
    n_
    (
        dict.getCheck<Pair<label>>
        (
            "n",
            [=](const Pair<label> x)
            {
                return min(x.first(), x.second()) > 0 ? true : false;
            }
        )
    ),
    delta_(initialiseDelta()),
    ranges_(initialiseRanges()),
    sliceRanges_(initialiseSliceRanges()),
    componentRanges_(initialiseComponentRanges()),
    lastSlice_(initialiseLastSlice()),
    rndGen_(fixSeed_ ? 1234 : time(0)),
    box_
    (
        continuous_
      ? dict.getOrDefault<List<List<scalar>>>
        (
            "box",
            initialiseBox() // First time-step
        )
      :
        initialiseBox()
    ),
    patchPoints_(initialisePatchPoints()),
    patchFaces_(initialisePatchFaces()),
    patch_(initialisePatch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::syntheticTurbulence::randomNumberBox::shift()
{
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        List<scalar>& slice = box_[dir];

        // Shift forward from the back to the front / First Out
        inplaceRotateList(slice, sliceRanges_[dir]);
    }
}


void Foam::syntheticTurbulence::randomNumberBox::refill()
{
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        List<scalar>& slice = box_[dir];

        // Refill the back with a new random-number set / First In
        for (label i = 0; i < sliceRanges_[dir]; ++i)
        {
            slice[i] = rndGen_.GaussNormal<scalar>();
        }
    }
}


Foam::Field<Foam::vector>
Foam::syntheticTurbulence::randomNumberBox::convolve
(
    const List<List<scalar>>& kernel
) const
{
    Field<vector> fld(n_.first()*n_.second(), Zero);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        const List<scalar>& in = box_[dir];
        Field<scalar> out(n_.first()*n_.second(), Zero);

        const List<scalar>& filter1 = kernel[dir];
        const List<scalar>& filter2 = kernel[3 + dir];
        const List<scalar>& filter3 = kernel[6 + dir];

        const label sz1 = ranges_[dir];
        const label sz2 = ranges_[3 + dir];
        const label sz3 = ranges_[6 + dir];
        const label szfilter1 = kernel[dir].size();
        const label szfilter2 = kernel[3 + dir].size();
        const label szfilter3 = kernel[6 + dir].size();
        const label sz23 = sliceRanges_[dir];
        const label sz123 = componentRanges_[dir];
        const label validSlice2 = sz2 - (szfilter2 - 1);
        const label validSlice3 = sz3 - (szfilter3 - 1);

        // Convolution summation - Along 1st direction
        scalarField tmp(sz123);
        label filterCentre = label(szfilter2/label(2));
        label endIndex = sz2 - filterCentre;
        label i0 = 0;
        label i1 = 0;
        label i2 = 0;

        for (label i = 0; i < sz1; ++i)
        {
            for (label j = 0; j < sz3; ++j)
            {
                i1 += filterCentre;

                for (label k = filterCentre; k < endIndex; ++k)
                {
                    tmp[i1] = 0.0;
                    label q = 0;

                    for (label p = szfilter2 - 1; p >= 0; --p, ++q)
                    {
                        tmp[i1] += in[i0 + q]*filter2[p];
                    }
                    ++i0;
                    ++i1;
                }
                i0 += (filterCentre + filterCentre);
                i1 += filterCentre;
            }
        }

        // Convolution summation - Along 2nd direction
        scalarField tmp2(tmp);
        filterCentre = label(szfilter3/label(2));
        endIndex = sz3 - filterCentre;
        i1 = 0;
        i2 = 0;

        for (label i = 0; i < sz1; ++i)
        {
           label sl = i*sz23;

            for (label j = 0; j < sz2; ++j)
            {
                i1 = j + sl;
                i2 = i1;

                for (label k = 0; k < endIndex - filterCentre; ++k)
                {
                    tmp[i1] = 0.0;
                    label q = 0;

                    for (label p = szfilter3 - 1; p >= 0; --p, ++q)
                    {
                        tmp[i1] += tmp2[i2 + q*sz2]*filter3[p];
                    }
                    i1 += sz2;
                    i2 += sz2;
                }
                i1 += (sz2 + filterCentre);
                i2 += (sz2 + filterCentre);
            }
        }

        // Convolution summation - Along 3rd direction
        filterCentre = label(szfilter1/label(2));
        endIndex = sz1 - filterCentre;
        i1 = (szfilter2 - label(1))/label(2);
        i2 = (szfilter2 - label(1))/label(2);
        label i3 = 0;

        for (label i = 0; i < validSlice3; ++i)
        {
            for (label j = 0; j < validSlice2; ++j)
            {
                scalar sum = 0.0;
                i1 = i2 + j;

                for (label k = szfilter1 - 1; k >= 0; --k)
                {
                    sum += tmp[i1]*filter1[k];
                    i1 += sz23;
                }
                out[i3] = sum;
                ++i3;
            }
            i2 += sz2;
        }

        fld.replace(dir, out);
    }

    return fld;
}


void Foam::syntheticTurbulence::randomNumberBox::write
(
    Ostream& os
) const
{
    os.writeEntry("L", L_);
    os.writeEntry("n", n_);
    os.writeEntryIfDifferent<bool>("fixSeed", true, fixSeed_);
    os.writeEntryIfDifferent<bool>("continuous", false, continuous_);

    if (continuous_ && Pstream::master())
    {
        os.writeEntry("box", box_);
    }
}


// ************************************************************************* //
