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
    Foam::syntheticTurbulence::randomNumberBox

Description

Usage

SourceFiles
    randomNumberBox.C

\*---------------------------------------------------------------------------*/

#ifndef randomNumberBox_H
#define randomNumberBox_H

#include "List.H"
#include "dictionary.H"
#include "Random.H"
#include "fvPatch.H"
#include "primitivePatch.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace syntheticTurbulence
{

/*---------------------------------------------------------------------------*\
                      Class randomNumberBox Declaration
\*---------------------------------------------------------------------------*/

class randomNumberBox
{
    // Private Data

        //- Constant reference to patch
        const fvPatch& p_;

        //- Flag to fix the random-number generator seed to 1234 or
        //- generate a new seed based on clock-time per simulation
        const bool fixSeed_;

        //- Flag to write random-number sets at output time, and to read them
        //- on restart. Otherwise, generate new random-number sets on restart
        const bool continuous_;

        //- Integral length-scale set in local coordinates [m]
        const tensor L_;

        //- Number of faces on random-number box inlet plane (<height> <width>)
        const Pair<label> n_;

        //- Uniform cell size on random-number box inlet plane [m]
        const Pair<scalar> delta_;

        //- Random-number box ranges to be filtered by corresponding kernels
        //  e.g. for U: (Lxu, Lxv, Lxw, Lyu, Lyv, Lyw, Lzu, Lzv, Lzw)
        const List<label> ranges_;

        //- Ranges of each inlet-normal slice of the random-number box
        //  e.g. for U: (Lyu*Lzu, Lyv*Lzv, Lyw*Lzw)
        const List<label> sliceRanges_;

        //- Ranges of each physical-variable component of the random-number box
        //  e.g. for U: (Lxu*Lyu*Lzu, Lxv*Lyv*Lzv, Lxw*Lyw*Lzw)
        const List<label> componentRanges_;

        //- Indices to the first elem of last plane of the random-number box
        const List<label> lastSlice_;

        //- Random number generator
        Random rndGen_;

        //- Random-number box
        //  First List: Components of physical variable, e.g. (u,v,w)
        //  Second List: Containing all random numbers in a continuous list
        //  for each component of Cartesian coordinate system, e.g. (x,y,z)
        //  e.g. for U, box_[0] contains sets for u-component in x-y-z order
        List<List<scalar>> box_; // Vector<scalarList> box_;

        //- Vertices of virtual patch of the random-number box
        pointField patchPoints_;

        //- Faces of virtual patch of the random-number box
        faceList patchFaces_;

        //- Virtual patch of the random-number box
        primitivePatch patch_;


    // Private Member Functions

        //- Initialise the bounding box in the local coordinate system
        boundBox initialiseBoundBox() const;

        //- Initialise the object: delta
        Pair<scalar> initialiseDelta() const;

        //- Initialise the object: ranges
        List<label> initialiseRanges() const;

        //- Initialise the object: sliceRanges
        List<label> initialiseSliceRanges() const;

        //- Initialise the object: componentRanges
        List<label> initialiseComponentRanges() const;

        //- Initialise the object: lastSlice
        List<label> initialiseLastSlice() const;

        //- Initialise the object: box
        List<List<scalar>> initialiseBox();

        //- Initialise the object: patchPoints
        pointField initialisePatchPoints() const;

        //- Initialise the object: patchFaces
        faceList initialisePatchFaces() const;

        //- Initialise the object: virtualPatch
        primitivePatch initialisePatch() const;

        //- Convert integral-length scales in meters
        //- to random-number box inlet cell-size units
        tensor meterToCell(const tensor& L) const;


        //- No copy assignment
        void operator=(const randomNumberBox&) = delete;


public:

    // Constructors

        //- Default construct
        randomNumberBox() noexcept;

        //- Construct from patch and dictionary
        randomNumberBox
        (
            const fvPatch& p,
            const dictionary& dict
        );


    //- Destructor
    virtual ~randomNumberBox() = default;


    // Public Data

        //- Flag to activate debug statements
        static int debug;


    // Member Functions

        // Access

            //- Return reference to random-number box inlet patch
            const primitivePatch& patch() const noexcept
            {
                return patch_;
            }

            //- Return the integral-length scale set
            tensor L() const
            {
                return meterToCell(L_);
            }


        // Evaluation

            //- Discard current time-step random-box slice
            //- (the closest to the patch) by shifting from the back to the front
            void shift();

            //- Add a new random-box slice to the rear of the box
            void refill();

            //- Embed two-point correlations, i.e. L
            //  Apply three-dimensional "valid"-type separable
            //  convolution summation algorithm
            //  (Based on Song Ho Ahn's two-dimensional "full"-type convolution)
            Field<vector> convolve
            (
                const List<List<scalar>>& kernel
            ) const;


        // I-O

            //- Write random-number box settings
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace syntheticTurbulence
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
