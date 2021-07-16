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
    Foam::solidBodyFvGeometryScheme

Description

SourceFiles
    solidBodyFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef solidBodyFvGeometryScheme_H
#define solidBodyFvGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class solidBodyFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class solidBodyFvGeometryScheme
:
    public basicFvGeometryScheme
{
    // Private Data

        //- Partial update flag
        bool partialUpdate_;

        //- Cache mesh motion flag
        bool cacheMotion_;

        //- Flag to indicate that the cache has been initialised
        bool cacheInitialised_;

        //- Changed face IDs
        labelList changedFaces_;

        //- Changed patch IDs
        labelList changedPatches_;


    // Private Member Functions

        //- Set the mesh motion data (point, face IDs)
        void setMeshMotionData();

        //- No copy construct
        solidBodyFvGeometryScheme(const solidBodyFvGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const solidBodyFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("solidBody");


    // Constructors

        //- Construct from mesh
        solidBodyFvGeometryScheme(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~solidBodyFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
