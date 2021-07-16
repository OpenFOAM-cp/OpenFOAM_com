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

#include "solidBodyFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyFvGeometryScheme, 0);
    addToRunTimeSelectionTable
    (
        fvGeometryScheme,
        solidBodyFvGeometryScheme,
        dict
    );
}


void Foam::solidBodyFvGeometryScheme::setMeshMotionData()
{
    if (!cacheInitialised_ || !cacheMotion_)
    {
        DebugInFunction << "Creating cache" << endl;

        const pointField& oldPoints = mesh_.oldPoints();
        const pointField& currPoints = mesh_.points();

        if (oldPoints.size() != currPoints.size())
        {
            FatalErrorInFunction
                << "Old and current points sizes must be the same. "
                << "Old points:" << oldPoints.size()
                << " Current points:" << currPoints.size()
                << abort(FatalError);
        }


        bitSet changedPoints(oldPoints.size());

        // Find the changed points
        forAll(changedPoints, pointi)
        {
            if (oldPoints[pointi] != currPoints[pointi])
            {
                changedPoints.set(pointi);
            }
        }

        // Quick return if no points have moved
        if (returnReduce(changedPoints.count(), sumOp<label>()) == 0)
        {
            return;
        }

        // Set changed face and patch IDs
        DynamicList<label> changedFaces(mesh_.faces().size());
        DynamicList<label> changedPatches(mesh_.faces().size());

        const faceList& faces = mesh_.faces();

        for (label facei = 0; facei < mesh_.nInternalFaces(); ++facei)
        {
            const face& f = faces[facei];

            for (label fpi : f)
            {
                if (changedPoints[fpi])
                {
                    changedFaces.append(facei);
                    changedPatches.append(-1);
                    break;
                }
            }
        }

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            forAll(pp, patchFacei)
            {
                const label facei = pp.start() + patchFacei;
                const face& f = faces[facei];

                for (label fpi : f)
                {
                    if (changedPoints[fpi])
                    {
                        changedFaces.append(facei);
                        changedPatches.append(patchi);
                        break;
                    }
                }
            }
        }

        changedFaces_.transfer(changedFaces);
        changedPatches_.transfer(changedPatches);
    }

    cacheInitialised_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyFvGeometryScheme::solidBodyFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicFvGeometryScheme(mesh, dict),
    partialUpdate_(dict.get<bool>("partialUpdate")),
    cacheMotion_(dict.get<bool>("cacheMotion")),
    cacheInitialised_(false),
    changedFaces_(),
    changedPatches_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidBodyFvGeometryScheme::movePoints()
{
    // Note: not calling fvGeometryScheme::movePoints since we want to perform
    // our own geometry manipulations

    bool haveGeometry =
        mesh_.hasCellCentres()
     && mesh_.hasFaceCentres()
     && mesh_.hasCellVolumes()
     && mesh_.hasFaceAreas();

    if (!haveGeometry)
    {
        DebugInFunction
            << "Creating initial geometry using primitiveMesh::updateGeom"
            << endl;

        const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
        return;
    }

    if (mesh_.moving())
    {
        setMeshMotionData();

        DebugInFunction << "Performing partial meshPhi construction" << endl;

        // Set the mesh flux
        auto& meshPhi = const_cast<fvMesh&>(mesh_).setPhi();
        auto& meshPhii = meshPhi.primitiveFieldRef();
        auto& meshPhiBf = meshPhi.boundaryFieldRef();

        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const faceList& faces = mesh_.faces();
        const pointField& oldPoints = mesh_.oldPoints();
        const pointField& currPoints = mesh_.points();

        const scalar rdt = 1.0/mesh_.time().deltaTValue();

        //meshPhi == dimensionedScalar(dimVolume/dimTime, Zero);
        meshPhii = Zero;
        meshPhiBf == Zero;

        forAll(changedFaces_, i)
        {
            const face& f = faces[changedFaces_[i]];

            if (changedPatches_[i] == -1)
            {
                const label facei = changedFaces_[i];
                meshPhii[facei] = f.sweptVol(oldPoints, currPoints)*rdt;
            }
            else
            {
                const label patchi = changedPatches_[i];
                const label patchFacei = changedFaces_[i] - pbm[patchi].start();
                meshPhiBf[patchi][patchFacei] =
                    f.sweptVol(oldPoints, currPoints)*rdt;
            }
        }

        // Note: ACMI can delete the cell centres and volumes independently
        // from the face centres and areas. This will lead to an infinite loop
        // when performing partial updates.
        if (partialUpdate_ && haveGeometry)
        {
            // Keep base geometry and update as needed
            DebugInFunction << "Performing partial geometry update" << endl;

            // Initialise geometry using the old/existing values
            vectorField faceCentres(mesh_.faceCentres());
            vectorField faceAreas(mesh_.faceAreas());

            primitiveMeshTools::updateFaceCentresAndAreas
            (
                mesh_,
                changedFaces_,
                mesh_.points(),
                faceCentres,
                faceAreas
            );

            vectorField cellCentres(mesh_.cellCentres());
            scalarField cellVolumes(mesh_.cellVolumes());

            primitiveMeshTools::updateCellCentresAndVols
            (
                mesh_,
                changedFaces_,
                faceCentres,
                faceAreas,
                cellCentres,
                cellVolumes
            );

            const_cast<fvMesh&>(mesh_).primitiveMesh::resetGeometry
            (
                std::move(faceCentres),
                std::move(faceAreas),
                std::move(cellCentres),
                std::move(cellVolumes)
            );
        }
        else
        {
            DebugInFunction
                << "Performing complete geometry clear and update" << endl;

            // Clear out old geometry
            // Note: this recreates the old primitiveMesh::movePoints behaviour
            const_cast<fvMesh&>(mesh_).primitiveMesh::clearGeom();

            // Use lower level to calculate the geometry
            const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
        }
    }
    else
    {
        DebugInFunction << "Performing complete geometry update" << endl;

        // Use lower level to calculate the geometry
        const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
    }
}


// ************************************************************************* //
