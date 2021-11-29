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

#include "softWall.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(softWall, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        softWall,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::softWall::softWall
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    anchor_(),
    refAttachmentPt_(),
    wallNormal_(),
    psi_(),
    C_()
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::softWall::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);
    restraintForce = Zero;
    restraintMoment = Zero;

    const vector r(restraintPosition - anchor_);
    // scalar magR = mag(r);
    // r /= (magR + VSMALL);

    const vector v(motion.velocity(restraintPosition));

    const scalar d = (wallNormal_/mag(wallNormal_)) & r;

    const vector rDir(r/(mag(r) + VSMALL));

    const scalar wn = 3.14/C_;
    const scalar damping = psi_*2*wn;
    const scalar stiffness = sqr(wn);

    if (d < 0)
    {
        restraintForce = (-damping*(rDir & v) + stiffness*d)*rDir;

        restraintMoment = (p ^ force);
    }

    if (motion.report())
    {
        /* Info<< " attachmentPt - anchor " << r*magR
            << " spring length " << magR
            << " force " << restraintForce
            << endl;*/
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::softWall::read
(
    const dictionary& sDoFRBMRDict
)
{
    if (sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict))
    {
        sDoFRBMRCoeffs_.readEntry("anchor", anchor_);
        sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
        sDoFRBMRCoeffs_.readEntry("wallNormal", wallNormal_);
        sDoFRBMRCoeffs_.readEntry("psi", psi_);
        sDoFRBMRCoeffs_.readEntry("C", C_);

        return true;
    }

    return false;
}


void Foam::sixDoFRigidBodyMotionRestraints::softWall::write
(
    Ostream& os
) const
{
    os.writeEntry("anchor", anchor_);
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("wallNormal", wallNormal_);
    os.writeEntry("psi", psi_);
    os.writeEntry("C", C_);
}

// ************************************************************************* //
