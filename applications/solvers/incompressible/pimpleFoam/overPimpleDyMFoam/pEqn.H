volScalarField rAU("rAU", 1.0/UEqn.A());

surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
volVectorField H("H", UEqn.H());

volVectorField HbyA("HbyA", U);
HbyA = constrainHbyA(rAU*H, U, p);

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

phiHbyA = fvc::flux(HbyA);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

// WIP: To adjust phi on fringe faces to help mass
// conservation
if (adjustFringe)
{
    fvc::makeRelative(phiHbyA, U);
    oversetAdjustPhi(phiHbyA, U);
    fvc::makeAbsolute(phiHbyA, U);
}

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
        pEqn.relax();
        U =
            cellMask*
            (
                HbyA -
                rAU*fvc::reconstruct((pEqn.flux())/rAUf)
            );
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

// Excludes error in interpolated/hole cells
#include "continuityErrs.H"

{
    Uf = fvc::interpolate(U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    Uf += n*(phi/mesh.magSf() - (n & Uf));
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
phi *= faceMask;
