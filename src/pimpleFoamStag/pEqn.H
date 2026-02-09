// Staggered pressure equation: pressure correction acts directly on
// face-normal velocity Un.  NO Rhie-Chow interpolation.

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

// Update persistent HbyA field for function object access.
HbyA_stored = HbyA;

// Project HbyA onto face normals to get predicted face-normal velocity
surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
surfaceScalarField UnStar
(
    "UnStar",
    (fvc::interpolate(HbyA) & mesh.Sf()) / mesh.magSf()
);

// Predicted face flux
surfaceScalarField phiStar("phiStar", UnStar * mesh.magSf());

// DDT correction for temporal consistency
if (pimple.ddtCorr())
{
    phiStar += MRF.zeroFilter(rAUf*fvc::ddtCorr(U, phi, Uf));
}

MRF.makeRelative(phiStar);

if (p.needReference())
{
    fvc::makeRelative(phiStar, U);
    adjustPhi(phiStar, U, p);
    fvc::makeAbsolute(phiStar, U);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update pressure BCs for flux consistency
constrainPressure(p, U, phiStar, rAUf, MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiStar)
    );

    pEqn.setReference(pRefCell, pRefValue);
    fvOptions.constrain(pEqn);

    pEqn.solve(p.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        // EXACT face-normal pressure correction (no Rhie-Chow)
        Un = UnStar - rAUf * fvc::snGrad(p);
        phi = Un * mesh.magSf();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Reconstruct cell-centered velocity from corrected pressure
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make fluxes relative to mesh motion
fvc::makeRelative(phi, U);
