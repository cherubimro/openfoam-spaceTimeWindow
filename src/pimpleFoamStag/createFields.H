#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Persistent HbyA field for function object access (e.g. spaceTimeWindowExtract).
// The intermediate velocity u* is updated each PIMPLE iteration in pEqn.H.
volVectorField HbyA_stored
(
    IOobject
    (
        "HbyA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::REGISTER
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"


// --- Staggered-specific fields ---

// Primary face-normal velocity (the staggered kinematic variable)
Info<< "Creating face-normal velocity field Un\n" << endl;
surfaceScalarField Un
(
    IOobject
    (
        "Un",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (fvc::interpolate(U) & (mesh.Sf()/mesh.magSf()))
);

// Unit face-normal vector (cached)
const surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());

// Reconstruction tensor: inv(surfaceSum(SfHat * Sf))
// Maps face-normal data back to cell-centered vectors via least-squares.
// On 2D meshes the empty-direction component is singular, so we add a
// small diagonal perturbation (same fix as fvc::reconstruct uses
// implicitly via its extrapolatedCalculated BC).
volTensorField reconTensor
(
    IOobject
    (
        "reconTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inv
    (
        fvc::surfaceSum(SfHat * mesh.Sf())
      + dimensionedTensor
        (
            "small",
            dimArea,
            symmTensor(SMALL, 0, 0, SMALL, 0, SMALL)
        )
    )
);
