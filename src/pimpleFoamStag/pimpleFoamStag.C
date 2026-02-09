/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    pimpleFoamStag

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    using a Mahesh-style hybrid staggered mesh approach.

    Face-normal velocity Un is the primary kinematic variable. The pressure
    correction acts directly on Un, eliminating Rhie-Chow interpolation
    and providing exact discrete mass conservation.

    Momentum is assembled at cell centers (reusing all standard fvm/fvc
    operators, turbulence models, and boundary conditions). Cell-centered
    velocity U is reconstructed from Un via least-squares each iteration.

    \heading Algorithm per time step, per PIMPLE iteration:
    1. RECONSTRUCT: U_cell from Un (faces) via least-squares
    2. MOMENTUM:    Cell-centered UEqn (standard fvm operators)
    3. PROJECT:     HbyA -> face-normal UnStar
    4. PRESSURE:    fvm::laplacian(rAUf, p) == fvc::div(phiStar)
    5. CORRECT:     Un = UnStar - rAUf * snGrad(p)  [EXACT, no Rhie-Chow]
    6. TURBULENCE:  turbulence->correct() on final iteration

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " using hybrid staggered mesh (Mahesh-style)."
        " Face-normal velocity is the primary kinematic variable."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // Ghost cell overwrite: check if ghostCells zone exists
    const bool useGhostCells =
        mesh.cellZones().findZoneID("ghostCells") >= 0;

    if (useGhostCells)
    {
        Info<< "Ghost cell mode enabled: "
            << mesh.cellZones()["ghostCells"].size()
            << " ghost cells will be overwritten each timestep" << endl;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    // Recompute reconTensor after mesh motion
                    reconTensor = inv
                    (
                        fvc::surfaceSum(SfHat * mesh.Sf())
                      + dimensionedTensor
                        (
                            "small",
                            dimArea,
                            symmTensor(SMALL, 0, 0, SMALL, 0, SMALL)
                        )
                    );
                }
            }

            // Step 0: Overwrite ghost cells from extracted data
            if (useGhostCells)
            {
                #include "overwriteGhostCells.H"
            }

            // Step 1: Reconstruct cell-centered U from face-normal Un
            #include "reconstructU.H"

            // Step 2: Assemble and solve cell-centered momentum equation
            #include "UEqn.H"

            // Step 3-5: Pressure correction loop (acts on Un directly)
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            // Re-overwrite ghost cells after pressure correction
            if (useGhostCells)
            {
                #include "overwriteGhostCells.H"
            }

            // Step 6: Turbulence correction
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
