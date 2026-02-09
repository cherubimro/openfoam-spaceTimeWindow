// Explicit ghost cell overwrite from extracted global simulation data.
// Called inside PIMPLE loop when useGhostCells is true.
// Reads cell-centre values from constant/boundaryData/ghostCells/<time>/
// and directly overwrites the ghost cell values in U and p.
{
    const labelList& ghostCells = mesh.cellZones()["ghostCells"];

    if (ghostCells.size() > 0)
    {
        fileName ghostDir =
            mesh.time().path()
          / "constant" / "boundaryData" / "ghostCells"
          / runTime.timeName();

        // Read and apply U
        {
            IFstream is(ghostDir / "U");
            if (is.good())
            {
                // Skip FoamFile header if present
                token firstToken(is);
                if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
                {
                    dictionary headerDict(is);
                }
                else
                {
                    is.putBack(firstToken);
                }

                vectorField ghostU(is);
                forAll(ghostCells, i)
                {
                    U[ghostCells[i]] = ghostU[i];
                }
            }
            else
            {
                WarningInFunction
                    << "Ghost cell U data not found at " << ghostDir / "U"
                    << endl;
            }
        }

        // Read and apply p
        {
            IFstream is(ghostDir / "p");
            if (is.good())
            {
                // Skip FoamFile header if present
                token firstToken(is);
                if (firstToken.isWord() && firstToken.wordToken() == "FoamFile")
                {
                    dictionary headerDict(is);
                }
                else
                {
                    is.putBack(firstToken);
                }

                scalarField ghostP(is);
                forAll(ghostCells, i)
                {
                    p[ghostCells[i]] = ghostP[i];
                }
            }
            else
            {
                WarningInFunction
                    << "Ghost cell p data not found at " << ghostDir / "p"
                    << endl;
            }
        }

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();
    }
}
