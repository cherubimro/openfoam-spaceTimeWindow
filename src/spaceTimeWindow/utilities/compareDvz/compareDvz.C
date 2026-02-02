/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Description
    Utility to compare two DVZ files and report max difference.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "deltaVarintCodec.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addArgument("file1", "First DVZ file");
    argList::addArgument("file2", "Second DVZ file");
    argList::noParallel();
    argList::noFunctionObjects();

    #include "setRootCase.H"

    fileName file1 = args.get<fileName>(1);
    fileName file2 = args.get<fileName>(2);

    Info<< "Comparing:" << nl
        << "  " << file1 << nl
        << "  " << file2 << nl << endl;

    // Read both files as vectors
    vectorField field1 = deltaVarintCodec::readVector(file1);
    vectorField field2 = deltaVarintCodec::readVector(file2);

    if (field1.size() != field2.size())
    {
        FatalErrorInFunction
            << "Field sizes differ: " << field1.size() << " vs " << field2.size()
            << exit(FatalError);
    }

    Info<< "Field size: " << field1.size() << nl << endl;

    // Compute statistics
    scalar maxDiff = 0.0;
    scalar sumDiff = 0.0;
    label maxDiffIdx = 0;

    for (label i = 0; i < field1.size(); ++i)
    {
        scalar diff = mag(field1[i] - field2[i]);
        sumDiff += diff;

        if (diff > maxDiff)
        {
            maxDiff = diff;
            maxDiffIdx = i;
        }
    }

    scalar avgDiff = sumDiff / field1.size();

    Info<< "Results:" << nl
        << "  Max difference:  " << maxDiff << " at index " << maxDiffIdx << nl
        << "  Avg difference:  " << avgDiff << nl;

    if (maxDiff > 0)
    {
        Info<< nl << "Value at max diff:" << nl
            << "  File1: " << field1[maxDiffIdx] << nl
            << "  File2: " << field2[maxDiffIdx] << nl;
    }

    if (maxDiff < 1e-5)
    {
        Info<< nl << "Files are effectively identical (max diff < 1e-5)" << endl;
    }
    else
    {
        Info<< nl << "WARNING: Significant differences found!" << endl;
    }

    return 0;
}


// ************************************************************************* //
