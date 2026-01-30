/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Application
    spaceTimeWindowKeygen

Description
    Generate a keypair for spaceTimeWindow encryption.

    Creates a public/private keypair using libsodium's sealed box scheme.
    The public key is used during extraction (in controlDict), and the
    private key is used during case initialization to decrypt boundary data.

    SECURITY: Store the private key securely and NEVER include it in
    any case files or version control. The private key should only be
    entered interactively when running spaceTimeWindowInitCase.

Usage
    spaceTimeWindowKeygen [OPTIONS]

    Options:
        -format <text|json>   Output format (default: text)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "sodiumCrypto.H"
#include <iostream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a keypair for spaceTimeWindow encryption"
    );

    argList::addOption
    (
        "format",
        "text|json",
        "Output format (default: text)"
    );

    argList::noParallel();
    argList::noFunctionObjects();
    argList::noBanner();

    #include "setRootCase.H"

#ifdef FOAM_USE_SODIUM
    if (!sodiumCrypto::available())
    {
        FatalErrorInFunction
            << "libsodium failed to initialize" << nl
            << exit(FatalError);
    }

    // Generate keypair
    std::vector<uint8_t> publicKey;
    std::vector<uint8_t> privateKey;

    if (!sodiumCrypto::generateKeypair(publicKey, privateKey))
    {
        FatalErrorInFunction
            << "Failed to generate keypair" << nl
            << exit(FatalError);
    }

    // Convert to base64
    std::string publicKeyB64 = sodiumCrypto::toBase64(publicKey);
    std::string privateKeyB64 = sodiumCrypto::toBase64(privateKey);

    // Output format
    word format = args.getOrDefault<word>("format", "text");

    if (format == "json")
    {
        std::cout << "{" << std::endl;
        std::cout << "  \"publicKey\": \"" << publicKeyB64 << "\"," << std::endl;
        std::cout << "  \"privateKey\": \"" << privateKeyB64 << "\"" << std::endl;
        std::cout << "}" << std::endl;
    }
    else
    {
        std::cout << "========================================" << std::endl;
        std::cout << "spaceTimeWindow Keypair Generator" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << std::endl;
        std::cout << "PUBLIC KEY (add to controlDict):" << std::endl;
        std::cout << publicKeyB64 << std::endl;
        std::cout << std::endl;
        std::cout << "PRIVATE KEY (keep SECRET - enter at init time):" << std::endl;
        std::cout << privateKeyB64 << std::endl;
        std::cout << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "SECURITY NOTES:" << std::endl;
        std::cout << "- Store the private key in a secure location" << std::endl;
        std::cout << "- NEVER include the private key in case files" << std::endl;
        std::cout << "- NEVER commit the private key to version control" << std::endl;
        std::cout << "- The private key is required for spaceTimeWindowInitCase" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << std::endl;
        std::cout << "Usage in controlDict:" << std::endl;
        std::cout << "    spaceTimeWindowExtract1" << std::endl;
        std::cout << "    {" << std::endl;
        std::cout << "        type            spaceTimeWindowExtract;" << std::endl;
        std::cout << "        ..." << std::endl;
        std::cout << "        encrypt         true;" << std::endl;
        std::cout << "        publicKey       \"" << publicKeyB64 << "\";" << std::endl;
        std::cout << "    }" << std::endl;
    }

    return 0;
#else
    FatalErrorInFunction
        << "spaceTimeWindowKeygen requires libsodium support" << nl
        << "Rebuild with FOAM_USE_SODIUM=1 to enable encryption" << nl
        << exit(FatalError);

    return 1;
#endif
}


// ************************************************************************* //
