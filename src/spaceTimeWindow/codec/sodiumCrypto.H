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

Class
    Foam::sodiumCrypto

Description
    Optional asymmetric encryption for boundary condition data using libsodium.

    Uses sealed boxes (crypto_box_seal) for anonymous encryption:
    - Encryption requires only the recipient's PUBLIC key
    - Decryption requires both PUBLIC and PRIVATE keys
    - Provides confidentiality and integrity

    Key format:
    - Keys are 32 bytes, represented as base64 strings for controlDict
    - Public key can be freely shared in controlDict
    - Private key is NEVER stored - entered via stdin at init time

    File format (.enc):
    - 4 bytes: magic "ENC1"
    - 4 bytes: flags (reserved)
    - 8 bytes: original data size (little-endian)
    - N bytes: sealed box ciphertext

    Workflow:
    1. Generate keypair: spaceTimeWindowKeygen
    2. Extraction: PUBLIC key in controlDict, BC data encrypted to .enc files
    3. Init: spaceTimeWindowInitCase prompts for PRIVATE key (stdin, no echo)
       and decrypts .enc files to plain format
    4. Reconstruction: runs with decrypted data

    Usage in controlDict (extraction):
    \verbatim
    spaceTimeWindowExtract1
    {
        ...
        encrypt         true;
        publicKey       "base64-encoded-32-byte-public-key";
    }
    \endverbatim

Note
    This class compiles to no-ops if FOAM_USE_SODIUM is not defined.
    Enable with: export FOAM_USE_SODIUM=1 before wmake

SourceFiles
    sodiumCrypto.C

\*---------------------------------------------------------------------------*/

#ifndef sodiumCrypto_H
#define sodiumCrypto_H

#include "fileName.H"
#include "word.H"
#include <vector>
#include <string>
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sodiumCrypto Declaration
\*---------------------------------------------------------------------------*/

class sodiumCrypto
{
public:

    // Public constants

        //- File extension for encrypted files
        static const word fileExtension;

        //- Magic number for encrypted files ("ENC1")
        static constexpr uint32_t MAGIC = 0x31434E45;

        //- Public key size in bytes
        static constexpr size_t PUBLIC_KEY_SIZE = 32;

        //- Private key size in bytes
        static constexpr size_t PRIVATE_KEY_SIZE = 32;

        //- Header size in bytes
        static constexpr size_t HEADER_SIZE = 16;


    // Static Member Functions

        //- Check if libsodium support is compiled in
        static bool available();

        //- Initialize libsodium (call once at startup)
        //  Returns true on success, false if not available or init failed
        static bool initialize();

        //- Generate a new keypair
        //  Returns true on success
        static bool generateKeypair
        (
            std::vector<uint8_t>& publicKey,
            std::vector<uint8_t>& privateKey
        );

        //- Encode binary data to base64 string
        static std::string toBase64(const std::vector<uint8_t>& data);

        //- Decode base64 string to binary data
        //  Returns empty vector on error
        static std::vector<uint8_t> fromBase64(const std::string& b64);

        //- Encrypt data using public key and write to file
        static void encryptToFile
        (
            const fileName& path,
            const std::vector<uint8_t>& plaintext,
            const std::vector<uint8_t>& publicKey
        );

        //- Encrypt data using public key (base64) and write to file
        static void encryptToFile
        (
            const fileName& path,
            const std::vector<uint8_t>& plaintext,
            const std::string& publicKeyBase64
        );

        //- Decrypt file using keypair, return plaintext
        static std::vector<uint8_t> decryptFromFile
        (
            const fileName& path,
            const std::vector<uint8_t>& publicKey,
            const std::vector<uint8_t>& privateKey
        );

        //- Decrypt file using keypair (base64), return plaintext
        static std::vector<uint8_t> decryptFromFile
        (
            const fileName& path,
            const std::string& publicKeyBase64,
            const std::string& privateKeyBase64
        );

        //- Check if file is encrypted (by extension or magic)
        static bool isEncryptedFile(const fileName& path);

        //- Read private key from stdin without echo
        //  Prompts user and reads key in base64 format
        static std::string readPrivateKeyFromStdin(const std::string& prompt);

        //- Derive public key from private key
        //  For X25519/sealed boxes, public key = crypto_scalarmult_base(private key)
        static std::vector<uint8_t> derivePublicKey
        (
            const std::vector<uint8_t>& privateKey
        );

        //- Derive public key from private key (base64 versions)
        static std::string derivePublicKeyBase64(const std::string& privateKeyBase64);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
