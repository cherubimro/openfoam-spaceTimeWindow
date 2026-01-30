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

\*---------------------------------------------------------------------------*/

#include "sodiumCrypto.H"
#include "error.H"
#include <fstream>
#include <cstring>

#ifdef FOAM_USE_SODIUM
#include <sodium.h>
#include <termios.h>
#include <unistd.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::sodiumCrypto::fileExtension("enc");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sodiumCrypto::available()
{
#ifdef FOAM_USE_SODIUM
    return true;
#else
    return false;
#endif
}


bool Foam::sodiumCrypto::initialize()
{
#ifdef FOAM_USE_SODIUM
    static bool initialized = false;
    static bool initResult = false;

    if (!initialized)
    {
        initResult = (sodium_init() >= 0);
        initialized = true;
    }

    return initResult;
#else
    return false;
#endif
}


bool Foam::sodiumCrypto::generateKeypair
(
    std::vector<uint8_t>& publicKey,
    std::vector<uint8_t>& privateKey
)
{
#ifdef FOAM_USE_SODIUM
    if (!initialize())
    {
        return false;
    }

    publicKey.resize(crypto_box_PUBLICKEYBYTES);
    privateKey.resize(crypto_box_SECRETKEYBYTES);

    crypto_box_keypair(publicKey.data(), privateKey.data());

    return true;
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return false;
#endif
}


std::string Foam::sodiumCrypto::toBase64(const std::vector<uint8_t>& data)
{
#ifdef FOAM_USE_SODIUM
    if (data.empty())
    {
        return "";
    }

    // Calculate required buffer size
    size_t b64MaxLen = sodium_base64_ENCODED_LEN(
        data.size(),
        sodium_base64_VARIANT_ORIGINAL
    );

    std::vector<char> b64(b64MaxLen);

    sodium_bin2base64(
        b64.data(),
        b64MaxLen,
        data.data(),
        data.size(),
        sodium_base64_VARIANT_ORIGINAL
    );

    return std::string(b64.data());
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return "";
#endif
}


std::vector<uint8_t> Foam::sodiumCrypto::fromBase64(const std::string& b64)
{
#ifdef FOAM_USE_SODIUM
    if (b64.empty())
    {
        return {};
    }

    // Output buffer (will be smaller than input)
    std::vector<uint8_t> bin(b64.size());
    size_t binLen = 0;

    int result = sodium_base642bin(
        bin.data(),
        bin.size(),
        b64.c_str(),
        b64.size(),
        nullptr,  // ignore characters
        &binLen,
        nullptr,  // end pointer
        sodium_base64_VARIANT_ORIGINAL
    );

    if (result != 0)
    {
        return {};
    }

    bin.resize(binLen);
    return bin;
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return {};
#endif
}


void Foam::sodiumCrypto::encryptToFile
(
    const fileName& path,
    const std::vector<uint8_t>& plaintext,
    const std::vector<uint8_t>& publicKey
)
{
#ifdef FOAM_USE_SODIUM
    if (!initialize())
    {
        FatalErrorInFunction
            << "Failed to initialize libsodium"
            << abort(FatalError);
    }

    if (publicKey.size() != crypto_box_PUBLICKEYBYTES)
    {
        FatalErrorInFunction
            << "Invalid public key size: " << publicKey.size()
            << " (expected " << crypto_box_PUBLICKEYBYTES << ")"
            << abort(FatalError);
    }

    // Sealed box ciphertext size
    size_t ciphertextLen = plaintext.size() + crypto_box_SEALBYTES;
    std::vector<uint8_t> ciphertext(ciphertextLen);

    // Encrypt
    int result = crypto_box_seal(
        ciphertext.data(),
        plaintext.data(),
        plaintext.size(),
        publicKey.data()
    );

    if (result != 0)
    {
        FatalErrorInFunction
            << "Encryption failed for file: " << path
            << abort(FatalError);
    }

    // Write to file with header
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    // Header: magic (4) + flags (4) + original size (8) = 16 bytes
    uint32_t magic = MAGIC;
    uint32_t flags = 0;
    uint64_t origSize = plaintext.size();

    ofs.write(reinterpret_cast<const char*>(&magic), 4);
    ofs.write(reinterpret_cast<const char*>(&flags), 4);
    ofs.write(reinterpret_cast<const char*>(&origSize), 8);

    // Ciphertext
    ofs.write(reinterpret_cast<const char*>(ciphertext.data()), ciphertext.size());

#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
#endif
}


void Foam::sodiumCrypto::encryptToFile
(
    const fileName& path,
    const std::vector<uint8_t>& plaintext,
    const std::string& publicKeyBase64
)
{
    std::vector<uint8_t> publicKey = fromBase64(publicKeyBase64);

    if (publicKey.empty())
    {
        FatalErrorInFunction
            << "Invalid base64 public key"
            << abort(FatalError);
    }

    encryptToFile(path, plaintext, publicKey);
}


std::vector<uint8_t> Foam::sodiumCrypto::decryptFromFile
(
    const fileName& path,
    const std::vector<uint8_t>& publicKey,
    const std::vector<uint8_t>& privateKey
)
{
#ifdef FOAM_USE_SODIUM
    if (!initialize())
    {
        FatalErrorInFunction
            << "Failed to initialize libsodium"
            << abort(FatalError);
    }

    if (publicKey.size() != crypto_box_PUBLICKEYBYTES)
    {
        FatalErrorInFunction
            << "Invalid public key size: " << publicKey.size()
            << abort(FatalError);
    }

    if (privateKey.size() != crypto_box_SECRETKEYBYTES)
    {
        FatalErrorInFunction
            << "Invalid private key size: " << privateKey.size()
            << abort(FatalError);
    }

    // Read file
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if (!ifs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for reading: " << path
            << abort(FatalError);
    }

    size_t fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    if (fileSize < HEADER_SIZE + crypto_box_SEALBYTES)
    {
        FatalErrorInFunction
            << "File too small to be encrypted: " << path
            << abort(FatalError);
    }

    // Read header
    uint32_t magic;
    uint32_t flags;
    uint64_t origSize;

    ifs.read(reinterpret_cast<char*>(&magic), 4);
    ifs.read(reinterpret_cast<char*>(&flags), 4);
    ifs.read(reinterpret_cast<char*>(&origSize), 8);

    if (magic != MAGIC)
    {
        FatalErrorInFunction
            << "Invalid encrypted file (bad magic): " << path
            << abort(FatalError);
    }

    // Read ciphertext
    size_t ciphertextLen = fileSize - HEADER_SIZE;
    std::vector<uint8_t> ciphertext(ciphertextLen);
    ifs.read(reinterpret_cast<char*>(ciphertext.data()), ciphertextLen);

    // Decrypt
    std::vector<uint8_t> plaintext(origSize);

    int result = crypto_box_seal_open(
        plaintext.data(),
        ciphertext.data(),
        ciphertext.size(),
        publicKey.data(),
        privateKey.data()
    );

    if (result != 0)
    {
        FatalErrorInFunction
            << "Decryption failed for file: " << path << nl
            << "    Possible causes:" << nl
            << "    - Wrong private key" << nl
            << "    - Corrupted file" << nl
            << "    - File encrypted with different public key"
            << abort(FatalError);
    }

    return plaintext;
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return {};
#endif
}


std::vector<uint8_t> Foam::sodiumCrypto::decryptFromFile
(
    const fileName& path,
    const std::string& publicKeyBase64,
    const std::string& privateKeyBase64
)
{
    std::vector<uint8_t> publicKey = fromBase64(publicKeyBase64);
    std::vector<uint8_t> privateKey = fromBase64(privateKeyBase64);

    if (publicKey.empty())
    {
        FatalErrorInFunction
            << "Invalid base64 public key"
            << abort(FatalError);
    }

    if (privateKey.empty())
    {
        FatalErrorInFunction
            << "Invalid base64 private key"
            << abort(FatalError);
    }

    return decryptFromFile(path, publicKey, privateKey);
}


bool Foam::sodiumCrypto::isEncryptedFile(const fileName& path)
{
    // Check by extension first
    if (path.ext() == fileExtension)
    {
        return true;
    }

    // Check by magic number
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.good())
    {
        return false;
    }

    uint32_t magic;
    ifs.read(reinterpret_cast<char*>(&magic), 4);

    if (!ifs.good() || ifs.gcount() != 4)
    {
        return false;
    }

    return magic == MAGIC;
}


std::string Foam::sodiumCrypto::readPrivateKeyFromStdin(const std::string& prompt)
{
#ifdef FOAM_USE_SODIUM
    // Print prompt to stderr (not stdout, to avoid capture)
    std::cerr << prompt << std::flush;

    // Disable echo
    struct termios oldTerm, newTerm;
    tcgetattr(STDIN_FILENO, &oldTerm);
    newTerm = oldTerm;
    newTerm.c_lflag &= ~(ECHO | ECHOE | ECHOK | ECHONL);
    tcsetattr(STDIN_FILENO, TCSANOW, &newTerm);

    // Read line
    std::string key;
    std::getline(std::cin, key);

    // Restore echo
    tcsetattr(STDIN_FILENO, TCSANOW, &oldTerm);

    // Print newline since user's Enter wasn't echoed
    std::cerr << std::endl;

    // Trim whitespace
    size_t start = key.find_first_not_of(" \t\r\n");
    size_t end = key.find_last_not_of(" \t\r\n");

    if (start == std::string::npos)
    {
        return "";
    }

    return key.substr(start, end - start + 1);
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return "";
#endif
}


std::vector<uint8_t> Foam::sodiumCrypto::derivePublicKey
(
    const std::vector<uint8_t>& privateKey
)
{
#ifdef FOAM_USE_SODIUM
    if (!initialize())
    {
        FatalErrorInFunction
            << "Failed to initialize libsodium"
            << abort(FatalError);
    }

    if (privateKey.size() != crypto_box_SECRETKEYBYTES)
    {
        FatalErrorInFunction
            << "Invalid private key size: " << privateKey.size()
            << " (expected " << crypto_box_SECRETKEYBYTES << ")"
            << abort(FatalError);
    }

    std::vector<uint8_t> publicKey(crypto_box_PUBLICKEYBYTES);

    // Derive public key from private key using scalar multiplication
    // For X25519: public = private * basepoint
    crypto_scalarmult_base(publicKey.data(), privateKey.data());

    return publicKey;
#else
    FatalErrorInFunction
        << "Encryption not available. Rebuild with FOAM_USE_SODIUM=1"
        << abort(FatalError);
    return {};
#endif
}


std::string Foam::sodiumCrypto::derivePublicKeyBase64
(
    const std::string& privateKeyBase64
)
{
    std::vector<uint8_t> privateKey = fromBase64(privateKeyBase64);

    if (privateKey.empty())
    {
        FatalErrorInFunction
            << "Invalid base64 private key"
            << abort(FatalError);
    }

    std::vector<uint8_t> publicKey = derivePublicKey(privateKey);
    return toBase64(publicKey);
}


// ************************************************************************* //
