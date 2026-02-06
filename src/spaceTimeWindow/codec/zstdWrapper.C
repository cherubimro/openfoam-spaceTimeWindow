/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024-2025 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "zstdWrapper.H"
#include "error.H"
#include <fstream>
#include <cstring>

#ifdef FOAM_USE_ZSTD
#include <zstd.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::zstdWrapper::fileExtension("zstd");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::zstdWrapper::available()
{
#ifdef FOAM_USE_ZSTD
    return true;
#else
    return false;
#endif
}


std::vector<uint8_t> Foam::zstdWrapper::compress
(
    const std::vector<uint8_t>& data,
    int level
)
{
#ifdef FOAM_USE_ZSTD
    size_t const maxDstSize = ZSTD_compressBound(data.size());
    std::vector<uint8_t> compressed(maxDstSize);

    size_t const compressedSize = ZSTD_compress
    (
        compressed.data(),
        maxDstSize,
        data.data(),
        data.size(),
        level
    );

    if (ZSTD_isError(compressedSize))
    {
        FatalErrorInFunction
            << "Zstd compression failed: " << ZSTD_getErrorName(compressedSize)
            << abort(FatalError);
    }

    compressed.resize(compressedSize);
    return compressed;
#else
    FatalErrorInFunction
        << "Zstd compression not available. Rebuild with FOAM_USE_ZSTD=1"
        << abort(FatalError);
    return {};
#endif
}


std::vector<uint8_t> Foam::zstdWrapper::decompress
(
    const std::vector<uint8_t>& data
)
{
#ifdef FOAM_USE_ZSTD
    unsigned long long const decompressedSize =
        ZSTD_getFrameContentSize(data.data(), data.size());

    if (decompressedSize == ZSTD_CONTENTSIZE_ERROR)
    {
        FatalErrorInFunction
            << "Data is not valid zstd compressed"
            << abort(FatalError);
    }

    if (decompressedSize == ZSTD_CONTENTSIZE_UNKNOWN)
    {
        FatalErrorInFunction
            << "Zstd frame does not contain original size"
            << abort(FatalError);
    }

    std::vector<uint8_t> decompressed(decompressedSize);

    size_t const actualSize = ZSTD_decompress
    (
        decompressed.data(),
        decompressedSize,
        data.data(),
        data.size()
    );

    if (ZSTD_isError(actualSize))
    {
        FatalErrorInFunction
            << "Zstd decompression failed: " << ZSTD_getErrorName(actualSize)
            << abort(FatalError);
    }

    decompressed.resize(actualSize);
    return decompressed;
#else
    FatalErrorInFunction
        << "Zstd compression not available. Rebuild with FOAM_USE_ZSTD=1"
        << abort(FatalError);
    return {};
#endif
}


void Foam::zstdWrapper::compressToFile
(
    const fileName& path,
    const std::vector<uint8_t>& data,
    int level
)
{
    std::vector<uint8_t> compressed = compress(data, level);

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    ofs.write
    (
        reinterpret_cast<const char*>(compressed.data()),
        compressed.size()
    );
}


std::vector<uint8_t> Foam::zstdWrapper::decompressFromFile
(
    const fileName& path
)
{
    std::ifstream ifs(path, std::ios::binary | std::ios::ate);
    if (!ifs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for reading: " << path
            << abort(FatalError);
    }

    size_t fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);

    std::vector<uint8_t> compressed(fileSize);
    ifs.read(reinterpret_cast<char*>(compressed.data()), fileSize);

    return decompress(compressed);
}


bool Foam::zstdWrapper::isZstdFile(const fileName& path)
{
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


// ************************************************************************* //
