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

Class
    Foam::zstdWrapper

Description
    Optional zstd compression layer for boundary condition data.

    Wraps DVZ/DVZT encoded data with zstd compression at level 3 (best
    speed/ratio balance for this data type: ~10-14% ratio at 300-400 MB/s).

    Extension stacking: fieldName.codec[.zstd][.enc]
    - Examples: U.dvz.zstd, p.dvzt.zstd, U.dvz.zstd.enc

Note
    This class compiles to no-ops if FOAM_USE_ZSTD is not defined.
    Enable with: export FOAM_USE_ZSTD=1 before wmake

SourceFiles
    zstdWrapper.C

\*---------------------------------------------------------------------------*/

#ifndef zstdWrapper_H
#define zstdWrapper_H

#include "fileName.H"
#include "word.H"
#include <vector>
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class zstdWrapper Declaration
\*---------------------------------------------------------------------------*/

class zstdWrapper
{
public:

    // Public constants

        //- File extension for zstd-compressed files
        static const word fileExtension;

        //- Zstd frame magic number
        static constexpr uint32_t MAGIC = 0xFD2FB528;


    // Static Member Functions

        //- Check if zstd support is compiled in
        static bool available();

        //- Compress data buffer with zstd
        static std::vector<uint8_t> compress
        (
            const std::vector<uint8_t>& data,
            int level = 3
        );

        //- Decompress zstd-compressed data buffer
        static std::vector<uint8_t> decompress
        (
            const std::vector<uint8_t>& data
        );

        //- Compress data and write to file
        static void compressToFile
        (
            const fileName& path,
            const std::vector<uint8_t>& data,
            int level = 3
        );

        //- Read and decompress file
        static std::vector<uint8_t> decompressFromFile
        (
            const fileName& path
        );

        //- Check if file starts with zstd magic number
        static bool isZstdFile(const fileName& path);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
