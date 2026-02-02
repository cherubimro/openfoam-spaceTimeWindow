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
    Foam::deltaVarintTemporalCodec

Description
    High-efficiency codec for time-series CFD field data using hybrid
    spatial-temporal prediction:

    - Component-major data ordering (better compression)
    - Hybrid prediction: α * spatial_neighbor + β * temporal_neighbor
    - Quantization to configurable precision (default: 6 significant digits)
    - Zigzag encoding for signed integers
    - Variable-length integer encoding (varint)
    - Keyframe support for random access and error recovery

    This codec achieves better compression than DVZ by exploiting temporal
    correlation in addition to spatial correlation. Expected improvement
    is 30-50% over DVZ for turbulent flow data.

    File format (.dvzt):
    - 4 bytes: magic number "DVT1"
    - 4 bytes: number of faces (uint32)
    - 4 bytes: number of components (uint32)
    - 4 bytes: precision exponent (uint32, e.g., 6 for 1e6)
    - 4 bytes: flags (bit 0 = isKeyframe)
    - For keyframes (same as DVZ):
        - For each component:
            - 8 bytes: first value (double)
            - Variable: spatial delta-encoded remaining values (varint)
    - For delta frames:
        - For each component:
            - Variable: hybrid-predicted residuals (varint)

    Decompression requires:
    - Keyframes: self-contained, no dependencies
    - Delta frames: previous timestep's decoded data

    Typical usage: keyframe every 10-20 timesteps for error recovery.

Usage
    Encoding (during extraction):
    \verbatim
    vectorField prevU = ...;  // Previous timestep (nullptr for keyframe)
    vectorField currU = ...;  // Current timestep
    deltaVarintTemporalCodec::write("U.dvzt", currU, &prevU, isKeyframe, 6);
    \endverbatim

    Decoding (in spaceTimeWindowInitCase):
    \verbatim
    vectorField prevU = ...;  // Previous decoded timestep
    vectorField currU = deltaVarintTemporalCodec::readVector("U.dvzt", &prevU);
    \endverbatim

SourceFiles
    deltaVarintTemporalCodec.C

\*---------------------------------------------------------------------------*/

#ifndef deltaVarintTemporalCodec_H
#define deltaVarintTemporalCodec_H

#include "scalarField.H"
#include "vectorField.H"
#include "fileName.H"
#include <fstream>
#include <vector>
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class deltaVarintTemporalCodec Declaration
\*---------------------------------------------------------------------------*/

class deltaVarintTemporalCodec
{
    // Private static data

        //- Magic number for file identification ("DVT1")
        static constexpr uint32_t MAGIC = 0x31545644;  // "DVT1" in little-endian

        //- Default spatial prediction weight (for hybrid prediction)
        static constexpr double ALPHA_SPATIAL = 0.3;

        //- Default temporal prediction weight (for hybrid prediction)
        static constexpr double BETA_TEMPORAL = 0.7;


    // Private enums

        //- Flag bits
        enum Flags
        {
            FLAG_KEYFRAME = 0x01   // Bit 0: is keyframe
        };


    // Private Member Functions

        //- Zigzag encode a signed 64-bit integer
        static inline uint64_t zigzagEncode(int64_t val)
        {
            return (static_cast<uint64_t>(val) << 1) ^ (val >> 63);
        }

        //- Zigzag decode to signed 64-bit integer
        static inline int64_t zigzagDecode(uint64_t val)
        {
            return static_cast<int64_t>((val >> 1) ^ -(val & 1));
        }

        //- Write a varint to buffer
        static void writeVarint(std::vector<uint8_t>& buf, uint64_t val);

        //- Read a varint from buffer, advancing position
        static uint64_t readVarint(const uint8_t* data, size_t& pos, size_t size);

        //- Encode a keyframe component (spatial-only, same as DVZ)
        static void encodeKeyframeComponent
        (
            std::vector<uint8_t>& buf,
            const scalar* data,
            label size,
            label stride,
            label offset,
            uint32_t precision
        );

        //- Decode a keyframe component
        static void decodeKeyframeComponent
        (
            const uint8_t* buf,
            size_t& pos,
            size_t bufSize,
            scalar* data,
            label size,
            label stride,
            label offset,
            uint32_t precision
        );

        //- Encode a delta frame component (hybrid spatial-temporal prediction)
        static void encodeDeltaFrameComponent
        (
            std::vector<uint8_t>& buf,
            const scalar* currData,      // Current timestep data
            const scalar* prevData,      // Previous timestep data (already decoded)
            label size,
            label stride,
            label offset,
            uint32_t precision
        );

        //- Decode a delta frame component
        static void decodeDeltaFrameComponent
        (
            const uint8_t* buf,
            size_t& pos,
            size_t bufSize,
            scalar* currData,            // Output: current timestep
            const scalar* prevData,      // Input: previous timestep (already decoded)
            label size,
            label stride,
            label offset,
            uint32_t precision
        );


public:

    // Static Member Functions

        // Writing

            //- Write scalar field to compressed file
            static void write
            (
                const fileName& path,
                const scalarField& field,
                const scalarField* prevField,   // nullptr for keyframe
                bool forceKeyframe = false,
                uint32_t precision = 6
            );

            //- Write vector field to compressed file
            static void write
            (
                const fileName& path,
                const vectorField& field,
                const vectorField* prevField,   // nullptr for keyframe
                bool forceKeyframe = false,
                uint32_t precision = 6
            );

        // Encoding (to buffer)

            //- Encode scalar field to buffer
            static std::vector<uint8_t> encode
            (
                const scalarField& field,
                const scalarField* prevField,
                bool forceKeyframe = false,
                uint32_t precision = 6
            );

            //- Encode vector field to buffer
            static std::vector<uint8_t> encode
            (
                const vectorField& field,
                const vectorField* prevField,
                bool forceKeyframe = false,
                uint32_t precision = 6
            );

        // Reading

            //- Read scalar field from compressed file
            static scalarField readScalar
            (
                const fileName& path,
                const scalarField* prevField = nullptr  // Required for delta frames
            );

            //- Read vector field from compressed file
            static vectorField readVector
            (
                const fileName& path,
                const vectorField* prevField = nullptr  // Required for delta frames
            );

        // Decoding (from buffer)

            //- Decode scalar field from buffer
            static scalarField decodeScalar
            (
                const std::vector<uint8_t>& buf,
                const scalarField* prevField = nullptr
            );

            //- Decode vector field from buffer
            static vectorField decodeVector
            (
                const std::vector<uint8_t>& buf,
                const vectorField* prevField = nullptr
            );

        // Utilities

            //- Check if file is in delta-varint-temporal format (by magic number)
            static bool isDvztFile(const fileName& path);

            //- Check if buffer is in delta-varint-temporal format
            static bool isDvztBuffer(const std::vector<uint8_t>& buf);

            //- Check if file/buffer is a keyframe (can be decoded standalone)
            static bool isKeyframe(const fileName& path);
            static bool isKeyframeBuffer(const std::vector<uint8_t>& buf);

            //- Get file extension for compressed files
            static word fileExtension()
            {
                return "dvzt";
            }

            //- Get magic number (for format detection)
            static uint32_t magic()
            {
                return MAGIC;
            }

            //- Get prediction weights
            static double alphaSpatial() { return ALPHA_SPATIAL; }
            static double betaTemporal() { return BETA_TEMPORAL; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
