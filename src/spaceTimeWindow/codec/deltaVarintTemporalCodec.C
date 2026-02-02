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

#include "deltaVarintTemporalCodec.H"
#include "error.H"
#include <cmath>
#include <cstring>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::deltaVarintTemporalCodec::writeVarint
(
    std::vector<uint8_t>& buf,
    uint64_t val
)
{
    while (val >= 0x80)
    {
        buf.push_back(static_cast<uint8_t>((val & 0x7F) | 0x80));
        val >>= 7;
    }
    buf.push_back(static_cast<uint8_t>(val & 0x7F));
}


uint64_t Foam::deltaVarintTemporalCodec::readVarint
(
    const uint8_t* data,
    size_t& pos,
    size_t size
)
{
    uint64_t result = 0;
    int shift = 0;

    while (pos < size)
    {
        uint8_t byte = data[pos++];
        result |= static_cast<uint64_t>(byte & 0x7F) << shift;

        if ((byte & 0x80) == 0)
        {
            return result;
        }

        shift += 7;

        if (shift >= 64)
        {
            FatalErrorInFunction
                << "Varint overflow while reading compressed data"
                << abort(FatalError);
        }
    }

    FatalErrorInFunction
        << "Unexpected end of data while reading varint"
        << abort(FatalError);

    return 0;
}


void Foam::deltaVarintTemporalCodec::encodeKeyframeComponent
(
    std::vector<uint8_t>& buf,
    const scalar* data,
    label size,
    label stride,
    label offset,
    uint32_t precision
)
{
    if (size <= 0) return;

    const double precisionMult = std::pow(10.0, precision);

    // Write first value as raw double
    double firstVal = data[offset];
    const uint8_t* firstBytes = reinterpret_cast<const uint8_t*>(&firstVal);
    for (int i = 0; i < 8; ++i)
    {
        buf.push_back(firstBytes[i]);
    }

    // Spatial delta-encode remaining values (same as DVZ)
    double prev = firstVal;
    for (label i = 1; i < size; ++i)
    {
        double val = data[offset + i * stride];
        double delta = val - prev;

        // Quantize delta
        int64_t quantized = static_cast<int64_t>(std::round(delta * precisionMult));

        // Zigzag encode and write as varint
        writeVarint(buf, zigzagEncode(quantized));

        // Reconstruct value for next delta (avoids accumulating error)
        prev = prev + static_cast<double>(quantized) / precisionMult;
    }
}


void Foam::deltaVarintTemporalCodec::decodeKeyframeComponent
(
    const uint8_t* buf,
    size_t& pos,
    size_t bufSize,
    scalar* data,
    label size,
    label stride,
    label offset,
    uint32_t precision
)
{
    if (size <= 0) return;

    const double precisionDiv = 1.0 / std::pow(10.0, precision);

    // Read first value as raw double
    if (pos + 8 > bufSize)
    {
        FatalErrorInFunction
            << "Unexpected end of data while reading first value"
            << abort(FatalError);
    }

    double firstVal;
    std::memcpy(&firstVal, buf + pos, 8);
    pos += 8;

    data[offset] = firstVal;

    // Spatial delta-decode remaining values
    double prev = firstVal;
    for (label i = 1; i < size; ++i)
    {
        uint64_t encoded = readVarint(buf, pos, bufSize);
        int64_t quantized = zigzagDecode(encoded);

        double val = prev + static_cast<double>(quantized) * precisionDiv;
        data[offset + i * stride] = val;
        prev = val;
    }
}


void Foam::deltaVarintTemporalCodec::encodeDeltaFrameComponent
(
    std::vector<uint8_t>& buf,
    const scalar* currData,
    const scalar* prevData,
    label size,
    label stride,
    label offset,
    uint32_t precision
)
{
    if (size <= 0) return;

    const double precisionMult = std::pow(10.0, precision);

    // For tracking reconstructed values (to avoid error accumulation)
    double prevReconstructed = 0.0;

    for (label i = 0; i < size; ++i)
    {
        double currVal = currData[offset + i * stride];
        double prevTimeVal = prevData[offset + i * stride];

        // Compute prediction
        double predicted;
        if (i == 0)
        {
            // First face: pure temporal prediction
            predicted = prevTimeVal;
        }
        else
        {
            // Hybrid prediction: spatial + temporal
            // prevReconstructed is the reconstructed value of face i-1 at current time
            // prevTimeVal is the exact value of face i at previous time
            predicted = ALPHA_SPATIAL * prevReconstructed + BETA_TEMPORAL * prevTimeVal;
        }

        // Compute residual
        double residual = currVal - predicted;

        // Quantize residual
        int64_t quantized = static_cast<int64_t>(std::round(residual * precisionMult));

        // Zigzag encode and write as varint
        writeVarint(buf, zigzagEncode(quantized));

        // Reconstruct value for next face's spatial prediction
        prevReconstructed = predicted + static_cast<double>(quantized) / precisionMult;
    }
}


void Foam::deltaVarintTemporalCodec::decodeDeltaFrameComponent
(
    const uint8_t* buf,
    size_t& pos,
    size_t bufSize,
    scalar* currData,
    const scalar* prevData,
    label size,
    label stride,
    label offset,
    uint32_t precision
)
{
    if (size <= 0) return;

    const double precisionDiv = 1.0 / std::pow(10.0, precision);

    double prevReconstructed = 0.0;

    for (label i = 0; i < size; ++i)
    {
        double prevTimeVal = prevData[offset + i * stride];

        // Compute prediction (same as encoder)
        double predicted;
        if (i == 0)
        {
            predicted = prevTimeVal;
        }
        else
        {
            predicted = ALPHA_SPATIAL * prevReconstructed + BETA_TEMPORAL * prevTimeVal;
        }

        // Read and decode residual
        uint64_t encoded = readVarint(buf, pos, bufSize);
        int64_t quantized = zigzagDecode(encoded);

        // Reconstruct value
        double val = predicted + static_cast<double>(quantized) * precisionDiv;
        currData[offset + i * stride] = val;

        // Update for next iteration
        prevReconstructed = val;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deltaVarintTemporalCodec::write
(
    const fileName& path,
    const scalarField& field,
    const scalarField* prevField,
    bool forceKeyframe,
    uint32_t precision
)
{
    std::vector<uint8_t> buf = encode(field, prevField, forceKeyframe, precision);

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    ofs.write(reinterpret_cast<const char*>(buf.data()), buf.size());
}


void Foam::deltaVarintTemporalCodec::write
(
    const fileName& path,
    const vectorField& field,
    const vectorField* prevField,
    bool forceKeyframe,
    uint32_t precision
)
{
    std::vector<uint8_t> buf = encode(field, prevField, forceKeyframe, precision);

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    ofs.write(reinterpret_cast<const char*>(buf.data()), buf.size());
}


std::vector<uint8_t> Foam::deltaVarintTemporalCodec::encode
(
    const scalarField& field,
    const scalarField* prevField,
    bool forceKeyframe,
    uint32_t precision
)
{
    // Determine if this is a keyframe
    bool isKey = forceKeyframe || (prevField == nullptr);

    std::vector<uint8_t> buf;
    buf.reserve(field.size() * 2);

    // Write header (20 bytes)
    // Magic number
    buf.push_back((MAGIC >> 0) & 0xFF);
    buf.push_back((MAGIC >> 8) & 0xFF);
    buf.push_back((MAGIC >> 16) & 0xFF);
    buf.push_back((MAGIC >> 24) & 0xFF);

    // Number of elements
    uint32_t nElements = static_cast<uint32_t>(field.size());
    buf.push_back((nElements >> 0) & 0xFF);
    buf.push_back((nElements >> 8) & 0xFF);
    buf.push_back((nElements >> 16) & 0xFF);
    buf.push_back((nElements >> 24) & 0xFF);

    // Number of components (1 for scalar)
    uint32_t nComponents = 1;
    buf.push_back((nComponents >> 0) & 0xFF);
    buf.push_back((nComponents >> 8) & 0xFF);
    buf.push_back((nComponents >> 16) & 0xFF);
    buf.push_back((nComponents >> 24) & 0xFF);

    // Precision exponent
    buf.push_back((precision >> 0) & 0xFF);
    buf.push_back((precision >> 8) & 0xFF);
    buf.push_back((precision >> 16) & 0xFF);
    buf.push_back((precision >> 24) & 0xFF);

    // Flags
    uint32_t flags = isKey ? FLAG_KEYFRAME : 0;
    buf.push_back((flags >> 0) & 0xFF);
    buf.push_back((flags >> 8) & 0xFF);
    buf.push_back((flags >> 16) & 0xFF);
    buf.push_back((flags >> 24) & 0xFF);

    // Encode data
    if (isKey)
    {
        encodeKeyframeComponent(buf, field.cdata(), field.size(), 1, 0, precision);
    }
    else
    {
        if (prevField->size() != field.size())
        {
            FatalErrorInFunction
                << "Field size mismatch: current=" << field.size()
                << " previous=" << prevField->size()
                << abort(FatalError);
        }
        encodeDeltaFrameComponent(buf, field.cdata(), prevField->cdata(),
                                  field.size(), 1, 0, precision);
    }

    return buf;
}


std::vector<uint8_t> Foam::deltaVarintTemporalCodec::encode
(
    const vectorField& field,
    const vectorField* prevField,
    bool forceKeyframe,
    uint32_t precision
)
{
    bool isKey = forceKeyframe || (prevField == nullptr);

    std::vector<uint8_t> buf;
    buf.reserve(field.size() * 6);

    // Write header
    buf.push_back((MAGIC >> 0) & 0xFF);
    buf.push_back((MAGIC >> 8) & 0xFF);
    buf.push_back((MAGIC >> 16) & 0xFF);
    buf.push_back((MAGIC >> 24) & 0xFF);

    uint32_t nElements = static_cast<uint32_t>(field.size());
    buf.push_back((nElements >> 0) & 0xFF);
    buf.push_back((nElements >> 8) & 0xFF);
    buf.push_back((nElements >> 16) & 0xFF);
    buf.push_back((nElements >> 24) & 0xFF);

    uint32_t nComponents = 3;
    buf.push_back((nComponents >> 0) & 0xFF);
    buf.push_back((nComponents >> 8) & 0xFF);
    buf.push_back((nComponents >> 16) & 0xFF);
    buf.push_back((nComponents >> 24) & 0xFF);

    buf.push_back((precision >> 0) & 0xFF);
    buf.push_back((precision >> 8) & 0xFF);
    buf.push_back((precision >> 16) & 0xFF);
    buf.push_back((precision >> 24) & 0xFF);

    uint32_t flags = isKey ? FLAG_KEYFRAME : 0;
    buf.push_back((flags >> 0) & 0xFF);
    buf.push_back((flags >> 8) & 0xFF);
    buf.push_back((flags >> 16) & 0xFF);
    buf.push_back((flags >> 24) & 0xFF);

    // Component-major ordering: all X, then all Y, then all Z
    const scalar* currData = reinterpret_cast<const scalar*>(field.cdata());

    if (isKey)
    {
        for (uint32_t comp = 0; comp < 3; ++comp)
        {
            encodeKeyframeComponent(buf, currData, field.size(), 3, comp, precision);
        }
    }
    else
    {
        if (prevField->size() != field.size())
        {
            FatalErrorInFunction
                << "Field size mismatch: current=" << field.size()
                << " previous=" << prevField->size()
                << abort(FatalError);
        }

        const scalar* prevData = reinterpret_cast<const scalar*>(prevField->cdata());

        for (uint32_t comp = 0; comp < 3; ++comp)
        {
            encodeDeltaFrameComponent(buf, currData, prevData,
                                      field.size(), 3, comp, precision);
        }
    }

    return buf;
}


Foam::scalarField Foam::deltaVarintTemporalCodec::readScalar
(
    const fileName& path,
    const scalarField* prevField
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

    std::vector<uint8_t> buf(fileSize);
    ifs.read(reinterpret_cast<char*>(buf.data()), fileSize);

    return decodeScalar(buf, prevField);
}


Foam::vectorField Foam::deltaVarintTemporalCodec::readVector
(
    const fileName& path,
    const vectorField* prevField
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

    std::vector<uint8_t> buf(fileSize);
    ifs.read(reinterpret_cast<char*>(buf.data()), fileSize);

    return decodeVector(buf, prevField);
}


Foam::scalarField Foam::deltaVarintTemporalCodec::decodeScalar
(
    const std::vector<uint8_t>& buf,
    const scalarField* prevField
)
{
    if (buf.size() < 20)
    {
        FatalErrorInFunction
            << "Buffer too small to contain header"
            << abort(FatalError);
    }

    size_t pos = 0;

    uint32_t magic = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (magic != MAGIC)
    {
        FatalErrorInFunction
            << "Invalid buffer format (bad magic number, expected DVT1)"
            << abort(FatalError);
    }

    uint32_t nElements = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    uint32_t nComponents = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (nComponents != 1)
    {
        FatalErrorInFunction
            << "Expected scalar field (1 component) but found " << nComponents
            << abort(FatalError);
    }

    uint32_t precision = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    uint32_t flags = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    bool isKey = (flags & FLAG_KEYFRAME) != 0;

    scalarField result(nElements);

    if (isKey)
    {
        decodeKeyframeComponent(buf.data(), pos, buf.size(),
                                result.data(), nElements, 1, 0, precision);
    }
    else
    {
        if (!prevField)
        {
            FatalErrorInFunction
                << "Delta frame requires previous field for decoding"
                << abort(FatalError);
        }
        if (prevField->size() != label(nElements))
        {
            FatalErrorInFunction
                << "Previous field size mismatch"
                << abort(FatalError);
        }
        decodeDeltaFrameComponent(buf.data(), pos, buf.size(),
                                  result.data(), prevField->cdata(),
                                  nElements, 1, 0, precision);
    }

    return result;
}


Foam::vectorField Foam::deltaVarintTemporalCodec::decodeVector
(
    const std::vector<uint8_t>& buf,
    const vectorField* prevField
)
{
    if (buf.size() < 20)
    {
        FatalErrorInFunction
            << "Buffer too small to contain header"
            << abort(FatalError);
    }

    size_t pos = 0;

    uint32_t magic = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (magic != MAGIC)
    {
        FatalErrorInFunction
            << "Invalid buffer format (bad magic number, expected DVT1)"
            << abort(FatalError);
    }

    uint32_t nElements = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    uint32_t nComponents = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (nComponents != 3)
    {
        FatalErrorInFunction
            << "Expected vector field (3 components) but found " << nComponents
            << abort(FatalError);
    }

    uint32_t precision = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    uint32_t flags = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    bool isKey = (flags & FLAG_KEYFRAME) != 0;

    vectorField result(nElements);
    scalar* currData = reinterpret_cast<scalar*>(result.data());

    if (isKey)
    {
        for (uint32_t comp = 0; comp < 3; ++comp)
        {
            decodeKeyframeComponent(buf.data(), pos, buf.size(),
                                    currData, nElements, 3, comp, precision);
        }
    }
    else
    {
        if (!prevField)
        {
            FatalErrorInFunction
                << "Delta frame requires previous field for decoding"
                << abort(FatalError);
        }
        if (prevField->size() != label(nElements))
        {
            FatalErrorInFunction
                << "Previous field size mismatch"
                << abort(FatalError);
        }

        const scalar* prevData = reinterpret_cast<const scalar*>(prevField->cdata());

        for (uint32_t comp = 0; comp < 3; ++comp)
        {
            decodeDeltaFrameComponent(buf.data(), pos, buf.size(),
                                      currData, prevData,
                                      nElements, 3, comp, precision);
        }
    }

    return result;
}


bool Foam::deltaVarintTemporalCodec::isDvztFile(const fileName& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.good())
    {
        return false;
    }

    uint8_t header[4];
    ifs.read(reinterpret_cast<char*>(header), 4);

    if (!ifs.good() || ifs.gcount() != 4)
    {
        return false;
    }

    uint32_t magic = header[0] | (header[1] << 8) | (header[2] << 16) | (header[3] << 24);
    return magic == MAGIC;
}


bool Foam::deltaVarintTemporalCodec::isDvztBuffer(const std::vector<uint8_t>& buf)
{
    if (buf.size() < 4)
    {
        return false;
    }

    uint32_t magic = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
    return magic == MAGIC;
}


bool Foam::deltaVarintTemporalCodec::isKeyframe(const fileName& path)
{
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs.good())
    {
        return false;
    }

    uint8_t header[20];
    ifs.read(reinterpret_cast<char*>(header), 20);

    if (!ifs.good() || ifs.gcount() != 20)
    {
        return false;
    }

    uint32_t magic = header[0] | (header[1] << 8) | (header[2] << 16) | (header[3] << 24);
    if (magic != MAGIC)
    {
        return false;
    }

    // Flags are at offset 16
    uint32_t flags = header[16] | (header[17] << 8) | (header[18] << 16) | (header[19] << 24);
    return (flags & FLAG_KEYFRAME) != 0;
}


bool Foam::deltaVarintTemporalCodec::isKeyframeBuffer(const std::vector<uint8_t>& buf)
{
    if (buf.size() < 20)
    {
        return false;
    }

    uint32_t magic = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
    if (magic != MAGIC)
    {
        return false;
    }

    uint32_t flags = buf[16] | (buf[17] << 8) | (buf[18] << 16) | (buf[19] << 24);
    return (flags & FLAG_KEYFRAME) != 0;
}


// ************************************************************************* //
