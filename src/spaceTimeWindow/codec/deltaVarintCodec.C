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

#include "deltaVarintCodec.H"
#include "error.H"
#include <cmath>
#include <cstring>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::deltaVarintCodec::writeVarint
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


uint64_t Foam::deltaVarintCodec::readVarint
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


void Foam::deltaVarintCodec::encodeComponent
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

    // Delta-encode remaining values
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


void Foam::deltaVarintCodec::decodeComponent
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

    // Delta-decode remaining values
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::deltaVarintCodec::write
(
    const fileName& path,
    const scalarField& field,
    uint32_t precision
)
{
    std::vector<uint8_t> buf;
    buf.reserve(field.size() * 2);  // Rough estimate

    // Write header
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

    // Encode data (single component, stride=1, offset=0)
    encodeComponent(buf, field.cdata(), field.size(), 1, 0, precision);

    // Write to file
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    ofs.write(reinterpret_cast<const char*>(buf.data()), buf.size());
}


void Foam::deltaVarintCodec::write
(
    const fileName& path,
    const vectorField& field,
    uint32_t precision
)
{
    std::vector<uint8_t> buf;
    buf.reserve(field.size() * 6);  // Rough estimate (3 components, ~2 bytes each)

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

    // Number of components (3 for vector)
    uint32_t nComponents = 3;
    buf.push_back((nComponents >> 0) & 0xFF);
    buf.push_back((nComponents >> 8) & 0xFF);
    buf.push_back((nComponents >> 16) & 0xFF);
    buf.push_back((nComponents >> 24) & 0xFF);

    buf.push_back((precision >> 0) & 0xFF);
    buf.push_back((precision >> 8) & 0xFF);
    buf.push_back((precision >> 16) & 0xFF);
    buf.push_back((precision >> 24) & 0xFF);

    // Component-major ordering: all X, then all Y, then all Z
    // This gives better compression as similar values are grouped
    const scalar* data = reinterpret_cast<const scalar*>(field.cdata());

    for (uint32_t comp = 0; comp < 3; ++comp)
    {
        // For vector field laid out as [x0,y0,z0, x1,y1,z1, ...]
        // Component c has stride=3, offset=c
        encodeComponent(buf, data, field.size(), 3, comp, precision);
    }

    // Write to file
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs.good())
    {
        FatalErrorInFunction
            << "Cannot open file for writing: " << path
            << abort(FatalError);
    }

    ofs.write(reinterpret_cast<const char*>(buf.data()), buf.size());
}


Foam::scalarField Foam::deltaVarintCodec::readScalar(const fileName& path)
{
    // Read file into buffer
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

    // Read header
    if (fileSize < 16)
    {
        FatalErrorInFunction
            << "File too small to contain header: " << path
            << abort(FatalError);
    }

    size_t pos = 0;

    // Check magic number
    uint32_t magic = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (magic != MAGIC)
    {
        FatalErrorInFunction
            << "Invalid file format (bad magic number): " << path
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
            << " components in file: " << path
            << abort(FatalError);
    }

    uint32_t precision = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    // Decode data
    scalarField result(nElements);
    decodeComponent(buf.data(), pos, fileSize, result.data(), nElements, 1, 0, precision);

    return result;
}


Foam::vectorField Foam::deltaVarintCodec::readVector(const fileName& path)
{
    // Read file into buffer
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

    // Read header
    if (fileSize < 16)
    {
        FatalErrorInFunction
            << "File too small to contain header: " << path
            << abort(FatalError);
    }

    size_t pos = 0;

    uint32_t magic = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    if (magic != MAGIC)
    {
        FatalErrorInFunction
            << "Invalid file format (bad magic number): " << path
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
            << " components in file: " << path
            << abort(FatalError);
    }

    uint32_t precision = buf[pos] | (buf[pos+1] << 8) | (buf[pos+2] << 16) | (buf[pos+3] << 24);
    pos += 4;

    // Decode data (component-major to interleaved)
    vectorField result(nElements);
    scalar* data = reinterpret_cast<scalar*>(result.data());

    for (uint32_t comp = 0; comp < 3; ++comp)
    {
        decodeComponent(buf.data(), pos, fileSize, data, nElements, 3, comp, precision);
    }

    return result;
}


std::vector<uint8_t> Foam::deltaVarintCodec::encode
(
    const scalarField& field,
    uint32_t precision
)
{
    std::vector<uint8_t> buf;
    buf.reserve(field.size() * 2);

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

    uint32_t nComponents = 1;
    buf.push_back((nComponents >> 0) & 0xFF);
    buf.push_back((nComponents >> 8) & 0xFF);
    buf.push_back((nComponents >> 16) & 0xFF);
    buf.push_back((nComponents >> 24) & 0xFF);

    buf.push_back((precision >> 0) & 0xFF);
    buf.push_back((precision >> 8) & 0xFF);
    buf.push_back((precision >> 16) & 0xFF);
    buf.push_back((precision >> 24) & 0xFF);

    encodeComponent(buf, field.cdata(), field.size(), 1, 0, precision);

    return buf;
}


std::vector<uint8_t> Foam::deltaVarintCodec::encode
(
    const vectorField& field,
    uint32_t precision
)
{
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

    const scalar* data = reinterpret_cast<const scalar*>(field.cdata());

    for (uint32_t comp = 0; comp < 3; ++comp)
    {
        encodeComponent(buf, data, field.size(), 3, comp, precision);
    }

    return buf;
}


Foam::scalarField Foam::deltaVarintCodec::decodeScalar
(
    const std::vector<uint8_t>& buf
)
{
    if (buf.size() < 16)
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
            << "Invalid buffer format (bad magic number)"
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

    scalarField result(nElements);
    decodeComponent(buf.data(), pos, buf.size(), result.data(), nElements, 1, 0, precision);

    return result;
}


Foam::vectorField Foam::deltaVarintCodec::decodeVector
(
    const std::vector<uint8_t>& buf
)
{
    if (buf.size() < 16)
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
            << "Invalid buffer format (bad magic number)"
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

    vectorField result(nElements);
    scalar* data = reinterpret_cast<scalar*>(result.data());

    for (uint32_t comp = 0; comp < 3; ++comp)
    {
        decodeComponent(buf.data(), pos, buf.size(), data, nElements, 3, comp, precision);
    }

    return result;
}


bool Foam::deltaVarintCodec::isDeltaVarintFile(const fileName& path)
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


bool Foam::deltaVarintCodec::isDeltaVarintBuffer(const std::vector<uint8_t>& buf)
{
    if (buf.size() < 4)
    {
        return false;
    }

    uint32_t magic = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
    return magic == MAGIC;
}


// ************************************************************************* //
