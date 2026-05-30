// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class String;

  /**
    @brief Compresses and uncompresses arbitrary byte buffers using zlib.

    Static utility class. All methods treat the @c std::string / @c String
    arguments as raw byte containers; the data is not assumed to be NUL-
    terminated and may contain embedded zeros.

    @note Only the raw zlib format is supported on the decompression side;
          gzip-framed streams are rejected.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ZlibCompression
  {
public:

    /**
      @brief Compress @p raw_data into @p compressed_data using zlib.

      @param[in]  raw_data        Data to compress (treated as a raw byte buffer; the reference is read-only despite the non-const type).
      @param[out] compressed_data Receives the compressed payload; any previous contents are replaced.

      @throws Exception::OutOfMemory     if zlib cannot allocate an output buffer of the worst-case size.
      @throws Exception::InvalidValue    if the destination buffer turns out to be too small for the compressed output.
      @throws Exception::ConversionError if zlib reports any other failure during compression.
    */
    static void compressString(std::string& raw_data, std::string& compressed_data);

    /**
      @brief Compress the @p in_length bytes pointed to by @p raw_data into @p compressed_data using zlib.

      @param[in]  raw_data        Pointer to the bytes to compress.
      @param[in]  in_length       Length of @p raw_data in bytes.
      @param[out] compressed_data Receives the compressed payload; any previous contents are replaced.

      @throws Exception::OutOfMemory     if zlib cannot allocate an output buffer of the worst-case size.
      @throws Exception::InvalidValue    if the destination buffer turns out to be too small for the compressed output.
      @throws Exception::ConversionError if zlib reports any other failure during compression.
    */
    static void compressData(const void* raw_data, const size_t in_length, std::string& compressed_data);

    /**
      @brief Uncompress @p compressed_data using zlib.

      When the uncompressed size is known up front, pass it as @p output_size
      to skip the iterative chunked path and decompress straight into a
      single preallocated buffer. Pass @c 0 (the default) when the size is
      not known; decompression then proceeds chunk by chunk.

      @note Only zlib-format data is accepted; gzip-framed streams are
            rejected with @c Exception::InternalToolError.

      @param[in]  compressed_data Pointer to the zlib-compressed bytes.
      @param[in]  nr_bytes        Length of @p compressed_data in bytes; passing @c 0 with a non-zero @p output_size results in an empty @p out.
      @param[out] out             Receives the decompressed bytes; any previous contents are replaced.
      @param[in]  output_size     Expected uncompressed size (optional). When @c 0, an iterative path is taken; when positive, @p out is sized to @p output_size up front. If the actual decompressed size turns out to be smaller, @p out is shrunk to match and an informational message is logged.

      @throws Exception::InvalidValue      if @p output_size was specified (>0) but is smaller than the actual uncompressed size.
      @throws Exception::InternalToolError if zlib cannot decompress the data (corrupted input, unsupported gzip framing, or any other zlib failure).
    */
    static void uncompressData(const void* compressed_data, size_t nr_bytes, std::string& out, size_t output_size = 0);

    /**
      @brief Convenience wrapper around @ref uncompressData that takes a @c String input.

      @param[in]  in          Compressed input.
      @param[out] out         Receives the decompressed bytes; any previous contents are replaced.
      @param[in]  output_size Expected uncompressed size; see @ref uncompressData.

      @throws Exception::InvalidValue      see @ref uncompressData.
      @throws Exception::InternalToolError see @ref uncompressData.
    */
    static void uncompressString(const String& in, std::string& out, size_t output_size = 0);

  };

} // namespace OpenMS


