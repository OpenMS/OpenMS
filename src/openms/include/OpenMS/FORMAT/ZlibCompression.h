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
    * @brief Compresses and uncompresses data using zlib
    *
    * @note The 'strings' here are not really null-terminated but rather
    * containers of data. If you want safe conversions, use QtByteArray.
    * 
  */
  class OPENMS_DLLAPI ZlibCompression
  {
public:

    /**
      * @brief Compresses data using zlib directly
      *
      * @param raw_data Data to be compressed
      * @param compressed_data Compressed result data
      * 
    */
    static void compressString(std::string& raw_data, std::string& compressed_data);

    /**
     * @brief Compresses data using zlib directly
     *
     * @param raw_data Data to be compressed
     * @param in_length Length of @p raw_data in bytes
     * @param compressed_data Compressed result data
     *
     */
    static void compressData(const void* raw_data, const size_t in_length, std::string& compressed_data);

    /**
      * @brief Uncompresses data using zlib
        
        If available, provide the size of the uncompressed data in @p output_size for a small performance gain.

        @note Does not support gzip format decompression (only zlib format).
       
        @param[in] compressed_data The zlib compressed data
        @param[in] nr_bytes Number of bytes in @p compressed data
        @param[out] out Uncompressed result data
        @param[in] output_size [optional] If known, provide the size of the uncompressed data
        
        @throws Exception::InvalidValue if output_size specified turns out to be smaller than actual size of uncompressed data.
        @throws Exception::InternalToolError if zlib cannot decompress the data (e.g. due to data corruption or unsupported gzip format)
      * 
    */
    static void uncompressData(const void* compressed_data, size_t nr_bytes, std::string& out, size_t output_size = 0);

    /// Convencience function calling @p uncompressData 
    static void uncompressString(const String& in, std::string& out, size_t output_size = 0);

  };

} // namespace OpenMS


