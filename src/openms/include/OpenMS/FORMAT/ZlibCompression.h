// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <vector>

class QByteArray;

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
      * @brief Compresses data using Qt
      *
      * @param raw_data Data to be compressed
      * @param compressed_data Compressed result data
      * 
    */
    static void compressString(const QByteArray& raw_data, QByteArray& compressed_data);

    /**
      * @brief Uncompresses data using Qt (wrapper around Qt function)
      *
      * @param compressed_data Compressed data
      * @param nr_bytes Number of bytes in compressed data
      * @param raw_data Uncompressed result data
      * 
    */
    static void uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data);

    /**
      * @brief Uncompresses data using Qt
      *
      * @param compressed_data Compressed data
      * @param raw_data Uncompressed result data
      * 
    */
    static void uncompressString(const QByteArray& compressed_data, QByteArray& raw_data);

  };

} // namespace OpenMS


