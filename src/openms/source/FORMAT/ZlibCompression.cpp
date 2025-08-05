// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZlibCompression.h>

#include <QtCore/QByteArray>

#include <zlib.h>

using namespace std;

namespace OpenMS
{

  void ZlibCompression::compressString(std::string& str, std::string& compressed)
  {
    compressData(reinterpret_cast<Bytef*>(&str[0]), str.size(), compressed);
  }

  void ZlibCompression::compressData(const void* raw_data, const size_t in_length, std::string& compressed)
  {
    compressed.clear();

    const unsigned long sourceLen = (unsigned long)in_length;
    unsigned long compressed_length =                         // compressBound((unsigned long)str.size());
      sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*

    int zlib_error;
    do
    {
      compressed.resize(compressed_length); // reserve enough space -- we may not need all of it
      zlib_error = compress(reinterpret_cast<Bytef*>(&compressed[0]), &compressed_length, (Bytef*)raw_data, sourceLen);

      switch (zlib_error)
      {
        case Z_MEM_ERROR:
          throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, compressed_length);

        case Z_BUF_ERROR:
          compressed_length *= 2;
      }
    } while (zlib_error == Z_BUF_ERROR);

    if (zlib_error != Z_OK)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compression error?");
    }
    compressed.resize(compressed_length); // cut down to the actual data
  }

  void ZlibCompression::compressString(const QByteArray& raw_data, QByteArray& compressed_data)
  {
    compressed_data = qCompress(raw_data);
    compressed_data.remove(0, 4);
  }

  void ZlibCompression::uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data)
  {
    // Use automatic size detection with streaming decompression
    raw_data.clear();
    
    z_stream strm = {};
    strm.next_in = (Bytef*)compressed_data;
    strm.avail_in = nr_bytes;

    // Initialize zlib for raw deflate format
    if (inflateInit(&strm) != Z_OK)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Failed to initialize zlib decompression");
    }

    // Decompress in chunks
    const size_t chunk_size = 16384;
    std::vector<char> buffer(chunk_size);
    int ret;

    do
    {
      strm.avail_out = chunk_size;
      strm.next_out = (Bytef*)buffer.data();

      ret = inflate(&strm, Z_NO_FLUSH);
      if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR)
      {
        inflateEnd(&strm);
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Zlib decompression error");
      }

      size_t have = chunk_size - strm.avail_out;
      raw_data.append(buffer.data(), have);

    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
  }

  void ZlibCompression::uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data, size_t expected_size)
  {
    if (expected_size > 0)
    {
      // Use size hint for more efficient decompression
      raw_data.resize(expected_size);
      uLongf dest_len = expected_size;
      
      int ret = uncompress((Bytef*)raw_data.data(), &dest_len, (Bytef*)compressed_data, nr_bytes);
      
      if (ret == Z_OK)
      {
        raw_data.resize(dest_len);
      }
      else if (ret == Z_BUF_ERROR)
      {
        // Buffer too small, fall back to streaming method
        uncompressString(compressed_data, nr_bytes, raw_data);
      }
      else
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Zlib decompression error with size hint");
      }
    }
    else
    {
      // No size hint, use streaming method
      uncompressString(compressed_data, nr_bytes, raw_data);
    }
  }

  void ZlibCompression::uncompressString(const QByteArray& compressed_data, QByteArray& raw_data)
  {
    std::string uncompressed_std;
    uncompressString(compressed_data.constData(), compressed_data.size(), uncompressed_std);
    
    raw_data.clear();
    raw_data.append(uncompressed_std.data(), uncompressed_std.size());
  }

}

