// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZlibCompression.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <array>
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
    unsigned long compressed_length = compressBound(in_length);
    compressed.resize(compressed_length); // reserve enough space -- we may not need all of it

    int zlib_error = compress(reinterpret_cast<Bytef*>(&compressed[0]), &compressed_length, (Bytef*)raw_data, sourceLen);

    switch (zlib_error)
    {
      case Z_OK:         // Compression successful
        break;
      case Z_MEM_ERROR:
        throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, compressed_length);
      case Z_BUF_ERROR:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                      "Destination buffer too small for compressed output",
                                      String(compressed_length));
      default:
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compression error!");
    }

    compressed.resize(compressed_length); // cut down to the actual data
  }

  void uncompressStringSingleShot___(const void* compressed_data, size_t nr_bytes, std::string& raw_data, size_t output_size)
  {
    if (nr_bytes == 0)
    {
      raw_data.clear();
      return;
    }
    raw_data.resize(output_size);
    uLongf uncompressedSize = output_size;
    int ret = uncompress((Bytef*)raw_data.data(), &uncompressedSize, (Bytef*)compressed_data, nr_bytes);

    if (ret == Z_OK)
    {
      if (uncompressedSize != raw_data.size())
      {
        OPENMS_LOG_INFO << "ZlibCompression::uncompressString: Warning: decompressed data was smaller (" << std::to_string(uncompressedSize) << ") than anticipated: " << std::to_string(output_size) << std::endl;
        raw_data.resize(uncompressedSize);
      }
    }
    else if (ret == Z_BUF_ERROR)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression failed because specified output_size was too small. Size of data after decompression is larger than anticipated: " + std::to_string(uncompressedSize), std::to_string(output_size));
    }
    else
    {
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "zlib::inflate failed with code " + std::to_string(ret) + " .");
    }

  }

  void uncompressStringIterative___(const void* compressed_data, size_t nr_bytes, std::string& uncompressed)
  {
    uncompressed.clear();

    // Setup input
    z_stream strm = {};
    strm.next_in = (Bytef*)(compressed_data);
    strm.avail_in = nr_bytes;

    // Initialize zlib (use inflateInit2 for gzip or raw deflate)
    int ret = inflateInit(&strm);
    if (ret != Z_OK)
    {
      throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "zlib::inflateInit failed with code " + std::to_string(ret) + " .");
    }

    // Decompress loop
    const size_t CHUNK_SIZE = 16384;
    std::array<char, CHUNK_SIZE> buffer;
    do
    {
      strm.avail_out = CHUNK_SIZE;
      strm.next_out = (Bytef*)buffer.data();

      ret = inflate(&strm, Z_NO_FLUSH);
      if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR || ret == Z_BUF_ERROR)
      {
        inflateEnd(&strm);
        throw Exception::InternalToolError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "zlib::inflate failed with code " + std::to_string(ret) + " .");
      }

      size_t bytes_decompressed = CHUNK_SIZE - strm.avail_out;
      uncompressed.insert(uncompressed.end(), buffer.begin(), buffer.begin() + bytes_decompressed);

    } while (ret != Z_STREAM_END);

    inflateEnd(&strm);
  }

  
  void ZlibCompression::uncompressData(const void* compressed_data, size_t nr_bytes, std::string& raw_data, size_t output_size)
  {
    if (output_size == 0)
    { // unknown output size - decompress chunk by chunk
      uncompressStringIterative___(compressed_data, nr_bytes, raw_data);
    }
    else
    { // known output size - decompress into a preallocated string
      uncompressStringSingleShot___(compressed_data, nr_bytes, raw_data, output_size);
    }
  }

  void ZlibCompression::uncompressString(const String& in, std::string& out, size_t output_size)
  {
    uncompressData(in.data(), in.size(), out, output_size);
  }
  
} // namespace OpenMS
