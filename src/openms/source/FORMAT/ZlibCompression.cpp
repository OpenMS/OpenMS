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

  void ZlibCompression::uncompressString(const void * tt, size_t blob_bytes, std::string& uncompressed)
  {
    // take a leap of faith and assume the input is valid
    QByteArray compressed_data = QByteArray::fromRawData((const char*)tt, (int)blob_bytes);
    QByteArray raw_data;

    ZlibCompression::uncompressString(compressed_data, raw_data);

    // Note that we may have zero bytes in the string, so we cannot use QString
    uncompressed.clear();
    uncompressed = std::string(raw_data.data(), raw_data.size());
  }

  void ZlibCompression::uncompressString(const QByteArray& compressed_data, QByteArray& raw_data)
  {
    QByteArray czip;
    czip.resize(4);
    czip[0] = (compressed_data.size() & 0xff000000) >> 24;
    czip[1] = (compressed_data.size() & 0x00ff0000) >> 16;
    czip[2] = (compressed_data.size() & 0x0000ff00) >> 8;
    czip[3] = (compressed_data.size() & 0x000000ff);
    czip += compressed_data;
    raw_data = qUncompress(czip);

    if (raw_data.isEmpty())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression error?");
    }
  }

}

