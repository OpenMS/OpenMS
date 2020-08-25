// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZlibCompression.h>

#include <zlib.h>

using namespace std;

namespace OpenMS
{

  void ZlibCompression::compressString(std::string& str, std::string& compressed)
  {
    compressed.clear();

    unsigned long sourceLen =   (unsigned long)str.size();
    unsigned long compressed_length = //compressBound((unsigned long)str.size());
                                      sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*

    int zlib_error;
    do
    {
      compressed.resize(compressed_length);
      zlib_error = compress(reinterpret_cast<Bytef*>(&compressed[0]), &compressed_length, reinterpret_cast<Bytef*>(&str[0]), (unsigned long) str.size());

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
    compressed.resize(compressed_length);
  }

  void ZlibCompression::compressString(const QByteArray& raw_data, QByteArray& compressed_data)
  {
    compressed_data = qCompress(raw_data);
    compressed_data.remove(0, 4);
  }

  void ZlibCompression::uncompressString(const void * tt, size_t blob_bytes, std::string& uncompressed)
  {
    // take a leap of faith and assume the input is valid
    QByteArray compressed_data = QByteArray::fromRawData((const char*)tt, blob_bytes);
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

