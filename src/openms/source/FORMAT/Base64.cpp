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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Base64.h>

#include <QtCore/QList>
#include <QtCore/QString>

using namespace std;

namespace OpenMS
{

  /*

   Background in the following two encoding / decoding mapping arrays.

   While encoding we have to map a binary value to its character value using
   the base 64 mapping:

    binary  ->    char = val

       0    ->     A   = 65  
                   ...     
      25    ->     Z   = 90  
      26    ->     a   = 97  
                   ...     
      51    ->     z   = 122 
      52    ->     0   = 48  
                   ...     
      61    ->     9   = 57  
      62    ->     +   = 43  
      63    ->     /   = 47  
      
   
   While decoding we have to map a character to its base 64 target using the
   base 64 mapping:

    char = val      ->      target 
   
     A   = 65       ->        0
     ...     
     Z   = 90       ->       25
     a   = 97       ->       26
     ...     
     z   = 122      ->       51
     0   = 48       ->       52
     ...     
     9   = 57       ->       61
     +   = 43       ->       62
     /   = 47       ->       63


  this can be done by first subtracting 43 from each value
  and then looking up a value in a table as found in the code,
  e.g. lookup[char - 43] - 62

  The following string can be produced by this Python snippet:

offset = 62
s = ""
s += chr(62 + offset)
s += '$' * 3
s += "".join([ chr(52 + i + offset) for i in range(10) ])
s += '$' * 7
s += "".join([ chr(i + offset) for i in range(26) ])
s += '$' * 6
s += "".join([ chr(26 + i + offset) for i in range(26) ])
print s

"|$$$rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\]^_`abcdefghijklmnopq"

  however, this is basically just convenience to produce printable characters.
  We could also go for a direct mapping:

s = '$' * 42
s += hex(62)
s += '$' * 3
s += "".join([ hex(52 + i) for i in range(10) ])
s += '$' * 7
s += "".join([ hex(i) for i in range(26) ])
s += '$' * 6
s += "".join([ hex(26 + i) for i in range(26) ])
print s

"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$0x3e$$$0x340x350x360x370x380x390x3a0x3b0x3c0x3d$$$$$$$0x00x10x20x30x40x50x60x70x80x90xa0xb0xc0xd0xe0xf0x100x110x120x130x140x150x160x170x180x19$$$$$$0x1a0x1b0x1c0x1d0x1e0x1f0x200x210x220x230x240x250x260x270x280x290x2a0x2b0x2c0x2d0x2e0x2f0x300x310x320x33"

  which would allow lookup[char] without any adding or subtraction step (but
  needs a 42 byte extra memory in the local cache).

  */

  const char Base64::encoder_[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  const char Base64::decoder_[] = "|$$$}rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\\]^_`abcdefghijklmnopq";

  void Base64::encodeStrings(const std::vector<String>& in, String& out, bool zlib_compression, bool append_null_byte)
  {
    out.clear();
    if (in.empty())
      return;

    std::string str;
    std::string compressed;
    Byte* it;
    Byte* end;
    for (Size i = 0; i < in.size(); ++i)
    {
      str = str.append(in[i]);
      if (append_null_byte) str.push_back('\0');
    }

    if (zlib_compression)
    {
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

      it = reinterpret_cast<Byte*>(&compressed[0]);
      end = it + compressed_length;
      // TODO check integer overflow
      out.resize((Size)ceil(compressed_length / 3.) * 4); //resize output array in order to have enough space for all characters
    }
    else
    {
      // TODO check integer overflow
      out.resize((Size)ceil(str.size() / 3.) * 4); //resize output array in order to have enough space for all characters
      it = reinterpret_cast<Byte*>(&str[0]);
      end = it + str.size();
    }
    Byte* to = reinterpret_cast<Byte*>(&out[0]);
    Size written = 0;

    while (it != end)
    {
      Int int_24bit = 0;
      Int padding_count = 0;

      // construct 24-bit integer from 3 bytes
      for (Size i = 0; i < 3; i++)
      {
        if (it != end)
        {
          int_24bit |= *it++ << ((2 - i) * 8);
        }
        else
        {
          padding_count++;
        }
      }

      // write out 4 characters
      for (Int i = 3; i >= 0; i--)
      {
        to[i] = encoder_[int_24bit & 0x3F];
        int_24bit >>= 6;
      }

      // fixup for padding
      if (padding_count > 0)
        to[3] = '=';
      if (padding_count > 1)
        to[2] = '=';

      to += 4;
      written += 4;
    }

    out.resize(written); //no more space is needed
  }

  void Base64::decodeStrings(const String& in, std::vector<String>& out, bool zlib_compression)
  {
    out.clear();

    // The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }

    QByteArray base64_uncompressed;
    decodeSingleString(in, base64_uncompressed, zlib_compression);
    QList<QByteArray> null_strings = base64_uncompressed.split('\0');
    for (QList<QByteArray>::iterator it = null_strings.begin(); it < null_strings.end(); ++it)
    {
      if (!it->isEmpty())
      {
        out.push_back(QString(*it).toStdString());
      }
    }
  }

  void Base64::decodeSingleString(const String& in, QByteArray& base64_uncompressed, bool zlib_compression)
  {
    // The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }

    QByteArray herewego = QByteArray::fromRawData(in.c_str(), (int) in.size());
    base64_uncompressed = QByteArray::fromBase64(herewego);
    if (zlib_compression)
    {
      QByteArray czip;
      czip.resize(4);
      czip[0] = (base64_uncompressed.size() & 0xff000000) >> 24;
      czip[1] = (base64_uncompressed.size() & 0x00ff0000) >> 16;
      czip[2] = (base64_uncompressed.size() & 0x0000ff00) >> 8;
      czip[3] = (base64_uncompressed.size() & 0x000000ff);
      czip += base64_uncompressed;
      base64_uncompressed = qUncompress(czip);

      if (base64_uncompressed.isEmpty())
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression error?");
      }
    }
  }

} //end OpenMS
