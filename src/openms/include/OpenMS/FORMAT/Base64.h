// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_FORMAT_BASE64_H
#define OPENMS_FORMAT_BASE64_H

#ifndef OPENMS_IS_BIG_ENDIAN
#if defined OPENMS_BIG_ENDIAN
#define OPENMS_IS_BIG_ENDIAN true
#else
#define OPENMS_IS_BIG_ENDIAN false
#endif
#endif

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <vector>

#include <QByteArray>
#include <zlib.h>

namespace OpenMS
{
  /**
    @brief Class to encode and decode Base64

    Base64 supports two precisions: 32 bit (float) and 64 bit (double).
  */
  class OPENMS_DLLAPI Base64
  {

public:

    /// default constructor
    Base64();

    /// Destructor
    virtual ~Base64();

    /// Byte order type
    enum ByteOrder
    {
      BYTEORDER_BIGENDIAN,                  ///< Big endian type
      BYTEORDER_LITTLEENDIAN            ///< Little endian type
    };

    /**
        @brief Encodes a vector of floating point numbers to a Base64 string

        You can specify the byte order of the output and if it is zlib-compressed.

        @note @p in will be empty after this method
    */
    template <typename FromType>
    void encode(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a vector of floating point numbers

        You have to specify the byte order of the input and if it is zlib-compressed.
    */
    template <typename ToType>
    void decode(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression = false);

    /**
        @brief Encodes a vector of integer point numbers to a Base64 string

        You can specify the byte order of the output and if it is zlib-compressed.

        @note @p in will be empty after this method
    */
    template <typename FromType>
    void encodeIntegers(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a vector of integer numbers

        You have to specify the byte order of the input and if it is zlib-compressed.
    */
    template <typename ToType>
    void decodeIntegers(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression = false);

    /**
        @brief Encodes a vector of strings to a Base64 string

        You can specify zlib-compression.

        @param in A vector of data to be encoded (String)
        @param out A String containing the Base64 encoded data
        @param zlib_compression Whether the data should be compressed with zlib before encoding in Base64
        @param append_null_byte Whether a null-byte should be appended after each of the Strings contained in the in vector
      
        @note Unless append_null_byte is false, will add a null byte ("\0") at the end of each input
    */
    void encodeStrings(const std::vector<String> & in, String & out, bool zlib_compression = false, bool append_null_byte = true);

    /**
        @brief Decodes a Base64 string to a vector of (null-terminated) strings

        You have to specify whether the Base64 string is zlib-compressed.

        @param in A String containing the Base64 encoded data
        @param out A vector containing the decoded data (split at null "\0") bytes
        @param zlib_compression Whether the data should be decompressed with zlib after decoding in Base64
    */
    void decodeStrings(const String & in, std::vector<String> & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a QByteArray

        @param in A String containing the Base64 encoded data
        @param out A ByteArray containing the decoded data
        @param zlib_compression Whether the data should be decompressed with zlib after decoding in Base64
    */
    void decodeSingleString(const String & in, QByteArray & base64_uncompressed, bool zlib_compression);

private:

    ///Internal class needed for type-punning
    union Reinterpreter64_
    {
      double f;
      UInt64 i;
    };

    ///Internal class needed for type-punning
    union Reinterpreter32_
    {
      float f;
      UInt32 i;
    };

    static const char encoder_[];
    static const char decoder_[];
    /// Decodes a Base64 string to a vector of floating point numbers
    template <typename ToType>
    void decodeUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    ///Decodes a compressed Base64 string to a vector of floating point numbers
    template <typename ToType>
    void decodeCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    /// Decodes a Base64 string to a vector of integer numbers
    template <typename ToType>
    void decodeIntegersUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    ///Decodes a compressed Base64 string to a vector of integer numbers
    template <typename ToType>
    void decodeIntegersCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);
  };

  /// Endianizes a 32 bit type from big endian to little endian and vice versa
  inline UInt32 endianize32(const UInt32& n)
  {
    return ((n & 0x000000ff) << 24) | 
           ((n & 0x0000ff00) <<  8) |
           ((n & 0x00ff0000) >>  8) |
           ((n & 0xff000000) >> 24);
  }

  /// Endianizes a 64 bit type from  big endian to little endian and vice versa
  inline UInt64 endianize64(const UInt64& n)
  {
    return ((n >> 56) & 0x00000000000000FF) |
           ((n >> 40) & 0x000000000000FF00) | 
           ((n >> 24) & 0x0000000000FF0000) | 
           ((n >>  8) & 0x00000000FF000000) |
           ((n <<  8) & 0x000000FF00000000) |
           ((n << 24) & 0x0000FF0000000000) |
           ((n << 40) & 0x00FF000000000000) | 
           ((n << 56) & 0xFF00000000000000);
  }

  template <typename FromType>
  void Base64::encode(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression)
  {
    out.clear();
    if (in.empty())
      return;

    //initialize
    const Size element_size = sizeof(FromType);
    const Size input_bytes = element_size * in.size();
    String compressed;
    Byte * it;
    Byte * end;
    //Change endianness if necessary
    if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      if (element_size == 4)
      {
        for (Size i = 0; i < in.size(); ++i)
        {
          Reinterpreter32_ tmp;
          tmp.f = in[i];
          tmp.i = endianize32(tmp.i);
          in[i] = tmp.f;
        }
      }
      else
      {
        for (Size i = 0; i < in.size(); ++i)
        {
          Reinterpreter64_ tmp;
          tmp.f = static_cast<double>(in[i]);
          tmp.i = endianize64(tmp.i);
          in[i] = tmp.f;
        }
      }
    }

    //encode with compression
    if (zlib_compression)
    {
      unsigned long sourceLen =   (unsigned long)in.size();
      unsigned long compressed_length =       //compressBound((unsigned long)in.size());
                                        sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*
      //
      // (*) compressBound is not defined in the QtCore lib, which forces the linker under windows to link in our zlib.
      //     This leads to multiply defined symbols as compress() is then defined twice.

      int zlib_error;
      do
      {
        compressed.resize(compressed_length);
        zlib_error = compress(reinterpret_cast<Bytef *>(&compressed[0]), &compressed_length, reinterpret_cast<Bytef *>(&in[0]), (unsigned long)input_bytes);

        switch (zlib_error)
        {
        case Z_MEM_ERROR:
          throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, compressed_length);
          break;

        case Z_BUF_ERROR:
          compressed_length *= 2;
        }
      }
      while (zlib_error == Z_BUF_ERROR);

      if (zlib_error != Z_OK)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compression error?");
      }

      String(compressed).swap(compressed);
      it = reinterpret_cast<Byte *>(&compressed[0]);
      end = it + compressed_length;
      out.resize((Size)ceil(compressed_length / 3.) * 4);     //resize output array in order to have enough space for all characters
    }
    //encode without compression
    else
    {
      out.resize((Size)ceil(input_bytes / 3.) * 4);     //resize output array in order to have enough space for all characters
      it = reinterpret_cast<Byte *>(&in[0]);
      end = it + input_bytes;
    }

    Byte * to = reinterpret_cast<Byte *>(&out[0]);


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

    out.resize(written);         //no more space is needed
  }

  template <typename ToType>
  void Base64::decode(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression)
  {
    if (zlib_compression)
    {
      decodeCompressed_(in, from_byte_order, out);
    }
    else
    {
      decodeUncompressed_(in, from_byte_order, out);
    }
  }

  template <typename ToType>
  void Base64::decodeCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();
    if (in == "") return;

    const Size element_size = sizeof(ToType);

    String decompressed;

    QByteArray qt_byte_array = QByteArray::fromRawData(in.c_str(), (int) in.size());
    QByteArray bazip = QByteArray::fromBase64(qt_byte_array);
    QByteArray czip;
    czip.resize(4);
    czip[0] = (bazip.size() & 0xff000000) >> 24;
    czip[1] = (bazip.size() & 0x00ff0000) >> 16;
    czip[2] = (bazip.size() & 0x0000ff00) >> 8;
    czip[3] = (bazip.size() & 0x000000ff);
    czip += bazip;
    QByteArray base64_uncompressed = qUncompress(czip);

    if (base64_uncompressed.isEmpty())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression error?");
    }
    decompressed.resize(base64_uncompressed.size());

    std::copy(base64_uncompressed.begin(), base64_uncompressed.end(), decompressed.begin());

    void* byte_buffer = reinterpret_cast<void *>(&decompressed[0]);
    Size buffer_size = decompressed.size();

    const ToType * float_buffer = reinterpret_cast<const ToType *>(byte_buffer);
    if (buffer_size % element_size != 0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount?");
    }
    
    Size float_count = buffer_size / element_size;
    
    // change endianness if necessary
    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      if (element_size == 4) // 32 bit
      {
        UInt32 * p = reinterpret_cast<UInt32 *>(byte_buffer);
        std::transform(p, p + float_count, p, endianize32);
      }
      else // 64 bit
      {
        UInt64 * p = reinterpret_cast<UInt64 *>(byte_buffer);
        std::transform(p, p + float_count, p, endianize64);
      }
    }

    // copy values
    out.assign(float_buffer, float_buffer + float_count);
  }

  template <typename ToType>
  void Base64::decodeUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();

    // The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }
    if (in.size() % 4 != 0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Malformed base64 input, length is not a multiple of 4.");
    }

    Size src_size = in.size();
    // last one or two '=' are skipped if contained
    int padding = 0;
    if (in[src_size - 1] == '=') padding++;
    if (in[src_size - 2] == '=') padding++;

    src_size -= padding;

    register UInt a;
    register UInt b;

    UInt offset = 0;
    int inc = 1;
    UInt written = 0;

    const Size element_size = sizeof(ToType);

    // enough for either float or double
    char element[8] = "\x00\x00\x00\x00\x00\x00\x00";

    // Parse little endian data in big endian OpenMS (or other way round)
    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || 
       (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      offset = (element_size - 1);              // other endian
      inc = -1;
    }
    else
    {
      offset = 0;
      inc = 1;
    }

    //reserve enough space in the output vector
    out.reserve((UInt)(std::ceil((4.0 * src_size) / 3.0) + 6.0));

    // sort all read bytes correctly into a char[4] (double) or
    // char[8] (float) and push_back when necessary.
    for (Size i = 0; i < src_size; i += 4)
    {
      // decode 4 Base64-Chars to 3 Byte
      // -------------------------------

      // decode the first two chars
      a = decoder_[(int)in[i] - 43] - 62;
      b = decoder_[(int)in[i + 1] - 43] - 62;
      if (i + 1 >= src_size)
      {
        b = 0;
      }
      // write first byte (6 bits from a and 2 highest bits from b)
      element[offset] = (unsigned char) ((a << 2) | (b >> 4));
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType * to_type = reinterpret_cast<ToType *>(&element[0]);
        out.push_back((*to_type));
        strcpy(element, "");
      }

      // decode the third char
      a = decoder_[(int)in[i + 2] - 43] - 62;
      if (i + 2 >= src_size)
      {
        a = 0;
      }
      // write second byte (4 lowest bits from b and 4 highest bits from a)
      element[offset] = (unsigned char) (((b & 15) << 4) | (a >> 2));
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType * to_type = reinterpret_cast<ToType *>(&element[0]);
        out.push_back((*to_type));
        strcpy(element, "");
      }

      // decode the fourth char
      b = decoder_[(int)in[i + 3] - 43] - 62;
      if (i + 3 >= src_size)
      {
        b = 0;
      }
      // write third byte (2 lowest bits from a and 6 bits from b)
      element[offset] = (unsigned char) (((a & 3) << 6) | b);
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType * to_type = reinterpret_cast<ToType *>(&element[0]);
        out.push_back((*to_type));
        strcpy(element, "");
      }
    }
  }

  template <typename FromType>
  void Base64::encodeIntegers(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression)
  {
    out.clear();
    if (in.empty())
      return;

    //initialize
    const Size element_size = sizeof(FromType);
    const Size input_bytes = element_size * in.size();
    String compressed;
    Byte * it;
    Byte * end;
    //Change endianness if necessary
    if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      if (element_size == 4)
      {
        for (Size i = 0; i < in.size(); ++i)
        {
          UInt32 tmp = in[i];
          tmp = endianize32(tmp);
          in[i] = tmp;
        }
      }
      else
      {
        for (Size i = 0; i < in.size(); ++i)
        {
          UInt64 tmp = in[i];
          tmp = endianize64(tmp);
          in[i] = tmp;
        }
      }
    }

    //encode with compression (use Qt because of zlib support)
    if (zlib_compression)
    {
      unsigned long sourceLen =   (unsigned long)input_bytes;
      unsigned long compressed_length =       //compressBound((unsigned long)in.size());
                                        sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*

      compressed.resize(compressed_length);
      while (compress(reinterpret_cast<Bytef *>(&compressed[0]), &compressed_length, reinterpret_cast<Bytef *>(&in[0]), (unsigned long)input_bytes) != Z_OK)
      {
        compressed_length *= 2;
        compressed.reserve(compressed_length);
      }


      String(compressed).swap(compressed);
      it = reinterpret_cast<Byte *>(&compressed[0]);
      end = it + compressed_length;
      out.resize((Size)ceil(compressed_length / 3.) * 4);     //resize output array in order to have enough space for all characters
    }
    //encode without compression
    else
    {
      out.resize((Size)ceil(input_bytes / 3.) * 4);     //resize output array in order to have enough space for all characters
      it = reinterpret_cast<Byte *>(&in[0]);
      end = it + input_bytes;
    }

    Byte * to = reinterpret_cast<Byte *>(&out[0]);


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

    out.resize(written);         //no more space is needed
  }

  template <typename ToType>
  void Base64::decodeIntegers(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression)
  {
    if (zlib_compression)
    {
      decodeIntegersCompressed_(in, from_byte_order, out);
    }
    else
    {
      decodeIntegersUncompressed_(in, from_byte_order, out);
    }
  }

  template <typename ToType>
  void Base64::decodeIntegersCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();
    if (in == "")
      return;

    void * byte_buffer;
    Size buffer_size;
    const Size element_size = sizeof(ToType);

    String decompressed;

    QByteArray qt_byte_array = QByteArray::fromRawData(in.c_str(), (int) in.size());
    QByteArray bazip = QByteArray::fromBase64(qt_byte_array);
    QByteArray czip;
    czip.resize(4);
    czip[0] = (bazip.size() & 0xff000000) >> 24;
    czip[1] = (bazip.size() & 0x00ff0000) >> 16;
    czip[2] = (bazip.size() & 0x0000ff00) >> 8;
    czip[3] = (bazip.size() & 0x000000ff);
    czip += bazip;
    QByteArray base64_uncompressed = qUncompress(czip);
    if (base64_uncompressed.isEmpty())
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression error?");
    }
    decompressed.resize(base64_uncompressed.size());

    std::copy(base64_uncompressed.begin(), base64_uncompressed.end(), decompressed.begin());

    byte_buffer = reinterpret_cast<void *>(&decompressed[0]);
    buffer_size = decompressed.size();

    //change endianness if necessary
    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      if (element_size == 4)
      {
        const Int32 * float_buffer = reinterpret_cast<const Int32 *>(byte_buffer);
        if (buffer_size % element_size != 0)
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount?");
        Size float_count = buffer_size / element_size;
        UInt32 * p = reinterpret_cast<UInt32 *>(byte_buffer);
        std::transform(p, p + float_count, p, endianize32);

        out.resize(float_count);
        // do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
        for (Size i = 0; i < float_count; ++i)
        {
          out[i] = (ToType) * float_buffer;
          ++float_buffer;
        }
      }
      else
      {
        const Int64 * float_buffer = reinterpret_cast<const Int64 *>(byte_buffer);

        if (buffer_size % element_size != 0)
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount?");

        Size float_count = buffer_size / element_size;

        UInt64 * p = reinterpret_cast<UInt64 *>(byte_buffer);
        std::transform(p, p + float_count, p, endianize64);

        out.resize(float_count);
        // do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
        for (Size i = 0; i < float_count; ++i)
        {
          out[i] = (ToType) * float_buffer;
          ++float_buffer;
        }
      }
    }
    else
    {
      if (element_size == 4)
      {
        const Int * float_buffer = reinterpret_cast<const Int *>(byte_buffer);
        if (buffer_size % element_size != 0)
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount while decoding?");

        Size float_count = buffer_size / element_size;
        out.resize(float_count);
        // do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
        for (Size i = 0; i < float_count; ++i)
        {
          out[i] = (ToType) * float_buffer;
          ++float_buffer;
        }
      }
      else
      {
        const Int64 * float_buffer = reinterpret_cast<const Int64 *>(byte_buffer);

        if (buffer_size % element_size != 0)
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Bad BufferCount while decoding?");

        Size float_count = buffer_size / element_size;
        out.resize(float_count);
        // do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
        for (Size i = 0; i < float_count; ++i)
        {
          out[i] = (ToType) * float_buffer;
          ++float_buffer;
        }
      }
    }

  }

  template <typename ToType>
  void Base64::decodeIntegersUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();

    // The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }

    Size src_size = in.size();
    // last one or two '=' are skipped if contained
    int padding = 0;
    if (in[src_size - 1] == '=') padding++;
    if (in[src_size - 2] == '=') padding++;

    src_size -= padding;

    register UInt a;
    register UInt b;

    UInt offset = 0;
    int inc = 1;
    UInt written = 0;

    const Size element_size = sizeof(ToType);

    // enough for either float or double
    char element[8] = "\x00\x00\x00\x00\x00\x00\x00";

    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      offset = (element_size - 1);              // other endian
      inc = -1;
    }
    else
    {
      offset = 0;
      inc = 1;
    }

    //reserve enough space in the output vector
    out.reserve((UInt)(std::ceil((4.0 * src_size) / 3.0) + 6.0));

    // sort all read bytes correctly into a char[4] (double) or
    // char[8] (float) and push_back when necessary.
    for (Size i = 0; i < src_size; i += 4)
    {

      // decode 4 Base64-Chars to 3 Byte
      // -------------------------------

      // decode the first two chars
      a = decoder_[(int)in[i] - 43] - 62;
      b = decoder_[(int)in[i + 1] - 43] - 62;
      if (i + 1 >= src_size)
      {
        b = 0;
      }
      // write first byte (6 bits from a and 2 highest bits from b)
      element[offset] = (unsigned char) ((a << 2) | (b >> 4));
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType float_value;
        if (element_size == 4)
        {
          Int32 * value = reinterpret_cast<Int32 *>(&element[0]);
          float_value = (ToType) * value;
        }
        else
        {
          Int64 * value = reinterpret_cast<Int64 *>(&element[0]);
          float_value = (ToType) * value;
        }
        out.push_back(float_value);
        strcpy(element, "");
      }

      // decode the third char
      a = decoder_[(int)in[i + 2] - 43] - 62;
      if (i + 2 >= src_size)
      {
        a = 0;
      }
      // write second byte (4 lowest bits from b and 4 highest bits from a)
      element[offset] = (unsigned char) (((b & 15) << 4) | (a >> 2));
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType float_value;
        if (element_size == 4)
        {
          Int32 * value = reinterpret_cast<Int32 *>(&element[0]);
          float_value = (ToType) * value;
        }
        else
        {
          Int64 * value = reinterpret_cast<Int64 *>(&element[0]);
          float_value = (ToType) * value;
        }
        out.push_back(float_value);
        strcpy(element, "");
      }

      // decode the fourth char
      b = decoder_[(int)in[i + 3] - 43] - 62;
      if (i + 3 >= src_size)
      {
        b = 0;
      }
      // write third byte (2 lowest bits from a and 6 bits from b)
      element[offset] = (unsigned char) (((a & 3) << 6) | b);
      written++;
      offset = (offset + inc) % element_size;

      if (written % element_size == 0)
      {
        ToType float_value;
        if (element_size == 4)
        {
          Int32 * value = reinterpret_cast<Int32 *>(&element[0]);
          float_value = (ToType) * value;
        }
        else
        {
          Int64 * value = reinterpret_cast<Int64 *>(&element[0]);
          float_value = (ToType) * value;
        }
        out.push_back(float_value);
        strcpy(element, "");
      }
    }
  }

} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
