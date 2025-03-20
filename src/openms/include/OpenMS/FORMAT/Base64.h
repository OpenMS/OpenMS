// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Moritz Aubermann $
// --------------------------------------------------------------------------

#pragma once

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
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <QtCore/QByteArray>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#ifdef OPENMS_COMPILER_MSVC
#pragma comment(linker, "/export:compress")
#endif

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
    Base64() = default;

    /// Byte order type
    enum ByteOrder
    {
      BYTEORDER_BIGENDIAN,              ///< Big endian type
      BYTEORDER_LITTLEENDIAN            ///< Little endian type
    };
	
    /**
        @brief Encodes a vector of floating point numbers to a Base64 string

        You can specify the byte order of the output and if it is zlib-compressed.

        @note @p in will be empty after this method
    */
    template <typename FromType>
    static void encode(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a vector of floating point numbers

        You have to specify the byte order of the input and if it is zlib-compressed.
    */
    template <typename ToType>
    static void decode(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression = false);

    /**
        @brief Encodes a vector of integer point numbers to a Base64 string

        You can specify the byte order of the output and if it is zlib-compressed.

        @note @p in will be empty after this method
    */
    template <typename FromType>
    static void encodeIntegers(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a vector of integer numbers

        You have to specify the byte order of the input and if it is zlib-compressed.
    */
    template <typename ToType>
    static void decodeIntegers(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out, bool zlib_compression = false);

    /**
        @brief Encodes a vector of strings to a Base64 string

        You can specify zlib-compression.

        @param in A vector of data to be encoded (String)
        @param out A String containing the Base64 encoded data
        @param zlib_compression Whether the data should be compressed with zlib before encoding in Base64
        @param append_null_byte Whether a null-byte should be appended after each of the Strings contained in the in vector
      
        @note Unless append_null_byte is false, will add a null byte ("\0") at the end of each input
    */
    static void encodeStrings(const std::vector<String> & in, String & out, bool zlib_compression = false, bool append_null_byte = true);

    /**
        @brief Decodes a Base64 string to a vector of (null-terminated) strings

        You have to specify whether the Base64 string is zlib-compressed.

        @param in A String containing the Base64 encoded data
        @param out A vector containing the decoded data (split at null "\0") bytes
        @param zlib_compression Whether the data should be decompressed with zlib after decoding in Base64
    */
    static void decodeStrings(const String & in, std::vector<String> & out, bool zlib_compression = false);

    /**
        @brief Decodes a Base64 string to a QByteArray

        @param in A String containing the Base64 encoded data
        @param base64_uncompressed A ByteArray containing the decoded data
        @param zlib_compression Whether the data should be decompressed with zlib after decoding in Base64
    */
    static void decodeSingleString(const String& in, QByteArray& base64_uncompressed, bool zlib_compression);

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
    static void decodeUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    ///Decodes a compressed Base64 string to a vector of floating point numbers
    template <typename ToType>
    static void decodeCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    /// Decodes a Base64 string to a vector of integer numbers
    template <typename ToType>
    static void decodeIntegersUncompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    ///Decodes a compressed Base64 string to a vector of integer numbers
    template <typename ToType>
    static void decodeIntegersCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out);

    static void stringSimdEncoder_(std::string& in, std::string& out);

    static void stringSimdDecoder_(const std::string& in, std::string& out);
  };

  // Possible optimization: add simd registerwise endianizer (this will only be beneficial for ARM, since mzML + x64 CPU does not need to convert since both use LITTLE_ENDIAN).
  // mzXML(!), which is outdated uses BIG_ENDIAN, i.e. "network", in its base64 encoding, so there x64 will benefit, but not ARM.
  // However: the code below gets optimized to the bswap instruction by most compilers, which is very fast (1 cycle latency + 1 ops)
  // and it is doubtful that SSE4's _mm_shuffle_epi8 will do better, see https://dev.to/wunk/fast-array-reversal-with-simd-j3p
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
    {
      return;
    }

    // initialize
    const Size element_size = sizeof(FromType);
    const Size input_bytes = element_size * in.size();
    // change endianness if necessary
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
      /////////////////////////endianize registerwise
    }

    // encode with compression
    if (zlib_compression)
    {
      String compressed;
      ZlibCompression::compressData((void*)in.data(), input_bytes, compressed);
      stringSimdEncoder_(compressed, out);
    }
    else // encode without compression
    {
      String str((char*)in.data(), input_bytes);
      stringSimdEncoder_(str, out);
    }

  }

  template <typename ToType>  ////////////////////////////////////////////nothing to change here, magic happenes elsewhere
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

  template <int type_size>
  inline void invertEndianess(void* byte_buffer, const size_t element_count);
  template<>
  inline void invertEndianess<4>(void* byte_buffer, const size_t element_count)
  {
    UInt32* p = reinterpret_cast<UInt32*>(byte_buffer);
    std::transform(p, p + element_count, p, endianize32);
  }
  template<>
  inline void invertEndianess<8>(void* byte_buffer, const size_t element_count)
  {
    UInt64* p = reinterpret_cast<UInt64*>(byte_buffer);
    std::transform(p, p + element_count, p, endianize64);
  }


  template <typename ToType>
  void Base64::decodeCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();
    if (in.empty()) return;

    constexpr Size element_size = sizeof(ToType);

    String decompressed;

    String s;
    stringSimdDecoder_(in, s);
    QByteArray bazip = QByteArray::fromRawData(s.c_str(), (int) s.size());

   /////////////////////////////////////////////////////////////////////////////////////if faster: first encode then call fromRawData
   // QByteArray qt_byte_array = QByteArray::fromRawData(in.c_str(), (int) in.size());
   // QByteArray bazip = QByteArray::fromBase64(qt_byte_array);
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
      invertEndianess<element_size>(byte_buffer, float_count);
    }

    // copy values
    out.assign(float_buffer, float_buffer + float_count);
  }

  template <typename ToType>
  void Base64::decodeUncompressed_(const String& in, ByteOrder from_byte_order , std::vector<ToType>& out)
  {
    out.clear();

    // The length of a base64 string is always a multiple of 4 (always 3
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

    constexpr Size element_size = sizeof(ToType);
    String s;
    stringSimdDecoder_(in,s);

    // change endianness if necessary (mzML is always LITTLE_ENDIAN; x64 is LITTLE_ENDIAN)
    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      invertEndianess<element_size>((void*)s.data(), s.size() / element_size);
    }

    const char* cptr = s.data();
    const ToType * fptr = reinterpret_cast<const ToType*>(cptr);
    out.assign(fptr,fptr + s.size()/element_size);
  }

  template <typename FromType>
  void Base64::encodeIntegers(std::vector<FromType> & in, ByteOrder to_byte_order, String & out, bool zlib_compression)
  {
    out.clear();
    if (in.empty())
      return;

    // initialize
    const Size element_size = sizeof(FromType);
    const Size input_bytes = element_size * in.size();

    // change endianness if necessary
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

    // encode with compression (use Qt because of zlib support)
    if (zlib_compression)
    {
      String compressed;
      ZlibCompression::compressData((void*)in.data(), input_bytes, compressed);
      stringSimdEncoder_(compressed, out);
    }
    else // encode without compression
    {
      String str((char*)in.data(), input_bytes);
      stringSimdEncoder_(str, out);
    }
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
    if (in.empty())
      return;

    void * byte_buffer;
    Size buffer_size;
    constexpr Size element_size = sizeof(ToType);

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

    // change endianness if necessary
    if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
    {
      if constexpr(element_size == 4)
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
      if constexpr(element_size == 4)
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

    UInt a;
    UInt b;

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

    // reserve enough space in the output vector
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

