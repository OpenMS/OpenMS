// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
// $Authors: Marc Sturm, Moritz Aubermann $
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

#include <QByteArray>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>


#include <simde/x86/ssse3.h>
#include <simde/x86/sse2.h>
#include <simde/x86/sse.h>

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
      BYTEORDER_BIGENDIAN,                  ///< Big endian type
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

    bool testfunction(std::string s){
      std::string t;
      std::string a;
      Base64 unit;
      unit.stringSimdEncoder_(s,t);
      unit.stringSimdDecoder_(t,a);
      std::cout << a << std::endl;
      if(s==a)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

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

    
    const simde__m128i mask1_ = simde_mm_set1_epi32(0x3F000000);//00111111 00000000 00000000 00000000
    const simde__m128i mask2_ = simde_mm_set1_epi32(0x003F0000);//00000000 00111111 00000000 00000000
    const simde__m128i mask3_ = simde_mm_set1_epi32(0x00003F00);//00000000 00000000 00111111 00000000
    const simde__m128i mask4_ = simde_mm_set1_epi32(0x0000003F);//00000000 00000000 00000000 00111111

    const simde__m128i mask1d_ = simde_mm_set1_epi32(0xFF000000); // 11111111 00000000 00000000 00000000
    const simde__m128i mask2d_ = simde_mm_set1_epi32(0x00FF0000); // 00000000 11111111 00000000 00000000
    const simde__m128i mask3d_ = simde_mm_set1_epi32(0x0000FF00); // 00000000 00000000 11111111 00000000
    const simde__m128i mask4d_ = simde_mm_set1_epi32(0x000000FF); // 00000000 00000000 00000000 11111111

    //difference between base64 encoding and ascii encoding, used to cast from base64 binaries to characters
    
    const simde__m128i difference_A_ = simde_mm_set1_epi8('A');
    const simde__m128i difference_a_ = simde_mm_set1_epi8('a' - 26);
    const simde__m128i difference_0_ = simde_mm_set1_epi8('0' - 52);
    const simde__m128i difference_plus_ = simde_mm_set1_epi8('+');
    const simde__m128i difference_slash_ = simde_mm_set1_epi8('/');



    const simde__m128i shuffle_mask_1_ = simde_mm_setr_epi8(2, 2, 1, 0, 5, 5, 4, 3, 8, 8, 7, 6, 11, 11, 10, 9);
    const simde__m128i shuffle_mask_2_ = simde_mm_setr_epi8(3,2,1,0,7,6,5,4,11,10,9,8,15,14,13,12);

    const simde__m128i shuffle_mask_big_endian_ = simde_mm_setr_epi8(0,1,2,2,3,4,5,5,6,7,8,8,9,10,11,11);
    //second shuffle doesnt need to happen

    //decoding shuffle masks:
    //shuffle_mask_2 gets used
    const simde__m128i shuffle_mask_d_2_ = simde_mm_setr_epi8(3,2,1,7,6,5,11,10,9,15,14,13,0,4,8,12);
    

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

    inline void registerEncoder_(simde__m128i & data);

    inline void registerDecoder_(simde__m128i & data);

    inline void stringSimdEncoder_(std::string & in, std::string & out);

    inline void stringSimdDecoder_(const std::string & in, std::string & out);
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
  //TODO: add simd registerwise endianizer

  // these operators are defined for GCC/clang, but not in MSVC (TODO: maybe use SFINAE, but that is overkill for the moment)
  #ifdef _MSC_VER 
  inline simde__m128i operator|(const simde__m128i& left, const simde__m128i& right)
  {
    return simde_mm_or_si128(left, right);
  }
  inline simde__m128i& operator|=(simde__m128i& left, const simde__m128i& right)
  {
    left = simde_mm_or_si128(left, right);
    return left;
  }
  inline simde__m128i operator&(const simde__m128i left, const simde__m128i& right)
  {
    return simde_mm_and_si128(left, right);
  }
  #endif

  ///encode the first 12 bytes of a 128 bit simde integer type to base64
  void Base64::registerEncoder_(simde__m128i &data)
  {
    
    if(! OPENMS_IS_BIG_ENDIAN)
    {
      data = simde_mm_shuffle_epi8(data, shuffle_mask_1_);
    //by shuffling every 3 8bit ASCII Letters now take up 4 bytes, "ABC" gets shuffled to "CCBA" to match the 4 bytes of the Base64 Encoding, and deal with little Endianness.
    }
    else
    {
      data = simde_mm_shuffle_epi8(data, shuffle_mask_big_endian_);
    }
    //shifting and masking data, so now every 6 bit of the ASCII encoding have their own byte. 
    //shifting over 32 bits takes endianness into account, which needs to be accounted for when shuffeling
    data =   (simde_mm_srli_epi32(data, 2 ) & mask1_) |
             (simde_mm_srli_epi32(data, 4 ) & mask2_) |
             (simde_mm_srli_epi32(data, 6 ) & mask3_) |
             (data & mask4_);

    if(! OPENMS_IS_BIG_ENDIAN)  //otherwise the data already is ordered correctly
    {      
      data = simde_mm_shuffle_epi8(data, shuffle_mask_2_);
    }

    //masking data and adding/substracting to match base64 codes to fitting characters
    simde__m128i capital_mask = simde_mm_cmplt_epi8(data, _mm_set1_epi8(26)); // (a < b) ? 0xFF : 0x00
    simde__m128i all_mask= capital_mask;
    simde__m128i lower_case_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data,simde_mm_set1_epi8(52))); //not allMask and  b where b is 0xFF if binaries are smaller than 52
    all_mask |= lower_case_mask;
    simde__m128i number_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(62)));
    all_mask |= number_mask;
    simde__m128i plus_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(63)));
    all_mask |= plus_mask;
    simde__m128i & slash_negative_mask = all_mask;

    data =  (capital_mask   & simde_mm_add_epi8(data, difference_A_))|
            (lower_case_mask & simde_mm_add_epi8(data, difference_a_))|
            (number_mask      & simde_mm_add_epi8(data, difference_0_))|
            (plus_mask        & difference_plus_                      )|
            (simde_mm_andnot_si128(slash_negative_mask, difference_slash_));

  } 

  void Base64::registerDecoder_(simde__m128i & data)
  {
    //ASCII letters must be translated over to base64. This cannot be achieved by just adding/substracting a single value
    //from each letter since in ASCII Alphabet capital letters aren't followed up by small Letters, and small Letters not by numbers (..).

    //Therefore certain kinds of characters must be masked out for further processing:

    //plusMask equals 0xFF for each corresponding plus, otherwise its 0
    simde__m128i plusMask = simde_mm_cmpeq_epi8(data, difference_plus_);
    simde__m128i allMask = plusMask;
    //slashMask similar to plusMask
    simde__m128i slashMask = simde_mm_cmpeq_epi8(data, difference_slash_);
    allMask |= slashMask;
    //for number mask: all characters less than '9' plus 1 must be numbers, '+' or '/' because input is Base64
    //therefore "not allMask and less than '9' + 1 (see allNumbers in header) " applied on data sets all bytes corresponding to numbers in the mask to 0xFF
    simde__m128i numberMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('9' + 1)));
    allMask |= numberMask;
    simde__m128i bigLetterMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('Z'+1)));
    allMask |= bigLetterMask;
    simde__m128i smallLetterMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('z'+ 1)));

    //match ASCII characters with coresponding Base64 codes:
    data =  plusMask    & simde_mm_set1_epi8(62)|
            slashMask   & simde_mm_set1_epi8(63)|
            numberMask  & simde_mm_add_epi8( data, simde_mm_set1_epi8(4))|
            bigLetterMask & simde_mm_sub_epi8( data, simde_mm_set1_epi8(65))| //ASCII 'A' is 65, Base64 'A' is 0
            smallLetterMask & simde_mm_sub_epi8(data, simde_mm_set1_epi8(71));

    //convert little endian to big endian:
    data = simde_mm_shuffle_epi8(data, shuffle_mask_2_);

    //the actual magic (conversion base64 to ASCII) happens here by shifting and masking:
    data =  simde_mm_slli_epi32((data & mask1d_),2)   |
            simde_mm_slli_epi32((data & mask2d_),4)   |
            simde_mm_slli_epi32((data & mask3d_),6)   |
            simde_mm_slli_epi32((data & mask4d_),8);

    //convert big endian to little endian
    data= simde_mm_shuffle_epi8(data, shuffle_mask_d_2_);

  }

  void Base64::stringSimdEncoder_(std::string & in, std::string & out)
  {
          // TODO check integer overflow
      out.resize((Size)(in.size() / 3) * 4 + 16); //resize output array, so the register encoder doesnt write memory to unallocated memory
      uint8_t padding = (3 - in.size() % 3 ) % 3;
      const int loop = in.size() / 12;
     
      in.resize(in.size() + 4, '\0');
      //otherwise there are cases where register encoder isnt allowed to access last bytes   
          
      Base64 unit;
      simde__m128i data{};
      //loop  through input as long as it's safe to access memory
      for(int i = 0; i < loop; i++)
      {
          //each time the last 4 out of 16 byte string data get lost through processing, therefore jumps of 12 bytes (/characters)
          data = simde_mm_lddqu_si128((simde__m128i*) & in[12*i] );
          unit.registerEncoder_(data);
          simde_mm_storeu_si128((simde__m128*) & out[i*16], data);
      }
      
      size_t read = loop *12;
      size_t written = loop * 16;

      //create buffer to translate last bytes without accessing memory that hasn't been allocated 
      std::array<char,16> buffer{};
      memcpy(& buffer[0],& in[read],in.size()-read -4); //minus 4 because of 4 appended null bytes
      data = simde_mm_lddqu_si128((simde__m128i*) & buffer[0]);
      unit.registerEncoder_(data);
      simde_mm_storeu_si128((simde__m128*) & out[written], data);
      
      in.resize(in.size() -4); //remove null bytes

      //resizing out and add padding if necessary
      if(padding) 
      {
        size_t newsize = ceil((double)in.size() / 3.) * 4; 
        out.resize(newsize );
        for(size_t j = newsize - 1; j >= newsize -padding; j--)
        {
            out[j] = '=';
        }      
      }
      else
      {        
          out.resize((in.size() / 3) * 4);
      }
    
  }

  void Base64::stringSimdDecoder_(const std::string & in, std::string & out)
  {
    out.clear();
    //out.resize((in.size() / 4) * 3);
    const char* inPtr = &in[0];

    //padding count:
    uint8_t g = 0;
    if( in[in.size() -1] == '=')
        g++;
    if( in[in.size() -2] == '=')
        g++;

    unsigned outsize = (in.size() / 16) * 12 + 16;
    //not final size (final rezize later to cutoff unwanted characters)
    out.resize(outsize);
    char* outPtr = &out[0];
    int loop= in.size()  / 16;

    for(int i = 0; i < loop; i++)
    {
        simde__m128i data = simde_mm_lddqu_si128((simde__m128i*) (inPtr + i*16) );
        registerDecoder_(data);
        simde_mm_storeu_si128((simde__m128 *) (outPtr + i * 12), data);
    }

    size_t read = loop * 16;
    std::array<char, 16> rest;
    std::fill(rest.begin(), rest.end(), 'x');
    std::copy(in.begin() + read, in.end(), rest.begin());

    simde__m128i data = simde_mm_lddqu_si128((simde__m128i*) &rest[0] );
    registerDecoder_(data);
    size_t written = loop * 12;
    simde_mm_storeu_si128((simde__m128*)(outPtr + written), data);

    //cutting off decoding of appendix
    outsize = (in.size() / 4) * 3 - g;
    out.resize(outsize);
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
      /////////////////////////endianize registerwise
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
      Base64 unit;
      unit.stringSimdEncoder_(compressed, out);   //resize output array in order to have enough space for all characters
    }
    //encode without compression
    else
    {
      String str((char*)in.data(), input_bytes);
      Base64 unit;
      unit.stringSimdEncoder_(str, out);

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

  template <typename ToType>
  void Base64::decodeCompressed_(const String & in, ByteOrder from_byte_order, std::vector<ToType> & out)
  {
    out.clear();
    if (in.empty()) return;

    const Size element_size = sizeof(ToType);

    String decompressed;

//  
    String s;
    Base64 unit;
    unit.stringSimdDecoder_(in,s);
    QByteArray bazip = QByteArray::fromRawData(s.c_str(), (int) s.size());
//
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

      ////////////////////////////////use faster register endianizer
    }

    // copy values
    out.assign(float_buffer, float_buffer + float_count);
  }

  template <typename ToType>
  void Base64::decodeUncompressed_(const String& in, ByteOrder /*from_byte_order*/ , std::vector<ToType>& out)  // TODO byte order not needed?
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

    const Size element_size = sizeof(ToType);
    Base64 unit;
    String s;
    unit.stringSimdDecoder_(in,s);
    const char* cptr = s.data();
    const ToType * fptr = reinterpret_cast<const ToType*>(cptr);
    out.assign(fptr,fptr+s.size()/element_size);
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

      Base64 unit;
      unit.stringSimdEncoder_(compressed,out);      
    }
    //encode without compression
    else
    {
      String str((char*)in.data(), input_bytes);
      Base64 unit;
      unit.stringSimdEncoder_(str,out);

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

