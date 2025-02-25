// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Moritz Aubermann $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Base64.h>

#include <QtCore/QList>
#include <QtCore/QString>

#include <OpenMS/SYSTEM/SIMDe.h>

using namespace std;

namespace OpenMS
{
  const simde__m128i mask1_ = simde_mm_set1_epi32(0x3F000000); // 00111111 00000000 00000000 00000000
  const simde__m128i mask2_ = simde_mm_set1_epi32(0x003F0000); // 00000000 00111111 00000000 00000000
  const simde__m128i mask3_ = simde_mm_set1_epi32(0x00003F00); // 00000000 00000000 00111111 00000000
  const simde__m128i mask4_ = simde_mm_set1_epi32(0x0000003F); // 00000000 00000000 00000000 00111111

  const simde__m128i mask1d_ = simde_mm_set1_epi32(0xFF000000); // 11111111 00000000 00000000 00000000
  const simde__m128i mask2d_ = simde_mm_set1_epi32(0x00FF0000); // 00000000 11111111 00000000 00000000
  const simde__m128i mask3d_ = simde_mm_set1_epi32(0x0000FF00); // 00000000 00000000 11111111 00000000
  const simde__m128i mask4d_ = simde_mm_set1_epi32(0x000000FF); // 00000000 00000000 00000000 11111111

  // difference between base64 encoding and ascii encoding, used to cast from base64 binaries to characters

  const simde__m128i difference_A_ = simde_mm_set1_epi8('A');
  const simde__m128i difference_a_ = simde_mm_set1_epi8('a' - 26);
  const simde__m128i difference_0_ = simde_mm_set1_epi8('0' - 52);
  const simde__m128i difference_plus_ = simde_mm_set1_epi8('+');
  const simde__m128i difference_slash_ = simde_mm_set1_epi8('/');

  const simde__m128i shuffle_mask_1_ = simde_mm_setr_epi8(2, 2, 1, 0, 5, 5, 4, 3, 8, 8, 7, 6, 11, 11, 10, 9);
  const simde__m128i shuffle_mask_2_ = simde_mm_setr_epi8(3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12);

  const simde__m128i shuffle_mask_big_endian_ = simde_mm_setr_epi8(0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11);
  // second shuffle doesnt need to happen

  // decoding shuffle masks:
  // shuffle_mask_2 gets used
  const simde__m128i shuffle_mask_d_2_ = simde_mm_setr_epi8(3, 2, 1, 7, 6, 5, 11, 10, 9, 15, 14, 13, 0, 4, 8, 12);

  /// Encode the first 12 bytes of a 128 bit simde integer type to base64
  void registerEncoder_(simde__m128i& data)
  {
    if constexpr (!OPENMS_IS_BIG_ENDIAN)
    {
      data = simde_mm_shuffle_epi8(data, shuffle_mask_1_);
      // by shuffling every 3 8bit ASCII Letters now take up 4 bytes, "ABC" gets shuffled to "CCBA" to match the 4 bytes of the Base64 Encoding, and deal with little Endianness.
    }
    else
    {
      data = simde_mm_shuffle_epi8(data, shuffle_mask_big_endian_);
    }
    // shifting and masking data, so now every 6 bit of the ASCII encoding have their own byte.
    // shifting over 32 bits takes endianness into account, which needs to be accounted for when shuffeling
    data = (simde_mm_srli_epi32(data, 2) & mask1_) | (simde_mm_srli_epi32(data, 4) & mask2_) | (simde_mm_srli_epi32(data, 6) & mask3_) | (data & mask4_);

    if constexpr (!OPENMS_IS_BIG_ENDIAN) // otherwise the data already is ordered correctly
    {
      data = simde_mm_shuffle_epi8(data, shuffle_mask_2_);
    }

    // masking data and adding/substracting to match base64 codes to fitting characters
    simde__m128i capital_mask = simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(26)); // (a < b) ? 0xFF : 0x00
    simde__m128i all_mask = capital_mask;
    simde__m128i lower_case_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(52))); // not allMask and  b where b is 0xFF if binaries are smaller than 52
    all_mask |= lower_case_mask;
    simde__m128i number_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(62)));
    all_mask |= number_mask;
    simde__m128i plus_mask = simde_mm_andnot_si128(all_mask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8(63)));
    all_mask |= plus_mask;
    simde__m128i& slash_negative_mask = all_mask;

    data = (capital_mask & simde_mm_add_epi8(data, difference_A_)) | (lower_case_mask & simde_mm_add_epi8(data, difference_a_)) | (number_mask & simde_mm_add_epi8(data, difference_0_)) |
           (plus_mask & difference_plus_) | (simde_mm_andnot_si128(slash_negative_mask, difference_slash_));
  }

  void registerDecoder_(simde__m128i& data)
  {
    // ASCII letters must be translated over to base64. This cannot be achieved by just adding/substracting a single value
    // from each letter since in ASCII Alphabet capital letters aren't followed up by small Letters, and small Letters not by numbers (..).

    // Therefore certain kinds of characters must be masked out for further processing:

    // plusMask equals 0xFF for each corresponding plus, otherwise its 0
    simde__m128i plusMask = simde_mm_cmpeq_epi8(data, difference_plus_);
    simde__m128i allMask = plusMask;
    // slashMask similar to plusMask
    simde__m128i slashMask = simde_mm_cmpeq_epi8(data, difference_slash_);
    allMask |= slashMask;
    // for number mask: all characters less than '9' plus 1 must be numbers, '+' or '/' because input is Base64
    // therefore "not allMask and less than '9' + 1 (see allNumbers in header) " applied on data sets all bytes corresponding to numbers in the mask to 0xFF
    simde__m128i numberMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('9' + 1)));
    allMask |= numberMask;
    simde__m128i bigLetterMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('Z' + 1)));
    allMask |= bigLetterMask;
    simde__m128i smallLetterMask = simde_mm_andnot_si128(allMask, simde_mm_cmplt_epi8(data, simde_mm_set1_epi8('z' + 1)));

    // match ASCII characters with coresponding Base64 codes:
    data = (plusMask & simde_mm_set1_epi8(62)) | (slashMask & simde_mm_set1_epi8(63)) | (numberMask & simde_mm_add_epi8(data, simde_mm_set1_epi8(4))) |
           (bigLetterMask & simde_mm_sub_epi8(data, simde_mm_set1_epi8(65))) | // ASCII 'A' is 65, Base64 'A' is 0
           (smallLetterMask & simde_mm_sub_epi8(data, simde_mm_set1_epi8(71)));

    // convert little endian to big endian:
    data = simde_mm_shuffle_epi8(data, shuffle_mask_2_);

    // the actual magic (conversion base64 to ASCII) happens here by shifting and masking:
    data = simde_mm_slli_epi32((data & mask1d_), 2) | simde_mm_slli_epi32((data & mask2d_), 4) | simde_mm_slli_epi32((data & mask3d_), 6) | simde_mm_slli_epi32((data & mask4d_), 8);

    // convert big endian to little endian
    data = simde_mm_shuffle_epi8(data, shuffle_mask_d_2_);
  }

  void Base64::stringSimdEncoder_(std::string& in, std::string& out)
  {
    out.resize((Size)(in.size() / 3) * 4 + 16); // resize output array, so the register encoder doesnt write memory to unallocated memory
    uint8_t padding = (3 - in.size() % 3) % 3;
    const int loop = in.size() / 12;

    in.resize(in.size() + 4, '\0');
    // otherwise there are cases where register encoder isnt allowed to access last bytes

    simde__m128i data {};
    // loop  through input as long as it's safe to access memory
    for (int i = 0; i < loop; i++)
    {
      // each time the last 4 out of 16 byte string data get lost through processing, therefore jumps of 12 bytes (/characters)
      data = simde_mm_lddqu_si128((simde__m128i*)&in[12 * i]);
      registerEncoder_(data);
      simde_mm_storeu_si128((simde__m128*)&out[i * 16], data);
    }

    size_t read = loop * 12;
    size_t written = loop * 16;

    // create buffer to translate last bytes without accessing memory that hasn't been allocated
    std::array<char, 16> buffer {};
    memcpy(&buffer[0], &in[read], in.size() - read - 4); // minus 4 because of 4 appended null bytes
    data = simde_mm_lddqu_si128((simde__m128i*)&buffer[0]);
    registerEncoder_(data);
    simde_mm_storeu_si128((simde__m128*)&out[written], data);

    in.resize(in.size() - 4); // remove null bytes

    // resizing out and add padding if necessary
    if (padding)
    {
      size_t newsize = ceil((double)in.size() / 3.) * 4;
      out.resize(newsize);
      for (size_t j = newsize - 1; j >= newsize - padding; j--)
      {
        out[j] = '=';
      }
    }
    else
    {
      out.resize((in.size() / 3) * 4);
    }
  }

  void Base64::stringSimdDecoder_(const std::string& in, std::string& out)
  {
    out.clear();
    const char* inPtr = &in[0];

    // padding count:
    uint8_t g = 0;
    if (in[in.size() - 1] == '=')
      g++;
    if (in[in.size() - 2] == '=')
      g++;

    unsigned outsize = (in.size() / 16) * 12 + 16;
    // not final size (final rezize later to cutoff unwanted characters)
    out.resize(outsize);
    char* outPtr = &out[0];
    int loop = in.size() / 16;

    for (int i = 0; i < loop; i++)
    {
      simde__m128i data = simde_mm_lddqu_si128((simde__m128i*)(inPtr + i * 16));
      registerDecoder_(data);
      simde_mm_storeu_si128((simde__m128*)(outPtr + i * 12), data);
    }

    size_t read = loop * 16;
    std::array<char, 16> rest;
    std::fill(rest.begin(), rest.end(), 'x');
    std::copy(in.begin() + read, in.end(), rest.begin());

    simde__m128i data = simde_mm_lddqu_si128((simde__m128i*)&rest[0]);
    registerDecoder_(data);
    size_t written = loop * 12;
    simde_mm_storeu_si128((simde__m128*)(outPtr + written), data);

    // cutting off decoding of appendix
    outsize = (in.size() / 4) * 3 - g;
    out.resize(outsize);
  }

  const char Base64::encoder_[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  const char Base64::decoder_[] = "|$$$}rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\\]^_`abcdefghijklmnopq";

  void Base64::encodeStrings(const std::vector<String>& in, String& out, bool zlib_compression, bool append_null_byte)
  {
    out.clear();
    if (in.empty())
    {
      return;
    }
    std::string str;
    for (Size i = 0; i < in.size(); ++i)
    {          
      str.append(in[i]);
      if (append_null_byte)
      {
        str.push_back('\0');
      }
    }

    if (zlib_compression)
    {
      String compressed;
      ZlibCompression::compressString(str, compressed);
      Base64::stringSimdEncoder_(compressed, out);
    }
    else
    {
      Base64::stringSimdEncoder_(str, out);
    }
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
    decodeSingleString(in, base64_uncompressed, zlib_compression);    //////////////////////////////////////////////the magic happenes here
    QList<QByteArray> null_strings = base64_uncompressed.split('\0');
    for (QList<QByteArray>::iterator it = null_strings.begin(); it < null_strings.end(); ++it)
    {
      if (!it->isEmpty())
      {
        out.emplace_back(QString(*it).toStdString());
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
    ////////////////////compare our decoding to QT decoding, and possibly decode first using simde, then copy into QByte Array
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
