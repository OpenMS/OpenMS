// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/SYSTEM/SIMDe.h>

#include <QtCore/QList>
#include <QtCore/QString>
#include <OpenMS/FORMAT/Vec.h>

using namespace std;

namespace OpenMS
{
  const Vec128 mask1_ = simde_mm_set1_epi32(0x3F000000);
  const Vec128 mask2_ = simde_mm_set1_epi32(0x003F0000);
  const Vec128 mask3_ = simde_mm_set1_epi32(0x00003F00);
  const Vec128 mask4_ = simde_mm_set1_epi32(0x0000003F);

  const Vec128 mask1d_ = simde_mm_set1_epi32(0xFF000000);
  const Vec128 mask2d_ = simde_mm_set1_epi32(0x00FF0000);
  const Vec128 mask3d_ = simde_mm_set1_epi32(0x0000FF00);
  const Vec128 mask4d_ = simde_mm_set1_epi32(0x000000FF);

  const Vec128 difference_A_ = simde_mm_set1_epi8('A');
  const Vec128 difference_a_ = simde_mm_set1_epi8('a' - 26);
  const Vec128 difference_0_ = simde_mm_set1_epi8('0' - 52);
  const Vec128 difference_plus_ = simde_mm_set1_epi8('+');
  const Vec128 difference_slash_ = simde_mm_set1_epi8('/');

  const Vec128 shuffle_mask_1_ = simde_mm_setr_epi8(2, 2, 1, 0, 5, 5, 4, 3, 8, 8, 7, 6, 11, 11, 10, 9);
  const Vec128 shuffle_mask_2_ = simde_mm_setr_epi8(3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12);
  const Vec128 shuffle_mask_big_endian_ = simde_mm_setr_epi8(0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11);
  const Vec128 shuffle_mask_d_2_ = simde_mm_setr_epi8(3, 2, 1, 7, 6, 5, 11, 10, 9, 15, 14, 13, 0, 4, 8, 12);

  void registerEncoder_(Vec128& data_raw)
  {
    Vec128 data(data_raw);

    data = (!OPENMS_IS_BIG_ENDIAN ? data.shuffle_epi8(shuffle_mask_1_) : data.shuffle_epi8(shuffle_mask_big_endian_));

    data = (data.srli_epi32(2) & mask1_) | (data.srli_epi32(4) & mask2_) |
           (data.srli_epi32(6) & mask3_) | (data & mask4_);

    if (!OPENMS_IS_BIG_ENDIAN)
      data = data.shuffle_epi8(shuffle_mask_2_);

    Vec128 capital_mask = Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8(26)));
    Vec128 all_mask = capital_mask;

    Vec128 lower_case_mask = Vec128::andnot(all_mask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8(52))));
    all_mask |= lower_case_mask;

    Vec128 number_mask = Vec128::andnot(all_mask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8(62))));
    all_mask |= number_mask;

    Vec128 plus_mask = Vec128::andnot(all_mask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8(63))));
    all_mask |= plus_mask;

    Vec128 slash_mask = all_mask;

    data = (capital_mask & (data + difference_A_)) |
           (lower_case_mask & (data + difference_a_)) |
           (number_mask & (data + difference_0_)) |
           (plus_mask & difference_plus_) |
           (Vec128::andnot(slash_mask, difference_slash_));

    data_raw = data;
  }

  void registerDecoder_(Vec128& data_raw)
  {
    Vec128 data(data_raw);

    Vec128 plusMask = Vec128::cmpeq_epi8(data, difference_plus_);
    Vec128 allMask = plusMask;
    Vec128 slashMask = Vec128::cmpeq_epi8(data, difference_slash_);
    allMask |= slashMask;
    Vec128 numberMask = Vec128::andnot(allMask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8('9' + 1))));
    allMask |= numberMask;
    Vec128 bigLetterMask = Vec128::andnot(allMask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8('Z' + 1))));
    allMask |= bigLetterMask;
    Vec128 smallLetterMask = Vec128::andnot(allMask, Vec128::cmplt_epi8(data, Vec128(simde_mm_set1_epi8('z' + 1))));

    data = (plusMask & Vec128(simde_mm_set1_epi8(62))) |
           (slashMask & Vec128(simde_mm_set1_epi8(63))) |
           (numberMask & (data + Vec128(simde_mm_set1_epi8(4)))) |
           (bigLetterMask & (data - Vec128(simde_mm_set1_epi8(65)))) |
           (smallLetterMask & (data - Vec128(simde_mm_set1_epi8(71))));

    data = data.shuffle_epi8(shuffle_mask_2_);

    data = ((data & mask1d_).slli_epi32(2)) |
           ((data & mask2d_).slli_epi32(4)) |
           ((data & mask3d_).slli_epi32(6)) |
           ((data & mask4d_).slli_epi32(8));

    data = data.shuffle_epi8(shuffle_mask_d_2_);

    data_raw = data;
  }
   void Base64::stringSimdEncoder_(std::string& in, std::string& out)
  {
    out.resize((Size)(in.size() / 3) * 4 + 16); // resize output array, so the register encoder doesnt write memory to unallocated memory
    uint8_t padding = (3 - in.size() % 3) % 3;
    const int loop = in.size() / 12;

    in.resize(in.size() + 4, '\0');
    // otherwise there are cases where register encoder isnt allowed to access last bytes

    Vec128 data {};
    // loop  through input as long as it's safe to access memory
    for (int i = 0; i < loop; i++)
    {
      // each time the last 4 out of 16 byte string data get lost through processing, therefore jumps of 12 bytes (/characters)
      const simde__m128i* ptr = reinterpret_cast<const simde__m128i*>(&in[12*i]);
      data = Vec128(simde_mm_lddqu_si128(ptr));
      registerEncoder_(data);
      simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(&out[i * 16]), data);
    }

    size_t read = loop * 12;
    size_t written = loop * 16;

    // create buffer to translate last bytes without accessing memory that hasn't been allocated
    std::array<char, 16> buffer {};
    memcpy(&buffer[0], &in[read], in.size() - read - 4); // minus 4 because of 4 appended null bytes

    const simde__m128i* buffer_ptr = reinterpret_cast<const simde__m128i*>(&buffer[0]);
    data = Vec128(simde_mm_lddqu_si128(buffer_ptr));
    registerEncoder_(data);
    simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(&out[written]), data);

    in.resize(in.size() - 4); // remove null bytes

    // resizing out and add padding if necessary
    if (padding)
    {
      // size_t newsize = ceil((double)in.size() / 3.) * 4;
      size_t newsize = static_cast<size_t>(ceil(static_cast<double>(in.size()) / 3.)) * 4;
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
      const simde__m128i raw = simde_mm_lddqu_si128(reinterpret_cast<const simde__m128i*>(inPtr + i * 16));
      Vec128 data(raw);
      registerDecoder_(data);
      simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(outPtr + i * 12), data.raw());
    }
    size_t read = loop * 16;
    std::array<char, 16> rest;
    std::fill(rest.begin(), rest.end(), 'x');
    std::copy(in.begin() + read, in.end(), rest.begin());

    const simde__m128i raw = simde_mm_lddqu_si128(reinterpret_cast<const simde__m128i*>(&rest[0]));
    Vec128 data(raw);
    registerDecoder_(data);
    size_t written = loop * 12;
    simde_mm_storeu_si128(reinterpret_cast<simde__m128i*>(outPtr + written), data.raw());
    // cutting off decoding of appendix
    outsize = (in.size() / 4) * 3 - g;
    out.resize(outsize); 
    }
<<<<<<< Updated upstream
=======
  

>>>>>>> Stashed changes
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
  // stringSimdEncoder_ and stringSimdDecoder_ don't need modification unless you want to use Vec128::load/store
  // For now, the core logic is adapted above. Let me know if you want to refactor those too.
} // namespace OpenMS
