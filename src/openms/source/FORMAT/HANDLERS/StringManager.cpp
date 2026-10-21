// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/StringManager.h>

#include <cmath> // MSVC: SIMDe's simde-math.h needs std::isnormal/isnan/... in scope
#include <OpenMS/SYSTEM/SIMDe.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>

namespace OpenMS::Internal
{
  namespace
  {
    [[noreturn]] void invalidUnicode()
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid UTF-8 or UTF-16 sequence");
    }
  }

  std::u16string StringManager::convert(const char* str)
  {
    std::u16string result;
    if (!str) return result;
    const auto* input = reinterpret_cast<const unsigned char*>(str);
    const auto* end = input + std::strlen(str);
    result.reserve(end - input);
    while (input != end)
    {
      const unsigned lead = *input++;
      unsigned code = lead;
      unsigned continuation = 0;
      unsigned minimum = 0;
      if (lead < 0x80) {}
      else if (lead >= 0xc2 && lead <= 0xdf) { code &= 0x1f; continuation = 1; minimum = 0x80; }
      else if (lead >= 0xe0 && lead <= 0xef) { code &= 0x0f; continuation = 2; minimum = 0x800; }
      else if (lead >= 0xf0 && lead <= 0xf4) { code &= 0x07; continuation = 3; minimum = 0x10000; }
      else invalidUnicode();
      if (static_cast<Size>(end - input) < continuation) invalidUnicode();
      for (unsigned j = 0; j < continuation; ++j)
      {
        if ((*input & 0xc0) != 0x80) invalidUnicode();
        code = (code << 6) | (*input++ & 0x3f);
      }
      if (code < minimum || code > 0x10ffff || (code >= 0xd800 && code <= 0xdfff)) invalidUnicode();
      if (code <= 0xffff) result.push_back(static_cast<char16_t>(code));
      else
      {
        code -= 0x10000;
        result.push_back(static_cast<char16_t>(0xd800 + (code >> 10)));
        result.push_back(static_cast<char16_t>(0xdc00 + (code & 0x3ff)));
      }
    }
    return result;
  }

  std::u16string StringManager::convert(const std::string& str)
  {
    return convert(str.c_str());
  }

  std::unique_ptr<char16_t[]> StringManager::convertPtr(const char* str)
  {
    const auto converted = convert(str);
    auto result = std::make_unique<char16_t[]>(converted.size() + 1);
    std::copy(converted.begin(), converted.end(), result.get());
    result[converted.size()] = 0;
    return result;
  }

  std::unique_ptr<char16_t[]> StringManager::convertPtr(const std::string& str)
  {
    return convertPtr(str.c_str());
  }

  std::string StringManager::convert(const char16_t* str)
  {
    std::string result;
    if (!str) return result;
    const Size length = std::char_traits<char16_t>::length(str);
    if (isASCII(str, length))
    {
      appendASCII(str, length, result);
      return result;
    }
    for (Size i = 0; i < length; ++i)
    {
      unsigned code = str[i];
      if (code >= 0xd800 && code <= 0xdbff)
      {
        if (i + 1 == length || str[i + 1] < 0xdc00 || str[i + 1] > 0xdfff) invalidUnicode();
        code = 0x10000 + ((code - 0xd800) << 10) + (str[++i] - 0xdc00);
      }
      else if (code >= 0xdc00 && code <= 0xdfff) invalidUnicode();
      if (code < 0x80) result.push_back(static_cast<char>(code));
      else if (code < 0x800)
      {
        result.push_back(static_cast<char>(0xc0 | (code >> 6)));
        result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
      }
      else if (code < 0x10000)
      {
        result.push_back(static_cast<char>(0xe0 | (code >> 12)));
        result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
      }
      else
      {
        result.push_back(static_cast<char>(0xf0 | (code >> 18)));
        result.push_back(static_cast<char>(0x80 | ((code >> 12) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3f)));
        result.push_back(static_cast<char>(0x80 | (code & 0x3f)));
      }
    }
    return result;
  }

  Int StringManager::parseInt(const char16_t* str)
  {
    return StringUtils::toInt32(convert(str));
  }

  // https://github.com/OpenMS/OpenMS/issues/8122
  #if defined(__GNUC__)
  __attribute__((no_sanitize("address")))
  #elif defined(_MSC_VER)
  __declspec(no_sanitize_address)
  #endif
  Size StringManager::strLength(const char16_t* input_ptr)
  {
    if (input_ptr == nullptr)
    {
      return 0;
    }

    Size processed_chars = 0;
    const char16_t* pos_ptr = input_ptr;

    // find number of characters until we reach a 16-byte aligned address
    uintptr_t ptr_value = reinterpret_cast<uintptr_t>(pos_ptr);
    size_t misalignment = ptr_value & 15; // modulo 16 to find misalignment
    int chars_to_align = misalignment ? (16 - misalignment) / sizeof(char16_t) : 0;

    // process one by one until we reach a 16-byte aligned address (where a SIMD load can be used)
    for (int i = 0; i < chars_to_align; ++i)
    {
      if (*pos_ptr == 0)
      {
        return i;
      }
      ++pos_ptr;
    }
    processed_chars = chars_to_align;

    // SIMD loop: find the first zero character in blocks of 8 UTF-16 characters (16 bytes each)
    const simde__m128i zero = simde_mm_setzero_si128();
    while (true)
    {
      simde__m128i bits = simde_mm_load_si128(reinterpret_cast<const simde__m128i*>(pos_ptr));
      simde__m128i cmp_zero = simde_mm_cmpeq_epi16(bits, zero); // sets bits to 0xFFFF (2 bytes) for each character that is zero
      uint16_t zero_mask = simde_mm_movemask_epi8(cmp_zero);    // extracts MSB from each byte

      if (zero_mask != 0)
      { // Found a zero character
        auto byte_pos_zero = std::countr_zero(zero_mask); // count trailing zeros to find the first zero character
        auto char_pos_zero = byte_pos_zero / 2;           // each UTF-16 character is 2 bytes, so divide by 2 to get character position
        return processed_chars + char_pos_zero;
      }

      // 8 chars (each 2 bytes) had no zero
      pos_ptr += 8;
      processed_chars += 8;
    }

    // never reached
    return processed_chars;
  }

  void StringManager::compress64_(const char16_t* inputIt, char* outputIt)
  {
    simde__m128i bits = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(inputIt));

    // Select every second byte (little-endian lower byte of each UTF-16 character)
    const simde__m128i shuffleMask = simde_mm_setr_epi8(
      0, 2, 4, 6, 8, 10, 12, 14,
      -1, -1, -1, -1, -1, -1, -1, -1
    );

    simde__m128i compressed = simde_mm_shuffle_epi8(bits, shuffleMask);

    // Store the lower 64 bits (8 ASCII characters)
    simde_mm_storel_epi64(reinterpret_cast<simde__m128i*>(outputIt), compressed);
  }

  bool StringManager::isASCII(const char16_t* chars, const Size length)
  {
    if (length == 0)
    {
      return true;
    }

    Size fullBlocks = length / 8;
    Size remainder = length % 8;

    const char16_t* inputPtr = chars;
    simde__m128i mask = simde_mm_set1_epi16(0xFF80);

    // Process blocks of 8 UTF-16 characters using SIMD
    for (Size i = 0; i < fullBlocks; ++i)
    {
      simde__m128i bits = simde_mm_loadu_si128(reinterpret_cast<const simde__m128i*>(inputPtr));
      simde__m128i zero = simde_mm_setzero_si128();
      simde__m128i andOp = simde_mm_and_si128(bits, mask);
      simde__m128i cmp = simde_mm_cmpeq_epi16(andOp, zero);

      if (simde_mm_movemask_epi8(cmp) != 0xFFFF)
      {
        return false;
      }

      inputPtr += 8;
    }

    // Check remaining characters individually
    for (Size i = 0; i < remainder; ++i)
    {
      if (inputPtr[i] & 0xFF80)
      {
        return false;
      }
    }

    return true;
  }

  void StringManager::appendASCII(const char16_t* chars, const Size length, std::string& result)
  {
    // char16_t are characters in UTF16.
    // We know that the Base64 string here can only contain plain ASCII
    // and all bytes except the least significant one will be zero. Thus
    // we can convert to char directly (only keeping the least
    // significant byte).

    Size fullBlocks = length / 8;
    Size remainder = length % 8;

    const char16_t* inputPtr = chars;

    Size currentSize = result.size();
    result.resize(currentSize + length);
    char* outputPtr = &result[currentSize];

    // Copy blocks of 8 characters at a time
    for (Size i = 0; i < fullBlocks; ++i)
    {
      compress64_(inputPtr, outputPtr);
      inputPtr += 8;
      outputPtr += 8;
    }

    // Copy any remaining characters individually
    for (Size i = 0; i < remainder; ++i)
    {
      outputPtr[i] = static_cast<char>(inputPtr[i] & 0xFF);
    }
  }

  StringManager::StringManager() = default;

  StringManager::~StringManager() = default;

} // namespace OpenMS::Internal
