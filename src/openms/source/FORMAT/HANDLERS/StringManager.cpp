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

#include <xercesc/util/XMLString.hpp>

#include <bit>
#include <cstdint>

using namespace xercesc;

namespace OpenMS::Internal
{
  // XMLCh is the Xerces UTF-16 code unit; it is guaranteed to match char16_t in
  // size and, on the supported toolchains, is char16_t. The StringManager public
  // API is expressed in char16_t so its header stays Xerces-free.
  static_assert(sizeof(::XMLCh) == sizeof(char16_t), "XMLCh is not sized correctly for UTF-16.");

  // Specializations for character types, released by XMLString::release
  template<> void shared_xerces_ptr<char>::doRelease_(char* item)
  {
    xercesc::XMLString::release(&item);
  }

  template<> void shared_xerces_ptr<char16_t>::doRelease_(char16_t* item)
  {
    xercesc::XMLString::release(reinterpret_cast<XMLCh**>(&item));
  }

  // Specializations for character types, which needs to be released by XMLString::release
  template <>
  void unique_xerces_ptr<char>::doRelease_(char*& item)
  {
    xercesc::XMLString::release(&item);
  }

  template <>
  void unique_xerces_ptr<char16_t>::doRelease_(char16_t*& item)
  {
    xercesc::XMLString::release(reinterpret_cast<XMLCh**>(&item));
  }

  unique_xerces_ptr<char16_t> StringManager::fromNative_(const char* str)
  {
    return unique_xerces_ptr<char16_t>(reinterpret_cast<char16_t*>(xercesc::XMLString::transcode(str)));
  }

  unique_xerces_ptr<char16_t> StringManager::fromNative_(const std::string& str)
  {
    return fromNative_(str.c_str());
  }

  std::string StringManager::toNative_(const char16_t* str)
  {
    std::string r;
    Size l = strLength(str);
    if (isASCII(str, l))
    {
      appendASCII(str, l, r);
    }
    else
    {
      r = (unique_xerces_ptr<char>(xercesc::XMLString::transcode(reinterpret_cast<const XMLCh*>(str))).get());
    }
    return r;
  }

  std::string StringManager::toNative_(const unique_xerces_ptr<char16_t>& str)
  {
    return toNative_(str.get());
  }

  StringManager::XercesString StringManager::convert(const char* str)
  {
    return fromNative_(str).get();
  }

  StringManager::XercesString StringManager::convert(const std::string& str)
  {
    return fromNative_(str.c_str()).get();
  }

  unique_xerces_ptr<char16_t> StringManager::convertPtr(const char* str)
  {
    return fromNative_(str);
  }

  unique_xerces_ptr<char16_t> StringManager::convertPtr(const std::string& str)
  {
    return fromNative_(str.c_str());
  }

  std::string StringManager::convert(const char16_t* str)
  {
    return toNative_(str);
  }

  Int StringManager::parseInt(const char16_t* str)
  {
    return xercesc::XMLString::parseInt(reinterpret_cast<const XMLCh*>(str));
  }

  // https://github.com/OpenMS/OpenMS/issues/8122
  #if defined(__GNUC__)
  __attribute__((no_sanitize("address")))
  #elif defined(_MSC_VER)
  __declspec(no_sanitize_address)
  #endif
  Size StringManager::strLength(const char16_t* input_ptr)
  {
    if (input_ptr == nullptr) {
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
        if (*pos_ptr == 0) {
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
    simde__m128i mask = simde_mm_set1_epi16(0xFF00);

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
      if (inputPtr[i] & 0xFF00)
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
