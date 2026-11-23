// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <memory>
#include <string>

namespace OpenMS::Internal
{
  /**
    @brief Strict UTF-8/UTF-16 conversion and ASCII helpers for XML handlers.

    Returned strings and arrays own their storage through the standard library.
    Invalid Unicode sequences throw Exception::ConversionError.
  */
  class OPENMS_DLLAPI StringManager
  {
protected:
    /// Compresses eight 8x16bit Chars in char16_t* to 8x8bit Chars by cutting upper byte
    static void compress64_ (const char16_t * input_it, char* output_it);

public:
    /// Constructor
    StringManager();

    /// Destructor
    ~StringManager();

    /// Calculates the length of a char16_t* string using SIMDe
    // https://github.com/OpenMS/OpenMS/issues/8122
    static Size strLength(const char16_t* input_ptr);

    /// Transcode the supplied C string to a UTF-16 string
    static std::u16string convert(const char * str);

    /**
      @brief Converts the UTF-8 sequence [@p str, @p str + @p length) to UTF-16 into @p result

      @p result is overwritten; its capacity is reused. Throws Exception::ConversionError for invalid UTF-8.
    */
    static void convert(const char * str, Size length, std::u16string & result);

    /// Transcode the supplied C++ string to a UTF-16 string
    static std::u16string convert(const std::string & str);

    /// Transcode the supplied C string to a UTF-16 string pointer
    static std::unique_ptr<char16_t[]> convertPtr(const char * str);

    /// Transcode the supplied C++ string to a UTF-16 string pointer
    static std::unique_ptr<char16_t[]> convertPtr(const std::string & str);

    /// Transcode the supplied char16_t* to a std::string
    static std::string convert(const char16_t * str);

    /// Parse an integer from a UTF-16 numeric string.
    static Int parseInt(const char16_t * str);

    /// Checks if supplied chars in char16_t* can be encoded with ASCII (i.e. the upper byte of each char is 0)
    static bool isASCII(const char16_t * chars, const Size length);

    /**
     * @brief Transcodes the supplied char16_t* and appends it to the OpenMS String
     *
     * @note Assumes that the char16_t* only contains ASCII characters
     *
    */
    static void appendASCII(const char16_t * str, const Size length, std::string & result);
  };

} // namespace OpenMS::Internal
