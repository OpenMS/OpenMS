// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <string>

namespace OpenMS
{
  /// Convert a UTF-8 std::string to std::filesystem::path safely on all platforms.
  /// On Windows, std::filesystem::path(std::string) uses the current code page, not UTF-8,
  /// so we explicitly construct from u8string. If the bytes are not valid UTF-8 (e.g., a
  /// filename from argv in the ANSI code page), fall back to the native code page.
  inline std::filesystem::path to_path(const std::string& s)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    try
    {
      return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
    }
    catch (const std::system_error&)
    {
      // Not valid UTF-8 (e.g., ANSI-encoded filename from argv); let Windows use the current code page.
      return std::filesystem::path(s);
    }
#else
    return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
#endif
  }
} // namespace OpenMS
