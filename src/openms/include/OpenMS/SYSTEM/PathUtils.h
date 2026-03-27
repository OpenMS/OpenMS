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
  /// On Windows, std::filesystem::path(std::string) uses the current code page, not UTF-8.
  inline std::filesystem::path to_path(const std::string& s)
  {
    return std::filesystem::path(
      std::u8string(reinterpret_cast<const char8_t*>(s.data()), s.size()));
  }
} // namespace OpenMS
