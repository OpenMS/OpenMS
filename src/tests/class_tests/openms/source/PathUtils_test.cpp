// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/SYSTEM/PathUtils.h>

using namespace OpenMS;

START_TEST(PathUtils, "$Id$")

START_SECTION(std::filesystem::path to_path(const std::string& s))
{
  // ASCII: identical on all platforms
  TEST_EQUAL(to_path("ASCII_only.txt").string(), "ASCII_only.txt")

  // UTF-8 encoded 'ä' (U+00E4 → 0xC3 0xA4): valid UTF-8, works on all platforms
  const std::string utf8_ae{'\xC3', '\xA4'};
  std::filesystem::path p_utf8 = to_path(utf8_ae + ".mzML");
  TEST_FALSE(p_utf8.empty())

#ifdef _WIN32
  // Windows-1252 encoded 'ä' (0xE4) — previously crashed with std::system_error
  // because the lone byte 0xE4 is not valid UTF-8. The fixed to_path() must fall
  // back to the ANSI code page and succeed.
  const std::string ansi_ae{'\xE4'};
  bool threw = false;
  std::filesystem::path p_ansi;
  try
  {
    p_ansi = to_path(ansi_ae + ".mzML");
  }
  catch (...)
  {
    threw = true;
  }
  TEST_FALSE(threw)
  TEST_FALSE(p_ansi.empty())
  // Verify the wide representation contains 'ä' (U+00E4)
  TEST_TRUE(p_ansi.wstring().find(L'\xE4') != std::wstring::npos)
#endif
}
END_SECTION

END_TEST
