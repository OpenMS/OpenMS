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

namespace
{
  // Read a filesystem::path back as UTF-8 bytes. u8string() is UTF-8 on every
  // platform (unlike string(), which is lossy on the Windows ANSI code page),
  // so it is the correct cross-platform way to assert a byte-exact round-trip.
  std::string asUtf8(const std::filesystem::path& p)
  {
    const std::u8string u8 = p.u8string();
    return std::string(reinterpret_cast<const char*>(u8.data()), u8.size());
  }
}

START_TEST(PathUtils, "$Id$")

START_SECTION(std::filesystem::path to_path(const std::string& s))
{
  // ASCII: identical on all platforms
  TEST_EQUAL(to_path("ASCII_only.txt").string(), "ASCII_only.txt")

  // UTF-8 encoded 'ä' (U+00E4 → 0xC3 0xA4): valid UTF-8, works on all platforms
  const std::string utf8_ae{'\xC3', '\xA4'};
  std::filesystem::path p_utf8 = to_path(utf8_ae + ".mzML");
  TEST_FALSE(p_utf8.empty())
  TEST_EQUAL(asUtf8(p_utf8), utf8_ae + ".mzML")

  // Round-trip invariant for valid UTF-8: to_path() -> u8string() must reproduce
  // the original bytes on every platform. All file I/O now flows through to_path()
  // after the Qt removal, so these guard the non-ASCII path handling end to end.

  // CJK (Japanese "日本語")
  const std::string cjk{'\xE6', '\x97', '\xA5', '\xE6', '\x9C', '\xAC', '\xE8', '\xAA', '\x9E'};
  const std::string cjk_file = cjk + ".mzML";
  std::filesystem::path p_cjk = to_path(cjk_file);
  TEST_FALSE(p_cjk.empty())
  TEST_EQUAL(asUtf8(p_cjk), cjk_file)

  // Embedded spaces and parentheses
  const std::string spaced = "my data file (run 1).mzML";
  std::filesystem::path p_spaced = to_path(spaced);
  TEST_FALSE(p_spaced.empty())
  TEST_EQUAL(asUtf8(p_spaced), spaced)

  // Long name (well beyond the historical Windows MAX_PATH of 260): constructing
  // the path object must not truncate or throw (the 260 limit applies to filesystem
  // operations, not to the path object itself).
  const std::string longname = std::string(300, 'a') + ".mzML";
  std::filesystem::path p_long = to_path(longname);
  TEST_FALSE(p_long.empty())
  TEST_EQUAL(asUtf8(p_long), longname)
  TEST_EQUAL(asUtf8(p_long).size(), longname.size())

  // Mixed CJK + Latin-1 diacritics with a space ("测试 äö")
  const std::string mixed{'\xE6', '\xB5', '\x8B', '\xE8', '\xAF', '\x95', ' ', '\xC3', '\xA4', '\xC3', '\xB6'};
  const std::string mixed_file = mixed + ".featureXML";
  std::filesystem::path p_mixed = to_path(mixed_file);
  TEST_FALSE(p_mixed.empty())
  TEST_EQUAL(asUtf8(p_mixed), mixed_file)

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
