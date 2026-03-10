// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/SYSTEM/SHA.h>

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

static std::string digest_to_hex(const SHA::digest_type& d)
{
  std::ostringstream oss;
  oss << std::hex << std::nouppercase << std::setfill('0');
  for (int i = 0; i < 5; ++i)
  {
    oss << std::setw(8) << static_cast<unsigned int>(d[i]);
  }
  return oss.str();
}

START_TEST(SHA, "$Id$")

START_SECTION((SHA()))
{
  SHA sha;
  // Empty message SHA-1
  // Expected: da39a3ee5e6b4b0d3255bfef95601890afd80709
  SHA::digest_type dg = {0,0,0,0,0};
  sha.get_digest(dg);

  std::string hex = digest_to_hex(dg);
  TEST_EQUAL(hex, "da39a3ee5e6b4b0d3255bfef95601890afd80709");
}
END_SECTION

START_SECTION((void process_bytes(void const* buffer, std::size_t byte_count)))
{
  // "abc"
  // Expected: a9993e364706816aba3e25717850c26c9cd0d89d
  const std::string s = "abc";

  // one-shot
  SHA sha1;
  sha1.process_bytes(s.data(), s.size());
  SHA::digest_type dg1 = {0,0,0,0,0};
  sha1.get_digest(dg1);
  std::string hex1 = digest_to_hex(dg1);
  TEST_EQUAL(hex1, "a9993e364706816aba3e25717850c26c9cd0d89d");

  // byte-by-byte
  SHA sha2;
  for (unsigned char c : s)
  {
    sha2.process_byte(c);
  }
  SHA::digest_type dg2 = {0,0,0,0,0};
  sha2.get_digest(dg2);
  std::string hex2 = digest_to_hex(dg2);
  TEST_EQUAL(hex2, "a9993e364706816aba3e25717850c26c9cd0d89d");
}
END_SECTION

START_SECTION((get_digest correctness on chunked input))
{
  // "The quick brown fox jumps over the lazy dog"
  // Expected: 2fd4e1c67a2d28fced849ee1bb76e7391b93eb12
  const std::string msg = "The quick brown fox jumps over the lazy dog";
  const std::vector<size_t> chunks = {1, 2, 3, 7, 13, 64, 8192};

  for (size_t chunk : chunks)
  {
    SHA sha;
    for (size_t pos = 0; pos < msg.size(); pos += chunk)
    {
      const size_t len = std::min(chunk, msg.size() - pos);
      sha.process_bytes(msg.data() + pos, len);
    }
    SHA::digest_type dg = {0,0,0,0,0};
    sha.get_digest(dg);
    std::string hex = digest_to_hex(dg);
    TEST_EQUAL(hex, "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
  }
}
END_SECTION

START_SECTION((reset()))
{
  SHA sha;
  const std::string s = "abc";
  sha.process_bytes(s.data(), s.size());
  SHA::digest_type dg1 = {0,0,0,0,0};
  sha.get_digest(dg1);
  TEST_EQUAL(digest_to_hex(dg1), "a9993e364706816aba3e25717850c26c9cd0d89d");

  // reuse after reset should produce same result
  sha.reset();
  sha.process_bytes(s.data(), s.size());
  SHA::digest_type dg2 = {0,0,0,0,0};
  sha.get_digest(dg2);
  TEST_EQUAL(digest_to_hex(dg2), "a9993e364706816aba3e25717850c26c9cd0d89d");
}
END_SECTION

END_TEST