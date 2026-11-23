// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/PortableMath.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
///////////////////////////

using namespace OpenMS;

namespace
{
  std::uint64_t bits(double x)
  {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof(b));
    return b;
  }
}

START_TEST(PortableMath, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// The point of these functions is that they give the same bits on every platform, so the expected values are bit
// patterns (fdlibm's results), not tolerances. A failure on one platform means a compiler contracted or reordered
// the arithmetic.

START_SECTION((double exp(double x)))
{
  const std::pair<double, std::uint64_t> cases[] = {
    {-745, 0x0000000000000001ULL},
    {-20.5, 0x3e157a3afeed00abULL},
    {-1, 0x3fd78b56362cef38ULL},
    {-1e-10, 0x3feffffffff24190ULL},
    {0, 0x3ff0000000000000ULL},
    {1e-10, 0x3ff000000006df38ULL},
    {0.5, 0x3ffa61298e1e069cULL},
    {1, 0x4005bf0a8b14576aULL},
    {2.3025850929940459, 0x4024000000000001ULL},
    {10, 0x40d5829dcf950560ULL},
    {88.700000000000003, 0x47ef4705bbffae5cULL},
    {709.70000000000005, 0x7fed75ae7a50ee14ULL},
  };
  for (const auto& [x, expected] : cases)
  {
    TEST_EQUAL(bits(PortableMath::exp(x)), expected)
    TEST_REAL_SIMILAR(PortableMath::exp(x), std::exp(x))
  }
  TEST_EQUAL(PortableMath::exp(1000.0), std::numeric_limits<double>::infinity())
  TEST_EQUAL(PortableMath::exp(-1000.0), 0.0)
  TEST_EQUAL(std::isnan(PortableMath::exp(std::numeric_limits<double>::quiet_NaN())), true)
}
END_SECTION

START_SECTION((double log(double x)))
{
  const std::pair<double, std::uint64_t> cases[] = {
    {1e-300, 0xc085963447f87fb5ULL},
    {1e-10, 0xc037069e2aa2aa5bULL},
    {0.10000000000000001, 0xc0026bb1bbb55515ULL},
    {0.5, 0xbfe62e42fefa39efULL},
    {1, 0x0000000000000000ULL},
    {1.5, 0x3fd9f323ecbf984cULL},
    {2.7182818284590451, 0x3ff0000000000000ULL},
    {10, 0x40026bb1bbb55516ULL},
    {1013.4567, 0x401baf3aab3f9e74ULL},
    {1.0000000000000001e+300, 0x4085963447f87fb5ULL},
  };
  for (const auto& [x, expected] : cases)
  {
    TEST_EQUAL(bits(PortableMath::log(x)), expected)
    TEST_REAL_SIMILAR(PortableMath::log(x), std::log(x))
  }
  TEST_EQUAL(PortableMath::log(0.0), -std::numeric_limits<double>::infinity())
  TEST_EQUAL(std::isnan(PortableMath::log(-1.0)), true)
  TEST_EQUAL(PortableMath::log(std::numeric_limits<double>::infinity()), std::numeric_limits<double>::infinity())
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
