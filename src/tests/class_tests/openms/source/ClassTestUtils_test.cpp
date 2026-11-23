// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/ClassTestUtils.h>

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

///////////////////////////

/*
  Differential test: the class-test framework against libOpenMS.

  OpenMSTestFramework depends on the C++ standard library only, so the numeric
  formatting and number parsing it needs are *copies* of libOpenMS code:

    * ClassTestUtils.h  detail::appendNumeric  <- StringUtils.cpp  appendNumeric
    * FuzzyStringComparator.cpp  tryParseNaN / fromCharsFloat / parseFloat /
      extractDouble                           <- StringUtils.cpp / StringUtils.h

  Those copies exist to produce *identical* results, and nothing in the build
  enforces that: the two sides compile independently, so a change to one is
  silently accepted. The consequence of drift is not a build error but a
  behavioural one -- TEST_EQUAL(std::string, <number>) would stop agreeing with
  StringUtils::toStr(), and FuzzyStringComparator would stop reading back the
  numbers the library writes into the TOPP reference files.

  This test is the guard. It is deliberately a *comparison* against libOpenMS
  rather than a table of expected strings: hard-coded expectations would have to
  be updated whenever the library's format legitimately changes, and would then
  no longer say anything about the two implementations agreeing.

  It lives in the class-test project (not in the framework) because only here
  are both implementations visible: the framework must not link libOpenMS.

  Verified by mutation: dropping one digit of precision in the framework's
  appendNumeric, making its extractDouble always fail, and making its parseFloat
  reject exponents are each caught here.
*/

namespace
{
  // Both sides of each comparison are already std::string, so TEST_EQUAL here is
  // a plain string comparison and does not itself depend on the code under test.
  //
  // STATUS() names the input, because the failure message only shows the two
  // formatted results -- which for a formatting bug can look confusingly similar.
#define TEST_SAME_FORMAT(value)                                                          \
  {                                                                                      \
    const std::string lib_ = OpenMS::StringUtils::toStr(value);                          \
    const std::string fw_ = OpenMS::Internal::ClassTest::detail::toString(value);        \
    if (fw_ != lib_) STATUS("input: " # value)                                           \
    TEST_EQUAL(fw_, lib_)                                                                \
  }

  // Values chosen to exercise every branch of appendNumeric: the to_chars call,
  // the trailing-zero trim, the "keep at least one fractional digit" rule, the
  // NaN/infinity short-circuits, and the extremes where precision runs out.
  const std::vector<double> doubles = {
    0.0, -0.0, 1.0, -1.0, 114.0, 0.5, -0.5,
    0.1, 1.0 / 3.0, 2.0 / 3.0, 1.5, 2.10, 100.5,
    1e-5, 1e5, 1e15, 1e16, 1e-15, 1e100, 1e-100, 1e300, 1e-300,
    3.14159265358979, 2.718281828459045,
    std::numeric_limits<double>::min(),
    std::numeric_limits<double>::max(),
    std::numeric_limits<double>::lowest(),
    std::numeric_limits<double>::epsilon(),
    std::numeric_limits<double>::denorm_min(),
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity()
  };

  const std::vector<float> floats = {
    0.0f, -0.0f, 1.0f, -1.0f, 114.0f, 0.5f, 0.1f, 1.5f, 2.10f,
    1e-5f, 1e5f, 1e10f, 1e-10f, 3.14159f,
    std::numeric_limits<float>::min(),
    std::numeric_limits<float>::max(),
    std::numeric_limits<float>::lowest(),
    std::numeric_limits<float>::epsilon(),
    std::numeric_limits<float>::denorm_min(),
    std::numeric_limits<float>::quiet_NaN(),
    std::numeric_limits<float>::infinity(),
    -std::numeric_limits<float>::infinity()
  };

  const std::vector<long double> long_doubles = {
    0.0L, -1.0L, 114.0L, 0.1L, 1.0L / 3.0L, 1e300L, 1e-300L,
    std::numeric_limits<long double>::min(),
    std::numeric_limits<long double>::max(),
    std::numeric_limits<long double>::epsilon(),
    std::numeric_limits<long double>::quiet_NaN(),
    std::numeric_limits<long double>::infinity()
  };
}

START_TEST(ClassTestUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


START_SECTION(([EXTRA] detail::toString(double) agrees with StringUtils::toStr(double)))
{
  for (double d : doubles) TEST_SAME_FORMAT(d)
}
END_SECTION

START_SECTION(([EXTRA] detail::toString(float) agrees with StringUtils::toStr(float)))
{
  for (float f : floats) TEST_SAME_FORMAT(f)
}
END_SECTION

START_SECTION(([EXTRA] detail::toString(long double) agrees with StringUtils::toStr(long double)))
{
  for (long double ld : long_doubles) TEST_SAME_FORMAT(ld)
}
END_SECTION

START_SECTION(([EXTRA] detail::toString(integral) agrees with StringUtils::toStr(integral)))
{
  // The framework formats integers with std::to_chars directly rather than
  // through appendNumeric; this checks that shortcut still matches the library.
  for (int i : {0, 1, -1, 114, -114, 2147483647, -2147483647}) TEST_SAME_FORMAT(i)
  for (unsigned int u : {0u, 1u, 114u, 4294967295u}) TEST_SAME_FORMAT(u)
  for (short s : {short(0), short(-1), short(32767), short(-32768)}) TEST_SAME_FORMAT(s)
  for (long l : {0L, -1L, 114L}) TEST_SAME_FORMAT(l)
  for (unsigned long ul : {0uL, 114uL}) TEST_SAME_FORMAT(ul)
  for (long long ll : {0LL, -1LL, 9223372036854775807LL}) TEST_SAME_FORMAT(ll)
  for (unsigned long long ull : {0uLL, 18446744073709551615uLL}) TEST_SAME_FORMAT(ull)

  TEST_SAME_FORMAT('x')
  TEST_SAME_FORMAT(std::string("already a string"))
  TEST_SAME_FORMAT("a char pointer")
}
END_SECTION

START_SECTION(([EXTRA] ClassTest::writtenDigits agrees with OpenMS::writtenDigits))
{
  // The precision fed to appendNumeric. Only the floating point types have to
  // agree -- they decide the format of every number the library writes.
  TEST_EQUAL(OpenMS::Internal::ClassTest::writtenDigits<float>(), OpenMS::writtenDigits<float>())
  TEST_EQUAL(OpenMS::Internal::ClassTest::writtenDigits<double>(), OpenMS::writtenDigits<double>())
  TEST_EQUAL(OpenMS::Internal::ClassTest::writtenDigits<long double>(), OpenMS::writtenDigits<long double>())

  // Documented, deliberate divergences (report precision only, never the
  // pass/fail decision of TEST_REAL_SIMILAR). Pinned here so that they stay
  // deliberate: if someone "fixes" one of them, this is where it shows up.
  TEST_EQUAL(OpenMS::Internal::ClassTest::writtenDigits<long>(), std::numeric_limits<long>::digits10)
  TEST_EQUAL(OpenMS::Internal::ClassTest::writtenDigits<unsigned long>(), std::numeric_limits<unsigned long>::digits10)
  // ... whereas the library clamps both to int's digits10.
  TEST_EQUAL(OpenMS::writtenDigits<long>(), std::numeric_limits<int>::digits10)
}
END_SECTION

START_SECTION(([EXTRA] FuzzyStringComparator parses the numbers StringUtils writes))
{
  // The discriminating property: two *textually different* spellings of nearly
  // the same value compare equal only if the comparator really parsed them as
  // numbers. Were its copy of the parser to stop recognising the library's
  // output, the comparison would fall back to text and report a difference.
  FuzzyStringComparator fsc;
  std::ostringstream sink; // keep the failure report out of the test log
  fsc.setLogDestination(sink);
  fsc.setAcceptableRelative(1.001); // 0.1 %
  fsc.setAcceptableAbsolute(0.0);

  for (double d : doubles)
  {
    if (!std::isfinite(d) || d == 0.0) continue; // perturbing these is meaningless
    const double nudged_value = d * 1.0001; // 0.01 %, inside the tolerance
    if (!std::isfinite(nudged_value)) continue; // nudging the extremes overflows to infinity

    const std::string written = StringUtils::toStr(d);
    const std::string nudged = StringUtils::toStr(nudged_value);

    if (nudged == written) continue; // precision exhausted, nothing to tell apart
    if (!fsc.compareStrings(written, nudged)) STATUS("input: " + written + " vs " + nudged)
    TEST_EQUAL(fsc.compareStrings(written, nudged), true)
  }

  // ... and a difference outside the tolerance must still be caught, so the
  // check above cannot pass by the comparator simply accepting everything.
  TEST_EQUAL(fsc.compareStrings(StringUtils::toStr(1.0), StringUtils::toStr(2.0)), false)

  // The spelling the library emits ("NaN") must read back as the same value as
  // the lower-case spelling. Note this is a weaker check than the loop above:
  // the copied tryParseNaN is not the only path, because std::from_chars accepts
  // "nan" too -- so it pins the accepted spellings, not that helper.
  TEST_EQUAL(fsc.compareStrings("NaN", "nan"), true)
  TEST_EQUAL(StringUtils::toStr(std::numeric_limits<double>::quiet_NaN()), "NaN")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
