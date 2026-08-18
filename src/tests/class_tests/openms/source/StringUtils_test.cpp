// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <cmath>
#include <limits>

using namespace OpenMS;
using namespace std;
using namespace StringUtils;

START_TEST(StringUtils, "$Id$")

string whitespaces = "\t\r\n "; ///< all whitespaces we need to test

START_SECTION(inline const char* skipWhitespace(const char* p, const char* p_end))
{
  // postfix with 16x, to enable SIMD on the prefix
  #define x16 "xxxxxxxxxxxxxxxx"
  #define s16 "                "
  for (const char whitespace : whitespaces)
  {
    std::string at1 = "0 2  3456789101112" x16;
    StringUtils::substitute(at1, ' ', whitespace);
    TEST_EQUAL(skipWhitespace(at1), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 1)), 1);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 2)), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 3)), 2);
    std::string at2 = s16 s16 "1" x16;
    StringUtils::substitute(at2, ' ', whitespace);
    TEST_EQUAL(skipWhitespace(std::string_view(at2.data())), 32);
    TEST_EQUAL(skipWhitespace(std::string_view(at2.data() + 2)), 30);
    std::string at1_noSSE = "0 2  34";
    StringUtils::substitute(at1_noSSE, ' ', whitespace);
    TEST_EQUAL(skipWhitespace(std::string_view(at1_noSSE.data())), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1_noSSE.data() + 1)), 1);
    TEST_EQUAL(skipWhitespace(std::string_view(at1_noSSE.data() + 2)), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1_noSSE.data() + 3)), 2);
  }
}
END_SECTION

START_SECTION(inline const char* skipNonWhitespace(const char* p, const char* p_end))
{
  // postfix with 16x, to enable SIMD on the prefix
  #define x16 "xxxxxxxxxxxxxxxx"
  #define s16 "                "
  for (const char whitespace : whitespaces)
  {
    std::string at1 = "0 2  3456789101112" x16;
    StringUtils::substitute(at1, ' ', whitespace);
    TEST_EQUAL(skipNonWhitespace(at1), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 1)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 2)), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 3)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 5)), 13 + 16);
    std::string at2 = x16 x16 " " x16;
    StringUtils::substitute(at2, ' ', whitespace);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data())), 32);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data() + 31)), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data() + 33)), 16);
    std::string at1_noSSE = "0 2  34";
    StringUtils::substitute(at1_noSSE, ' ', whitespace);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1_noSSE.data())), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1_noSSE.data() + 1)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1_noSSE.data() + 2)), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1_noSSE.data() + 3)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1_noSSE.data() + 5)), 2);
  }
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StringUtilsHelper* ptr = nullptr;
StringUtilsHelper* null_ptr = nullptr;



START_SECTION(StringUtilsHelper())
{
	ptr = new StringUtilsHelper();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~StringUtilsHelper())
{
	delete ptr;
}
END_SECTION

START_SECTION((static std::string numberLength(double d, UInt n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string number(double d, UInt n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& fillLeft(String &this_s, char c, UInt size)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& fillRight(String &this_s, char c, UInt size)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasPrefix(const std::string &this_s, const std::string &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasSuffix(const std::string &this_s, const std::string &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasSubstring(const std::string &this_s, const std::string &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool has(const std::string &this_s, Byte byte)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string prefix(const std::string &this_s, size_t length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string suffix(const std::string &this_s, size_t length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string prefix(const std::string &this_s, Int length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string suffix(const std::string &this_s, Int length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string prefix(const std::string &this_s, char delim)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string suffix(const std::string &this_s, char delim)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string substr(const std::string &this_s, size_t pos, size_t n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string chop(const std::string &this_s, Size n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& trim(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& quote(String &this_s, char q, OpenMS::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& unquote(String &this_s, char q, OpenMS::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& simplify(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string random(UInt length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& reverse(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split(const std::string &this_s, const char splitter, std::vector<std::string> &substrings, bool quote_protect)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split(const std::string &this_s, const std::string &splitter, std::vector<std::string> &substrings)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split_quoted(const std::string &this_s, const std::string &splitter, std::vector<std::string> &substrings, char q, OpenMS::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static QString toQString(const std::string &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION


START_SECTION((static Int32 toInt32(const std::string &this_s)))
{
  // easy case
  TEST_EQUAL(StringUtils::toInt32("2147483647"), 2147483647)
  // with spaces (allowed)
  TEST_EQUAL(StringUtils::toInt32("  2147483647"), 2147483647)
  TEST_EQUAL(StringUtils::toInt32("2147483647 "), 2147483647)
  TEST_EQUAL(StringUtils::toInt32("   2147483647  "), 2147483647)
  //
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32("2147483648")) // +1 too large
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32("-2147483649")) // -1 too small

  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32("1234  moreText")) // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32(" 1234 911.0"))    // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32(" abc "))
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32(" 123.45 "))
}
END_SECTION



START_SECTION((static Int64 toInt64(const std::string &this_s)))
{
  // easy case
  TEST_EQUAL(StringUtils::toInt64("9223372036854775807"), 9223372036854775807)
  // with spaces (allowed)
  TEST_EQUAL(StringUtils::toInt64("  9223372036854775807"), 9223372036854775807)
  TEST_EQUAL(StringUtils::toInt64("9223372036854775807 "), 9223372036854775807)
  TEST_EQUAL(StringUtils::toInt64("   9223372036854775807  "), 9223372036854775807)
  //
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64("9223372036854775808")) // +1 too large
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64("-9223372036854775809")) // -1 too small

  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64("1234  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64(" 1234 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64(" abc "))
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt64(" 123.45 "))
}
END_SECTION

START_SECTION((static float toFloat(const std::string &this_s)))
{
  // easy case
  TEST_REAL_SIMILAR(StringUtils::toFloat("1234.45"), 1234.45)
  // with spaces (allowed)
  TEST_REAL_SIMILAR(StringUtils::toFloat("  1234.45"), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toFloat("1234.45 "), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toFloat("   1234.45  "), 1234.45)
  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat("1234.45  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat(" 1234.45 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat(" abc "))
  // subnormal value: underflows below the smallest *normal* float but is still representable --> must parse (no out-of-range)
  TEST_TRUE(StringUtils::toFloat("1e-40") > 0.0f)
  // overflow above FLT_MAX --> error
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat("1e40"))
}
END_SECTION

START_SECTION((static double toDouble(const std::string &this_s)))
{
  // easy case
  TEST_REAL_SIMILAR(StringUtils::toDouble("1234.45"), 1234.45)
  // with spaces (allowed)
  TEST_REAL_SIMILAR(StringUtils::toDouble("  1234.45"), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toDouble("1234.45 "), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toDouble("   1234.45  "), 1234.45)
  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble("1234.45  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble(" 1234.45 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble(" abc "))
  // subnormal value: underflows below the smallest *normal* double but is still representable --> must parse (no out-of-range)
  TEST_EQUAL(StringUtils::toDouble("2.17388884170148e-321"), 2.17388884170148e-321)
  TEST_TRUE(StringUtils::toDouble("4.9e-324") > 0.0) // ~smallest positive subnormal double
  // overflow above DBL_MAX --> error
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble("1e400"))
}
END_SECTION


START_SECTION((template <typename IteratorT> static bool extractDouble(IteratorT& begin, const IteratorT& end, double& target)))
{
  double d;
  {
    std::string ss("12345.45  ");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 12345.45)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 8); // was the iterator advanced?
  }

  {
    std::string ss("+1234.45!");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 1234.45)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 8); // was the iterator advanced?
  }
  {
    d = 0;
    std::string ss("  -123.45");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), false);
    TEST_REAL_SIMILAR(d, 0)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 0); // was the iterator advanced?
  }
  {
    std::string ss("15.0e6");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 15.0e6)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 6); // was the iterator advanced?
  }
  {
    // try two doubles in a single stream (should stop after the first)
    std::string ss("-5.0	9.1");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, -5.0)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 4); // was the iterator advanced?
    auto it2 = ss.begin() + 5;
    TEST_EQUAL(StringUtils::extractDouble(it2, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 9.1)
    TEST_EQUAL((int)std::distance(ss.begin(), it2), 8); // was the iterator advanced?
  }
  {
    // explicitly test X.FeY vs XeY since some compilers implementation of the native operator>> stop reading at 'e' if no '.F' was seen
    std::string ss("15.0e6 x");   
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 15.0e6)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 6); // was the iterator advanced?
  }
  {
    std::string ss("16e6!");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_REAL_SIMILAR(d, 16e+06)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 4); // was the iterator advanced?
  }
  {
    std::string ss("!noNumber");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), false);
    TEST_EQUAL((int)std::distance(ss.begin(), it), 0); // was the iterator advanced?
  }
  {
    // subnormal value (underflow): must be parsed, not rejected as out-of-range
    std::string ss("2.17388884170148e-321 x");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractDouble(it, ss.end(), d), true);
    TEST_EQUAL(d, 2.17388884170148e-321);
    TEST_EQUAL((int)std::distance(ss.begin(), it), 21); // was the iterator advanced past the full number?
  }
}
END_SECTION

START_SECTION((template <typename IteratorT> static bool extractInt(IteratorT& begin, const IteratorT& end, int& target)))
{
  int i;
  {
    std::string ss("12345  ");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), true);
    TEST_EQUAL(i, 12345)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 5); // was the iterator advanced?
  }

  {
    std::string ss("+1234!");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), true);
    TEST_EQUAL(i, 1234)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 5); // was the iterator advanced?
  }
  {
    i = 0;
    std::string ss("  -123");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), false);
    TEST_EQUAL(i, 0)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 0); // was the iterator advanced?
  }
  {
    std::string ss("y15++");
    auto it = ss.begin() + 1;
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), true);
    TEST_EQUAL(i, 15)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 3); // was the iterator advanced?
  }
  {
    // try two ints in a single stream (should stop after the first)
    std::string ss("-5	9");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), true);
    TEST_EQUAL(i, -5)
    TEST_EQUAL((int)std::distance(ss.begin(), it), 2); // was the iterator advanced?
    auto it2 = ss.begin() + 3;
    TEST_EQUAL(StringUtils::extractInt(it2, ss.end(), i), true);
    TEST_EQUAL(i, 9)
    TEST_EQUAL((int)std::distance(ss.begin(), it2), 4); // was the iterator advanced?
  }
  {
    std::string ss("!noNumber");
    auto it = ss.begin();
    TEST_EQUAL(StringUtils::extractInt(it, ss.end(), i), false);
    TEST_EQUAL((int)std::distance(ss.begin(), it), 0); // was the iterator advanced?
  }
}
END_SECTION


START_SECTION((static std::string& toUpper(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& firstToUpper(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& toLower(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& substitute(String &this_s, char from, char to)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& substitute(String &this_s, const std::string &from, const std::string &to)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& remove(String &this_s, char what)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& ensureLastChar(String &this_s, char end)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static std::string& removeWhitespaces(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION


START_SECTION([EXTRA] non-finite values round-trip through toStr/toDouble/toFloat)
{
  // std::to_chars writes "inf", which has neither '.' nor 'e', so the "keep at least one digit after
  // the decimal point" rule used to append ".0" and emit "inf.0" - a token the parsers reject. Any
  // file carrying an infinity then failed to load in full.
  const double dinf = std::numeric_limits<double>::infinity();
  const float finf = std::numeric_limits<float>::infinity();

  TEST_EQUAL(StringUtils::toStr(dinf), "inf")
  TEST_EQUAL(StringUtils::toStr(-dinf), "-inf")
  TEST_EQUAL(StringUtils::toStr(finf), "inf")
  TEST_EQUAL(StringUtils::toStr(-finf), "-inf")
  TEST_EQUAL(StringUtils::toStr(std::numeric_limits<double>::quiet_NaN()), "NaN")

  // and back again
  TEST_EQUAL(std::isinf(StringUtils::toDouble(StringUtils::toStr(dinf))), true)
  TEST_EQUAL(StringUtils::toDouble(StringUtils::toStr(dinf)) > 0, true)
  TEST_EQUAL(std::isinf(StringUtils::toDouble(StringUtils::toStr(-dinf))), true)
  TEST_EQUAL(StringUtils::toDouble(StringUtils::toStr(-dinf)) < 0, true)
  TEST_EQUAL(std::isinf(StringUtils::toFloat(StringUtils::toStr(finf))), true)
  TEST_EQUAL(StringUtils::toFloat(StringUtils::toStr(-finf)) < 0, true)
  TEST_EQUAL(std::isnan(StringUtils::toDouble(StringUtils::toStr(std::numeric_limits<double>::quiet_NaN()))), true)

  // finite values are untouched by the added branch
  TEST_EQUAL(StringUtils::toStr(5.0), "5.0")
  TEST_EQUAL(StringUtils::toStr(-0.5), "-0.5")
  TEST_EQUAL(StringUtils::toStr(0.0), "0.0")

  // a list of them survives the DataValue encoding used by the XML writers
  const DoubleList mixed {1.5, dinf, -dinf, std::numeric_limits<double>::quiet_NaN(), 2.5};
  const StringList parts = ListUtils::create<std::string>(
    StringUtils::substr(DataValue(mixed).toString(), 1, DataValue(mixed).toString().size() - 2), ',');
  TEST_EQUAL(parts.size(), 5)
  const DoubleList back = ListUtils::create<double>(parts);
  TEST_REAL_SIMILAR(back[0], 1.5)
  TEST_EQUAL(std::isinf(back[1]) && back[1] > 0, true)
  TEST_EQUAL(std::isinf(back[2]) && back[2] < 0, true)
  TEST_EQUAL(std::isnan(back[3]), true)
  TEST_REAL_SIMILAR(back[4], 2.5)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



