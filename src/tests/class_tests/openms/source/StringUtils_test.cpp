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
#include <OpenMS/DATASTRUCTURES/StringUtilsSimple.h>
///////////////////////////

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
    String at1 = "0 2  3456789101112" x16;
    at1.substitute(' ', whitespace);
    TEST_EQUAL(skipWhitespace(at1), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 1)), 1);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 2)), 0);
    TEST_EQUAL(skipWhitespace(std::string_view(at1.data() + 3)), 2);
    String at2 = s16 s16 "1" x16;
    at2.substitute(' ', whitespace);
    TEST_EQUAL(skipWhitespace(std::string_view(at2.data())), 32);
    TEST_EQUAL(skipWhitespace(std::string_view(at2.data() + 2)), 30);
    String at1_noSSE = "0 2  34";
    at1_noSSE.substitute(' ', whitespace);
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
    String at1 = "0 2  3456789101112" x16;
    at1.substitute(' ', whitespace);
    TEST_EQUAL(skipNonWhitespace(at1), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 1)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 2)), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 3)), 0);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at1.data() + 5)), 13 + 16);
    String at2 = x16 x16 " " x16;
    at2.substitute(' ', whitespace);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data())), 32);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data() + 31)), 1);
    TEST_EQUAL(skipNonWhitespace(std::string_view(at2.data() + 33)), 16);
    String at1_noSSE = "0 2  34";
    at1_noSSE.substitute(' ', whitespace);
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

START_SECTION((static String numberLength(double d, UInt n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String number(double d, UInt n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& fillLeft(String &this_s, char c, UInt size)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& fillRight(String &this_s, char c, UInt size)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasPrefix(const String &this_s, const String &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasSuffix(const String &this_s, const String &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool hasSubstring(const String &this_s, const String &string)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool has(const String &this_s, Byte byte)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, size_t length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, size_t length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, Int length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, Int length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, char delim)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, char delim)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String substr(const String &this_s, size_t pos, size_t n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String chop(const String &this_s, Size n)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& trim(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& quote(String &this_s, char q, String::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& unquote(String &this_s, char q, String::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& simplify(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String random(UInt length)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& reverse(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split(const String &this_s, const char splitter, std::vector< String > &substrings, bool quote_protect)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split(const String &this_s, const String &splitter, std::vector< String > &substrings)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static bool split_quoted(const String &this_s, const String &splitter, std::vector< String > &substrings, char q, String::QuotingMethod method)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static QString toQString(const String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION


START_SECTION((static Int32 toInt32(const String &this_s)))
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



START_SECTION((static Int64 toInt64(const String &this_s)))
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

START_SECTION((static float toFloat(const String &this_s)))
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
}
END_SECTION

START_SECTION((static double toDouble(const String &this_s)))
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
}
END_SECTION

START_SECTION((static String& toUpper(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& firstToUpper(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& toLower(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& substitute(String &this_s, char from, char to)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& substitute(String &this_s, const String &from, const String &to)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& remove(String &this_s, char what)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& ensureLastChar(String &this_s, char end)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION

START_SECTION((static String& removeWhitespaces(String &this_s)))
{
  NOT_TESTABLE // tested in String_test.cpp
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



