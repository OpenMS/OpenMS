// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------
//
// Tests for std::string + StringUtils free functions.
// OpenMS::String is now a type alias for std::string; all former member
// functions are available as StringUtils:: free functions.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/DATASTRUCTURES/StringUtils.h>    // brings in StringUtils + using String = std::string
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

using namespace OpenMS;          // pulls in StringUtils::toStr(= std::string), QuotingMethod, StringUtils
using namespace std;

START_TEST(String, "$Id$")

// -------------------------------------------------------------------------
// Construction from std::string, char*, char, size
// -------------------------------------------------------------------------

std::string* s_ptr = nullptr;
std::string* s_nullPointer = nullptr;
START_SECTION((std::string()))
  s_ptr = new std::string;
  TEST_NOT_EQUAL(s_ptr, s_nullPointer)
END_SECTION

START_SECTION(([EXTRA] ~std::string()))
  delete s_ptr;
END_SECTION

START_SECTION((std::string(const char* s, size_t length)))
  std::string s("abcdedfg", 5);
  TEST_EQUAL(s, "abcde")

  std::string s2("abcdedfg", 0);
  TEST_EQUAL(s2, "")

  std::string s3("abcdedfg", 8);
  TEST_EQUAL(s3, "abcdedfg")

  char test_utf8[]   = "T\xc3\xbc\x62ingen";
  char test_utf16[] {'\xff','\xfe','T','\x00','\xfc','\x00','b','\x00','i','\x00','n','\x00','g','\x00','e','\x00','n','\x00'};
  char test_iso8859[] = "T\xfc\x62ingen";

  std::string s_utf8(test_utf8, 9);
  std::string s_utf16(test_utf16, 18);
  std::string s_iso8859(test_iso8859, 8);

  TEST_EQUAL(s_utf16.size(), 18U)
  TEST_EQUAL(s_utf8.size(), 9U)
  TEST_EQUAL(s_iso8859.size(), 8U)
  TEST_EQUAL(s_iso8859[1], '\xfc')
  TEST_EQUAL(s_utf8[1], '\xc3')
  TEST_EQUAL(s_utf8[2], '\xbc')
END_SECTION

START_SECTION((std::string(const std::string&)))
  std::string s(std::string("blablabla"));
  TEST_EQUAL(s, "blablabla")
END_SECTION

START_SECTION((std::string(const char*)))
  std::string s("blablabla");
  TEST_EQUAL(s, "blablabla")
END_SECTION

START_SECTION((std::string(size_t, char)))
  std::string s(17, 'b');
  TEST_EQUAL(s, "bbbbbbbbbbbbbbbbb")
END_SECTION

// -------------------------------------------------------------------------
// toStr() — numeric to string
// -------------------------------------------------------------------------

START_SECTION((StringUtils::toStr(int)))
  TEST_EQUAL(StringUtils::toStr(int(-17)), "-17")
END_SECTION

START_SECTION((StringUtils::toStr(unsigned int)))
  TEST_EQUAL(StringUtils::toStr((unsigned int)(17)), "17")
END_SECTION

START_SECTION((StringUtils::toStr(long int)))
  TEST_EQUAL(StringUtils::toStr((long int)(-17)), "-17")
END_SECTION

START_SECTION((StringUtils::toStr(long unsigned int)))
  TEST_EQUAL(StringUtils::toStr((long unsigned int)(17)), "17")
END_SECTION

START_SECTION((StringUtils::toStr(short int)))
  TEST_EQUAL(StringUtils::toStr((short int)(-17)), "-17")
END_SECTION

START_SECTION((StringUtils::toStr(short unsigned int)))
  TEST_EQUAL(StringUtils::toStr((short unsigned int)(17)), "17")
END_SECTION

START_SECTION((StringUtils::toStr(long long unsigned int)))
  TEST_EQUAL(StringUtils::toStr((long long unsigned int)(12345678)), "12345678")
END_SECTION

START_SECTION((StringUtils::toStr(long long signed int)))
  TEST_EQUAL(StringUtils::toStr((long long signed int)(-12345678)), "-12345678")
END_SECTION

START_SECTION((StringUtils::toStr(float, bool full_precision)))
  TEST_EQUAL(StringUtils::toStr(17.0123456f),         "17.012346")
  TEST_EQUAL(StringUtils::toStr(17.0123f, false),     "17.012")
  TEST_EQUAL(StringUtils::toStr(50254.199219f),       "5.02542e04")

  constexpr float denorm = std::numeric_limits<float>::min() / 10;
  TEST_EQUAL(StringUtils::toStr(denorm, true),  "1.175495e-39")
  TEST_EQUAL(StringUtils::toStr(denorm, false), "1.175e-39")

  constexpr float nan_f = std::numeric_limits<float>::quiet_NaN();
  TEST_EQUAL(StringUtils::toStr(nan_f, true),  "NaN")
  TEST_EQUAL(StringUtils::toStr(nan_f, false), "NaN")
END_SECTION

START_SECTION((StringUtils::toStr(double, bool full_precision)))
  TEST_EQUAL(StringUtils::toStr(double(17.012345)),         "17.012345")
  TEST_EQUAL(StringUtils::toStr(double(17.012345), false),  "17.012")

  constexpr double denorm = std::numeric_limits<double>::min() / 10;
  TEST_EQUAL(StringUtils::toStr(denorm, true),  "2.225073858507203e-309")
  TEST_EQUAL(StringUtils::toStr(denorm, false), "2.225e-309")

  constexpr double nan_d = std::numeric_limits<double>::quiet_NaN();
  TEST_EQUAL(StringUtils::toStr(nan_d, true),  "NaN")
  TEST_EQUAL(StringUtils::toStr(nan_d, false), "NaN")
END_SECTION

START_SECTION((StringUtils::toStr(long double, bool full_precision)))
  TEST_EQUAL(StringUtils::toStr(17.012345L),        "17.012345")
  TEST_EQUAL(StringUtils::toStr(17.012345L, false), "17.012")

  constexpr long double nan_ld = std::numeric_limits<long double>::quiet_NaN();
  TEST_EQUAL(StringUtils::toStr(nan_ld, true),  "NaN")
  TEST_EQUAL(StringUtils::toStr(nan_ld, false), "NaN")
END_SECTION

START_SECTION((StringUtils::toStr(DataValue)))
  TEST_EQUAL(StringUtils::toStr(DataValue(17.012345)),              "17.012345")
  TEST_EQUAL(StringUtils::toStr(DataValue(17.012345), false),       "17.012")
  TEST_EQUAL(StringUtils::toStr(DataValue(DoubleList({17.012345, 2.0}))),        "[17.012345, 2.0]")
  TEST_EQUAL(StringUtils::toStr(DataValue(DoubleList({17.012345, 2.0})), false), "[17.012, 2.0]")
  TEST_EQUAL(StringUtils::toStr(DataValue("bla")),  "bla")
  TEST_EQUAL(StringUtils::toStr(DataValue(4711)),   "4711")
END_SECTION

// -------------------------------------------------------------------------
// StringUtils::numberLength / number / random
// -------------------------------------------------------------------------

START_SECTION((StringUtils::numberLength(double, UInt)))
  TEST_EQUAL(StringUtils::numberLength(12345678.9123,  11), "12345678.91")
  TEST_EQUAL(StringUtils::numberLength(-12345678.9123, 11), "-12345678.9")
  TEST_EQUAL(StringUtils::numberLength(12345678.9123,  10), "12345678.9")
  TEST_EQUAL(StringUtils::numberLength(-12345678.9123, 10), "-1234.5e04")
  TEST_EQUAL(StringUtils::numberLength(12345678.9123,   9), "1234.5e04")
  TEST_EQUAL(StringUtils::numberLength(-12345678.9123,  9), "-123.4e05")
END_SECTION

START_SECTION((StringUtils::number(double, UInt)))
  TEST_EQUAL(StringUtils::number(123.1234, 0), "123")
  TEST_EQUAL(StringUtils::number(123.1234, 1), "123.1")
  TEST_EQUAL(StringUtils::number(123.1234, 2), "123.12")
  TEST_EQUAL(StringUtils::number(123.1234, 3), "123.123")
  TEST_EQUAL(StringUtils::number(123.1234, 4), "123.1234")
  TEST_EQUAL(StringUtils::number(123.1234, 5), "123.12340")
  TEST_EQUAL(StringUtils::number(0.0,      5), "0.00000")
END_SECTION

START_SECTION((StringUtils::random(UInt)))
  std::string r = StringUtils::random(10);
  TEST_EQUAL(r.size(), 10U)
END_SECTION

// -------------------------------------------------------------------------
// Predicates
// -------------------------------------------------------------------------

std::string s("ACDEFGHIKLMNPQRSTVWY");

START_SECTION((StringUtils::hasPrefix))
  TEST_EQUAL(StringUtils::hasPrefix(s, ""),                        true)
  TEST_EQUAL(StringUtils::hasPrefix(s, "ACDEF"),                   true)
  TEST_EQUAL(StringUtils::hasPrefix(s, "ACDEFGHIKLMNPQRSTVWY"),    true)
  TEST_EQUAL(StringUtils::hasPrefix(s, "ABCDEF"),                  false)
  TEST_EQUAL(StringUtils::hasPrefix(s, "ACDEFGHIKLMNPQRSTVWYACDEF"), false)
END_SECTION

START_SECTION((StringUtils::hasSuffix))
  TEST_EQUAL(StringUtils::hasSuffix(s, ""),                        true)
  TEST_EQUAL(StringUtils::hasSuffix(s, "TVWY"),                    true)
  TEST_EQUAL(StringUtils::hasSuffix(s, "ACDEFGHIKLMNPQRSTVWY"),    true)
  TEST_EQUAL(StringUtils::hasSuffix(s, "WXYZ"),                    false)
  TEST_EQUAL(StringUtils::hasSuffix(s, "ACDEFACDEFGHIKLMNPQRSTVWY"), false)
END_SECTION

START_SECTION((StringUtils::hasSubstring))
  TEST_EQUAL(StringUtils::hasSubstring(s, ""),                       true)
  TEST_EQUAL(StringUtils::hasSubstring(s, "GHIKLM"),                 true)
  TEST_EQUAL(StringUtils::hasSubstring(s, "ACDEFGHIKLMNPQRSTVWY"),   true)
  TEST_EQUAL(StringUtils::hasSubstring(s, "MLKIGH"),                 false)
  TEST_EQUAL(StringUtils::hasSubstring(s, "ACDEFGHIKLMNPQRSTVWYACDEF"), false)
END_SECTION

START_SECTION((StringUtils::has))
  TEST_EQUAL(StringUtils::has(s, 'A'), true)
  TEST_EQUAL(StringUtils::has(s, 'O'), false)
END_SECTION

// -------------------------------------------------------------------------
// Accessors
// -------------------------------------------------------------------------

START_SECTION((StringUtils::prefix(Int)))
  TEST_EQUAL(StringUtils::prefix(s, (Int)4), "ACDE")
  TEST_EQUAL(StringUtils::prefix(s, (Int)0), "")
  TEST_EXCEPTION(Exception::IndexOverflow, StringUtils::prefix(s, (Int)(s.size()+1)))
  TEST_EXCEPTION(Exception::IndexUnderflow, StringUtils::prefix(s, (Int)(-1)))
END_SECTION

START_SECTION((StringUtils::suffix(Int)))
  TEST_EQUAL(StringUtils::suffix(s, (Int)4), "TVWY")
  TEST_EQUAL(StringUtils::suffix(s, (Int)0), "")
  TEST_EXCEPTION(Exception::IndexOverflow, StringUtils::suffix(s, (Int)(s.size()+1)))
  TEST_EXCEPTION(Exception::IndexUnderflow, StringUtils::suffix(s, (Int)(-1)))
END_SECTION

START_SECTION((StringUtils::prefix(size_t)))
  TEST_EQUAL(StringUtils::prefix(s, (size_t)4), "ACDE")
  TEST_EQUAL(StringUtils::prefix(s, (size_t)0), "")
  TEST_EXCEPTION(Exception::IndexOverflow, StringUtils::prefix(s, s.size()+1))
END_SECTION

START_SECTION((StringUtils::suffix(size_t)))
  TEST_EQUAL(StringUtils::suffix(s, (size_t)4), "TVWY")
  TEST_EQUAL(StringUtils::suffix(s, (size_t)0), "")
  TEST_EXCEPTION(Exception::IndexOverflow, StringUtils::suffix(s, s.size()+1))
END_SECTION

START_SECTION((StringUtils::prefix(char delim)))
  TEST_EQUAL(StringUtils::prefix(s, 'F'), "ACDE")
  TEST_EQUAL(StringUtils::prefix(s, 'A'), "")
  TEST_EXCEPTION(Exception::ElementNotFound, StringUtils::prefix(s, 'Z'))
END_SECTION

START_SECTION((StringUtils::suffix(char delim)))
  TEST_EQUAL(StringUtils::suffix(s, 'S'), "TVWY")
  TEST_EQUAL(StringUtils::suffix(s, 'Y'), "")
  TEST_EXCEPTION(Exception::ElementNotFound, StringUtils::suffix(s, 'Z'))
END_SECTION

START_SECTION((StringUtils::substr))
  std::string s2("abcdef");
  TEST_EQUAL(StringUtils::substr(s2, 0, 4), "abcd")
  TEST_EQUAL(StringUtils::substr(s2, 1, 1), "b")
  TEST_EQUAL(StringUtils::substr(s2, 1, 3), "bcd")
  TEST_EQUAL(StringUtils::substr(s2, 0, 6), "abcdef")
  TEST_EQUAL(StringUtils::substr(s2, 5, 1), "f")
  TEST_EQUAL(StringUtils::substr(s2, 6, 1), "")
  TEST_EQUAL(StringUtils::substr(s2, 0, 7), "abcdef")
  TEST_EQUAL(StringUtils::substr(s2, 0, std::string::npos), "abcdef")
  TEST_EQUAL(StringUtils::substr(s2, 0), "abcdef")
  TEST_EQUAL(StringUtils::substr(s2, 1), "bcdef")
  TEST_EQUAL(StringUtils::substr(s2, 5), "f")
  TEST_EQUAL(StringUtils::substr(s2, 6), "")
END_SECTION

START_SECTION((StringUtils::chop))
  std::string s2("abcdef");
  TEST_EQUAL(StringUtils::chop(s2, 0), "abcdef")
  TEST_EQUAL(StringUtils::chop(s2, 1), "abcde")
  TEST_EQUAL(StringUtils::chop(s2, 2), "abcd")
  TEST_EQUAL(StringUtils::chop(s2, 3), "abc")
  TEST_EQUAL(StringUtils::chop(s2, 4), "ab")
  TEST_EQUAL(StringUtils::chop(s2, 5), "a")
  TEST_EQUAL(StringUtils::chop(s2, 6), "")
  TEST_EQUAL(StringUtils::chop(s2, 9), "")
  TEST_EQUAL(StringUtils::chop(s2, (Size)(-1)), "") // wraps to huge -> empty
END_SECTION

// -------------------------------------------------------------------------
// Mutators
// -------------------------------------------------------------------------

START_SECTION((StringUtils::reverse))
  std::string r("ACDEFGHIKLMNPQRSTVWY");
  StringUtils::reverse(r);
  TEST_EQUAL(r, "YWVTSRQPNMLKIHGFEDCA")
  r.clear();
  StringUtils::reverse(r);
  TEST_EQUAL(r, "")
END_SECTION

START_SECTION((StringUtils::trim))
  std::string t("\n\r\t test \n\r\t");
  StringUtils::trim(t);
  TEST_EQUAL(t, "test")
  StringUtils::trim(t);
  TEST_EQUAL(t, "test")
  t = "";
  StringUtils::trim(t);
  TEST_EQUAL(t, "")
  t = " t";
  StringUtils::trim(t);
  TEST_EQUAL(t, "t")
  t = "t ";
  StringUtils::trim(t);
  TEST_EQUAL(t, "t")
  t = "\t\r\n ";
  StringUtils::trim(t);
  TEST_EQUAL(t, "")
END_SECTION

START_SECTION((StringUtils::quote / StringUtils::unquote))
  std::string q;
  StringUtils::quote(q, '\'', QuotingMethod::NONE);
  TEST_EQUAL(q, "''")
  StringUtils::quote(q, '\'', QuotingMethod::ESCAPE);
  TEST_EQUAL(q, "'\\'\\''")

  q = "ab\"cd\\ef";
  StringUtils::quote(q, '"', QuotingMethod::NONE);
  TEST_EQUAL(q, "\"ab\"cd\\ef\"")
  StringUtils::quote(q, '"', QuotingMethod::ESCAPE);
  TEST_EQUAL(q, "\"\\\"ab\\\"cd\\\\ef\\\"\"")

  q = "ab\"cd\\ef";
  StringUtils::quote(q, '"', QuotingMethod::DOUBLE);
  TEST_EQUAL(q, "\"ab\"\"cd\\ef\"")

  std::string u;
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::unquote(u))
  u = "''";
  StringUtils::unquote(u, '\'', QuotingMethod::NONE);
  TEST_EQUAL(u, "")
  u = "'\\'\\''";
  StringUtils::unquote(u, '\'', QuotingMethod::ESCAPE);
  TEST_EQUAL(u, "''")
  u = R"("ab"cd\ef")";
  StringUtils::unquote(u, '"', QuotingMethod::NONE);
  TEST_EQUAL(u, "ab\"cd\\ef")
  u = R"("\"ab\"cd\\ef\"")";
  StringUtils::unquote(u, '"', QuotingMethod::ESCAPE);
  TEST_EQUAL(u, "\"ab\"cd\\ef\"")
  u = R"("ab""cd\ef")";
  StringUtils::unquote(u, '"', QuotingMethod::DOUBLE);
  TEST_EQUAL(u, "ab\"cd\\ef")
END_SECTION

START_SECTION((StringUtils::simplify))
  std::string si("\n\r\t te\tst \n\r\t");
  StringUtils::simplify(si);
  TEST_EQUAL(si, " te st ")
  StringUtils::simplify(si);
  TEST_EQUAL(si, " te st ")
  si = "";
  StringUtils::simplify(si);
  TEST_EQUAL(si, "")
  si = " t";
  StringUtils::simplify(si);
  TEST_EQUAL(si, " t")
  si = "t ";
  StringUtils::simplify(si);
  TEST_EQUAL(si, "t ")
  si = "\t\r\n ";
  StringUtils::simplify(si);
  TEST_EQUAL(si, " ")
END_SECTION

START_SECTION((StringUtils::fillLeft))
  std::string fl("TEST");
  StringUtils::fillLeft(fl, 'x', 4);
  TEST_EQUAL(fl, "TEST")
  StringUtils::fillLeft(fl, 'y', 5);
  TEST_EQUAL(fl, "yTEST")
  StringUtils::fillLeft(fl, 'z', 7);
  TEST_EQUAL(fl, "zzyTEST")
END_SECTION

START_SECTION((StringUtils::fillRight))
  std::string fr("TEST");
  StringUtils::fillRight(fr, 'x', 4);
  TEST_EQUAL(fr, "TEST")
  StringUtils::fillRight(fr, 'y', 5);
  TEST_EQUAL(fr, "TESTy")
  StringUtils::fillRight(fr, 'z', 7);
  TEST_EQUAL(fr, "TESTyzz")
END_SECTION

// -------------------------------------------------------------------------
// Converters
// -------------------------------------------------------------------------

START_SECTION((StringUtils::toInt32))
  std::string sv;
  sv = "123";
  TEST_EQUAL(StringUtils::toInt32(sv), 123)
  sv = "-123";
  TEST_EQUAL(StringUtils::toInt32(sv), -123)
  sv = "  -123  ";
  TEST_EQUAL(StringUtils::toInt32(sv), -123)
  sv = "524 starts with an int";
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt32(sv))
  sv = "not an int";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, StringUtils::toInt32(sv),
    "Could not convert string 'not an int' to an integer value")
END_SECTION

START_SECTION((StringUtils::toFloat))
  std::string sv;
  sv = "123.456";
  TEST_REAL_SIMILAR(StringUtils::toFloat(sv), 123.456f)
  sv = "-123.456";
  TEST_REAL_SIMILAR(StringUtils::toFloat(sv), -123.456f)
  sv = "73629.9";
  TEST_EQUAL(StringUtils::toStr(StringUtils::toFloat(sv)), "7.36299e04")
  sv = "nan";
  TEST_EQUAL(std::isnan(StringUtils::toFloat(sv)), true)
  sv = "NaN";
  TEST_EQUAL(std::isnan(StringUtils::toFloat(sv)), true)
  sv = "not a number";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, StringUtils::toFloat(sv),
    "Could not convert string 'not a number' to a float value")
END_SECTION

START_SECTION((StringUtils::toDouble))
  std::string sv;
  sv = "123.456";
  TEST_REAL_SIMILAR(StringUtils::toDouble(sv), 123.456)
  sv = "-123.4567890123";
  TEST_REAL_SIMILAR(StringUtils::toDouble(sv), -123.4567890123)
  sv = "73629.980123";
  TEST_EQUAL(StringUtils::toStr(StringUtils::toDouble(sv)), "7.3629980123e04")
  sv = "nan";
  TEST_TRUE(std::isnan(StringUtils::toDouble(sv)))
  sv = "NaN";
  TEST_TRUE(std::isnan(StringUtils::toDouble(sv)))
  sv = "not a number";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, StringUtils::toDouble(sv),
    "Could not convert string 'not a number' to a double value")
END_SECTION

// -------------------------------------------------------------------------
// substitute / remove / ensureLastChar / removeWhitespaces / toUpper / toLower
// -------------------------------------------------------------------------

START_SECTION((StringUtils::toUpper))
  std::string tu = "test45%#.,";
  StringUtils::toUpper(tu);
  TEST_EQUAL(tu, "TEST45%#.,")
  tu = "";
  StringUtils::toUpper(tu);
  TEST_EQUAL(tu, "")
END_SECTION

START_SECTION((StringUtils::toLower))
  std::string tl = "TEST45%#.,";
  StringUtils::toLower(tl);
  TEST_EQUAL(tl, "test45%#.,")
  tl = "";
  StringUtils::toLower(tl);
  TEST_EQUAL(tl, "")
END_SECTION

START_SECTION((StringUtils::firstToUpper))
  std::string fu = "test45%#.,";
  StringUtils::firstToUpper(fu);
  TEST_EQUAL(fu, "Test45%#.,")
  fu = " ";
  StringUtils::firstToUpper(fu);
  TEST_EQUAL(fu, " ")
  fu = "";
  StringUtils::firstToUpper(fu);
  TEST_EQUAL(fu, "")
END_SECTION

START_SECTION((StringUtils::substitute char->char))
  std::string sub = "abcdefg";
  StringUtils::substitute(sub, 'a', 'x');
  TEST_EQUAL(sub, "xbcdefg")
  StringUtils::substitute(sub, 'g', 'y');
  TEST_EQUAL(sub, "xbcdefy")
  StringUtils::substitute(sub, 'c', '-');
  TEST_EQUAL(sub, "xb-defy")
  sub = ".....";
  StringUtils::substitute(sub, '.', ',');
  TEST_EQUAL(sub, ",,,,,")
END_SECTION

START_SECTION((StringUtils::substitute string->string))
  std::string sub = "abcdefg";
  StringUtils::substitute(sub, std::string("a"), std::string("x"));
  TEST_EQUAL(sub, "xbcdefg")
  StringUtils::substitute(sub, std::string("bcd"), std::string("y"));
  TEST_EQUAL(sub, "xyefg")
  StringUtils::substitute(sub, std::string("fg"), std::string(""));
  TEST_EQUAL(sub, "xye")
  StringUtils::substitute(sub, std::string("e"), std::string("z!"));
  TEST_EQUAL(sub, "xyz!")
  StringUtils::substitute(sub, std::string("u"), std::string("blblblblbl"));
  TEST_EQUAL(sub, "xyz!")
  StringUtils::substitute(sub, std::string(""), std::string("blblblblbl"));
  TEST_EQUAL(sub, "xyz!")
  sub = "abcdefgabcdefgabcdefgab";
  StringUtils::substitute(sub, std::string("ab"), std::string("x"));
  TEST_EQUAL(sub, "xcdefgxcdefgxcdefgx")
  StringUtils::substitute(sub, std::string("x"), std::string(""));
  TEST_EQUAL(sub, "cdefgcdefgcdefg")
END_SECTION

START_SECTION((StringUtils::remove))
  std::string rm = "abcabc";
  StringUtils::remove(rm, 'a');
  TEST_EQUAL(rm, "bcbc")
  StringUtils::remove(rm, 'c');
  TEST_EQUAL(rm, "bb")
  StringUtils::remove(rm, 'b');
  TEST_EQUAL(rm, "")
END_SECTION

START_SECTION((StringUtils::ensureLastChar))
  std::string elc = "/";
  StringUtils::ensureLastChar(elc, '/');
  TEST_EQUAL(elc, "/")
  StringUtils::ensureLastChar(elc, '\\');
  TEST_EQUAL(elc, "/\\")
  StringUtils::ensureLastChar(elc, '\\');
  TEST_EQUAL(elc, "/\\")
  StringUtils::ensureLastChar(elc, '/');
  TEST_EQUAL(elc, "/\\/")
END_SECTION

START_SECTION((StringUtils::removeWhitespaces))
{
  std::string rw;
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, "")

  rw = " \t \n ";
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, "")

  rw = "test";
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, "test")

  rw = "\n\r\t test \n\r\t";
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, "test")

  rw = "\n\r\t t\ne \ts\rt \n\r\t";
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, "test")

  const std::string test(16 * 1024 + 1, 'x');
  rw = test + std::string(100, ' ');
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(rw.begin(), rw.end(), g);
  StringUtils::removeWhitespaces(rw);
  TEST_EQUAL(rw, test)
}
END_SECTION

// -------------------------------------------------------------------------
// Global operator+ and operator+= (std::string + numeric types)
// -------------------------------------------------------------------------

const std::string fixed("test");

START_SECTION((operator+(std::string, int)))
  TEST_EQUAL(fixed + (int)(4),                    "test4")
END_SECTION

START_SECTION((operator+(std::string, unsigned int)))
  TEST_EQUAL(fixed + (unsigned int)(4),            "test4")
END_SECTION

START_SECTION((operator+(std::string, short int)))
  TEST_EQUAL(fixed + (short int)(4),               "test4")
END_SECTION

START_SECTION((operator+(std::string, short unsigned int)))
  TEST_EQUAL(fixed + (short unsigned int)(4),      "test4")
END_SECTION

START_SECTION((operator+(std::string, long int)))
  TEST_EQUAL(fixed + (long int)(4),                "test4")
END_SECTION

START_SECTION((operator+(std::string, long unsigned int)))
  TEST_EQUAL(fixed + (long unsigned int)(4),       "test4")
END_SECTION

START_SECTION((operator+(std::string, long long unsigned int)))
  TEST_EQUAL(fixed + (long long unsigned int)(4),  "test4")
END_SECTION

START_SECTION((operator+(std::string, float)))
  TEST_EQUAL(fixed + (float)(4),                   "test4.0")
END_SECTION

START_SECTION((operator+(std::string, double)))
  TEST_EQUAL(fixed + (double)(4),                  "test4.0")
END_SECTION

START_SECTION((operator+(std::string, long double)))
  TEST_EQUAL(fixed + (long double)(4),             "test4.0")
END_SECTION

START_SECTION((operator+=(std::string&, int)))
  std::string s2 = "test";
  s2 += (int)(7);
  TEST_EQUAL(s2, "test7")
END_SECTION

START_SECTION((operator+=(std::string&, unsigned int)))
  std::string s2 = "test";
  s2 += (unsigned int)(7);
  TEST_EQUAL(s2, "test7")
END_SECTION

START_SECTION((operator+=(std::string&, float)))
  std::string s2 = "test";
  s2 += (float)(7.4);
  TEST_EQUAL(s2, "test7.4")
END_SECTION

START_SECTION((operator+=(std::string&, double)))
  std::string s2 = "test";
  s2 += (double)(7.4);
  TEST_EQUAL(s2, "test7.4")
END_SECTION

START_SECTION((operator+=(std::string&, long double)))
{
  std::string s2 = "test";
  long double x = 7.4L;
  s2 += x;
  TEST_EQUAL(s2, "test7.4")
}
END_SECTION

// -------------------------------------------------------------------------
// split / split_quoted / concatenate
// -------------------------------------------------------------------------

START_SECTION((StringUtils::split char))
  std::string sp(";1;2;3;4;5;");
  vector<std::string> parts;
  bool result = StringUtils::split(sp, ';', parts);
  TEST_EQUAL(result, true)
  TEST_EQUAL(parts.size(), 7U)
  TEST_EQUAL(parts[0], "")
  TEST_EQUAL(parts[1], "1")
  TEST_EQUAL(parts[6], "")

  sp = "1;2;3;4;5";
  result = StringUtils::split(sp, ';', parts);
  TEST_EQUAL(result, true)
  TEST_EQUAL(parts.size(), 5U)

  sp = "";
  result = StringUtils::split(sp, ',', parts);
  TEST_EQUAL(result, false)
  TEST_EQUAL(parts.size(), 0U)

  sp = ";";
  result = StringUtils::split(sp, ';', parts);
  TEST_EQUAL(result, true)
  TEST_EQUAL(parts.size(), 2U)

  sp = "nodelim";
  result = StringUtils::split(sp, ';', parts);
  TEST_EQUAL(result, false)
  TEST_EQUAL(parts.size(), 1U)

  // quote_protect
  sp = " \"hello\", world, 23.3";
  result = StringUtils::split(sp, ',', parts, true);
  TEST_EQUAL(result, true)
  TEST_EQUAL(parts.size(), 3U)
  TEST_EQUAL(parts[0], "hello")
  TEST_EQUAL(parts[1], "world")
  TEST_EQUAL(parts[2], "23.3")

  sp = R"( "hello", " donot,splitthis ", "23.4 " )";
  result = StringUtils::split(sp, ',', parts, true);
  TEST_EQUAL(result, true)
  TEST_EQUAL(parts.size(), 3U)
  TEST_EQUAL(parts[1], " donot,splitthis ")

  sp = R"( "first", "seconds"<thisshouldnotbehere>, third)";
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::split(sp, ',', parts, true))
END_SECTION

START_SECTION((StringUtils::split string))
  std::string sp2 = "abcdebcfghbc";
  vector<std::string> subs;
  bool result = StringUtils::split(sp2, std::string("bc"), subs);
  TEST_EQUAL(result, true)
  TEST_EQUAL(subs.size(), 4U)
  TEST_EQUAL(subs[0], "a")
  TEST_EQUAL(subs[1], "de")
  TEST_EQUAL(subs[2], "fgh")
  TEST_EQUAL(subs[3], "")

  result = StringUtils::split(sp2, std::string("xy"), subs);
  TEST_EQUAL(result, false)
  TEST_EQUAL(subs.size(), 1U)

  result = StringUtils::split(sp2, std::string(""), subs);
  TEST_EQUAL(result, true)
  TEST_EQUAL(subs.size(), sp2.size())

  result = StringUtils::split(std::string(""), std::string(","), subs);
  TEST_EQUAL(result, false)
  TEST_EQUAL(subs.size(), 0U)
END_SECTION

START_SECTION((StringUtils::split_quoted))
  std::string sq = "abcdebcfghbc";
  vector<std::string> subs;
  bool result = StringUtils::split_quoted(sq, "bc", subs);
  TEST_EQUAL(result, true)
  TEST_EQUAL(subs.size(), 4U)

  sq = R"("a,b,c","d,\",f","")";
  result = StringUtils::split_quoted(sq, ",", subs, '"', QuotingMethod::ESCAPE);
  TEST_EQUAL(result, true)
  TEST_EQUAL(subs.size(), 3U)
  TEST_EQUAL(subs[0], "\"a,b,c\"")
  TEST_EQUAL(subs[1], "\"d,\\\",f\"")
  TEST_EQUAL(subs[2], "\"\"")

  sq = R"("a,"b")";
  TEST_EXCEPTION(Exception::ConversionError,
    StringUtils::split_quoted(sq, ",", subs, '"', QuotingMethod::ESCAPE))

  sq = R"("ab"___"cd""ef")";
  result = StringUtils::split_quoted(sq, "___", subs, '"', QuotingMethod::DOUBLE);
  TEST_EQUAL(result, true)
  TEST_EQUAL(subs.size(), 2U)
  TEST_EQUAL(subs[0], "\"ab\"")
  TEST_EQUAL(subs[1], "\"cd\"\"ef\"")
END_SECTION

START_SECTION((StringUtils::concatenate))
  vector<std::string> parts;
  StringUtils::split(std::string("1;2;3;4;5"), ';', parts);
  std::string joined = StringUtils::concatenate(parts, "g");
  TEST_EQUAL(joined, "1g2g3g4g5")

  StringUtils::split(std::string("1;2;3;4;5"), ';', parts);
  joined = StringUtils::concatenate(parts);
  TEST_EQUAL(joined, "12345")

  StringUtils::split(std::string(""), ';', parts);
  joined = StringUtils::concatenate(parts);
  TEST_EQUAL(joined, "")
END_SECTION

END_TEST
