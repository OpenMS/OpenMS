// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

#include <algorithm>
#include <cmath>

#include <limits>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

String* s_ptr = nullptr;
String* s_nullPointer = nullptr;
START_SECTION((String()))
  s_ptr = new String;
  TEST_NOT_EQUAL(s_ptr, s_nullPointer)
END_SECTION

START_SECTION(([EXTRA] ~String()))
  delete s_ptr;
END_SECTION

START_SECTION((String(const QString &s)))
  QString qs("bla");
  String s(qs);
  TEST_EQUAL(s=="bla",true)
END_SECTION

START_SECTION((QString toQString() const))
  QString qs("bla");
  String s("bla");
  TEST_EQUAL(s.toQString()==qs,true)
END_SECTION

START_SECTION((String(const char* s, SizeType length)))
  String s("abcdedfg",5);
  TEST_EQUAL(s,"abcde")

  String s2("abcdedfg",0);
  TEST_EQUAL(s2,"")

  String s3("abcdedfg",8);
  TEST_EQUAL(s3,"abcdedfg")

  // This function does *not* check for null bytes as these may be part of a
  // valid string! If you provide a length that is inaccurate, you are to blame
  // (we do not test this since we would cause undefined behavior).
  // String s4("abcdedfg", 15);
  // TEST_NOT_EQUAL(s4,"abcdedfg")

  // Unicode test
  // print(b"T\xc3\xbc\x62ingen".decode("utf8"))
  // print(b"\xff\xfeT\x00\xfc\x00b\x00i\x00n\x00g\x00e\x00n\x00".decode("utf16"))
  // print(b"T\xfc\x62ingen".decode("iso8859"))
  char test_utf8[] =  "T\xc3\xbc\x62ingen";
  char test_utf16[] {'\xff','\xfe','T','\x00','\xfc','\x00','b','\x00','i','\x00','n','\x00','g','\x00','e','\x00','n','\x00'};
  char test_iso8859[] = "T\xfc\x62ingen";

  String s_utf8(test_utf8, 9);
  String s_utf16(test_utf16, 18);
  String s_iso8859(test_iso8859, 8);

  TEST_EQUAL(s_utf16.size(), 18);
  TEST_EQUAL(s_utf8.size(), 9);
  TEST_EQUAL(s_iso8859.size(), 8);

  TEST_EQUAL(s_iso8859[0], 'T')
  TEST_EQUAL(s_iso8859[1], '\xfc') // single byte for u
  TEST_EQUAL(s_iso8859[2], 'b')
  TEST_EQUAL(s_utf8[0], 'T')
  TEST_EQUAL(s_utf8[1], '\xc3')
  TEST_EQUAL(s_utf8[2], '\xbc') // two bytes for u
  TEST_EQUAL(s_utf8[3], 'b')

  // Its very easy to produce nonsense with unicode (e.g. chop string in half
  // in the middle of a character)
  String nonsense(test_utf8, 2);
  TEST_EQUAL(nonsense.size(), 2);
  TEST_EQUAL(nonsense[1], '\xc3')
END_SECTION


START_SECTION((String(const std::string& s)))
  String s(string("blablabla"));
  TEST_EQUAL(s,"blablabla")

  // Unicode test
  // print(b"T\xc3\xbc\x62ingen".decode("utf8"))
  // print(b"\xff\xfeT\x00\xfc\x00b\x00i\x00n\x00g\x00e\x00n\x00".decode("utf16"))
  // print(b"T\xfc\x62ingen".decode("iso8859"))
  char test_utf8[] =  "T\xc3\xbc\x62ingen";
  char test_utf16[] {'\xff','\xfe','T','\x00','\xfc','\x00','b','\x00','i','\x00','n','\x00','g','\x00','e','\x00','n','\x00'};
  char test_iso8859[] = "T\xfc\x62ingen";

  std::string std_s_utf8(test_utf8, 9);
  std::string std_s_utf16(test_utf16, 18);
  std::string std_s_iso8859(test_iso8859, 8);

  String s_utf8(std_s_utf8);
  String s_utf16(std_s_utf16);
  String s_iso8859(std_s_iso8859);

  TEST_EQUAL(s_utf16.size(), 18);
  TEST_EQUAL(s_utf8.size(), 9);
  TEST_EQUAL(s_iso8859.size(), 8);

  TEST_EQUAL(s_iso8859[0], 'T')
  TEST_EQUAL(s_iso8859[1], '\xfc') // single byte for u
  TEST_EQUAL(s_iso8859[2], 'b')
  TEST_EQUAL(s_utf8[0], 'T')
  TEST_EQUAL(s_utf8[1], '\xc3')
  TEST_EQUAL(s_utf8[2], '\xbc') // two bytes for u
  TEST_EQUAL(s_utf8[3], 'b')
END_SECTION

START_SECTION((String(const char* s)))
  String s("blablabla");
  TEST_EQUAL(s,"blablabla")
END_SECTION

START_SECTION((String(size_t len, char c)))
  String s(17,'b');
  TEST_EQUAL(s,"bbbbbbbbbbbbbbbbb")
END_SECTION

START_SECTION((String(const char c)))
  String s('v');
  TEST_EQUAL(s,"v")
END_SECTION

START_SECTION((String(int i)))
  String s(int (-17));
  TEST_EQUAL(s,"-17")
END_SECTION

START_SECTION((String(unsigned int i)))
  String s((unsigned int) (17));
  TEST_EQUAL(s,"17")
END_SECTION

START_SECTION((String(long int i)))
  String s((long int)(-17));
  TEST_EQUAL(s,"-17")
END_SECTION

START_SECTION((String(long unsigned int i)))
  String s((long unsigned int)(17));
  TEST_EQUAL(s,"17")
END_SECTION

START_SECTION((String(short int i)))
  String s((short int)(-17));
  TEST_EQUAL(s,"-17")
END_SECTION

START_SECTION((String(short unsigned int i)))
  String s((short unsigned int)(17));
  TEST_EQUAL(s,"17")
END_SECTION

START_SECTION((String(float f, bool full_precision = true)))
  String s(17.0123456f);
  TEST_EQUAL(s,"17.012346")
  String s2(17.0123f, false);
  TEST_EQUAL(s2, "17.012")
END_SECTION

START_SECTION((String(float f)))
  float f = 50254.199219;
  double d = f;
  String s(f);
  TEST_EQUAL(s,"5.02542e04")
  String s2(d);
  TEST_EQUAL(s2, "5.025419921875e04")

  // test denormals
  constexpr float denorm = std::numeric_limits<float>::min() / 10;
  //assert(std::fpclassify(denorm) == FP_SUBNORMAL);
  TEST_EQUAL(String(denorm, true), "1.175495e-39")
  TEST_EQUAL(String(denorm, false), "1.175e-39")

  // we need this special 'NaN' since the default 'nan' is not recognized by downstream tools
  // such as any Java-based tool (e.g. KNIME) trying to parse our output files
  constexpr float nan = std::numeric_limits<float>::quiet_NaN();
  assert(std::isnan(nan));
  TEST_EQUAL(String(nan, true), "NaN")
  TEST_EQUAL(String(nan, false), "NaN")

END_SECTION

START_SECTION((String(double d, bool full_precision = true)))
  String s(double(17.012345));
  TEST_EQUAL(s,"17.012345")
  String s2(double(17.012345), false);
  TEST_EQUAL(s2, "17.012")
  // test denormals
  constexpr double denorm = std::numeric_limits<double>::min() / 10;
  assert(std::fpclassify(denorm) == FP_SUBNORMAL);
  TEST_EQUAL(String(denorm, true), "2.225073858507203e-309")
  TEST_EQUAL(String(denorm, false), "2.225e-309")

  // we need this special 'NaN' since the default 'nan' is not recognized by downstream tools 
  // such as any Java-based tool (e.g. KNIME) trying to parse our output files
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  assert(std::isnan(nan));
  TEST_EQUAL(String(nan, true), "NaN")
  TEST_EQUAL(String(nan, false), "NaN")

END_SECTION

START_SECTION((String(long double ld, bool full_precision = true)))
  String s(17.012345L); // suffix L indicates long double
  TEST_EQUAL(s,"17.012345")
  String s2(17.012345L, false); // suffix L indicates long double
  TEST_EQUAL(s2, "17.012")
  // test denormals
  long double denorm = std::numeric_limits<long double>::min() / 10;
  assert(std::fpclassify(denorm) == FP_SUBNORMAL);
  // we cannot test `long double` since it's size is very platform dependent
  // TEST_EQUAL(String(denorm, true), "9.999888671826829e-321")
  //TEST_EQUAL(String(denorm, false), "1.0e-320")
  
  // we need this special 'NaN' since the default 'nan' is not recognized by downstream tools
  // such as any Java-based tool (e.g. KNIME) trying to parse our output files
  constexpr long double nan = std::numeric_limits<long double>::quiet_NaN();
  assert(std::isnan(nan));
  TEST_EQUAL(String(nan, true), "NaN")
  TEST_EQUAL(String(nan, false), "NaN")
END_SECTION

START_SECTION((String(const DataValue& d, bool full_precision = true)))
  TEST_EQUAL(String(DataValue(17.012345)), "17.012345")
  TEST_EQUAL(String(DataValue(17.012345), false), "17.012")
  TEST_EQUAL(String(DataValue(DoubleList({17.012345, 2.0}))), "[17.012345, 2.0]")
  TEST_EQUAL(String(DataValue(DoubleList({17.012345, 2.0})), false), "[17.012, 2.0]")
  TEST_EQUAL(String(DataValue("bla")), "bla")
  TEST_EQUAL(String(DataValue(4711)), "4711")
END_SECTION

START_SECTION((String(long long unsigned int i)))
  String s((long long unsigned int)(12345678));
  TEST_EQUAL(s,"12345678")
END_SECTION

START_SECTION((String(long long signed int i)))
  String s((long long signed int)(-12345678));
  TEST_EQUAL(s,"-12345678")
END_SECTION


START_SECTION((static String numberLength(double d, UInt n)))
  TEST_EQUAL(String::numberLength(12345678.9123,11),"12345678.91")
  TEST_EQUAL(String::numberLength(-12345678.9123,11),"-12345678.9")
  TEST_EQUAL(String::numberLength(12345678.9123,10),"12345678.9")
  TEST_EQUAL(String::numberLength(-12345678.9123,10),"-1234.5e04")
  TEST_EQUAL(String::numberLength(12345678.9123,9),"1234.5e04")
  TEST_EQUAL(String::numberLength(-12345678.9123,9),"-123.4e05")
END_SECTION


START_SECTION((static String number(double d, UInt n)))
  TEST_EQUAL(String::number(123.1234,0),"123")
  TEST_EQUAL(String::number(123.1234,1),"123.1")
  TEST_EQUAL(String::number(123.1234,2),"123.12")
  TEST_EQUAL(String::number(123.1234,3),"123.123")
  TEST_EQUAL(String::number(123.1234,4),"123.1234")
  TEST_EQUAL(String::number(123.1234,5),"123.12340")
  TEST_EQUAL(String::number(0.0,5),"0.00000")
END_SECTION


START_SECTION((template<class InputIterator> String(InputIterator first, InputIterator last)))
  String s("ABCDEFGHIJKLMNOP");
  String::Iterator i = s.begin();
  String::Iterator j = s.end();
  String s2(i,j);
  TEST_EQUAL(s,s2)
  ++i;++i;
  --j;--j;
  s2 = String(i,j);
  TEST_EQUAL(s2,"CDEFGHIJKLMN")

  //test cases where the begin is equal to the end
  i = s.begin();
  j = s.begin();
  s2 = String(i,j);
  TEST_EQUAL(s2,"")
  TEST_EQUAL(s2.size(),0U)

  i = s.end();
  j = s.end();
  s2 = String(i,j);
  TEST_EQUAL(s2,"")
  TEST_EQUAL(s2.size(),0U)


  // Unicode test
  // print(b"T\xc3\xbc\x62ingen".decode("utf8"))
  // print(b"\xff\xfeT\x00\xfc\x00b\x00i\x00n\x00g\x00e\x00n\x00".decode("utf16"))
  // print(b"T\xfc\x62ingen".decode("iso8859"))
  char test_utf8[] =  "T\xc3\xbc\x62ingen";
  char test_utf16[] {'\xff','\xfe','T','\x00','\xfc','\x00','b','\x00','i','\x00','n','\x00','g','\x00','e','\x00','n','\x00'};
  char test_iso8859[] = "T\xfc\x62ingen";

  std::string std_s_utf8(test_utf8, 9);
  std::string std_s_utf16(test_utf16, 18);
  std::string std_s_iso8859(test_iso8859, 8);

  String s_utf8(std_s_utf8.begin(), std_s_utf8.end());
  String s_utf16(std_s_utf16.begin(), std_s_utf16.end());
  String s_iso8859(std_s_iso8859.begin(), std_s_iso8859.end());

  TEST_EQUAL(s_utf16.size(), 18);
  TEST_EQUAL(s_utf8.size(), 9);
  TEST_EQUAL(s_iso8859.size(), 8);

  TEST_EQUAL(s_iso8859[0], 'T')
  TEST_EQUAL(s_iso8859[1], '\xfc') // single byte for u
  TEST_EQUAL(s_iso8859[2], 'b')
  TEST_EQUAL(s_utf8[0], 'T')
  TEST_EQUAL(s_utf8[1], '\xc3')
  TEST_EQUAL(s_utf8[2], '\xbc') // two bytes for u
  TEST_EQUAL(s_utf8[3], 'b')

  // Its very easy to produce nonsense with unicode (e.g. chop string in half
  // in the middle of a character)
  String nonsense(std_s_utf8.begin(), std_s_utf8.begin() + 2);
  TEST_EQUAL(nonsense.size(), 2);
  TEST_EQUAL(nonsense[1], '\xc3')

END_SECTION

String s("ACDEFGHIKLMNPQRSTVWY");
START_SECTION((bool hasPrefix(const String& string) const))
  TEST_EQUAL(s.hasPrefix(""), true);
  TEST_EQUAL(s.hasPrefix("ACDEF"), true);
  TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWY"), true);
  TEST_EQUAL(s.hasPrefix("ABCDEF"), false);
  TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
END_SECTION

START_SECTION((bool hasSuffix(const String& string) const))
  TEST_EQUAL(s.hasSuffix(""), true);
  TEST_EQUAL(s.hasSuffix("TVWY"), true);
  TEST_EQUAL(s.hasSuffix("ACDEFGHIKLMNPQRSTVWY"), true);
  TEST_EQUAL(s.hasSuffix("WXYZ"), false);
  TEST_EQUAL(s.hasSuffix("ACDEFACDEFGHIKLMNPQRSTVWY"), false);
END_SECTION

START_SECTION((bool hasSubstring(const String& string) const))
  TEST_EQUAL(s.hasSubstring(""), true);
  TEST_EQUAL(s.hasSubstring("GHIKLM"), true);
  TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWY"), true);
  TEST_EQUAL(s.hasSubstring("MLKIGH"), false);
  TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
END_SECTION

START_SECTION((bool has(Byte byte) const))
  TEST_EQUAL(s.has('A'), true);
  TEST_EQUAL(s.has('O'), false);
END_SECTION

START_SECTION((String prefix(Int length) const))
  TEST_EQUAL(s.prefix((Int)4), "ACDE");
  TEST_EQUAL(s.prefix((Int)0), "");
  TEST_EXCEPTION(Exception::IndexOverflow, s.prefix(s.size()+1));
  TEST_EXCEPTION(Exception::IndexUnderflow, s.prefix(-1));
END_SECTION

START_SECTION((String suffix(Int length) const))
  TEST_EQUAL(s.suffix((Int)4), "TVWY");
  TEST_EQUAL(s.suffix((Int)0), "");
  TEST_EXCEPTION(Exception::IndexOverflow, s.suffix(s.size()+1));
  TEST_EXCEPTION(Exception::IndexUnderflow, s.suffix(-1));
END_SECTION

START_SECTION((String prefix(SizeType length) const))
  TEST_EQUAL(s.prefix((String::SizeType)4), "ACDE");
  TEST_EQUAL(s.prefix((String::SizeType)0), "");
  TEST_EXCEPTION(Exception::IndexOverflow, s.prefix(s.size()+1));
END_SECTION

START_SECTION((String suffix(SizeType length) const))
  TEST_EQUAL(s.suffix((String::SizeType)4), "TVWY");
  TEST_EQUAL(s.suffix((String::SizeType)0), "");
  TEST_EXCEPTION(Exception::IndexOverflow, s.suffix(s.size()+1));
END_SECTION

START_SECTION((String prefix(char delim) const))
  TEST_EQUAL(s.prefix('F'), "ACDE");
  TEST_EQUAL(s.prefix('A'), "");
  TEST_EXCEPTION(Exception::ElementNotFound, s.prefix('Z'));
END_SECTION

START_SECTION((String suffix(char delim) const))
  TEST_EQUAL(s.suffix('S'), "TVWY");
  TEST_EQUAL(s.suffix('Y'), "");
  TEST_EXCEPTION(Exception::ElementNotFound, s.suffix('Z'));
END_SECTION

START_SECTION((String substr(size_t pos=0, size_t n=npos) const))
  String s("abcdef");
  //std::string functionality
  TEST_EQUAL(s.substr(0,4),"abcd");
  TEST_EQUAL(s.substr(1,1),"b")
  TEST_EQUAL(s.substr(1,3),"bcd")
  TEST_EQUAL(s.substr(0,4),"abcd")
  TEST_EQUAL(s.substr(0,6),"abcdef")
  TEST_EQUAL(s.substr(5,1),"f")
  TEST_EQUAL(s.substr(6,1),"")
  TEST_EQUAL(s.substr(0,7),"abcdef")

  TEST_EQUAL(s.substr(0,String::npos), "abcdef")

  // check with defaults
  TEST_EQUAL(s.substr(0),"abcdef");
  TEST_EQUAL(s.substr(1),"bcdef")
  TEST_EQUAL(s.substr(5),"f")
  TEST_EQUAL(s.substr(6),"")
END_SECTION

START_SECTION((String chop(Size n) const))
  String s("abcdef");
  TEST_EQUAL(s.chop(0), "abcdef")
  TEST_EQUAL(s.chop(1), "abcde")
  TEST_EQUAL(s.chop(2), "abcd")
  TEST_EQUAL(s.chop(3), "abc")
  TEST_EQUAL(s.chop(4), "ab")
  TEST_EQUAL(s.chop(5), "a")
  TEST_EQUAL(s.chop(6), "")
  TEST_EQUAL(s.chop(9), "")

  TEST_EQUAL(s.chop(-1), "")
END_SECTION

START_SECTION((String& reverse()))
  s.reverse();
  TEST_EQUAL(s, "YWVTSRQPNMLKIHGFEDCA");
  s = "";
  s.reverse();
  TEST_EQUAL(s, "");
END_SECTION

START_SECTION((String& trim()))
  String s("\n\r\t test \n\r\t");
  s.trim();
  TEST_EQUAL(s,"test");
  s.trim();
  TEST_EQUAL(s,"test");
  s = "";
  s.trim();
  TEST_EQUAL(s,"");
  s = " t";
  s.trim();
  TEST_EQUAL(s,"t");
  s = "t ";
  s.trim();
  TEST_EQUAL(s,"t");
  s = "\t\r\n ";
  s.trim();
  TEST_EQUAL(s,"");
END_SECTION

START_SECTION((String& quote(char q = '"', QuotingMethod method = ESCAPE)))
  String s;
  s.quote('\'', String::NONE);
  TEST_EQUAL(s, "''");
  s.quote('\'', String::ESCAPE);
  TEST_EQUAL(s, "'\\'\\''");
  s = "ab\"cd\\ef";
  s.quote('"', String::NONE);
  TEST_EQUAL(s, "\"ab\"cd\\ef\"");
  s.quote('"', String::ESCAPE);
  TEST_EQUAL(s, "\"\\\"ab\\\"cd\\\\ef\\\"\"");
  s = "ab\"cd\\ef";
  s.quote('"', String::DOUBLE);
  TEST_EQUAL(s, "\"ab\"\"cd\\ef\"");
END_SECTION

START_SECTION((String& unquote(char q = '"', QuotingMethod method = ESCAPE)))
  String s;
  TEST_EXCEPTION(Exception::ConversionError, s.unquote());
  s = "''";
  s.unquote('\'', String::NONE);
  TEST_EQUAL(s, "");
  s = "'\\'\\''";
  s.unquote('\'', String::ESCAPE);
  TEST_EQUAL(s, "''");
  s = R"("ab"cd\ef")";
  s.unquote('"', String::NONE);
  TEST_EQUAL(s, "ab\"cd\\ef");
  s = R"("\"ab\"cd\\ef\"")";
  s.unquote('"', String::ESCAPE);
  TEST_EQUAL(s, "\"ab\"cd\\ef\"");
  s = R"("ab""cd\ef")";
  s.unquote('"', String::DOUBLE);
  TEST_EQUAL(s, "ab\"cd\\ef");
END_SECTION

START_SECTION((String& simplify()))
  String s("\n\r\t te\tst \n\r\t");
  s.simplify();
  TEST_EQUAL(s," te st ");
  s.simplify();
  TEST_EQUAL(s," te st ");
  s = "";
  s.simplify();
  TEST_EQUAL(s,"");
  s = " t";
  s.simplify();
  TEST_EQUAL(s," t");
  s = "t ";
  s.simplify();
  TEST_EQUAL(s,"t ");
  s = "\t\r\n ";
  s.simplify();
  TEST_EQUAL(s," ");
END_SECTION

START_SECTION((String& fillLeft(char c, UInt size)))
  String s("TEST");
  s.fillLeft('x',4);
  TEST_EQUAL(s,"TEST")
  s.fillLeft('y',5);
  TEST_EQUAL(s,"yTEST")
  s.fillLeft('z',7);
  TEST_EQUAL(s,"zzyTEST")
END_SECTION

START_SECTION((String& fillRight(char c, UInt size)))
  String s("TEST");
  s.fillRight('x',4);
  TEST_EQUAL(s,"TEST")
  s.fillRight('y',5);
  TEST_EQUAL(s,"TESTy")
  s.fillRight('z',7);
  TEST_EQUAL(s,"TESTyzz")
END_SECTION

START_SECTION((Int toInt() const))
  String s;
  s = "123";
  TEST_EQUAL(s.toInt(),123);
  s = "-123";
  TEST_EQUAL(s.toInt(),-123);
  s = "  -123";
  TEST_EQUAL(s.toInt(),-123);
  s = "-123  ";
  TEST_EQUAL(s.toInt(),-123);
  s = "  -123  ";
  TEST_EQUAL(s.toInt(),-123);
  // expect errors:
  s = "524 starts with an int";
  TEST_EXCEPTION(Exception::ConversionError, s.toInt())
  s = "not an int";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, s.toInt(), String("Could not convert string '") + s + "' to an integer value")
  s = "contains an 13135 int";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, s.toInt(), String("Could not convert string '") + s + "' to an integer value")
  s = "ends with an int 525";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, s.toInt(), String("Could not convert string '") + s + "' to an integer value")
END_SECTION

START_SECTION((float toFloat() const))
  String s;
  s = "123.456";
  TEST_REAL_SIMILAR(s.toFloat(),123.456);
  s = "-123.456";
  TEST_REAL_SIMILAR(s.toFloat(),-123.456);
  s = "123.9";
  TEST_REAL_SIMILAR(s.toFloat(),123.9);
  s = "73629.9";
  TEST_EQUAL(String(s.toFloat()),"7.36299e04");
  s = "47218.8";
  TEST_EQUAL(String(s.toFloat()),"4.72188e04");
  s = String("nan");
  TEST_EQUAL(std::isnan(s.toFloat()),true);
  s = "NaN";
  TEST_EQUAL(std::isnan(s.toFloat()),true);
  s = "not a number";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, s.toFloat(), String("Could not convert string '") + s + "' to a float value")
END_SECTION

START_SECTION((double toDouble() const))
  String s;
  s = "123.456";
  TEST_REAL_SIMILAR(s.toDouble(),123.456);
  s = "-123.4567890123";
  TEST_REAL_SIMILAR(s.toDouble(),-123.4567890123);
  s = "123.99999";
  TEST_REAL_SIMILAR(s.toDouble(),123.99999);
  s = "73629.980123";
  TEST_EQUAL(String(s.toDouble()),"7.3629980123e04");
  s = "47218.890000001";
  TEST_EQUAL(String(s.toDouble()),"4.7218890000001e04");
  s = "nan";
  TEST_TRUE(std::isnan(s.toDouble()));
  s = "NaN"; // used in INI files, for Java compatibility
  TEST_TRUE(std::isnan(s.toDouble()));
  s = "not a number";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ConversionError, s.toDouble(), String("Could not convert string '") + s + "' to a double value")
END_SECTION

START_SECTION((static String random(UInt length)))
  String s;
  String s2 = s.random(10);
  TEST_EQUAL(s2.size(),10);
END_SECTION

START_SECTION((bool split(const char splitter, std::vector<String>& substrings, bool quote_protect=false) const))
  String s(";1;2;3;4;5;");
  vector<String> split;
  bool result = s.split(';',split);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),7);
  TEST_EQUAL(split[0],String(""));
  TEST_EQUAL(split[1],String("1"));
  TEST_EQUAL(split[2],String("2"));
  TEST_EQUAL(split[3],String("3"));
  TEST_EQUAL(split[4],String("4"));
  TEST_EQUAL(split[5],String("5"));
  TEST_EQUAL(split[6],String(""));


  s = "1;2;3;4;5";
  result = s.split(';', split);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),5);
  TEST_EQUAL(split[0],String("1"));
  TEST_EQUAL(split[1],String("2"));
  TEST_EQUAL(split[2],String("3"));
  TEST_EQUAL(split[3],String("4"));
  TEST_EQUAL(split[4],String("5"));

  s = "";
  result = s.split(',', split);
  TEST_EQUAL(result, false);
  TEST_EQUAL(split.size(), 0);

  s = ";";
  result = s.split(';', split);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),2);
  TEST_EQUAL(split[0],"");
  TEST_EQUAL(split[1],"");

  result = s.split(',', split);
  TEST_EQUAL(result,false);
  TEST_EQUAL(split.size(),1);

  s = "nodelim";
  result = s.split(';', split);
  TEST_EQUAL(result,false);
  TEST_EQUAL(split.size(),1);

  // testing quoting behaviour
  s=" \"hello\", world, 23.3";
  result = s.split(',', split, true);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),3);
  TEST_EQUAL(split[0],"hello");
  TEST_EQUAL(split[1],"world");
  TEST_EQUAL(split[2],"23.3");

  s=R"( "hello", " donot,splitthis ", "23.4 " )";
  result = s.split(',', split, true);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),3);
  TEST_EQUAL(split[0],"hello");
  TEST_EQUAL(split[1]," donot,splitthis ");
  TEST_EQUAL(split[2],"23.4 ");

  s=R"( "hello", " donot,splitthis ", "23.5 " )";
  result = s.split(',', split, true);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),3);
  TEST_EQUAL(split[0],"hello");
  TEST_EQUAL(split[1]," donot,splitthis ");
  TEST_EQUAL(split[2],"23.5 ");

  s=R"( "hello", " donot,splitthis ", "23.6 " )";
  result = s.split(',', split, true);
  TEST_EQUAL(result,true);
  TEST_EQUAL(split.size(),3);
  TEST_EQUAL(split[0],"hello");
  TEST_EQUAL(split[1]," donot,splitthis ");
  TEST_EQUAL(split[2],"23.6 ");

  s = " \"nodelim \"";
  result = s.split(';', split, true);
  TEST_EQUAL(result,false);
  TEST_EQUAL(split.size(),1);

  // testing invalid quoting...
  s = R"( "first", "seconds"<thisshouldnotbehere>, third)";
  TEST_EXCEPTION(Exception::ConversionError, s.split(',', split, true));
END_SECTION

START_SECTION((bool split(const String& splitter, std::vector<String>& substrings) const))
String s = "abcdebcfghbc";
vector<String> substrings;
bool result = s.split("bc", substrings);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 4);
TEST_EQUAL(substrings[0], "a");
TEST_EQUAL(substrings[1], "de");
TEST_EQUAL(substrings[2], "fgh");
TEST_EQUAL(substrings[3], "");

s = "abcdabcdabcd";
result = s.split("abcd", substrings);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 4);
TEST_EQUAL(substrings[0], "");
TEST_EQUAL(substrings[1], "");
TEST_EQUAL(substrings[2], "");
TEST_EQUAL(substrings[3], "");

result = s.split("xy", substrings);
TEST_EQUAL(result, false);
TEST_EQUAL(substrings.size(), 1);
TEST_EQUAL(substrings[0], s);

result = s.split("", substrings);
TEST_EQUAL(result, true);
TEST_EQUAL(s.size(), substrings.size());
TEST_EQUAL(substrings[0], "a");
TEST_EQUAL(substrings[substrings.size() - 1], "d");

result = String("").split(",", substrings);
TEST_EQUAL(result, false);
TEST_EQUAL(substrings.size(), 0);
END_SECTION

START_SECTION((bool split_quoted(const String& splitter,  std::vector<String>& substrings, char q = '"', QuotingMethod method = ESCAPE) const))
String s = "abcdebcfghbc";
vector<String> substrings;
bool result = s.split_quoted("bc", substrings);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 4);
TEST_EQUAL(substrings[0], "a");
TEST_EQUAL(substrings[1], "de");
TEST_EQUAL(substrings[2], "fgh");
TEST_EQUAL(substrings[3], "");

s = "abcdabcdabcd";
result = s.split_quoted("abcd", substrings);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 4);
TEST_EQUAL(substrings[0], "");
TEST_EQUAL(substrings[1], "");
TEST_EQUAL(substrings[2], "");
TEST_EQUAL(substrings[3], "");

result = s.split_quoted("xy", substrings);
TEST_EQUAL(result, false);
TEST_EQUAL(substrings.size(), 1);
TEST_EQUAL(substrings[0], s);

s = R"("a,b,c","d,\",f","")";
result = s.split_quoted(",", substrings, '"', String::ESCAPE);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 3);
TEST_EQUAL(substrings[0], "\"a,b,c\"");
TEST_EQUAL(substrings[1], "\"d,\\\",f\"");
TEST_EQUAL(substrings[2], "\"\"");

s = R"("a,"b")";
TEST_EXCEPTION(Exception::ConversionError, s.split_quoted(",", substrings, '"',
                                                          String::ESCAPE));
s = R"("ab"___"cd""ef")";
result = s.split_quoted("___", substrings, '"', String::DOUBLE);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 2);
TEST_EQUAL(substrings[0], "\"ab\"");
TEST_EQUAL(substrings[1], "\"cd\"\"ef\"");
END_SECTION

START_SECTION((template<class StringIterator> void concatenate(StringIterator first, StringIterator last, const String& glue = "")))
  vector<String> split;
  String("1;2;3;4;5").split(';',split);
  String s;
  s.concatenate(split.begin(),split.end(),"g");
  TEST_EQUAL(s,"1g2g3g4g5");

  String("1;2;3;4;5").split(';',split);
  s.concatenate(split.begin(),split.end());
  TEST_EQUAL(s,"12345");

  String("").split(';',split);
  s.concatenate(split.begin(),split.end());
  TEST_EQUAL(s,"");

  s.concatenate(split.begin(),split.end(),"_");
  TEST_EQUAL(s,"");
END_SECTION

START_SECTION((String& toUpper()))
  String s;
  s = "test45%#.,";
  s.toUpper();
  TEST_EQUAL(s,"TEST45%#.,");
  s = "";
  s.toUpper();
  TEST_EQUAL(s,"");
END_SECTION

START_SECTION((String& toLower()))
  String s;
  s = "TEST45%#.,";
  s.toLower();
  TEST_EQUAL(s,"test45%#.,");
  s = "";
  s.toLower();
  TEST_EQUAL(s,"");
END_SECTION

START_SECTION((String& firstToUpper()))
  String s;
  s = "test45%#.,";
  s.firstToUpper();
  TEST_EQUAL(s,"Test45%#.,");
  s = " ";
  s.firstToUpper();
  TEST_EQUAL(s," ");
  s = "";
  s.firstToUpper();
  TEST_EQUAL(s,"");
END_SECTION

START_SECTION((String& substitute(char from, char to)))
  String s = "abcdefg";

  s.substitute('a','x');
  TEST_EQUAL(s,"xbcdefg")

  s.substitute('g','y');
  TEST_EQUAL(s,"xbcdefy")

  s.substitute('c','-');
  TEST_EQUAL(s,"xb-defy")

  s = ".....";
  s.substitute('.',',');
  TEST_EQUAL(s,",,,,,")

  s = ".....";
  s.substitute(',','.');
  TEST_EQUAL(s,".....")
END_SECTION

START_SECTION((String& substitute(const String& from, const String& to)))
  //single occurence
  String s = "abcdefg";

  s.substitute("a","x");
  TEST_EQUAL(s,"xbcdefg")

  s.substitute("bcd","y");
  TEST_EQUAL(s,"xyefg")

  s.substitute("fg","");
  TEST_EQUAL(s,"xye")

  s.substitute("e","z!");
  TEST_EQUAL(s,"xyz!")

  s.substitute("u","blblblblbl");
  TEST_EQUAL(s,"xyz!")

  s.substitute("","blblblblbl");
  TEST_EQUAL(s,"xyz!")

  //mutiple occurrences
  s = "abcdefgabcdefgabcdefgab";
  s.substitute("ab","x");
  TEST_EQUAL(s,"xcdefgxcdefgxcdefgx")

  s.substitute("x","");
  TEST_EQUAL(s,"cdefgcdefgcdefg")
END_SECTION

START_SECTION((String& remove(char what)))
  String s = "abcabc";

  s.remove('a');
  TEST_EQUAL(s, "bcbc");

  s.remove('c');
  TEST_EQUAL(s, "bb");

  s.remove('b');
  TEST_EQUAL(s, "");
END_SECTION

START_SECTION((String& ensureLastChar(char end)))
  String s = "/";
  s.ensureLastChar('/');
  TEST_EQUAL(s, "/")

  s.ensureLastChar('\\');
  TEST_EQUAL(s, "/\\")

  s.ensureLastChar('\\');
  TEST_EQUAL(s, "/\\")

  s.ensureLastChar('/');
  TEST_EQUAL(s, "/\\/")
END_SECTION

START_SECTION((String& removeWhitespaces()))
{
  String s;

  s.removeWhitespaces();
  TEST_EQUAL(s,"");
  
  s = " \t \n ";
  s.removeWhitespaces();
  TEST_EQUAL(s, "");

  s = "test";
  s.removeWhitespaces();
  TEST_EQUAL(s,"test");

  s = "\n\r\t test \n\r\t";
  s.removeWhitespaces();
  TEST_EQUAL(s,"test");

  s = "\n\r\t t\ne \ts\rt \n\r\t";
  s.removeWhitespaces();
  TEST_EQUAL(s,"test");

  const std::string test(16 * 1024 + 1, 'x'); // not a multiple of 16, so any SSE code needs to deal with a remainder
  s = test + std::string(100, ' ');
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(s.begin(), s.end(), g);
  s.removeWhitespaces();
  TEST_EQUAL(s, test);
}
END_SECTION

const String fixed("test");

START_SECTION((String operator+ (int i) const))
  TEST_EQUAL(fixed + (int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (unsigned int i) const))
  TEST_EQUAL(fixed + (unsigned int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (short int i) const))
  TEST_EQUAL(fixed + (short int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (short unsigned int i) const))
  TEST_EQUAL(fixed + (short unsigned int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (long int i) const))
  TEST_EQUAL(fixed + (long int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (long unsigned int i) const))
  TEST_EQUAL(fixed + (long unsigned int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (long long unsigned int i) const))
  TEST_EQUAL(fixed + (long long unsigned int)(4), "test4")
END_SECTION

START_SECTION((String operator+ (float f) const))
  TEST_EQUAL(fixed + (float)(4), "test4.0")
END_SECTION

START_SECTION((String operator+ (double d) const))
  TEST_EQUAL(fixed + (double)(4), "test4.0")
END_SECTION

START_SECTION((String operator+(long double ld) const ))
  TEST_EQUAL(fixed + (long double)(4), "test4.0")
END_SECTION

START_SECTION((String operator+ (char c) const))
  TEST_EQUAL(fixed + '4', "test4")
END_SECTION

START_SECTION((String operator+ (const char* s) const))
  TEST_EQUAL(fixed + "bla4", "testbla4")
END_SECTION

START_SECTION((String operator+ (const String& s) const))
  TEST_EQUAL(fixed + String("bla4"), "testbla4")
END_SECTION

START_SECTION((String operator+ (const std::string& s) const))
TEST_EQUAL(fixed.operator+(std::string("bla4")), "testbla4")
END_SECTION


START_SECTION((String& operator+= (int i)))
  String s = "test";
  s += (int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (unsigned int i)))
  String s = "test";
  s += (unsigned int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (short int i)))
  String s = "test";
  s += (short int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (short unsigned int i)))
  String s = "test";
  s += (short unsigned int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (long int i)))
  String s = "test";
  s += (long int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (long unsigned int i)))
  String s = "test";
  s += (long unsigned int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (long long unsigned int i)))
  String s = "test";
  s += (long long unsigned int)(7);
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (float f)))
  String s = "test";
  s += (float)(7.4);
  TEST_EQUAL(s, "test7.4")
END_SECTION

START_SECTION((String& operator+= (double d)))
  String s = "test";
  s += (double)(7.4);
  TEST_EQUAL(s, "test7.4")
END_SECTION

START_SECTION((String& operator+= (long double d)))
{
  String s = "test";
  // long double x = 7.4; // implictly double (not long double!)  =>  7.40000000000000036
  long double x = 7.4L; // explictly long double  =>  7.4
  s += x;
  TEST_EQUAL(s, "test7.4");
}
END_SECTION

START_SECTION((String& operator+= (char c)))
  String s = "test";
  s += 'x';
  TEST_EQUAL(s, "testx")
END_SECTION

START_SECTION((String& operator+= (const char* s)))
  String s = "test";
  s += 7;
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (const String& s)))
  String s = "test";
  s += String("7");
  TEST_EQUAL(s, "test7")
END_SECTION

START_SECTION((String& operator+= (const std::string& s)))
  String s = "test";
  s += std::string("7");
  TEST_EQUAL(s, "test7")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

