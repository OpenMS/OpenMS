// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include <QtCore/QString>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

String* s_ptr = 0;
String* s_nullPointer = 0;
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

	String s4("abcdedfg",15);
	TEST_EQUAL(s4,"abcdedfg")
END_SECTION

START_SECTION((String(const DataValue& d)))
	TEST_EQUAL(String(DataValue(1.4)),"1.4")
	TEST_EQUAL(String(DataValue("bla")),"bla")
	TEST_EQUAL(String(DataValue(4711)),"4711")
END_SECTION

START_SECTION((String(const std::string& s)))
	String s(string("blablabla"));
	TEST_EQUAL(s,"blablabla")
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

START_SECTION((String(float f)))
	String s(float(17.0123));
	TEST_EQUAL(s,"17.0123")
END_SECTION

START_SECTION((String(double d)))
	String s(double(17.012345));
	TEST_EQUAL(s,"17.012345")
END_SECTION

START_SECTION((String(long double ld)))
	String s(17.012345L); // suffix L indicates long double
	TEST_EQUAL(s,"17.012345")
END_SECTION

START_SECTION((String(long long unsigned int i)))
	String s((long long unsigned int)(12345678));
	TEST_EQUAL(s,"12345678")
END_SECTION

START_SECTION((String(long long signed int i)))
	String s((long long signed int)(-12345678));
	TEST_EQUAL(s,"-12345678")
END_SECTION


START_SECTION((static String numberLength(DoubleReal d, UInt n)))
	TEST_EQUAL(String::numberLength(12345678.9123,11),"12345678.91")
	TEST_EQUAL(String::numberLength(-12345678.9123,11),"-12345678.9")
	TEST_EQUAL(String::numberLength(12345678.9123,10),"12345678.9")
	TEST_EQUAL(String::numberLength(-12345678.9123,10),"-1234.5e04")
	TEST_EQUAL(String::numberLength(12345678.9123,9),"1234.5e04")
	TEST_EQUAL(String::numberLength(-12345678.9123,9),"-123.4e05")
END_SECTION


START_SECTION((static String number(DoubleReal d, UInt n)))
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
  s = "\"ab\"cd\\ef\"";
  s.unquote('"', String::NONE);
  TEST_EQUAL(s, "ab\"cd\\ef");
  s = "\"\\\"ab\\\"cd\\\\ef\\\"\"";
  s.unquote('"', String::ESCAPE);
  TEST_EQUAL(s, "\"ab\"cd\\ef\"");
  s = "\"ab\"\"cd\\ef\"";
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
	s = "123.456";
	TEST_EQUAL(s.toInt(),123);
	s = "-123.456";
	TEST_EQUAL(s.toInt(),-123);
	s = "123.9";
	TEST_EQUAL(s.toInt(),123);
	s = "73629.00";
	TEST_EQUAL(s.toInt(),73629);
	s = "73629.50";
	TEST_EQUAL(s.toInt(),73629);
	s = "73629.99";
	TEST_EQUAL(s.toInt(),73629);
END_SECTION

START_SECTION((Real toFloat() const))
	String s;
	s = "123.456";
	TEST_REAL_SIMILAR(s.toFloat(),123.456);
	s = "-123.456";
	TEST_REAL_SIMILAR(s.toFloat(),-123.456);
	s = "123.9";
	TEST_REAL_SIMILAR(s.toFloat(),123.9);
	s = "73629.9";
	TEST_EQUAL(String(s.toFloat()),"73629.9");
	s = "47218.8";
	TEST_EQUAL(String(s.toFloat()),"47218.8");
END_SECTION

START_SECTION((DoubleReal toDouble() const))
	String s;
	s = "123.456";
	TEST_REAL_SIMILAR(s.toDouble(),123.456);
	s = "-123.4567890123";
	TEST_REAL_SIMILAR(s.toDouble(),-123.4567890123);
	s = "123.99999";
	TEST_REAL_SIMILAR(s.toDouble(),123.99999);
	s = "73629.980123";
	TEST_EQUAL(String(s.toDouble()),"73629.980123");
	s = "47218.890000001";
	TEST_EQUAL(String(s.toDouble()),"47218.890000001");
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

	s=" \"hello\", \" donot,splitthis \", \"23.4 \" ";
	result = s.split(',', split, true);
	TEST_EQUAL(result,true);
	TEST_EQUAL(split.size(),3);
	TEST_EQUAL(split[0],"hello");
	TEST_EQUAL(split[1]," donot,splitthis ");
	TEST_EQUAL(split[2],"23.4 ");

	s=" \"hello\", \" donot,splitthis \", \"23.5 \" ";
	result = s.split(',', split, true);
	TEST_EQUAL(result,true);
	TEST_EQUAL(split.size(),3);
	TEST_EQUAL(split[0],"hello");
	TEST_EQUAL(split[1]," donot,splitthis ");
	TEST_EQUAL(split[2],"23.5 ");

	s=" \"hello\", \" donot,splitthis \", \"23.6 \" ";
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
	s = " \"first\", \"seconds\"<thisshouldnotbehere>, third";
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

START_SECTION((bool split_quoted(const String& splitter,	std::vector<String>& substrings, char q = '"', QuotingMethod method = ESCAPE) const))
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

s = "\"a,b,c\",\"d,\\\",f\",\"\"";
result = s.split_quoted(",", substrings, '"', String::ESCAPE);
TEST_EQUAL(result, true);
TEST_EQUAL(substrings.size(), 3);
TEST_EQUAL(substrings[0], "\"a,b,c\"");
TEST_EQUAL(substrings[1], "\"d,\\\",f\"");
TEST_EQUAL(substrings[2], "\"\"");

s = "\"a,\"b\"";
TEST_EXCEPTION(Exception::ConversionError, s.split_quoted(",", substrings, '"',
																													String::ESCAPE));
s = "\"ab\"___\"cd\"\"ef\"";
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

	//mutiple occurences
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
	String s;

	s.removeWhitespaces();
	TEST_EQUAL(s,"");
	
	s = "test";
	s.removeWhitespaces();
	TEST_EQUAL(s,"test");

	s = "\n\r\t test \n\r\t";
	s.removeWhitespaces();
	TEST_EQUAL(s,"test");

	s = "\n\r\t t\ne \ts\rt \n\r\t";
	s.removeWhitespaces();
	TEST_EQUAL(s,"test");
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
	TEST_EQUAL(fixed + (float)(4), "test4")
END_SECTION

START_SECTION((String operator+ (double d) const))
	TEST_EQUAL(fixed + (double)(4), "test4")
END_SECTION

START_SECTION((String operator+(long double ld) const ))
	TEST_EQUAL(fixed + (long double)(4), "test4")
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
	/*
	NOTE (by Clemens): Windows platforms do not really support long double.  See
	<CONCEPT/Types.h>.  I am leaving this code here because it will help to
	clarify how to set writtenDigits on new platforms.
	*/
#if 0
#define ECHO_AND_DO(bla) STATUS(""#bla); bla
#define ECHO_AND_VALUE(bla) STATUS(""#bla << ": " << bla);
  typedef long double longdouble;

	ECHO_AND_VALUE(sizeof(double));
	ECHO_AND_VALUE(std::numeric_limits<double>::digits);
	ECHO_AND_VALUE(std::numeric_limits<double>::digits10);
	ECHO_AND_VALUE(writtenDigits<double>());

	ECHO_AND_DO(std::cout.precision(std::numeric_limits<double>::digits10));
	ECHO_AND_VALUE(typeAsString(7.4) << ": " << 7.4);
	ECHO_AND_VALUE(typeAsString(7.4L) << ": " << 7.4L);

	ECHO_AND_DO(std::cout.precision( writtenDigits<>(double()) ));
	ECHO_AND_VALUE(typeAsString(7.4) << ": " << 7.4);
	ECHO_AND_VALUE(typeAsString(7.4L) << ": " << 7.4L);

	ECHO_AND_VALUE(sizeof(long double));
	ECHO_AND_VALUE(std::numeric_limits<long double>::digits);
	ECHO_AND_VALUE(std::numeric_limits<long double>::digits10);
	ECHO_AND_VALUE( writtenDigits<>( longdouble() ) );

	ECHO_AND_DO(std::cout.precision(std::numeric_limits<long double>::digits10));
	STATUS(typeAsString(7.4) << ": " << 7.4);
	STATUS(typeAsString(7.4L) << ": " << 7.4L);

	ECHO_AND_DO(std::cout.precision(writtenDigits<>( longdouble() )));
	STATUS(typeAsString(7.4) << ": " << 7.4);
	STATUS(typeAsString(7.4L) << ": " << 7.4L);

	const UInt save_prec  = std::cout.precision();
	for ( UInt prec = 10; prec <= 30; ++prec)
	{
		std::cout.precision(prec);
		STATUS("prec: " << prec << "   7.4: " << 7.4 << "   7.4L: " << 7.4L);
	}
	std::cout.precision(save_prec);

#undef ECHO_AND_DO
#undef ECHO_AND_VALUE
#endif /* End of funny stuff by Clemens */

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
