// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

String* s_ptr = 0;
CHECK((String()))
	s_ptr = new String;
	TEST_NOT_EQUAL(s_ptr, 0)
RESULT

CHECK(([EXTRA] ~String()))
	delete s_ptr;
RESULT

CHECK((String(const char* s, SizeType length)))
	String s("abcdedfg",5);
	TEST_EQUAL(s,"abcde")

	String s2("abcdedfg",0);
	TEST_EQUAL(s2,"")

	String s3("abcdedfg",8);
	TEST_EQUAL(s3,"abcdedfg")

	String s4("abcdedfg",15);
	TEST_EQUAL(s4,"abcdedfg")
RESULT

CHECK((String(const DataValue& d)))
	TEST_EQUAL(String(DataValue(1.4f)),"1.4")
	TEST_EQUAL(String(DataValue("bla")),"bla")
	TEST_EQUAL(String(DataValue(4711)),"4711")
RESULT

CHECK((String(const std::string& s)))
	String s(string("blablabla"));
	TEST_EQUAL(s,"blablabla")
RESULT

CHECK((String(const char* s)))
	String s("blablabla");
	TEST_EQUAL(s,"blablabla")
RESULT

CHECK((String(size_t len, char c)))
	String s(17,'b');
	TEST_EQUAL(s,"bbbbbbbbbbbbbbbbb")
RESULT

CHECK((String(const char c)))
	String s('v');
	TEST_EQUAL(s,"v")
RESULT

CHECK((String(int i)))
	String s(int (-17));
	TEST_EQUAL(s,"-17")
RESULT

CHECK((String(unsigned int i)))
	String s((unsigned int) (17));
	TEST_EQUAL(s,"17")
RESULT

CHECK((String(long int i)))
	String s((long int)(-17));
	TEST_EQUAL(s,"-17")
RESULT

CHECK((String(long unsigned int i)))
	String s((long unsigned int)(17));
	TEST_EQUAL(s,"17")
RESULT

CHECK((String(short int i)))
	String s((short int)(-17));
	TEST_EQUAL(s,"-17")
RESULT

CHECK((String(short unsigned int i)))
	String s((short unsigned int)(17));
	TEST_EQUAL(s,"17")
RESULT

CHECK((String(float f)))
	String s(float(17.0123));
	TEST_EQUAL(s,"17.0123")
RESULT

CHECK((String(double d)))
	String s(double(17.012345));
	TEST_EQUAL(s,"17.012345")
RESULT

CHECK((String(long double d)))
	String s((long double)(17.012345));
	TEST_EQUAL(s,"17.012345")
RESULT

CHECK((String(long long unsigned int i)))
	String s((long long unsigned int)(12345678));
	TEST_EQUAL(s,"12345678")
RESULT

CHECK((String(double d, UnsignedInt size)))
	String s;
	s = String(12345678.9123,11);
	TEST_EQUAL(s,"12345678.91")
	s = String(-12345678.9123,11);
	TEST_EQUAL(s,"-12345678.9")
	
	s = String(12345678.9123,10);
	TEST_EQUAL(s,"12345678.9")
	s = String(-12345678.9123,10);
	TEST_EQUAL(s,"-1234.5e04")
	
	s = String(12345678.9123,9);
	TEST_EQUAL(s,"1234.5e04")
	s = String(-12345678.9123,9);
	TEST_EQUAL(s,"-123.4e05")
RESULT

CHECK((template<class InputIterator> String(InputIterator first, InputIterator last)))
	String s("ABCDEFGHIJKLMNOP");
	String::Iterator i = s.begin();
	String::Iterator j = s.end();
	String s2(i,j);
	TEST_EQUAL(s,s2)
	++i;++i;
	--j;--j;
	s2 = String(i,j);
	TEST_EQUAL("CDEFGHIJKLMN",s2)
RESULT

String s("ACDEFGHIKLMNPQRSTVWY");
CHECK((bool hasPrefix(const String& string) const))
	TEST_EQUAL(s.hasPrefix(""), true);
	TEST_EQUAL(s.hasPrefix("ACDEF"), true);
	TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasPrefix("ABCDEF"), false);
	TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
RESULT

CHECK((bool hasSuffix(const String& string) const))
	TEST_EQUAL(s.hasSuffix(""), true);
	TEST_EQUAL(s.hasSuffix("TVWY"), true);
	TEST_EQUAL(s.hasSuffix("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasSuffix("WXYZ"), false);
	TEST_EQUAL(s.hasSuffix("ACDEFACDEFGHIKLMNPQRSTVWY"), false);
RESULT

CHECK((bool hasSubstring(const String& string) const))
	TEST_EQUAL(s.hasSubstring(""), true);
	TEST_EQUAL(s.hasSubstring("GHIKLM"), true);
	TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasSubstring("MLKIGH"), false);
	TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
RESULT

CHECK((bool has(Byte byte) const))
	TEST_EQUAL(s.has('A'), true);
	TEST_EQUAL(s.has('O'), false);
RESULT

CHECK((String prefix(SignedInt length) const throw(Exception::IndexUnderflow, Exception::IndexOverflow)))
	TEST_EQUAL(s.prefix((SignedInt)4), "ACDE");
	TEST_EQUAL(s.prefix((SignedInt)0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.prefix(s.size()+1));
	TEST_EXCEPTION(Exception::IndexUnderflow, s.prefix(-1));
RESULT

CHECK((String suffix(SignedInt length) const throw(Exception::IndexUnderflow, Exception::IndexOverflow)))
	TEST_EQUAL(s.suffix((SignedInt)4), "TVWY");
	TEST_EQUAL(s.suffix((SignedInt)0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.suffix(s.size()+1));
	TEST_EXCEPTION(Exception::IndexUnderflow, s.suffix(-1));
RESULT

CHECK((String prefix(SizeType length) const throw(Exception::IndexOverflow)))
	TEST_EQUAL(s.prefix((String::SizeType)4), "ACDE");
	TEST_EQUAL(s.prefix((String::SizeType)0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.prefix(s.size()+1));
RESULT

CHECK((String suffix(SizeType length) const throw(Exception::IndexOverflow)))
	TEST_EQUAL(s.suffix((String::SizeType)4), "TVWY");
	TEST_EQUAL(s.suffix((String::SizeType)0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.suffix(s.size()+1));
RESULT

CHECK((String prefix(char delim) const throw (Exception::ElementNotFound)))
	TEST_EQUAL(s.prefix('F'), "ACDE");
	TEST_EQUAL(s.prefix('A'), "");
	TEST_EXCEPTION(Exception::ElementNotFound<char>, s.suffix('Z'));
RESULT

CHECK((String suffix(char delim) const throw (Exception::ElementNotFound)))
	TEST_EQUAL(s.suffix('S'), "TVWY");
	TEST_EQUAL(s.suffix('Y'), "");
	TEST_EXCEPTION(Exception::ElementNotFound<char>, s.suffix('Z'));
RESULT

CHECK((String substr(SignedInt start=0, SignedInt n=NPOS) const))
	String s("abcdef");
	//std::string functionality
	TEST_EQUAL(s.substr(0),"abcdef");
	TEST_EQUAL(s.substr(0,4),"abcd");
	TEST_EQUAL(s.substr(1,1),"b")
	TEST_EQUAL(s.substr(1),"bcdef")
	TEST_EQUAL(s.substr(1,3),"bcd")
	TEST_EQUAL(s.substr(0,4),"abcd")
	TEST_EQUAL(s.substr(0,8),"abcdef")
	//start negative
	TEST_EQUAL(s.substr(-1),"f")
	TEST_EQUAL(s.substr(-2),"ef")
	TEST_EQUAL(s.substr(-3),"def");
	TEST_EQUAL(s.substr(-3,1),"d")
	//n negative
	TEST_EQUAL(s.substr(0,-2),"abcd")
	TEST_EQUAL(s.substr(0,-1),"abcde")
	TEST_EQUAL(s.substr(2,-1),"cde")
	TEST_EQUAL(s.substr(4,-4),"")
	TEST_EQUAL(s.substr(1,-1),"bcde")
	TEST_EQUAL(s.substr(4,-3),"")
	//both negative
	TEST_EQUAL(s.substr(-4,-2),"cd")
	TEST_EQUAL(s.substr(-1,-2),"")
	TEST_EQUAL(s.substr(-3,-2),"d")
	TEST_EQUAL(s.substr(-4,-1),"cde")
	TEST_EQUAL(s.substr(-1,-1),"")
	TEST_EQUAL(s.substr(-3,-1),"de")
RESULT

CHECK((String& reverse()))
	s.reverse();
	TEST_EQUAL(s, "YWVTSRQPNMLKIHGFEDCA");
	s = "";
	s.reverse();
	TEST_EQUAL(s, "");
RESULT

CHECK((String& trim()))
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
RESULT

CHECK((String& fillLeft(char c, UnsignedInt size)))
	String s("TEST");
	s.fillLeft('x',4);
	TEST_EQUAL(s,"TEST")
	s.fillLeft('y',5);
	TEST_EQUAL(s,"yTEST")
	s.fillLeft('z',7);
	TEST_EQUAL(s,"zzyTEST")
RESULT

CHECK((String& fillRight(char c, UnsignedInt size)))
	String s("TEST");
	s.fillRight('x',4);
	TEST_EQUAL(s,"TEST")
	s.fillRight('y',5);
	TEST_EQUAL(s,"TESTy")
	s.fillRight('z',7);
	TEST_EQUAL(s,"TESTyzz")	
RESULT

CHECK((int toInt() const throw(Exception::ConversionError)))
	String s;
	s = "123.456";
	TEST_EQUAL(s.toInt(),123);
	s = "-123.456";
	TEST_EQUAL(s.toInt(),-123);
	s = "123.9";
	TEST_EQUAL(s.toInt(),123);
	s = "73629.00";
	TEST_REAL_EQUAL(s.toInt(),73629);
	s = "73629.50";
	TEST_REAL_EQUAL(s.toInt(),73629);
	s = "73629.99";
	TEST_REAL_EQUAL(s.toInt(),73629);
RESULT

CHECK((float toFloat() const throw(Exception::ConversionError)))
	String s;
	s = "123.456";
	TEST_REAL_EQUAL(s.toFloat(),123.456);
	s = "-123.456";
	TEST_REAL_EQUAL(s.toFloat(),-123.456);
	s = "123.9";
	TEST_REAL_EQUAL(s.toFloat(),123.9);
	s = "73629.98";
	TEST_EQUAL(String(s.toFloat()),"73629.98");
	s = "47218.89";
	TEST_EQUAL(String(s.toFloat()),"47218.89");
RESULT

CHECK((double toDouble() const throw(Exception::ConversionError)))
	String s;
	s = "123.456";
	TEST_REAL_EQUAL(s.toDouble(),123.456);
	s = "-123.456";
	TEST_REAL_EQUAL(s.toDouble(),-123.456);
	s = "123.9";
	TEST_REAL_EQUAL(s.toDouble(),123.9);
	s = "73629.98";
	TEST_EQUAL(String(s.toDouble()),"73629.98");
	s = "47218.89";
	TEST_EQUAL(String(s.toDouble()),"47218.89");
RESULT

CHECK((String random(Size length)))
	String s;
	String s2 = s.random(10);
	TEST_EQUAL(s2.size(),10);
RESULT

CHECK((bool split(char splitter, std::vector<String>& substrings) const))
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

	result = s.split(',', split);
	TEST_EQUAL(result,false);
	TEST_EQUAL(split.size(),0);
RESULT

CHECK((void implode(std::vector<String>::const_iterator first, std::vector<String>::const_iterator last, const std::string& glue = "")))
	vector<String> split;
	String("1;2;3;4;5").split(';',split);
	String s;
	s.implode(split.begin(),split.end(),"g");
	TEST_EQUAL(s,"1g2g3g4g5");
	
	String("1;2;3;4;5").split(';',split);
	s.implode(split.begin(),split.end());
	TEST_EQUAL(s,"12345");

	String("").split(';',split);
	s.implode(split.begin(),split.end());
	TEST_EQUAL(s,"");
	
	s.implode(split.begin(),split.end(),"_");
	TEST_EQUAL(s,"");
RESULT

CHECK((String& toUpper()))
	String s;
	s = "test45%#.,";
	s.toUpper();
	TEST_EQUAL(s,"TEST45%#.,");
	s = "";
	s.toUpper();
	TEST_EQUAL(s,"");
RESULT

CHECK((String& toLower()))
	String s;
	s = "TEST45%#.,";
	s.toLower();
	TEST_EQUAL(s,"test45%#.,");
	s = "";
	s.toLower();
	TEST_EQUAL(s,"");
RESULT

CHECK((String& firstToUpper()))
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
RESULT

CHECK((String& substitute(char from, char to)))
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
RESULT

CHECK((String& remove(char what)))
	String s = "abcabc";

	s.remove('a');
	TEST_EQUAL(s, "bcbc");
	
	s.remove('c');
	TEST_EQUAL(s, "bb");
	
	s.remove('b');
	TEST_EQUAL(s, "");
RESULT

CHECK((String& ensureLastChar(char end)))
	String s = "/";
	s.ensureLastChar('/');
	TEST_EQUAL("/", s)

	s.ensureLastChar('\\');
	TEST_EQUAL("/\\", s)

	s.ensureLastChar('\\');
	TEST_EQUAL("/\\", s)

	s.ensureLastChar('/');
	TEST_EQUAL("/\\/", s)
RESULT

CHECK((String& removeWhitespaces()))
	String s;
	
	s.removeWhitespaces();	
	TEST_EQUAL(s,"");
	
	s = "\n\r\t test \n\r\t";
	s.removeWhitespaces();	
	TEST_EQUAL(s,"test");

	s = "\n\r\t te \n\r\tst \n\r\t";
	s.removeWhitespaces();	
	TEST_EQUAL(s,"test");
RESULT


const String fixed("test");

CHECK((String operator+ (int i) const))
	TEST_EQUAL(fixed + (int)(4), "test4")
RESULT

CHECK((String operator+ (unsigned int i) const))
	TEST_EQUAL(fixed + (unsigned int)(4), "test4")
RESULT

CHECK((String operator+ (short int i) const))
	TEST_EQUAL(fixed + (short int)(4), "test4")
RESULT

CHECK((String operator+ (short unsigned int i) const))
	TEST_EQUAL(fixed + (short unsigned int)(4), "test4")
RESULT

CHECK((String operator+ (long int i) const))
	TEST_EQUAL(fixed + (long int)(4), "test4")
RESULT

CHECK((String operator+ (long unsigned int i) const))
	TEST_EQUAL(fixed + (long unsigned int)(4), "test4")
RESULT

CHECK((String operator+ (long long unsigned int i) const))
	TEST_EQUAL(fixed + (long long unsigned int)(4), "test4")
RESULT

CHECK((String operator+ (float f) const))
	TEST_EQUAL(fixed + (float)(4), "test4")
RESULT

CHECK((String operator+ (double d) const))
	TEST_EQUAL(fixed + (double)(4), "test4")
RESULT

CHECK((String operator+ (long double d) const))
	TEST_EQUAL(fixed + (long double)(4), "test4")
RESULT

CHECK((String operator+ (char c) const))
	TEST_EQUAL(fixed + '4', "test4")
RESULT

CHECK((String operator+ (const char* s) const))
	TEST_EQUAL(fixed + "bla4", "testbla4")
RESULT

CHECK((String operator+ (const String& s) const))
	TEST_EQUAL(fixed + String("bla4"), "testbla4")
RESULT

CHECK((String operator+ (const std::string& s) const))
	TEST_EQUAL(fixed + std::string("bla4"), "testbla4")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
