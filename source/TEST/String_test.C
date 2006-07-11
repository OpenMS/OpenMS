// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id: String_test.C,v 1.7 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

String* s_ptr = 0;
CHECK(String())
	s_ptr = new String;
	TEST_NOT_EQUAL(s_ptr, 0)
RESULT

CHECK(~String())
	delete s_ptr;
RESULT

CHECK(String(const std::string& s))
	String s(string("blablabla"));
	TEST_EQUAL(s,"blablabla")
RESULT

CHECK(String(const char* s))
	String s("blablabla");
	TEST_EQUAL(s,"blablabla")
RESULT

CHECK(String(size_t len, char c))
	String s(17,'b');
	TEST_EQUAL(s,"bbbbbbbbbbbbbbbbb")
RESULT

CHECK(String(char c))
	String s('v');
	TEST_EQUAL(s,"v")
RESULT

CHECK(String(SignedInt i))
	String s(SignedInt(-17));
	TEST_EQUAL(s,"-17")
RESULT

CHECK(String(UnsignedInt i))
	String s(UnsignedInt(17));
	TEST_EQUAL(s,"17")
RESULT

CHECK(String(float f))
	String s(float(17.0123));
	TEST_EQUAL(s,"17.0123")
RESULT

CHECK(String(double d))
	String s(double(17.012345));
	TEST_EQUAL(s,"17.012345")
RESULT

CHECK(String(double d))
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

CHECK(String(InputIterator first,InputIterator last))
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
CHECK( hasPrefix(String) )
	TEST_EQUAL(s.hasPrefix(""), true);
	TEST_EQUAL(s.hasPrefix("ACDEF"), true);
	TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasPrefix("ABCDEF"), false);
	TEST_EQUAL(s.hasPrefix("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
RESULT

CHECK ( hasSuffix(String) )
	TEST_EQUAL(s.hasSuffix(""), true);
	TEST_EQUAL(s.hasSuffix("TVWY"), true);
	TEST_EQUAL(s.hasSuffix("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasSuffix("WXYZ"), false);
	TEST_EQUAL(s.hasSuffix("ACDEFACDEFGHIKLMNPQRSTVWY"), false);
RESULT

CHECK ( hasSubstring(String) )
	TEST_EQUAL(s.hasSubstring(""), true);
	TEST_EQUAL(s.hasSubstring("GHIKLM"), true);
	TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWY"), true);
	TEST_EQUAL(s.hasSubstring("MLKIGH"), false);
	TEST_EQUAL(s.hasSubstring("ACDEFGHIKLMNPQRSTVWYACDEF"), false);
RESULT

CHECK ( has(Byte) )
	TEST_EQUAL(s.has('A'), true);
	TEST_EQUAL(s.has('O'), false);
RESULT

CHECK ( prefix(Index) )
	TEST_EQUAL(s.prefix(4), "ACDE");
	TEST_EQUAL(s.prefix(0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.prefix(s.size()+1));
RESULT

CHECK ( suffix(Index) )
	TEST_EQUAL(s.suffix(4), "TVWY");
	TEST_EQUAL(s.suffix(0), "");
	TEST_EXCEPTION(Exception::IndexOverflow, s.suffix(s.size()+1));
RESULT

CHECK ( prefix(char) )
	TEST_EQUAL(s.prefix('F'), "ACDE");
	TEST_EQUAL(s.prefix('A'), "");
	TEST_EXCEPTION(Exception::ElementNotFound<char>, s.suffix('Z'));
RESULT

CHECK ( suffix(char) )
	TEST_EQUAL(s.suffix('S'), "TVWY");
	TEST_EQUAL(s.suffix('Y'), "");
	TEST_EXCEPTION(Exception::ElementNotFound<char>, s.suffix('Z'));
RESULT

CHECK ( String substr(start, n) )
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

CHECK ( reverse() )
	s.reverse();
	TEST_EQUAL(s, "YWVTSRQPNMLKIHGFEDCA");
	s = "";
	s.reverse();
	TEST_EQUAL(s, "");
RESULT

CHECK ( trim() )
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

CHECK(String& fillLeft(char c, UnsignedInt size))
	String s("TEST");
	s.fillLeft('x',4);
	TEST_EQUAL(s,"TEST")
	s.fillLeft('y',5);
	TEST_EQUAL(s,"yTEST")
	s.fillLeft('z',7);
	TEST_EQUAL(s,"zzyTEST")
RESULT

CHECK(String& fillRight(char c, UnsignedInt size))
	String s("TEST");
	s.fillRight('x',4);
	TEST_EQUAL(s,"TEST")
	s.fillRight('y',5);
	TEST_EQUAL(s,"TESTy")
	s.fillRight('z',7);
	TEST_EQUAL(s,"TESTyzz")	
RESULT

CHECK ( toInt() )
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

CHECK ( toFloat() )
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

CHECK ( toDouble() )
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

CHECK( random(int) )
	String s;
	String s2 = s.random(10);
	TEST_EQUAL(s2.size(),10);
RESULT

CHECK( vector<String> split() )
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

CHECK( void implode(InputIterator, InputIterator, glue))
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

CHECK( void toUpper())
	String s;
	s = "test45%#.,";
	s.toUpper();
	TEST_EQUAL(s,"TEST45%#.,");
RESULT

CHECK( void toLower())
	String s;
	s = "TEST45%#.,";
	s.toLower();
	TEST_EQUAL(s,"test45%#.,");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
