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

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <sstream>

///////////////////////////

START_TEST(DataValue, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

// default ctor
DataValue* dv_ptr = 0;
CHECK((DataValue()))
	dv_ptr = new DataValue;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK((~DataValue()))
	delete dv_ptr;
RESULT

// ctor for all supported types a DataValue object can hold

CHECK((DataValue(double)))
	double x = -3.0;
	DataValue d(x);
	TEST_REAL_EQUAL((double)d, -3.0)
RESULT

CHECK((DataValue(float)))
	float x = 3.0;
	DataValue d(x);
	TEST_REAL_EQUAL((float)d, 3.0)
RESULT

CHECK((DataValue(long)))
	long n = 100000000;
	DataValue d(n);
	TEST_EQUAL((long)d, 100000000)
RESULT

CHECK((DataValue(int)))
	int n = -3000;
	DataValue d(n);
	TEST_EQUAL((int)d, -3000)
RESULT

CHECK((DataValue(short)))
	short n = -3;
	DataValue d(n);
	TEST_EQUAL((short)d, -3)
RESULT

CHECK((DataValue(char*)))
	char* s = "test char";
	DataValue d(s);
	TEST_EQUAL((std::string)d, "test char")
RESULT

CHECK((DataValue(std::string)))
	std::string s = "test string";
	DataValue d(s);
	TEST_EQUAL((std::string)d, "test string")
RESULT

// copy ctor

CHECK((DataValue(const DataValue&)))
	DataValue p1((double) 1.23);
	DataValue p2((short) 3);
	DataValue p3((float) 1.23);
	DataValue p4((int) 3);
	DataValue p5((long) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8;
	DataValue copy_of_p1(p1);
	DataValue copy_of_p2(p2);
	DataValue copy_of_p3(p3);
	DataValue copy_of_p4(p4);
	DataValue copy_of_p5(p5);
	DataValue copy_of_p6(p6);
	DataValue copy_of_p7(p7);
	DataValue copy_of_p8(p8);
	TEST_REAL_EQUAL( (double) copy_of_p1, 1.23)
	TEST_EQUAL( (short) copy_of_p2, 3)
	TEST_REAL_EQUAL( (float) copy_of_p3, 1.23)
	TEST_EQUAL( (int) copy_of_p4, 3)
	TEST_REAL_EQUAL( (long) copy_of_p5, 123)
	TEST_EQUAL( (std::string) copy_of_p6, "test char")
	TEST_EQUAL( (std::string) copy_of_p7, "test string")
	TEST_EQUAL( (copy_of_p8.isEmpty()),true)
RESULT

// assignment operator

CHECK((DataValue& operator = (const DataValue&)))
	DataValue p1((double) 1.23);
	DataValue p2((short) 3);
	DataValue p3((float) 1.23);
	DataValue p4((int) 3);
	DataValue p5((long) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8;
	DataValue copy_of_p;
	copy_of_p = p1;
	TEST_REAL_EQUAL( (double) copy_of_p, 1.23)
	copy_of_p = p2;
	TEST_EQUAL( (short) copy_of_p, 3)
	copy_of_p = p3;
	TEST_REAL_EQUAL( (float) copy_of_p, 1.23)
	copy_of_p = p4;
	TEST_EQUAL( (int) copy_of_p, 3)
	copy_of_p = p5;
	TEST_REAL_EQUAL( (long) copy_of_p, 123)
	copy_of_p = p6;
	TEST_EQUAL( (std::string) copy_of_p, "test char")
	copy_of_p = p7;
	TEST_EQUAL( (std::string) copy_of_p, "test string")
	copy_of_p = p8;
	TEST_EQUAL( (copy_of_p.isEmpty()),true)
RESULT

// Is DataValue object empty?

CHECK((bool isEmpty() const))
	DataValue p1;
	bool res1 =  p1.isEmpty();
	TEST_NOT_EQUAL( res1, 0)
	DataValue p2((float)1.2);
	bool res2 =  p2.isEmpty();
	TEST_EQUAL( res2, 0)
	TEST_REAL_EQUAL( (float) p2, 1.2)
	DataValue p3((short)2);
	bool res3 =  p3.isEmpty();
	TEST_EQUAL( res3, 0)
	TEST_EQUAL( (short) p3, 2)
	DataValue p4("2");
	bool res4 =  p4.isEmpty();
	TEST_EQUAL( res4, 0)
	TEST_EQUAL( (std::string) p4, "2")
RESULT

// conversion operators

CHECK((operator std::string() const throw(Exception::ConversionError)))
	DataValue d((std::string) "test string");
	std::string k = d;
	TEST_EQUAL(k,"test string")
RESULT

CHECK((operator double() const throw(Exception::ConversionError)))
	DataValue d((double) 5.5);
	double k = d;
	TEST_REAL_EQUAL(k,5.5)
RESULT

CHECK((operator float() const throw(Exception::ConversionError)))
	DataValue d((float) 5.45);
	float k = d;
	TEST_REAL_EQUAL(k,5.45)
RESULT

CHECK((operator int() const throw(Exception::ConversionError)))
	DataValue d((int) 55);
	int k = d;
	TEST_EQUAL(k,55)
RESULT

CHECK((operator unsigned int() const throw(Exception::ConversionError)))
	DataValue d((int) 55);
	unsigned int k = d;
	TEST_EQUAL(k,55)
RESULT

CHECK((operator short() const throw(Exception::ConversionError)))
	DataValue d((short) 5);
	short k = d;
	TEST_EQUAL(k,5)
RESULT

CHECK((operator long() const throw(Exception::ConversionError)))
	DataValue d((long) 555);
	long k = d;
	TEST_EQUAL(k,555)
RESULT

CHECK(([EXTRA] (int/unsigned int/float/double/short/long)DataValue("-123.0354")))
	std::string s = "-123.0354";
	DataValue d(s);
	TEST_EQUAL((int)d,-123)
	TEST_EQUAL((unsigned int)d,123)
	TEST_REAL_EQUAL((float)d,float(-123.0354))
	TEST_REAL_EQUAL((double)d,double(-123.0354))
	TEST_EQUAL((short)d,-123)
	TEST_EQUAL((long)d,-123)
RESULT

CHECK(([EXTRA] (int/unsigned int/float/double/short/long)DataValue("-123.0354 Km")))
	std::string s = "-123.0354 Km";
	DataValue d(s);
	TEST_EQUAL((int)d,-123)
	TEST_EQUAL((unsigned int)d,123)
	TEST_REAL_EQUAL((float)d,float(-123.0354))
	TEST_REAL_EQUAL((double)d,double(-123.0354))
	TEST_EQUAL((short)d,-123)
	TEST_EQUAL((long)d,-123)
RESULT

CHECK(([EXTRA] (int/unsigned int/float/double/short/long)DataValue(double(-222.234))))
	DataValue d((double)(-222.234));
	TEST_EQUAL((int)d,-222)
	TEST_EQUAL((unsigned int)d,222)
	TEST_REAL_EQUAL((float)d,float(-222.234))
	TEST_REAL_EQUAL((double)d,double(-222.234))
	TEST_EQUAL((short)d,-222)
	TEST_EQUAL((long)d,-222)
RESULT

CHECK(([EXTRA] (int/unsigned int/float/double/short/long)DataValue(float(-222.234))))
	float val = -222.234;
	DataValue d(val);
	TEST_EQUAL((int)d,-222)
	TEST_EQUAL((unsigned int)d,222)
	TEST_REAL_EQUAL((float)d,(float)(-222.234))
	TEST_REAL_EQUAL((double)d,(double)(-222.234))
	TEST_EQUAL((short)d,-222)
	TEST_EQUAL((long)d,-222)
RESULT

CHECK(([EXTRA] (int/unsigned int/float/double/short/long)DataValue(int(-222))))
	int val = -222;
	DataValue d(val);
	TEST_EQUAL((int)d,-222)
	TEST_EQUAL((unsigned int)d,222)
	TEST_REAL_EQUAL((float)d,(float)(-222.0))
	TEST_REAL_EQUAL((double)d,(double)(-222.0))
	TEST_EQUAL((short)d,-222)
	TEST_EQUAL((long)d,-222)
RESULT

CHECK(([EXTRA] friend bool operator==(const DataValue&, const DataValue&)))
  DataValue a(5.0);
  DataValue b(5.0);
  TEST_EQUAL(a==b,true);
  a = DataValue((double)15.13);
  b = DataValue((double)15.13);
  TEST_EQUAL(a==b,true);
  a = DataValue((float)15.13);
  b = DataValue((float)(17-1.87));
  TEST_EQUAL(a==b,true);
  a = DataValue((int)5);
  b = DataValue((int)5);
  TEST_EQUAL(a==b,true);
  a = DataValue((long)5000);
  b = DataValue((long)5000);
  TEST_EQUAL(a==b,true);
  a = DataValue("hello");
  b = DataValue(std::string("hello"));
  TEST_EQUAL(a==b,true);
  a = DataValue((float)15.13);
  b = DataValue((float)(15.13001));
  TEST_EQUAL(a==b,false);
RESULT

CHECK(([EXTRA] friend bool operator!=(const DataValue&, const DataValue&)))
  DataValue a(5.0);
  DataValue b(5.1);
  TEST_EQUAL(a!=b,true);
  a = DataValue((double)15.13001);
  b = DataValue((double)15.13);
  TEST_EQUAL(a!=b,true);

RESULT

CHECK((const char* toChar() const throw(Exception::ConversionError)))
	DataValue a;
  TEST_EQUAL(a.toChar() == NULL, true)  
  a = DataValue("hello");
  TEST_EQUAL(a.toChar(),std::string("hello"))
	a = DataValue(5);
  TEST_EXCEPTION(Exception::ConversionError, a.toChar() )
RESULT

CHECK((std::string toString() const))
	DataValue a;
  TEST_EQUAL(a.toString(), "")  
  a = DataValue("hello");
  TEST_EQUAL(a.toString(),"hello")
	a = DataValue(5);
  TEST_EQUAL(a.toString(), "5")
  a = DataValue(47.11);
  TEST_EQUAL(a.toString(), "47.11")
  a = DataValue(-23456.78);
  TEST_EQUAL(a.toString(), "-23456.78")
  
RESULT

CHECK(([EXTRA] friend std::ostream& operator<<(std::ostream&, const DataValue&)))
	DataValue a((int)5), b((long)100), c((double)1.111), d((float)1.1), e("hello "), f(std::string("world")), g;
	std::ostringstream os;
  os << a << b << c << d << e << f << g;
  TEST_EQUAL(os.str(),"51001.1111.1hello world")
RESULT

CHECK((DataType valueType() const))
	DataValue a;
	TEST_EQUAL(a.valueType(), DataValue::EMPTYVALUE);

	DataValue a1(1.45);
	TEST_EQUAL(a1.valueType(), DataValue::DOUVALUE);

	DataValue a2(1.34f);
	TEST_EQUAL(a2.valueType(), DataValue::FLOVALUE);

	DataValue a3(123);
	TEST_EQUAL(a3.valueType(), DataValue::INTVALUE);

	DataValue a4("bla");
	TEST_EQUAL(a4.valueType(), DataValue::STRVALUE);

	DataValue a5(short(2));
	TEST_EQUAL(a5.valueType(), DataValue::SHOVALUE);

	DataValue a6(long(2));
	TEST_EQUAL(a6.valueType(), DataValue::LONVALUE);	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
