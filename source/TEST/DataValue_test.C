// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
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
using namespace std;

// default ctor
DataValue* dv_ptr = 0;
START_SECTION((DataValue()))
	dv_ptr = new DataValue;
	TEST_NOT_EQUAL(dv_ptr, 0)
END_SECTION

// destructor
START_SECTION((virtual ~DataValue()))
	delete dv_ptr;
END_SECTION

// ctor for all supported types a DataValue object can hold

START_SECTION((DataValue(long double)))
	long double x = -3.4L;
	DataValue d(x);
	// Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((DoubleReal)d, -3.4L)
END_SECTION

START_SECTION((DataValue(double)))
	double x = -3.0;
	DataValue d(x);
	// Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((DoubleReal)d, -3.0);
END_SECTION

START_SECTION((DataValue(float)))
	float x = 3.0;
	DataValue d(x);
	// Note: The implementation uses typedef DoubleReal (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((DoubleReal)d, 3.0);
END_SECTION


START_SECTION((DataValue(short int)))
	short int n = -3000;
	DataValue d(n);
	TEST_EQUAL((short int)d, -3000)
END_SECTION

START_SECTION((DataValue(unsigned short int)))
	unsigned short int n = 3000u;
	DataValue d(n);
	TEST_EQUAL((unsigned short)d, 3000u)
END_SECTION

START_SECTION((DataValue(int)))
  int n = -3000;
  DataValue d(n);
  TEST_EQUAL((int)d, -3000)
END_SECTION

START_SECTION((DataValue(unsigned)))
  unsigned int n = 3000u;
  DataValue d(n);
  TEST_EQUAL((unsigned int)d, 3000u)
END_SECTION

START_SECTION((DataValue(long int)))
  long int n = -3000;
  DataValue d(n);
  TEST_EQUAL((long int)d, -3000)
END_SECTION

START_SECTION((DataValue(unsigned long)))
  unsigned long int n = 3000u;
  DataValue d(n);
  TEST_EQUAL((unsigned long int)d, 3000u)
END_SECTION

START_SECTION((DataValue(const char*)))
	const char* s = "test char";
	DataValue d(s);
	TEST_EQUAL((std::string)d, "test char")
END_SECTION

START_SECTION((DataValue(const std::string&)))
	string s = "test string";
	DataValue d(s);
	TEST_EQUAL((String)d, "test string")
END_SECTION

START_SECTION((DataValue(const QString&)))
	QString s = "test string";
	DataValue d(s);
	TEST_EQUAL((String)d, "test string")
END_SECTION

START_SECTION((DataValue(const String&)))
	String s = "test string";
	DataValue d(s);
	TEST_EQUAL((String)d, "test string")
END_SECTION

START_SECTION((DataValue(const StringList &)))
	StringList sl;
	sl << "test string" << "test String 2";
	DataValue d(sl);
	TEST_EQUAL((StringList)d, sl)
END_SECTION

START_SECTION((DataValue(const IntList &)))
	IntList il;
	il << 1 <<2 ;
	DataValue d(il);
	TEST_EQUAL((IntList)d,il)
END_SECTION

START_SECTION((DataValue(const DoubleList &)))
	DoubleList dl;
	dl << 1.2 << 22.3333;
	DataValue d(dl);
	TEST_EQUAL((DoubleList)d,dl);
END_SECTION
// copy ctor

START_SECTION((DataValue(const DataValue&)))
	DataValue p1((DoubleReal) 1.23);
	DataValue p3((Real) 1.23);
	DataValue p4((Int) -3);
	DataValue p5((UInt) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8(StringList::create("test string,string2,last string"));
	DataValue p9;
	DataValue p10(IntList::create("1,2,3,4,5"));
	DataValue p11(DoubleList::create("1.2,2.3,3.4"));
	DataValue copy_of_p1(p1);
	DataValue copy_of_p3(p3);
	DataValue copy_of_p4(p4);
	DataValue copy_of_p5(p5);
	DataValue copy_of_p6(p6);
	DataValue copy_of_p7(p7);
	DataValue copy_of_p8(p8);
	DataValue copy_of_p9(p9);
	DataValue copy_of_p10(p10);
	DataValue copy_of_p11(p11);
	TEST_REAL_SIMILAR( (DoubleReal) copy_of_p1, 1.23)
	TEST_REAL_SIMILAR( (Real) copy_of_p3, 1.23)
	TEST_EQUAL( (Int) copy_of_p4, -3)
	TEST_EQUAL( (UInt) copy_of_p5, 123)
	TEST_EQUAL( (std::string) copy_of_p6, "test char")
	TEST_EQUAL( (std::string) copy_of_p7, "test string")
	TEST_EQUAL( (StringList) copy_of_p8, StringList::create("test string,string2,last string"))
	TEST_EQUAL( (copy_of_p9.isEmpty()),true)
	TEST_EQUAL((IntList)copy_of_p10,IntList::create("1,2,3,4,5"))
	TEST_EQUAL((DoubleList)copy_of_p11,DoubleList::create("1.2,2.3,3.4"))
END_SECTION

// assignment operator

START_SECTION((DataValue& operator = (const DataValue&)))
	DataValue p1((DoubleReal) 1.23);
	DataValue p3((Real) 1.23);
	DataValue p4((Int) -3);
	DataValue p5((UInt) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8(StringList::create("test string,string2,last string"));
	DataValue p9;
	DataValue p10(IntList::create("1,2,3,4,5"));
	DataValue p11(DoubleList::create("1.2,2.3,3.4"));
	DataValue copy_of_p;
	copy_of_p = p1;
	TEST_REAL_SIMILAR( (DoubleReal) copy_of_p, 1.23)
	copy_of_p = p3;
	TEST_REAL_SIMILAR( (Real) copy_of_p, 1.23)
	copy_of_p = p4;
	TEST_EQUAL( (Int) copy_of_p, -3)
	copy_of_p = p5;
	TEST_EQUAL( (UInt) copy_of_p, 123)
	copy_of_p = p6;
	TEST_EQUAL( (std::string) copy_of_p, "test char")
	copy_of_p = p7;
	TEST_EQUAL( (std::string) copy_of_p, "test string")
	copy_of_p = p8;
	TEST_EQUAL( (StringList) copy_of_p, StringList::create("test string,string2,last string"))
	copy_of_p = p9;
	TEST_EQUAL( (copy_of_p.isEmpty()),true)
	copy_of_p = p10;
	TEST_EQUAL((IntList)copy_of_p,IntList::create("1,2,3,4,5"))
	copy_of_p = p11;
	TEST_EQUAL((DoubleList)copy_of_p,DoubleList::create("1.2,2.3,3.4"))
END_SECTION

// Is DataValue object empty?

START_SECTION((bool isEmpty() const))
	DataValue p1;
	bool res1 =  p1.isEmpty();
	TEST_NOT_EQUAL( res1, false)
	DataValue p2((Real)1.2);
	bool res2 =  p2.isEmpty();
	TEST_EQUAL( res2, false)
	TEST_REAL_SIMILAR( (Real) p2, 1.2)
	DataValue p4("2");
	bool res4 =  p4.isEmpty();
	TEST_EQUAL( res4, false)
	TEST_EQUAL( (std::string) p4, "2")
END_SECTION

// conversion operators

START_SECTION((operator std::string() const))
	DataValue d((std::string) "test string");
	std::string k = d;
	TEST_EQUAL(k,"test string")
END_SECTION

START_SECTION((operator StringList() const))
	StringList sl;
	sl << "test string list";
	DataValue d(sl);
	StringList sl_op = d;
	TEST_EQUAL(sl_op, d)
END_SECTION

START_SECTION((operator IntList() const))
	IntList il;
	il << 1<< 2;
	DataValue d(il);
	IntList il_op = d;
	TEST_EQUAL(il_op,d);


  TEST_EXCEPTION(Exception::ConversionError, (StringList)DataValue("abc,ab"))
END_SECTION

START_SECTION((operator DoubleList() const))
	DoubleList dl;
	dl<< 1.2 <<22.34455;
	DataValue d(dl);
	DoubleList dl_op = d;
	TEST_EQUAL(dl_op,d);
END_SECTION

START_SECTION((operator long double() const))
	DataValue d(5.4L);
	long double k = d;
	TEST_REAL_SIMILAR(k,5.4L)
END_SECTION

START_SECTION((operator double() const))
	DataValue d(5.4);
	double k = d;
	TEST_REAL_SIMILAR(k,5.4)
END_SECTION

START_SECTION((operator float() const))
	DataValue d(5.4f);
	float k = d;
	TEST_REAL_SIMILAR(k,5.4f)
END_SECTION

START_SECTION((operator int() const ))
	DataValue d((Int) -55);
	int k = d;
	TEST_EQUAL(k,-55)

  TEST_EXCEPTION(Exception::ConversionError, (int)DataValue(55.4))
END_SECTION

START_SECTION((operator unsigned int() const ))
	DataValue d((Int) 55);
	unsigned int k = d;
	TEST_EQUAL(k,55)

  TEST_EXCEPTION(Exception::ConversionError, (unsigned int)DataValue(-55))
  TEST_EXCEPTION(Exception::ConversionError, (unsigned int)DataValue(55.4))
END_SECTION

START_SECTION((operator short int() const))
	DataValue d((short int) -55);
	short int k = d;
	TEST_EQUAL(k,-55)

  TEST_EXCEPTION(Exception::ConversionError, (short int)DataValue(55.4))
END_SECTION

START_SECTION((operator unsigned short int() const))
	DataValue d((short int) 55);
	unsigned short int k = d;
	TEST_EQUAL(k,55)

  TEST_EXCEPTION(Exception::ConversionError, (unsigned short int)DataValue(-55))
  TEST_EXCEPTION(Exception::ConversionError, (unsigned short int)DataValue(55.4))
END_SECTION

START_SECTION((operator long int() const))
	DataValue d((long int) -55);
	long int k = d;
	TEST_EQUAL(k,-55)

  TEST_EXCEPTION(Exception::ConversionError, (long int)DataValue(55.4))
END_SECTION

START_SECTION((operator unsigned long int() const))
	DataValue d((long int) 55);
	unsigned long int k = d;
	TEST_EQUAL(k,55)

  TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)DataValue(-55))
  TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)DataValue(55.4))
END_SECTION

START_SECTION((operator long long() const))
	{
	DataValue d((long long) 55);
	long long k = d;
	TEST_EQUAL(k,55)
	}
	{
	DataValue d((long long) -1);
	long long k = d;
	TEST_EQUAL(k,-1)
	}
	{
	DataValue d((SignedSize) -55);
	SignedSize k = d;
	TEST_EQUAL(k,-55)
	}
	
	TEST_EXCEPTION(Exception::ConversionError, (long int)DataValue(55.4))
END_SECTION

START_SECTION((operator unsigned long long() const))
	{
	DataValue d((unsigned long long) 55);
	unsigned long long k = d;
	TEST_EQUAL(k,55)
	}
	{
	DataValue d((Size) 55);
	Size k = d;
	TEST_EQUAL(k,55)
	}

  TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)DataValue(-55))
  TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)DataValue(55.4))
END_SECTION

START_SECTION(([EXTRA] friend bool operator==(const DataValue&, const DataValue&)))
  DataValue a(5.0);
  DataValue b(5.0);
  TEST_EQUAL(a==b,true);
  a = DataValue((DoubleReal)15.13);
  b = DataValue((DoubleReal)15.13);
  TEST_EQUAL(a==b,true);
  a = DataValue((Real)15.13);
  b = DataValue((Real)(17-1.87));
  TEST_EQUAL(a==b,true);
  a = DataValue((Int)5);
  b = DataValue((Int)5);
  TEST_EQUAL(a==b,true);
  a = DataValue((UInt)5000);
  b = DataValue((UInt)5000);
  TEST_EQUAL(a==b,true);
  a = DataValue("hello");
  b = DataValue(std::string("hello"));
  TEST_EQUAL(a==b,true);
  a = DataValue((Real)15.13);
  b = DataValue((Real)(15.13001));
  TEST_EQUAL(a==b,false);
END_SECTION

START_SECTION(([EXTRA] friend bool operator!=(const DataValue&, const DataValue&)))
  DataValue a(5.0);
  DataValue b(5.1);
  TEST_EQUAL(a!=b,true);
  a = DataValue((DoubleReal)15.13001);
  b = DataValue((DoubleReal)15.13);
  TEST_EQUAL(a!=b,true);

END_SECTION

START_SECTION((const char* toChar() const))
	DataValue a;
  TEST_EQUAL(a.toChar() == NULL, true)
  a = DataValue("hello");
  TEST_STRING_EQUAL(a.toChar(),"hello")
	a = DataValue(5);
  TEST_EXCEPTION(Exception::ConversionError, a.toChar() )
END_SECTION

START_SECTION((String toString() const))
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
  a = DataValue(StringList::create("test string,string2,last string"));
  TEST_EQUAL(a.toString(), "[test string, string2, last string]")
  a = DataValue(IntList::create("1,2,3,4,5"));
  TEST_EQUAL(a.toString(),"[1, 2, 3, 4, 5]")
  a= DataValue(DoubleList::create("1.2,23.3333"));
  TEST_EQUAL(a.toString(),"[1.2, 23.3333]")
END_SECTION

START_SECTION((bool toBool() const))
	//valid cases
	DataValue a("true");
	TEST_EQUAL(a.toBool(),true)
	a = DataValue("false");
	TEST_EQUAL(a.toBool(),false)

	//invalid cases
	a = DataValue();
	TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
	a = DataValue("bla");
	TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
	a = DataValue(12);
	TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
	a = DataValue(34.45);
	TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
END_SECTION

START_SECTION((QString toQString() const))
	DataValue a;
  TEST_EQUAL(a.toQString().toStdString(), "")
  a = DataValue("hello");
  TEST_EQUAL(a.toQString().toStdString(),"hello")
	a = DataValue(5);
  TEST_EQUAL(a.toQString().toStdString(), "5")
  a = DataValue(47.11);
  TEST_EQUAL(a.toQString().toStdString(), "47.110000")
  a = DataValue(-23456.78);
  TEST_EQUAL(a.toQString().toStdString(), "-23456.780000")
  a = DataValue(StringList::create("test string,string2,last string"));
  TEST_EQUAL(a.toQString().toStdString(), "[test string, string2, last string]")
  a =DataValue(IntList::create("1,2,3"));
  TEST_EQUAL(a.toQString().toStdString(), "[1, 2, 3]")
  a = DataValue(DoubleList::create("1.22,43.23232"));
  TEST_EQUAL(a.toQString().toStdString(),"[1.22, 43.23232]")
END_SECTION

START_SECTION(([EXTRA] friend std::ostream& operator<<(std::ostream&, const DataValue&)))
	DataValue a((Int)5), b((UInt)100), c((DoubleReal)1.111), d((DoubleReal)1.1), e("hello "), f(std::string("world")), g;
	std::ostringstream os;
  os << a << b << c << d << e << f << g;
  TEST_EQUAL(os.str(),"51001.1111.1hello world")
END_SECTION

START_SECTION((DataType valueType() const))
	DataValue a;
	TEST_EQUAL(a.valueType(), DataValue::EMPTY_VALUE);

	DataValue a1(1.45);
	TEST_EQUAL(a1.valueType(), DataValue::DOUBLE_VALUE);

	DataValue a2(1.34f);
	TEST_EQUAL(a2.valueType(), DataValue::DOUBLE_VALUE);

	DataValue a3(123);
	TEST_EQUAL(a3.valueType(), DataValue::INT_VALUE);

	DataValue a4("bla");
	TEST_EQUAL(a4.valueType(), DataValue::STRING_VALUE);

	DataValue a5(StringList::create("test string,string2,last string"));
	TEST_EQUAL(a5.valueType(), DataValue::STRING_LIST)

	DataValue a6(UInt(2));
	TEST_EQUAL(a6.valueType(), DataValue::INT_VALUE);

	DataValue a7(IntList::create("1,2,3"));
	TEST_EQUAL(a7.valueType(),DataValue::INT_LIST)

	DataValue a8(DoubleList::create("1.2,32.4567"));
	TEST_EQUAL(a8.valueType(),DataValue::DOUBLE_LIST);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

