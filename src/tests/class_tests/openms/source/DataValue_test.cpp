// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <QString>

#include <sstream>

// we ignore the -Wunused-value warning here, since we do not want the compiler
// to report problems like
// DataValue_test.cpp:285:3: warning: expression result unused [-Wunused-value]
//   TEST_EXCEPTION(Exception::ConversionError, (StringList)DataValue("abc,ab"))

#ifdef __clang__
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunused-value"
#endif

///////////////////////////

START_TEST(DataValue, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// default ctor
DataValue* dv_ptr = nullptr;
DataValue* dv_nullPointer = nullptr;
START_SECTION((DataValue()))
	dv_ptr = new DataValue;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~DataValue()))
	delete dv_ptr;
END_SECTION

// ctor for all supported types a DataValue object can hold

START_SECTION((DataValue(long double)))
	long double x = -3.4L;
	DataValue d(x);
	// Note: The implementation uses typedef double (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((double)d, -3.4L)
END_SECTION

START_SECTION((DataValue(double)))
	double x = -3.0;
	DataValue d(x);
	// Note: The implementation uses typedef double (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((double)d, -3.0);
END_SECTION

START_SECTION((DataValue(float)))
	float x = 3.0;
	DataValue d(x);
	// Note: The implementation uses typedef double (as opposed to float, double, long double.)
	TEST_REAL_SIMILAR((double)d, 3.0);
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

START_SECTION((DataValue(long long)))
	long long n = -3000;
	DataValue d(n);
	TEST_EQUAL((long long) d, -3000)
END_SECTION

START_SECTION((DataValue(unsigned long long)))
	unsigned long long n = 3000;
	DataValue d(n);
	TEST_EQUAL((unsigned long long) d, 3000)
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
	TEST_EQUAL(d == sl, true)
END_SECTION

START_SECTION((DataValue(const IntList &)))
	IntList il;
	il.push_back(1);
  il.push_back(2);
	DataValue d(il);
	TEST_EQUAL(d == il, true)
END_SECTION

START_SECTION((DataValue(const DoubleList &)))
	DoubleList dl;
	dl.push_back(1.2);
  dl.push_back(22.3333);
	DataValue d(dl);
  DoubleList dldv = d;
	TEST_EQUAL(dldv == dl, true);
END_SECTION
// copy ctor

START_SECTION((DataValue(const DataValue&)))
	DataValue p1((double) 1.23);
	DataValue p3((float) 1.23);
	DataValue p4((Int) -3);
	DataValue p5((UInt) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8(ListUtils::create<String>("test string,string2,last string"));
	DataValue p9;
	DataValue p10(ListUtils::create<Int>("1,2,3,4,5"));
	DataValue p11(ListUtils::create<double>("1.2,2.3,3.4"));
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
	TEST_REAL_SIMILAR( (double) copy_of_p1, 1.23)
	TEST_REAL_SIMILAR( (float) copy_of_p3, 1.23)
	TEST_EQUAL( (Int) copy_of_p4, -3)
	TEST_EQUAL( (UInt) copy_of_p5, 123)
	TEST_EQUAL( (std::string) copy_of_p6, "test char")
	TEST_EQUAL( (std::string) copy_of_p7, "test string")
	TEST_EQUAL( copy_of_p8 == ListUtils::create<String>("test string,string2,last string"), true)
	TEST_EQUAL( (copy_of_p9.isEmpty()), true)
	TEST_EQUAL( copy_of_p10 == ListUtils::create<Int>("1,2,3,4,5"), true)
	TEST_EQUAL( copy_of_p11 == ListUtils::create<double>("1.2,2.3,3.4"), true)
END_SECTION

// assignment operator

START_SECTION((DataValue& operator=(const DataValue&)))
	DataValue p1((double) 1.23);
	DataValue p3((float) 1.23);
	DataValue p4((Int) -3);
	DataValue p5((UInt) 123);
	DataValue p6("test char");
	DataValue p7(std::string("test string"));
	DataValue p8(ListUtils::create<String>("test string,string2,last string"));
	DataValue p9;
	DataValue p10(ListUtils::create<Int>("1,2,3,4,5"));
	DataValue p11(ListUtils::create<double>("1.2,2.3,3.4"));
	DataValue copy_of_p;
	copy_of_p = p1;
	TEST_REAL_SIMILAR( (double) copy_of_p, 1.23)
	copy_of_p = p3;
	TEST_REAL_SIMILAR( (float) copy_of_p, 1.23)
	copy_of_p = p4;
	TEST_EQUAL( (Int) copy_of_p, -3)
	copy_of_p = p5;
	TEST_EQUAL( (UInt) copy_of_p, 123)
	copy_of_p = p6;
	TEST_EQUAL( (std::string) copy_of_p, "test char")
	copy_of_p = p7;
	TEST_EQUAL( (std::string) copy_of_p, "test string")
	copy_of_p = p8;
	TEST_EQUAL( copy_of_p == ListUtils::create<String>("test string,string2,last string"), true)
	copy_of_p = p9;
	TEST_EQUAL( (copy_of_p.isEmpty()), true)
	copy_of_p = p10;
	TEST_EQUAL(copy_of_p == ListUtils::create<Int>("1,2,3,4,5"), true)
	copy_of_p = p11;
	TEST_EQUAL(copy_of_p == ListUtils::create<double>("1.2,2.3,3.4"), true)
END_SECTION

// Is DataValue object empty?

START_SECTION((bool isEmpty() const))
{
	DataValue p1;
  TEST_EQUAL(p1.isEmpty(), true);
	DataValue p2((float)1.2);
	TEST_EQUAL(p2.isEmpty(), false);
	TEST_REAL_SIMILAR((float) p2, 1.2);
  DataValue p3("");
  TEST_EQUAL(p3.isEmpty(), false); // empty string does not count as empty!
	DataValue p4("2");
	TEST_EQUAL(p4.isEmpty(), false)
  TEST_EQUAL((std::string) p4, "2");
}
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
	TEST_EQUAL(sl_op == d, true)
END_SECTION

START_SECTION((StringList toStringList() const))
	StringList sl;
	sl << "test string list";
	DataValue d(sl);
	StringList sl_op = d.toStringList();
	TEST_EQUAL(sl_op == d, true)
END_SECTION

START_SECTION((operator IntList() const))
	IntList il;
	il.push_back(1);
  il.push_back(2);
	DataValue d(il);
	IntList il_op = d;
	TEST_EQUAL(il_op == il, true)
  TEST_EXCEPTION(Exception::ConversionError, StringList sl = DataValue("abc,ab");)
END_SECTION

START_SECTION((IntList toIntList() const))
	IntList il;
	il.push_back(1);
  il.push_back(2);
	DataValue d(il);
	IntList il_op = d.toIntList();
	TEST_EQUAL(il_op == il, true)
  TEST_EXCEPTION(Exception::ConversionError, StringList sl = DataValue("abc,ab").toStringList();)
END_SECTION

START_SECTION((operator DoubleList() const))
	DoubleList dl;
	dl.push_back(1.2);
  dl.push_back(22.34455);
	DataValue d(dl);
	DoubleList dl_op = d;
	TEST_EQUAL(dl_op == d, true);
END_SECTION

START_SECTION((DoubleList toDoubleList() const))
	DoubleList dl;
	dl.push_back(1.2);
  dl.push_back(22.34455);
	DataValue d(dl);
	DoubleList dl_op = d.toDoubleList();
	TEST_EQUAL(dl_op == d, true);
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
  a = DataValue((double)15.13);
  b = DataValue((double)15.13);
  TEST_EQUAL(a==b,true);
  a = DataValue((float)15.13);
  b = DataValue((float)(17-1.87));
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
  a = DataValue((float)15.13);
  b = DataValue((float)(15.13001));
  TEST_EQUAL(a==b,false);
END_SECTION

START_SECTION(([EXTRA] friend bool operator!=(const DataValue&, const DataValue&)))
  DataValue a(5.0);
  DataValue b(5.1);
  TEST_EQUAL(a!=b,true);
  a = DataValue((double)15.13001);
  b = DataValue((double)15.13);
  TEST_EQUAL(a!=b,true);

END_SECTION

START_SECTION((const char* toChar() const))
	DataValue a;
  TEST_EQUAL(a.toChar() == nullptr, true)
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
  a = DataValue(ListUtils::create<String>("test string,string2,last string"));
  TEST_EQUAL(a.toString(), "[test string, string2, last string]")
  a = DataValue(ListUtils::create<Int>("1,2,3,4,5"));
  TEST_EQUAL(a.toString(),"[1, 2, 3, 4, 5]")
  a= DataValue(ListUtils::create<double>("1.2,23.3333"));
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
  a = DataValue(ListUtils::create<String>("test string,string2,last string"));
  TEST_EQUAL(a.toQString().toStdString(), "[test string, string2, last string]")
  a =DataValue(ListUtils::create<Int>("1,2,3"));
  TEST_EQUAL(a.toQString().toStdString(), "[1, 2, 3]")
  a = DataValue(ListUtils::create<double>("1.22,43.23232"));
  TEST_EQUAL(a.toQString().toStdString(),"[1.22, 43.23232]")
END_SECTION

START_SECTION(([EXTRA] friend std::ostream& operator<<(std::ostream&, const DataValue&)))
	DataValue a((Int)5), b((UInt)100), c((double)1.111), d((double)1.1), e("hello "), f(std::string("world")), g;
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

	DataValue a5(ListUtils::create<String>("test string,string2,last string"));
	TEST_EQUAL(a5.valueType(), DataValue::STRING_LIST)

	DataValue a6(UInt(2));
	TEST_EQUAL(a6.valueType(), DataValue::INT_VALUE);

	DataValue a7(ListUtils::create<Int>("1,2,3"));
	TEST_EQUAL(a7.valueType(),DataValue::INT_LIST)

	DataValue a8(ListUtils::create<double>("1.2,32.4567"));
	TEST_EQUAL(a8.valueType(),DataValue::DOUBLE_LIST);
END_SECTION

START_SECTION((bool hasUnit() const))
{
  DataValue a;
  TEST_EQUAL(a.hasUnit(), false)

	DataValue a1("bla");
  TEST_EQUAL(a1.hasUnit(), false)

	DataValue a2(1.45);
  TEST_EQUAL(a2.hasUnit(), false)

  a2.setUnit("millimeters");
  TEST_EQUAL(a2.hasUnit(), true)
}
END_SECTION

START_SECTION((const String& getUnit() const))
{
  DataValue a;
  TEST_EQUAL(a.getUnit(), "")

	DataValue a1(2.2);
  TEST_EQUAL(a1.getUnit(), "")

  a1.setUnit("ppm");
  TEST_EQUAL(a1.getUnit(), "ppm")
}
END_SECTION

START_SECTION((void setUnit(const String& unit)))
{
	DataValue a1(2.2);
  TEST_EQUAL(a1.getUnit(), "")

  a1.setUnit("ppm");
  TEST_EQUAL(a1.getUnit(), "ppm")

  a1.setUnit("kg");
  TEST_EQUAL(a1.getUnit(), "kg")

}
END_SECTION

START_SECTION((DataValue& operator=(const char*)))
{
  const char * v = "value";
  DataValue a("v");
  a = v;
  TEST_EQUAL((String)a, "value")
}
END_SECTION

START_SECTION((DataValue& operator=(const std::string&)))
{
  std::string v = "value";
  DataValue a("v");
  a = v;
  TEST_EQUAL((String)a, "value")
}
END_SECTION

START_SECTION((DataValue& operator=(const String&)))
{
  String v = "value";
  DataValue a("v");
  a = v;
  TEST_EQUAL((String)a, "value")
}
END_SECTION

START_SECTION((DataValue& operator=(const QString&)))
{
  QString v = "value";
  DataValue a("v");
  a = v;
  TEST_EQUAL((String)a, "value")
}
END_SECTION

START_SECTION((DataValue& operator=(const StringList&)))
{
  StringList v = ListUtils::create<String>("value,value2");
  DataValue a("v");
  a = v;
  StringList sla = a;
  TEST_EQUAL(sla.size(), 2)
  ABORT_IF(sla.size() != 2)
  TEST_EQUAL(sla[0], "value")
  TEST_EQUAL(sla[1], "value2")
}
END_SECTION

START_SECTION((DataValue& operator=(const IntList&)))
{
  IntList v = ListUtils::create<Int>("2,-3");
  DataValue a("v");
  a = v;
  IntList dv = a;
  TEST_EQUAL(dv.size(), 2)
  ABORT_IF(dv.size() != 2)
  TEST_EQUAL(dv[0], 2)
  TEST_EQUAL(dv[1], -3)
}
END_SECTION

START_SECTION((DataValue& operator=(const DoubleList&)))
{
  DoubleList v = ListUtils::create<double>("2.14,-3.45");
  DataValue a("v");
  a = v;
  DoubleList adl = a;
  TEST_EQUAL(adl.size(), 2)
  ABORT_IF(adl.size() != 2)
  TEST_EQUAL(adl[0], 2.14)
  TEST_EQUAL(adl[1], -3.45)
}
END_SECTION

START_SECTION((DataValue& operator=(const long double)))
{
  const long double v = 2.44;
  DataValue a("v");
  a = v;
  TEST_EQUAL((long double)a, 2.44)
}
END_SECTION

START_SECTION((DataValue& operator=(const double)))
{
  const double v = 2.44;
  DataValue a("v");
  a = v;
  TEST_EQUAL((double)a, 2.44)
}
END_SECTION

START_SECTION((DataValue& operator=(const float)))
{
  const float v = 2.44f;
  DataValue a("v");
  a = v;
  TEST_EQUAL((float)a, 2.44f)
}
END_SECTION

START_SECTION((DataValue& operator=(const short int)))
{
  const short int v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((short int)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const unsigned short int)))
{
  const unsigned short int v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((unsigned short int)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const int)))
{
  const int v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((int)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const unsigned)))
{
  const unsigned v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((unsigned)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const long int)))
{
  const long int v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((long int)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const unsigned long)))
{
  const unsigned long v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((unsigned long)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const long long)))
{
  const long long v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((long long)a, 2)
}
END_SECTION

START_SECTION((DataValue& operator=(const unsigned long long)))
{
  const unsigned long long v = 2;
  DataValue a("v");
  a = v;
  TEST_EQUAL((unsigned long long)a, 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#ifdef __clang__
	#pragma clang diagnostic pop
#endif
