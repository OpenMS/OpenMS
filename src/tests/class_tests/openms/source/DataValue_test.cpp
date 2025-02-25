// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ParamValue.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <QString>

#include <sstream>
#include <iostream>

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

  // Just as a sanity check, the size of DataValue should be exactly 16 bytes
  // on a 64 bit system:
  // - 1 byte for the data type
  // - 1 byte for the unit type
  // - 4 bytes for the unit identifier (32bit integer)
  // - 1 byte padding
  // - 8 bytes for the actual data / pointers to data
  std::cout << "\n\n --- Size of DataValue " << sizeof(DataValue) << std::endl;

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
  TEST_TRUE(d == sl)
END_SECTION

START_SECTION((DataValue(const IntList &)))
  IntList il;
  il.push_back(1);
  il.push_back(2);
  DataValue d(il);
  TEST_TRUE(d == il)
END_SECTION

START_SECTION((DataValue(const DoubleList &)))
  DoubleList dl;
  dl.push_back(1.2);
  dl.push_back(22.3333);
  DataValue d(dl);
  DoubleList dldv = d;
  TEST_TRUE(dldv == dl);
END_SECTION

// copy ctor
START_SECTION((DataValue(const DataValue&)))
{
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

  DataValue val;
  {
    DataValue p1((double) 1.23);
    p1.setUnit(8);
    val = DataValue(p1);
  }
  DataValue val2(val);

  TEST_REAL_SIMILAR( (double) val, 1.23)
  TEST_EQUAL( val.getUnit(), 8)
  TEST_REAL_SIMILAR( (double) val2, 1.23)
  TEST_EQUAL( val2.getUnit(), 8)
}
END_SECTION

// move ctor
START_SECTION((DataValue(DataValue&&) noexcept))
{
  // Ensure that DataValue has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(DataValue(std::declval<DataValue&&>())), true)

  DataValue empty;
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
  DataValue copy_of_p1(std::move(p1));
  DataValue copy_of_p3(std::move(p3));
  DataValue copy_of_p4(std::move(p4));
  DataValue copy_of_p5(std::move(p5));
  DataValue copy_of_p6(std::move(p6));
  DataValue copy_of_p7(std::move(p7));
  DataValue copy_of_p8(std::move(p8));
  DataValue copy_of_p9(std::move(p9));
  DataValue copy_of_p10(std::move(p10));
  DataValue copy_of_p11(std::move(p11));
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

  TEST_TRUE(p1 == empty)
  TEST_TRUE(p3 == empty)
  TEST_TRUE(p4 == empty)
  TEST_TRUE(p5 == empty)
  TEST_TRUE(p6 == empty)
  TEST_TRUE(p7 == empty)
  TEST_TRUE(p8 == empty)
  TEST_TRUE(p9 == empty)
  TEST_TRUE(p10 == empty)
  TEST_TRUE(p11 == empty)

  DataValue val;
  {
    DataValue p1((double) 1.23);
    p1.setUnit(8);
    val = DataValue(p1);
  }
  DataValue val2(std::move(val));

  TEST_TRUE(val == empty)
  TEST_REAL_SIMILAR( (double) val2, 1.23)
  TEST_EQUAL( val2.getUnit(), 8)
}
END_SECTION

// assignment operator
START_SECTION((DataValue& operator=(const DataValue&)))
{
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

  DataValue val;
  {
    DataValue p1((double) 1.23);
    p1.setUnit(9);
    val = p1;
  }
  DataValue val2 = val;

  TEST_REAL_SIMILAR( (double) val, 1.23)
  TEST_EQUAL( val.getUnit(), 9)
  TEST_REAL_SIMILAR( (double) val2, 1.23)
  TEST_EQUAL( val2.getUnit(), 9)

}
END_SECTION

// move assignment operator
START_SECTION(( DataValue& operator=(DataValue&&) noexcept ))
{
  // Ensure that DataValue has a no-except move assignment operator.
  TEST_EQUAL(noexcept(declval<DataValue&>() = declval<DataValue &&>()), true)

  DataValue empty;
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
  copy_of_p = std::move(p1);
  TEST_REAL_SIMILAR( (double) copy_of_p, 1.23)
  copy_of_p = std::move(p3);
  TEST_REAL_SIMILAR( (float) copy_of_p, 1.23)
  copy_of_p = std::move(p4);
  TEST_EQUAL( (Int) copy_of_p, -3)
  copy_of_p = std::move(p5);
  TEST_EQUAL( (UInt) copy_of_p, 123)
  copy_of_p = std::move(p6);
  TEST_EQUAL( (std::string) copy_of_p, "test char")
  copy_of_p = std::move(p7);
  TEST_EQUAL( (std::string) copy_of_p, "test string")
  copy_of_p = std::move(p8);
  TEST_EQUAL( copy_of_p == ListUtils::create<String>("test string,string2,last string"), true)
  copy_of_p = std::move(p9);
  TEST_EQUAL( (copy_of_p.isEmpty()), true)
  copy_of_p = std::move(p10);
  TEST_EQUAL(copy_of_p == ListUtils::create<Int>("1,2,3,4,5"), true)
  copy_of_p = std::move(p11);
  TEST_EQUAL(copy_of_p == ListUtils::create<double>("1.2,2.3,3.4"), true)

  TEST_TRUE(p1 == empty)
  TEST_TRUE(p3 == empty)
  TEST_TRUE(p4 == empty)
  TEST_TRUE(p5 == empty)
  TEST_TRUE(p6 == empty)
  TEST_TRUE(p7 == empty)
  TEST_TRUE(p8 == empty)
  TEST_TRUE(p9 == empty)
  TEST_TRUE(p10 == empty)
  TEST_TRUE(p11 == empty)

  DataValue val;
  {
    DataValue p1((double) 1.23);
    p1.setUnit(8);
    val = p1;
  }
  DataValue val2 = std::move(val);

  TEST_TRUE(val == empty)
  TEST_REAL_SIMILAR( (double) val2, 1.23)
  TEST_EQUAL( val2.getUnit(), 8)
}
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
START_SECTION((operator ParamValue() const))
{
  int i = 12;
  double d = 3.41;
  String s = "test";
  IntList i_l = {1, 2};
  DoubleList d_l = {2.71, 3.41};
  StringList s_l = {"test", "list"};
  vector<std::string> std_s_l = {"test", "list"};

  DataValue d_i(i);
  ParamValue p_i = d_i;
  TEST_EQUAL(p_i, ParamValue(i))

  DataValue d_d(d);
  ParamValue p_d = d_d;
  TEST_EQUAL(p_d, ParamValue(d))

  DataValue d_s(s);
  ParamValue p_s = d_s;
  TEST_EQUAL(p_s, ParamValue(s))

  DataValue d_i_l(i_l);
  ParamValue p_i_l = d_i_l;
  TEST_EQUAL(p_i_l, ParamValue(i_l))
  
  DataValue d_d_l(d_l);
  ParamValue p_d_l = d_d_l;
  TEST_EQUAL(p_d_l, ParamValue(d_l))

  DataValue d_s_l(s_l);
  ParamValue p_s_l = d_s_l;
  TEST_EQUAL(p_s_l, ParamValue(std_s_l))
}
END_SECTION

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
  TEST_TRUE(sl_op == d)
END_SECTION

START_SECTION((StringList toStringList() const))
  StringList sl;
  sl << "test string list";
  DataValue d(sl);
  StringList sl_op = d.toStringList();
  TEST_TRUE(sl_op == d)
END_SECTION

START_SECTION((operator IntList() const))
  IntList il;
  il.push_back(1);
  il.push_back(2);
  DataValue d(il);
  IntList il_op = d;
  TEST_TRUE(il_op == il)
  TEST_EXCEPTION(Exception::ConversionError, StringList sl = DataValue("abc,ab");)
END_SECTION

START_SECTION((IntList toIntList() const))
  IntList il;
  il.push_back(1);
  il.push_back(2);
  DataValue d(il);
  IntList il_op = d.toIntList();
  TEST_TRUE(il_op == il)
  TEST_EXCEPTION(Exception::ConversionError, StringList sl = DataValue("abc,ab").toStringList();)
END_SECTION

START_SECTION((operator DoubleList() const))
  DoubleList dl;
  dl.push_back(1.2);
  dl.push_back(22.34455);
  DataValue d(dl);
  DoubleList dl_op = d;
  TEST_TRUE(dl_op == d);
END_SECTION

START_SECTION((DoubleList toDoubleList() const))
  DoubleList dl;
  dl.push_back(1.2);
  dl.push_back(22.34455);
  DataValue d(dl);
  DoubleList dl_op = d.toDoubleList();
  TEST_TRUE(dl_op == d);
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
}
END_SECTION

START_SECTION((operator unsigned long long() const))
{
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
}
END_SECTION

START_SECTION(([EXTRA] friend bool operator==(const DataValue&, const DataValue&)))
{
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

  a = DataValue("hello");
  b = DataValue(std::string("hello"));
  TEST_EQUAL(a==b,true);
  a.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
  TEST_EQUAL(a==b,false);
  b.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
  TEST_EQUAL(a==b,true);
  a.setUnit(1);
  TEST_EQUAL(a==b,false);
  b.setUnit(1);
  TEST_EQUAL(a==b,true);
}
END_SECTION

START_SECTION(([EXTRA] friend bool operator!=(const DataValue&, const DataValue&)))
{
  DataValue a(5.0);
  DataValue b(5.1);
  TEST_EQUAL(a!=b,true);
  a = DataValue((double)15.13001);
  b = DataValue((double)15.13);
  TEST_EQUAL(a!=b,true);

  a = DataValue("hello");
  b = DataValue(std::string("hello"));
  TEST_EQUAL(a!=b,false);
}
END_SECTION

START_SECTION((const char* toChar() const))
  DataValue a;
  TEST_EQUAL(a.toChar() == nullptr, true)
  a = DataValue("hello");
  TEST_STRING_EQUAL(a.toChar(),"hello")
  a = DataValue(5);
  TEST_EXCEPTION(Exception::ConversionError, a.toChar() )
END_SECTION

START_SECTION((String toString(bool full_precision) const))
  DataValue a;
  TEST_EQUAL(a.toString(), "")
  a = DataValue("hello");
  TEST_EQUAL(a.toString(),"hello")
  a = DataValue(5);
  TEST_EQUAL(a.toString(), "5")
  a = DataValue(47.11);
  TEST_EQUAL(a.toString(), "47.109999999999999")
  TEST_EQUAL(a.toString(false), "47.11")
  a = DataValue(-23456.78);
  TEST_EQUAL(a.toString(), "-2.345678e04")
  a = DataValue(ListUtils::create<String>("test string,string2,last string"));
  TEST_EQUAL(a.toString(), "[test string, string2, last string]")
  a = DataValue(ListUtils::create<Int>("1,2,3,4,5"));
  TEST_EQUAL(a.toString(),"[1, 2, 3, 4, 5]")
  a = DataValue(ListUtils::create<double>("1.2,47.11,1.2345678e05"));
  TEST_EQUAL(a.toString(),"[1.2, 47.109999999999999, 1.2345678e05]")
  TEST_EQUAL(a.toString(false), "[1.2, 47.11, 1.235e05]")
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
  TEST_EQUAL(a.toQString().toStdString(), "47.109999999999999")
  a = DataValue(-23456.78);
  TEST_EQUAL(a.toQString().toStdString(), "-2.345678e04")
  a = DataValue(ListUtils::create<String>("test string,string2,last string"));
  TEST_EQUAL(a.toQString().toStdString(), "[test string, string2, last string]")
  a =DataValue(ListUtils::create<Int>("1,2,3"));
  TEST_EQUAL(a.toQString().toStdString(), "[1, 2, 3]")
  a = DataValue(ListUtils::create<double>("1.22,43.23232"));
  TEST_EQUAL(a.toQString().toStdString(),"[1.22, 43.232320000000001]")
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

  a2.setUnit(16); // "millimeters"
  TEST_EQUAL(a2.hasUnit(), true)

}
END_SECTION

START_SECTION((const String& getUnit() const))
{
  DataValue a;
  TEST_EQUAL(a.getUnit(), -1)

  DataValue a1(2.2);
  TEST_EQUAL(a1.getUnit(), -1)

  a1.setUnit(169); // id: UO:0000169 name: parts per million
  TEST_EQUAL(a1.getUnit(), 169)
}
END_SECTION

START_SECTION((void setUnit(const String& unit)))
{
  DataValue a1(2.2);
  TEST_EQUAL(a1.getUnit(), -1)

  a1.setUnit(169); // id: UO:0000169 name: parts per million
  TEST_EQUAL(a1.getUnit(), 169)

  a1.setUnit(9); // id: UO:0000009 name: kilogram
  TEST_EQUAL(a1.getUnit(), 9)

  a1.setUnit(a1.getUnit());
  TEST_EQUAL(a1.hasUnit(), true)
  TEST_EQUAL(a1.getUnit(), 9)
}
END_SECTION

START_SECTION((inline UnitType getUnitType() const))
{
  DataValue a("v");
  TEST_EQUAL(a.getUnitType(), DataValue::UnitType::OTHER)
}
END_SECTION

START_SECTION((inline void setUnitType(const UnitType & u)))
{
  DataValue a("v");
  a.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
  TEST_EQUAL(a.getUnitType(), DataValue::UnitType::MS_ONTOLOGY)
  a.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
  TEST_EQUAL(a.getUnitType(), DataValue::UnitType::UNIT_ONTOLOGY)
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
