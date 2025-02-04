// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ParamValue.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

#include <sstream>
#include <iostream>

// we ignore the -Wunused-value warning here, since we do not want the compiler
// to report problems like
// ParamValue_test.cpp:285:3: warning: expression result unused [-Wunused-value]
//   TEST_EXCEPTION(Exception::ConversionError, (StringList)ParamValue("abc,ab"))

#ifdef __clang__
#pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-value"
#endif

///////////////////////////

START_TEST(ParamValue, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

        using namespace OpenMS;
        using namespace std;

// default ctor
        ParamValue* dv_ptr = nullptr;
        ParamValue* dv_nullPointer = nullptr;
        START_SECTION((ParamValue()))

                // Just as a sanity check, the size of ParamValue should be exactly 16 bytes
                // on a 64 bit system:
                // - 1 byte for the data type
                // - 7 byte padding
                // - 8 bytes for the actual data / pointers to data
                std::cout << "\n\n --- Size of ParamValue " << sizeof(ParamValue) << std::endl;

                dv_ptr = new ParamValue;
                TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
        END_SECTION

// destructor
        START_SECTION((virtual ~ParamValue()))
                delete dv_ptr;
        END_SECTION

// ctor for all supported types a ParamValue object can hold

        START_SECTION((ParamValue(long double)))
                long double x = -3.4L;
                ParamValue d(x);
                // Note: The implementation uses typedef double (as opposed to float, double, long double.)
                TEST_REAL_SIMILAR((double)d, -3.4L)
        END_SECTION

        START_SECTION((ParamValue(double)))
                double x = -3.0;
                ParamValue d(x);
                // Note: The implementation uses typedef double (as opposed to float, double, long double.)
                TEST_REAL_SIMILAR((double)d, -3.0);
        END_SECTION

        START_SECTION((ParamValue(float)))
                float x = 3.0;
                ParamValue d(x);
                // Note: The implementation uses typedef double (as opposed to float, double, long double.)
                TEST_REAL_SIMILAR((double)d, 3.0);
        END_SECTION


        START_SECTION((ParamValue(short int)))
                short int n = -3000;
                ParamValue d(n);
                TEST_EQUAL((short int)d, -3000)
        END_SECTION

        START_SECTION((ParamValue(unsigned short int)))
                unsigned short int n = 3000u;
                ParamValue d(n);
                TEST_EQUAL((unsigned short)d, 3000u)
        END_SECTION

        START_SECTION((ParamValue(int)))
                int n = -3000;
                ParamValue d(n);
                TEST_EQUAL((int)d, -3000)
        END_SECTION

        START_SECTION((ParamValue(unsigned)))
                unsigned int n = 3000u;
                ParamValue d(n);
                TEST_EQUAL((unsigned int)d, 3000u)
        END_SECTION

        START_SECTION((ParamValue(long int)))
                long int n = -3000;
                ParamValue d(n);
                TEST_EQUAL((long int)d, -3000)
        END_SECTION

        START_SECTION((ParamValue(unsigned long)))
                unsigned long int n = 3000u;
                ParamValue d(n);
                TEST_EQUAL((unsigned long int)d, 3000u)
        END_SECTION

        START_SECTION((ParamValue(long long)))
                long long n = -3000;
                ParamValue d(n);
                TEST_EQUAL((long long) d, -3000)
        END_SECTION

        START_SECTION((ParamValue(unsigned long long)))
                unsigned long long n = 3000;
                ParamValue d(n);
                TEST_EQUAL((unsigned long long) d, 3000)
        END_SECTION

        START_SECTION((ParamValue(const char*)))
                const char* s = "test char";
                ParamValue d(s);
                TEST_EQUAL((std::string)d, "test char")
        END_SECTION

        START_SECTION((ParamValue(const std::string&)))
                string s = "test string";
                ParamValue d(s);
                TEST_EQUAL(d, "test string")
        END_SECTION

        START_SECTION((ParamValue(const vector<string> &)))
                vector<string> sl = {"test string", "test String 2"};
                ParamValue d(sl);
                TEST_TRUE(d == sl)
        END_SECTION

        START_SECTION((ParamValue(const vector<int> &)))
                vector<int> il = {1, 2};
                ParamValue d(il);
                TEST_TRUE(d == il)
        END_SECTION

        START_SECTION((ParamValue(const vector<double> &)))
                vector<double> dl = {1.2, 22.3333};
                ParamValue d(dl);
                //vector<double> dldv = d;
                TEST_TRUE(d == dl);
        END_SECTION

// copy ctor
        START_SECTION((ParamValue(const ParamValue&)))
                {
                    ParamValue p1((double) 1.23);
                    ParamValue p3((float) 1.23);
                    ParamValue p4((Int) -3);
                    ParamValue p5((UInt) 123);
                    ParamValue p6("test char");
                    ParamValue p7(std::string("test string"));
                    ParamValue p8({"test string","string2","last string"});
                    ParamValue p9;
                    ParamValue p10(vector<int>{1,2,3,4,5});
                    ParamValue p11(vector<double>{1.2,2.3,3.4});
                    ParamValue copy_of_p1(p1);
                    ParamValue copy_of_p3(p3);
                    ParamValue copy_of_p4(p4);
                    ParamValue copy_of_p5(p5);
                    ParamValue copy_of_p6(p6);
                    ParamValue copy_of_p7(p7);
                    ParamValue copy_of_p8(p8);
                    ParamValue copy_of_p9(p9);
                    ParamValue copy_of_p10(p10);
                    ParamValue copy_of_p11(p11);
                    TEST_REAL_SIMILAR( (double) copy_of_p1, 1.23)
                    TEST_REAL_SIMILAR( (float) copy_of_p3, 1.23)
                    TEST_EQUAL( (Int) copy_of_p4, -3)
                    TEST_EQUAL( (UInt) copy_of_p5, 123)
                    TEST_EQUAL( (std::string) copy_of_p6, "test char")
                    TEST_EQUAL( (std::string) copy_of_p7, "test string")
                    TEST_EQUAL( copy_of_p8 == ListUtils::create<string>("test string,string2,last string"), true)
                    TEST_EQUAL( (copy_of_p9.isEmpty()), true)
                    TEST_EQUAL( copy_of_p10 == ListUtils::create<int>("1,2,3,4,5"), true)
                    TEST_EQUAL( copy_of_p11 == ListUtils::create<double>("1.2,2.3,3.4"), true)
                }
        END_SECTION

// move ctor
        START_SECTION((ParamValue(ParamValue&&) noexcept))
                {
                    // Ensure that ParamValue has a no-except move constructor (otherwise
                    // std::vector is inefficient and will copy instead of move).
                    TEST_EQUAL(noexcept(ParamValue(std::declval<ParamValue&&>())), true)

                    ParamValue empty;
                    ParamValue p1((double) 1.23);
                    ParamValue p3((float) 1.23);
                    ParamValue p4((Int) -3);
                    ParamValue p5((UInt) 123);
                    ParamValue p6("test char");
                    ParamValue p7(std::string("test string"));
                    ParamValue p8({"test string","string2","last string"});
                    ParamValue p9;
                    ParamValue p10(vector<int>{1,2,3,4,5});
                    ParamValue p11(vector<double>{1.2,2.3,3.4});
                    ParamValue copy_of_p1(std::move(p1));
                    ParamValue copy_of_p3(std::move(p3));
                    ParamValue copy_of_p4(std::move(p4));
                    ParamValue copy_of_p5(std::move(p5));
                    ParamValue copy_of_p6(std::move(p6));
                    ParamValue copy_of_p7(std::move(p7));
                    ParamValue copy_of_p8(std::move(p8));
                    ParamValue copy_of_p9(std::move(p9));
                    ParamValue copy_of_p10(std::move(p10));
                    ParamValue copy_of_p11(std::move(p11));
                    TEST_REAL_SIMILAR( (double) copy_of_p1, 1.23)
                    TEST_REAL_SIMILAR( (float) copy_of_p3, 1.23)
                    TEST_EQUAL( (Int) copy_of_p4, -3)
                    TEST_EQUAL( (UInt) copy_of_p5, 123)
                    TEST_EQUAL( (std::string) copy_of_p6, "test char")
                    TEST_EQUAL( (std::string) copy_of_p7, "test string")
                    TEST_EQUAL( copy_of_p8 == ListUtils::create<string>("test string,string2,last string"), true)
                    TEST_EQUAL( (copy_of_p9.isEmpty()), true)
                    TEST_EQUAL( copy_of_p10 == ListUtils::create<int>("1,2,3,4,5"), true)
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
                }
        END_SECTION

// assignment operator
        START_SECTION((ParamValue& operator=(const ParamValue&)))
                {
                    ParamValue p1((double) 1.23);
                    ParamValue p3((float) 1.23);
                    ParamValue p4((Int) -3);
                    ParamValue p5((UInt) 123);
                    ParamValue p6("test char");
                    ParamValue p7(std::string("test string"));
                    ParamValue p8({"test string","string2","last string"});
                    ParamValue p9;
                    ParamValue p10(vector<int>{1,2,3,4,5});
                    ParamValue p11(vector<double>{1.2,2.3,3.4});
                    ParamValue copy_of_p;
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
                    TEST_EQUAL( copy_of_p == ListUtils::create<string>("test string,string2,last string"), true)
                    copy_of_p = p9;
                    TEST_EQUAL( (copy_of_p.isEmpty()), true)
                    copy_of_p = p10;
                    TEST_EQUAL(copy_of_p == ListUtils::create<int>("1,2,3,4,5"), true)
                    copy_of_p = p11;
                    TEST_EQUAL(copy_of_p == ListUtils::create<double>("1.2,2.3,3.4"), true)
                }
        END_SECTION

// move assignment operator
        START_SECTION(( ParamValue& operator=(ParamValue&&) noexcept ))
                {
                    // Ensure that ParamValue has a no-except move assignment operator.
                    TEST_EQUAL(noexcept(declval<ParamValue&>() = declval<ParamValue &&>()), true)

                    ParamValue empty;
                    ParamValue p1((double) 1.23);
                    ParamValue p3((float) 1.23);
                    ParamValue p4((Int) -3);
                    ParamValue p5((UInt) 123);
                    ParamValue p6("test char");
                    ParamValue p7(std::string("test string"));
                    ParamValue p8({"test string","string2","last string"});
                    ParamValue p9;
                    ParamValue p10(vector<int>{1,2,3,4,5});
                    ParamValue p11(vector<double>{1.2,2.3,3.4});
                    ParamValue copy_of_p;
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
                    TEST_EQUAL( copy_of_p == ListUtils::create<string>("test string,string2,last string"), true)
                    copy_of_p = std::move(p9);
                    TEST_EQUAL( (copy_of_p.isEmpty()), true)
                    copy_of_p = std::move(p10);
                    TEST_EQUAL(copy_of_p == ListUtils::create<int>("1,2,3,4,5"), true)
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
                }
        END_SECTION

// Is ParamValue object empty?
        START_SECTION((bool isEmpty() const))
                {
                    ParamValue p1;
                    TEST_EQUAL(p1.isEmpty(), true);

                    ParamValue p2((float)1.2);
                    TEST_EQUAL(p2.isEmpty(), false);
                    TEST_REAL_SIMILAR((float) p2, 1.2);

                    ParamValue p3("");
                    TEST_EQUAL(p3.isEmpty(), false); // empty string does not count as empty!

                    ParamValue p4("2");
                    TEST_EQUAL(p4.isEmpty(), false)
                    TEST_EQUAL((std::string) p4, "2");
                }
        END_SECTION

// conversion operators
        START_SECTION((operator std::string() const))
                ParamValue d((std::string) "test string");
                std::string k = d;
                TEST_EQUAL(k,"test string")
        END_SECTION

        START_SECTION((operator vector<string>() const))
                vector<string> sl = {"test string list"};
                ParamValue d(sl);
                vector<string> sl_op = d;
                TEST_TRUE(sl_op == d)
        END_SECTION

        START_SECTION((vector<string> toStringVector() const))
                vector<string> sl = {"test string list"};
                ParamValue d(sl);
                vector<string> sl_op = d.toStringVector();
                TEST_TRUE(sl_op == d)
        END_SECTION

        START_SECTION((operator vector<int>() const))
                vector<int> il = {1, 2};
                ParamValue d(il);
                vector<int> il_op = d;
                TEST_TRUE(il_op == il)
                TEST_EXCEPTION(Exception::ConversionError, vector<string> sl = ParamValue("abc,ab");)
        END_SECTION

        START_SECTION((vector<int> toIntVector() const))
                vector<int> il = {1, 2};
                ParamValue d(il);
                vector<int> il_op = d.toIntVector();
                TEST_TRUE(il_op == il)
                TEST_EXCEPTION(Exception::ConversionError, vector<string> sl = ParamValue("abc,ab").toStringVector();)
        END_SECTION

        START_SECTION((operator vector<double>() const))
                vector<double> dl = {1.2, 22.34455};
                ParamValue d(dl);
                vector<double> dl_op = d;
                TEST_TRUE(dl_op == d);
        END_SECTION

        START_SECTION((DoubleList toDoubleVector() const))
                vector<double> dl = {1.2, 22.34455};
                ParamValue d(dl);
                vector<double> dl_op = d.toDoubleVector();
                TEST_TRUE(dl_op == d);
        END_SECTION

        START_SECTION((operator long double() const))
                ParamValue d(5.4L);
                long double k = d;
                TEST_REAL_SIMILAR(k,5.4L)
        END_SECTION

        START_SECTION((operator double() const))
                ParamValue d(5.4);
                double k = d;
                TEST_REAL_SIMILAR(k,5.4)
        END_SECTION

        START_SECTION((operator float() const))
                ParamValue d(5.4f);
                float k = d;
                TEST_REAL_SIMILAR(k,5.4f)
        END_SECTION

        START_SECTION((operator int() const ))
                ParamValue d((Int) -55);
                int k = d;
                TEST_EQUAL(k,-55)

                TEST_EXCEPTION(Exception::ConversionError, (int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator unsigned int() const ))
                ParamValue d((Int) 55);
                unsigned int k = d;
                TEST_EQUAL(k,55)

                TEST_EXCEPTION(Exception::ConversionError, (unsigned int)ParamValue(-55))
                TEST_EXCEPTION(Exception::ConversionError, (unsigned int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator short int() const))
                ParamValue d((short int) -55);
                short int k = d;
                TEST_EQUAL(k,-55)

                TEST_EXCEPTION(Exception::ConversionError, (short int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator unsigned short int() const))
                ParamValue d((short int) 55);
                unsigned short int k = d;
                TEST_EQUAL(k,55)

                TEST_EXCEPTION(Exception::ConversionError, (unsigned short int)ParamValue(-55))
                TEST_EXCEPTION(Exception::ConversionError, (unsigned short int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator long int() const))
                ParamValue d((long int) -55);
                long int k = d;
                TEST_EQUAL(k,-55)

                TEST_EXCEPTION(Exception::ConversionError, (long int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator unsigned long int() const))
                ParamValue d((long int) 55);
                unsigned long int k = d;
                TEST_EQUAL(k,55)

                TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)ParamValue(-55))
                TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)ParamValue(55.4))
        END_SECTION

        START_SECTION((operator long long() const))
                {
                    {
                        ParamValue d((long long) 55);
                        long long k = d;
                        TEST_EQUAL(k,55)
                    }
                    {
                        ParamValue d((long long) -1);
                        long long k = d;
                        TEST_EQUAL(k,-1)
                    }
                    {
                        ParamValue d((SignedSize) -55);
                        SignedSize k = d;
                        TEST_EQUAL(k,-55)
                    }

                    TEST_EXCEPTION(Exception::ConversionError, (long int)ParamValue(55.4))
                }
        END_SECTION

        START_SECTION((operator unsigned long long() const))
                {
                    {
                        ParamValue d((unsigned long long) 55);
                        unsigned long long k = d;
                        TEST_EQUAL(k,55)
                    }
                    {
                        ParamValue d((Size) 55);
                        Size k = d;
                        TEST_EQUAL(k,55)
                    }

                    TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)ParamValue(-55))
                    TEST_EXCEPTION(Exception::ConversionError, (unsigned long int)ParamValue(55.4))
                }
        END_SECTION

        START_SECTION(([EXTRA] friend bool operator==(const ParamValue&, const ParamValue&)))
                {
                    ParamValue a(5.0);
                    ParamValue b(5.0);
                    TEST_EQUAL(a==b,true);
                    a = ParamValue((double)15.13);
                    b = ParamValue((double)15.13);
                    TEST_EQUAL(a==b,true);
                    a = ParamValue((float)15.13);
                    b = ParamValue((float)(17-1.87));
                    TEST_EQUAL(a==b,true);
                    a = ParamValue((Int)5);
                    b = ParamValue((Int)5);
                    TEST_EQUAL(a==b,true);
                    a = ParamValue((UInt)5000);
                    b = ParamValue((UInt)5000);
                    TEST_EQUAL(a==b,true);
                    a = ParamValue("hello");
                    b = ParamValue(std::string("hello"));
                    TEST_EQUAL(a==b,true);
                    a = ParamValue((float)15.13);
                    b = ParamValue((float)(15.13001));
                    TEST_EQUAL(a==b,false);
                }
        END_SECTION

        START_SECTION(([EXTRA] friend bool operator!=(const ParamValue&, const ParamValue&)))
                {
                    ParamValue a(5.0);
                    ParamValue b(5.1);
                    TEST_EQUAL(a!=b,true);
                    a = ParamValue((double)15.13001);
                    b = ParamValue((double)15.13);
                    TEST_EQUAL(a!=b,true);

                    a = ParamValue("hello");
                    b = ParamValue(std::string("hello"));
                    TEST_EQUAL(a!=b,false);
                }
        END_SECTION

        START_SECTION((const char* toChar() const))
                ParamValue a;
                TEST_EQUAL(a.toChar() == nullptr, true)
                a = ParamValue("hello");
                TEST_STRING_EQUAL(a.toChar(),"hello")
                a = ParamValue(5);
                TEST_EXCEPTION(Exception::ConversionError, a.toChar() )
        END_SECTION

        START_SECTION((String toString(bool full_precision) const))
                ParamValue a;
                TEST_EQUAL(a.toString(), "")
                a = ParamValue("hello");
                TEST_EQUAL(a.toString(),"hello")
                a = ParamValue(5);
                TEST_EQUAL(a.toString(), "5")
                a = ParamValue(47.11);
                TEST_EQUAL(a.toString(), "47.109999999999999")
                TEST_EQUAL(a.toString(false), "47.11")
                a = ParamValue(-23456.78);
                TEST_EQUAL(a.toString(), "-2.345678e04")
                a = ParamValue(ListUtils::create<string>("test string,string2,last string"));
                TEST_EQUAL(a.toString(), "[test string, string2, last string]")
                a = ParamValue(ListUtils::create<int>("1,2,3,4,5"));
                TEST_EQUAL(a.toString(),"[1, 2, 3, 4, 5]")
                a = ParamValue(ListUtils::create<double>("1.2,47.11,1.2345678e05"));
                TEST_EQUAL(a.toString(),"[1.2, 47.109999999999999, 1.2345678e05]")
                TEST_EQUAL(a.toString(false), "[1.2, 47.11, 1.235e05]")
        END_SECTION

        START_SECTION((bool toBool() const))
                //valid cases
                ParamValue a("true");
                TEST_EQUAL(a.toBool(),true)
                a = ParamValue("false");
                TEST_EQUAL(a.toBool(),false)

                //invalid cases
                a = ParamValue();
                TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
                a = ParamValue("bla");
                TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
                a = ParamValue(12);
                TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
                a = ParamValue(34.45);
                TEST_EXCEPTION(Exception::ConversionError, a.toBool() )
        END_SECTION

        START_SECTION(([EXTRA] friend std::ostream& operator<<(std::ostream&, const ParamValue&)))
                ParamValue a((Int)5), b((UInt)100), c((double)1.111), d((double)1.1), e("hello "), f(std::string("world")), g;
                std::ostringstream os;
                os << a << b << c << d << e << f << g;
                TEST_EQUAL(os.str(),"51001.1111.1hello world")
        END_SECTION

        START_SECTION((DataType valueType() const))
                ParamValue a;
                TEST_EQUAL(a.valueType(), ParamValue::EMPTY_VALUE);

                ParamValue a1(1.45);
                TEST_EQUAL(a1.valueType(), ParamValue::DOUBLE_VALUE);

                ParamValue a2(1.34f);
                TEST_EQUAL(a2.valueType(), ParamValue::DOUBLE_VALUE);

                ParamValue a3(123);
                TEST_EQUAL(a3.valueType(), ParamValue::INT_VALUE);

                ParamValue a4("bla");
                TEST_EQUAL(a4.valueType(), ParamValue::STRING_VALUE);

                ParamValue a5(ListUtils::create<string>("test string,string2,last string"));
                TEST_EQUAL(a5.valueType(), ParamValue::STRING_LIST)

                ParamValue a6(UInt(2));
                TEST_EQUAL(a6.valueType(), ParamValue::INT_VALUE);

                ParamValue a7(ListUtils::create<int>("1,2,3"));
                TEST_EQUAL(a7.valueType(),ParamValue::INT_LIST)

                ParamValue a8(ListUtils::create<double>("1.2,32.4567"));
                TEST_EQUAL(a8.valueType(),ParamValue::DOUBLE_LIST);
        END_SECTION

        START_SECTION((ParamValue& operator=(const char*)))
                {
                    const char * v = "value";
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL(a, "value")
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const std::string&)))
                {
                    std::string v = "value";
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL(a, "value")
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const vector<string>&)))
                {
                    vector<string> v = {"value","value2"};
                    ParamValue a("v");
                    a = v;
                    vector<string> sla = a;
                    TEST_EQUAL(sla.size(), 2)
                    ABORT_IF(sla.size() != 2)
                    TEST_EQUAL(sla[0], "value")
                    TEST_EQUAL(sla[1], "value2")
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const vector<int>&)))
                {
                    vector<int> v = {2,-3};
                    ParamValue a("v");
                    a = v;
                    vector<int> dv = a;
                    TEST_EQUAL(dv.size(), 2)
                    ABORT_IF(dv.size() != 2)
                    TEST_EQUAL(dv[0], 2)
                    TEST_EQUAL(dv[1], -3)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const vector<double>&)))
                {
                    vector<double> v = {2.14,-3.45};
                    ParamValue a("v");
                    a = v;
                    vector<double> adl = a;
                    TEST_EQUAL(adl.size(), 2)
                    ABORT_IF(adl.size() != 2)
                    TEST_EQUAL(adl[0], 2.14)
                    TEST_EQUAL(adl[1], -3.45)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const long double)))
                {
                    const long double v = 2.44;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((long double)a, 2.44)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const double)))
                {
                    const double v = 2.44;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((double)a, 2.44)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const float)))
                {
                    const float v = 2.44f;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((float)a, 2.44f)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const short int)))
                {
                    const short int v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((short int)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const unsigned short int)))
                {
                    const unsigned short int v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((unsigned short int)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const int)))
                {
                    const int v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((int)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const unsigned)))
                {
                    const unsigned v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((unsigned)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const long int)))
                {
                    const long int v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((long int)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const unsigned long)))
                {
                    const unsigned long v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((unsigned long)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const long long)))
                {
                    const long long v = 2;
                    ParamValue a("v");
                    a = v;
                    TEST_EQUAL((long long)a, 2)
                }
        END_SECTION

        START_SECTION((ParamValue& operator=(const unsigned long long)))
                {
                    const unsigned long long v = 2;
                    ParamValue a("v");
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
