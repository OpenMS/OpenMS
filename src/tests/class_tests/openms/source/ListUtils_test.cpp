// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(ListUtils, "$Id$")

START_SECTION((template < typename T, typename E > static bool contains(const std::vector< T > &container, const E &elem)))
{
  // Int
  std::vector<Int> iv;
  iv.push_back(1);
  iv.push_back(2);
  iv.push_back(3);
  iv.push_back(4);

  TEST_EQUAL(ListUtils::contains(iv, 1),true)
  TEST_EQUAL(ListUtils::contains(iv, 2),true)
  TEST_EQUAL(ListUtils::contains(iv, 3),true)
  TEST_EQUAL(ListUtils::contains(iv, 4),true)
  TEST_EQUAL(ListUtils::contains(iv, 5),false)
  TEST_EQUAL(ListUtils::contains(iv, 1011),false)

  // String
  std::vector<String> sv;
  sv.push_back("yes");
  sv.push_back("no");
  TEST_EQUAL(ListUtils::contains(sv, "yes"),true)
  TEST_EQUAL(ListUtils::contains(sv, "no"),true)
  TEST_EQUAL(ListUtils::contains(sv, "jup"),false)
  TEST_EQUAL(ListUtils::contains(sv, ""),false)
  TEST_EQUAL(ListUtils::contains(sv, "noe"),false)
}
END_SECTION

START_SECTION((static bool contains(const std::vector< double > &container, const double &elem, double tolerance=0.00001)))
{
  //
  std::vector<double> dv;
  dv.push_back(1.2);
  dv.push_back(3.4);
  TEST_EQUAL(ListUtils::contains(dv, 1.2),true)
  TEST_EQUAL(ListUtils::contains(dv, 1.21),false)
  TEST_EQUAL(ListUtils::contains(dv, 1.19),false)
  TEST_EQUAL(ListUtils::contains(dv, 1.21,0.02),true)
  TEST_EQUAL(ListUtils::contains(dv, 1.19,0.02),true)
  TEST_EQUAL(ListUtils::contains(dv, 3.4),true)
  TEST_EQUAL(ListUtils::contains(dv, 4.2),false)
  TEST_EQUAL(ListUtils::contains(dv, 2.0),false)
  TEST_EQUAL(ListUtils::contains(dv, 0.0),false)
}
END_SECTION

START_SECTION((template < typename T > static std::vector<T> create(const std::vector< String > &s)))
{
  std::vector<String> iv;
  iv.push_back("1.2");
  iv.push_back("1.56");
  iv.push_back("10.4");

  std::vector<String> sv = ListUtils::create<String>(iv);
  TEST_EQUAL(sv.size(), 3)
  ABORT_IF(sv.size() != 3)
  TEST_EQUAL(sv[0], iv[0])
  TEST_EQUAL(sv[1], iv[1])
  TEST_EQUAL(sv[2], iv[2])

  // create double vector
  std::vector<double> dv = ListUtils::create<double>(iv);
  TEST_EQUAL(dv.size(), 3)
  ABORT_IF(dv.size() != 3)
  TEST_EQUAL(dv[0], 1.2)
  TEST_EQUAL(dv[1], 1.56)
  TEST_EQUAL(dv[2], 10.4)

  iv.push_back("a");
  std::vector<String> sv2 = ListUtils::create<String>(iv);
  TEST_EQUAL(sv2.size(), 4)
  ABORT_IF(sv2.size() != 4)
  TEST_EQUAL(sv2[3], iv[3])

  TEST_EXCEPTION(Exception::ConversionError, ListUtils::create<double>(iv))
}
END_SECTION

START_SECTION((template < typename T > static std::vector<T> create(const String &str, const char splitter= ',')))
{
  std::vector<String> sv = ListUtils::create<String>("yes,no, maybe");
  TEST_EQUAL(sv.size(), 3)
  ABORT_IF(sv.size() != 3)
  TEST_EQUAL(sv[0], "yes")
  TEST_EQUAL(sv[1], "no")
  TEST_EQUAL(sv[2], " maybe")

  std::vector<double> dv = ListUtils::create<double>("1.2,3.5");
  TEST_EQUAL(dv.size(), 2)
  ABORT_IF(dv.size() != 2)
  TEST_EQUAL(dv[0], 1.2)
  TEST_EQUAL(dv[1], 3.5)

  std::vector<Int> iv = ListUtils::create<Int>("1,5");
  TEST_EQUAL(iv.size(),2)
  ABORT_IF(iv.size() != 2)
  TEST_EQUAL(iv[0], 1)
  TEST_EQUAL(iv[1], 5)

  IntList iv2 = ListUtils::create<Int>("2");
  TEST_EQUAL(iv2.size(),1)
  TEST_EQUAL(iv2[0],2)

  IntList iv3 = ListUtils::create<Int>("");
  TEST_EQUAL(iv3.size(),0)

  StringList sl1 = ListUtils::create<String>("test string,string2,last string");
  TEST_EQUAL(sl1.size(),3)
  ABORT_IF(sl1.size() != 3)
  TEST_EQUAL(sl1[0], "test string")
  TEST_EQUAL(sl1[1], "string2")
  TEST_EQUAL(sl1[2], "last string")

  StringList list = ListUtils::create<String>("yes,no");
  TEST_EQUAL(list.size(),2)
  ABORT_IF(list.size() != 2)
  TEST_STRING_EQUAL(list[0],"yes")
  TEST_STRING_EQUAL(list[1],"no")

  StringList list2 = ListUtils::create<String>("no");
  TEST_EQUAL(list2.size(),1)
  ABORT_IF(list2.size() != 1)
  TEST_STRING_EQUAL(list2[0],"no")

  StringList list3 = ListUtils::create<String>("");
  TEST_EQUAL(list3.size(),0)

  StringList sl4 = ListUtils::create<String>("test string#string2#last string", '#');
  TEST_EQUAL(sl4.size(),3)
  ABORT_IF(sl4.size() != 3)
  TEST_EQUAL(sl4[0], "test string")
  TEST_EQUAL(sl4[1], "string2")
  TEST_EQUAL(sl4[2], "last string")
}
END_SECTION

START_SECTION((template < typename T > static String concatenate(const std::vector< T > &container, const String &glue="")))
{
  std::vector<String> list;
  list.push_back("1");
  list.push_back("2");
  list.push_back("3");
  list.push_back("4");
  list.push_back("5");
  TEST_STRING_EQUAL(ListUtils::concatenate(list, "g"),"1g2g3g4g5");
  TEST_STRING_EQUAL(ListUtils::concatenate(list, ""),"12345");

  list.clear();
  TEST_STRING_EQUAL(ListUtils::concatenate(list, "g"),"");
  TEST_STRING_EQUAL(ListUtils::concatenate(list, ""),"");

  //test2 (from StringList)
  std::vector<String> tmp;
  TEST_EQUAL(ListUtils::concatenate(tmp),"")
  tmp.push_back("1\n");
  tmp.push_back("2\n");
  tmp.push_back("3\n");
  TEST_EQUAL(ListUtils::concatenate(tmp),"1\n2\n3\n")
}
END_SECTION

START_SECTION((template <typename T> static Int getIndex(const std::vector<T>& container, const E& elem)))
{
  IntList ints;
  ints.push_back(4);
  ints.push_back(3);
  ints.push_back(1);
  ints.push_back(2);

  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 0), -1);
  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 1), 2);
  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 2), 3);
  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 3), 1);
  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 4), 0);
  TEST_EQUAL(ListUtils::getIndex<Int>(ints, 5), -1);

  StringList strings;
  strings.push_back("four");
  strings.push_back("three");
  strings.push_back("one");
  strings.push_back("two");

  TEST_EQUAL(ListUtils::getIndex<String>(strings, "zero"), -1);
  TEST_EQUAL(ListUtils::getIndex<String>(strings, "one"), 2);
  TEST_EQUAL(ListUtils::getIndex<String>(strings, "two"), 3);
  TEST_EQUAL(ListUtils::getIndex<String>(strings, "three"), 1);
  TEST_EQUAL(ListUtils::getIndex<String>(strings, "four"), 0);
  TEST_EQUAL(ListUtils::getIndex<String>(strings, "five"), -1);
}
END_SECTION

END_TEST
