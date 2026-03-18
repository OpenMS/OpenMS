// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Your Name $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/DataArrays.h>

using namespace OpenMS;
using namespace OpenMS::DataArrays;
using namespace std;

///////////////////////////

START_TEST(DataArrays, "$Id$")

/////////////////////////////////////////////////////////////

FloatDataArray* fda_ptr = nullptr;
FloatDataArray* fda_nullPointer = nullptr;

START_SECTION((FloatDataArray()))
  fda_ptr = new FloatDataArray;
  TEST_NOT_EQUAL(fda_ptr, fda_nullPointer)
  TEST_EQUAL(fda_ptr->size(), 0)
END_SECTION

START_SECTION(([EXTRA] ~FloatDataArray()))
  delete fda_ptr;
END_SECTION

START_SECTION((FloatDataArray copy constructor))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  fda1.setName("test_array");
  FloatDataArray fda2(fda1);
  TEST_EQUAL(fda2.size(), 3)
  TEST_EQUAL(fda2[0], 1.0f)
  TEST_EQUAL(fda2[1], 2.0f)
  TEST_EQUAL(fda2[2], 3.0f)
  TEST_EQUAL(fda2.getName(), "test_array")
END_SECTION

START_SECTION((FloatDataArray operator==))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  FloatDataArray fda4{1.0f, 2.0f};
  
  TEST_TRUE(fda1 == fda2)
  TEST_FALSE(fda1 == fda3)
  TEST_FALSE(fda1 == fda4)
  
  // Test with metadata
  fda1.setName("array1");
  TEST_FALSE(fda1 == fda2)
  fda2.setName("array1");
  TEST_TRUE(fda1 == fda2)
END_SECTION

START_SECTION((FloatDataArray operator!=))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  
  TEST_FALSE(fda1 != fda2)
  TEST_TRUE(fda1 != fda3)
END_SECTION

START_SECTION((FloatDataArray operator<))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  
  TEST_FALSE(fda1 < fda2)
  TEST_TRUE(fda1 < fda3)
  TEST_FALSE(fda3 < fda1)
END_SECTION

START_SECTION((FloatDataArray operator<=))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  
  TEST_TRUE(fda1 <= fda2)
  TEST_TRUE(fda1 <= fda3)
  TEST_FALSE(fda3 <= fda1)
END_SECTION

START_SECTION((FloatDataArray operator>))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  
  TEST_FALSE(fda1 > fda2)
  TEST_FALSE(fda1 > fda3)
  TEST_TRUE(fda3 > fda1)
END_SECTION

START_SECTION((FloatDataArray operator>=))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 4.0f};
  
  TEST_TRUE(fda1 >= fda2)
  TEST_FALSE(fda1 >= fda3)
  TEST_TRUE(fda3 >= fda1)
END_SECTION

// Test special float cases
START_SECTION((FloatDataArray with NaN))
  FloatDataArray fda1{1.0f, std::numeric_limits<float>::quiet_NaN(), 3.0f};
  FloatDataArray fda2{1.0f, std::numeric_limits<float>::quiet_NaN(), 3.0f};
  FloatDataArray fda3{1.0f, 2.0f, 3.0f};
  
  // NaN comparisons should follow IEEE 754 rules
  TEST_FALSE(fda1 == fda2)  // NaN != NaN
  TEST_TRUE(fda1 != fda2)
  
END_SECTION

/////////////////////////////////////////////////////////////
// IntegerDataArray tests

IntegerDataArray* ida_ptr = nullptr;
IntegerDataArray* ida_nullPointer = nullptr;

START_SECTION((IntegerDataArray()))
  ida_ptr = new IntegerDataArray;
  TEST_NOT_EQUAL(ida_ptr, ida_nullPointer)
  TEST_EQUAL(ida_ptr->size(), 0)
END_SECTION

START_SECTION(([EXTRA] ~IntegerDataArray()))
  delete ida_ptr;
END_SECTION

START_SECTION((IntegerDataArray copy constructor))
  IntegerDataArray ida1{1, 2, 3};
  ida1.setName("int_array");
  IntegerDataArray ida2(ida1);
  TEST_EQUAL(ida2.size(), 3)
  TEST_EQUAL(ida2[0], 1)
  TEST_EQUAL(ida2[1], 2)
  TEST_EQUAL(ida2[2], 3)
  TEST_EQUAL(ida2.getName(), "int_array")
END_SECTION

START_SECTION((IntegerDataArray operator==))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_TRUE(ida1 == ida2)
  TEST_FALSE(ida1 == ida3)
END_SECTION

START_SECTION((IntegerDataArray operator!=))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_FALSE(ida1 != ida2)
  TEST_TRUE(ida1 != ida3)
END_SECTION

START_SECTION((IntegerDataArray operator<))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_FALSE(ida1 < ida2)
  TEST_TRUE(ida1 < ida3)
  TEST_FALSE(ida3 < ida1)
END_SECTION

START_SECTION((IntegerDataArray operator<=))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_TRUE(ida1 <= ida2)
  TEST_TRUE(ida1 <= ida3)
  TEST_FALSE(ida3 <= ida1)
END_SECTION

START_SECTION((IntegerDataArray operator>))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_FALSE(ida1 > ida2)
  TEST_FALSE(ida1 > ida3)
  TEST_TRUE(ida3 > ida1)
END_SECTION

START_SECTION((IntegerDataArray operator>=))
  IntegerDataArray ida1{1, 2, 3};
  IntegerDataArray ida2{1, 2, 3};
  IntegerDataArray ida3{1, 2, 4};
  
  TEST_TRUE(ida1 >= ida2)
  TEST_FALSE(ida1 >= ida3)
  TEST_TRUE(ida3 >= ida1)
END_SECTION

/////////////////////////////////////////////////////////////
// StringDataArray tests

StringDataArray* sda_ptr = nullptr;
StringDataArray* sda_nullPointer = nullptr;

START_SECTION((StringDataArray()))
  sda_ptr = new StringDataArray;
  TEST_NOT_EQUAL(sda_ptr, sda_nullPointer)
  TEST_EQUAL(sda_ptr->size(), 0)
END_SECTION

START_SECTION(([EXTRA] ~StringDataArray()))
  delete sda_ptr;
END_SECTION

START_SECTION((StringDataArray copy constructor))
  StringDataArray sda1{"apple", "banana", "cherry"};
  sda1.setName("string_array");
  StringDataArray sda2(sda1);
  TEST_EQUAL(sda2.size(), 3)
  TEST_EQUAL(sda2[0], "apple")
  TEST_EQUAL(sda2[1], "banana")
  TEST_EQUAL(sda2[2], "cherry")
  TEST_EQUAL(sda2.getName(), "string_array")
END_SECTION

START_SECTION((StringDataArray operator==))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_TRUE(sda1 == sda2)
  TEST_FALSE(sda1 == sda3)
END_SECTION

START_SECTION((StringDataArray operator!=))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_FALSE(sda1 != sda2)
  TEST_TRUE(sda1 != sda3)
END_SECTION

START_SECTION((StringDataArray operator<))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_FALSE(sda1 < sda2)
  TEST_TRUE(sda1 < sda3)
  TEST_FALSE(sda3 < sda1)
END_SECTION

START_SECTION((StringDataArray operator<=))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_TRUE(sda1 <= sda2)
  TEST_TRUE(sda1 <= sda3)
  TEST_FALSE(sda3 <= sda1)
END_SECTION

START_SECTION((StringDataArray operator>))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_FALSE(sda1 > sda2)
  TEST_FALSE(sda1 > sda3)
  TEST_TRUE(sda3 > sda1)
END_SECTION

START_SECTION((StringDataArray operator>=))
  StringDataArray sda1{"apple", "banana", "cherry"};
  StringDataArray sda2{"apple", "banana", "cherry"};
  StringDataArray sda3{"apple", "banana", "date"};
  
  TEST_TRUE(sda1 >= sda2)
  TEST_FALSE(sda1 >= sda3)
  TEST_TRUE(sda3 >= sda1)
END_SECTION

// Test std::vector functionality inheritance
START_SECTION((DataArrays std::vector functionality))
  FloatDataArray fda;
  fda.push_back(1.0f);
  fda.push_back(2.0f);
  fda.push_back(3.0f);
  
  TEST_EQUAL(fda.size(), 3)
  TEST_EQUAL(fda.front(), 1.0f)
  TEST_EQUAL(fda.back(), 3.0f)
  
  fda.clear();
  TEST_EQUAL(fda.size(), 0)
  TEST_TRUE(fda.empty())
  
  // Test reserve and capacity
  fda.reserve(100);
  TEST_TRUE(fda.capacity() >= 100)
  
  // Test iteration
  IntegerDataArray ida{1, 2, 3, 4, 5};
  int sum = 0;
  for (const auto& val : ida)
  {
    sum += val;
  }
  TEST_EQUAL(sum, 15)
END_SECTION

// Test MetaInfoDescription functionality inheritance
START_SECTION((DataArrays MetaInfoDescription functionality))
  FloatDataArray fda;
  fda.setName("my_float_array");
  TEST_EQUAL(fda.getName(), "my_float_array")
  
  IntegerDataArray ida;
  ida.setName("my_int_array");
  TEST_EQUAL(ida.getName(), "my_int_array")
  
  StringDataArray sda;
  sda.setName("my_string_array");
  TEST_EQUAL(sda.getName(), "my_string_array")
  
  // Test that metadata affects comparison
  FloatDataArray fda1{1.0f, 2.0f};
  FloatDataArray fda2{1.0f, 2.0f};
  
  TEST_TRUE(fda1 == fda2)
  
  fda1.setName("array1");
  TEST_FALSE(fda1 == fda2)
  
  fda2.setName("array1");
  TEST_TRUE(fda1 == fda2)
  
  // Test MetaInfo functionality
  fda1.setMetaValue("label", "test_label");
  TEST_FALSE(fda1 == fda2)
  
  fda2.setMetaValue("label", "test_label");
  TEST_TRUE(fda1 == fda2)
END_SECTION

// Test mixed comparisons with different metadata
START_SECTION((DataArrays comparison with mixed metadata))
  FloatDataArray fda1{1.0f, 2.0f, 3.0f};
  FloatDataArray fda2{1.0f, 2.0f, 3.0f};
  
  // Same data, different names
  fda1.setName("A");
  fda2.setName("B");
  
  // Arrays with different metadata should compare consistently
  // The exact ordering doesn't matter as long as it's consistent
  bool less_than = fda1 < fda2;
  bool greater_than = fda1 > fda2;
  TEST_TRUE(less_than != greater_than) // One should be true, one false
  TEST_FALSE(fda1 == fda2)
  TEST_TRUE(fda1 != fda2)
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
