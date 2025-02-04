// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

START_TEST(DataFilters, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef BaseFeature::QualityType QualityType;

///constructor and destructor test
DataFilters* ptr;
DataFilters* nullPointer = nullptr;
START_SECTION((DataFilters()))
	ptr = new DataFilters();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(([EXTRA]~DataFilters()))
	delete ptr;	
END_SECTION

DataFilters::DataFilter* ptr2 = nullptr;
DataFilters::DataFilter* nullPointer2 = nullptr;

START_SECTION(([EXTRA]DataFilters::DataFilter()))
	ptr2 = new DataFilters::DataFilter();
  TEST_NOT_EQUAL(ptr2, nullPointer2)
END_SECTION

START_SECTION(([EXTRA]~DataFilters::DataFilter()))
	delete ptr2;	
END_SECTION


DataFilters::DataFilter filter_1;
DataFilters::DataFilter filter_2;
DataFilters::DataFilter filter_3;
DataFilters::DataFilter filter_4;
DataFilters::DataFilter filter_5;
DataFilters::DataFilter filter_6;
DataFilters::DataFilter filter_7;
DataFilters::DataFilter filter_8;
DataFilters::DataFilter filter_9;
DataFilters::DataFilter filter_10;
DataFilters::DataFilter filter_11;
DataFilters::DataFilter filter_12;


START_SECTION(([EXTRA]void DataFilter::fromString(const String& filter)))

	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString(""), "the value '' was used but is not valid; invalid filter format")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("not_enough_arguments"), "the value 'not_enough_arguments' was used but is not valid; invalid filter format")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("invalid_fieldname = 0"), "the value 'invalid_fieldname' was used but is not valid; invalid field name")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("Intensity invalid_operator 5"), "the value 'invalid_operator' was used but is not valid; invalid operator")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("Meta::test = string without enclosing quotation marks"), "the value 'string without enclosing quotation marks' was used but is not valid; invalid value")
	//second argument of binary relation missing:
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("Charge = "), "the value '=' was used but is not valid; invalid filter format")
	//string value and non-meta field:
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("Quality = \"a string\""), "the value 'a string' was used but is not valid; invalid value")
	//operation "exists" and non-meta field:
	TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidValue, filter_1.fromString("Intensity exists"), "the value 'exists' was used but is not valid; invalid operator")
	
	filter_1.fromString("Intensity <= 201.334");
	filter_2.fromString("Intensity >= 1000");
	filter_3.fromString("Charge = 4");
	filter_4.fromString("Quality <= 1.0");
	filter_5.fromString("Meta::test_int <= 0");
	filter_6.fromString("Meta::test_double = 0");
	filter_7.fromString("Meta::test_string = \"hello world 2\"");
	filter_8.fromString("Meta::test_dummy exists");
	//valid, but nonsense (nothing will pass this filter):
	filter_9.fromString("Meta::test_string >= \"a string\"");
	filter_10.fromString("Meta::test_string = \"hello world 2\"");
	filter_11.fromString("Meta::unknown_metavalue = 5");
	filter_12.fromString("Meta::test_dummy2 exists");
	
END_SECTION


START_SECTION(([EXTRA]String DataFilter::toString() const))
	
	TEST_STRING_EQUAL(filter_1.toString(), "Intensity <= 201.334000000000003")
	TEST_STRING_EQUAL(filter_2.toString(), "Intensity >= 1000.0")
	TEST_STRING_EQUAL(filter_3.toString(), "Charge = 4.0")
	TEST_STRING_EQUAL(filter_4.toString(), "Quality <= 1.0")
	TEST_STRING_EQUAL(filter_5.toString(), "Meta::test_int <= 0.0")
	TEST_STRING_EQUAL(filter_6.toString(), "Meta::test_double = 0.0")
	TEST_STRING_EQUAL(filter_7.toString(), "Meta::test_string = \"hello world 2\"")
	TEST_STRING_EQUAL(filter_8.toString(), "Meta::test_dummy exists")
	TEST_STRING_EQUAL(filter_9.toString(), "Meta::test_string >= \"a string\"")

END_SECTION


START_SECTION(([EXTRA]bool DataFilter::operator==(const DataFilter& rhs) const))

	TEST_TRUE(filter_10 == filter_7)
	TEST_EQUAL(filter_1 == filter_2, false)
	TEST_TRUE(filter_3 == filter_3)

END_SECTION


START_SECTION(([EXTRA]bool DataFilter::operator!=(const DataFilter& rhs) const))

	TEST_EQUAL(filter_10 != filter_7, false)
	TEST_FALSE(filter_3 == filter_4)
	TEST_EQUAL(filter_4 != filter_4, false)
	
END_SECTION

START_SECTION((bool isActive() const))
	DataFilters tmp;
	TEST_EQUAL(tmp.isActive(), false)
END_SECTION

START_SECTION((void setActive(bool is_active)))
	DataFilters tmp;
	tmp.setActive(true);
	TEST_EQUAL(tmp.isActive(), true)
END_SECTION

DataFilters filters;

START_SECTION((void add(const DataFilter& filter)))

	filters.add(filter_1);
	filters.add(filter_2);
	filters.add(filter_3);
	
	TEST_EQUAL(filters[0] == filter_1, true)
	TEST_EQUAL(filters[1] == filter_2, true)
	TEST_EQUAL(filters[2] == filter_3, true)
	
END_SECTION


START_SECTION((const DataFilter& operator[](Size index) const ))
	
	TEST_EXCEPTION(Exception::IndexOverflow, filters[3])
	filters.add(filter_1);
	TEST_EQUAL(filters[0] == filters[3], true)
	filters.remove(3);
	
END_SECTION


START_SECTION((Size size() const))

	TEST_EQUAL(filters.size(), 3)
	filters.add(filter_4);
	TEST_EQUAL(filters.size(), 4)
	filters.add(filter_5);
	filters.add(filter_6);
	filters.add(filter_7);
	filters.add(filter_8);
	filters.add(filter_9);
	TEST_EQUAL(filters.size(), 9)
	filters.remove(0);
	TEST_EQUAL(filters.size(), 8)
	filters.remove(0);
	TEST_EQUAL(filters.size(), 7)

END_SECTION


START_SECTION((void remove(Size index)))

	TEST_EXCEPTION(Exception::IndexOverflow, filters.remove(7))
	filters.remove(0);
	TEST_EQUAL(filters[0] == filter_4, true)
	filters.remove(0);
	TEST_EQUAL(filters[0] == filter_5, true)
	
END_SECTION


START_SECTION((void replace(Size index, const DataFilter &filter)))
	
	TEST_EXCEPTION(Exception::IndexOverflow, filters.replace(10, filter_1))
	//at the moment: filters[0] == filter_5, ..., filters[4] == filter_9
	filters.replace(0, filter_1);
	filters.replace(1, filter_2);
	filters.replace(2, filter_3);
	filters.replace(3, filter_4);
	filters.replace(4, filter_5);
	TEST_EQUAL(filters[0] == filter_1, true)
	TEST_EQUAL(filters[1] == filter_2, true)
	TEST_EQUAL(filters[2] == filter_3, true)
	TEST_EQUAL(filters[3] == filter_4, true)
	TEST_EQUAL(filters[4] == filter_5, true)
	TEST_EQUAL(filters.size(), 5)
	
END_SECTION


START_SECTION((void clear()))
	
	filters.clear();
	TEST_EQUAL(filters.size(), 0)

END_SECTION


///construct some test features
Feature feature_1;
feature_1.setIntensity(1000.00f);
feature_1.setCharge(4);
feature_1.setOverallQuality((QualityType)31.3334);
feature_1.setMetaValue(String("test_int"), 5);
feature_1.setMetaValue(String("test_double"), 23.42);
feature_1.setMetaValue(String("test_string"), String("hello world 1"));

Feature feature_2;
feature_2.setIntensity(122.01f);
feature_2.setCharge(3);
feature_2.setOverallQuality((QualityType)0.002);
feature_2.setMetaValue(String("test_int"), 10);
feature_2.setMetaValue(String("test_double"), 0.042);
feature_2.setMetaValue(String("test_string"), String("hello world 2"));

Feature feature_3;
feature_3.setIntensity(55.0f);
feature_3.setCharge(4);
feature_3.setOverallQuality((QualityType) 1.);
feature_3.setMetaValue(String("test_int"), 0);
feature_3.setMetaValue(String("test_double"), 100.01);
feature_3.setMetaValue(String("test_string"), String("hello world 3"));

///construct some test consensus features
ConsensusFeature c_feature_1;
c_feature_1.setIntensity(1000.00f);
c_feature_1.setCharge(4);
c_feature_1.setQuality((QualityType) 31.3334);

ConsensusFeature c_feature_2;
c_feature_2.setIntensity(122.01f);
c_feature_2.setCharge(3);
c_feature_2.setQuality((QualityType) 0.002);

ConsensusFeature c_feature_3;
c_feature_3.setIntensity(55.0f);
c_feature_3.setCharge(4);
c_feature_3.setQuality((QualityType) 1.);

///construct some test peaks
MSSpectrum spec;
Peak1D peak;
peak.setIntensity(201.334f);
spec.push_back(peak);
peak.setIntensity(2008.2f);
spec.push_back(peak);
peak.setIntensity(0.001f);
spec.push_back(peak);

MSSpectrum::FloatDataArrays& mdas = spec.getFloatDataArrays();
mdas.resize(3);

mdas[0].setName("test_int");
mdas[0].resize(3);
mdas[0][0] = 5;
mdas[0][1] = 10;
mdas[0][2] = 0;

mdas[1].setName("test_double");
mdas[1].resize(3);
mdas[1][0] =  23.42f;
mdas[1][1] = 0.0f;
mdas[1][2] = 100.01f;

mdas[2].setName("test_dummy");
mdas[2].resize(3);

START_SECTION((template < class PeakType > bool passes(const MSSpectrum &spectrum, Size peak_index) const ))

	filters.add(filter_1); // "Intensity <= 201.334"
	TEST_EQUAL(filters.passes(spec,0), true) // 201.334
	TEST_EQUAL(filters.passes(spec,1), false) // 2008.2
	TEST_EQUAL(filters.passes(spec,2), true) // 0.001
	
	filters.add(filter_2); // "Intensity <= 201.334" && "Intensity >= 1000"
	TEST_EQUAL(filters.passes(spec,0), false) // 201.334
	TEST_EQUAL(filters.passes(spec,1), false) // 2008.2
	TEST_EQUAL(filters.passes(spec,2), false) // 0.001
	
	filters.remove(0); // "Intensity >= 1000"
	TEST_EQUAL(filters.passes(spec,0), false) // 201.334
	TEST_EQUAL(filters.passes(spec,1), true) // 2008.2
	TEST_EQUAL(filters.passes(spec,2), false) // 0.001
	
	filters.clear();
	filters.add(filter_5); // "Meta::test_int <= 0"
	TEST_EQUAL(filters.passes(spec,0), false) // 5
	TEST_EQUAL(filters.passes(spec,1), false) // 10
	TEST_EQUAL(filters.passes(spec,2), true) // 0
	
	filters.clear();
	filters.add(filter_8); // Meta::test_dummy exists
	TEST_EQUAL(filters.passes(spec,0), true)
	TEST_EQUAL(filters.passes(spec,1), true)
	TEST_EQUAL(filters.passes(spec,2), true)

	filters.clear();
	filters.add(filter_12); // Meta::test_dummy2 exists
	TEST_EQUAL(filters.passes(spec,0), false)
	TEST_EQUAL(filters.passes(spec,1), false)
	TEST_EQUAL(filters.passes(spec,2), false)


	filters.clear();
	filters.add(filter_6); // Meta::test_double = 0
	TEST_EQUAL(filters.passes(spec,0), false)
	TEST_EQUAL(filters.passes(spec,1), true)
	TEST_EQUAL(filters.passes(spec,2), false)

END_SECTION


START_SECTION((bool passes(const Feature& feature) const))

	filters.clear();
	filters.add(filter_3); // "Charge = 4"
	TEST_EQUAL(filters.passes(feature_1), true) // 4
	TEST_EQUAL(filters.passes(feature_2), false) // 3
	TEST_EQUAL(filters.passes(feature_3), true) // 4
	
	filters.add(filter_4); // "Quality <= 1.0" && "Charge = 4"
	TEST_EQUAL(filters.passes(feature_1), false) // Quality = 31.3334; Charge = 4
	TEST_EQUAL(filters.passes(feature_2), false) // Quality = 0.002; Charge = 3
	TEST_EQUAL(filters.passes(feature_3), true) // Quality = 1; Charge = 4
	
	filters.remove(0); // "Quality <= 1.0"
	TEST_EQUAL(filters.passes(feature_1), false) // Quality = 31.3334
	TEST_EQUAL(filters.passes(feature_2), true) // Quality = 0.002
	TEST_EQUAL(filters.passes(feature_3), true) // Quality = 1
	
	filters.clear();
	filters.add(filter_2); // "Intensity >= 1000"
	TEST_EQUAL(filters.passes(feature_1), true) // 1000.00
	TEST_EQUAL(filters.passes(feature_2), false) // 122.01
	TEST_EQUAL(filters.passes(feature_3), false) // 55.0
	
	filters.clear();
	filters.add(filter_7); // "Meta::test_string = \"hello world 2\""
	TEST_EQUAL(filters.passes(feature_1), false)
	TEST_EQUAL(filters.passes(feature_2), true)
	TEST_EQUAL(filters.passes(feature_3), false)
	
	filters.add(filter_8); // "Meta::test_dummy exists"
	TEST_EQUAL(filters.passes(feature_1), false)
	TEST_EQUAL(filters.passes(feature_2), false)
	TEST_EQUAL(filters.passes(feature_3), false)
	
	filters.clear();
	filters.add(filter_5); // "Meta::test_int <= 0"
	TEST_EQUAL(filters.passes(feature_1), false) // 5
	TEST_EQUAL(filters.passes(feature_2), false) // 10
	TEST_EQUAL(filters.passes(feature_3), true) // 0

END_SECTION

START_SECTION((bool passes(const ConsensusFeature& consensus_feature) const))

	filters.clear();
	filters.add(filter_3); // "Charge = 4"
	TEST_EQUAL(filters.passes(c_feature_1), true) // 4
	TEST_EQUAL(filters.passes(c_feature_2), false) // 3
	TEST_EQUAL(filters.passes(c_feature_3), true) // 4
	
	filters.add(filter_4); // "Quality <= 1.0" && "Charge = 4"
	TEST_EQUAL(filters.passes(c_feature_1), false) // Quality = 31.3334; Charge = 4
	TEST_EQUAL(filters.passes(c_feature_2), false) // Quality = 0.002; Charge = 3
	TEST_EQUAL(filters.passes(c_feature_3), true) // Quality = 1; Charge = 4
	
	filters.remove(0); // "Quality <= 1.0"
	TEST_EQUAL(filters.passes(c_feature_1), false) // Quality = 31.3334
	TEST_EQUAL(filters.passes(c_feature_2), true) // Quality = 0.002
	TEST_EQUAL(filters.passes(c_feature_3), true) // Quality = 1
	
	filters.clear();
	filters.add(filter_2); // "Intensity >= 1000"
	TEST_EQUAL(filters.passes(c_feature_1), true) // 1000.00
	TEST_EQUAL(filters.passes(c_feature_2), false) // 122.01
	TEST_EQUAL(filters.passes(c_feature_3), false) // 55.0

END_SECTION

DataFilters::DataFilter* df_ptr;
START_SECTION(([DataFilters::DataFilter] DataFilter()))
{
  df_ptr = new DataFilters::DataFilter();
  TEST_NOT_EQUAL(df_ptr, nullPointer2)

  delete df_ptr;
}
END_SECTION

START_SECTION(([DataFilters::DataFilter] String toString() const ))
{
  DataFilters::DataFilter df1;
  df1.field = DataFilters::INTENSITY;
  df1.op = DataFilters::LESS_EQUAL;
  df1.value = 25.5;

  TEST_EQUAL(df1.toString(), "Intensity <= 25.5")

  df1.field = DataFilters::META_DATA;
  df1.meta_name = "meta-value";
  df1.op = DataFilters::EXISTS;
  df1.value_is_numerical = false;

  TEST_EQUAL(df1.toString(), "Meta::meta-value exists")

  df1.op = DataFilters::EQUAL;
  df1.value_string = "value";
  TEST_EQUAL(df1.toString(), "Meta::meta-value = \"value\"")
}
END_SECTION

START_SECTION(([DataFilters::DataFilter] void fromString(const String &filter)))
{
  DataFilters::DataFilter df1;
  df1.fromString("Intensity <= 25.3");
  TEST_EQUAL(df1.field, DataFilters::INTENSITY)
  TEST_EQUAL(df1.op, DataFilters::LESS_EQUAL)
  TEST_EQUAL(df1.value, 25.3)
  TEST_EQUAL(df1.value_is_numerical, true)

  DataFilters::DataFilter df2;
  df2.fromString("Meta::meta-value exists");
  TEST_EQUAL(df2.field, DataFilters::META_DATA)
  TEST_EQUAL(df2.op, DataFilters::EXISTS)
  TEST_EQUAL(df2.meta_name, "meta-value")

  DataFilters::DataFilter df3;
  df3.fromString("Meta::meta-value = \"value\"");
  TEST_EQUAL(df3.field, DataFilters::META_DATA)
  TEST_EQUAL(df3.op, DataFilters::EQUAL)
  TEST_EQUAL(df3.meta_name, "meta-value")
  TEST_EQUAL(df3.value_string, "value")
  TEST_EQUAL(df3.value_is_numerical, false)

  // test some wrong cases
  DataFilters::DataFilter exception_filter;
  TEST_EXCEPTION(Exception::InvalidValue, exception_filter.fromString("Intensity <> 24.5"))
  TEST_EXCEPTION(Exception::InvalidValue, exception_filter.fromString("Intensity < 24.5"))
  TEST_EXCEPTION(Exception::InvalidValue, exception_filter.fromString("Insenity = 2.0"))
  TEST_EXCEPTION(Exception::InvalidValue, exception_filter.fromString("Charge = text-value"))
}
END_SECTION

START_SECTION(([DataFilters::DataFilter] bool operator==(const DataFilter &rhs) const ))
{
  DataFilters::DataFilter df1,df2,df3;

  TEST_TRUE(df1 == df2)

  // field
  df1.field = DataFilters::CHARGE;
  df2.field = DataFilters::CHARGE;
  df3.field = DataFilters::INTENSITY;

  TEST_EQUAL(df1==df2,true)
  TEST_EQUAL(df1==df3,false)
  df3.field = DataFilters::CHARGE;

  // op
  df1.op = DataFilters::EQUAL;
  df2.op = DataFilters::EQUAL;
  df3.op = DataFilters::GREATER_EQUAL;

  TEST_EQUAL(df1==df2,true)
  TEST_EQUAL(df1==df3,false)
  df3.op = DataFilters::EQUAL;

  // value_is_numerical
  df1.value = 0.0;
  df2.value = 0.0;
  df3.value = 0.2;

  TEST_EQUAL(df1==df2,true)
  TEST_EQUAL(df1==df3,false)

  df1.meta_name = "df1";
  df2.meta_name = "df1";

  TEST_EQUAL(df1==df2,true)
  df2.meta_name = "df2";
  TEST_EQUAL(df1==df2,false)
  df2.meta_name = "df1";

  df1.value_string = "df1";
  df2.value_string = "df1";
  TEST_EQUAL(df1==df2,true)
  df2.value_string = "df2";
  TEST_EQUAL(df1==df2,false)
  df2.value_string = "df1";

  df1.value_is_numerical = true;
  df2.value_is_numerical = true;
  TEST_EQUAL(df1==df2,true)
  df2.value_is_numerical = false;
  TEST_EQUAL(df1==df2,false)
}
END_SECTION

START_SECTION(([DataFilters::DataFilter] bool operator!=(const DataFilter &rhs) const ))
{
  DataFilters::DataFilter df1,df2,df3;

  TEST_TRUE(df1 == df2)

  // field
  df1.field = DataFilters::CHARGE;
  df2.field = DataFilters::CHARGE;
  df3.field = DataFilters::INTENSITY;

  TEST_EQUAL(df1!=df2,false)
  TEST_EQUAL(df1!=df3,true)
  df3.field = DataFilters::CHARGE;

  // op
  df1.op = DataFilters::EQUAL;
  df2.op = DataFilters::EQUAL;
  df3.op = DataFilters::GREATER_EQUAL;

  TEST_EQUAL(df1!=df2,false)
  TEST_EQUAL(df1!=df3,true)
  df3.op = DataFilters::EQUAL;

  // value_is_numerical
  df1.value = 0.0;
  df2.value = 0.0;
  df3.value = 0.2;

  TEST_EQUAL(df1!=df2,false)
  TEST_EQUAL(df1!=df3,true)

  df1.meta_name = "df1";
  df2.meta_name = "df1";

  TEST_EQUAL(df1!=df2,false)
  df2.meta_name = "df2";
  TEST_EQUAL(df1!=df2,true)
  df2.meta_name = "df1";

  df1.value_string = "df1";
  df2.value_string = "df1";
  TEST_EQUAL(df1!=df2,false)
  df2.value_string = "df2";
  TEST_EQUAL(df1!=df2,true)
  df2.value_string = "df1";

  df1.value_is_numerical = true;
  df2.value_is_numerical = true;
  TEST_EQUAL(df1!=df2,false)
  df2.value_is_numerical = false;
  TEST_EQUAL(df1!=df2,true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
