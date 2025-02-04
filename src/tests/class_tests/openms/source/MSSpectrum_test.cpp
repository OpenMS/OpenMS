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
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

#include <OpenMS/IONMOBILITY/IMDataConverter.h>

#include <sstream>

using namespace OpenMS;
using namespace std;

static_assert(OpenMS::Test::fulfills_rule_of_5<MSSpectrum>(), "Must fulfill rule of 5");
static_assert(OpenMS::Test::fulfills_rule_of_6<MSSpectrum>(), "Must fulfill rule of 6");
static_assert(OpenMS::Test::fulfills_fast_vector<MSSpectrum>(), "Must have fast vector semantics");
static_assert(std::is_nothrow_move_constructible_v<MSSpectrum>, "Must have nothrow move constructible");

/// A spec with RT, m/z, intensity, and meta data arrays, marked as an IM spectrum (i.e. spec.containsIMData() == true)
MSSpectrum getPrefilledSpec()
{
  MSSpectrum ds;
  MSSpectrum::FloatDataArray float_array {56.0, 201.0, 31, 31, 31, 37, 29, 34, 60, 29};
  MSSpectrum::StringDataArray string_array {"56", "201", "31", "31", "31", "37", "29", "34", "60", "29"};
  MSSpectrum::IntegerDataArray int_array {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};
  std::vector<double> mzs {423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321};
  std::vector<double> intensities {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3, float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  IMDataConverter::setIMUnit(ds.getFloatDataArrays()[1], DriftTimeUnit::MILLISECOND);
  TEST_TRUE(ds.containsIMData())

  ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(2, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");
  
  ds.setRT(5.0);
  return ds;
}

START_TEST(MSSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
// Dummy peak data

Peak1D p1;
p1.setIntensity(1.0f);
p1.setMZ(2.0);

Peak1D p2;
p2.setIntensity(2.0f);
p2.setMZ(10.0);

Peak1D p3;
p3.setIntensity(3.0f);
p3.setMZ(30.0);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum* ptr = nullptr;
MSSpectrum* nullPointer = nullptr;
START_SECTION((MSSpectrum()))
{
  ptr = new MSSpectrum();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MSSpectrum()))
{
  delete ptr;
}
END_SECTION

START_SECTION(([EXTRA] MSSpectrum()))
{
  MSSpectrum tmp;
  Peak1D peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  TEST_EQUAL(tmp.size(),1);
  TEST_REAL_SIMILAR(tmp[0].getMZ(), 47.11);
}
END_SECTION

START_SECTION((MSSpectrum(const std::initializer_list<Peak1D>& init)))
{
  MSSpectrum tmp {{47.11, 2}, {500.0, 3}};
  TEST_EQUAL(tmp.size(), 2);
  TEST_REAL_SIMILAR(tmp[0].getMZ(), 47.11);
  TEST_REAL_SIMILAR(tmp[1].getMZ(), 500.0);
  TEST_REAL_SIMILAR(tmp[0].getIntensity(), 2)
  TEST_REAL_SIMILAR(tmp[1].getIntensity(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Member accessors

START_SECTION((UInt getMSLevel() const))
{
  MSSpectrum spec;
  TEST_EQUAL(spec.getMSLevel(), 1)
}
END_SECTION

START_SECTION((void setMSLevel(UInt ms_level)))
{
  MSSpectrum spec;
  spec.setMSLevel(17);
  TEST_EQUAL(spec.getMSLevel(), 17)
}
END_SECTION

START_SECTION((const String& getName() const))
{
  MSSpectrum s;
  TEST_STRING_EQUAL(s.getName(), "")
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  MSSpectrum s;
  s.setName("bla");
  TEST_STRING_EQUAL(s.getName(),"bla")
}
END_SECTION

START_SECTION((double getRT() const ))
{
  MSSpectrum s;
  TEST_REAL_SIMILAR(s.getRT(), -1.0)
}
END_SECTION

START_SECTION((void setRT(double rt)))
{
  MSSpectrum s;
  s.setRT(0.451);
  TEST_REAL_SIMILAR(s.getRT(), 0.451)
}
END_SECTION

START_SECTION((double getDriftTime() const ))
{
  MSSpectrum s;
  TEST_REAL_SIMILAR(s.getDriftTime(), -1.0)
}
END_SECTION

START_SECTION((void setDriftTime(double dt)))
{
  MSSpectrum s;
  s.setDriftTime(0.451);
  TEST_REAL_SIMILAR(s.getDriftTime(), 0.451)
}
END_SECTION

START_SECTION((double getDriftTimeUnit() const ))
{
  MSSpectrum s;
  TEST_EQUAL(s.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
}
END_SECTION

START_SECTION((double getDriftTimeUnitAsString() const))
{
  MSSpectrum s;
  TEST_EQUAL(s.getDriftTimeUnitAsString(), "<NONE>");
}
END_SECTION


START_SECTION((void setDriftTimeUnit(double dt)))
{
  MSSpectrum s;
  s.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(s.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(s.getDriftTimeUnitAsString(), "ms");
}
END_SECTION

START_SECTION((const FloatDataArrays& getFloatDataArrays() const))
{
  MSSpectrum s;
  TEST_EQUAL(s.getFloatDataArrays().size(),0)
}
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
{
  MSSpectrum s;
  s.getFloatDataArrays().resize(2);
  TEST_EQUAL(s.getFloatDataArrays().size(),2)
}
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const))
{
  MSSpectrum s;
  TEST_EQUAL(s.getStringDataArrays().size(),0)
}
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
{
  MSSpectrum s;
  s.getStringDataArrays().resize(2);
  TEST_EQUAL(s.getStringDataArrays().size(),2)
}
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const))
{
  MSSpectrum s;
  TEST_EQUAL(s.getIntegerDataArrays().size(),0)
}
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
{
  MSSpectrum s;
  s.getIntegerDataArrays().resize(2);
  TEST_EQUAL(s.getIntegerDataArrays().size(),2)
}
END_SECTION

START_SECTION((MSSpectrum& select(const std::vector<Size>& indices)))
{
  MSSpectrum s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p3);
  s.push_back(p3);
  s.push_back(p2);

  MSSpectrum::IntegerDataArray aia{1, 2, 3, 4, 5};
  MSSpectrum::FloatDataArray afa{1.0, 2.0, 3.0, 4.0, 5.0};
  MSSpectrum::StringDataArray asa{"1", "2", "3", "4", "5"};
  s.getFloatDataArrays().push_back(afa);
  s.getIntegerDataArrays().push_back(aia);
  s.getStringDataArrays().push_back(asa);
  s.getFloatDataArrays().push_back(afa);
  s.getIntegerDataArrays().push_back(aia);
  s.getStringDataArrays().push_back(asa);

  TEST_REAL_SIMILAR(s[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(s[4].getIntensity(), 2.0)
  TEST_EQUAL(s.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s.getFloatDataArrays()[0].size(), 5)
  TEST_EQUAL(s.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s.getIntegerDataArrays()[0].size(), 5)
  TEST_EQUAL(s.getStringDataArrays().size(), 2)
  TEST_EQUAL(s.getStringDataArrays()[0].size(), 5)

  // re-order
  MSSpectrum s2 = s;
  Size order[] = {4, 2, 3, 1, 0};
  s2.select(std::vector<Size>(&order[0], &order[5]));
  TEST_REAL_SIMILAR(s2[0].getIntensity(), 2.0)
  TEST_REAL_SIMILAR(s2[4].getIntensity(), 1.0)
  TEST_EQUAL(s2.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s2.getFloatDataArrays()[0].size(), 5)
  TEST_EQUAL(s2.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s2.getIntegerDataArrays()[0].size(), 5)
  TEST_EQUAL(s2.getStringDataArrays().size(), 2)
  TEST_EQUAL(s2.getStringDataArrays()[0].size(), 5)

  TEST_REAL_SIMILAR(s2.getFloatDataArrays()[0][1], 3.0)
  TEST_EQUAL(s2.getIntegerDataArrays()[0][1], 3)
  TEST_EQUAL(s2.getStringDataArrays()[0][1], "3")

  // subset
  s2 = s;
  Size subset[] = {4, 2, 3};
  // --> new values in Meta arrays are:
  //     5, 3, 4
  s2.select(std::vector<Size>(&subset[0], &subset[3]));
  TEST_REAL_SIMILAR(s2[0].getIntensity(), 2.0)
  TEST_REAL_SIMILAR(s2[1].getIntensity(), 3.0)
  TEST_REAL_SIMILAR(s2[2].getIntensity(), 3.0)
  TEST_EQUAL(s2.getFloatDataArrays().size(), 2)
  TEST_EQUAL(s2.getFloatDataArrays()[0].size(), 3)
  TEST_EQUAL(s2.getIntegerDataArrays().size(), 2)
  TEST_EQUAL(s2.getIntegerDataArrays()[0].size(), 3)
  TEST_EQUAL(s2.getStringDataArrays().size(), 2)
  TEST_EQUAL(s2.getStringDataArrays()[0].size(), 3)

  TEST_REAL_SIMILAR(s2.getFloatDataArrays()[0][1], 3.0)
  TEST_EQUAL(s2.getIntegerDataArrays()[0][1], 3)
  TEST_EQUAL(s2.getStringDataArrays()[0][1], "3")
}
END_SECTION

/////////////////////////////////////////////////////////////
// RangeManager

START_SECTION((virtual void updateRanges()))
{
  MSSpectrum s = getPrefilledSpec();

  for (int i = 0; i < 2; ++i) // second time to check the initialization
  {
    s.updateRanges();
    TEST_REAL_SIMILAR(s.getMinIntensity(), 29)
    TEST_REAL_SIMILAR(s.getMaxIntensity(), 201)
    TEST_REAL_SIMILAR(s.getMinMZ(), 412.321)
    TEST_REAL_SIMILAR(s.getMaxMZ(), 423.269)
    TEST_REAL_SIMILAR(s.getMinMobility(), 29)
    TEST_REAL_SIMILAR(s.getMaxMobility(), 201)
  }

  //test with only one peak
  s = MSSpectrum{};
  s.push_back(p1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMinIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMaxMZ(), 2)
  TEST_REAL_SIMILAR(s.getMinMZ(), 2)
  TEST_TRUE(s.RangeMobility::isEmpty())
}
END_SECTION

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((MSSpectrum(const MSSpectrum& source)))
{
  MSSpectrum tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  tmp.setMSLevel(17);
  tmp.setRT(7.0);
  tmp.setDriftTime(8.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  tmp.setName("bla");
  //peaks
  MSSpectrum::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  
  MSSpectrum tmp2(tmp);
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 17)
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 8.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(tmp2.getName(),"bla")
  //peaks
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
}
END_SECTION

START_SECTION((MSSpectrum(const MSSpectrum&& source)))
{
  // Ensure that MSSpectrum has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(MSSpectrum(std::declval<MSSpectrum&&>())), true)

  MSSpectrum tmp;
  tmp.setRT(9.0);
  tmp.setDriftTime(5.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  tmp.setMSLevel(18);
  tmp.setName("bla2");
  tmp.setMetaValue("label2",5.0);
  tmp.getInstrumentSettings().getScanWindows().resize(2);
  //peaks
  MSSpectrum::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  peak.getPosition()[0] = 48.11;
  tmp.push_back(peak);
  
  //copy tmp so we can move one of them
  MSSpectrum orig = tmp;
  MSSpectrum tmp2(std::move(tmp));

  TEST_EQUAL(tmp2, orig); // should be equal to the original

  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),2);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label2"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 18)
  TEST_REAL_SIMILAR(tmp2.getRT(), 9.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 5.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true);
  TEST_EQUAL(tmp2.getName(),"bla2")
  TEST_EQUAL(tmp2.size(),2);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
  TEST_REAL_SIMILAR(tmp2[1].getPosition()[0],48.11);

  // test move
  TEST_EQUAL(tmp.size(),0);
  TEST_EQUAL(tmp.metaValueExists("label2"), false);
}
END_SECTION

START_SECTION((MSSpectrum& operator= (const MSSpectrum& source)))
{
  MSSpectrum tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.setMetaValue("label",5.0);
  tmp.setMSLevel(17);
  tmp.setRT(7.0);
  tmp.setDriftTime(8.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  tmp.setName("bla");
  //peaks
  MSSpectrum::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  //normal assignment
  MSSpectrum tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(), 1);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 17)
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 8.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(tmp2.getName(), "bla")
  TEST_EQUAL(tmp2.size(), 1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0], 47.11);

  //Assignment of empty object
  //normal assignment
  tmp2 = MSSpectrum();
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(), 0);
  TEST_EQUAL(tmp2.metaValueExists("label"), false)
  TEST_EQUAL(tmp2.getMSLevel(),1)
  TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), -1.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
  TEST_EQUAL(tmp2.getName(), "")
  TEST_EQUAL(tmp2.size(), 0);
}
END_SECTION

START_SECTION((MSSpectrum& operator= (const MSSpectrum&& source)))
{
  MSSpectrum tmp {{47.11, 0}, {48.11, 0}};
  tmp.setRT(9.0);
  tmp.setDriftTime(5.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  tmp.setMSLevel(18);
  tmp.setName("bla2");
  tmp.setMetaValue("label2",5.0);
  tmp.getInstrumentSettings().getScanWindows().resize(2);

  //copy tmp so we can move one of them
  MSSpectrum orig = tmp;

  //move assignment
  MSSpectrum tmp2;
  tmp2 = std::move(tmp);

  TEST_EQUAL(tmp2, orig); // should be equal to the original

  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),2);
  TEST_REAL_SIMILAR(tmp2.getMetaValue("label2"), 5.0)
  TEST_EQUAL(tmp2.getMSLevel(), 18)
  TEST_REAL_SIMILAR(tmp2.getRT(), 9.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), 5.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true);
  TEST_EQUAL(tmp2.getName(),"bla2")
  TEST_EQUAL(tmp2.size(),2);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
  TEST_REAL_SIMILAR(tmp2[1].getPosition()[0],48.11);

  // test move
  TEST_EQUAL(tmp.size(),0);
  TEST_EQUAL(tmp.metaValueExists("label2"), false);

  //Assignment of empty object
  //normal assignment
#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic push
// Ignore -Wpessimizing-move, because we want to test the move assignment operator.
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
  tmp2 = std::move(MSSpectrum());
#ifndef OPENMS_WINDOWSPLATFORM
#pragma clang diagnostic pop
#endif
  TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),0);
  TEST_FALSE(tmp2.metaValueExists("label"))
  TEST_EQUAL(tmp2.getMSLevel(),1)
  TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
  TEST_REAL_SIMILAR(tmp2.getDriftTime(), -1.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
  TEST_EQUAL(tmp2.getName(),"")
  TEST_EQUAL(tmp2.size(),0);
}
END_SECTION

START_SECTION((bool operator== (const MSSpectrum& rhs) const))
{
  MSSpectrum edit, empty;
  
  TEST_TRUE(edit==empty);

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_FALSE(edit==empty);

  edit = empty;
  edit.resize(1);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setMetaValue("label",String("bla"));
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.setDriftTime(5);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.setRT(5);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.setMSLevel(5);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.getFloatDataArrays().resize(5);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.getStringDataArrays().resize(5);
  TEST_FALSE(empty == edit);

  edit = empty;
  edit.getIntegerDataArrays().resize(5);
  TEST_FALSE(empty == edit);

  //name is not checked => no change
  edit = empty;
  edit.setName("bla");
  TEST_TRUE(empty == edit);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear(false);
  TEST_TRUE(empty == edit);
}
END_SECTION

START_SECTION((bool operator!= (const MSSpectrum& rhs) const))
{
  MSSpectrum edit, empty;
  
  TEST_FALSE(edit != empty);

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.resize(1);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.setMetaValue("label",String("bla"));
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.setDriftTime(5);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.setRT(5);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.setMSLevel(5);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.getFloatDataArrays().resize(5);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.getIntegerDataArrays().resize(5);
  TEST_TRUE(edit != empty);

  edit = empty;
  edit.getStringDataArrays().resize(5);
  TEST_TRUE(edit != empty);

  //name is not checked => no change
  edit = empty;
  edit.setName("bla");
  TEST_FALSE(edit != empty);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear(false);
  TEST_TRUE(edit == empty);
}
END_SECTION


/////////////////////////////////////////////////////////////
// Sorting

START_SECTION((void sortByIntensity(bool reverse=false)))
{
  MSSpectrum ds;
  Peak1D p;
  MSSpectrum::FloatDataArray float_array { 420.13f, 412.824f, 423.269f, 415.287f, 413.8f, 419.113f, 416.293f, 418.232f, 414.301f, 412.321f };
  MSSpectrum::StringDataArray string_array {"420.13", "412.82", "423.27", "415.29", "413.80", "419.11", "416.29", "418.23", "414.30", "412.32"};
  MSSpectrum::IntegerDataArray int_array {420, 412, 423, 415, 413, 419, 416, 418, 414, 412};
  std::vector<double> mzs {420.130, 412.824, 423.269, 415.287, 413.800, 419.113, 416.293, 418.232, 414.301, 412.321};
  std::vector<double> intensities {201, 60, 56, 37, 34, 31, 31, 31, 29, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }
  ds.sortByIntensity();
  std::vector<double> intensities_copy(intensities);
  std::sort(intensities_copy.begin(), intensities_copy.end());
  MSSpectrum::iterator it_ds = ds.begin();
  ABORT_IF(ds.size() != intensities_copy.size())
  for(std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    TEST_EQUAL(it_ds->getIntensity(), *it);
    ++it_ds;
  }
  ds.clear(true);
  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
    ds.push_back(p);
  }

  ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(1, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByIntensity();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  MSSpectrum::iterator it1 = ds.begin();
  MSSpectrum::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSSpectrum::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSSpectrum::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  TOLERANCE_ABSOLUTE(0.0001)
    for (std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
    {
      if (it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end() && it4 != ds.getIntegerDataArrays()[0].end())
      {
        //metadataarray values == mz values
        TEST_REAL_SIMILAR(it1->getIntensity(), *it);
        TEST_REAL_SIMILAR(*it2 , it1->getMZ());
        TEST_STRING_EQUAL(*it3 , String::number(it1->getMZ(),2));
        TEST_EQUAL(*it4 , (Int)floor(it1->getMZ()));
        ++it1;
        ++it2;
        ++it3;
        ++it4;
      }
      else
      {
        TEST_EQUAL(true,false)
      }
    }
}
END_SECTION



START_SECTION((void sortByPosition()))
{
  MSSpectrum ds;
  MSSpectrum::FloatDataArray float_array {56.0, 201.0, 31, 31, 31, 37, 29, 34, 60, 29};
  MSSpectrum::StringDataArray string_array {"56", "201", "31", "31", "31", "37", "29", "34", "60", "29"};
  MSSpectrum::IntegerDataArray int_array {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};
  std::vector<double> mzs {423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321};
  std::vector<double> intensities {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.sortByPosition();
  MSSpectrum::iterator it = ds.begin();
  for (std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    if(it == ds.end())
    {
      TEST_EQUAL(true,false)
    }
    TEST_EQUAL(it->getIntensity(), *rit);
    ++it;
  }
  ds.clear(true);
  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(2, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByPosition();

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  Size size = intensities.size();
  ABORT_IF(ds.size() != size);
  ABORT_IF(ds.getFloatDataArrays()[1].size() != size);
  ABORT_IF(ds.getStringDataArrays()[0].size() != size);
  ABORT_IF(ds.getIntegerDataArrays()[0].size() != size);
  MSSpectrum::iterator it1 = ds.begin();
  MSSpectrum::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSSpectrum::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSSpectrum::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  for (std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    //metadataarray values == intensity values
    TEST_REAL_SIMILAR(it1->getIntensity(), *rit);
    TEST_REAL_SIMILAR(*it2 , *rit);
    TEST_STRING_EQUAL(*it3 , String::number(*rit,0));
    TEST_EQUAL(*it4 , (Int)floor(*rit));
    ++it1;
    ++it2;
    ++it3;
    ++it4;
  }
}
END_SECTION

START_SECTION(void sortByIonMobility())
{
  auto ds = getPrefilledSpec();

  TEST_FALSE(ds.isSortedByIM())
  ds.sortByIonMobility();
  TEST_TRUE(ds.isSortedByIM())
  auto [idx, unit] = ds.getIMData();
  TEST_EQUAL(idx, 1)
  const auto& im = ds.getFloatDataArrays()[idx];
  TEST_TRUE(std::is_sorted(im.begin(), im.end())) 
}
END_SECTION

START_SECTION(void isSortedByIM() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void sortByPositionPresorted()))
{
  MSSpectrum ds;
  MSSpectrum::FloatDataArray float_array {19, 20, 23, 15, 16, 18, 13, 14, 12, 12};
  MSSpectrum::StringDataArray string_array {"19", "20", "23", "15", "16", "18", "13", "14", "12", "12"};
  MSSpectrum::IntegerDataArray int_array {19, 20, 23, 15, 16, 18, 13, 14, 12, 12};
  std::vector<double> mzs {419.113, 420.130, 423.269, 415.287, 416.293, 418.232, 413.800, 414.301, 412.824, 412.321};
  std::vector<double> intensities {19, 20, 23, 15, 16, 18, 13, 14, 12, 12};

  MSSpectrum::Chunks chunks(ds);
  double last_added = 0;
  for (Size i = 0; i < mzs.size(); ++i)
  {
    if (mzs[i] < last_added) chunks.add(true);
    last_added = mzs[i];
    ds.emplace_back(mzs[i], intensities[i]);
  }
  chunks.add(true); // Add the last chunk

  ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3,float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getFloatDataArrays()[1].setName("f2");
  ds.getFloatDataArrays()[2].setName("f3");

  ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(2, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByPositionPresorted(chunks.getChunks());

  TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
  TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

  TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
  TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

  TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

  Size size = intensities.size();
  ABORT_IF(ds.size() != size);
  ABORT_IF(ds.getFloatDataArrays()[1].size() != size);
  ABORT_IF(ds.getStringDataArrays()[0].size() != size);
  ABORT_IF(ds.getIntegerDataArrays()[0].size() != size);
  MSSpectrum::iterator it1 = ds.begin();
  MSSpectrum::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
  MSSpectrum::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
  MSSpectrum::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
  std::sort(intensities.begin(), intensities.end());
  for (std::vector<double>::iterator it = intensities.begin(); it != intensities.end(); ++it)
  {
    //metadataarray values == intensity values
    TEST_REAL_SIMILAR(it1->getIntensity(), *it);
    TEST_REAL_SIMILAR(*it2 , *it);
    TEST_STRING_EQUAL(*it3 , String::number(*it,0));
    TEST_EQUAL(*it4 , (Int)floor(*it));
    ++it1; ++it2; ++it3; ++it4;
  }
}
END_SECTION

START_SECTION(bool isSorted() const)
{
  //make test dataset
  MSSpectrum spec {{1000.0, 3}, {1001, 5}, {1002, 1}};

  TEST_EQUAL(spec.isSorted(),true)

  reverse(spec.begin(), spec.end());
  TEST_EQUAL(spec.isSorted(),false)
}
END_SECTION

START_SECTION(template<class Predicate>
              bool isSorted(const Predicate& lamdba) const)
{
  MSSpectrum ds;
  MSSpectrum::FloatDataArray float_array {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};
  MSSpectrum::StringDataArray string_array {"56", "201", "31", "31", "31", "37", "29", "34", "60", "29"};
  MSSpectrum::IntegerDataArray int_array {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};
  std::vector<double> mzs{423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321};
  std::vector<double> intensities{56, 201, 31, 31, 31, 37, 29, 34, 60, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.getFloatDataArrays() = std::vector<MSSpectrum::FloatDataArray>(3, float_array);
  ds.getFloatDataArrays()[0].setName("f1");
  ds.getStringDataArrays() = std::vector<MSSpectrum::StringDataArray>(2, string_array);
  ds.getStringDataArrays()[0].setName("s1");
  ds.getStringDataArrays()[1].setName("s2");

  ds.getIntegerDataArrays() = std::vector<MSSpectrum::IntegerDataArray>(1, int_array);
  ds.getIntegerDataArrays()[0].setName("i1");

  ds.sortByPosition();

  // more expensive than isSorted(), but just to make sure
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) {return ds[a].getMZ() < ds[b].getMZ();}), true)
  TEST_EQUAL(ds.isSorted(), true) // call other method. Should give the same result

  ds.sortByIntensity();
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getIntensity() < ds[b].getIntensity(); }), true)
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getMZ() < ds[b].getMZ(); }), false)
  TEST_EQUAL(ds.isSorted(), false)// call other method. Should give the same result

  // sort by metadata array; float data is identical to intensities here, so we can easily check
  auto float_sort_func = [&ds](Size a, Size b) {
    return ds.getFloatDataArrays()[0][a] < ds.getFloatDataArrays()[0][b]; 
  };
  ds.sortByPosition();// make sure the order is wrong before calling .sort(...)
  ds.sort(float_sort_func);
  TEST_EQUAL(ds[0].getIntensity(), 29)
  TEST_EQUAL(ds.isSorted(), false) // not sorted by m/z
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getIntensity() < ds[b].getIntensity(); }), true)

}
END_SECTION

START_SECTION(template<class Predicate> void sort(const Predicate& lambda))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
// Finding peaks or peak ranges

const MSSpectrum spec_find{{1.0, 29.0f}, {2.0, 60.0f}, {3.0, 34.0f}, {4.0, 29.0f}, {5.0, 37.0f}, {6.0, 31.0f}};

START_SECTION((Iterator MZEnd(CoordinateType mz)))
{
  MSSpectrum::Iterator it;
  auto tmp = spec_find;
  it = tmp.MZEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = tmp.MZEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = tmp.MZEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
{
  MSSpectrum::Iterator it;
  auto tmp = spec_find;

  it = tmp.MZBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

START_SECTION((Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)))
{
  MSSpectrum::Iterator it;
  auto tmp = spec_find;

  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
{
  MSSpectrum::Iterator it;
  auto tmp = spec_find;

  it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)))
{
  MSSpectrum::Iterator it;
  auto tmp = spec_find;

  it = tmp.MZEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
    it = tmp.MZEnd(tmp.begin(), 5, tmp.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
    it = tmp.MZEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
{
  MSSpectrum::ConstIterator it;

  it = spec_find.MZEnd(spec_find.begin(), 4.5, spec_find.end());
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = spec_find.MZEnd(spec_find.begin(), 5, spec_find.end());
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = spec_find.MZEnd(spec_find.begin(), 4.5, spec_find.begin());
  TEST_EQUAL(it->getPosition()[0], spec_find.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const))
{
  MSSpectrum::ConstIterator it;

  it = spec_find.MZEnd(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = spec_find.MZEnd(5.0);
  TEST_EQUAL(it->getPosition()[0],6.0)
  it = spec_find.MZEnd(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const))
{
  MSSpectrum::ConstIterator it;

  it = spec_find.MZBegin(4.5);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = spec_find.MZBegin(5.0);
  TEST_EQUAL(it->getPosition()[0],5.0)
  it = spec_find.MZBegin(5.5);
  TEST_EQUAL(it->getPosition()[0],6.0)
}
END_SECTION

auto tmp = spec_find;

START_SECTION((Iterator PosBegin(CoordinateType mz)))
{
  MSSpectrum::Iterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosBegin(Iterator begin, CoordinateType mz, Iterator end)))
{
  MSSpectrum::Iterator it;
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 5.5, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosBegin(CoordinateType mz) const ))
{
  MSSpectrum::ConstIterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  MSSpectrum::ConstIterator it;
  it = tmp.PosBegin(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((Iterator PosEnd(CoordinateType mz)))
{
  MSSpectrum::Iterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosEnd(Iterator begin, CoordinateType mz, Iterator end)))
{
  MSSpectrum::Iterator it;
  it = tmp.PosEnd(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosEnd(tmp.begin(), 4.0, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosEnd(CoordinateType mz) const ))
{
  MSSpectrum::ConstIterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const ))
{
  MSSpectrum::ConstIterator it;
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 5.0, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it-1)->getPos(), (tmp.end()-1)->getPos())
}
END_SECTION

START_SECTION(bool containsIMData() const)
{
  auto ds = getPrefilledSpec();
  TEST_TRUE(ds.containsIMData())
}
END_SECTION

START_SECTION((std::pair<Size,DriftTimeUnit> getIMData() const))
{
  auto ds = getPrefilledSpec();
  auto [im_data_index, unit] = ds.getIMData();
  TEST_EQUAL(im_data_index, 1)
  TEST_TRUE(unit == DriftTimeUnit::MILLISECOND)
}
END_SECTION

const MSSpectrum spec_test {
  {412.321, 29.0f},
  {412.824, 60.0f},
  {413.8, 34.0f},
  {414.301, 29.0f},
  {415.287, 37.0f},
  {416.293, 31.0f},
  {418.232, 31.0f},
  {419.113, 31.0f},
  {420.13, 201.0f},
  {423.269, 56.0f},
  {426.292, 34.0f},
  {427.28, 82.0f},
  {428.322, 87.0f},
  {430.269, 30.0f},
  {431.246, 29.0f},
  {432.289, 42.0f},
  {436.161, 32.0f},
  {437.219, 54.0f},
  {439.186, 40.0f},
  {440.27, 40},
  {441.224, 23.0f}};

START_SECTION((Size findNearest(CoordinateType mz) const))
{
  MSSpectrum tmp = spec_test;

  //test outside mass range
  TEST_EQUAL(tmp.findNearest(400.0),0);
  TEST_EQUAL(tmp.findNearest(500.0),20);
  //test mass range borders
  TEST_EQUAL(tmp.findNearest(412.4),0);
  TEST_EQUAL(tmp.findNearest(441.224),20);
  //test inside scan
  TEST_EQUAL(tmp.findNearest(426.29),10);
  TEST_EQUAL(tmp.findNearest(426.3),10);
  TEST_EQUAL(tmp.findNearest(427.2),11);
  TEST_EQUAL(tmp.findNearest(427.3),11);

  //empty spectrum
  MSSpectrum tmp2;
  TEST_PRECONDITION_VIOLATED(tmp2.findNearest(427.3));
}
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz, CoordinateType tolerance) const))
{
  //test outside mass range
  TEST_EQUAL(spec_test.findNearest(400.0, 1.0), -1);
  TEST_EQUAL(spec_test.findNearest(500.0, 1.0), -1);

  //test mass range borders
  TEST_EQUAL(spec_test.findNearest(412.4, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(412.4, 0.1), 0);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.01),-1);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.1), 20);

  //test inside scan
  TEST_EQUAL(spec_test.findNearest(426.29, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(426.3, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(427.2, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001), -1);

  //empty spectrum
  MSSpectrum spec_test2;
  TEST_EQUAL(spec_test2.findNearest(427.3, 1.0, 1.0), -1);
}
END_SECTION
START_SECTION((Size findNearest(CoordinateType mz, CoordinateType left_tolerance, CoordinateType right_tolerance) const))
{
  //test outside mass range
  TEST_EQUAL(spec_test.findNearest(400.0, 1.0, 1.0), -1);
  TEST_EQUAL(spec_test.findNearest(500.0, 1.0, 1.0), -1);

  //test mass range borders
  TEST_EQUAL(spec_test.findNearest(412.4, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(412.4, 0.1, 0.1), 0);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.01, 0.01),-1);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.1, 0.1), 20);

  //test inside scan
  TEST_EQUAL(spec_test.findNearest(426.29, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(426.3, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(427.2, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 0.001), -1);

  TEST_EQUAL(spec_test.findNearest(427.3, 0.1, 0.001), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 1.01), -1);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 1.1), 12);

  //empty spectrum
  MSSpectrum spec_test2;
  TEST_EQUAL(spec_test2.findNearest(427.3, 1.0, 1.0), -1);
}
END_SECTION
START_SECTION((Size findHighestInWindow(CoordinateType mz, CoordinateType tolerance_left, CoordinateType tolerance_righ) const))
{
  //test outside mass range
  TEST_EQUAL(spec_test.findHighestInWindow(400.0, 1.0, 1.0), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(500.0, 1.0, 1.0), -1);

  //test mass range borders
  TEST_EQUAL(spec_test.findHighestInWindow(412.4, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(412.4, 0.1, 0.1), 0);
  TEST_EQUAL(spec_test.findHighestInWindow(441.3, 0.01, 0.01),-1);
  TEST_EQUAL(spec_test.findHighestInWindow(441.3, 0.1, 0.1), 20);

  //test inside scan
  TEST_EQUAL(spec_test.findHighestInWindow(426.29, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findHighestInWindow(426.3, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findHighestInWindow(427.2, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 0.001), -1);

  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.1, 0.001), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 1.01), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 1.1), 12);

  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 9.0, 4.0), 8);
  TEST_EQUAL(spec_test.findHighestInWindow(430.25, 1.9, 1.01), 13);

  //empty spectrum
  MSSpectrum spec_test2;
  TEST_EQUAL(spec_test2.findHighestInWindow(427.3, 1.0, 1.0), -1);
}
END_SECTION

START_SECTION( SpectrumSettings::SpectrumType MSSpectrum::getType(const bool query_data) const)
{
  // test empty spectrum
  MSSpectrum edit;
  TEST_EQUAL(edit.getType(false), SpectrumSettings::UNKNOWN);
  TEST_EQUAL(edit.getType(true), SpectrumSettings::UNKNOWN);

  // easiest: type is explicitly given
  edit.setType(SpectrumSettings::PROFILE);
  TEST_EQUAL(edit.getType(false), SpectrumSettings::PROFILE);
  TEST_EQUAL(edit.getType(true), SpectrumSettings::PROFILE);

  // second easiest: type is given in data processing
  DataProcessing dp;
  dp.setProcessingActions( { DataProcessing::PEAK_PICKING } );
  boost::shared_ptr< DataProcessing > dp_(new DataProcessing(dp));
  edit.getDataProcessing().push_back(dp_);
  // still profile, since DP is only checked when type is unknown
  TEST_EQUAL(edit.getType(false), SpectrumSettings::PROFILE);
  TEST_EQUAL(edit.getType(true), SpectrumSettings::PROFILE);
  edit.setType(SpectrumSettings::UNKNOWN);
  TEST_EQUAL(edit.getType(false), SpectrumSettings::CENTROID);
  TEST_EQUAL(edit.getType(true), SpectrumSettings::CENTROID);

  // third case: estimation from data
  edit.getDataProcessing().clear();
  // too few points
  edit.push_back( { 100.0, 1.0 } );
  edit.push_back( { 200.0, 1.0 } );
  edit.push_back( { 300.0, 1.0 } );
  edit.push_back( { 400.0, 1.0 } );
  TEST_EQUAL(edit.getType(false), SpectrumSettings::UNKNOWN);
  TEST_EQUAL(edit.getType(true), SpectrumSettings::UNKNOWN);
  edit.push_back( { 500.0, 1.0 } );
  edit.push_back( { 600.0, 1.0 } );
  TEST_EQUAL(edit.getType(false), SpectrumSettings::UNKNOWN); // data is not inspected
  TEST_EQUAL(edit.getType(true), SpectrumSettings::CENTROID);
}
END_SECTION


START_SECTION(ConstIterator getBasePeak() const)
{
  const auto it = spec_test.getBasePeak();
  TEST_REAL_SIMILAR(it->getIntensity(), 201.0)
  TEST_EQUAL(std::distance(spec_test.begin(), it), 8);
  MSSpectrum empty;
  TEST_EQUAL(empty.getBasePeak() == empty.end(), true);
}
END_SECTION


START_SECTION(Iterator getBasePeak())
{
  MSSpectrum test = spec_test;
  auto it = test.getBasePeak();
  it->setIntensity(it->getIntensity() + 0.0);
  TEST_REAL_SIMILAR(it->getIntensity(), 201.0)
  TEST_EQUAL(std::distance(test.begin(), it), 8);
}
END_SECTION


START_SECTION(PeakType::IntensityType calculateTIC() const)
{
  auto r = spec_test.calculateTIC();
  TEST_REAL_SIMILAR(r, 1032.0)
  TEST_EQUAL(MSSpectrum().calculateTIC(), 0.0);
}
END_SECTION


START_SECTION(void clear(bool clear_meta_data))
{
  MSSpectrum edit;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  edit.resize(1);
  edit.setMetaValue("label",String("bla"));
  edit.setRT(5);
  edit.setDriftTime(6);
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  edit.setMSLevel(5);
  edit.getFloatDataArrays().resize(5);
  edit.getIntegerDataArrays().resize(5);
  edit.getStringDataArrays().resize(5);

  edit.clear(false);
  TEST_EQUAL(edit.size(),0)
  TEST_EQUAL(edit == MSSpectrum(),false)
  TEST_EQUAL(edit.empty(),true)

  edit.clear(true);
  TEST_EQUAL(edit.empty(),true)
  TEST_EQUAL(edit == MSSpectrum(),true)
}
END_SECTION

START_SECTION(([MSSpectrum::RTLess] bool operator()(const MSSpectrum &a, const MSSpectrum &b) const))
{
  vector< MSSpectrum> v;

  MSSpectrum sp1;
  sp1.setRT(3.0f);
  v.push_back(sp1);

  MSSpectrum sp2;
  sp2.setRT(2.0f);
  v.push_back(sp2);

  MSSpectrum sp3;
  sp3.setRT(1.0f);
  v.push_back(sp3);

  std::sort(v.begin(),v.end(), MSSpectrum::RTLess());

  TEST_REAL_SIMILAR(v[0].getRT(), 1.0);
  TEST_REAL_SIMILAR(v[1].getRT(), 2.0);
  TEST_REAL_SIMILAR(v[2].getRT(), 3.0);

  ///
  MSSpectrum s1;
  s1.setRT(0.451);

  MSSpectrum s2;
  s2.setRT(0.5);

  TEST_EQUAL(MSSpectrum::RTLess()(s1,s2), true);
  TEST_EQUAL(MSSpectrum::RTLess()(s2,s1), false);
  TEST_EQUAL(MSSpectrum::RTLess()(s2,s2), false);
}
END_SECTION

START_SECTION((std::pair<DriftTimeUnit, std::vector<float>> maybeGetIMData() const))
{
  // Test successful retrieval of ion mobility data
  MSSpectrum spec;
  
  // Create a float data array with ion mobility data
  DataArrays::FloatDataArray im_array;
  im_array.setName("Ion Mobility");
  im_array.resize(3);
  im_array[0] = 1.0f;
  im_array[1] = 2.0f;
  im_array[2] = 3.0f;
  im_array.setMetaValue("unit", "millisecond");
  
  // Add the array to spectrum's float data arrays
  std::vector<DataArrays::FloatDataArray> fda;
  fda.push_back(im_array);
  spec.setFloatDataArrays(fda);
  
  // Test successful case
  std::pair<DriftTimeUnit, std::vector<float>> result = spec.maybeGetIMData();
  TEST_TRUE(result.first == DriftTimeUnit::MILLISECOND)
  TEST_EQUAL(result.second.size(), 3)
  TEST_REAL_SIMILAR(result.second[0], 1.0)
  TEST_REAL_SIMILAR(result.second[1], 2.0)
  TEST_REAL_SIMILAR(result.second[2], 3.0)
  
  // Test case with missing ion mobility data
  MSSpectrum spec_no_im;
  result = spec_no_im.maybeGetIMData();
  TEST_TRUE(result.first == DriftTimeUnit::NONE)
  TEST_EQUAL(result.second.empty(), true)
  
  // Test case with empty float arrays
  MSSpectrum spec_empty;
  spec_empty.getFloatDataArrays().clear();
  result = spec_empty.maybeGetIMData();
  TEST_TRUE(result.first == DriftTimeUnit::NONE)
  TEST_EQUAL(result.second.empty(), true)
  
  // Test case with wrong array name
  MSSpectrum spec_wrong_name;
  DataArrays::FloatDataArray wrong_array;
  wrong_array.setName("Wrong Name");
  wrong_array.resize(2);
  wrong_array[0] = 4.0f;
  wrong_array[1] = 5.0f;
  fda.clear();
  fda.push_back(wrong_array);
  spec_wrong_name.setFloatDataArrays(fda);
  result = spec_wrong_name.maybeGetIMData();
  TEST_TRUE(result.first == DriftTimeUnit::NONE)
  TEST_EQUAL(result.second.empty(), true)
}
END_SECTION

START_SECTION(([EXTRA] std::ostream& operator << (std::ostream& os, const MSSpectrum& spec)))
{
  MSSpectrum spec
  { {412.321, 29.0f}, //0
    {412.824, 60.0f}, //1
    {413.8, 34.0f},   //2
    {414.301, 29.0f}, //3
    {415.287, 37.0f}, //4
    {416.293, 31.0f}, //5
    {418.232, 31.0f}, //6
    {419.113, 31.0f}, //7
    {420.13, 201.0f}, //8
    {423.269, 56.0f}, //9
    {426.292, 34.0f}  //10
  };

  spec.getInstrumentSettings().getScanWindows().resize(1);
  spec.setMetaValue("label", 5.0);
  spec.setMSLevel(17);
  spec.setRT(7.0);
  spec.setName("bla");

  ostringstream test_stream;
  test_stream << spec;

  TEST_EQUAL(test_stream.str(), "-- MSSPECTRUM BEGIN --\n"
                                "-- SPECTRUMSETTINGS BEGIN --\n"
                                "-- SPECTRUMSETTINGS END --\n"
                                "POS: 412.321 INT: 29\n"
                                "POS: 412.824 INT: 60\n"
                                "POS: 413.8 INT: 34\n"
                                "POS: 414.301 INT: 29\n"
                                "POS: 415.287 INT: 37\n"
                                "POS: 416.293 INT: 31\n"
                                "POS: 418.232 INT: 31\n"
                                "POS: 419.113 INT: 31\n"
                                "POS: 420.13 INT: 201\n"
                                "POS: 423.269 INT: 56\n"
                                "POS: 426.292 INT: 34\n"
                                "-- MSSPECTRUM END --\n")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
