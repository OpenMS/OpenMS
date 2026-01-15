// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/CalibrationData.h>
///////////////////////////

#include <OpenMS/MATH/MathFunctions.h>

using namespace OpenMS;
using namespace std;

START_TEST(Adduct, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CalibrationData* ptr = nullptr;
CalibrationData* nullPointer = nullptr;
START_SECTION(CalibrationData())
{
	ptr = new CalibrationData();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CalibrationData())
{
	delete ptr;
}
END_SECTION


CalibrationData cd;
for (Size i = 0; i < 10; ++i)
{
  cd.insertCalibrationPoint(100.100 + i, 200.200 + i, 128.5 + i, 200.0 + i, 1, 66);
  cd.insertCalibrationPoint(120.100 + i + 0.5, 400.200 + i, 128.5 + i, 200.0 + i, 1, 77);
}

START_SECTION(CalDataType::CoordinateType getMZ(Size i) const)
  TEST_REAL_SIMILAR(cd.getMZ(0), 200.200 + 0)
  TEST_REAL_SIMILAR(cd.getMZ(3), 400.200 + 1)
END_SECTION

START_SECTION(CalDataType::CoordinateType getRT(Size i) const)
  TEST_REAL_SIMILAR(cd.getRT(0), 100.100 + 0)
  TEST_REAL_SIMILAR(cd.getRT(3), 120.100 + 1 + 0.5)
END_SECTION
      
START_SECTION(CalDataType::CoordinateType getIntensity(Size i) const)
  TEST_REAL_SIMILAR(cd.getIntensity(0), 128.5 + 0)
  TEST_REAL_SIMILAR(cd.getIntensity(3), 128.5 + 1)
END_SECTION

START_SECTION(const_iterator begin() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION
      
START_SECTION(const_iterator end() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION

START_SECTION(Size size() const)
  TEST_EQUAL(cd.size(), Size(cd.end() - cd.begin()))
END_SECTION


START_SECTION(bool empty() const)
  TEST_EQUAL(cd.empty(), false)
  TEST_EQUAL(CalibrationData().empty(), true)
END_SECTION

START_SECTION(void clear())
  CalibrationData cd2(cd);
  TEST_EQUAL(cd2.empty(), false)
  cd2.clear();
  TEST_EQUAL(cd2.empty(), true)
END_SECTION
      
START_SECTION(void setUsePPM(bool usePPM))
  cd.setUsePPM(false);
  TEST_EQUAL(cd.usePPM(), false)
  cd.setUsePPM(true);
  TEST_EQUAL(cd.usePPM(), true)
END_SECTION

      
START_SECTION(bool usePPM() const)
  NOT_TESTABLE // tested above
END_SECTION

      
START_SECTION(void insertCalibrationPoint(CalDataType::CoordinateType rt, CalDataType::CoordinateType mz_obs, CalDataType::IntensityType intensity, CalDataType::CoordinateType mz_ref, double weight, int group = -1))
  NOT_TESTABLE // tested above
END_SECTION

      
START_SECTION(Size getNrOfGroups() const)
  TEST_EQUAL(cd.getNrOfGroups(), 2);
  CalibrationData cd2;
  TEST_EQUAL(cd2.getNrOfGroups(), 0);
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getError(Size i) const)
  TEST_REAL_SIMILAR(cd.getError(0), Math::getPPM(200.200 + 0, 200.0))
  TEST_REAL_SIMILAR(cd.getError(3), Math::getPPM(400.200 + 1, 200.0 + 1))
  cd.setUsePPM(false); // use absolute error
  TEST_REAL_SIMILAR(cd.getError(0), (200.200 + 0 - 200.0))
  TEST_REAL_SIMILAR(cd.getError(3), (400.200 + 1 - (200.0 + 1)))
  cd.setUsePPM(true); // reset
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getRefMZ(Size i) const)
  TEST_REAL_SIMILAR(cd.getRefMZ(0), 200.0 + 0)
  TEST_REAL_SIMILAR(cd.getRefMZ(3), 200.0 + 1)
END_SECTION

      
START_SECTION(CalDataType::CoordinateType getWeight(Size i) const)
  TEST_REAL_SIMILAR(cd.getWeight(0), 1.0)
  TEST_REAL_SIMILAR(cd.getWeight(3), 1.0)
END_SECTION

      
START_SECTION(int getGroup(Size i) const)
  TEST_EQUAL(cd.getGroup(0), 66)
  TEST_EQUAL(cd.getGroup(3), 77)
END_SECTION

      
START_SECTION(static StringList getMetaValues())
  TEST_EQUAL(ListUtils::concatenate(CalibrationData::getMetaValues(), ","), "mz_ref,ppm_error,weight")
END_SECTION

      
START_SECTION(CalibrationData median(double rt_left, double rt_right) const)
  CalibrationData m = cd.median(0, 1e6);
  TEST_EQUAL(m.size(), 2); // two medians (of two groups)
  TEST_REAL_SIMILAR(m.getMZ(0),  200.200 + 9.0/2)
  TEST_REAL_SIMILAR(m.getMZ(1),  400.200 + 9.0/2)
END_SECTION

      
START_SECTION(void sortByRT())
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



