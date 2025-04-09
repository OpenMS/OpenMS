// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/MaxLikeliFitter1D.h>
///////////////////////////

///////////////////////////

using namespace OpenMS;
using namespace std;

class TestModel : public MaxLikeliFitter1D
{
  public: TestModel() : MaxLikeliFitter1D()
  {
    setName("TestModel");
    check_defaults_ = false;
    defaultsToParam_();
  }


  TestModel(const TestModel& source) : MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  ~TestModel() override
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    MaxLikeliFitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_() override
  {
     MaxLikeliFitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, std::unique_ptr<InterpolationModel>&  /*model*/) override
  {
//    double center = 0.0;
//    center = model->getCenter();

    return 1.0;
  }

  QualityType fitOffset_(InterpolationModel* /* model */, const RawDataArrayType& /*set*/ , const CoordinateType /* stdev1 */, const CoordinateType /* stdev2 */, const CoordinateType /* offset_step */)
  {
//    double center = 0.0;
//    center = model->getCenter();

//    double st_dev_1 = 0.0;
//    st_dev_1 = stdev1;
//    double st_dev_2 = 0.0;
//    st_dev_2 = stdev2;
//    double offset = 0.0;
//    offset = offset_step;

    return 1.0;
  }

};

/////////////////////////////////////////////////////////////

START_TEST(MaxLikeliFitter1D, "$Id$")

/////////////////////////////////////////////////////////////

TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
START_SECTION(MaxLikeliFitter1D())
{
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((MaxLikeliFitter1D(const  MaxLikeliFitter1D &source)))
	TestModel tm1;
  TestModel tm2(tm1);
END_SECTION

START_SECTION((virtual ~MaxLikeliFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual MaxLikeliFitter1D& operator=(const  MaxLikeliFitter1D &source)))
	TestModel tm1;
  TestModel tm2;
  tm2 = tm1;
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



