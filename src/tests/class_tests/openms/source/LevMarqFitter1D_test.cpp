// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/FEATUREFINDER/Fitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

///////////////////////////


class TestModel : public LevMarqFitter1D
{
  public: TestModel() : LevMarqFitter1D()
  {
    setName("TestModel");
    check_defaults_ = false;
    defaultsToParam_();
  }


  TestModel(const TestModel& source) : LevMarqFitter1D(source)
  {
    updateMembers_();
  }

  ~TestModel() override
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    LevMarqFitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_() override
  {
     LevMarqFitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, std::unique_ptr<InterpolationModel>&  /*model*/) override
  {
//    double center = 0.0;
//    center = model->getCenter();

    return 1.0;
  }

  void optimize_()
  {
  }
};

/////////////////////////////////////////////////////////////

START_TEST(LevMarqFitter1D, "$Id$")

///////////////////////////


/////////////////////////////////////////////////////////////


TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
START_SECTION((LevMarqFitter1D()))
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((LevMarqFitter1D(const  LevMarqFitter1D &source)))
	TestModel tm1;
  TestModel tm2(tm1);
END_SECTION

START_SECTION((virtual ~LevMarqFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual LevMarqFitter1D& operator=(const  LevMarqFitter1D &source)))
	TestModel tm1;
  TestModel tm2;
  tm2 = tm1;
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



