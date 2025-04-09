// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>


using namespace OpenMS;


class TestModel : public Fitter1D
{
  public:
  TestModel()
    : Fitter1D()
  {
    setName("TestModel");

    check_defaults_ = false;

    defaultsToParam_();
  }

  TestModel(const TestModel& source) : Fitter1D(source)
  {
    updateMembers_();
  }

  ~TestModel() override
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    Fitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_() override
  {
    Fitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, std::unique_ptr<InterpolationModel>&  /*model*/) override
  {
//    double center = 0.0;
//    center = model->getCenter();

    return 1.0;
  }

};


START_TEST(Fitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using std::stringstream;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
START_SECTION(Fitter1D())
{
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~Fitter1D()))
  delete ptr;
END_SECTION

START_SECTION((virtual QualityType fit1d(const  RawDataArrayType &, InterpolationModel *&)))
	Fitter1D f1d;
  Fitter1D::RawDataArrayType rft;
  std::unique_ptr<InterpolationModel> ipm;
	TEST_EXCEPTION(Exception::NotImplemented,f1d.fit1d(rft,ipm));
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



