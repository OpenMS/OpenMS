// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>


using namespace OpenMS;


class TestModel : public Fitter1D
{
  public:
  TestModel()
    : Fitter1D()
  {
    setName(getProductName());

    check_defaults_ = false;

    defaultsToParam_();
  }

  TestModel(const TestModel& source) : Fitter1D(source)
  {
    updateMembers_();
  }

  virtual ~TestModel()
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    Fitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_()
  {
    Fitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, InterpolationModel*& model)
  {
    DoubleReal center = 0.0;
    center = model->getCenter();

    return 1.0;
  }

  static const String getProductName()
  {
    return "TestModel";
  }

};


START_TEST(Fitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using std::stringstream;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


TestModel* ptr = 0;
TestModel* nullPointer = 0;
START_SECTION(Fitter1D())
{
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((Fitter1D(const  Fitter1D &source)))
	TestModel tm1;

  TestModel tm2(tm1);
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION

START_SECTION((virtual ~Fitter1D()))
  delete ptr;
END_SECTION

START_SECTION((virtual Fitter1D& operator=(const  Fitter1D &source)))
	TestModel tm1;
  TestModel tm2;

  tm2 = tm1;
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION

START_SECTION((virtual QualityType fit1d(const  RawDataArrayType &, InterpolationModel *&)))
	Fitter1D f1d;
  Fitter1D::RawDataArrayType rft;
  InterpolationModel *ipm = 0;
	TEST_EXCEPTION(Exception::NotImplemented,f1d.fit1d(rft,ipm));
END_SECTION

START_SECTION((void registerChildren()))
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



