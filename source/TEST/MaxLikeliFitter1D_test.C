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

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>
///////////////////////////

///////////////////////////

using namespace OpenMS;
using namespace std;

class TestModel : public MaxLikeliFitter1D
{
  public: TestModel() : MaxLikeliFitter1D()
  {
    setName(getProductName());
    check_defaults_ = false;
    defaultsToParam_();
  }


  TestModel(const TestModel& source) : MaxLikeliFitter1D(source)
  {
    updateMembers_();
  }

  virtual ~TestModel()
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    MaxLikeliFitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_()
  {
     MaxLikeliFitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, InterpolationModel*& model)
  {
    DoubleReal center = 0.0;
    center = model->getCenter();

    return 1.0;
  }

  QualityType fitOffset_(InterpolationModel* model, const RawDataArrayType& /*set*/ , const CoordinateType stdev1, const CoordinateType stdev2, const CoordinateType offset_step)
  {
    DoubleReal center = 0.0;
    center = model->getCenter();

    DoubleReal st_dev_1 = 0.0;
    st_dev_1 = stdev1;
    DoubleReal st_dev_2 = 0.0;
    st_dev_2 = stdev2;
    DoubleReal offset = 0.0;
    offset = offset_step;

    return 1.0;
  }

  static const String getProductName()
  {
    return "TestModel";
  }

};

/////////////////////////////////////////////////////////////

START_TEST(MaxLikeliFitter1D, "$Id$")

/////////////////////////////////////////////////////////////

TestModel* ptr = 0;
TestModel* nullPointer = 0;
START_SECTION(MaxLikeliFitter1D())
{
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((MaxLikeliFitter1D(const  MaxLikeliFitter1D &source)))
	TestModel tm1;

  TestModel tm2(tm1);
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION

START_SECTION((virtual ~MaxLikeliFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual MaxLikeliFitter1D& operator=(const  MaxLikeliFitter1D &source)))
	TestModel tm1;
  TestModel tm2;

  tm2 = tm1;
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



