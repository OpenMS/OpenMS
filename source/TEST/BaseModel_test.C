// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/CONCEPT/Exception.h>

///////////////////////////

using namespace OpenMS;
using std::stringstream;

class TestModel : public BaseModel<2>
{
  public:
	TestModel()
		: BaseModel<2>()
	{
		setName(getProductName());

		check_defaults_ = false;

		defaultsToParam_();
	}

	TestModel(const TestModel& source)
		: BaseModel<2>(source)
	{
		updateMembers_();
	}

	virtual ~TestModel()
	{
	}

	virtual TestModel& operator = (const TestModel& source)
	{
		if (&source == this) return *this;

		BaseModel<2>::operator = (source);
		updateMembers_();

		return *this;
	}

	void updateMembers_()
	{
		BaseModel<2>::updateMembers_();
	}

	IntensityType getIntensity(const PositionType& pos) const
	{
		return pos[0]+pos[1];
	}

	bool isContained(const PositionType& pos) const
	{
		return getIntensity(pos)>cut_off_;
	}

	void getSamples(SamplesType& /*cont*/) const
	{
	}

	static const String getProductName()
	{
		return "TestModel";
	}

};

START_TEST(BaseModel, "$Id$")

// default ctor
TestModel* ptr = 0;
TestModel* nullPointer = 0;
START_SECTION((BaseModel()))
	ptr = new TestModel();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~BaseModel()))
	delete ptr;
END_SECTION

// assignment operator
START_SECTION((virtual BaseModel& operator=(const BaseModel &source)))
	TestModel tm1;
  TestModel tm2;

  tm1.setCutOff(3.3);
  tm2 = tm1;
	TEST_REAL_SIMILAR(tm1.getCutOff(),tm2.getCutOff())
END_SECTION

// copy constructor
START_SECTION((BaseModel(const BaseModel &source)))
	TestModel tm1;
  tm1.setCutOff(0.1);

  TestModel tm2(tm1);
	TEST_REAL_SIMILAR(tm1.getCutOff(),tm2.getCutOff())
END_SECTION

START_SECTION(([EXTRA]IntensityType getCutOff() const))
  const TestModel s;
  TEST_REAL_SIMILAR(s.getCutOff(), TestModel::IntensityType(0))
END_SECTION

START_SECTION((virtual void setCutOff(IntensityType cut_off)))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_SIMILAR(s.getCutOff(), 4.4)
END_SECTION

START_SECTION(([EXTRA]const String& getName() const))
	TestModel s;
  TEST_EQUAL(s.getName(), "TestModel")
END_SECTION

START_SECTION((virtual IntensityType getIntensity(const PositionType &pos) const =0))
{
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  pos[1]=0.2;
	TEST_REAL_SIMILAR(s.getIntensity(pos), 0.3);
}
END_SECTION

START_SECTION((virtual bool isContained(const PositionType &pos) const))
	TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  pos[1]=0.2;
  const TestModel& t = s;
  TEST_EQUAL(t.isContained(pos), false)
END_SECTION

START_SECTION((template <typename PeakType> void fillIntensity(PeakType &peak) const))
	const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.getPosition()[1]=0.2;
  p.setIntensity(0.1f);
  t.fillIntensity(p);
  TEST_REAL_SIMILAR(p.getIntensity(), 0.3)
END_SECTION

START_SECTION((template <class PeakIterator> void fillIntensities(PeakIterator begin, PeakIterator end) const))
	const TestModel t;
  std::vector< TestModel::PeakType > vec(4);
  for (Size i=0; i<4; ++i)
  {
		vec[i].setIntensity(-0.5);
		vec[i].getPosition()[0] = i;
	}
  t.fillIntensities(vec.begin()+1, vec.end()-1);
  TEST_EQUAL(vec[0].getIntensity(), -0.5)
  TEST_EQUAL(vec[1].getIntensity(), 1.0)
  TEST_EQUAL(vec[2].getIntensity(), 2.0)
  TEST_EQUAL(vec[3].getIntensity(), -0.5)
END_SECTION

START_SECTION([EXTRA] DefaultParmHandler::setParameters(...))
	Param p;
	p.setValue("cutoff",17.0);
	TestModel m;
	m.setParameters(p);
	TEST_REAL_SIMILAR(m.getParameters().getValue("cutoff"), 17.0)
END_SECTION

START_SECTION((static void registerChildren()))
	// TODO
END_SECTION

START_SECTION((virtual IntensityType getCutOff() const))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_SIMILAR(s.getCutOff(), 4.4)
END_SECTION

START_SECTION((virtual void getSamples(SamplesType &cont) const =0))
  NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void getSamples(std::ostream &os)))
  NOT_TESTABLE;
END_SECTION

START_SECTION((template <class PeakIterator> void registerChildren()))
{
  NOT_TESTABLE;
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
