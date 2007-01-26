// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(InterpolationModel , "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

class TestModel : public InterpolationModel< >
{
  public:
	TestModel()
		: InterpolationModel< >()
	{
		setName(getProductName());
		
		check_defaults_ = false;
		
		defaultsToParam_();
	}


	TestModel(const TestModel& source)
		: InterpolationModel< >(source)
	{
		updateMembers_();
	}
	
	virtual ~TestModel()
	{
	}
	
	virtual TestModel& operator = (const TestModel& source)
	{
		if (&source == this) return *this;
		
		InterpolationModel< >::operator = (source);
		updateMembers_();
		
		return *this;
	}
	
	void updateMembers_()
	{
		 InterpolationModel< >::updateMembers_();
	}

	IntensityType getIntensity(const PositionType& pos) const
	{
		return pos[0]*3.0;
	}

	bool isContained(const PositionType& pos) const
	{
		return getIntensity(pos)>cut_off_;
	}

	void  fillIntensity(PeakType& peak) const
	{
		peak.setIntensity(getIntensity(peak.getPosition()));
	}

	void getSamples(SamplesType& /*cont*/) const
	{
	}
	
	void setSamples() 
	{
	}
	
	const CoordinateType getCenter() const
	{
		return 10.0;
	}

	static const String getProductName()
	{ 
		return "TestModel"; 
	}



};


// default ctor
TestModel* ptr = 0;
CHECK(TestModel())
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~TestModel())
	delete ptr;
RESULT

// assignment operator
CHECK(TestModel& operator = (const TestModel& source))
	TestModel tm1;
  TestModel tm2;
  
  tm1.setCutOff(3.3);
  tm2 = tm1;
	TEST_REAL_EQUAL(tm1.getCutOff(),tm2.getCutOff())
	TEST_REAL_EQUAL(tm1.getScalingFactor(),tm2.getScalingFactor())
RESULT

// copy constructor
CHECK(TestModel(const TestModel& source))
	TestModel fp1;	
  fp1.setCutOff(0.1);

  TestModel fp2(fp1);

  TestModel fp3;
  fp3.setCutOff(0.1);

  fp1 = TestModel();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK(IntensityType getCutOff() const)
  const TestModel s;
  TEST_REAL_EQUAL(s.getCutOff(), TestModel::IntensityType(0))
RESULT

CHECK(void setCutOff(IntensityType cut_off))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_EQUAL(s.getCutOff(), 4.4)
RESULT

CHECK(const String& getName() const)
	TestModel s;
  TEST_EQUAL(s.getName(), "TestModel")
RESULT

CHECK(IntensityType getIntensity(const PositionType& pos) const)
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  TEST_REAL_EQUAL(s.getIntensity(pos), 0.3)
RESULT

CHECK(bool isContained(const PositionType& pos) const)
	TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  const TestModel& t = s;
  TEST_REAL_EQUAL(t.isContained(pos), false)
RESULT

CHECK(void fillIntensity(PeakType& peak) const)
	const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.getIntensity() = 0.1;
  t.fillIntensity(p);
  TEST_REAL_EQUAL(p.getIntensity(), 0.3)
RESULT

CHECK(void  fillIntensities(PeakIterator beg, PeakIterator end) const)
	const TestModel t;
  std::vector< TestModel::PeakType > vec(4);
  for (UnsignedInt i=0; i<4; ++i)
  {
		vec[i].getIntensity() = -0.5;
		vec[i].getPosition()[0] = i;
	}
  t.fillIntensities(vec.begin()+1, vec.end()-1);
  TEST_EQUAL(vec[0].getIntensity(), -0.5)
  TEST_EQUAL(vec[1].getIntensity(), 3.0)
  TEST_EQUAL(vec[2].getIntensity(), 6.0)
  TEST_EQUAL(vec[3].getIntensity(), -0.5)
RESULT

CHECK(void  getCenter() const)
	const TestModel t;
 TEST_REAL_EQUAL(t.getCenter(),10.0);
RESULT

CHECK([EXTRA] DefaultParmHandler::setParameters(...))
	Param p;
	p.setValue("cutoff",17.0);
	TestModel m;
	m.setParameters(p);
	TEST_REAL_EQUAL(m.getParameters().getValue("cutoff"), 17.0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
