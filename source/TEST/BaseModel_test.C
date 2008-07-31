// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(BaseModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

class TestModel : public BaseModel<3>
{
  public:
	TestModel()
		: BaseModel<3>()
	{
		setName(getProductName());
		
		check_defaults_ = false;
		
		defaultsToParam_();
	}

	TestModel(const TestModel& source)
		: BaseModel<3>(source)
	{
		updateMembers_();
	}
	
	virtual ~TestModel()
	{
	}
	
	virtual TestModel& operator = (const TestModel& source)
	{
		if (&source == this) return *this;
		
		BaseModel<3>::operator = (source);
		updateMembers_();
		
		return *this;
	}
	
	void updateMembers_()
	{
		BaseModel<3>::updateMembers_();
	}

	IntensityType getIntensity(const PositionType& pos) const
	{
		return pos[0]+pos[1]+pos[2];
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

	static const String getProductName()
	{ 
		return "TestModel"; 
	}

};


// default ctor
TestModel* ptr = 0;
CHECK((BaseModel()))
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~BaseModel()))
	delete ptr;
RESULT

// assignment operator
CHECK((virtual BaseModel& operator=(const BaseModel &source)))
	TestModel tm1;
  TestModel tm2;
  
  tm1.setCutOff(3.3);
  tm2 = tm1;
	TEST_REAL_EQUAL(tm1.getCutOff(),tm2.getCutOff())
RESULT

// copy constructor
CHECK((BaseModel(const BaseModel &source)))
	TestModel tm1;	
  tm1.setCutOff(0.1);
	
  TestModel tm2(tm1);
	TEST_REAL_EQUAL(tm1.getCutOff(),tm2.getCutOff())
RESULT

CHECK(([EXTRA]IntensityType getCutOff() const))
  const TestModel s;
  TEST_REAL_EQUAL(s.getCutOff(), TestModel::IntensityType(0))
RESULT

CHECK((virtual void setCutOff(IntensityType cut_off)))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_EQUAL(s.getCutOff(), 4.4)
RESULT

CHECK(([EXTRA]const String& getName() const))
	TestModel s;
  TEST_EQUAL(s.getName(), "TestModel")
RESULT

CHECK((virtual IntensityType getIntensity(const PositionType &pos) const =0))
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  pos[1]=0.2;
  pos[2]=0.3;
  TEST_REAL_EQUAL(s.getIntensity(pos), 0.6)
RESULT

CHECK((virtual bool isContained(const PositionType &pos) const ))
	TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  pos[1]=0.2;
  pos[2]=0.3;
  const TestModel& t = s;
  TEST_REAL_EQUAL(t.isContained(pos), false)
RESULT

CHECK((virtual void fillIntensity(PeakType &peak) const ))
	const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.getPosition()[1]=0.2;
  p.getPosition()[2]=0.3;
  p.setIntensity(0.1);
  t.fillIntensity(p);
  TEST_REAL_EQUAL(p.getIntensity(), 0.6)
RESULT

CHECK((template <class PeakIterator> void fillIntensities(PeakIterator beg, PeakIterator end) const ))
	const TestModel t;
  std::vector< TestModel::PeakType > vec(4);
  for (UInt i=0; i<4; ++i)
  {
		vec[i].setIntensity(-0.5);
		vec[i].getPosition()[0] = i;
	}
  t.fillIntensities(vec.begin()+1, vec.end()-1);
  TEST_EQUAL(vec[0].getIntensity(), -0.5)
  TEST_EQUAL(vec[1].getIntensity(), 1.0)
  TEST_EQUAL(vec[2].getIntensity(), 2.0)
  TEST_EQUAL(vec[3].getIntensity(), -0.5)
RESULT
	
CHECK([EXTRA] DefaultParmHandler::setParameters(...))
	Param p;
	p.setValue("cutoff",17.0);
	TestModel m;
	m.setParameters(p);
	TEST_REAL_EQUAL(m.getParameters().getValue("cutoff"), 17.0)
RESULT

CHECK((static void registerChildren()))
	// TODO
RESULT

CHECK( virtual IntensityType getCutOff() const )
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_EQUAL(s.getCutOff(), 4.4)
RESULT

CHECK( virtual void getSamples(SamplesType &cont) const =0 )
	// TODO
RESULT

CHECK( virtual void getSamples(std::ostream &os) )
	// TODO
RESULT

CHECK((template <class PeakIterator> void registerChildren()))
{
	// TODO
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
