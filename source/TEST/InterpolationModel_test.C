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

class TestModel : public InterpolationModel
{
  public:
	TestModel()
		: InterpolationModel()
	{
		setName(getProductName());
		
		check_defaults_ = false;
		
		defaultsToParam_();
	}


	TestModel(const TestModel& source)
		: InterpolationModel(source)
	{
		updateMembers_();
	}
	
	virtual ~TestModel()
	{
	}
	
	virtual TestModel& operator = (const TestModel& source)
	{
		if (&source == this) return *this;
		
		InterpolationModel::operator = (source);
		updateMembers_();
		
		return *this;
	}
	
	void updateMembers_()
	{
		 InterpolationModel::updateMembers_();
	}

	IntensityType getIntensity(const PositionType& pos) const
	{
		return pos[0]*3.0;
	}
	
	IntensityType getIntensity(CoordinateType coord) const
	{
		return coord*3.0;
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
	
	CoordinateType getCenter() const
	{
		return 10.0;
	}

	static const String getProductName()
	{ 
		return "TestModel"; 
	}



};

//////////////////////////////////////

// default ctor
TestModel* ptr = 0;
CHECK((InterpolationModel()))
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~InterpolationModel()))
	delete ptr;
RESULT

// assignment operator
CHECK((virtual InterpolationModel& operator=(const InterpolationModel &source)))
	TestModel tm1;
  TestModel tm2;
  
  tm1.setCutOff(3.3);
  tm2 = tm1;
	TEST_REAL_EQUAL(tm1.getCutOff(),tm2.getCutOff())
	TEST_REAL_EQUAL(tm1.getScalingFactor(),tm2.getScalingFactor())
RESULT

// copy constructor
CHECK((InterpolationModel(const InterpolationModel &source)))
	TestModel fp1;	
  fp1.setCutOff(0.1);

  TestModel fp2(fp1);

  TestModel fp3;
  fp3.setCutOff(0.1);

  fp1 = TestModel();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK(([EXTRA]IntensityType getCutOff() const))
  const TestModel s;
  TEST_REAL_EQUAL(s.getCutOff(), TestModel::IntensityType(0))
RESULT

CHECK(([EXTRA]void setCutOff(IntensityType cut_off)))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_EQUAL(s.getCutOff(), 4.4)
RESULT

CHECK(([EXTRA]const String& getName() const))
	TestModel s;
  TEST_EQUAL(s.getName(), "TestModel")
RESULT

CHECK((IntensityType getIntensity(const PositionType& pos) const))
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  TEST_REAL_EQUAL(s.getIntensity(pos), 0.3)
RESULT

CHECK(([EXTRA]bool isContained(const PositionType& pos) const))
	TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  const TestModel& t = s;
  TEST_REAL_EQUAL(t.isContained(pos), false)
RESULT

CHECK(([EXTRA]void fillIntensity(PeakType& peak) const))
	const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.setIntensity(0.1);
  t.fillIntensity(p);
  TEST_REAL_EQUAL(p.getIntensity(), 0.3)
RESULT

CHECK(([EXTRA]void  fillIntensities(PeakIterator beg, PeakIterator end) const))
	const TestModel t;
  std::vector< TestModel::PeakType > vec(4);
  for (UInt i=0; i<4; ++i)
  {
		vec[i].setIntensity(-0.5);
		vec[i].getPosition()[0] = i;
	}
  t.fillIntensities(vec.begin()+1, vec.end()-1);
  TEST_EQUAL(vec[0].getIntensity(), -0.5)
  TEST_EQUAL(vec[1].getIntensity(), 3.0)
  TEST_EQUAL(vec[2].getIntensity(), 6.0)
  TEST_EQUAL(vec[3].getIntensity(), -0.5)
RESULT

CHECK( virtual CoordinateType getCenter() const) 
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

CHECK( void setScalingFactor(CoordinateType scaling) )
	TestModel tm;
	tm.setScalingFactor(2.0);
	
	TEST_REAL_EQUAL(tm.getParameters().getValue("intensity_scaling"),2.0)
	TEST_REAL_EQUAL(tm.getScalingFactor(),2.0)		
RESULT

CHECK( void setInterpolationStep(CoordinateType interpolation_step) )
	TestModel tm;
	tm.setInterpolationStep( 10.5 );
	
	TEST_REAL_EQUAL(tm.getParameters().getValue("interpolation_step"), 10.5 )
RESULT

CHECK( virtual void setSamples() )
	// not much to be tested here
RESULT

CHECK( void getSamples(SamplesType &cont) const )
	// not much to be tested here
RESULT

CHECK( virtual void setOffset(CoordinateType offset) )

RESULT

CHECK( CoordinateType getScalingFactor() const )
	TestModel tm;
	tm.setScalingFactor(666.0);
	
	TEST_REAL_EQUAL(tm.getParameters().getValue("intensity_scaling"),666.0)
	TEST_REAL_EQUAL(tm.getScalingFactor(),666.0)		
RESULT

CHECK( const LinearInterpolation& getInterpolation() const )
	TestModel tm;
	InterpolationModel::LinearInterpolation interpol1;
	InterpolationModel::LinearInterpolation interpol2 = tm.getInterpolation();
	
	// compare models
	TEST_REAL_EQUAL(interpol1.getScale(), interpol2.getScale());
	TEST_REAL_EQUAL(interpol1.getInsideReferencePoint(), interpol2.getInsideReferencePoint());
	TEST_REAL_EQUAL(interpol1.getOutsideReferencePoint(), interpol2.getOutsideReferencePoint() );
	
RESULT

CHECK( IntensityType getIntensity(CoordinateType coord) const )
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  TEST_REAL_EQUAL(s.getIntensity(pos), 0.3)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
