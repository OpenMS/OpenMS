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

#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

using namespace OpenMS;

class TestModel : public InterpolationModel
{
  public:
  TestModel()
    : InterpolationModel()
  {
    setName("TestModel");

    check_defaults_ = false;

    defaultsToParam_();
  }


  TestModel(const TestModel& source)
    : InterpolationModel(source)
  {
    updateMembers_();
  }

  ~TestModel() override
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    InterpolationModel::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_() override
  {
     InterpolationModel::updateMembers_();
  }

  IntensityType getIntensity(const PositionType& pos) const override
  {
    return pos[0]*3.0;
  }

  IntensityType getIntensity(CoordinateType coord) const
  {
    return coord*3.0;
  }

  bool isContained(const PositionType& pos) const override
  {
    return getIntensity(pos)>cut_off_;
  }

  void  fillIntensity(PeakType& peak) const
  {
    peak.setIntensity(getIntensity(peak.getPosition()));
  }

  void getSamples(SamplesType& /*cont*/) const override
  {
  }

  void setSamples() override
  {
  }

  CoordinateType getCenter() const override
  {
    return 10.0;
  }

};


START_TEST(InterpolationModel , "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using std::stringstream;

//////////////////////////////////////

// default ctor
TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
START_SECTION((InterpolationModel()))
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~InterpolationModel()))
	delete ptr;
END_SECTION

// assignment operator
START_SECTION((virtual InterpolationModel& operator=(const InterpolationModel &source)))
	TestModel tm1;
  TestModel tm2;

  tm1.setCutOff(3.3);
  tm2 = tm1;
	TEST_REAL_SIMILAR(tm1.getCutOff(),tm2.getCutOff())
	TEST_REAL_SIMILAR(tm1.getScalingFactor(),tm2.getScalingFactor())
END_SECTION

// copy constructor
START_SECTION((InterpolationModel(const InterpolationModel &source)))
	TestModel fp1;
  fp1.setCutOff(0.1);

  TestModel fp2(fp1);

  TestModel fp3;
  fp3.setCutOff(0.1);

  fp1 = TestModel();
	TEST_TRUE(fp2 == fp3)
END_SECTION

START_SECTION(([EXTRA]IntensityType getCutOff() const))
  const TestModel s;
  TEST_REAL_SIMILAR(s.getCutOff(), TestModel::IntensityType(0))
END_SECTION

START_SECTION(([EXTRA]void setCutOff(IntensityType cut_off)))
	TestModel s;
	s.setCutOff(4.4);
  TEST_REAL_SIMILAR(s.getCutOff(), 4.4)
END_SECTION

START_SECTION(([EXTRA]const String& getName() const))
	TestModel s;
  TEST_EQUAL(s.getName(), "TestModel")
END_SECTION

START_SECTION((IntensityType getIntensity(const PositionType& pos) const))
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  TEST_REAL_SIMILAR(s.getIntensity(pos), 0.3)
END_SECTION

START_SECTION(([EXTRA]bool isContained(const PositionType& pos) const))
	TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  const TestModel& t = s;
  TEST_EQUAL(t.isContained(pos), false)
END_SECTION

START_SECTION(([EXTRA]void fillIntensity(PeakType& peak) const))
	const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.setIntensity(0.1f);
  t.fillIntensity(p);
  TEST_REAL_SIMILAR(p.getIntensity(), 0.3)
END_SECTION

START_SECTION(([EXTRA]void  fillIntensities(PeakIterator beg, PeakIterator end) const))
	const TestModel t;
  std::vector< TestModel::PeakType > vec(4);
  for (Size i=0; i<4; ++i)
  {
		vec[i].setIntensity(-0.5);
		vec[i].getPosition()[0] = i;
	}
  t.fillIntensities(vec.begin()+1, vec.end()-1);
  TEST_EQUAL(vec[0].getIntensity(), -0.5)
  TEST_EQUAL(vec[1].getIntensity(), 3.0)
  TEST_EQUAL(vec[2].getIntensity(), 6.0)
  TEST_EQUAL(vec[3].getIntensity(), -0.5)
END_SECTION

START_SECTION( virtual CoordinateType getCenter() const)
	const TestModel t;
 TEST_REAL_SIMILAR(t.getCenter(),10.0);
END_SECTION

START_SECTION([EXTRA] DefaultParmHandler::setParameters(...))
	Param p;
	p.setValue("cutoff",17.0);
	TestModel m;
	m.setParameters(p);
	TEST_REAL_SIMILAR(m.getParameters().getValue("cutoff"), 17.0)
END_SECTION

START_SECTION( void setScalingFactor(CoordinateType scaling) )
	TestModel tm;
	tm.setScalingFactor(2.0);

	TEST_REAL_SIMILAR(tm.getParameters().getValue("intensity_scaling"),2.0)
	TEST_REAL_SIMILAR(tm.getScalingFactor(),2.0)
END_SECTION

START_SECTION( void setInterpolationStep(CoordinateType interpolation_step) )
	TestModel tm;
	tm.setInterpolationStep( 10.5 );

	TEST_REAL_SIMILAR(tm.getParameters().getValue("interpolation_step"), 10.5 )
END_SECTION

START_SECTION( virtual void setSamples() )
	// not much to be tested here
END_SECTION

START_SECTION( void getSamples(SamplesType &cont) const )
	// not much to be tested here
END_SECTION

START_SECTION( virtual void setOffset(CoordinateType offset) )

END_SECTION

START_SECTION( CoordinateType getScalingFactor() const )
	TestModel tm;
	tm.setScalingFactor(666.0);

	TEST_REAL_SIMILAR(tm.getParameters().getValue("intensity_scaling"),666.0)
	TEST_REAL_SIMILAR(tm.getScalingFactor(),666.0)
END_SECTION

START_SECTION( const LinearInterpolation& getInterpolation() const )
	TestModel tm;
	InterpolationModel::LinearInterpolation interpol1;
	InterpolationModel::LinearInterpolation interpol2 = tm.getInterpolation();

	// compare models
	TEST_REAL_SIMILAR(interpol1.getScale(), interpol2.getScale());
	TEST_REAL_SIMILAR(interpol1.getInsideReferencePoint(), interpol2.getInsideReferencePoint());
	TEST_REAL_SIMILAR(interpol1.getOutsideReferencePoint(), interpol2.getOutsideReferencePoint() );

END_SECTION

START_SECTION( IntensityType getIntensity(CoordinateType coord) const )
	const TestModel s;
  TestModel::PositionType pos;
  pos[0]=0.1;
  TEST_REAL_SIMILAR(s.getIntensity(pos), 0.3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
