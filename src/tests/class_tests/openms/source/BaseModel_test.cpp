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

#include <OpenMS/FEATUREFINDER/BaseModel.h>
#include <OpenMS/FEATUREFINDER/BaseModel_impl.h>
#include <OpenMS/CONCEPT/Exception.h>

///////////////////////////

using namespace OpenMS;
using std::stringstream;

START_TEST(BaseModel, "$Id$")

class TestModel : public BaseModel
{
  public:
	TestModel()
		: BaseModel()
	{
		setName("TestModel");
		check_defaults_ = false;
		defaultsToParam_();
	}

	TestModel(const TestModel& source)
		: BaseModel(source)
	{
		updateMembers_();
	}

	~TestModel() override	{}

	virtual TestModel& operator = (const TestModel& source)
	{
		if (&source == this) return *this;

		BaseModel::operator = (source);
		updateMembers_();

		return *this;
	}

	void updateMembers_() override
	{
		BaseModel::updateMembers_();
	}

	IntensityType getIntensity(const PositionType& pos) const override
	{
		return pos[0];
	}

	bool isContained(const PositionType& pos) const override
	{
		return getIntensity(pos)>cut_off_;
	}

	void getSamples(SamplesType& /*cont*/) const override
	{
	}

};


// default ctor
TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
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
  TEST_REAL_SIMILAR(s.getIntensity(pos), 0.1);
}
END_SECTION

START_SECTION((virtual bool isContained(const PositionType &pos) const))
  TestModel s;
  s.setCutOff(0.9);
  TestModel::PositionType pos;
  pos[0]=0.1;
  const TestModel& t = s;
  TEST_EQUAL(t.isContained(pos), false)
END_SECTION

START_SECTION((template <typename PeakType> void fillIntensity(PeakType &peak) const))
  const TestModel t;
  TestModel::PeakType p;
  p.getPosition()[0]=0.1;
  p.setIntensity(0.1f);
  t.fillIntensity(p);
  TEST_REAL_SIMILAR(p.getIntensity(), 0.1)
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
