// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/Sample.h>
#include <sstream>

///////////////////////////

START_TEST(Sample, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TOLERANCE_ABSOLUTE(0.001)

// default ctor
Sample* dv_ptr = nullptr;
Sample* dv_nullPointer = nullptr;
START_SECTION((Sample()))
	dv_ptr = new Sample;
  TEST_NOT_EQUAL(dv_ptr, dv_nullPointer)
END_SECTION

// destructor
START_SECTION((~Sample()))
	delete dv_ptr;
END_SECTION

START_SECTION((const String& getName() const))
	Sample s;
	TEST_EQUAL(s.getName(),"")
END_SECTION

START_SECTION((const String& getOrganism() const))
	Sample s;
	TEST_EQUAL(s.getOrganism(),"")
END_SECTION

START_SECTION((const String& getNumber() const))
	Sample s;
	TEST_EQUAL(s.getNumber(),"")
END_SECTION

START_SECTION((const String& getComment() const))
	Sample s;
	TEST_EQUAL(s.getComment(),"")
END_SECTION

START_SECTION((SampleState getState() const))
	Sample s;
	TEST_EQUAL(s.getState(),Sample::SAMPLENULL)
END_SECTION

START_SECTION((double getMass() const ))
	Sample s;
	TEST_REAL_SIMILAR(s.getMass(),0.0)
END_SECTION

START_SECTION((double getVolume() const ))
	Sample s;
	TEST_REAL_SIMILAR(s.getVolume(),0.0)
END_SECTION

START_SECTION((double getConcentration() const ))
	Sample s;
	TEST_REAL_SIMILAR(s.getConcentration(),0.0)
END_SECTION

START_SECTION((void setName(const String& name)))
	Sample s;
	s.setName("TTEST");
	TEST_EQUAL(s.getName(),"TTEST")
END_SECTION

START_SECTION((void setOrganism(const String& organism)))
	Sample s;
	s.setOrganism("TTEST");
	TEST_EQUAL(s.getOrganism(),"TTEST")
END_SECTION

START_SECTION((void setNumber(const String& number)))
	Sample s;
	s.setNumber("Sample4711");
	TEST_EQUAL(s.getNumber(),"Sample4711")
END_SECTION

START_SECTION((void setComment(const String& comment)))
	Sample s;
	s.setComment("Sample Description");
	TEST_EQUAL(s.getComment(),"Sample Description")
END_SECTION

START_SECTION((void setState(SampleState state)))
	Sample s;
	s.setState(Sample::LIQUID);
	TEST_EQUAL(s.getState(),Sample::LIQUID)
END_SECTION

START_SECTION((void setMass(double mass)))
	Sample s;
	s.setMass(4711.2);
	TEST_REAL_SIMILAR(s.getMass(),4711.2)
END_SECTION

START_SECTION((void setVolume(double volume)))
	Sample s;
	s.setVolume(4711.3);
	TEST_REAL_SIMILAR(s.getVolume(),4711.3)
END_SECTION

START_SECTION((void setConcentration(double concentration)))
	Sample s;
	s.setConcentration(4711.4);
	TEST_REAL_SIMILAR(s.getConcentration(),4711.4)
END_SECTION

START_SECTION((const std::vector<Sample>& getSubsamples() const))
	Sample s;
	TEST_EQUAL(s.getSubsamples().size(),0)
END_SECTION

START_SECTION((std::vector<Sample>& getSubsamples()))
	Sample s,s2;
	s.getSubsamples().push_back(s2);
	TEST_EQUAL(s.getSubsamples().size(),1)
END_SECTION

START_SECTION((void setSubsamples(const std::vector<Sample>& subsamples)))
	Sample s,s2,s3;
	vector<Sample> v;

	//size=2
	s2.setName("2");
	s3.setName("3");
	v.push_back(s2);
	v.push_back(s3);
	s.setSubsamples(v);
	TEST_EQUAL(s.getSubsamples().size(),2)
	TEST_EQUAL(s.getSubsamples()[0].getName(),"2")
	TEST_EQUAL(s.getSubsamples()[1].getName(),"3")
END_SECTION

//copy ctr
START_SECTION((Sample(const Sample& source)))
	Sample s;

	//basic stuff
	s.setOrganism("TTEST2");
	s.setName("TTEST");
	s.setNumber("Sample4711");
	s.setComment("Sample Description");
	s.setState(Sample::LIQUID);
	s.setMass(4711.2);
	s.setVolume(4711.3);
	s.setConcentration(4711.4);

	//meta info
	s.setMetaValue("label",String("horse"));

	//subsamples
	Sample ss;
	ss.setName("2");
	s.getSubsamples().push_back(ss);

	//-----------------
	//Copy construction
	//-----------------
	Sample s2(s);

	//basic stuff
	TEST_EQUAL(s2.getName(),"TTEST")
	TEST_EQUAL(s2.getNumber(),"Sample4711")
	TEST_EQUAL(s2.getComment(),"Sample Description")
	TEST_EQUAL(s2.getState(),Sample::LIQUID)
	TEST_REAL_SIMILAR(s2.getMass(),4711.2)
	TEST_REAL_SIMILAR(s2.getVolume(),4711.3)
	TEST_REAL_SIMILAR(s2.getConcentration(),4711.4)
	TEST_EQUAL(s2.getOrganism(),"TTEST2")

	//meta
	TEST_EQUAL(s.getMetaValue("label"),"horse")

	//subsamples
	TEST_EQUAL(s.getSubsamples()[0].getName(),"2")
END_SECTION

//assignment operator
START_SECTION((Sample& operator= (const Sample& source)))
	Sample s;

	//basic stuff
	s.setName("TTEST");
	s.setOrganism("TTEST2");
	s.setNumber("Sample4711");
	s.setComment("Sample Description");
	s.setState(Sample::LIQUID);
	s.setMass(4711.2);
	s.setVolume(4711.3);
	s.setConcentration(4711.4);

	//meta
	s.setMetaValue("label",String("horse"));

	//subsamples
	Sample ss;
	ss.setName("2");
	s.getSubsamples().push_back(ss);

	//-----------------
	//Copy construction
	//-----------------
	Sample s2;
	s2=s;

	//basic stuff
	TEST_EQUAL(s2.getName(),"TTEST")
	TEST_EQUAL(s2.getNumber(),"Sample4711")
	TEST_EQUAL(s2.getComment(),"Sample Description")
	TEST_EQUAL(s2.getOrganism(),"TTEST2")
	TEST_EQUAL(s2.getState(),Sample::LIQUID)
	TEST_REAL_SIMILAR(s2.getMass(),4711.2)
	TEST_REAL_SIMILAR(s2.getVolume(),4711.3)
	TEST_REAL_SIMILAR(s2.getConcentration(),4711.4)

	//meta
	TEST_EQUAL(s.getMetaValue("label"),"horse")

	//subsamples
	TEST_EQUAL(s.getSubsamples()[0].getName(),"2")
END_SECTION

START_SECTION((bool operator== (const Sample& rhs) const))
	const Sample empty;
	Sample edit;

	TEST_EQUAL(edit==empty,true)

	edit.setName("TTEST");
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setOrganism("TTEST2");
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setNumber("Sample4711");
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setComment("Sample Description");
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setState(Sample::LIQUID);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setMass(4711.2);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setVolume(4711.3);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setConcentration(4711.4);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.getSubsamples().push_back(empty);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)

	edit.setMetaValue("color",45);
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
