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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <sstream>

///////////////////////////

START_TEST(Sample, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PRECISION(0.001)

// default ctor
Sample* dv_ptr = 0;
CHECK((Sample()))
	dv_ptr = new Sample;
	TEST_NOT_EQUAL(dv_ptr, 0)
RESULT

// destructor
CHECK((~Sample()))
	delete dv_ptr;
RESULT

CHECK((const String& getName() const))
	Sample s;
	TEST_EQUAL(s.getName(),"")
RESULT

CHECK((const String& getOrganism() const))
	Sample s;
	TEST_EQUAL(s.getOrganism(),"")
RESULT

CHECK((const String& getNumber() const))
	Sample s;
	TEST_EQUAL(s.getNumber(),"")
RESULT

CHECK((const String& getComment() const))
	Sample s;
	TEST_EQUAL(s.getComment(),"")
RESULT

CHECK((SampleState getState() const))
	Sample s;
	TEST_EQUAL(s.getState(),Sample::SAMPLENULL)
RESULT

CHECK((float getMass() const))
	Sample s;
	TEST_REAL_EQUAL(s.getMass(),0.0)
RESULT

CHECK((float getVolume() const))
	Sample s;
	TEST_REAL_EQUAL(s.getVolume(),0.0)
RESULT

CHECK((float getConcentration() const))
	Sample s;
	TEST_REAL_EQUAL(s.getConcentration(),0.0)
RESULT

CHECK((void setName(const String& name)))
	Sample s;
	s.setName("TTEST");
	TEST_EQUAL(s.getName(),"TTEST")
RESULT

CHECK((void setOrganism(const String& organism)))
	Sample s;
	s.setOrganism("TTEST");
	TEST_EQUAL(s.getOrganism(),"TTEST")
RESULT

CHECK((void setNumber(const String& number)))
	Sample s;
	s.setNumber("Sample4711");
	TEST_EQUAL(s.getNumber(),"Sample4711")
RESULT

CHECK((void setComment(const String& comment)))
	Sample s;
	s.setComment("Sample Description");
	TEST_EQUAL(s.getComment(),"Sample Description")
RESULT

CHECK((void setState(SampleState state)))
	Sample s;
	s.setState(Sample::LIQUID);
	TEST_EQUAL(s.getState(),Sample::LIQUID)
RESULT

CHECK((void setMass(float mass)))
	Sample s;
	s.setMass(4711.2);
	TEST_REAL_EQUAL(s.getMass(),4711.2)
RESULT

CHECK((void setVolume(float volume)))
	Sample s;
	s.setVolume(4711.3);
	TEST_REAL_EQUAL(s.getVolume(),4711.3)
RESULT

CHECK((void setConcentration(float concentration)))
	Sample s;
	s.setConcentration(4711.4);
	TEST_REAL_EQUAL(s.getConcentration(),4711.4)
RESULT

CHECK((const std::vector<Sample>& getSubsamples() const))
	Sample s;
	TEST_EQUAL(s.getSubsamples().size(),0)
RESULT

CHECK((std::vector<Sample>& getSubsamples()))
	Sample s,s2;
	s.getSubsamples().push_back(s2);
	TEST_EQUAL(s.getSubsamples().size(),1)
RESULT

CHECK((void setSubsamples(const std::vector<Sample>& subsamples)))
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
RESULT

//treatments

CHECK((Int countTreatments() const))
	Sample s;
	TEST_EQUAL(s.countTreatments(),0)
	Digestion d;
	s.addTreatment(d);
	TEST_EQUAL(s.countTreatments(),1)
RESULT

CHECK((const SampleTreatment& getTreatment(UInt position) const ))
	Sample s;
	TEST_EXCEPTION(Exception::IndexOverflow, s.getTreatment(0))
RESULT

CHECK((void addTreatment(const SampleTreatment& treatment, Int before_position=-1) ))
	Sample s;
	Digestion d;
	Modification m,m2,m3;
	Tagging t;
	
	//different treatments
	d.setEnzyme("D");
	m.setReagentName("m");
	t.setMassShift(5.0);
	s.addTreatment(d);
	s.addTreatment(m);
	s.addTreatment(t);
	TEST_EQUAL(s.countTreatments(),3)
	TEST_EQUAL(s.getTreatment(0).getType(),"Digestion")
	TEST_EQUAL(s.getTreatment(1).getType(),"Modification")
	TEST_EQUAL(s.getTreatment(2).getType(),"Tagging")
	
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(0))).getEnzyme(),"D")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(1))).getReagentName(),"m")
	TEST_REAL_EQUAL((dynamic_cast<const Tagging&>(s.getTreatment(2))).getMassShift(),5.0)
RESULT

CHECK((SampleTreatment& getTreatment(UInt position) ))
	Sample s;
	Digestion d;
	s.addTreatment(d);
	(dynamic_cast<Digestion&>(s.getTreatment(0))).setEnzyme("bluff");
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(0))).getEnzyme(),"bluff")
RESULT

CHECK((void removeTreatment(UInt position) ))
	Sample s;
	Digestion d;
	Modification m,m2,m3;
	Tagging t;
	
	//different treatments
	d.setEnzyme("D");
	m.setReagentName("m");
	t.setMassShift(5.0);
	s.addTreatment(d);
	s.addTreatment(m);
	s.addTreatment(t);
	
	//removeTreatment
	m2.setReagentName("m2");
	m3.setReagentName("m3");
	s.addTreatment(m2,0);
	s.addTreatment(m3,3);
	TEST_EQUAL(s.countTreatments(),5)
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(0))).getReagentName(),"m2")
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(1))).getEnzyme(),"D")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(2))).getReagentName(),"m")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(3))).getReagentName(),"m3")
	TEST_REAL_EQUAL((dynamic_cast<const Tagging&>(s.getTreatment(4))).getMassShift(),5.0)	
	
	s.removeTreatment(4);
	TEST_EQUAL(s.countTreatments(),4)
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(0))).getReagentName(),"m2")
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(1))).getEnzyme(),"D")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(2))).getReagentName(),"m")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(3))).getReagentName(),"m3")

	s.removeTreatment(2);
	TEST_EQUAL(s.countTreatments(),3)
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(0))).getReagentName(),"m2")
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(1))).getEnzyme(),"D")
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(2))).getReagentName(),"m3")

	s.removeTreatment(2);
	TEST_EQUAL(s.countTreatments(),2)
	TEST_EQUAL((dynamic_cast<const Modification&>(s.getTreatment(0))).getReagentName(),"m2")
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(1))).getEnzyme(),"D")

	//exceptions: index overflow
	TEST_EXCEPTION(Exception::IndexOverflow,s.removeTreatment(2))
	TEST_EXCEPTION(Exception::IndexOverflow,s.getTreatment(2))
	TEST_EXCEPTION(Exception::IndexOverflow,s.addTreatment(m,3))
RESULT

//copy ctr
CHECK((Sample(const Sample& source)))
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
	
	//treatments
	Digestion d;
	d.setEnzyme("D");
	s.addTreatment(d);
	
	//-----------------
	//Copy construction	
	//-----------------
	Sample s2(s);
	
	//basic stuff
	TEST_EQUAL(s2.getName(),"TTEST")
	TEST_EQUAL(s2.getNumber(),"Sample4711")
	TEST_EQUAL(s2.getComment(),"Sample Description")
	TEST_EQUAL(s2.getState(),Sample::LIQUID)
	TEST_REAL_EQUAL(s2.getMass(),4711.2)
	TEST_REAL_EQUAL(s2.getVolume(),4711.3)
	TEST_REAL_EQUAL(s2.getConcentration(),4711.4)	
	TEST_EQUAL(s2.getOrganism(),"TTEST2")
	
	//meta
	TEST_EQUAL("horse",s.getMetaValue("label"))
	
	//subsamples
	TEST_EQUAL("2",s.getSubsamples()[0].getName())
	
	//treatments
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(0))).getEnzyme(),"D")
RESULT

//assignment operator
CHECK((Sample& operator= (const Sample& source)))
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
	
	//treatments
	Digestion d;
	d.setEnzyme("D");
	s.addTreatment(d);
	
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
	TEST_REAL_EQUAL(s2.getMass(),4711.2)
	TEST_REAL_EQUAL(s2.getVolume(),4711.3)
	TEST_REAL_EQUAL(s2.getConcentration(),4711.4)	

	//meta
	TEST_EQUAL("horse",s.getMetaValue("label"))
	
	//subsamples
	TEST_EQUAL("2",s.getSubsamples()[0].getName())
	
	//treatments
	TEST_EQUAL((dynamic_cast<const Digestion&>(s.getTreatment(0))).getEnzyme(),"D")
RESULT

CHECK((bool operator== (const Sample& rhs) const))
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

	edit.addTreatment(Modification());
	TEST_EQUAL(edit==empty,false)
	edit = empty;
	TEST_EQUAL(edit==empty,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
