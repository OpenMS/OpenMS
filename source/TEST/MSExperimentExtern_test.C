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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////

START_TEST(MSExperimentExtern, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MSExperimentExtern<>* ptr = 0;
CHECK((MSExperimentExtern()))
	ptr = new MSExperimentExtern<>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~MSExperimentExtern()))
	delete ptr;
RESULT

CHECK((MSExperimentExtern(const MSExperimentExtern& source)))
  MSExperimentExtern<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts().at(0).setFirstName("Name");
  tmp.setBufferSize(1);
  tmp.updateBuffer();
  
  MSExperimentExtern<> tmp2(tmp);
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.getBufferSize(),1);
RESULT

CHECK( MSExperimentExtern& operator= (const MSExperimentExtern& source) )
  MSExperimentExtern<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.setBufferSize(1);
  tmp.updateBuffer();
  MSSpectrum<> spec;
  Peak1D p;
  p.setMZ(5.0);
  spec.push_back(p);
  p.setMZ(10.0);
  spec.push_back(p);
  tmp.push_back(spec);
  tmp.updateRanges();

  MSExperimentExtern<> tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.getBufferSize(),1);
	TEST_REAL_EQUAL(tmp2.getMin()[1],5.0);
	TEST_REAL_EQUAL(tmp2.getMax()[1],10.0);
	
  tmp2 = MSExperimentExtern<>();
  TEST_EQUAL(tmp2.getContacts().size(),0);  
  TEST_EQUAL(tmp2.getBufferSize(),100);
RESULT

CHECK((bool operator== (const MSExperimentExtern& rhs) const))
  	MSExperimentExtern<> edit,empty;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.getContacts().resize(1);
  	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
  	edit.setBufferSize(1);
  	edit.updateBuffer();
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK((bool operator!= (const MSExperimentExtern& rhs) const))
  	MSExperimentExtern<> edit,empty;
	
	TEST_EQUAL(edit!=empty, false);
	
	edit.getContacts().resize(1);
  	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setBufferSize(1);
  	edit.updateBuffer();
	TEST_EQUAL(edit!=empty, true);
RESULT

CHECK((template<class Container> void get2DData(Container& cont) const))
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	MSExperimentExtern<>::PeakType peak;

	// first spectrum (MS)
	spec.setRT(11.1);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 5;
	peak.setIntensity(47.11);
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 10;
	peak.setIntensity(48.11);
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 15;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);

	// second spectrum (MS/MS)
	spec.getContainer().clear();
	spec.setRT(11.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 6;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 11;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// third spectrum (MS)
	spec.getContainer().clear();
	spec.setRT(12.2);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 20;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 25;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// forth spectrum (MS/MS)
	spec.getContainer().clear();
	spec.setRT(12.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 21;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 26;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 31;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	
	
	//Convert
	DPeakArray< Peak2D > a;
	exp.get2DData(a);

	//Tests
	TEST_REAL_EQUAL(a.size(),5);
	TEST_REAL_EQUAL(a[0].getRT(),11.1);
	TEST_REAL_EQUAL(a[0].getMZ(),5);
	TEST_REAL_EQUAL(a[0].getIntensity(),47.11);
	TEST_REAL_EQUAL(a[1].getRT(),11.1);
	TEST_REAL_EQUAL(a[1].getMZ(),10);
	TEST_REAL_EQUAL(a[1].getIntensity(),48.11);
	TEST_REAL_EQUAL(a[2].getRT(),11.1);
	TEST_REAL_EQUAL(a[2].getMZ(),15);
	TEST_REAL_EQUAL(a[3].getRT(),12.2);
	TEST_REAL_EQUAL(a[3].getMZ(),20);
	TEST_REAL_EQUAL(a[4].getRT(),12.2);
	TEST_REAL_EQUAL(a[4].getMZ(),25);

	//Convert
	DPeakArray< Peak2D> list;
	exp.get2DData(list);

	//Tests
	TEST_REAL_EQUAL(list.size(),5);
	DPeakArray< Peak2D>::const_iterator it = list.begin();
	TEST_REAL_EQUAL(it->getRT(),11.1);
	TEST_REAL_EQUAL(it->getMZ(),5);
	TEST_REAL_EQUAL(it->getIntensity(),47.11);
	++it;
	TEST_REAL_EQUAL(it->getRT(),11.1);
	TEST_REAL_EQUAL(it->getMZ(),10);
	TEST_REAL_EQUAL(it->getIntensity(),48.11);
	++it;
	TEST_REAL_EQUAL(it->getRT(),11.1);
	TEST_REAL_EQUAL(it->getMZ(),15);
	++it;
	TEST_REAL_EQUAL(it->getRT(),12.2);
	TEST_REAL_EQUAL(it->getMZ(),20);
	++it;
	TEST_REAL_EQUAL(it->getRT(),12.2);
	TEST_REAL_EQUAL(it->getMZ(),25);
RESULT

CHECK((template<class Container> void set2DData(Container& cont)))
	MSExperimentExtern<> exp;
	
	// create sample data
	DPeakArray< Peak2D> input;
	
	Peak2D p1;
	p1.setIntensity(1.0);
	p1.setRT(2.0);
	p1.setMZ(3.0);
	input.push_back(p1);
	
	Peak2D p2;
	p2.setIntensity(4.0);
	p2.setRT(5.0);
	p2.setMZ(6.0);
	input.push_back(p2);
	
	Peak2D p3;
	p3.setIntensity(7.5);
	p3.setRT(8.5);
	p3.setMZ(9.5);
	input.push_back(p3);
	
	exp.set2DData(input);
	
	// retrieve data again and check for changes
	DPeakArray< Peak2D> output;
	
	exp.get2DData(output);
	TEST_EQUAL(output==input,true);
RESULT


CHECK((UInt getSize() const))
	MSExperimentExtern<Peak1D > tmp;
	TEST_EQUAL(tmp.getSize(),0);
	
	Peak1D p1;
	MSSpectrum<Peak1D > spec;
	spec.push_back(p1);
	spec.push_back(p1);
	spec.push_back(p1);
	
	tmp.push_back(spec);
	tmp.updateRanges();
	TEST_EQUAL(tmp.getSize(),3);
		
RESULT

CHECK((CoordinateType getMinMZ() const))
	MSExperiment<Peak1D > tmp;
	TEST_REAL_EQUAL(tmp.getMinMZ(),numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMaxMZ() const))
	MSExperiment<Peak1D > tmp;
	TEST_REAL_EQUAL(tmp.getMaxMZ(),-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMinRT() const))
	MSExperiment<Peak1D > tmp;
	TEST_REAL_EQUAL(tmp.getMinRT(),numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMaxRT() const))
	MSExperiment<Peak1D > tmp;
	TEST_REAL_EQUAL(tmp.getMaxRT(),-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK(([EXTRA]const std::vector<UInt>& getMSLevels() const))
	MSExperiment<Peak1D > tmp;
	TEST_EQUAL(tmp.getMSLevels().size(),0)
	TEST_REAL_EQUAL(tmp.getDataRange().min()[1],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[1],-numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().min()[0],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[0],-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK(([EXTRA]const AreaType& getDataRange() const))
	MSExperiment<Peak1D > tmp;
	TEST_REAL_EQUAL(tmp.getDataRange().min()[1],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[1],-numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().min()[0],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[0],-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((virtual void updateRanges()))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setMSLevel(1);
	s.setRT(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);
	
	s.getContainer().clear();
	s.setMSLevel(1);
	s.setRT(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRT(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRT(50.0);
	p.getPosition()[0] = 10.0;
	p.setIntensity(-9.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	
	
	tmp.updateRanges();
	tmp.updateRanges(); //second time to check the initialization
	
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),10.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-10.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),50.0)
	TEST_EQUAL(tmp.getMSLevels().size(),2)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getMSLevels()[1],3)
	TEST_EQUAL(tmp.getSize(),4)
	tmp.updateRanges();
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),10.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-10.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),50.0)
	
	TEST_REAL_EQUAL(tmp.getDataRange().min()[1],5.0)
	TEST_REAL_EQUAL(tmp.getDataRange().max()[1],10.0)
	TEST_REAL_EQUAL(tmp.getDataRange().min()[0],30.0)
	TEST_REAL_EQUAL(tmp.getDataRange().max()[0],50.0)
	
	TEST_EQUAL(tmp.getMSLevels().size(),2)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getMSLevels()[1],3)
	
	TEST_EQUAL(tmp.getSize(),4)
	
	//Update for MS level 1
	
	tmp.updateRanges(1);
	tmp.updateRanges(1);
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),7.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-7.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)
	tmp.updateRanges(1);
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),7.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-7.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)

	//test with only one peak

	MSExperiment< Peak1D > tmp2;
	MSSpectrum< Peak1D > s2;
	Peak1D p2;
	
	s2.setRT(30.0);
	p2.getPosition()[0] = 5.0;
	p2.setIntensity(-5.0);
	s2.getContainer().push_back(p2);
	tmp2.push_back(s2);

	tmp2.updateRanges();
	TEST_REAL_EQUAL(tmp2.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMaxMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMinInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp2.getMaxRT(),30.0)

	tmp2.updateRanges(1);
	TEST_REAL_EQUAL(tmp2.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMaxMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMinInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp2.getMaxRT(),30.0)
	
RESULT

CHECK((void updateRanges(Int ms_level)))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setMSLevel(1);
	s.setRT(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);
	
	s.getContainer().clear();
	s.setMSLevel(1);
	s.setRT(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRT(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRT(50.0);
	p.getPosition()[0] = 10.0;
	p.setIntensity(-9.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	
	
	//Update for MS level 1
	
	tmp.updateRanges(1);
	tmp.updateRanges(1);
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),7.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-7.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)
	tmp.updateRanges(1);
	TEST_REAL_EQUAL(tmp.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp.getMaxMZ(),7.0)
	TEST_REAL_EQUAL(tmp.getMinInt(),-7.0)
	TEST_REAL_EQUAL(tmp.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp.getMaxRT(),40.0)
	TEST_EQUAL(tmp.getMSLevels().size(),1)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getSize(),2)

	//test with only one peak

	MSExperiment< Peak1D > tmp2;
	MSSpectrum< Peak1D > s2;
	Peak1D p2;
	
	s2.setRT(30.0);
	p2.getPosition()[0] = 5.0;
	p2.setIntensity(-5.0);
	s2.getContainer().push_back(p2);
	tmp2.push_back(s2);

	tmp2.updateRanges(1);
	TEST_REAL_EQUAL(tmp2.getMinMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMaxMZ(),5.0)
	TEST_REAL_EQUAL(tmp2.getMinInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMaxInt(),-5.0)
	TEST_REAL_EQUAL(tmp2.getMinRT(),30.0)
	TEST_REAL_EQUAL(tmp2.getMaxRT(),30.0)
	
RESULT

CHECK((Iterator RTBegin(double rt)))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setRT(30.0);
	tmp.push_back(s);	
	s.setRT(40.0);
	tmp.push_back(s);	
	s.setRT(45.0);
	tmp.push_back(s);	
	s.setRT(50.0);
	tmp.push_back(s);
	
	MSExperiment< Peak1D >::Iterator it;
	
	it = tmp.RTBegin(20.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((Iterator RTEnd(double rt)))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setRT(30.0);
	tmp.push_back(s);	
	s.setRT(40.0);
	tmp.push_back(s);	
	s.setRT(45.0);
	tmp.push_back(s);	
	s.setRT(50.0);
	tmp.push_back(s);
	
	MSExperiment< Peak1D >::Iterator it;
	
	it = tmp.RTEnd(20.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((ConstIterator RTBegin(double rt) const))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setRT(30.0);
	tmp.push_back(s);	
	s.setRT(40.0);
	tmp.push_back(s);	
	s.setRT(45.0);
	tmp.push_back(s);	
	s.setRT(50.0);
	tmp.push_back(s);
	
	MSExperiment< Peak1D >::Iterator it;
	
	it = tmp.RTBegin(20.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((ConstIterator RTEnd(double rt) const))
	MSExperiment< Peak1D > tmp;
	MSSpectrum< Peak1D > s;
	Peak1D p;
	
	s.setRT(30.0);
	tmp.push_back(s);	
	s.setRT(40.0);
	tmp.push_back(s);	
	s.setRT(45.0);
	tmp.push_back(s);	
	s.setRT(50.0);
	tmp.push_back(s);
	
	MSExperiment< Peak1D >::Iterator it;
	
	it = tmp.RTEnd(20.0);
	TEST_REAL_EQUAL(it->getRT(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_EQUAL(it->getRT(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((void sortSpectra(bool sort_mz = true)))
	DPeakArray< Peak2D> plist;
	
	Peak2D p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 5.0;
	plist.push_back(p1);
		
	Peak2D p2;
	p2.getPosition()[0] = 1.0;
	p2.getPosition()[1] = 3.0;
	plist.push_back(p2);
		
	Peak2D p3;
	p3.getPosition()[0] = 2.0;
	p3.getPosition()[1] = 14.0;
	plist.push_back(p3);
	
	Peak2D p4;
	p4.getPosition()[0] = 2.0;
	p4.getPosition()[1] = 11.0;
	plist.push_back(p4);
	
	MSExperimentExtern<> exp;
	exp.set2DData(plist);
	
	exp.sortSpectra(true);
 	
	TEST_REAL_EQUAL(exp[0][0].getMZ(),3.0);
	TEST_REAL_EQUAL(exp[0][1].getMZ(),5.0);
	TEST_REAL_EQUAL(exp[1][0].getMZ(),11.0);
	TEST_REAL_EQUAL(exp[1][1].getMZ(),14.0);
RESULT

CHECK((void reset()))
	DPeakArray< Peak2D> plist;
	
	Peak2D p;
	p.getPosition()[0] = 1.0;
	p.getPosition()[1] = 5.0;
	plist.push_back(p);
	p.getPosition()[0] = 2.0;
	p.getPosition()[1] = 3.0;
	plist.push_back(p);
		
	MSExperimentExtern<> exp;
	exp.set2DData(plist);
	exp.updateRanges();
	
	exp.reset();
	
	TEST_EQUAL( (exp.size()==0),true);
RESULT

CHECK((const ExperimentalSettings& getExperimentalSettings() const))
	MSExperimentExtern<> exp;
	exp.setComment("test");
	TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
RESULT

CHECK((ExperimentalSettings& getExperimentalSettings()))
	MSExperimentExtern<> exp;
	exp.getExperimentalSettings().setComment("test");
	TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
RESULT

CHECK((AIterator areaBegin(const AreaType &area)))
  // write tests for areaBegin(), areaEnd() etc (whoever wrote these methods)
RESULT

CHECK((AIterator areaEnd()))
  // write tests for areaBegin(), areaEnd() etc (whoever wrote these methods)
RESULT

CHECK((AConstIterator areaBegin(const AreaType &area) const))
  // write tests for areaBegin(), areaEnd() etc (whoever wrote these methods)
RESULT

CHECK((AConstIterator areaEnd() const))
  // write tests for areaBegin(), areaEnd() etc (whoever wrote these methods)
RESULT

CHECK((ExperimentType& getType() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getType(),ExperimentalSettings::UNKNOWN);
RESULT

CHECK((const Date& getDate() const))
  ExperimentalSettings tmp;
  String s;
  tmp.getDate().get(s);
  TEST_EQUAL(s,"0000-00-00");
RESULT

CHECK((const HPLC& getHPLC() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getHPLC()==HPLC(),true);
RESULT

CHECK((const Instrument& getInstrument() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getInstrument()==Instrument(),true);
RESULT

CHECK((const ProcessingMethod& getProcessingMethod() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getProcessingMethod()==ProcessingMethod(),true);
RESULT

CHECK((const Sample& getSample() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSample()==Sample(),true);
RESULT

CHECK((const Software& getSoftware() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSoftware()==Software(),true);
RESULT

CHECK((const SourceFile& getSourceFile() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(),true);
RESULT

CHECK((const String& getComment() const))
	ExperimentalSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
RESULT

CHECK((void setComment(const String& comment)))
	ExperimentalSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
RESULT

CHECK((const std::vector<ContactPerson>& getContacts() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getContacts().size(),0);
RESULT

CHECK((void setContacts(const std::vector<ContactPerson>& contacts)))
  ExperimentalSettings tmp;
  std::vector<ContactPerson> dummy;
  ContactPerson c;
  c.setFirstName("bla17");
  c.setLastName("blubb17");
  dummy.push_back(c);
  c.setFirstName("bla18");
  c.setLastName("blubb18");
  dummy.push_back(c);
  tmp.setContacts(dummy);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getFirstName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getFirstName(),"bla18");
  TEST_EQUAL(tmp.getContacts()[0].getLastName(),"blubb17");
  TEST_EQUAL(tmp.getContacts()[1].getLastName(),"blubb18");
RESULT

CHECK((void setDate(const Date& date)))
  ExperimentalSettings tmp;
  Date dummy;
  String s;
  dummy.set("02/07/2006");
  tmp.setDate(dummy);
  tmp.getDate().get(s);
  TEST_EQUAL(s,"2006-02-07");
RESULT

CHECK((void setHPLC(const HPLC& hplc)))
  ExperimentalSettings tmp;
  HPLC dummy;
  dummy.setFlux(5);
  tmp.setHPLC(dummy);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
RESULT

CHECK((void setInstrument(const Instrument& instrument)))
  ExperimentalSettings tmp;
  Instrument dummy;
  dummy.setName("bla");
  tmp.setInstrument(dummy);
  TEST_EQUAL(tmp.getInstrument().getName(),"bla");
RESULT

CHECK((void setProcessingMethod(const ProcessingMethod& processing_method)))
  ExperimentalSettings tmp;
  ProcessingMethod dummy;
  dummy.setDeisotoping(true);
  tmp.setProcessingMethod(dummy);
  TEST_EQUAL(tmp.getProcessingMethod().getDeisotoping(),true);
RESULT

CHECK((void setSample(const Sample& sample)))
  ExperimentalSettings tmp;
  Sample dummy;
  dummy.setName("bla3");
  tmp.setSample(dummy);
  TEST_EQUAL(tmp.getSample().getName(),"bla3");
RESULT

CHECK((void setSoftware(const Software& software)))
  ExperimentalSettings tmp;
  Software dummy;
  dummy.setName("bla33");
  tmp.setSoftware(dummy);
  TEST_EQUAL(tmp.getSoftware().getName(),"bla33");
RESULT

CHECK((void setSourceFile(const SourceFile& source_file)))
  ExperimentalSettings tmp;
  SourceFile dummy;
  dummy.setNameOfFile("bla4");
  tmp.setSourceFile(dummy);
  TEST_EQUAL(tmp.getSourceFile().getNameOfFile(),"bla4");
RESULT

CHECK((void setType(ExperimentType type)))
  ExperimentalSettings tmp;
  tmp.setType(ExperimentalSettings::HPLC_MS);
  TEST_EQUAL(tmp.getType(),ExperimentalSettings::HPLC_MS);
RESULT

CHECK((HPLC& getHPLC()))
  ExperimentalSettings tmp;
  tmp.getHPLC().setFlux(5);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
RESULT

CHECK((Instrument& getInstrument()))
  ExperimentalSettings tmp;
  tmp.getInstrument().setName("bla55");
  TEST_EQUAL(tmp.getInstrument().getName(),"bla55");
RESULT

CHECK((ProcessingMethod& getProcessingMethod()))
  ExperimentalSettings tmp;
  tmp.getProcessingMethod().setDeisotoping(true);
  TEST_EQUAL(tmp.getProcessingMethod().getDeisotoping(),true);
RESULT

CHECK((Sample& getSample()))
  ExperimentalSettings tmp;
  tmp.getSample().setName("bla2");
  TEST_EQUAL(tmp.getSample().getName(),"bla2");
RESULT

CHECK((Software& getSoftware()))
  ExperimentalSettings tmp;
  tmp.getSoftware().setName("bla3");
  TEST_EQUAL(tmp.getSoftware().getName(),"bla3");
RESULT

CHECK((SourceFile& getSourceFile()))
  ExperimentalSettings tmp;
  tmp.getSourceFile().setNameOfFile("bla4");
  TEST_EQUAL(tmp.getSourceFile().getNameOfFile(),"bla4");
RESULT

CHECK((std::vector<ContactPerson>& getContacts()))
  ExperimentalSettings tmp;
  ContactPerson c;
  c.setFirstName("bla17");
  tmp.getContacts().push_back(c);
  c.setFirstName("bla18");
  tmp.getContacts().push_back(c);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getFirstName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getFirstName(),"bla18");
RESULT


CHECK((void reserve(size_type n)))
	MSExperimentExtern<> exp;
	exp.reserve(10000);
RESULT

CHECK( void clear() )
	MSExperimentExtern<> exp;
	
	// create sample data
	DPeakArray< Peak2D> input;
	
	Peak2D p1;
	p1.setIntensity(1.0);
	p1.setRT(2.0);
	p1.setMZ(3.0);
	input.push_back(p1);
	
	Peak2D p2;
	p2.setIntensity(4.0);
	p2.setRT(5.0);
	p2.setMZ(6.0);
	input.push_back(p2);
	
	Peak2D p3;
	p3.setIntensity(7.5);
	p3.setRT(8.5);
	p3.setMZ(9.5);
	input.push_back(p3);
	
	exp.set2DData(input);

	// empty data structure
	exp.clear();
	TEST_EQUAL(exp.size(),0)
RESULT	

CHECK( size_type size() const )
	
	MSExperimentExtern<> exp;
	DPeakArray< Peak2D> input;
	
	Peak2D p1;
	p1.setIntensity(1.0);
	p1.setRT(2.0);
	p1.setMZ(3.0);
	input.push_back(p1);
	
	Peak2D p2;
	p2.setIntensity(4.0);
	p2.setRT(5.0);
	p2.setMZ(6.0);
	input.push_back(p2);
	
	Peak2D p3;
	p3.setIntensity(7.5);
	p3.setRT(5.0);
	p3.setMZ(9.5);
	input.push_back(p3);
	
	exp.set2DData(input);
	TEST_EQUAL(exp.size(),2) // two scans
RESULT

CHECK( UInt getBufferSize() const )
	// will be tested below
RESULT

CHECK( UInt& getBufferSize() )
	// will be tested below	
RESULT

CHECK( void setBufferSize(UInt sz) )
	// will be tested below	
RESULT

CHECK( void updateBuffer() )
	MSExperimentExtern<> exp;
	exp.setBufferSize(999);
	exp.updateBuffer();
	
	TEST_EQUAL(exp.getBufferSize(),999)

RESULT

CHECK( void push_back(const SpectrumType &spec) )
	MSExperimentExtern<> exp;
	
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=1;i<=10;++i)
	{
		for (UInt j=0;j<4;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
	
	TEST_EQUAL(exp.size(),10)	
RESULT

CHECK(Iterator end())
	// will be tested below	
RESULT

CHECK(Iterator begin())
	MSExperimentExtern<> exp;
	exp.setBufferSize(2);
	exp.updateBuffer();
	
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=1;i<5;++i)
	{
		for (UInt j=0;j<4;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
	
	UInt rt = 1;
		
	for (MSExperimentExtern<>::Iterator it = exp.begin();
				it != exp.end();
				++it)
	{
		for (UInt mz = 0;mz <4;++mz)
		{
			TEST_REAL_EQUAL(it->getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(it->getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(it->getRT(),rt)
		++rt;
	}		
RESULT

CHECK( ConstIterator end() const )
	// test below
RESULT

CHECK( ConstIterator begin() const )

	MSExperimentExtern<> exp;
	exp.setBufferSize(2);	// check
	exp.updateBuffer();
	
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=1;i<5;++i)
	{
		for (UInt j=0;j<4;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
	
	UInt rt = 1;
	
	const MSExperimentExtern<> cexp(exp);
		
	for (MSExperimentExtern<>::ConstIterator cit = cexp.begin();
				cit != cexp.end();
				++cit)
	{
		for (UInt mz = 0;mz <4;++mz)
		{
			TEST_REAL_EQUAL(cit->getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(cit->getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(cit->getRT(),rt)
		++rt;
	}	
RESULT

CHECK( void resize(UInt new_size) )
	MSExperimentExtern<> exp;
	exp.resize(5000);
	
	TEST_EQUAL(exp.capacity(),5000)	
RESULT

CHECK( size_type capacity() const)
 	// test above
RESULT

CHECK( reference operator[](size_type n) )
	MSExperimentExtern<> exp;

	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=0;i<5;++i)
	{
		for (UInt j=0;j<5;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
		
	for (UInt rt=0; rt <5;++rt)
	{
		for (UInt mz = 0;mz <5;++mz)
		{
			TEST_REAL_EQUAL(exp[rt].getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(exp[rt].getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(exp[rt].getRT(),rt)
	}	
RESULT

CHECK( const_reference operator[](size_type n) const )
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=0;i<5;++i)
	{
		for (UInt j=0;j<5;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
	
	const MSExperimentExtern<> cexp(exp);
		
	for (UInt rt=0; rt <5;++rt)
	{
		for (UInt mz = 0;mz <5;++mz)
		{
			TEST_REAL_EQUAL(cexp[rt].getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(cexp[rt].getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(cexp[rt].getRT(),rt)
	}	
RESULT

CHECK(reference at(size_type n))
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=0;i<5;++i)
	{
		for (UInt j=0;j<5;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
		
	for (UInt rt=0; rt <5;++rt)
	{
		for (UInt mz = 0;mz <5;++mz)
		{
			TEST_REAL_EQUAL(exp.at(rt).getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(exp.at(rt).getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(exp.at(rt).getRT(),rt)
	}	
RESULT

CHECK( const_reference at(size_type n) const )
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=0;i<5;++i)
	{
		for (UInt j=0;j<5;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		spec.clear();
	}
	
	const MSExperimentExtern<> cexp(exp);
		
	for (UInt rt=0; rt <5;++rt)
	{
		for (UInt mz = 0;mz <5;++mz)
		{
			TEST_REAL_EQUAL(cexp.at(rt).getContainer()[mz].getMZ(),mz)
			TEST_REAL_EQUAL(cexp.at(rt).getContainer()[mz].getIntensity(), (mz*rt) )
		}
		TEST_REAL_EQUAL(cexp.at(rt).getRT(),rt)
	}	
RESULT

CHECK( reference back() )
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	Peak1D p;	
	for (UInt i=0;i<5;++i)
	{
		for (UInt j=0;j<5;++j)
		{			
			p.setMZ(j);
			p.setIntensity(i*j);
			spec.push_back(p);
		}
		spec.setRT(i);
		exp.push_back(spec);
		TEST_REAL_EQUAL(exp.back().getRT(),spec.getRT() )
		spec.clear();
	}
	
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
