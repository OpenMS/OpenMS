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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(MSExperiment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

const int MZ = DimensionDescription < LCMS_Tag >::MZ;
const int RT = DimensionDescription < LCMS_Tag >::RT;

MSExperiment<>* ptr = 0;
CHECK((MSExperiment()))
	ptr = new MSExperiment<>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(([EXTRA]~MSExperiment()))
	delete ptr;
RESULT

CHECK((MSExperiment(const MSExperiment& source)))
  MSExperiment<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  
  MSExperiment<> tmp2(tmp);
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
RESULT

CHECK((MSExperiment& operator= (const MSExperiment& source)))
  MSExperiment<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  DPeak<1> p;
  p.setPos(5.0);
  tmp[0].push_back(p);
  p.setPos(10.0);
  tmp[0].push_back(p);
  tmp.updateRanges();
  
  MSExperiment<> tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_EQUAL(tmp2.getMinMZ(),5.0);
	TEST_REAL_EQUAL(tmp2.getMaxMZ(),10.0);
	
  tmp2 = MSExperiment<>();
  TEST_EQUAL(tmp2.getContacts().size(),0);  
  TEST_EQUAL(tmp2.size(),0);
RESULT

CHECK((bool operator== (const MSExperiment& rhs) const))
  MSExperiment<> edit,empty;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty, false);
	
	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK((bool operator!= (const MSExperiment& rhs) const))
  MSExperiment<> edit,empty;
	
	TEST_EQUAL(edit!=empty, false);
	
	edit.getContacts().resize(1);
  TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit!=empty, true);
RESULT

CHECK((template<class Container> void get2DData(Container& cont) const))
	MSExperiment<> exp;
	MSExperiment<>::SpectrumType spec;
	MSExperiment<>::PeakType peak;

	// first spectrum (MS)
	spec.setRetentionTime(11.1);
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
	spec.setRetentionTime(11.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 6;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 11;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// third spectrum (MS)
	spec.getContainer().clear();
	spec.setRetentionTime(12.2);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 20;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 25;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// forth spectrum (MS/MS)
	spec.getContainer().clear();
	spec.setRetentionTime(12.5);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 21;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 26;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 31;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	
	
	//Convert
	DPeakArray<2, DRawDataPoint<2> > a;
	exp.get2DData(a);

	//Tests
	TEST_REAL_EQUAL(a.size(),5);
	TEST_REAL_EQUAL(a[0].getPosition()[RT],11.1);
	TEST_REAL_EQUAL(a[0].getPosition()[MZ],5);
	TEST_REAL_EQUAL(a[0].getIntensity(),47.11);
	TEST_REAL_EQUAL(a[1].getPosition()[RT],11.1);
	TEST_REAL_EQUAL(a[1].getPosition()[MZ],10);
	TEST_REAL_EQUAL(a[1].getIntensity(),48.11);
	TEST_REAL_EQUAL(a[2].getPosition()[RT],11.1);
	TEST_REAL_EQUAL(a[2].getPosition()[MZ],15);
	TEST_REAL_EQUAL(a[3].getPosition()[RT],12.2);
	TEST_REAL_EQUAL(a[3].getPosition()[MZ],20);
	TEST_REAL_EQUAL(a[4].getPosition()[RT],12.2);
	TEST_REAL_EQUAL(a[4].getPosition()[MZ],25);

	//Convert
	DPeakArray<2> list;
	exp.get2DData(list);

	//Tests
	TEST_REAL_EQUAL(list.size(),5);
	DPeakArray<2>::const_iterator it = list.begin();
	TEST_REAL_EQUAL(it->getPosition()[RT],11.1);
	TEST_REAL_EQUAL(it->getPosition()[MZ],5);
	TEST_REAL_EQUAL(it->getIntensity(),47.11);
	++it;
	TEST_REAL_EQUAL(it->getPosition()[RT],11.1);
	TEST_REAL_EQUAL(it->getPosition()[MZ],10);
	TEST_REAL_EQUAL(it->getIntensity(),48.11);
	++it;
	TEST_REAL_EQUAL(it->getPosition()[RT],11.1);
	TEST_REAL_EQUAL(it->getPosition()[MZ],15);
	++it;
	TEST_REAL_EQUAL(it->getPosition()[RT],12.2);
	TEST_REAL_EQUAL(it->getPosition()[MZ],20);
	++it;
	TEST_REAL_EQUAL(it->getPosition()[RT],12.2);
	TEST_REAL_EQUAL(it->getPosition()[MZ],25);
RESULT

CHECK((template<class Container> void set2DData(const Container& cont) throw(Exception::Precondition)))
	MSExperiment<> exp;
	
	// create sample data
	DPeakArray<2> input;
	
	DPeak<2> p1;
	p1.getIntensity()    = 1.0;
	p1.getPosition()[RT] = 2.0;
	p1.getPosition()[MZ] = 3.0;
	input.push_back(p1);
	
	DPeak<2> p2;
	p2.getIntensity()    = 4.0;
	p2.getPosition()[RT] = 5.0;
	p2.getPosition()[MZ] = 6.0;
	input.push_back(p2);
	
	DPeak<2> p3;
	p3.getIntensity()    = 7.5;
	p3.getPosition()[RT] = 8.5;
	p3.getPosition()[MZ] = 9.5;
	input.push_back(p3);
	
	exp.set2DData(input);
	
	// retrieve data again and check for changes
	DPeakArray<2> output;
	
	exp.get2DData(output);
	TEST_EQUAL(output==input,true);
	
	//test precondition
	input.push_back(p1);
	TEST_EXCEPTION(Exception::Precondition,exp.set2DData(input));
	
RESULT

CHECK(([EXTRA] MSExperiment<DRawDataPoint<1> >()))
	MSExperiment<DRawDataPoint<1> > tmp;
	tmp.resize(1);
	tmp[0].getContainer().resize(1);
	tmp[0].getContainer()[0].getPosition()[0] = 47.11;
	TEST_REAL_EQUAL(tmp[0].getContainer()[0].getPosition()[0],47.11)
RESULT

CHECK((CoordinateType getMinMZ() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMinMZ(),numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMaxMZ() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMaxMZ(),-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMinRT() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMinRT(),numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((CoordinateType getMaxRT() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMaxRT(),-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((const std::vector<UnsignedInt>& getMSLevels() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_EQUAL(tmp.getMSLevels().size(),0)
	TEST_REAL_EQUAL(tmp.getDataRange().min()[1],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[1],-numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().min()[0],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[0],-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((UnsignedInt getSize() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_EQUAL(tmp.getSize(),0)
RESULT

CHECK((const AreaType& getDataRange() const))
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getDataRange().min()[1],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[1],-numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().min()[0],numeric_limits<DPosition<2>::CoordinateType>::max())
	TEST_REAL_EQUAL(tmp.getDataRange().max()[0],-numeric_limits<DPosition<2>::CoordinateType>::max())
RESULT

CHECK((virtual void updateRanges()))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setMSLevel(1);
	s.setRetentionTime(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);
	
	s.getContainer().clear();
	s.setMSLevel(1);
	s.setRetentionTime(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRetentionTime(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRetentionTime(50.0);
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

	MSExperiment< DRawDataPoint<1> > tmp2;
	MSSpectrum< DRawDataPoint<1> > s2;
	DRawDataPoint<1> p2;
	
	s2.setRetentionTime(30.0);
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

CHECK((void updateRanges(SignedInt ms_level)))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setMSLevel(1);
	s.setRetentionTime(30.0);
	p.getPosition()[0] = 5.0;
	p.setIntensity(-5.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);
	
	s.getContainer().clear();
	s.setMSLevel(1);
	s.setRetentionTime(40.0);
	p.getPosition()[0] = 7.0;
	p.setIntensity(-7.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRetentionTime(45.0);
	p.getPosition()[0] = 9.0;
	p.setIntensity(-10.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRetentionTime(50.0);
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

	MSExperiment< DRawDataPoint<1> > tmp2;
	MSSpectrum< DRawDataPoint<1> > s2;
	DRawDataPoint<1> p2;
	
	s2.setRetentionTime(30.0);
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

CHECK(ConstAreaIterator areaEndConst() const)
  //Implicitly tested in next test
RESULT

CHECK( ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const )
	DPeakArray<2> plist;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 3.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 10.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 11.0;
	plist.push_back(p1);
	
	MSExperiment<> exp;
	exp.set2DData(plist);

	MSExperiment<>::ConstAreaIterator it = exp.areaBeginConst(0,15,0,15);
	
	TEST_EQUAL(it->getPosition()[0],2.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],3.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],10.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],11.0);
	it++;
	TEST_EQUAL(it==exp.areaEndConst(),true);
RESULT

CHECK(AreaIterator areaEnd())
  //Implicitly tested in next test
RESULT

CHECK( AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) )
	DPeakArray<2> plist;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 3.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 10.0;
	plist.push_back(p1);
	p1.getPosition()[0] = 2.0;
	p1.getPosition()[1] = 11.0;
	plist.push_back(p1);
	
	MSExperiment<> exp;
	exp.set2DData(plist);

	MSExperiment<>::AreaIterator it = exp.areaBegin(0,15,0,15);
	
	TEST_EQUAL(it->getPosition()[0],2.0);
	it->getPosition()[0] = 4711.0;
	TEST_EQUAL(it->getPosition()[0],4711.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],3.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],10.0);
	it++;
	TEST_EQUAL(it->getPosition()[0],11.0);
	it++;
	TEST_EQUAL(it==exp.areaEnd(),true);
RESULT

CHECK((Iterator RTBegin(double rt)))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setRetentionTime(30.0);
	tmp.push_back(s);	
	s.setRetentionTime(40.0);
	tmp.push_back(s);	
	s.setRetentionTime(45.0);
	tmp.push_back(s);	
	s.setRetentionTime(50.0);
	tmp.push_back(s);
	
	MSExperiment< DRawDataPoint<1> >::Iterator it;
	
	it = tmp.RTBegin(20.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((Iterator RTEnd(double rt)))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setRetentionTime(30.0);
	tmp.push_back(s);	
	s.setRetentionTime(40.0);
	tmp.push_back(s);	
	s.setRetentionTime(45.0);
	tmp.push_back(s);	
	s.setRetentionTime(50.0);
	tmp.push_back(s);
	
	MSExperiment< DRawDataPoint<1> >::Iterator it;
	
	it = tmp.RTEnd(20.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((ConstIterator RTBegin(double rt) const))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setRetentionTime(30.0);
	tmp.push_back(s);	
	s.setRetentionTime(40.0);
	tmp.push_back(s);	
	s.setRetentionTime(45.0);
	tmp.push_back(s);	
	s.setRetentionTime(50.0);
	tmp.push_back(s);
	
	MSExperiment< DRawDataPoint<1> >::Iterator it;
	
	it = tmp.RTBegin(20.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTBegin(30.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTBegin(31.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((ConstIterator RTEnd(double rt) const))
	MSExperiment< DRawDataPoint<1> > tmp;
	MSSpectrum< DRawDataPoint<1> > s;
	DRawDataPoint<1> p;
	
	s.setRetentionTime(30.0);
	tmp.push_back(s);	
	s.setRetentionTime(40.0);
	tmp.push_back(s);	
	s.setRetentionTime(45.0);
	tmp.push_back(s);	
	s.setRetentionTime(50.0);
	tmp.push_back(s);
	
	MSExperiment< DRawDataPoint<1> >::Iterator it;
	
	it = tmp.RTEnd(20.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),30.0)
	it = tmp.RTEnd(30.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	it = tmp.RTEnd(31.0);
	TEST_REAL_EQUAL(it->getRetentionTime(),40.0)
	TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true) 
RESULT

CHECK((void sortSpectra(bool sort_mz = true)))
	DPeakArray<2> plist;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 5.0;
	plist.push_back(p1);
		
	DPeak<2> p2;
	p2.getPosition()[0] = 1.0;
	p2.getPosition()[1] = 3.0;
	plist.push_back(p2);
		
	DPeak<2> p3;
	p3.getPosition()[0] = 2.0;
	p3.getPosition()[1] = 14.0;
	plist.push_back(p3);
	
	DPeak<2> p4;
	p4.getPosition()[0] = 2.0;
	p4.getPosition()[1] = 11.0;
	plist.push_back(p4);
	
	MSExperiment<> exp;
	exp.set2DData(plist);
	
	exp.sortSpectra(true);
 	
	TEST_REAL_EQUAL(exp[0][0].getPos(),3.0);
	TEST_REAL_EQUAL(exp[0][1].getPos(),5.0);
	TEST_REAL_EQUAL(exp[1][0].getPos(),11.0);
	TEST_REAL_EQUAL(exp[1][1].getPos(),14.0);
RESULT

CHECK((void reset()))
	DPeakArray<2> plist;
	
	DPeak<2> p;
	p.getPosition()[0] = 1.0;
	p.getPosition()[1] = 5.0;
	plist.push_back(p);
	p.getPosition()[0] = 2.0;
	p.getPosition()[1] = 3.0;
	plist.push_back(p);
		
	MSExperiment<> exp;
	exp.set2DData(plist);
	exp.updateRanges();
	
	exp.reset();
	
	TEST_EQUAL(exp==MSExperiment<>(),true);
RESULT

CHECK((const ExperimentalSettings& getExperimentalSettings() const))
	MSExperiment<> exp;
	exp.setComment("test");
	TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
