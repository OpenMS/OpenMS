// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MSExperiment_test.C,v 1.18 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPeakList.h>
#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(MSExperiment, "$Id: MSExperiment_test.C,v 1.18 2006/06/09 14:46:55 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

const int MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ;
const int RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT;

MSExperiment<>* ptr = 0;
CHECK(MSExperiment())
	ptr = new MSExperiment<>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MSExperiment())
	delete ptr;
RESULT

CHECK(MSExperiment(const MSExperiment& source))
  MSExperiment<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setName("Name");
  tmp.resize(1);
  tmp[0].setName("bla4711");
  
  MSExperiment<> tmp2(tmp);
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
  TEST_EQUAL(tmp2[0].getName(),"bla4711");
RESULT

CHECK(MSExperiment& operator= (const MSExperiment& source))
  MSExperiment<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setName("Name");
  tmp.resize(1);
  tmp[0].setName("bla4711");
  
  MSExperiment<> tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getContacts().size(),1);  
  TEST_EQUAL(tmp2.getContacts()[0].getName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
  TEST_EQUAL(tmp2[0].getName(),"bla4711");

  tmp2 = MSExperiment<>();
  TEST_EQUAL(tmp2.getContacts().size(),0);  
  TEST_EQUAL(tmp2.size(),0);
RESULT

CHECK(bool operator== (const MSExperiment& rhs) const)
  MSExperiment<> edit,empty;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty, false);
	
	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK(bool operator!= (const MSExperiment& rhs) const)
  MSExperiment<> edit,empty;
	
	TEST_EQUAL(edit!=empty, false);
	
	edit.getContacts().resize(1);
  TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
  edit.resize(1);
	TEST_EQUAL(edit!=empty, true);
RESULT

CHECK(template<class Container> void get2DData(Container& cont) const)
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
	DPeakArrayNonPolymorphic<2, DRawDataPoint<2> > a;
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
	DPeakList<2> list;
	exp.get2DData(list);

	//Tests
	TEST_REAL_EQUAL(list.size(),5);
	DPeakList<2>::const_iterator it = list.begin();
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

CHECK(template<class Container> void set2DData(Container& cont))
	MSExperiment<> exp;
	
	// create sample data
	DPeakArrayNonPolymorphic<2> input;
	
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
	DPeakArrayNonPolymorphic<2> output;
	
	exp.get2DData(output);
	TEST_EQUAL(output,input);
	
	//test precondition
	input.push_back(p1);
	TEST_EXCEPTION(Exception::Precondition,exp.set2DData(input));
	
RESULT

CHECK(MSExperiment<DRawDataPoint<1> >())
	MSExperiment<DRawDataPoint<1> > tmp;
	tmp.resize(1);
	tmp[0].getContainer().resize(1);
	tmp[0].getContainer()[0].getPosition()[0] = 47.11;
	TEST_REAL_EQUAL(tmp[0].getContainer()[0].getPosition()[0],47.11)
RESULT

CHECK(typename SpectrumType::PeakType::CoordinateType getMinMZ() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMinMZ(),0.0)
RESULT

CHECK(typename SpectrumType::PeakType::CoordinateType getMaxMZ() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMaxMZ(),0.0)
RESULT

CHECK(double getMinInt() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMinInt(),0.0)
RESULT

CHECK(double getMaxInt() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMaxInt(),0.0)
RESULT

CHECK(double getMinRT() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMinRT(),0.0)
RESULT

CHECK(double getMaxRT() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_REAL_EQUAL(tmp.getMaxRT(),0.0)
RESULT

CHECK(const std::vector<UnsignedInt>& getMSLevels() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_EQUAL(tmp.getMSLevels().size(),0)
RESULT

CHECK(UnsignedInt getSize() const)
	MSExperiment<DRawDataPoint<1> > tmp;
	TEST_EQUAL(tmp.getSize(),0)
RESULT

CHECK(void updateRanges())
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
	p.setIntensity(-9.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	

	s.getContainer().clear();
	s.setMSLevel(3);
	s.setRetentionTime(50.0);
	p.getPosition()[0] = 10.0;
	p.setIntensity(-10.0);
	s.getContainer().push_back(p);
	tmp.push_back(s);	
	
	tmp.updateRanges();
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
	TEST_EQUAL(tmp.getMSLevels().size(),2)
	TEST_EQUAL(tmp.getMSLevels()[0],1)
	TEST_EQUAL(tmp.getMSLevels()[1],3)
	TEST_EQUAL(tmp.getSize(),4)
RESULT

CHECK(getPeak())
	DPeakList<2> plist;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
		
	DPeak<2> p2;
	p2.getPosition()[0] = 1.0;
	p2.getPosition()[1] = 3.0;
	plist.push_back(p2);
		
	DPeak<2> p3;
	p3.getPosition()[0] = 2.0;
	p3.getPosition()[1] = 10.0;
	plist.push_back(p3);
	
	DPeak<2> p4;
	p4.getPosition()[0] = 2.0;
	p4.getPosition()[1] = 11.0;
	plist.push_back(p4);
	
	MSExperiment<> exp;
	exp.set2DData(plist);
	
	exp.updateRanges();
	 	
	TEST_EQUAL(exp.getPeak(0).getPosition()[1],2.0);
	
	TEST_EQUAL(exp.getPeak(1).getPosition()[1],3.0);
	
	TEST_EQUAL(exp.getPeak(2).getPosition()[1],10.0);
	
	TEST_EQUAL(exp.getPeak(3).getPosition()[1],11.0);
	

RESULT

CHECK(PeakIterator())
	DPeakList<2> plist;
	
	DPeak<2> p1;
	p1.getPosition()[0] = 1.0;
	p1.getPosition()[1] = 2.0;
	plist.push_back(p1);
		
	DPeak<2> p2;
	p2.getPosition()[0] = 1.0;
	p2.getPosition()[1] = 3.0;
	plist.push_back(p2);
		
	DPeak<2> p3;
	p3.getPosition()[0] = 2.0;
	p3.getPosition()[1] = 10.0;
	plist.push_back(p3);
	
	DPeak<2> p4;
	p4.getPosition()[0] = 2.0;
	p4.getPosition()[1] = 11.0;
	plist.push_back(p4);
	
	MSExperiment<> exp;
	exp.set2DData(plist);
		
	MSExperiment<>::PeakIterator rt_it = exp.RTBegin();
 	
	TEST_EQUAL(rt_it->getPosition()[0],2.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],3.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],10.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],11.0);
	rt_it++;

RESULT

CHECK(sortSpectra())
	DPeakList<2> plist;
	
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
	
	MSExperiment<>::PeakIterator rt_it = exp.RTBegin();
 	
	TEST_EQUAL(rt_it->getPosition()[0],3.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],5.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],11.0);
	rt_it++;
	TEST_EQUAL(rt_it->getPosition()[0],14.0);
	rt_it++;

RESULT

CHECK(const String& getName() const)
	MSExperiment<> exp;
	TEST_EQUAL("",exp.getName());
RESULT

CHECK(void setName(const String& name))
	MSExperiment<> exp;
	exp.setName("bla");
	TEST_EQUAL("bla",exp.getName());
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
