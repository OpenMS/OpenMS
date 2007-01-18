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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(MSExperimentExtern, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

const int MZ = DimensionDescription < LCMS_Tag >::MZ;
const int RT = DimensionDescription < LCMS_Tag >::RT;

MSExperimentExtern<>* ptr = 0;
CHECK(MSExperimentExternExtern())
	ptr = new MSExperimentExtern<>;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MSExperimentExtern())
	delete ptr;
RESULT

CHECK(MSExperimentExtern(const MSExperimentExtern& source))
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

CHECK(MSExperimentExtern& operator= (const MSExperimentExtern& source))
  MSExperimentExtern<> tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.setBufferSize(1);
  tmp.updateBuffer();
  MSSpectrum<> spec;
  DPeak<1> p;
  p.setPos(5.0);
  spec.push_back(p);
  p.setPos(10.0);
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

CHECK(bool operator== (const MSExperimentExtern& rhs) const)
  	MSExperimentExtern<> edit,empty;
	
	TEST_EQUAL(edit==empty, true);
	
	edit.getContacts().resize(1);
  	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
  	edit.setBufferSize(1);
  	edit.updateBuffer();
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK(bool operator!= (const MSExperimentExtern& rhs) const)
  	MSExperimentExtern<> edit,empty;
	
	TEST_EQUAL(edit!=empty, false);
	
	edit.getContacts().resize(1);
  	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setBufferSize(1);
  	edit.updateBuffer();
	TEST_EQUAL(edit!=empty, true);
RESULT

CHECK(template<class Container> void get2DData(Container& cont) const)
	MSExperimentExtern<> exp;
	MSExperimentExtern<>::SpectrumType spec;
	MSExperimentExtern<>::PeakType peak;

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

CHECK(template<class Container> void set2DData(Container& cont))
	MSExperimentExtern<> exp;
	
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
RESULT


CHECK(UnsignedInt getSize() const)
	MSExperimentExtern<DRawDataPoint<1> > tmp;
	TEST_EQUAL(tmp.getSize(),0);
	
	DRawDataPoint<1> p1;
	MSSpectrum<DRawDataPoint<1> > spec;
	spec.push_back(p1);
	spec.push_back(p1);
	spec.push_back(p1);
	
	tmp.push_back(spec);
	tmp.updateRanges();
	TEST_EQUAL(tmp.getSize(),3);
		
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
