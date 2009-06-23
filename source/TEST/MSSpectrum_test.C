// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
// Dummy peak data

Peak1D p1;
p1.setIntensity(1.0f);
p1.setMZ(2.0);

Peak1D p2;
p2.setIntensity(2.0f);
p2.setMZ(10.0);

Peak1D p3;
p3.setIntensity(3.0f);
p3.setMZ(30.0);


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum<>* ptr = 0;
START_SECTION((MSSpectrum()))
	ptr = new MSSpectrum<>();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~MSSpectrum()))
	delete ptr;
END_SECTION

START_SECTION(([EXTRA] MSSpectrum<RichPeak1D>()))
	MSSpectrum<RichPeak1D > tmp;
	RichPeak1D peak;
	peak.getPosition()[0] = 47.11;
	tmp.push_back(peak);
	TEST_EQUAL(tmp.size(),1);
	TEST_REAL_SIMILAR(tmp[0].getMZ(),47.11);	
END_SECTION

/////////////////////////////////////////////////////////////
// Member accessors

START_SECTION((UInt getMSLevel() const))
	MSSpectrum<> spec;
	TEST_EQUAL(spec.getMSLevel(),1)
END_SECTION

START_SECTION((void setMSLevel(UInt ms_level)))
	MSSpectrum<> spec;
	spec.setMSLevel(17);
	TEST_EQUAL(spec.getMSLevel(),17)
END_SECTION

START_SECTION((const String& getName() const))
	MSSpectrum<> s;
	TEST_STRING_EQUAL(s.getName(),"")
END_SECTION

START_SECTION((void setName(const String &name)))
	MSSpectrum<> s;
	s.setName("bla");
	TEST_STRING_EQUAL(s.getName(),"bla")
END_SECTION

START_SECTION((DoubleReal getRT() const ))
  MSSpectrum<> s;
  TEST_REAL_SIMILAR(s.getRT(),-1.0)
END_SECTION

START_SECTION((void setRT(DoubleReal rt)))
  MSSpectrum<> s;
  s.setRT(0.451);
  TEST_REAL_SIMILAR(s.getRT(),0.451)
END_SECTION

START_SECTION((const MetaDataArrays& getMetaDataArrays() const))
   MSSpectrum<> s;
  TEST_EQUAL(s.getMetaDataArrays().size(),0)
END_SECTION

START_SECTION((MetaDataArrays& getMetaDataArrays()))
  MSSpectrum<> s;
  s.getMetaDataArrays().resize(2);
  TEST_EQUAL(s.getMetaDataArrays().size(),2)
END_SECTION

/////////////////////////////////////////////////////////////
// RangeManager

START_SECTION((virtual void updateRanges()))
  MSSpectrum<> s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p1);

  s.updateRanges();
  s.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxInt(),2)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],10)
  TEST_REAL_SIMILAR(s.getMin()[0],2)

  //test with only one peak

	s.clear();
  s.push_back(p1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxInt(),1)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],2)
  TEST_REAL_SIMILAR(s.getMin()[0],2)
END_SECTION


/////////////////////////////////////////////////////////////
// Copy constructor, assignement operator, equality

START_SECTION((MSSpectrum(const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.setMetaValue("label",5.0);
	tmp.setMSLevel(17);
	tmp.setRT(7.0);
	tmp.setName("bla");
	//peaks
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	tmp.push_back(peak);
	
	MSSpectrum<> tmp2(tmp);
	TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
	TEST_EQUAL(tmp2.getMSLevel(), 17)
	TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
	TEST_EQUAL(tmp2.getName(),"bla")
	//peaks
	TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
END_SECTION


START_SECTION((MSSpectrum& operator= (const MSSpectrum& source)))
  MSSpectrum<> tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.setMetaValue("label",5.0);
	tmp.setMSLevel(17);
	tmp.setRT(7.0);
	tmp.setName("bla");
	//peaks
	MSSpectrum<>::PeakType peak;
	peak.getPosition()[0] = 47.11;
	tmp.push_back(peak);
	
	//normal assignment
	MSSpectrum<> tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),1);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("label"), 5.0)
	TEST_EQUAL(tmp2.getMSLevel(), 17)
	TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
	TEST_EQUAL(tmp2.getName(),"bla")
	TEST_EQUAL(tmp2.size(),1);
	TEST_REAL_SIMILAR(tmp2[0].getPosition()[0],47.11);
	
	//Assignment of empty object
	//normal assignment
	tmp2 = MSSpectrum<>();
	TEST_EQUAL(tmp2.getInstrumentSettings().getScanWindows().size(),0);
	TEST_EQUAL(tmp2.metaValueExists("label"), false)
	TEST_EQUAL(tmp2.getMSLevel(),1)
	TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
	TEST_EQUAL(tmp2.getName(),"")
	TEST_EQUAL(tmp2.size(),0);
END_SECTION

START_SECTION((bool operator== (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
	edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.resize(1);
	TEST_EQUAL(edit==empty,false);

	edit = empty;
	edit.setMetaValue("label",String("bla"));
	TEST_EQUAL(empty==edit, false);

	edit.setRT(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getMetaDataArrays().resize(5);
	TEST_EQUAL(empty==edit, false);

	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(empty==edit, true);

	edit = empty;
	edit.push_back(p1);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION((bool operator!= (const MSSpectrum& rhs) const))
  MSSpectrum<> edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
	edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(edit!=empty,true);
	
	edit = empty;
	edit.resize(1);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("bla"));
	TEST_EQUAL(edit!=empty,true);

	edit.setRT(5);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.getMetaDataArrays().resize(5);
	TEST_EQUAL(edit!=empty,true);


	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(edit!=empty,false);

	edit = empty;
	edit.push_back(p1);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(edit!=empty,true);

END_SECTION


/////////////////////////////////////////////////////////////
// Sorting


START_SECTION((void sortByIntensity(bool reverse=false)))
	MSSpectrum<> ds;
	Peak1D p;
	MSSpectrum<>::MetaDataArray tmp;
	tmp.resize(10);
	std::vector<DoubleReal> mzs, intensities;
	intensities.push_back(201); tmp[0] = 420.130f; mzs.push_back(420.130);
	intensities.push_back(60);  tmp[1] = 412.824f; mzs.push_back(412.824);
	intensities.push_back(56);  tmp[2] = 423.269f; mzs.push_back(423.269);
	intensities.push_back(37);  tmp[3] = 415.287f; mzs.push_back(415.287);
	intensities.push_back(34);  tmp[4] = 413.800f; mzs.push_back(413.800);
	intensities.push_back(31);  tmp[5] = 419.113f; mzs.push_back(419.113);
	intensities.push_back(31);  tmp[6] = 416.293f; mzs.push_back(416.293);
	intensities.push_back(31);  tmp[7] = 418.232f; mzs.push_back(418.232);
	intensities.push_back(29);  tmp[8] = 414.301f; mzs.push_back(414.301);
	intensities.push_back(29);  tmp[9] = 412.321f; mzs.push_back(412.321);

	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.sortByIntensity();
	std::vector<DoubleReal> intensities_copy(intensities);
	std::sort(intensities_copy.begin(),intensities_copy.end());
	MSSpectrum<>::iterator it_ds = ds.begin();
	for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
	{
		if(it_ds == ds.end()){ /* fail */ TEST_EQUAL(true,false) }
		TEST_EQUAL(it_ds->getIntensity(), *it);
		++it_ds;
	}
	ds.clear();
	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	intensities_copy = intensities;
	std::sort(intensities_copy.begin(),intensities_copy.end());

	ds.getMetaDataArrays() = std::vector<MSSpectrum<>::MetaDataArray>(3,tmp);
	ds.getMetaDataArrays()[0].setName("a1");
	ds.getMetaDataArrays()[1].setName("a2");
	ds.getMetaDataArrays()[2].setName("a3");

	ds.sortByIntensity();

	TEST_STRING_EQUAL(ds.getMetaDataArrays()[0].getName(),"a1")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[1].getName(),"a2")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[2].getName(),"a3")
	MSSpectrum<>::iterator it1 = ds.begin();
	MSSpectrum<>::MetaDataArray::iterator it2 = ds.getMetaDataArrays()[1].begin();
	TOLERANCE_ABSOLUTE(0.0001)
	for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
	{
		if(it1 != ds.end() && it2 != ds.getMetaDataArrays()[1].end())
		{
			//metadataarray values == mz values
			TEST_EQUAL(it1->getIntensity(), *it);
			TEST_REAL_SIMILAR(*it2 , it1->getMZ());
			++it1;
			++it2;
		}
		else{ /* fail */ TEST_EQUAL(true,false) }
	}
END_SECTION

START_SECTION((void sortByPosition()))
	MSSpectrum<> ds;
	Peak1D p; MSSpectrum<>::MetaDataArray tmp; tmp.resize(10);
	std::vector<DoubleReal> mzs, intensities;
	intensities.push_back(56);  tmp[0] = 56;  mzs.push_back(423.269);
	intensities.push_back(201); tmp[1] = 201; mzs.push_back(420.130);
	intensities.push_back(31);  tmp[2] = 31;  mzs.push_back(419.113);
	intensities.push_back(31);  tmp[3] = 31;  mzs.push_back(418.232);
	intensities.push_back(31);  tmp[4] = 31;  mzs.push_back(416.293);
	intensities.push_back(37);  tmp[5] = 37;  mzs.push_back(415.287);
	intensities.push_back(29);  tmp[6] = 29;  mzs.push_back(414.301);
	intensities.push_back(34);  tmp[7] = 34;  mzs.push_back(413.800);
	intensities.push_back(60);  tmp[8] = 60;  mzs.push_back(412.824);
	intensities.push_back(29);  tmp[9] = 29;  mzs.push_back(412.321);

	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.sortByPosition();
	MSSpectrum<>::iterator it = ds.begin();
	for(std::vector<DoubleReal>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
	{
		if(it == ds.end()){ /* fail */ TEST_EQUAL(true,false) }
		TEST_EQUAL(it->getIntensity(), *rit);
		++it;
	}
	ds.clear();
	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.getMetaDataArrays() = std::vector<MSSpectrum<>::MetaDataArray>(3,tmp);
	ds.getMetaDataArrays()[0].setName("a1");
	ds.getMetaDataArrays()[1].setName("a2");
	ds.getMetaDataArrays()[2].setName("a3");

	ds.sortByPosition();

	TEST_STRING_EQUAL(ds.getMetaDataArrays()[0].getName(),"a1")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[1].getName(),"a2")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[2].getName(),"a3")
	MSSpectrum<>::iterator it1 = ds.begin();
	MSSpectrum<>::MetaDataArray::iterator it2 = ds.getMetaDataArrays()[1].begin();
	for(std::vector<DoubleReal>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
	{
		if(it1 != ds.end() && it2 != ds.getMetaDataArrays()[1].end())
		{
			//metadataarray values == intensity values
			TEST_EQUAL(it1->getIntensity(), *rit);
			TEST_EQUAL(*it2 , *rit);
			++it1;
			++it2;
		}
		else{ /* fail */ TEST_EQUAL(true,false) }
	}

END_SECTION

START_SECTION(bool isSorted() const)
	//make test dataset
	MSSpectrum<> spec;
	Peak1D p;
	p.setIntensity(1.0);
	p.setMZ(1000.0);
	spec.push_back(p);
	
	p.setIntensity(1.0);
	p.setMZ(1001.0);
	spec.push_back(p);
	
	p.setIntensity(1.0);
	p.setMZ(1002.0);
	spec.push_back(p);
	
	TEST_EQUAL(spec.isSorted(),true)
	
	reverse(spec.begin(), spec.end());
	TEST_EQUAL(spec.isSorted(),false)
END_SECTION

/////////////////////////////////////////////////////////////
// Finding peaks or peak ranges

START_SECTION((Iterator MZEnd(CoordinateType mz)))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::Iterator it;

	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::Iterator it;

	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::Iterator it;

	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
	TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::ConstIterator it;

	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(tmp.begin(), 4.5, tmp.begin());
	TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::Iterator it;

	it = tmp.MZEnd(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(tmp.begin(), 5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(tmp.begin(), 4.5, tmp.begin());
	TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::ConstIterator it;

	it = tmp.MZEnd(tmp.begin(), 4.5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(tmp.begin(), 5, tmp.end());
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(tmp.begin(), 4.5, tmp.begin());
	TEST_EQUAL(it->getPosition()[0],tmp.begin()->getPosition()[0])
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::ConstIterator it;

	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const))
	MSSpectrum<> tmp;
	MSSpectrum<>::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.push_back(rdp);

	MSSpectrum<>::ConstIterator it;

	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Size findNearest(CoordinateType mz) const))
	MSSpectrum<> tmp;
	Peak1D p;
	p.setIntensity(29.0f); p.setMZ(412.321); tmp.push_back(p); //0
	p.setIntensity(60.0f); p.setMZ(412.824); tmp.push_back(p); //1
	p.setIntensity(34.0f); p.setMZ(413.8); tmp.push_back(p); //2
	p.setIntensity(29.0f); p.setMZ(414.301); tmp.push_back(p); //3
	p.setIntensity(37.0f); p.setMZ(415.287); tmp.push_back(p); //4
	p.setIntensity(31.0f); p.setMZ(416.293); tmp.push_back(p); //5
	p.setIntensity(31.0f); p.setMZ(418.232); tmp.push_back(p); //6
	p.setIntensity(31.0f); p.setMZ(419.113); tmp.push_back(p); //7
	p.setIntensity(201.0f); p.setMZ(420.13); tmp.push_back(p); //8
	p.setIntensity(56.0f); p.setMZ(423.269); tmp.push_back(p); //9
	p.setIntensity(34.0f); p.setMZ(426.292); tmp.push_back(p); //10
	p.setIntensity(82.0f); p.setMZ(427.28); tmp.push_back(p); //11
	p.setIntensity(87.0f); p.setMZ(428.322); tmp.push_back(p); //12
	p.setIntensity(30.0f); p.setMZ(430.269); tmp.push_back(p); //13
	p.setIntensity(29.0f); p.setMZ(431.246); tmp.push_back(p); //14
	p.setIntensity(42.0f); p.setMZ(432.289); tmp.push_back(p); //15
	p.setIntensity(32.0f); p.setMZ(436.161); tmp.push_back(p); //16
	p.setIntensity(54.0f); p.setMZ(437.219); tmp.push_back(p); //17
	p.setIntensity(40.0f); p.setMZ(439.186); tmp.push_back(p); //18
	p.setIntensity(40); p.setMZ(440.27); tmp.push_back(p); //19
	p.setIntensity(23.0f); p.setMZ(441.224); tmp.push_back(p); //20

	//test outside mass range
	TEST_EQUAL(tmp.findNearest(400.0),0);
	TEST_EQUAL(tmp.findNearest(500.0),20);
	//test mass range borders
	TEST_EQUAL(tmp.findNearest(412.4),0);
	TEST_EQUAL(tmp.findNearest(441.224),20);
	//test inside scan
	TEST_EQUAL(tmp.findNearest(426.29),10);
	TEST_EQUAL(tmp.findNearest(426.3),10);
	TEST_EQUAL(tmp.findNearest(427.2),11);
	TEST_EQUAL(tmp.findNearest(427.3),11);

	//empty spectrum
	MSSpectrum<> tmp2;
	TEST_PRECONDITION_VIOLATED(tmp2.findNearest(427.3));
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



