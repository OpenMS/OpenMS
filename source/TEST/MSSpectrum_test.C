// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
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
MSSpectrum<>* nullPointer = 0;
START_SECTION((MSSpectrum()))
	ptr = new MSSpectrum<>();
	TEST_NOT_EQUAL(ptr, nullPointer)
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

START_SECTION((const FloatDataArrays& getFloatDataArrays() const))
	MSSpectrum<> s;
  TEST_EQUAL(s.getFloatDataArrays().size(),0)
END_SECTION

START_SECTION((FloatDataArrays& getFloatDataArrays()))
  MSSpectrum<> s;
  s.getFloatDataArrays().resize(2);
  TEST_EQUAL(s.getFloatDataArrays().size(),2)
END_SECTION

START_SECTION((const StringDataArrays& getStringDataArrays() const))
	MSSpectrum<> s;
  TEST_EQUAL(s.getStringDataArrays().size(),0)
END_SECTION

START_SECTION((StringDataArrays& getStringDataArrays()))
  MSSpectrum<> s;
  s.getStringDataArrays().resize(2);
  TEST_EQUAL(s.getStringDataArrays().size(),2)
END_SECTION

START_SECTION((const IntegerDataArrays& getIntegerDataArrays() const))
	MSSpectrum<> s;
  TEST_EQUAL(s.getIntegerDataArrays().size(),0)
END_SECTION

START_SECTION((IntegerDataArrays& getIntegerDataArrays()))
  MSSpectrum<> s;
  s.getIntegerDataArrays().resize(2);
  TEST_EQUAL(s.getIntegerDataArrays().size(),2)
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

	s.clear(true);
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
	edit.getFloatDataArrays().resize(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getStringDataArrays().resize(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getIntegerDataArrays().resize(5);
	TEST_EQUAL(empty==edit, false);

	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(empty==edit, true);

	edit = empty;
	edit.push_back(p1);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear(false);
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
	edit.getFloatDataArrays().resize(5);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.getIntegerDataArrays().resize(5);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.getStringDataArrays().resize(5);
	TEST_EQUAL(edit!=empty,true);

	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(edit!=empty,false);

	edit = empty;
	edit.push_back(p1);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear(false);
	TEST_EQUAL(edit!=empty,true);

END_SECTION


/////////////////////////////////////////////////////////////
// Sorting


START_SECTION((void sortByIntensity(bool reverse=false)))
	MSSpectrum<> ds;
	Peak1D p;
	MSSpectrum<>::FloatDataArray float_array;
	MSSpectrum<>::StringDataArray string_array;
	MSSpectrum<>::IntegerDataArray int_array;
	std::vector<DoubleReal> mzs, intensities;
	MSSpectrum<>::IntegerDataArray in_array;
	intensities.push_back(201); mzs.push_back(420.130); float_array.push_back(420.130f); string_array.push_back("420.13"); int_array.push_back(420);
	intensities.push_back(60);  mzs.push_back(412.824); float_array.push_back(412.824f); string_array.push_back("412.82"); int_array.push_back(412);
	intensities.push_back(56);  mzs.push_back(423.269); float_array.push_back(423.269f); string_array.push_back("423.27"); int_array.push_back(423);
	intensities.push_back(37);  mzs.push_back(415.287); float_array.push_back(415.287f); string_array.push_back("415.29"); int_array.push_back(415);
	intensities.push_back(34);  mzs.push_back(413.800); float_array.push_back(413.800f); string_array.push_back("413.80"); int_array.push_back(413);
	intensities.push_back(31);  mzs.push_back(419.113); float_array.push_back(419.113f); string_array.push_back("419.11"); int_array.push_back(419);
	intensities.push_back(31);  mzs.push_back(416.293); float_array.push_back(416.293f); string_array.push_back("416.29"); int_array.push_back(416);
	intensities.push_back(31);  mzs.push_back(418.232); float_array.push_back(418.232f); string_array.push_back("418.23"); int_array.push_back(418);
	intensities.push_back(29);  mzs.push_back(414.301); float_array.push_back(414.301f); string_array.push_back("414.30"); int_array.push_back(414);
	intensities.push_back(29);  mzs.push_back(412.321); float_array.push_back(412.321f); string_array.push_back("412.32"); int_array.push_back(412);

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
		if(it_ds == ds.end())
		{
			TEST_EQUAL(true,false) 
		}
		TEST_EQUAL(it_ds->getIntensity(), *it);
		++it_ds;
	}
	ds.clear(true);
	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	intensities_copy = intensities;
	std::sort(intensities_copy.begin(),intensities_copy.end());

	ds.getFloatDataArrays() = std::vector<MSSpectrum<>::FloatDataArray>(3,float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");

	ds.getStringDataArrays() = std::vector<MSSpectrum<>::StringDataArray>(2, string_array);
	ds.getStringDataArrays()[0].setName("s1");
	ds.getStringDataArrays()[1].setName("s2");

	ds.getIntegerDataArrays() = std::vector<MSSpectrum<>::IntegerDataArray>(1, int_array);
	ds.getIntegerDataArrays()[0].setName("i1");

	ds.sortByIntensity();

	TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
	TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
	TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

	TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
	TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")

	TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")
	
	MSSpectrum<>::iterator it1 = ds.begin();
	MSSpectrum<>::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
	MSSpectrum<>::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
	MSSpectrum<>::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
	TOLERANCE_ABSOLUTE(0.0001)
	for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
	{
		if(it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end() && it4 != ds.getIntegerDataArrays()[0].end())
		{
			//metadataarray values == mz values
			TEST_REAL_SIMILAR(it1->getIntensity(), *it);
			TEST_REAL_SIMILAR(*it2 , it1->getMZ());
			TEST_STRING_EQUAL(*it3 , String::number(it1->getMZ(),2));
			TEST_EQUAL(*it4 , (Int)floor(it1->getMZ()));
			++it1;
			++it2;
			++it3;
			++it4;
		}
		else
		{
			TEST_EQUAL(true,false)
		}
	}
END_SECTION

START_SECTION((void sortByPosition()))
	MSSpectrum<> ds;
	Peak1D p;
	MSSpectrum<>::FloatDataArray float_array;
	MSSpectrum<>::StringDataArray string_array;
	MSSpectrum<>::IntegerDataArray int_array;
	std::vector<DoubleReal> mzs, intensities;
	intensities.push_back(56);  mzs.push_back(423.269); float_array.push_back(56);  string_array.push_back("56");  int_array.push_back(56);  
	intensities.push_back(201); mzs.push_back(420.130); float_array.push_back(201); string_array.push_back("201"); int_array.push_back(201); 
	intensities.push_back(31);  mzs.push_back(419.113); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
	intensities.push_back(31);  mzs.push_back(418.232); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
	intensities.push_back(31);  mzs.push_back(416.293); float_array.push_back(31);  string_array.push_back("31");  int_array.push_back(31);  
	intensities.push_back(37);  mzs.push_back(415.287); float_array.push_back(37);  string_array.push_back("37");  int_array.push_back(37);  
	intensities.push_back(29);  mzs.push_back(414.301); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);  
	intensities.push_back(34);  mzs.push_back(413.800); float_array.push_back(34);  string_array.push_back("34");  int_array.push_back(34);  
	intensities.push_back(60);  mzs.push_back(412.824); float_array.push_back(60);  string_array.push_back("60");  int_array.push_back(60);  
	intensities.push_back(29);  mzs.push_back(412.321); float_array.push_back(29);  string_array.push_back("29");  int_array.push_back(29);  

	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.sortByPosition();
	MSSpectrum<>::iterator it = ds.begin();
	for(std::vector<DoubleReal>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
	{
		if(it == ds.end())
		{
			TEST_EQUAL(true,false)
		}
		TEST_EQUAL(it->getIntensity(), *rit);
		++it;
	}
	ds.clear(true);
	for (Size i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.getFloatDataArrays() = std::vector<MSSpectrum<>::FloatDataArray>(3,float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");
	
	ds.getStringDataArrays() = std::vector<MSSpectrum<>::StringDataArray>(2, string_array);
	ds.getStringDataArrays()[0].setName("s1");
	ds.getStringDataArrays()[1].setName("s2");

	ds.getIntegerDataArrays() = std::vector<MSSpectrum<>::IntegerDataArray>(2, int_array);
	ds.getIntegerDataArrays()[0].setName("i1");
	
	ds.sortByPosition();

	TEST_STRING_EQUAL(ds.getFloatDataArrays()[0].getName(),"f1")
	TEST_STRING_EQUAL(ds.getFloatDataArrays()[1].getName(),"f2")
	TEST_STRING_EQUAL(ds.getFloatDataArrays()[2].getName(),"f3")

	TEST_STRING_EQUAL(ds.getStringDataArrays()[0].getName(),"s1")
	TEST_STRING_EQUAL(ds.getStringDataArrays()[1].getName(),"s2")
	
	TEST_STRING_EQUAL(ds.getIntegerDataArrays()[0].getName(),"i1")

	MSSpectrum<>::iterator it1 = ds.begin();
	MSSpectrum<>::FloatDataArray::iterator it2 = ds.getFloatDataArrays()[1].begin();
	MSSpectrum<>::StringDataArray::iterator it3 = ds.getStringDataArrays()[0].begin();
	MSSpectrum<>::IntegerDataArray::iterator it4 = ds.getIntegerDataArrays()[0].begin();
	for(std::vector<DoubleReal>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
	{
		if(it1 != ds.end() && it2 != ds.getFloatDataArrays()[1].end() && it3 != ds.getStringDataArrays()[0].end())
		{
			//metadataarray values == intensity values
			TEST_REAL_SIMILAR(it1->getIntensity(), *rit);
			TEST_REAL_SIMILAR(*it2 , *rit);
			TEST_STRING_EQUAL(*it3 , String::number(*rit,0));
			TEST_EQUAL(*it4 , (Int)floor(*rit));
			++it1;
			++it2;
			++it3;
			++it4;
		}
		else
		{
			TEST_EQUAL(true,false)
		}
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

START_SECTION(void clear(bool clear_meta_data))
  MSSpectrum<> edit;
  edit.getInstrumentSettings().getScanWindows().resize(1);
	edit.resize(1);
	edit.setMetaValue("label",String("bla"));
	edit.setRT(5);
	edit.setMSLevel(5);
	edit.getFloatDataArrays().resize(5);
	edit.getIntegerDataArrays().resize(5);
	edit.getStringDataArrays().resize(5);

	edit.clear(false);
	TEST_EQUAL(edit.size(),0)
	TEST_EQUAL(edit==MSSpectrum<>(),false)

	edit.clear(true);
	TEST_EQUAL(edit==MSSpectrum<>(),true)
END_SECTION

START_SECTION(([MSSpectrum::RTLess] bool operator()(const MSSpectrum &a, const MSSpectrum &b) const))
  vector< MSSpectrum<> > v;

	MSSpectrum<> sp1;
	sp1.setRT(3.0f);
	v.push_back(sp1);

  MSSpectrum<> sp2;
  sp2.setRT(2.0f);
  v.push_back(sp2);

  MSSpectrum<> sp3;
  sp3.setRT(1.0f);
  v.push_back(sp3);

  std::sort(v.begin(),v.end(), MSSpectrum<>::RTLess());

  TEST_REAL_SIMILAR(v[0].getRT(), 1.0);
  TEST_REAL_SIMILAR(v[1].getRT(), 2.0);
  TEST_REAL_SIMILAR(v[2].getRT(), 3.0);

  ///
  MSSpectrum<> s1;
  s1.setRT(0.451);

  MSSpectrum<> s2;
  s2.setRT(0.5);

  TEST_EQUAL(MSSpectrum<>::RTLess()(s1,s2), true);
  TEST_EQUAL(MSSpectrum<>::RTLess()(s2,s1), false);
  TEST_EQUAL(MSSpectrum<>::RTLess()(s2,s2), false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



