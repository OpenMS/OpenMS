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

#include <string>

#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>
///////////////////////////

START_TEST(DSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef DSpectrum<Peak1D> DSpectrum1;
typedef DSpectrum<Peak2D> DSpectrum2;
typedef DSpectrum<Peak1D> DSpectrum3;
typedef DSpectrum<Peak2D> DSpectrum4;

Peak2D dp2_1;
dp2_1.setIntensity(1);
dp2_1.getPosition()[0] = 2;
dp2_1.getPosition()[1] = 3;

Peak2D dp2_2;
dp2_2.setIntensity(2);
dp2_2.getPosition()[0] = 10;
dp2_2.getPosition()[1] = 12;

Peak2D dp2_3;
dp2_3.setIntensity(3);
dp2_3.getPosition()[0] = 30;
dp2_3.getPosition()[1] = 32;

DSpectrum3* ptr = 0;
START_SECTION((DSpectrum()))
	ptr = new DSpectrum3;
	TEST_NOT_EQUAL(ptr, 0)

	TEST_EQUAL(ptr->getMin(), DSpectrum3::PositionType::max)
	TEST_EQUAL(ptr->getMax(), DSpectrum3::PositionType::min_negative)
	TEST_REAL_SIMILAR(ptr->getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(ptr->getMaxInt(), -numeric_limits<DoubleReal>::max())
END_SECTION

START_SECTION((~DSpectrum()))
	delete ptr;
END_SECTION

START_SECTION((UInt getMSLevel() const))
	DSpectrum2 ds;
	TEST_EQUAL(ds.getMSLevel(),1)
END_SECTION

START_SECTION((void setMSLevel(UInt ms_level)))
	DSpectrum2 ds;
	ds.setMSLevel(17);
	TEST_EQUAL(ds.getMSLevel(),17)
END_SECTION

START_SECTION((virtual void updateRanges()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  s.push_back(dp2_1);

  s.updateRanges();
  s.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxInt(),2)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],10)
  TEST_REAL_SIMILAR(s.getMax()[1],12)
  TEST_REAL_SIMILAR(s.getMin()[0],2)
  TEST_REAL_SIMILAR(s.getMin()[1],3)

  //test with only one peak

	s.clear();
  s.push_back(dp2_1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxInt(),1)
  TEST_REAL_SIMILAR(s.getMinInt(),1)
  TEST_REAL_SIMILAR(s.getMax()[0],2)
  TEST_REAL_SIMILAR(s.getMax()[1],3)
  TEST_REAL_SIMILAR(s.getMin()[0],2)
  TEST_REAL_SIMILAR(s.getMin()[1],3)
END_SECTION

Peak1D p;
p.setIntensity(0.0);
p.getPosition()[0] = 500.0;
Peak1D p2;
p2.setIntensity(100.0);
p2.getPosition()[0] = 1300.0;

START_SECTION((DSpectrum(const DSpectrum& rhs)))
	DSpectrum1 s;
	s.setMetaValue("label",5.0);
	s.push_back(p);
	s.push_back(p2);
	s.getPrecursorPeak().setIntensity(200.0);
	s.setMSLevel(17);
	s.setRT(7.0);
	s.setName("bla");
	s.updateRanges();
	DSpectrum1 s2(s);

	TEST_REAL_SIMILAR(5.0 , (float)s2.getMetaValue("label"))
	TEST_EQUAL(2 , s2.size())
	TEST_REAL_SIMILAR(200.0 , s2.getPrecursorPeak().getIntensity())
	TEST_REAL_SIMILAR(500.0 , s2.getMin()[0])
	TEST_REAL_SIMILAR(1300.0 , s2.getMax()[0])
	TEST_REAL_SIMILAR(0.0 , s2.getMinInt())
	TEST_REAL_SIMILAR(100.0 , s2.getMaxInt())
	TEST_EQUAL(17 , s2.getMSLevel())
	TEST_REAL_SIMILAR(7.0 , s2.getRT())
	TEST_EQUAL("bla",s2.getName())
END_SECTION

START_SECTION((DSpectrum& operator = (const DSpectrum& rhs)))
	DSpectrum1 s;
	s.setMetaValue("label",5.0);
	s.push_back(p);
	s.push_back(p2);
	s.getPrecursorPeak().setIntensity(200.0);
	s.setMSLevel(17);
	s.setRT(7.0);
	s.setName("bla");
	s.updateRanges();
	DSpectrum1 s2;
	s2 = s;

	TEST_REAL_SIMILAR(5.0 , (float)s2.getMetaValue("label"))
	TEST_EQUAL(2 , s2.size())
	TEST_REAL_SIMILAR(200.0 , s2.getPrecursorPeak().getIntensity())
	TEST_REAL_SIMILAR(500.0 , s2.getMin()[0])
	TEST_REAL_SIMILAR(1300.0 , s2.getMax()[0])
	TEST_REAL_SIMILAR(0.0 , s2.getMinInt())
	TEST_REAL_SIMILAR(100.0 , s2.getMaxInt())
	TEST_EQUAL(17 , s2.getMSLevel())
	TEST_REAL_SIMILAR(7.0 , s2.getRT())
	TEST_EQUAL("bla",s2.getName())
END_SECTION

START_SECTION((bool operator == (const DSpectrum& rhs) const))
	DSpectrum1 empty,edit;

	TEST_EQUAL(empty==edit, true);

	edit.setMetaValue("label",String("DSpectrum"));
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(DSpectrum1::ContainerType::value_type());
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getPrecursorPeak().setIntensity(5.5);
	TEST_EQUAL(empty==edit, false);

	edit.setRT(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setRT(-2);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.getPrecursorPeak().getPosition()[0] = 1.5;
	TEST_EQUAL(empty==edit, false);

	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(empty==edit, true);

	edit = empty;
	edit.push_back(p);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty==edit, false);
END_SECTION

START_SECTION((bool operator != (const DSpectrum& rhs) const))
	DSpectrum1 empty,edit;

	TEST_EQUAL(empty!=edit, false);

	edit.setRT(5);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.setRT(-4);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.getPrecursorPeak().getPosition()[0] = 1.5;
	TEST_EQUAL(empty!=edit, true);

	edit.setMetaValue("label",String("DSpectrum"));
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(DSpectrum1::ContainerType::value_type());
	TEST_EQUAL(empty!=edit, true);

	//name is not checked => no change
	edit = empty;
	edit.setName("bla");
	TEST_EQUAL(empty!=edit, false);

	edit = empty;
	edit.push_back(p);
	edit.push_back(p2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty!=edit, true);
END_SECTION

START_SECTION((String getName() const))
	DSpectrum1 s;
	TEST_EQUAL(s.getName(),"")
END_SECTION

START_SECTION((void setName(const String &name)))
	DSpectrum1 s;
	s.setName("bla");
	TEST_EQUAL(s.getName(),"bla")
END_SECTION

//*************************** Tests for MetaInfoInterface ****************************************

START_SECTION(([EXTRA] setMetaValue(index,string)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("type",String("theoretical"));
	TEST_EQUAL(String(ds.getMetaValue("type")),String("theoretical"))
END_SECTION

START_SECTION(([EXTRA] setMetaValue(index,float)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("score",123.234);
	TEST_REAL_SIMILAR(float(ds.getMetaValue("score")),123.234)
END_SECTION

START_SECTION(([EXTRA] setMetaValue(index,int)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("score",123);
	TEST_EQUAL(float(ds.getMetaValue("score")),123)
END_SECTION

START_SECTION(([EXTRA] meta info with copy constructor))
	DSpectrum1 p;
	p.setMetaValue(2,String("bla"));
 	DSpectrum1 p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA] meta info with assignment))
	DSpectrum1 p;
	p.setMetaValue(2,String("bla"));
 	DSpectrum1 p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

//*************************** Tests for Metadata ****************************************


START_SECTION((PrecursorPeakType& getPrecursorPeak()))
  DSpectrum1 sdi;
  DSpectrum1::PrecursorPeakType p1;
  DSpectrum1::PrecursorPeakType p2 = sdi.getPrecursorPeak();
	TEST_EQUAL(p1==p2,true)
  sdi.getPrecursorPeak().setPosition(177);
	p2 = sdi.getPrecursorPeak();
	TEST_REAL_SIMILAR(p2.getPosition()[0],177.)
END_SECTION

START_SECTION((const PrecursorPeakType& getPrecursorPeak() const))
  DSpectrum1 sdi;
  DSpectrum1::PrecursorPeakType p1;
  DSpectrum1::PrecursorPeakType const p2= sdi.getPrecursorPeak();
  TEST_EQUAL(p1==p2,true)
END_SECTION

START_SECTION((void setPrecursorPeak(const PrecursorPeakType& peak)))
  DSpectrum1 sdi;
	DSpectrum1::PrecursorPeakType	p;
  p.setIntensity(47.11);
  sdi.setPrecursorPeak(p);
  TEST_EQUAL(p==sdi.getPrecursorPeak(),true)
  TEST_REAL_SIMILAR(sdi.getPrecursorPeak().getIntensity(),47.11)
END_SECTION

START_SECTION((CoordinateType getRT() const))
  DSpectrum1 sdi;
  TEST_REAL_SIMILAR(sdi.getRT(),-1)
END_SECTION

START_SECTION((void setRT(CoordinateType rt)))
  DSpectrum1 sdi;
  sdi.setRT(0.451);
  TEST_REAL_SIMILAR(sdi.getRT(),0.451)
END_SECTION

START_SECTION((Iterator MZEnd(CoordinateType mz)))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
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

	DSpectrum1::Iterator it;

	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((Iterator MZBegin(CoordinateType mz)))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
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

	DSpectrum1::Iterator it;

	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((ConstIterator MZEnd(CoordinateType mz) const))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
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

	DSpectrum1::ConstIterator it;

	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION((ConstIterator MZBegin(CoordinateType mz) const))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
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

	DSpectrum1::ConstIterator it;

	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
END_SECTION

START_SECTION(( void sortByIntensity() ))
	DSpectrum<Peak1D> ds;
	Peak1D p; DSpectrum<Peak1D>::MetaDataArray tmp; tmp.resize(10);
	std::vector<DoubleReal> mzs, intensities;
	intensities.push_back(201); tmp[0] = 420.130; mzs.push_back(420.130);
	intensities.push_back(60);  tmp[1] = 412.824; mzs.push_back(412.824);
	intensities.push_back(56);  tmp[2] = 423.269; mzs.push_back(423.269);
	intensities.push_back(37);  tmp[3] = 415.287; mzs.push_back(415.287);
	intensities.push_back(34);  tmp[4] = 413.800; mzs.push_back(413.800);
	intensities.push_back(31);  tmp[5] = 419.113; mzs.push_back(419.113);
	intensities.push_back(31);  tmp[6] = 416.293; mzs.push_back(416.293);
	intensities.push_back(31);  tmp[7] = 418.232; mzs.push_back(418.232);
	intensities.push_back(29);  tmp[8] = 414.301; mzs.push_back(414.301);
	intensities.push_back(29);  tmp[9] = 412.321; mzs.push_back(412.321);

	for(UInt i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.sortByIntensity();
	std::vector<DoubleReal> intensities_copy(intensities);
	std::sort(intensities_copy.begin(),intensities_copy.end());
	DSpectrum<Peak1D>::iterator it_ds = ds.begin();
	for(std::vector<DoubleReal>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
	{
		if(it_ds == ds.end()){ /* fail */ TEST_EQUAL(true,false) }
		TEST_EQUAL(it_ds->getIntensity(), *it);
		++it_ds;
	}
	ds.clear();
	for(UInt i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	intensities_copy = intensities;
	std::sort(intensities_copy.begin(),intensities_copy.end());

	ds.getMetaDataArrays() = std::vector<DSpectrum<Peak1D>::MetaDataArray>(3,tmp);
	ds.getMetaDataArrays()[0].setName("a1");
	ds.getMetaDataArrays()[1].setName("a2");
	ds.getMetaDataArrays()[2].setName("a3");
	
	ds.sortByIntensity();
	
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[0].getName(),"a1")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[1].getName(),"a2")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[2].getName(),"a3")
	DSpectrum<Peak1D>::iterator it1 = ds.begin();
	DSpectrum<Peak1D>::MetaDataArray::iterator it2 = ds.getMetaDataArrays()[1].begin();
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

START_SECTION(( void sortByPosition() ))
	DSpectrum<Peak1D> ds;
	Peak1D p; DSpectrum<Peak1D>::MetaDataArray tmp; tmp.resize(10);
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

	for(UInt i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.sortByPosition();
	DSpectrum<Peak1D>::iterator it = ds.begin();
	for(std::vector<DoubleReal>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
	{
		if(it == ds.end()){ /* fail */ TEST_EQUAL(true,false) }
		TEST_EQUAL(it->getIntensity(), *rit);
		++it;
	}
	ds.clear();
	for(UInt i = 0; i < mzs.size(); ++i)
	{
		p.setIntensity(intensities[i]); p.setMZ(mzs[i]);
		ds.push_back(p);
	}
	ds.getMetaDataArrays() = std::vector<DSpectrum<Peak1D>::MetaDataArray>(3,tmp);
	ds.getMetaDataArrays()[0].setName("a1");
	ds.getMetaDataArrays()[1].setName("a2");
	ds.getMetaDataArrays()[2].setName("a3");
	
	ds.sortByPosition();
	
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[0].getName(),"a1")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[1].getName(),"a2")
	TEST_STRING_EQUAL(ds.getMetaDataArrays()[2].getName(),"a3")
	DSpectrum<Peak1D>::iterator it1 = ds.begin();
	DSpectrum<Peak1D>::MetaDataArray::iterator it2 = ds.getMetaDataArrays()[1].begin();
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

	//Peak2D part
	DSpectrum<Peak2D> ds2, ds3;
	DSpectrum<Peak2D>::MetaDataArray tmp2; tmp2.resize(10);
	ds2.push_back(dp2_3);tmp2[0] = 3;
	ds2.push_back(dp2_2);tmp2[1] = 2;
	ds2.push_back(dp2_3);tmp2[2] = 3;
	ds2.push_back(dp2_1);tmp2[3] = 1;
	ds2.push_back(dp2_2);tmp2[4] = 2;
	ds2.push_back(dp2_3);tmp2[5] = 3;
	ds2.push_back(dp2_2);tmp2[6] = 2;
	ds2.push_back(dp2_3);tmp2[7] = 3;
	ds2.push_back(dp2_1);tmp2[8] = 1;
	ds2.push_back(dp2_3);tmp2[9] = 3;
	ds3 = ds2;
	ds2.sortByPosition();
	ds3.getMetaDataArrays() = std::vector<DSpectrum<Peak2D>::MetaDataArray> (3,tmp2);
	ds3.sortByPosition();
	for(UInt i = 0; i < ds2.size(); ++i)
	{
		TEST_EQUAL(ds2[i].getPosition() , ds3[i].getPosition());
		TEST_EQUAL(ds3[i].getIntensity() , ds3.getMetaDataArrays()[1][i]);
	}

	END_SECTION

START_SECTION((UInt findNearest(CoordinateType mz) const  ))
	DSpectrum<Peak1D> tmp;
	Peak1D p;
	p.setIntensity(29); p.setMZ(412.321); tmp.push_back(p); //0
	p.setIntensity(60); p.setMZ(412.824); tmp.push_back(p); //1
	p.setIntensity(34); p.setMZ(413.8); tmp.push_back(p); //2
	p.setIntensity(29); p.setMZ(414.301); tmp.push_back(p); //3
	p.setIntensity(37); p.setMZ(415.287); tmp.push_back(p); //4
	p.setIntensity(31); p.setMZ(416.293); tmp.push_back(p); //5
	p.setIntensity(31); p.setMZ(418.232); tmp.push_back(p); //6
	p.setIntensity(31); p.setMZ(419.113); tmp.push_back(p); //7
	p.setIntensity(201); p.setMZ(420.13); tmp.push_back(p); //8
	p.setIntensity(56); p.setMZ(423.269); tmp.push_back(p); //9
	p.setIntensity(34); p.setMZ(426.292); tmp.push_back(p); //10
	p.setIntensity(82); p.setMZ(427.28); tmp.push_back(p); //11
	p.setIntensity(87); p.setMZ(428.322); tmp.push_back(p); //12
	p.setIntensity(30); p.setMZ(430.269); tmp.push_back(p); //13
	p.setIntensity(29); p.setMZ(431.246); tmp.push_back(p); //14
	p.setIntensity(42); p.setMZ(432.289); tmp.push_back(p); //15
	p.setIntensity(32); p.setMZ(436.161); tmp.push_back(p); //16
	p.setIntensity(54); p.setMZ(437.219); tmp.push_back(p); //17
	p.setIntensity(40); p.setMZ(439.186); tmp.push_back(p); //18
	p.setIntensity(40); p.setMZ(440.27); tmp.push_back(p); //19
	p.setIntensity(23); p.setMZ(441.224); tmp.push_back(p); //20

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
	DSpectrum<Peak1D> tmp2;
	TEST_EXCEPTION(Exception::Precondition, tmp2.findNearest(427.3));
END_SECTION


START_SECTION(const MetaDataArrays& getMetaDataArrays() const)
  DSpectrum<> s1;
  TEST_EQUAL(s1.getMetaDataArrays().size(),0)
END_SECTION

START_SECTION(MetaDataArrays& getMetaDataArrays())
  DSpectrum<> s1;
  s1.getMetaDataArrays().resize(2);
  TEST_EQUAL(s1.getMetaDataArrays().size(),2)
END_SECTION

START_SECTION(([EXTRA] template <typename Container> std::ostream& operator << (std::ostream& os, const DSpectrum<Container>& rhs)))
	NOT_TESTABLE
	DSpectrum<> s1;
	s1.resize(2);
	std::cout << s1 << endl;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
