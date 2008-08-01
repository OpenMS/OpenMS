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
///////////////////////////

START_TEST(DSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef DSpectrum<DPeakArray<DPeak<1> > > DSpectrum1;
typedef DSpectrum<DPeakArray<DPeak<2> > > DSpectrum2;
typedef DSpectrum<DPeakArray<DPeak<3> > > DSpectrum3;
typedef DSpectrum<DPeakArray<DPeak<4> > > DSpectrum4;

DPeak<2> dp2_1;
dp2_1.setIntensity(1);
dp2_1.getPosition()[0] = 2;
dp2_1.getPosition()[1] = 3;

DPeak<2> dp2_2;
dp2_2.setIntensity(2);
dp2_2.getPosition()[0] = 10;
dp2_2.getPosition()[1] = 12;

DPeak<2> dp2_3;
dp2_3.setIntensity(3);
dp2_3.getPosition()[0] = 30;
dp2_3.getPosition()[1] = 32;

DSpectrum3* ptr = 0;
CHECK((DSpectrum()))
	ptr = new DSpectrum3;
	TEST_NOT_EQUAL(ptr, 0)

	TEST_EQUAL(ptr->getMin(), DSpectrum3::PositionType::max)
	TEST_EQUAL(ptr->getMax(), DSpectrum3::PositionType::min_negative)
	TEST_REAL_EQUAL(ptr->getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_EQUAL(ptr->getMaxInt(), -numeric_limits<DoubleReal>::max())
RESULT

CHECK((~DSpectrum()))
	delete ptr;
RESULT

CHECK((UInt getMSLevel() const))
	DSpectrum2 ds;
	TEST_EQUAL(ds.getMSLevel(),1)
RESULT

CHECK((void setMSLevel(UInt ms_level)))
	DSpectrum2 ds;
	ds.setMSLevel(17);
	TEST_EQUAL(ds.getMSLevel(),17)
RESULT

CHECK((virtual void updateRanges()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  s.push_back(dp2_1);
  
  s.updateRanges();
  s.updateRanges(); //second time to check the initialization
    
  TEST_REAL_EQUAL(s.getMaxInt(),2)
  TEST_REAL_EQUAL(s.getMinInt(),1)
  TEST_REAL_EQUAL(s.getMax()[0],10)
  TEST_REAL_EQUAL(s.getMax()[1],12)
  TEST_REAL_EQUAL(s.getMin()[0],2)
  TEST_REAL_EQUAL(s.getMin()[1],3)
  
  //test with only one peak

	s.clear();
  s.push_back(dp2_1);  
  s.updateRanges();
  TEST_REAL_EQUAL(s.getMaxInt(),1)
  TEST_REAL_EQUAL(s.getMinInt(),1)
  TEST_REAL_EQUAL(s.getMax()[0],2)
  TEST_REAL_EQUAL(s.getMax()[1],3)
  TEST_REAL_EQUAL(s.getMin()[0],2)
  TEST_REAL_EQUAL(s.getMin()[1],3)  
RESULT

DPeak<1> p;
p.setIntensity(0.0);
p.getPosition()[0] = 500.0;
DPeak<1> p2;
p2.setIntensity(100.0);
p2.getPosition()[0] = 1300.0;

CHECK((DSpectrum(const DSpectrum& rhs)))
	DSpectrum1 s;
	s.setMetaValue("label",5.0);
	s.getContainer().push_back(p);
	s.getContainer().push_back(p2);
	s.getPrecursorPeak().setIntensity(200.0);
	s.setMSLevel(17);
	s.setRT(7.0);
	s.setName("bla");
	s.updateRanges();
	DSpectrum1 s2(s);

	TEST_REAL_EQUAL(5.0 , (float)s2.getMetaValue("label"))
	TEST_EQUAL(2 , s2.getContainer().size())
	TEST_REAL_EQUAL(200.0 , s2.getPrecursorPeak().getIntensity())
	TEST_REAL_EQUAL(500.0 , s2.getMin()[0])
	TEST_REAL_EQUAL(1300.0 , s2.getMax()[0])
	TEST_REAL_EQUAL(0.0 , s2.getMinInt())
	TEST_REAL_EQUAL(100.0 , s2.getMaxInt())
	TEST_EQUAL(17 , s2.getMSLevel())
	TEST_REAL_EQUAL(7.0 , s2.getRT())
	TEST_EQUAL("bla",s2.getName())
RESULT

CHECK((DSpectrum& operator = (const DSpectrum& rhs)))
	DSpectrum1 s;
	s.setMetaValue("label",5.0);
	s.getContainer().push_back(p);
	s.getContainer().push_back(p2);
	s.getPrecursorPeak().setIntensity(200.0);
	s.setMSLevel(17);
	s.setRT(7.0);
	s.setName("bla");
	s.updateRanges();
	DSpectrum1 s2;
	s2 = s;

	TEST_REAL_EQUAL(5.0 , (float)s2.getMetaValue("label"))
	TEST_EQUAL(2 , s2.getContainer().size())
	TEST_REAL_EQUAL(200.0 , s2.getPrecursorPeak().getIntensity())
	TEST_REAL_EQUAL(500.0 , s2.getMin()[0])
	TEST_REAL_EQUAL(1300.0 , s2.getMax()[0])
	TEST_REAL_EQUAL(0.0 , s2.getMinInt())
	TEST_REAL_EQUAL(100.0 , s2.getMaxInt())
	TEST_EQUAL(17 , s2.getMSLevel())
	TEST_REAL_EQUAL(7.0 , s2.getRT())
	TEST_EQUAL("bla",s2.getName())
RESULT

CHECK((bool operator == (const DSpectrum& rhs) const))
	DSpectrum1 empty,edit;
	
	TEST_EQUAL(empty==edit, true);
	
	edit.setMetaValue("label",String("DSpectrum"));
	TEST_EQUAL(empty==edit, false);
	
	edit = empty;
	edit.setMSLevel(5);
	TEST_EQUAL(empty==edit, false);
	
	edit = empty;
	edit.getContainer().push_back(DSpectrum1::ContainerType::value_type());
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
RESULT

CHECK((bool operator != (const DSpectrum& rhs) const))
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
	edit.getContainer().push_back(DSpectrum1::ContainerType::value_type());
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
RESULT

CHECK((String getName() const))
	DSpectrum1 s;
	TEST_EQUAL(s.getName(),"")
RESULT

CHECK((void setName(const String &name)))
	DSpectrum1 s;
	s.setName("bla");
	TEST_EQUAL(s.getName(),"bla")
RESULT

//*************************** Container interface tests *********************************

CHECK((const ContainerType& getContainer() const))
	DSpectrum4 ds;
	TEST_EQUAL(ds.getContainer().empty(), true)
	TEST_EQUAL(ds.getContainer().size(), 0)
RESULT

CHECK((reference operator[] (size_type n)))
	DSpectrum2 ds;
	ds.getContainer().resize(1);
	ds[0].setIntensity(5.0);
	TEST_EQUAL(ds.getContainer()[0].getIntensity(),5.0);
RESULT

CHECK((const_reference operator[] (size_type n) const))
	DSpectrum2 ds;
	ds.getContainer().resize(1);
	ds.getContainer()[0].setIntensity(5.0);
	TEST_EQUAL(ds[0].getIntensity(),5.0);
RESULT

CHECK((ContainerType& getContainer()))
	DSpectrum4 ds;
	TEST_EQUAL(ds.getContainer().empty(), true)
	TEST_EQUAL(ds.getContainer().size(), 0)
	TEST_EQUAL(ds.empty(), true)
	TEST_EQUAL(ds.size(), 0)

	DPeak<4> p;
	p.getPosition()[0] = 0.0;
	p.getPosition()[1] = 1.1;
	p.getPosition()[2] = 2.2;
	p.getPosition()[3] = 3.3;
	p.setIntensity(15.0);
	ds.getContainer().push_back(p);
	TEST_EQUAL(ds.getContainer().empty(), false)
	TEST_EQUAL(ds.getContainer().size(), 1)
	TEST_EQUAL(ds.empty(), false)
	TEST_EQUAL(ds.size(), 1)

	DPeak<4> q;
	q.getPosition()[0] = 0.0;
	q.getPosition()[1] = 1.1;
	q.getPosition()[2] = 2.2;
	q.getPosition()[3] = 3.3;
	q.setIntensity(15.0);
	ds.getContainer().push_back(q);
	TEST_EQUAL(ds.getContainer().empty(), false)
	TEST_EQUAL(ds.getContainer().size(), 2)
	TEST_EQUAL(ds.empty(), false)
	TEST_EQUAL(ds.size(), 2)

	ds.getContainer().clear();
	TEST_EQUAL(ds.getContainer().empty(), true)
	TEST_EQUAL(ds.getContainer().size(), 0)	
	TEST_EQUAL(ds.empty(), true)
	TEST_EQUAL(ds.size(), 0)
RESULT

CHECK((void setContainer(const ContainerType& container)))
  DSpectrum2::ContainerType c;
  c.push_back(dp2_2);
  c.push_back(dp2_1);
  
  DSpectrum2 s;
  s.getContainer().push_back(dp2_1);
  s.setContainer(c);
  TEST_EQUAL(s.getContainer().size(),2);
  TEST_REAL_EQUAL(s.getContainer()[0].getIntensity(),2);
	TEST_REAL_EQUAL(s.getContainer()[1].getIntensity(),1);
RESULT

CHECK((UInt size() const))
	DSpectrum2 s;
	TEST_EQUAL(s.size(), 0)
RESULT

CHECK((void push_back( const value_type& val )))
  DSpectrum2 s;
  s.push_back(dp2_1);
	TEST_EQUAL(s.size(), 1)
  s.push_back(dp2_1);
	TEST_EQUAL(s.size(), 2)
  s.push_back(dp2_1);
	TEST_EQUAL(s.size(), 3)
RESULT

CHECK((void clear()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.clear();
	TEST_EQUAL(s.size(), 0)
RESULT

CHECK((bool empty() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
	TEST_EQUAL(s.empty(), false)
  s.clear();
	TEST_EQUAL(s.empty(), true)
RESULT

CHECK((ConstIterator begin() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	TEST_REAL_EQUAL(s.begin()->getIntensity(),1.0)
	TEST_REAL_EQUAL((++(s.begin()))->getIntensity(),2.0)
RESULT

CHECK((ConstIterator end() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	TEST_REAL_EQUAL((--(s.end()))->getIntensity(),2.0)
RESULT

CHECK((Iterator begin()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	s.begin()->setIntensity(4711);
	TEST_REAL_EQUAL(s.begin()->getIntensity(),4711)
RESULT

CHECK((Iterator end()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	(--(s.end()))->setIntensity(4711);
	TEST_REAL_EQUAL((--(s.end()))->getIntensity(),4711)
RESULT

CHECK((size_type max_size() const))
  DSpectrum2 s;
  TEST_EQUAL(s.max_size()>0,true)
RESULT

CHECK((void swap(ContainerType& rhs)))
  DSpectrum2::ContainerType c;
  c.push_back(dp2_2);
  c.push_back(dp2_1);
  
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.swap(c);
  TEST_EQUAL(s.getContainer().size(),2);
  TEST_REAL_EQUAL(s.getContainer()[0].getIntensity(),2);
	TEST_REAL_EQUAL(s.getContainer()[1].getIntensity(),1);
  TEST_EQUAL(c.size(),1);
  TEST_REAL_EQUAL(c[0].getIntensity(),1);
RESULT

CHECK((ConstReverseIterator rbegin() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	TEST_REAL_EQUAL(s.rbegin()->getIntensity(),2)
RESULT

CHECK((ConstReverseIterator rend() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	TEST_REAL_EQUAL((--(s.rend()))->getIntensity(),1)
RESULT

CHECK((ReverseIterator rbegin()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	s.rbegin()->setIntensity(4711);
	TEST_REAL_EQUAL(s.rbegin()->getIntensity(),4711)
RESULT

CHECK((ReverseIterator rend()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	(--(s.rend()))->setIntensity(4711);
	TEST_REAL_EQUAL((--(s.rend()))->getIntensity(),4711)
RESULT

CHECK((Iterator insert( Iterator loc, const value_type& val )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	s.insert(++(s.begin()),dp2_3);
	TEST_REAL_EQUAL(s.begin()->getIntensity(),1)
	TEST_REAL_EQUAL((++(s.begin()))->getIntensity(),3)
	TEST_REAL_EQUAL((++(++(s.begin())))->getIntensity(),2)
RESULT

CHECK((void insert( iterator loc, size_type num, const value_type& val )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
	s.insert(++(s.begin()),4,dp2_3);
	TEST_REAL_EQUAL(s.begin()->getIntensity(),1)
	TEST_REAL_EQUAL((++(s.begin()))->getIntensity(),3)
	TEST_REAL_EQUAL((--(--(s.end())))->getIntensity(),3)
	TEST_REAL_EQUAL((--(s.end()))->getIntensity(),2)
RESULT

CHECK((template<class InputIterator> void insert( iterator loc, InputIterator start, InputIterator end )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  
  DSpectrum2::ContainerType c;
  c.push_back(dp2_3);
  c.push_back(dp2_1);
  
  s.insert(++(s.begin()),c.begin(),c.end());
	TEST_REAL_EQUAL(s.begin()->getIntensity(),1)
	TEST_REAL_EQUAL((++(s.begin()))->getIntensity(),3)
	TEST_REAL_EQUAL((--(--(s.end())))->getIntensity(),1)
	TEST_REAL_EQUAL((--(s.end()))->getIntensity(),2)  
RESULT

CHECK((value_type& front()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
  s.front().setIntensity(5);
  TEST_REAL_EQUAL(s.front().getIntensity(),5)
RESULT

CHECK((const value_type& front() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
RESULT

CHECK((value_type& back()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  TEST_REAL_EQUAL(s.back().getIntensity(),2)
  s.back().setIntensity(5);
  TEST_REAL_EQUAL(s.back().getIntensity(),5)
RESULT

CHECK((const value_type& back() const))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  TEST_REAL_EQUAL(s.back().getIntensity(),2)
RESULT

CHECK((Iterator erase( iterator loc )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  s.push_back(dp2_3);
  s.erase(++(s.begin()));
  TEST_EQUAL(s.size(),2);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
  TEST_REAL_EQUAL(s.back().getIntensity(),3)
RESULT

CHECK((Iterator erase( iterator start, iterator end )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  s.push_back(dp2_2);
  s.push_back(dp2_2);
  s.push_back(dp2_3);
  s.erase(++(s.begin()), --(s.end()));
  TEST_EQUAL(s.size(),2);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
  TEST_REAL_EQUAL(s.back().getIntensity(),3)
RESULT

CHECK((void pop_back()))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  s.push_back(dp2_3);
  s.pop_back();
  TEST_REAL_EQUAL(s.back().getIntensity(),2)
  s.pop_back();
  TEST_REAL_EQUAL(s.back().getIntensity(),1)
RESULT

CHECK((void assign( size_type num, const value_type& val )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_1);
  s.push_back(dp2_3);
  s.assign(2, dp2_2);
  TEST_EQUAL(s.size(),2);
  TEST_REAL_EQUAL(s.front().getIntensity(),2)
  TEST_REAL_EQUAL(s.back().getIntensity(),2)  
RESULT

CHECK((template<class InputIterator> void assign( InputIterator start, InputIterator end )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  
  DSpectrum2::ContainerType c;
  c.push_back(dp2_3);
  c.push_back(dp2_1);
  
  s.assign(c.begin(),c.end());
  TEST_EQUAL(s.size(),2);
  TEST_REAL_EQUAL(s.front().getIntensity(),3)
  TEST_REAL_EQUAL(s.back().getIntensity(),1)  
RESULT

CHECK((void resize( size_type num, const value_type& val = value_type() )))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  
  s.resize(1);
  TEST_EQUAL(s.size(),1);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
  TEST_REAL_EQUAL(s.back().getIntensity(),1)  

  s.resize(3,dp2_3);
  TEST_EQUAL(s.size(),3);
  TEST_REAL_EQUAL(s.front().getIntensity(),1)
  TEST_REAL_EQUAL((++(s.begin()))->getIntensity(),3)  
  TEST_REAL_EQUAL(s.back().getIntensity(),3)  
RESULT

CHECK((bool operator<(const DSpectrum& rhs)))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  DSpectrum2 s2(s);

	TEST_EQUAL(s < s2, false)
	s2.push_back(dp2_2);
	TEST_EQUAL(s < s2 , true)
RESULT

CHECK((bool operator>(const DSpectrum& rhs)))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  DSpectrum2 s2(s);

	TEST_EQUAL(s > s2, false)
	s.push_back(dp2_2);
	TEST_EQUAL(s > s2 , true)
RESULT

CHECK((bool operator<=(const DSpectrum& rhs)))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  DSpectrum2 s2(s);

	TEST_EQUAL(s <= s2, true)
	s.push_back(dp2_2);
	TEST_EQUAL(s <= s2 , false)
RESULT

CHECK((bool operator>=(const DSpectrum& rhs)))
  DSpectrum2 s;
  s.push_back(dp2_1);
  s.push_back(dp2_2);
  DSpectrum2 s2(s);

	TEST_EQUAL(s >= s2, true)
	s2.push_back(dp2_2);
	TEST_EQUAL(s >= s2 , false)
RESULT

//*************************** Tests for MetaInfoInterface ****************************************

CHECK(([EXTRA] setMetaValue(index,string)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("type",String("theoretical"));
	TEST_EQUAL(String(ds.getMetaValue("type")),String("theoretical"))
RESULT

CHECK(([EXTRA] setMetaValue(index,float)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("score",123.234);
	TEST_REAL_EQUAL(float(ds.getMetaValue("score")),123.234)
RESULT

CHECK(([EXTRA] setMetaValue(index,int)/getMetaValue(index)))
	DSpectrum1 ds;
	ds.setMetaValue("score",123);
	TEST_EQUAL(float(ds.getMetaValue("score")),123)
RESULT

CHECK(([EXTRA] meta info with copy constructor))
	DSpectrum1 p;
	p.setMetaValue(2,String("bla"));
 	DSpectrum1 p2(p);
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

CHECK(([EXTRA] meta info with assignment))
	DSpectrum1 p;
	p.setMetaValue(2,String("bla"));
 	DSpectrum1 p2 = p;
	TEST_EQUAL(p.getMetaValue(2), "bla")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
 	p.setMetaValue(2,String("bluff"));
	TEST_EQUAL(p.getMetaValue(2), "bluff")
	TEST_EQUAL(p2.getMetaValue(2), "bla")
RESULT

//*************************** Tests for Metadata ****************************************


CHECK((PrecursorPeakType& getPrecursorPeak()))
  DSpectrum1 sdi;
  DSpectrum1::PrecursorPeakType p1;
  DSpectrum1::PrecursorPeakType p2 = sdi.getPrecursorPeak();
	TEST_EQUAL(p1==p2,true)
  sdi.getPrecursorPeak().setPosition(177);
	p2 = sdi.getPrecursorPeak();
	TEST_REAL_EQUAL(p2.getPosition()[0],177.)
RESULT

CHECK((const PrecursorPeakType& getPrecursorPeak() const))
  DSpectrum1 sdi;
  DSpectrum1::PrecursorPeakType p1;
  DSpectrum1::PrecursorPeakType const p2= sdi.getPrecursorPeak();
  TEST_EQUAL(p1==p2,true)
RESULT

CHECK((void setPrecursorPeak(const PrecursorPeakType& peak)))
  DSpectrum1 sdi;
	DSpectrum1::PrecursorPeakType	p;
  p.setIntensity(47.11);
  sdi.setPrecursorPeak(p);
  TEST_EQUAL(p==sdi.getPrecursorPeak(),true)
  TEST_REAL_EQUAL(sdi.getPrecursorPeak().getIntensity(),47.11)
RESULT

CHECK((CoordinateType getRT() const))
  DSpectrum1 sdi;
  TEST_REAL_EQUAL(sdi.getRT(),-1)
RESULT

CHECK((void setRT(CoordinateType rt)))
  DSpectrum1 sdi;
  sdi.setRT(0.451);
  TEST_REAL_EQUAL(sdi.getRT(),0.451)
RESULT

CHECK((Iterator MZEnd(CoordinateType mz)))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	DSpectrum1::Iterator it;
	
	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((Iterator MZBegin(CoordinateType mz)))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	DSpectrum1::Iterator it;
	
	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((ConstIterator MZEnd(CoordinateType mz) const))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	DSpectrum1::ConstIterator it;
	
	it = tmp.MZBegin(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.0);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZBegin(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((ConstIterator MZBegin(CoordinateType mz) const))
	DSpectrum1 tmp;
	DSpectrum1::PeakType rdp;
	rdp.getPosition()[0] = 1.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 2.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 3.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 4.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 5.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 6.0;
	tmp.getContainer().push_back(rdp);
	rdp.getPosition()[0] = 7.0;
	tmp.getContainer().push_back(rdp);
	
	DSpectrum1::ConstIterator it;
	
	it = tmp.MZEnd(4.5);
	TEST_EQUAL(it->getPosition()[0],5.0)
	it = tmp.MZEnd(5.0);
	TEST_EQUAL(it->getPosition()[0],6.0)
	it = tmp.MZEnd(5.5);
	TEST_EQUAL(it->getPosition()[0],6.0)
RESULT

CHECK((UInt findNearest(CoordinateType mz) const  ))
	DSpectrum< DPeakArray<Peak1D> > tmp;
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
	DSpectrum< DPeakArray<Peak1D> > tmp2;
	TEST_EXCEPTION(Exception::Precondition, tmp2.findNearest(427.3));
RESULT


CHECK(const MetaDataArrays& getMetaDataArrays() const)
  DSpectrum<> s1;
  TEST_EQUAL(s1.getMetaDataArrays().size(),0)
RESULT

CHECK(MetaDataArrays& getMetaDataArrays())
  DSpectrum<> s1;
  s1.getMetaDataArrays().resize(2);
  TEST_EQUAL(s1.getMetaDataArrays().size(),2)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
