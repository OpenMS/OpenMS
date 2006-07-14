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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/DPeakList.h>
#include <OpenMS/KERNEL/DPickedPeak.h>
///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

class Labeled1DPeak: public DPickedPeak<1>
{
	public:
	Labeled1DPeak(): label_() {}
	Labeled1DPeak(string label): label_(label) {}
	~Labeled1DPeak() {};
	string& getLabel() {return label_;}
	void setLabel(const string& label) {label_ = label;}
	virtual DPeak<1>* clone() const
	{
		DPeak<1>* tmp = new Labeled1DPeak(*this);
		return tmp;
	}
	
	virtual bool operator == (const DPeak<1>& rhs) const
	{
		if (typeid(*this) != typeid(rhs))
		{
			return false;
		}
		return (this->label_ == dynamic_cast<const Labeled1DPeak&>(rhs).label_) && (DPickedPeak<1>::operator==(rhs));
	}
	
	virtual bool operator != (const DPeak<1>& rhs) const
	{
		return !(operator==(rhs));
	}		

	protected:
	string label_;
};

class Labeled2DPeak: public DPickedPeak<2>
{
	public:
	Labeled2DPeak(): label_() {}
	Labeled2DPeak(string label): label_(label) {}
	~Labeled2DPeak() {};
	string& label() {return label_;}
	virtual DPeak<2>* clone() const
	{
		DPeak<2>* tmp = new Labeled2DPeak(*this);
		return tmp;
	}
	virtual bool operator == (const DPeak<2>& rhs) const
	{
		if (typeid(*this) != typeid(rhs))
		{
			return false;
		}
		return (this->label_ == dynamic_cast<const Labeled2DPeak&>(rhs).label_) && (DPickedPeak<2>::operator==(rhs));
	}		
	virtual bool operator != (const DPeak<2>& rhs) const
	{
		return !(operator==(rhs));
	}		

	protected:
	string label_;
};

class SameType
{
	public:
	SameType(){}
	~SameType() {}
	
	bool operator () (const  DPickedPeak<2>& a, const  DPickedPeak<2>& b)
	{
		return (typeid(a) ==typeid(b));
	}
};

START_TEST(DPeakList<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DPeakList<2, DPickedPeak<2> >* pl_ptr = 0;
CHECK(DPeakList())
	pl_ptr = new DPeakList<2, DPickedPeak<2> >;
	TEST_NOT_EQUAL(pl_ptr, 0)
	TEST_EQUAL(pl_ptr->size(), 0)
RESULT

CHECK(~DPeakList())
	delete pl_ptr;
RESULT


CHECK(DPeakList())
	std::vector<DPeak<2>*> l;
	l.push_back(0);
	l.push_back(0);
	l.push_back(0);
	std::list<DPeak<2>*> l2;
	l2.insert(l2.begin(),l.begin(),l.end());
	TEST_EQUAL(l2.size(),3)
RESULT

CHECK(DPeakList(const DPeakList& p))
	DPeakList<4 ,DPickedPeak<4> > pl;
	DPickedPeak<4> peak;
	peak.getIntensity() = 1.0;
  pl.push_back(peak);
	peak.getIntensity() = 2.0;
  pl.push_back(peak);
  
  DPeakList<4 ,DPickedPeak<4> > pl2(pl);
	TEST_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2.begin()->getIntensity(), 1.0)
	TEST_REAL_EQUAL((++pl2.begin())->getIntensity(), 2.0)
RESULT

DPeakList<2, DPickedPeak<2> > pl;

CHECK(bool empty() const)
	TEST_EQUAL(pl.empty(), true)
RESULT

DPickedPeak<2> peak1;
peak1.getPosition()[0] = 2.0;
peak1.getPosition()[1] = 3.0;
peak1.getIntensity() = 1.0;

DPickedPeak<2> peak2;
peak2.getPosition()[0] = 0.0;
peak2.getPosition()[1] = 2.5;
peak2.getIntensity() = 0.5;

DPickedPeak<2> peak3;
peak3.getPosition()[0] = 10.5;
peak3.getPosition()[1] = 0.0;
peak3.getIntensity() = 0.01;

CHECK(size_type size() const)
	TEST_EQUAL(pl.size(), 0)
	
	pl.push_back(peak1);
	TEST_EQUAL(pl.size(), 1)

	pl.push_back(peak2);
	TEST_EQUAL(pl.size(), 2)

	pl.push_back(peak3);
	TEST_EQUAL(pl.size(), 3)
RESULT

CHECK(bool empty() const)
	TEST_EQUAL(pl.empty(), false)
RESULT

CHECK([EXTRA] ConstIterator begin() const)
	const DPeakList<2, DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	TEST_REAL_EQUAL(c_pl.begin()->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[1], peak1.getPosition()[1])
RESULT

CHECK([EXTRA] ConstIterator end() const)
	const DPeakList<2, DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	bool result = (c_pl.begin() == c_pl.end());
	TEST_EQUAL(result, false)
	const DPeakList<2, DPickedPeak<2> > empty;
	result = (empty.begin() == empty.end());
	TEST_EQUAL(result, true)
	std::vector< DPickedPeak<2> > v(c_pl.size());
	std::copy(c_pl.begin(), c_pl.end(), v.begin());
	TEST_EQUAL(v.size(), 3)
	ABORT_IF(v.size() != 3)
	TEST_REAL_EQUAL(v[0].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(DPeakList& operator = (const DPeakList& rhs))
	DPeakList<2, DPickedPeak<2> > copy_of_pl;
	TEST_EQUAL(copy_of_pl.size(), 0)
	copy_of_pl = pl;
	TEST_EQUAL(copy_of_pl.size(), 3)
	
	std::vector< DPickedPeak<2> > v(copy_of_pl.size());
	std::copy(copy_of_pl.begin(), copy_of_pl.end(), v.begin());
	TEST_EQUAL(v.size(), 3)
	ABORT_IF(v.size() != 3)
	TEST_REAL_EQUAL(v[0].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(void sortByIntensity())
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByIntensity();
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector< DPickedPeak<2> > v(pl2.size());
	std::copy(pl2.begin(), pl2.end(), v.begin());
	TEST_EQUAL(v.size(), 3)
	ABORT_IF(v.size() != 3)
	TEST_REAL_EQUAL(v[2].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[0].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(void sortByNthPosition(UnsignedInt i) throw(Exception::NotImplemented))
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByNthPosition(0);
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector< DPickedPeak<2> > v(pl2.size());
	std::copy(pl2.begin(), pl2.end(), v.begin());
	TEST_EQUAL(v.size(), 3)
	ABORT_IF(v.size() != 3)
	TEST_REAL_EQUAL(v[1].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[0].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak3.getPosition()[1])


	pl2.sortByNthPosition(1);
	TEST_EQUAL(pl2.size(), 3)
	
	std::copy(pl2.begin(), pl2.end(), v.begin());
	TEST_EQUAL(v.size(), 3)
	ABORT_IF(v.size() != 3)
	TEST_REAL_EQUAL(v[2].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[0].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(Iterator begin())
	DPeakList<2, DPickedPeak<2> >::Iterator it = pl.begin();
	(*it).getIntensity()=1.4;
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(Iterator end())
	DPeakList<2, DPickedPeak<2> >::Iterator it = pl.end();
	--it;
	(*it).getIntensity()=4.1;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstIterator begin() const)
	DPeakList<2, DPickedPeak<2> >::ConstIterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstIterator end() const)
	DPeakList<2, DPickedPeak<2> >::ConstIterator it = pl.end();
	--it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rbegin())
	DPeakList<2, DPickedPeak<2> >::ReverseIterator it = pl.rbegin();
	(*it).getIntensity()=1.5;
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rend())
	DPeakList<2, DPickedPeak<2> >::ReverseIterator it = pl.rend();
	--it;
	(*it).getIntensity()=4.2;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstReverseIterator rbegin() const)
	DPeakList<2, DPickedPeak<2> >::ConstReverseIterator it;
	it = pl.rbegin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstReverseIterator rend() const)
	DPeakList<2, DPickedPeak<2> >::ConstReverseIterator it = pl.rend();
	--it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(DPeakList(size_type n))
	DPeakList<1> pl2(2);
	TEST_REAL_EQUAL(pl2.size(), 2)
	DPeakList<1>::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0)
RESULT

CHECK(DPeakList(size_type n, const PeakType& peak))
	 DPickedPeak<2> peak5;
	peak5.getPosition()[0] = 1.1;
	peak5.getIntensity() = 5.1;
	DPeakList<2, DPickedPeak<2> > pl2(3, peak5);
	TEST_REAL_EQUAL(pl2.size(), 3)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 5.1)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 5.1)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 5.1)
	++it;
RESULT

CHECK(const_reference front() const)
	 DPickedPeak<2> peak6(pl.front());
	TEST_REAL_EQUAL(peak6.getIntensity(), 4.2)
	TEST_REAL_EQUAL(peak6.getPosition()[0], 2.0)
	TEST_REAL_EQUAL(peak6.getPosition()[1], 3.0)
RESULT

CHECK(const_reference back() const)
	TEST_REAL_EQUAL(pl.back().getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl.back().getPosition()[0], 10.5)
	TEST_REAL_EQUAL(pl.back().getPosition()[1], 0.0)
RESULT

CHECK(reference front())
	pl.front().getIntensity()=4711.0;
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
RESULT

CHECK(reference back())
	pl.back().getIntensity()=4711.1;
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.end();
	--it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(void pop_back())
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_REAL_EQUAL(pl2.size(), 3)
	pl2.pop_back();
	TEST_REAL_EQUAL(pl2.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
RESULT

CHECK(void pop_front())
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_REAL_EQUAL(pl2.size(), 3)
	pl2.pop_front();
	TEST_REAL_EQUAL(pl2.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(void push_front(const PeakType& x))
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl2.size(), 3)
	
	pl2.push_front(peak1);
	TEST_EQUAL(pl2.size(), 4)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 1)
	pl2.push_front(peak2);
	TEST_EQUAL(pl2.size(), 5)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 0.5)
	pl2.push_front(peak3);
	TEST_EQUAL(pl2.size(), 6)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 0.01)
RESULT

CHECK(void swap(DPeakList& list))
	DPeakList<2, DPickedPeak<2> > pl2;
	
	 DPickedPeak<2> peak1;
	peak1.getPosition()[0] = 2.0;
	peak1.getPosition()[1] = 3.0;
	peak1.getIntensity() = 1.0;
	pl2.push_back(peak1);
	
	 DPickedPeak<2> peak2;
	peak2.getPosition()[0] = 0.0;
	peak2.getPosition()[1] = 2.5;
	peak2.getIntensity() = 2.5;
	pl2.push_back(peak2);


	TEST_REAL_EQUAL(pl2.front().getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2.back().getIntensity(), 2.5)	
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl.size(), 3)
	
	pl.swap(pl2);
	
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl.size(), 2)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl2.back().getIntensity(), 4711.1)
	TEST_REAL_EQUAL(pl.front().getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl.back().getIntensity(), 2.5)
	
	swap(pl,pl2);
	
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl.front().getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl.back().getIntensity(), 4711.1)
	TEST_REAL_EQUAL(pl2.front().getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2.back().getIntensity(), 2.5)
RESULT

CHECK(Iterator insert(Iterator pos, const PeakType& peak))
	 DPickedPeak<2> peak1;
	peak1.getIntensity() = 4712.0;
	TEST_REAL_EQUAL(pl.size(), 3)
	pl.insert(pl.end(),peak1);
	TEST_REAL_EQUAL(pl.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4712.0)
RESULT

CHECK(Iterator erase(Iterator pos))
	TEST_REAL_EQUAL(pl.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.end();
	--it;
	pl.erase(it);
	TEST_REAL_EQUAL(pl.size(), 3)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(void insert(Iterator pos, size_type n, const PeakType& peak))
	DPickedPeak<2> peak1;
	peak1.getIntensity() = 4711.0;
	pl.clear();
	pl.push_back(peak1);
	pl.push_back(peak2);
	peak1.getIntensity() = 4711.1;
	pl.push_back(peak1);
	
	TEST_REAL_EQUAL(pl.size(), 3)
	peak1.getIntensity() = 4714.0;
	pl.insert(pl.begin(),3,peak1);
	TEST_REAL_EQUAL(pl.size(), 6)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4714.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4714.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4714.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(Iterator erase(Iterator pos))
	TEST_REAL_EQUAL(pl.size(), 6)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	++it; ++it; ++it;
	pl.erase(pl.begin(),it);
	TEST_REAL_EQUAL(pl.size(), 3)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(template<class InputIterator> void insert(Iterator pos, InputIterator f, InputIterator l))
	TEST_REAL_EQUAL(pl.size(), 3)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	++it;
	pl.insert(pl.begin(),it,pl.end());
	TEST_REAL_EQUAL(pl.size(), 5)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(template<class InputIterator> DPeakList(InputIterator f, InputIterator l))
	DPeakList<2, DPickedPeak<2> > pl2(pl.begin(),pl.end());
	TEST_REAL_EQUAL(pl2.size(), 5)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(bool operator == (const DPeakList& rhs) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl == pl2 , true)
	pl2.front().getIntensity()=4.345;
	TEST_EQUAL(pl == pl2 , false)
RESULT

CHECK(bool operator !=(const DPeakList& list) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl != pl2 , false)
	pl2.front().getIntensity()=4.345;
	TEST_EQUAL(pl != pl2 , true)
RESULT

CHECK(bool operator < (const DPeakList& list) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl < pl2, false)
	pl2.push_back( DPickedPeak<2>());
	TEST_EQUAL(pl < pl2 , true)
RESULT

CHECK(bool operator > (const DPeakList& list) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl > pl2, false)
	pl2.erase(pl2.begin());
	TEST_EQUAL(pl > pl2 , true)
RESULT

CHECK(bool operator <= (const DPeakList& list) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl <= pl2, true)
	pl2.push_back( DPickedPeak<2>());
	TEST_EQUAL(pl <= pl2 , true)
	pl2.erase(pl2.begin());
	pl2.erase(pl2.begin());
	TEST_EQUAL(pl <= pl2 , false)
RESULT

CHECK(bool operator >= (const DPeakList& list) const)
	DPeakList<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl >= pl2, true)
	pl2.erase(pl2.begin());
	TEST_EQUAL(pl >= pl2 , true)
	pl2.insert(pl2.end(),2,pl2.front());
	TEST_EQUAL(pl >= pl2 , false)
RESULT

CHECK(void resize(size_type new_size, const PeakType& t=PeakType()))
	TEST_REAL_EQUAL(pl.size(), 5)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4711.1)
RESULT

CHECK(void clear())
	TEST_REAL_EQUAL(pl.size(), 2)
	pl.clear();
	TEST_REAL_EQUAL(pl.size(), 0)
RESULT

CHECK(void resize(size_type new_size, const PeakType& t=PeakType()))
	TEST_REAL_EQUAL(pl.size(), 0)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
RESULT

CHECK(void resize(size_type new_size, const PeakType& t=PeakType()))
	TEST_REAL_EQUAL(pl.size(), 2)
	 DPickedPeak<2> peak;
	peak.getIntensity()=4713.0;
	pl.resize(4,peak);
	TEST_EQUAL(pl.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void splice(iterator position, DPeakList& x))
	DPeakList<2, DPickedPeak<2> > dpl2;
	dpl2.push_back(peak1);
	dpl2.push_back(peak2);
	dpl2.push_back(peak3);
	TEST_EQUAL(dpl2.size(), 3)
	pl.splice(pl.begin(),dpl2);
	TEST_EQUAL(dpl2.size(), 0)
	TEST_EQUAL(pl.size(), 7)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void splice(iterator position, DPeakList& x, iterator i))
	DPeakList<2, DPickedPeak<2> > dpl2;
	dpl2.push_back(peak1);
	dpl2.push_back(peak2);
	dpl2.push_back(peak3);
	TEST_EQUAL(dpl2.size(), 3)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	++it; ++it; ++it; ++it;
	DPeakList<2, DPickedPeak<2> >::iterator it2 = dpl2.begin();
	++it2; 

	pl.splice(it,dpl2,it2);
	
	TEST_EQUAL(dpl2.size(), 2)
	it = dpl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	TEST_EQUAL(pl.size(), 8)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void splice(iterator position, DPeakList& x, iterator f, iterator l))
	DPeakList<2, DPickedPeak<2> > dpl2;
	dpl2.push_back(peak1);
	dpl2.push_back(peak2);
	dpl2.push_back(peak3);
	dpl2.push_back(peak3);
	TEST_EQUAL(dpl2.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it2 = dpl2.begin();
	++it2; 
	DPeakList<2, DPickedPeak<2> >::iterator it3 = dpl2.begin();
	++it3; ++it3; ++it3;
	
	pl.splice(pl.end(),dpl2,it2,it3);
	
	TEST_EQUAL(dpl2.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = dpl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	TEST_EQUAL(pl.size(), 10)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
RESULT

CHECK(void remove(const PeakType& p))
	pl.remove(peak1);
	TEST_EQUAL(pl.size(), 9)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	
	pl.remove(peak2);
	TEST_EQUAL(pl.size(), 6)
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
RESULT


CHECK(template<class Predicate> void remove_if(Predicate p))
	pl.remove_if(bind2nd(equal_to<  DPickedPeak<2> >() , peak3));
	TEST_EQUAL(pl.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void unique())
	pl.unique();
	TEST_EQUAL(pl.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(template<class BinaryPredicate> void unique(BinaryPredicate p))
	pl.insert(pl.begin(),1,peak1);
	pl.insert(pl.begin(),2,peak2);
	pl.insert(pl.begin(),3,peak3);
	TEST_EQUAL(pl.size(), 8)
	pl.unique(equal_to<  DPickedPeak<2> >());
	TEST_EQUAL(pl.size(), 5)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void sort())
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	it->getPosition()[0] = 1.0;
	it->getPosition()[1] = 3.0;
	++it;
	it->getPosition()[0] = 1.0;
	it->getPosition()[1] = 2.0;
	++it;
	it->getPosition()[0] = 0.0;
	it->getPosition()[1] = 0.0;
	++it;
	it->getPosition()[0] = 0.0;
	it->getPosition()[1] = 2.0;
	++it;
	it->getPosition()[0] = 0.0;
	it->getPosition()[1] = 1.0;
	pl.sort();
	it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
RESULT

CHECK(void merge(DPeakList& list))
	DPeakList<2, DPickedPeak<2> > dpl2;
	dpl2.push_back(peak2);
	dpl2.push_back(peak1);
	dpl2.push_back(peak3);
	pl.merge(dpl2);
	TEST_EQUAL(dpl2.size(), 0)
	TEST_EQUAL(pl.size(), 8)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
RESULT

CHECK(template<class BinaryPredicate> void sort(BinaryPredicate p))
	pl.sort( DPickedPeak<2>::IntensityLess());
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(template<class BinaryPredicate> void merge(DPeakList& list, BinaryPredicate p))
	DPeakList<2, DPickedPeak<2> > dpl2;
	dpl2.push_back(peak3);
	dpl2.push_back(peak2);
	dpl2.push_back(peak1);
	pl.merge(dpl2, DPickedPeak<2>::IntensityLess());
	TEST_EQUAL(dpl2.size(), 0)
	TEST_EQUAL(pl.size(), 11)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
RESULT

CHECK(void reverse())
	pl.reverse();
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4713.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.0)
RESULT

CHECK(template<class InputIterator> void assign(InputIterator f, InputIterator l))
	DPeakList<2, DPickedPeak<2> > dpa2;
	dpa2.push_back(peak1);
	dpa2.push_back(peak2);
	dpa2.push_back(peak3);
	TEST_EQUAL(pl.size(), 11)
	pl.assign(dpa2.begin(),dpa2.end());
	TEST_EQUAL(pl.size(), 3)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.5)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 0.01)
RESULT

CHECK(void assign(size_type n, const PeakType& x))
	pl.assign(5,peak3);
	TEST_EQUAL(pl.size(), 5)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl.begin();
	while (it!=pl.end())
	{
		TEST_REAL_EQUAL(it->getIntensity(), 0.01)
		++it;
	}
RESULT

CHECK(void sortByPosition())
DPeakList<2, DPickedPeak<2> > dpa2;
 DPickedPeak<2> p1(peak1);
p1.getIntensity()=1;
 DPickedPeak<2> p2(peak2);
p2.getIntensity()=2;
 DPickedPeak<2> p3(peak3);
p3.getIntensity()=3;
 DPickedPeak<2> p4(peak1);
p4.getPosition()[1]=4711;
p4.getIntensity()=4;
 DPickedPeak<2> p5(peak2);
p5.getPosition()[1]=4711;
p5.getIntensity()=5;
 DPickedPeak<2> p6(peak3);
p6.getPosition()[1]=4711;
p6.getIntensity()=6;
dpa2.push_back(p1);
dpa2.push_back(p2);
dpa2.push_back(p3);
dpa2.push_back(p4);
dpa2.push_back(p5);
dpa2.push_back(p6);
dpa2.sortByPosition();
DPeakList<2, DPickedPeak<2> >::iterator it = dpa2.begin();
TEST_REAL_EQUAL(it->getIntensity(), 2.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 5.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 1.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 4.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 3.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 6.0)
RESULT

CHECK(void sortByPosition())
DPeakList<2, DPickedPeak<2> > dpa2;
 DPickedPeak<2> p1(peak1);
p1.getIntensity()=1;
 DPickedPeak<2> p2(peak2);
p2.getIntensity()=2;
 DPickedPeak<2> p3(peak3);
p3.getIntensity()=3;
 DPickedPeak<2> p4(peak1);
p4.getPosition()[1]=4711;
p4.getIntensity()=4;
 DPickedPeak<2> p5(peak2);
p5.getPosition()[1]=4711;
p5.getIntensity()=5;
 DPickedPeak<2> p6(peak3);
p6.getPosition()[1]=4711;
p6.getIntensity()=6;
dpa2.push_back(p1);
dpa2.push_back(p2);
dpa2.push_back(p3);
dpa2.push_back(p4);
dpa2.push_back(p5);
dpa2.push_back(p6);
dpa2.sortByPosition();
DPeakList<2, DPickedPeak<2> >::iterator it=dpa2.begin();
TEST_REAL_EQUAL(it->getIntensity(), 2.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 5.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 1.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 4.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 3.0)
++it;
TEST_REAL_EQUAL(it->getIntensity(), 6.0)
RESULT

CHECK(template< typename ComparatorType > void sortByComparator())
	DPeakList<2, DPickedPeak<2> > pl2;
	pl2.push_back(peak1);
	pl2.push_back(peak2);
	pl2.push_back(peak3);
	
	pl2.sortByComparator<DPickedPeak<2>::PositionLess>();
	TEST_EQUAL(pl2.size(), 3)
	
	DPeakList<2, DPickedPeak<2> >::Iterator it=pl2.begin();

	TEST_REAL_EQUAL(it->getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak2.getPosition()[1])	
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak1.getPosition()[1])

	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak3.getPosition()[1])

	swap(*(pl2.begin()),*it);
	pl2.sortByComparator<DPickedPeak<2>::PositionLess>();
	
	it=pl2.begin();

	TEST_REAL_EQUAL(it->getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak2.getPosition()[1])	
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak1.getPosition()[1])

	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak3.getPosition()[1])
	
RESULT

CHECK(template< typename ComparatorType > void sortByComparator( ComparatorType const & comparator ))
	DPeakList<2, DPickedPeak<2> > pl2;
	pl2.push_back(peak1);
	pl2.push_back(peak2);
	pl2.push_back(peak3);
	pl2.sortByComparator<DPickedPeak<2>::NthPositionLess<1> >(DPickedPeak<2>::NthPositionLess<1>());
	TEST_EQUAL(pl2.size(), 3)
	
	DPeakList<2, DPickedPeak<2> >::Iterator it=pl2.begin();

	TEST_REAL_EQUAL(it->getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak3.getPosition()[1])
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak2.getPosition()[1])
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak1.getPosition()[1])

	swap(*(pl2.begin()),*it);
	pl2.sortByComparator<DPickedPeak<2>::NthPositionLess<0> >(DPickedPeak<2>::NthPositionLess<0>());
	
	it=pl2.begin();

	TEST_REAL_EQUAL(it->getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak2.getPosition()[1])
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak1.getPosition()[1])
	
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(it->getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(it->getPosition()[1], peak3.getPosition()[1])
RESULT

////////// Tests with an inhomogenous DPeakList ////////////
/////////////////////////////////////////////////////////////

DPeakList<1> dpa;
 DPickedPeak<1> p1,p3;
p1.getIntensity()=1;
p3.getIntensity()=3;
Labeled1DPeak p2,p4;
p2.getIntensity()=2;
p2.setLabel("L2");
p4.getIntensity()=4;
p4.setLabel("L4");

CHECK(push_back(const PeakType&) (inhomogenous list))
	TEST_EQUAL(dpa.size(), 0)
	dpa.push_back(p3);
	dpa.push_back(p4);
	dpa.push_front(p2);
	dpa.push_front(p1);
	TEST_EQUAL(dpa.size(), 4)
	DPeakList<1>::iterator it = dpa.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(back() (inhomogenous list))
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa.back()).getLabel(), "L4")
RESULT

CHECK(front() (inhomogenous list))
	TEST_EQUAL(typeid(dpa.front()) == typeid(Labeled1DPeak),false)
RESULT

CHECK(DPeakList(size_type n, const PrakType& p) (inhomogenous list))
	DPeakList<1> dpa2(4,dpa.back());
	DPeakList<1>::iterator it = dpa2.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(DPeakList( const PeakType& p) (inhomogenous list))
	DPeakList<1> dpa2(dpa);
	DPeakList<1>::iterator it = dpa2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(DPeakList(InputIterator f, InputIterator l) (inhomogenous list))
	DPeakList<1> dpa2(dpa.begin(),dpa.end());
	DPeakList<1>::iterator it = dpa2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(operator = (inhomogenous list))
	DPeakList<1> dpa2;
	dpa2 = dpa;
	DPeakList<1>::iterator it = dpa2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3)
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(swap(DPeakList&) (inhomogenous list))
	DPeakList<1> dpa2(2,dpa.back());
	dynamic_cast<Labeled1DPeak&>(dpa2.front()).setLabel("dpa2L1");
	dynamic_cast<Labeled1DPeak&>(dpa2.back()).setLabel("dpa2L2");
	dpa2.swap(dpa);
	DPeakList<1>::iterator it = dpa2.begin();
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	it = dpa.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "dpa2L1")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "dpa2L2")
	swap(dpa,dpa2);
	it = dpa.begin();
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	it = dpa2.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "dpa2L1")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "dpa2L2")
RESULT

CHECK(resize(size_type n, DPeakList& p) (inhomogenous list))
	TEST_EQUAL(dpa.size(), 4)
	dpa.resize(2);
	TEST_EQUAL(dpa.size(), 2)
	dpa.resize(4,p4);
	TEST_EQUAL(dpa.size(), 4)
	DPeakList<1>::iterator it = dpa.begin();
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(insert(Iterator pos, DPeakList& p) (inhomogenous list))
	dpa.insert(dpa.begin(),p4);
	DPeakList<1>::iterator it = dpa.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(insert(Iterator pos, size_type n, DPeakList& p) (inhomogenous list))
	DPeakList<1>::iterator it = dpa.begin();
	++it;
	dpa.insert(it,2,p2);
	it = dpa.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(insert(Iterator pos, InputIterator f, InputIterator l) (inhomogenous list))
	dpa.insert(dpa.end(),dpa.begin(),dpa.end());
	DPeakList<1>::iterator it = dpa.begin();
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(typeid(*it) == typeid(Labeled1DPeak),false)
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L2")
	++it;
	TEST_EQUAL(it->getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
	++it;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(*it).getLabel(), "L4")
RESULT

CHECK(operator == (inhomogenous list))
	DPeakList<1> dpa2(dpa);
	TEST_EQUAL( dpa == dpa2 ,true)
	dynamic_cast<Labeled1DPeak&>(dpa2.front()).getIntensity()=1.234234;
	TEST_EQUAL( dpa == dpa2 ,false)
	dpa2 = dpa;
	TEST_EQUAL( dpa == dpa2 ,true)
	dynamic_cast<Labeled1DPeak&>(dpa2.front()).setLabel("test");
	TEST_EQUAL( dpa == dpa2 ,false)
RESULT

CHECK(operator != (inhomogenous list))
	DPeakList<1> dpa2(dpa);
	TEST_EQUAL( dpa != dpa2 ,false)
	dpa2.front().getIntensity()=1.234234;
	TEST_EQUAL( dpa != dpa2 ,true)
	dpa2 = dpa;
	TEST_EQUAL( dpa != dpa2 ,false)
	dynamic_cast<Labeled1DPeak&>(dpa2.front()).setLabel("test");
	TEST_EQUAL( dpa != dpa2 ,true)
RESULT

 DPickedPeak<2> p1_2,p3_2;
p1_2.getIntensity()=1;
p1_2.getPosition()[0]=132;
p1_2.getPosition()[1]=12;
p3_2.getIntensity()=3;
p3_2.getPosition()[0]=9;
p3_2.getPosition()[1]=34;

Labeled2DPeak p2_2,p4_2;
p2_2.getIntensity()=2;
p2_2.getPosition()[0]=11;
p2_2.getPosition()[1]=3;
p2_2.label()="L2";
p4_2.getIntensity()=4;
p4_2.getPosition()[0]=1;
p4_2.getPosition()[1]=17;
p4_2.label()="L4";

Labeled2DPeak p2_2_mutant(p2_2);
p2_2_mutant.label()="L2_mutiert";

 DPickedPeak<2> p2_2_mutant2(p2_2);

DPeakList<2, DPickedPeak<2> > pl2;
pl2.push_back(p1_2);
pl2.push_back(p2_2);
pl2.push_back(p3_2);
pl2.push_back(p4_2);

CHECK(sort() / sorting by intensity/width/position (inhomogenous list))
	DPeakList<2, DPickedPeak<2> >::iterator it;
	DPeakList<2, DPickedPeak<2> > dpa3(pl2);
	it = dpa3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")

	dpa3.sortByNthPosition(0);
	it=dpa3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)

	dpa3.sortByNthPosition(1);
	it=dpa3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)

	dpa3.sortByIntensity();
	it = dpa3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")

	dpa3.sort();
	it=dpa3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
RESULT

DPeakList<2, DPickedPeak<2> > pl3=pl2;


CHECK(splice(interator pos, DPeakList list) (inhomogenous))
	TEST_EQUAL(pl2.size(), 4)
	TEST_EQUAL(pl3.size(), 4)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	pl2.splice(pl2.begin(),pl3);
	TEST_EQUAL(pl2.size(), 8)
	TEST_EQUAL(pl3.size(), 0)
	it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
RESULT

CHECK(splice(interator pos, DPeakList list, iterator i) (inhomogenous))
	pl3.splice(pl3.begin(),pl2,pl2.begin());
	pl3.splice(pl3.begin(),pl2,pl2.begin());
	TEST_EQUAL(pl2.size(), 6)
	TEST_EQUAL(pl3.size(), 2)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
RESULT

CHECK(splice(interator pos, DPeakList list, iterator f, iterator l) (inhomogenous))
	DPeakList<2, DPickedPeak<2> >::iterator it = pl2.begin();
	++it; 
	++it;
	pl3.splice(pl3.end(),pl2,pl2.begin(),it);
	TEST_EQUAL(pl2.size(), 4)
	TEST_EQUAL(pl3.size(), 4)
	it = pl2.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
RESULT

CHECK(remove(DPeak p) (inhomogenous))
	pl3.push_front(p2_2_mutant);
	pl3.push_front(p2_2_mutant2);
	TEST_EQUAL(pl3.size(), 6)
	pl3.remove(p2_2);
	TEST_EQUAL(pl3.size(), 5)
	pl3.remove(p2_2_mutant);
	TEST_EQUAL(pl3.size(), 4)
	pl3.remove(p2_2_mutant2);
	TEST_EQUAL(pl3.size(), 3)	
RESULT

CHECK(merge(DPeakList l) (inhomogenous))
	pl3=pl2;
	TEST_EQUAL(pl3.size(), 4)
	TEST_EQUAL(pl2.size(), 4)
	pl3.sort();
	pl2.sort();
	pl3.merge(pl2);
	TEST_EQUAL(pl3.size(), 8)
	TEST_EQUAL(pl2.size(), 0)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
RESULT

CHECK(unique() (inhomogenous))
	pl3.unique();
	TEST_EQUAL(pl3.size(), 4)

	DPeakList<2, DPickedPeak<2> >::iterator it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
RESULT

CHECK(merge(DPeakList l, BinaryPredicate comp) (inhomogenous))
	pl2=pl3;
	TEST_EQUAL(pl3.size(), 4)
	TEST_EQUAL(pl2.size(), 4)
	pl3.sortByIntensity();
	pl2.sortByIntensity();
	pl3.merge(pl2, DPickedPeak<2>::IntensityLess());
	TEST_EQUAL(pl3.size(), 8)
	TEST_EQUAL(pl2.size(), 0)
	DPeakList<2, DPickedPeak<2> >::iterator it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
RESULT

CHECK(unique(BinaryPredicate b) (inhomogenous))
	pl3.unique(SameType());
	TEST_EQUAL(pl3.size(), 4)

	DPeakList<2, DPickedPeak<2> >::iterator it = pl3.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 2.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L2")
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 3.0)
	TEST_EQUAL(typeid(*it) == typeid(Labeled2DPeak),false)
	++it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.0)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(*it).label(), "L4")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
