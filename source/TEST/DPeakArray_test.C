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
#include <OpenMS/KERNEL/DPickedPeak.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <string>

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

/////////////////////////////////////////////////////////////

START_TEST(DPeakArray<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PRECISION(0.0001)

DPeakArray<2,DPickedPeak<2> >* pl_ptr = 0;
CHECK(DPeakArray())
	pl_ptr = new DPeakArray<2 ,DPickedPeak<2> >;
	TEST_NOT_EQUAL(pl_ptr, 0)
	TEST_EQUAL(pl_ptr->size(), 0)
RESULT

CHECK(~DPeakArray())
	delete pl_ptr;
RESULT


CHECK(void push_back(const PeakType& x))
	DPeakArray<4 ,DPickedPeak<4> > pl;
	DPickedPeak<4> peak;
	peak.getIntensity() = 1.0;
  pl.push_back(peak);
	peak.getIntensity() = 2.0;
  pl.push_back(peak);
  
  DPeakArray<4 ,DPickedPeak<4> > pl2(pl);
	TEST_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.0)
RESULT

DPeakArray<2 ,DPickedPeak<2> > pl;

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

CHECK([EXTRA] bool empty() const)
	TEST_EQUAL(pl.empty(), false)
RESULT


CHECK([EXTRA] ConstIterator begin() const)
	const DPeakArray<2 ,DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	TEST_REAL_EQUAL(c_pl.begin()->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[1], peak1.getPosition()[1])
RESULT

CHECK([EXTRA] ConstIterator end() const)
	const DPeakArray<2 ,DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	bool result = (c_pl.begin() == c_pl.end());
	TEST_EQUAL(result, false)
	const DPeakArray<2 ,DPickedPeak<2> > empty;
	result = (empty.begin() == empty.end());
	TEST_EQUAL(result, true)
	std::vector<DPickedPeak<2> > v(c_pl.size());
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

CHECK(DPeakArray& operator = (const DPeakArray& rhs))
	DPeakArray<2 ,DPickedPeak<2> > copy_of_pl;
	TEST_EQUAL(copy_of_pl.size(), 0)
	copy_of_pl = pl;
	TEST_EQUAL(copy_of_pl.size(), 3)
	copy_of_pl = pl;
	TEST_EQUAL(copy_of_pl.size(), 3)
		
	std::vector<DPickedPeak<2> > v(copy_of_pl.size());
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

CHECK(DPeakArray(const DPeakArray& p))
	DPeakArray<2 ,DPickedPeak<2> > copy_of_pl(pl);
	TEST_EQUAL(copy_of_pl.size(), 3)
		
	std::vector<DPickedPeak<2> > v(copy_of_pl.size());
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
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByIntensity();
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector<DPickedPeak<2> > v(pl2.size());
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
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByNthPosition(0);
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector<DPickedPeak<2> > v(pl2.size());
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

CHECK(template< typename ComparatorType > void sortByComparator())
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByComparator<DPickedPeak<2>::PositionLess>();
	TEST_EQUAL(pl2.size(), 3)
	
	TEST_REAL_EQUAL(pl2[1].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(pl2[1].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(pl2[0].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(pl2[0].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(pl2[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak3.getPosition()[1])

	swap(pl2[0],pl2[2]);
	pl2.sortByComparator<DPickedPeak<2>::PositionLess>();
	
	TEST_REAL_EQUAL(pl2[1].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(pl2[1].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(pl2[0].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(pl2[0].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(pl2[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(template< typename ComparatorType > void sortByComparator( ComparatorType const & comparator ))
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	pl2.sortByComparator<DPickedPeak<2>::NthPositionLess<1> >(DPickedPeak<2>::NthPositionLess<1>());
	TEST_EQUAL(pl2.size(), 3)
	
	TEST_REAL_EQUAL(pl2[2].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(pl2[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(pl2[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(pl2[0].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(pl2[0].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak3.getPosition()[1])

	swap(pl2[0],pl2[2]);
	pl2.sortByComparator<DPickedPeak<2>::NthPositionLess<0> >(DPickedPeak<2>::NthPositionLess<0>());
	
	TEST_REAL_EQUAL(pl2[1].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(pl2[1].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(pl2[0].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(pl2[0].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(pl2[2].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak3.getPosition()[1])
RESULT

//
//DPickedPeak<2> peak1;
//peak1.getPosition()[0] = 2.0;
//peak1.getPosition()[1] = 3.0;
//peak1.getWidth()[0] = 1.0;
//peak1.getWidth()[1] = 0.5;
//peak1.getIntensity() = 1.0;
//
//DPickedPeak<2> peak2;
//peak2.getPosition()[0] = 0.0;
//peak2.getPosition()[1] = 2.5;
//peak2.getWidth()[0] = 0.2;
//peak2.getWidth()[1] = 0.2;
//peak2.getIntensity() = 0.5;
//
//DPickedPeak<2> peak3;
//peak3.getPosition()[0] = 10.5;
//peak3.getPosition()[1] = 0.0;
//peak3.getWidth()[0] = 0.5;
//peak3.getWidth()[1] = 0.01;
//peak3.getIntensity() = 0.01;


CHECK(Iterator begin())
	DPeakArray<2, DPickedPeak<2> >::Iterator it = pl.begin();
	(*it).getIntensity()=1.4;
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(Iterator end())
	DPeakArray<2, DPickedPeak<2> >::Iterator it = pl.end()-1;
	(*it).getIntensity()=4.1;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstIterator begin() const)
	DPeakArray<2, DPickedPeak<2> >::ConstIterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstIterator end() const)
	DPeakArray<2, DPickedPeak<2> >::ConstIterator it = pl.end();
	--it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rbegin())
	DPeakArray<2, DPickedPeak<2> >::ReverseIterator it = pl.rbegin();
	(*it).getIntensity()=1.5;
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rend())
	DPeakArray<2, DPickedPeak<2> >::ReverseIterator it = pl.rend()-1;
	(*it).getIntensity()=4.2;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstReverseIterator rbegin() const)
	DPeakArray<2, DPickedPeak<2> >::ConstReverseIterator it = pl.rbegin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstReverseIterator rend() const)
	DPeakArray<2, DPickedPeak<2> >::ConstReverseIterator it = pl.rend()-1;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(size_type capacity() const)
	TEST_EQUAL(pl.capacity(), 3)
	TEST_EQUAL(pl.size(), 3)
RESULT

CHECK(void reserve(size_type n))
	pl.reserve(4);
	TEST_EQUAL(pl.size(), 3)
	TEST_EQUAL(pl.capacity(), 4)

	DPickedPeak<2> peak4;
	peak4.getPosition()[0] = 1.1;
	peak4.getPosition()[1] = 1.1;
	peak4.getIntensity() = 1.1;
	pl.push_back(peak4);
	TEST_EQUAL(pl.size(), 4)
	TEST_EQUAL(pl.capacity(), 4)
RESULT

CHECK(const_reference operator [](size_type n) const)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[2].getPosition()[0], 10.5)
	TEST_REAL_EQUAL(pl[2].getPosition()[1], 0.0)
	
	TEST_REAL_EQUAL(pl[3].getIntensity(), 1.1)
	TEST_REAL_EQUAL(pl[3].getPosition()[0], 1.1)
	TEST_REAL_EQUAL(pl[3].getPosition()[1], 1.1)
RESULT

CHECK(reference operator [](size_type n))
	pl[3].getIntensity() = 1.2;
	pl[3].getPosition()[0] = 1.5;
	pl[3].getPosition()[1] = 1.6;

	TEST_REAL_EQUAL(pl[3].getIntensity(), 1.2)
	TEST_REAL_EQUAL(pl[3].getPosition()[0], 1.5)
	TEST_REAL_EQUAL(pl[3].getPosition()[1], 1.6)
RESULT

CHECK(DPeakArray(size_type n))
	DPeakArray<1> pl2(2);
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 0)
RESULT

CHECK(DPeakArray(size_type n, const PeakType& peak))
	DPickedPeak<2> peak5;
	peak5.getPosition()[0] = 1.1;
	peak5.getIntensity() = 5.1;
	DPeakArray<2, DPickedPeak<2> > pl2(3, peak5);
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 5.1)
RESULT

CHECK(const_reference front() const)
	DPickedPeak<2> peak6(pl.front());
	TEST_REAL_EQUAL(peak6.getIntensity(), 4.2)
	TEST_REAL_EQUAL(peak6.getPosition()[0], 2.0)
	TEST_REAL_EQUAL(peak6.getPosition()[1], 3.0)
RESULT

CHECK(const_reference back() const)
	TEST_REAL_EQUAL(pl.back().getIntensity(), 1.2)
	TEST_REAL_EQUAL(pl.back().getPosition()[0], 1.5)
	TEST_REAL_EQUAL(pl.back().getPosition()[1], 1.6)
RESULT

CHECK(reference front())
	pl.front().getIntensity()=4711.0;
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
RESULT

CHECK(reference back())
	pl.back().getIntensity()=4711.1;
	TEST_REAL_EQUAL(pl[3].getIntensity(), 4711.1)
RESULT

CHECK(void pop_back())
	TEST_REAL_EQUAL(pl.size(), 4)
	pl.pop_back();
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
RESULT

CHECK(void swap(DPeakArray& array))
	DPeakArray<2, DPickedPeak<2> > pl2;
	
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


	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.5)	
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl.size(), 3)
	
	pl.swap(pl2);

	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 2.5)
	
	pl.swap(pl2);
	
RESULT

CHECK(friend void swap(DPeakArray& a1, DPeakArray& a2))
	DPeakArray<2, DPickedPeak<2> > pkl, pkl2;
	
	DPickedPeak<2> peak1;
	peak1.getIntensity() = 1.0;
	
	DPickedPeak<2> peak2;
	peak2.getIntensity() = 2.5;
	
	pkl.push_back(peak1);
	pkl.push_back(peak2);
	pkl2.push_back(peak2);
	
	swap(pkl,pkl2);
	
	TEST_REAL_EQUAL(pkl.size(), 1)
	TEST_REAL_EQUAL(pkl2.size(), 2)
	TEST_REAL_EQUAL(pkl.front().getIntensity(), 2.5)
	TEST_REAL_EQUAL(pkl2.front().getIntensity(), 1.0)
	TEST_REAL_EQUAL(pkl2.back().getIntensity(), 2.5)
RESULT

CHECK(Iterator insert(Iterator pos, const PeakType& peak))
	DPickedPeak<2> peak1;
	peak1.getIntensity() = 4712.0;
	TEST_REAL_EQUAL(pl.size(), 3)
	pl.insert(pl.end(),peak1);
	TEST_REAL_EQUAL(pl.size(), 4)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 4712.0)
RESULT

CHECK(Iterator erase(Iterator pos))
	TEST_REAL_EQUAL(pl.size(), 4)
	pl.erase(pl.end()-1);
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
RESULT

CHECK(void insert(Iterator pos, size_type n, const PeakType& peak))
	DPickedPeak<2> peak1;
	peak1.getIntensity() = 4714.0;
	TEST_REAL_EQUAL(pl.size(), 3)
	pl.insert(pl.begin(),3,peak1);
	TEST_REAL_EQUAL(pl.size(), 6)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4714.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 4714.0)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 4714.0)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[4].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[5].getIntensity(), 1.5)
RESULT

CHECK(Iterator erase(Iterator first, Iterator last))
	TEST_REAL_EQUAL(pl.size(), 6)
	pl.erase(pl.begin(),pl.begin()+3);
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
RESULT

CHECK(template<class InputIterator> void insert(Iterator pos, InputIterator f, InputIterator l))
	TEST_REAL_EQUAL(pl.size(), 3)
	pl.insert(pl.begin(),pl.begin()+1,pl.end());
	TEST_REAL_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[4].getIntensity(), 1.5)
RESULT

CHECK(template<class InputIterator> DPeakArray(InputIterator f, InputIterator l))
	DPeakArray<2, DPickedPeak<2> > pl2(pl.begin()+1,pl.end()-1);
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 0.5)
RESULT

CHECK(bool operator == (const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl == pl2 , true)
	pl2[0].getIntensity()=4.345;
	TEST_EQUAL(pl == pl2 , false)
RESULT

CHECK(bool operator !=(const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl != pl2 , false)
	pl2[0].getIntensity()=4.345;
	TEST_EQUAL(pl != pl2 , true)
RESULT

CHECK(bool operator < (const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl < pl2, false)
	pl2.push_back(DPickedPeak<2>());
	TEST_EQUAL(pl < pl2 , true)
RESULT

CHECK(bool operator > (const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl > pl2, false)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl > pl2 , true)
RESULT

CHECK(bool operator <= (const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl <= pl2, true)
	pl2.push_back(DPickedPeak<2>());
	TEST_EQUAL(pl <= pl2 , true)
	pl2.erase(pl2.begin()+1,pl2.end()-2);
	TEST_EQUAL(pl <= pl2 , false)
RESULT

CHECK(bool operator >= (const DPeakArray& array) const)
	DPeakArray<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl >= pl2, true)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl >= pl2 , true)
	pl2.insert(pl2.end(),2,pl2.front());
	TEST_EQUAL(pl >= pl2 , false)
RESULT

CHECK(void resize(size_type new_size, const PeakType& t=PeakType()))
	TEST_REAL_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
RESULT

CHECK(void clear())
	TEST_REAL_EQUAL(pl.size(), 2)
	pl.clear();
	TEST_REAL_EQUAL(pl.size(), 0)
RESULT

CHECK([EXTRA] void resize(size_type new_size, const PeakType& t=PeakType()))
	TEST_REAL_EQUAL(pl.size(), 0)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)

	TEST_REAL_EQUAL(pl.size(), 2)
	DPickedPeak<2> peak;
	peak.getIntensity()=4713.0;
	pl.resize(4,peak);
	TEST_EQUAL(pl.size(), 4)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.0)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 4713.0)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 4713.0)
RESULT

CHECK(template<class InputIterator> void assign(InputIterator f, InputIterator l))
	DPeakArray<2, DPickedPeak<2> > dpa2;
	dpa2.push_back(peak1);
	dpa2.push_back(peak2);
	dpa2.push_back(peak3);
	TEST_EQUAL(pl.size(), 4)
	pl.assign(dpa2.begin(),dpa2.end());
	TEST_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
RESULT

CHECK(void assign(size_type n, const PeakType& x))
	pl.assign(5,peak3);
	TEST_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[4].getIntensity(), 0.01)
RESULT

CHECK(void sortByPosition())
DPeakArray<2, DPickedPeak<2> > dpa2;
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
TEST_REAL_EQUAL(dpa2[0].getIntensity(), 2.0)
TEST_REAL_EQUAL(dpa2[1].getIntensity(), 5.0)
TEST_REAL_EQUAL(dpa2[2].getIntensity(), 1.0)
TEST_REAL_EQUAL(dpa2[3].getIntensity(), 4.0)
TEST_REAL_EQUAL(dpa2[4].getIntensity(), 3.0)
TEST_REAL_EQUAL(dpa2[5].getIntensity(), 6.0)
RESULT

////////// Tests with an inhomogenous DPeakArray ////////////
/////////////////////////////////////////////////////////////

DPeakArray<1> dpa;
DPickedPeak<1> p1,p3;
p1.getIntensity()=1;
p3.getIntensity()=3;
Labeled1DPeak p2,p4;
p2.getIntensity()=2;
p2.setLabel("L2");
p4.getIntensity()=4;
p4.setLabel("L4");

CHECK([EXTRA] push_back(const PeakType&) / operator[](size_type n) (inhomogenous array))
	TEST_EQUAL(dpa.size(), 0)
	dpa.push_back(p1);
	dpa.push_back(p2);
	dpa.push_back(p3);
	dpa.push_back(p4);
	TEST_EQUAL(dpa.size(), 4)
	TEST_REAL_EQUAL(dpa[0].getIntensity(), 1)
	TEST_REAL_EQUAL(dpa[1].getIntensity(), 2)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "L2")
	TEST_REAL_EQUAL(dpa[2].getIntensity(), 3)
	TEST_REAL_EQUAL(dpa[3].getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] back() (inhomogenous array))
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa.back()).getLabel(), "L4")
RESULT

CHECK([EXTRA] DPeakArray(size_type n, const PrakType& p) (inhomogenous array))
	DPeakArray<1> dpa2(4,dpa.back());
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[0]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[2]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] DPeakArray( const PeakType& p) (inhomogenous array))
	DPeakArray<1> dpa2(dpa);
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] DPeakArray(InputIterator f, InputIterator l) (inhomogenous array))
	DPeakArray<1> dpa2(dpa.begin(),dpa.end());
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] operator = (inhomogenous array))
	DPeakArray<1> dpa2;
	dpa2 = dpa;
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] swap(DPeakArray&) (inhomogenous array))
	DPeakArray<1> dpa2(2,dpa.back());
	dynamic_cast<Labeled1DPeak&>(dpa2[0]).setLabel("dpa2L1");
	dynamic_cast<Labeled1DPeak&>(dpa2[1]).setLabel("dpa2L2");
	dpa2.swap(dpa);
	TEST_EQUAL(typeid(dpa2[0]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "L2")
	TEST_EQUAL(typeid(dpa2[2]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[3]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[0]).getLabel(), "dpa2L1")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "dpa2L2")
	swap(dpa,dpa2);
	TEST_EQUAL(typeid(dpa[0]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "L2")
	TEST_EQUAL(typeid(dpa[2]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[3]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[0]).getLabel(), "dpa2L1")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa2[1]).getLabel(), "dpa2L2")
RESULT

CHECK([EXTRA] resize(size_type n, DPeakArray& p) (inhomogenous array))
	TEST_EQUAL(dpa.size(), 4)
	dpa.resize(2);
	TEST_EQUAL(dpa.size(), 2)
	dpa.resize(4,p4);
	TEST_EQUAL(dpa.size(), 4)
	TEST_EQUAL(typeid(dpa[0]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "L2")
	TEST_EQUAL(dpa[2].getIntensity(), 4)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[2]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[3]).getLabel(), "L4")
RESULT

CHECK([EXTRA] insert(Iterator pos, DPeakArray& p) (inhomogenous array))
	dpa.insert(dpa.begin(),p4);
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[0]).getLabel(), "L4")
	TEST_EQUAL(typeid(dpa[1]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[2]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[3]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[4]).getLabel(), "L4")	
RESULT

CHECK([EXTRA] insert(Iterator pos, size_type n, DPeakArray& p) (inhomogenous array))
	dpa.insert(dpa.begin()+1,2,p2);
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[0]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[2]).getLabel(), "L2")
	TEST_EQUAL(typeid(dpa[3]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[4]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[5]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[6]).getLabel(), "L4")	
RESULT

CHECK([EXTRA] insert(Iterator pos, InputIterator f, InputIterator l) (inhomogenous array))
	dpa.insert(dpa.end(),dpa.begin(),dpa.end());
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[0]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[1]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[2]).getLabel(), "L2")
	TEST_EQUAL(typeid(dpa[3]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[4]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[5]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[6]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[7]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[8]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[9]).getLabel(), "L2")
	TEST_EQUAL(typeid(dpa[10]) == typeid(Labeled1DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[11]).getLabel(), "L2")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[12]).getLabel(), "L4")
	TEST_EQUAL(dynamic_cast<Labeled1DPeak&>(dpa[13]).getLabel(), "L4")	
RESULT

CHECK([EXTRA] bool operator == (const DPeakArray& array) const)
	DPeakArray<1> dpa2(dpa);
	TEST_EQUAL( dpa == dpa2 ,true)
	dynamic_cast<Labeled1DPeak&>(dpa2[0]).getIntensity()=1.234234;
	TEST_EQUAL( dpa == dpa2 ,false)
	dpa2 = dpa;
	TEST_EQUAL( dpa == dpa2 ,true)
	dynamic_cast<Labeled1DPeak&>(dpa2[0]).setLabel("test");
	TEST_EQUAL( dpa == dpa2 ,false)
RESULT

CHECK([EXTRA] bool operator !=(const DPeakArray& array) const)
	DPeakArray<1> dpa2(dpa);
	TEST_EQUAL( dpa != dpa2 ,false)
	dynamic_cast<Labeled1DPeak&>(dpa2[0]).getIntensity()=1.234234;
	TEST_EQUAL( dpa != dpa2 ,true)
	dpa2 = dpa;
	TEST_EQUAL( dpa != dpa2 ,false)
	dynamic_cast<Labeled1DPeak&>(dpa2[0]).setLabel("test");
	TEST_EQUAL( dpa != dpa2 ,true)
RESULT

CHECK([EXTRA] sorting by intensity/width/position (inhomogenous array))
	DPeakArray<2, DPickedPeak<2> > dpa2;
	DPickedPeak<2> p1,p3;
	p1.getIntensity()=1;
	p1.getPosition()[0]=132;
	p1.getPosition()[1]=12;
	p3.getIntensity()=3;
	p3.getPosition()[0]=9;
	p3.getPosition()[1]=34;

	Labeled2DPeak p2,p4;
	p2.getIntensity()=2;
	p2.getPosition()[0]=11;
	p2.getPosition()[1]=3;
	p2.label()="L2";
	p4.getIntensity()=4;
	p4.getPosition()[0]=1;
	p4.getPosition()[1]=17;
	p4.label()="L4";
	
	dpa2.push_back(p1);
	dpa2.push_back(p2);
	dpa2.push_back(p3);
	dpa2.push_back(p4);
	
	TEST_REAL_EQUAL(dpa2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(dpa2[1].getIntensity(), 2.0)
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 3.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 4.0)
	TEST_EQUAL(typeid(dpa2[0]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[1]).label(), "L2")
	TEST_EQUAL(typeid(dpa2[2]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[3]).label(), "L4")
	
	dpa2.sortByNthPosition(0);
	TEST_REAL_EQUAL(dpa2[0].getIntensity(), 4.0)
	TEST_REAL_EQUAL(dpa2[1].getIntensity(), 3.0)
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 2.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 1.0)	
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[0]).label(), "L4")
	TEST_EQUAL(typeid(dpa2[1]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[2]).label(), "L2")	
	TEST_EQUAL(typeid(dpa2[3]) == typeid(Labeled2DPeak),false)
	
	dpa2.sortByNthPosition(1);
	TEST_REAL_EQUAL(dpa2[0].getIntensity(), 2.0)
	TEST_REAL_EQUAL(dpa2[1].getIntensity(), 1.0)
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 4.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 3.0)	
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[0]).label(), "L2")
	TEST_EQUAL(typeid(dpa2[1]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[2]).label(), "L4")	
	TEST_EQUAL(typeid(dpa2[3]) == typeid(Labeled2DPeak),false)	

	dpa2.sortByIntensity();
	TEST_REAL_EQUAL(dpa2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(dpa2[1].getIntensity(), 2.0)
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 3.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 4.0)
	TEST_EQUAL(typeid(dpa2[0]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[1]).label(), "L2")
	TEST_EQUAL(typeid(dpa2[2]) == typeid(Labeled2DPeak),false)
	TEST_EQUAL(dynamic_cast<Labeled2DPeak&>(dpa2[3]).label(), "L4")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
