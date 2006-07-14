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
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/DPickedPeak.h>
#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DPeakArrayNonPolymorphic<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DPeakArrayNonPolymorphic<2, DPickedPeak<2> >* pl_ptr = 0;
CHECK(DPeakArrayNonPolymorphic())
	pl_ptr = new DPeakArrayNonPolymorphic<2, DPickedPeak<2> >;
	TEST_NOT_EQUAL(pl_ptr, 0)
	TEST_EQUAL(pl_ptr->size(), 0)
RESULT

CHECK(~DPeakArrayNonPolymorphic())
	delete pl_ptr;
RESULT

CHECK(DPeakArrayNonPolymorphic(const DPeakArrayNonPolymorphic& p))
	DPeakArrayNonPolymorphic<4 ,DPickedPeak<4> > pl;
	DPickedPeak<4> peak;
	peak.getIntensity() = 1.0;
  pl.push_back(peak);
	peak.getIntensity() = 2.0;
  pl.push_back(peak);
  
  DPeakArrayNonPolymorphic<4 ,DPickedPeak<4> > pl2(pl);
	TEST_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.0)
RESULT

CHECK(DPeakArrayNonPolymorphic(typename std::vector<PeakType>::size_type n, const PeakType& peak))
	pl_ptr = new DPeakArrayNonPolymorphic<2, DPickedPeak<2> >;

	DPickedPeak<2> peak;
	peak.getPosition()[0] = 1.0;
	peak.getPosition()[1] = 2.0;
	peak.setIntensity(4.123);
	
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > dpanp(3,peak);
	
	TEST_EQUAL(dpanp.size(), 3)
	TEST_EQUAL(dpanp[0] == peak, true)
	TEST_EQUAL(dpanp[1] == peak, true)
	TEST_EQUAL(dpanp[2] == peak, true)
RESULT

CHECK(template<class InputIterator> DPeakArrayNonPolymorphic(InputIterator f, InputIterator l))
	DPickedPeak<1> peak;
	peak.getPosition()[0] = 1.0;
	peak.setIntensity(1.01);
	
	DPeakArrayNonPolymorphic<1> dpanp;
	dpanp.push_back(peak);
	peak.setIntensity(2.02);
	dpanp.push_back(peak);
	peak.setIntensity(3.03);
	dpanp.push_back(peak);
	peak.setIntensity(4.04);
	dpanp.push_back(peak);
	
	DPeakArrayNonPolymorphic<1> dpanp2(dpanp.begin(),dpanp.end());
	
	TEST_EQUAL(dpanp.size(), dpanp2.size())
	TEST_EQUAL(dpanp[0] == dpanp2[0], true)
	TEST_EQUAL(dpanp[1] == dpanp2[1], true)
	TEST_EQUAL(dpanp[2] == dpanp2[2], true)
	TEST_EQUAL(dpanp[3] == dpanp2[3], true)
RESULT

DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl;

CHECK(empty() const)
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

CHECK(size() const)
	TEST_EQUAL(pl.size(), 0)
	
	pl.push_back(peak1);
	TEST_EQUAL(pl.size(), 1)

	pl.push_back(peak2);
	TEST_EQUAL(pl.size(), 2)

	pl.push_back(peak3);
	TEST_EQUAL(pl.size(), 3)
RESULT

CHECK(empty() const)
	TEST_EQUAL(pl.empty(), false)
RESULT


CHECK([EXTRA] ConstIterator begin() const)
	const DPeakArrayNonPolymorphic<2, DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	TEST_REAL_EQUAL(c_pl.begin()->getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(c_pl.begin()->getPosition()[1], peak1.getPosition()[1])
RESULT

CHECK([EXTRA] ConstIterator end() const)
	const DPeakArrayNonPolymorphic<2, DPickedPeak<2> >& c_pl(pl);
	TEST_EQUAL(c_pl.size(), 3)
	ABORT_IF(c_pl.size() != 3)
	bool result = (c_pl.begin() == c_pl.end());
	TEST_EQUAL(result, false)
	const DPeakArrayNonPolymorphic<2, DPickedPeak<2> > empty;
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

CHECK(DPeakArrayNonPolymorphic& operator = (const DPeakArrayNonPolymorphic& rhs))
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > copy_of_pl;
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


CHECK(void sortByIntensity())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
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
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
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

	pl2.sortByNthPosition(0);
	pl2[0].getPosition()[0] = 2.0;
	pl2[1].getPosition()[0] = 2.0;
	pl2.sortByPosition();
	
	TEST_REAL_EQUAL(pl2[0].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(pl2[1].getPosition()[0], 2.0)
	TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak1.getPosition()[1])
	
	TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak3.getPosition()[1])
RESULT

CHECK(template< typename ComparatorType > void sortByComparator())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
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
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
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

CHECK(Iterator begin())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::Iterator it = pl.begin();
	(*it).getIntensity()=1.4;
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(Iterator end())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::Iterator it = pl.end()-1;
	(*it).getIntensity()=4.1;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstIterator begin())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ConstIterator it = pl.begin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.4)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstIterator end())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ConstIterator it = pl.end();
	--it;
	TEST_REAL_EQUAL(it->getIntensity(), 4.1)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rbegin())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ReverseIterator it = pl.rbegin();
	(*it).getIntensity()=1.5;
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ReverseIterator rend())
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ReverseIterator it = pl.rend()-1;
	(*it).getIntensity()=4.2;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(ConstReverseIterator rbegin() const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ConstReverseIterator it = pl.rbegin();
	TEST_REAL_EQUAL(it->getIntensity(), 1.5)
	TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
	TEST_REAL_EQUAL(it->getPosition()[1], 0.0)
RESULT

CHECK(ConstReverseIterator rend() const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> >::ConstReverseIterator it = pl.rend()-1;
	TEST_REAL_EQUAL(it->getIntensity(), 4.2)
	TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
	TEST_REAL_EQUAL(it->getPosition()[1], 3.0)
RESULT

CHECK(void reserve(size_type))
	pl.reserve(4);
	TEST_EQUAL(pl.size(), 3)

	DPickedPeak<2> peak4;
	peak4.getPosition()[0] = 1.1;
	peak4.getPosition()[1] = 1.1;
	peak4.getIntensity() = 1.1;
	pl.push_back(peak4);
	TEST_EQUAL(pl.size(), 4)
RESULT

CHECK(DPeakArrayNonPolymorphic& operator[] const)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[2].getPosition()[0], 10.5)
	TEST_REAL_EQUAL(pl[2].getPosition()[1], 0.0)
	
	TEST_REAL_EQUAL(pl[3].getIntensity(), 1.1)
	TEST_REAL_EQUAL(pl[3].getPosition()[0], 1.1)
	TEST_REAL_EQUAL(pl[3].getPosition()[1], 1.1)
RESULT

CHECK(DPeakArrayNonPolymorphic& operator[])
	pl[3].getIntensity() = 1.2;
	pl[3].getPosition()[0] = 1.5;
	pl[3].getPosition()[1] = 1.6;

	TEST_REAL_EQUAL(pl[3].getIntensity(), 1.2)
	TEST_REAL_EQUAL(pl[3].getPosition()[0], 1.5)
	TEST_REAL_EQUAL(pl[3].getPosition()[1], 1.6)
RESULT

CHECK(DPeakArrayNonPolymorphic(typename std::vector<PeakType>::size_type n))
	DPeakArrayNonPolymorphic<1> pl2(2);
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 0)
RESULT

CHECK(DPeakArrayNonPolymorphic(typename std::vector<PeakType>::size_type n, const PeakType& peak))
	DPickedPeak<2> peak5;
	peak5.getPosition()[0] = 1.1;
	peak5.getIntensity() = 5.1;
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(3, peak5);
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 5.1)
RESULT

CHECK(reference front() const)
	DPickedPeak<2> peak6(pl.front());
	TEST_REAL_EQUAL(peak6.getIntensity(), 4.2)
	TEST_REAL_EQUAL(peak6.getPosition()[0], 2.0)
	TEST_REAL_EQUAL(peak6.getPosition()[1], 3.0)
RESULT

CHECK(reference back() const)
	TEST_REAL_EQUAL(pl.back().getIntensity(), 1.2)
	TEST_REAL_EQUAL(pl.back().getPosition()[0], 1.5)
	TEST_REAL_EQUAL(pl.back().getPosition()[1], 1.6)
RESULT

CHECK(reference front())
	pl.front().getIntensity()=4711.0;
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
RESULT

CHECK(reference back())
	PRECISION(0.01)
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

CHECK(void swap(DPeakArrayNonPolymorphic))
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2;
	
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
	
	swap(pl,pl2);
	
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.5)
RESULT

CHECK(iterator insert(iterator pos, const DPickedPeak<D>&))
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

CHECK(iterator erase(iterator pos))
	TEST_REAL_EQUAL(pl.size(), 4)
	pl.erase(pl.end()-1);
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
RESULT

CHECK(iterator insert(iterator pos, size_type n, const DPickedPeak<D>&))
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

CHECK(iterator erase(iterator pos))
	TEST_REAL_EQUAL(pl.size(), 6)
	pl.erase(pl.begin(),pl.begin()+3);
	TEST_REAL_EQUAL(pl.size(), 3)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 1.5)
RESULT

CHECK(iterator insert(iterator pos, InputIterator f, InputIterator l))
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > tmp(pl);
	TEST_REAL_EQUAL(pl.size(), 3)
	pl.insert(pl.begin(),tmp.begin()+1,tmp.end());
	TEST_REAL_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[4].getIntensity(), 1.5)
RESULT

CHECK(template<class InputIterator> DPeakArrayNonPolymorphic(InputIterator f, InputIterator l))
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl.begin()+1,pl.end()-1);
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.5)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 4711.0)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 0.5)
RESULT

CHECK(bool operator == (const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl == pl2 , true)
	pl2[0].getIntensity()=4.345;
	TEST_EQUAL(pl == pl2 , false)
RESULT

CHECK(bool operator !=(const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl != pl2 , false)
	pl2[0].getIntensity()=4.345;
	TEST_EQUAL(pl != pl2 , true)
RESULT

CHECK(bool operator < (const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl < pl2, false)
	pl2.push_back(DPickedPeak<2>());
	TEST_EQUAL(pl < pl2 , true)
RESULT

CHECK(bool operator > (const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl > pl2, false)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl > pl2 , true)
RESULT

CHECK(bool operator <= (const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl <= pl2, true)
	pl2.push_back(DPickedPeak<2>());
	TEST_EQUAL(pl <= pl2 , true)
	pl2.erase(pl2.begin()+1,pl2.end()-2);
	TEST_EQUAL(pl <= pl2 , false)
RESULT

CHECK(bool operator >= (const DPeakArrayNonPolymorphic& array) const)
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > pl2(pl);
	TEST_EQUAL(pl >= pl2, true)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl >= pl2 , true)
	pl2.insert(pl2.end(),2,pl2.front());
	TEST_EQUAL(pl >= pl2 , false)
RESULT

CHECK(resize() (shrink))
	TEST_REAL_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 1.5)
RESULT

CHECK(clear())
	TEST_REAL_EQUAL(pl.size(), 2)
	pl.clear();
	TEST_REAL_EQUAL(pl.size(), 0)
RESULT

CHECK(resize() (expand))
	TEST_REAL_EQUAL(pl.size(), 0)
	pl.resize(2);
	TEST_REAL_EQUAL(pl.size(), 2)
RESULT

CHECK(resize() (expand))
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

CHECK(template <class InputIterator> void assign(InputIterator f , InputIterator l))
	DPeakArrayNonPolymorphic<2, DPickedPeak<2> > dpa2;
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

CHECK(void assign(size_type n , const PeakType& x))
	pl.assign(5,peak3);
	TEST_EQUAL(pl.size(), 5)
	TEST_REAL_EQUAL(pl[0].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[1].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[3].getIntensity(), 0.01)
	TEST_REAL_EQUAL(pl[4].getIntensity(), 0.01)
RESULT

CHECK(void sortByPosition())
DPeakArrayNonPolymorphic<2, DPickedPeak<2> > dpa2;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
