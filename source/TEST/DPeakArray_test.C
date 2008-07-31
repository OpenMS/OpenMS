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
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DPeakArray<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DPeakArray<Peak2D >* pl_ptr = 0;
CHECK((DPeakArray()))
	pl_ptr = new DPeakArray<Peak2D >;
	TEST_NOT_EQUAL(pl_ptr, 0)
	TEST_EQUAL(pl_ptr->size(), 0)
RESULT

CHECK((~DPeakArray()))
	delete pl_ptr;
RESULT

CHECK((DPeakArray(const DPeakArray& p)))
	DPeakArray<DPeak<4> > pl;
	DPeak<4> peak;
	peak.setIntensity(1.0);
  pl.push_back(peak);
	peak.setIntensity(2.0);
  pl.push_back(peak);
  
  DPeakArray<DPeak<4> > pl2(pl);
	TEST_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.0)
RESULT

CHECK((template<class InputIterator> DPeakArray(InputIterator f, InputIterator l)))
	Peak1D peak;
	peak.getPosition() = 1.0;
	peak.setIntensity(1.01);
	
	DPeakArray<Peak1D> dpanp;
	dpanp.push_back(peak);
	peak.setIntensity(2.02);
	dpanp.push_back(peak);
	peak.setIntensity(3.03);
	dpanp.push_back(peak);
	peak.setIntensity(4.04);
	dpanp.push_back(peak);
	
	DPeakArray<Peak1D> dpanp2(dpanp.begin(),dpanp.end());
	
	TEST_EQUAL(dpanp.size(), dpanp2.size())
	TEST_EQUAL(dpanp[0] == dpanp2[0], true)
	TEST_EQUAL(dpanp[1] == dpanp2[1], true)
	TEST_EQUAL(dpanp[2] == dpanp2[2], true)
	TEST_EQUAL(dpanp[3] == dpanp2[3], true)
RESULT

Peak2D peak1;
peak1.getPosition()[0] = 2.0;
peak1.getPosition()[1] = 3.0;
peak1.setIntensity(1.0);

Peak2D peak2;
peak2.getPosition()[0] = 0.0;
peak2.getPosition()[1] = 2.5;
peak2.setIntensity(0.5);

Peak2D peak3;
peak3.getPosition()[0] = 10.5;
peak3.getPosition()[1] = 0.0;
peak3.setIntensity(0.01);

DPeakArray<Peak2D > pl;
pl.push_back(peak1);
pl.push_back(peak2);
pl.push_back(peak3);

CHECK((DPeakArray& operator = (const DPeakArray& rhs)))
	DPeakArray<Peak2D> copy_of_pl;
	TEST_EQUAL(copy_of_pl.size(), 0)
	copy_of_pl = pl;
	TEST_EQUAL(copy_of_pl.size(), 3)
	copy_of_pl = pl;
	TEST_EQUAL(copy_of_pl.size(), 3)
		
	std::vector<Peak2D > v(copy_of_pl.size());
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


CHECK((void sortByIntensity(bool reverse=false)))
	DPeakArray<Peak2D> pl2(pl);
	pl2.sortByIntensity();
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector<Peak2D > v(pl2.size());
	std::copy(pl2.begin(), pl2.end(), v.begin());
	TEST_REAL_EQUAL(v[2].getIntensity(), peak1.getIntensity())
	TEST_REAL_EQUAL(v[2].getPosition()[0], peak1.getPosition()[0])
	TEST_REAL_EQUAL(v[2].getPosition()[1], peak1.getPosition()[1])

	TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
	TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])
	TEST_REAL_EQUAL(v[1].getPosition()[1], peak2.getPosition()[1])

	TEST_REAL_EQUAL(v[0].getIntensity(), peak3.getIntensity())
	TEST_REAL_EQUAL(v[0].getPosition()[0], peak3.getPosition()[0])
	TEST_REAL_EQUAL(v[0].getPosition()[1], peak3.getPosition()[1])

	pl2.sortByIntensity(true);
	std::copy(pl2.begin(), pl2.end(), v.begin());
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

CHECK((void sortByNthPosition(UInt i) ))
	DPeakArray<Peak2D > pl2(pl);
	pl2.sortByNthPosition(0);
	TEST_EQUAL(pl2.size(), 3)
	
	std::vector<Peak2D > v(pl2.size());
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

CHECK((template< typename ComparatorType > void sortByComparator()))
	DPeakArray<Peak2D > pl2(pl);
	pl2.sortByComparator<Peak2D::PositionLess>();
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
	pl2.sortByComparator<Peak2D::PositionLess>();
	
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

CHECK((template< typename ComparatorType > void sortByComparator( ComparatorType const & comparator )))
	DPeakArray<Peak2D > pl2(pl);
	pl2.sortByComparator<Peak2D::NthPositionLess<1> >(Peak2D::NthPositionLess<1>());
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
	pl2.sortByComparator<Peak2D::NthPositionLess<0> >(Peak2D::NthPositionLess<0>());
	
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

CHECK((DPeakArray(typename std::vector< PeakType, AllocT >::size_type n)))
	DPeakArray<Peak1D> pl2(2);
	TEST_REAL_EQUAL(pl2.size(), 2)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 0)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 0)
RESULT

CHECK((DPeakArray(typename std::vector< PeakType, AllocT >::size_type n, const PeakType &peak)))
	Peak2D peak5;
	peak5.getPosition()[0] = 1.1;
	peak5.setIntensity(5.1);
	DPeakArray<Peak2D > pl2(3, peak5);
	TEST_REAL_EQUAL(pl2.size(), 3)
	TEST_REAL_EQUAL(pl2[0].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[1].getIntensity(), 5.1)
	TEST_REAL_EQUAL(pl2[2].getIntensity(), 5.1)
RESULT

CHECK((bool operator == (const DPeakArray& array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl == pl2 , true)
	pl2[0].setIntensity(4.345);
	TEST_EQUAL(pl == pl2 , false)
RESULT

CHECK((bool operator !=(const DPeakArray& array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl.size(), pl2.size())
	TEST_EQUAL(pl != pl2 , false)
	pl2[0].setIntensity(4.345);
	TEST_EQUAL(pl != pl2 , true)
RESULT

CHECK((template <typename AllocT2> bool operator<(const DPeakArray< PeakT, AllocT2 > &array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl < pl2, false)
	pl2.push_back(Peak2D());
	TEST_EQUAL(pl < pl2 , true)
RESULT

CHECK((template <typename AllocT2> bool operator>(const DPeakArray< PeakT, AllocT2 > &array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl > pl2, false)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl > pl2 , true)
RESULT

CHECK((bool operator <= (const DPeakArray& array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl <= pl2, true)
	pl2.push_back(Peak2D());
	TEST_EQUAL(pl <= pl2 , true)
	pl2.erase(pl2.begin()+1,pl2.end()-2);
	TEST_EQUAL(pl <= pl2 , false)
RESULT

CHECK((bool operator >= (const DPeakArray& array) const))
	DPeakArray<Peak2D > pl2(pl);
	TEST_EQUAL(pl >= pl2, true)
	pl2.erase(pl2.end()-1);
	TEST_EQUAL(pl >= pl2 , true)
	pl2.insert(pl2.end(),2,pl2.front());
	TEST_EQUAL(pl >= pl2 , false)
RESULT

CHECK((void sortByPosition()))
DPeakArray<Peak2D > dpa2;
Peak2D p1(peak1);
p1.setIntensity(1);
Peak2D p2(peak2);
p2.setIntensity(2);
Peak2D p3(peak3);
p3.setIntensity(3);
Peak2D p4(peak1);
p4.getPosition()[1]=4711;
p4.setIntensity(4);
Peak2D p5(peak2);
p5.getPosition()[1]=4711;
p5.setIntensity(5);
Peak2D p6(peak3);
p6.getPosition()[1]=4711;
p6.setIntensity(6);
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
