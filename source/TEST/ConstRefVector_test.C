// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef std::vector< Peak1D > PeakArrayType;
typedef std::vector< Peak2D > PeakArray2DType;

START_TEST(ConstRefVector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConstRefVector<PeakArrayType>* ptr = 0;
START_SECTION((ConstRefVector()))
	ptr = new ConstRefVector<PeakArrayType>();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~ConstRefVector()))
	delete ptr;
END_SECTION

START_SECTION((ConstRefVector(const ConstRefVector& p)))
  ConstRefVector<PeakArrayType> pl;
  Peak1D peak1;
  Peak1D peak2;
  peak1.setIntensity(1.0f);
  pl.push_back(peak1);
  peak2.setIntensity(2.0f);
  pl.push_back(peak2);

  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl2.size(), 2)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 2.0)
END_SECTION

START_SECTION((ConstRefVector& operator=(const ConstRefVector &rhs)))
  ConstRefVector<PeakArrayType> pl;
  Peak1D peak1;
  Peak1D peak2;
  peak1.setIntensity(1.0f);
  pl.push_back(peak1);
  peak2.setIntensity(2.0f);
  pl.push_back(peak2);

  ConstRefVector<PeakArrayType> pl2;
  pl2 = pl;
  TEST_EQUAL(pl2.size(), 2)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 2.0)
END_SECTION

ConstRefVector<PeakArrayType> pl;

Peak1D peak1;
peak1.setPosition(2.0);
peak1.setIntensity(1.0f);

Peak1D peak2;
peak2.setPosition(0.0);
peak2.setIntensity(0.5f);

Peak1D peak3;
peak3.setPosition(10.5);
peak3.setIntensity(0.01f);

START_SECTION((size_type size() const))
  TEST_EQUAL(pl.size(), 0)

  pl.push_back(peak1);
  TEST_EQUAL(pl.size(), 1)
END_SECTION

START_SECTION((void push_back(const ValueType &x)))
  pl.push_back(peak2);
  TEST_EQUAL(pl.size(), 2)
END_SECTION

START_SECTION((size_type max_size() const))
  ConstRefVector<PeakArrayType>::size_type max = pl.max_size();
  pl.push_back(peak3);
  TEST_EQUAL(pl.max_size() == max, true)
END_SECTION

START_SECTION((bool empty() const))
  TEST_EQUAL(pl.empty(), false)
END_SECTION

START_SECTION([EXTRA] ConstIterator begin() const)
  const ConstRefVector<PeakArrayType>& c_pl(pl);
  TEST_EQUAL(c_pl.size(), 3)
  ABORT_IF(c_pl.size() != 3)
  TEST_REAL_SIMILAR(c_pl.begin()->getIntensity(), peak1.getIntensity())
  TEST_REAL_SIMILAR(c_pl.begin()->getPosition()[0], peak1.getPosition()[0])
END_SECTION

START_SECTION([EXTRA] ConstIterator end() const)
  const ConstRefVector<PeakArrayType>& c_pl(pl);
  TEST_EQUAL(c_pl.size(), 3)
  ABORT_IF(c_pl.size() != 3)
  bool result = (c_pl.begin() == c_pl.end());
  TEST_EQUAL(result, false)
  const ConstRefVector<PeakArrayType> empty;
  result = (empty.begin() == empty.end());
  TEST_EQUAL(result, true)
  std::vector<Peak1D> v(c_pl.size());
  std::copy(c_pl.begin(), c_pl.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_SIMILAR(v[0].getIntensity(), peak1.getIntensity())
  TEST_REAL_SIMILAR(v[0].getPosition()[0], peak1.getPosition()[0])

  TEST_REAL_SIMILAR(v[1].getIntensity(), peak2.getIntensity())
  TEST_REAL_SIMILAR(v[1].getPosition()[0], peak2.getPosition()[0])

  TEST_REAL_SIMILAR(v[2].getIntensity(), peak3.getIntensity())
  TEST_REAL_SIMILAR(v[2].getPosition()[0], peak3.getPosition()[0])
END_SECTION

START_SECTION((void sortByIntensity(bool reverse=false)))
  ConstRefVector<PeakArrayType> pl2(pl);
  pl2.sortByIntensity();
  TEST_EQUAL(pl2.size(), 3)

  std::vector<Peak1D> v(pl2.size());
  std::copy(pl2.begin(), pl2.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_SIMILAR(v[2].getIntensity(), peak1.getIntensity())
  TEST_REAL_SIMILAR(v[2].getPosition()[0], peak1.getPosition()[0])

  TEST_REAL_SIMILAR(v[1].getIntensity(), peak2.getIntensity())
  TEST_REAL_SIMILAR(v[1].getPosition()[0], peak2.getPosition()[0])

  TEST_REAL_SIMILAR(v[0].getIntensity(), peak3.getIntensity())
  TEST_REAL_SIMILAR(v[0].getPosition()[0], peak3.getPosition()[0])
END_SECTION


ConstRefVector<PeakArray2DType> pl2;

Peak2D peak4;
peak4.getPosition()[0] = 2.0;
peak4.getPosition()[1] = 3.0;
peak4.setIntensity(1.0f);
pl2.push_back(peak4);

Peak2D peak5;
peak5.getPosition()[0] = 0.0;
peak5.getPosition()[1] = 2.5;
peak5.setIntensity(0.5f);
pl2.push_back(peak5);

Peak2D peak6;
peak6.getPosition()[0] = 10.5;
peak6.getPosition()[1] = 0.0;
peak6.setIntensity(0.01f);
pl2.push_back(peak6);

START_SECTION((Iterator begin()))
  ConstRefVector<PeakArrayType>::Iterator it = pl.begin();
  TEST_REAL_SIMILAR(it->getIntensity(), 1.0)
  TEST_REAL_SIMILAR(it->getPosition()[0], 2.0)
END_SECTION

START_SECTION((Iterator end()))
  ConstRefVector<PeakArrayType>::Iterator it = pl.end()-1;
  TEST_REAL_SIMILAR(it->getIntensity(), 0.01)
  TEST_REAL_SIMILAR(it->getPosition()[0], 10.5)
END_SECTION

START_SECTION((ConstIterator begin() const))
  ConstRefVector<PeakArrayType>::ConstIterator it = pl.begin();
  TEST_REAL_SIMILAR(it->getIntensity(), 1.0)
  TEST_REAL_SIMILAR(it->getPosition()[0], 2.0)
END_SECTION

START_SECTION((ConstIterator end() const))
  ConstRefVector<PeakArrayType>::ConstIterator it = pl.end();
  --it;
  TEST_REAL_SIMILAR(it->getIntensity(), 0.01)
  TEST_REAL_SIMILAR(it->getPosition()[0], 10.5)
END_SECTION

START_SECTION((ReverseIterator rbegin()))
  ConstRefVector<PeakArrayType>::ReverseIterator it = pl.rbegin();
  TEST_REAL_SIMILAR(it->getIntensity(), 0.01)
  TEST_REAL_SIMILAR(it->getPosition()[0], 10.5)
END_SECTION

START_SECTION((ReverseIterator rend()))
  ConstRefVector<PeakArrayType>::ReverseIterator it = pl.rend()-1;
  TEST_REAL_SIMILAR(it->getIntensity(), 1.0)
  TEST_REAL_SIMILAR(it->getPosition()[0], 2.0)
END_SECTION

START_SECTION((ConstReverseIterator rbegin() const))
  ConstRefVector<PeakArrayType>::ConstReverseIterator it = pl.rbegin();
  TEST_REAL_SIMILAR(it->getIntensity(), 0.01)
  TEST_REAL_SIMILAR(it->getPosition()[0], 10.5)
END_SECTION

START_SECTION((ConstReverseIterator rend() const))
  ConstRefVector<PeakArrayType>::ConstReverseIterator it = pl.rend()-1;
  TEST_REAL_SIMILAR(it->getIntensity(), 1.0)
  TEST_REAL_SIMILAR(it->getPosition()[0], 2.0)
END_SECTION

START_SECTION((size_type capacity() const))
  TEST_EQUAL(pl.capacity(), 3)
  TEST_EQUAL(pl.size(), 3)
END_SECTION

Peak1D peak7;
peak7.getPosition()[0] = 1.1;
peak7.setIntensity(1.1f);

START_SECTION((void reserve(size_type n)))
  pl.reserve(4);
  TEST_EQUAL(pl.size(), 3)
  TEST_EQUAL(pl.capacity(), 4)

  pl.push_back(peak7);

  TEST_EQUAL(pl.size(), 4)
  TEST_EQUAL(pl.capacity(), 4)
END_SECTION

START_SECTION((const_reference operator [](size_type n) const))
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[2].getPosition()[0], 10.5)

  TEST_REAL_SIMILAR(pl[3].getIntensity(), 1.1)
  TEST_REAL_SIMILAR(pl[3].getPosition()[0], 1.1)
END_SECTION

START_SECTION((ConstRefVector(size_type n)))
  ConstRefVector<PeakArrayType> pl2(2);

  TEST_EQUAL(pl2.size(), 2)
END_SECTION

START_SECTION((ConstRefVector(size_type n, const ValueType &element)))
  Peak2D peak;
  peak.getPosition()[0] = 1.1;
  peak.setIntensity(5.1f);
  ConstRefVector<PeakArray2DType> pl2(3, peak);
  TEST_EQUAL(pl2.size(), 3)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 5.1)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 5.1)
  TEST_REAL_SIMILAR(pl2[2].getIntensity(), 5.1)
END_SECTION

START_SECTION((const_reference front() const))
  Peak1D peak;
  peak = pl.front();

  TEST_REAL_SIMILAR(peak.getIntensity(), 1.0)
  TEST_REAL_SIMILAR(peak.getPosition()[0], 2)
END_SECTION

START_SECTION((const_reference back() const))
 	Peak1D peak;
	peak = pl.back();

  TEST_REAL_SIMILAR(peak.getIntensity(), 1.1)
  TEST_REAL_SIMILAR(peak.getPosition()[0], 1.1)
END_SECTION

START_SECTION((void pop_back()))
  TEST_EQUAL(pl.size(), 4)
  pl.pop_back();
  TEST_EQUAL(pl.size(), 3)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
END_SECTION

Peak1D peak8;
peak8.getPosition()[0] = 2.0;
peak8.setIntensity(1.0f);

Peak1D peak9;
peak9.getPosition()[0] = 0.0;
peak9.setIntensity(2.5f);

START_SECTION((void swap(ConstRefVector &array)))
  ConstRefVector<PeakArrayType> pl2;

  pl2.push_back(peak8);
  pl2.push_back(peak9);

  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 2.5)
  TEST_EQUAL(pl2.size(), 2)
  TEST_EQUAL(pl.size(), 3)

  pl.swap(pl2);

  TEST_EQUAL(pl2.size(), 3)
  TEST_EQUAL(pl.size(), 2)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl2[2].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 2.5)

  swap(pl,pl2);

  TEST_EQUAL(pl.size(), 3)
  TEST_EQUAL(pl2.size(), 2)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 2.5)
END_SECTION

Peak1D peak10;
peak10.setIntensity(4712.0);
START_SECTION((Iterator insert(Iterator pos, const ValueType &element)))
  TEST_EQUAL(pl.size(), 3)
  pl.insert(pl.end(),peak10);

  TEST_EQUAL(pl.size(), 4)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[3].getIntensity(), 4712.0)
END_SECTION

START_SECTION((Iterator erase(Iterator pos)))
  TEST_EQUAL(pl.size(), 4)
  pl.erase(pl.end()-1);

  TEST_EQUAL(pl.size(), 3)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
END_SECTION

START_SECTION((void insert(Iterator pos, size_type n, const ValueType &element)))
  peak10.setIntensity(4714.0);
  TEST_EQUAL(pl.size(), 3)
  pl.insert(pl.begin(),3,peak10);

  TEST_EQUAL(pl.size(), 6)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 4714.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 4714.0)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 4714.0)
  TEST_REAL_SIMILAR(pl[3].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[4].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[5].getIntensity(), 0.01)
END_SECTION

START_SECTION((template <class InputIterator> void insert(Iterator pos, InputIterator f, InputIterator l)))
  pl.erase(pl.begin(),pl.begin()+3);
  TEST_EQUAL(pl.size(), 3)
  pl.insert(pl.begin(),pl.begin()+1,pl.end());

  TEST_EQUAL(pl.size(), 5)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[3].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[4].getIntensity(), 0.01)
END_SECTION

START_SECTION((template <class InputIterator> ConstRefVector(InputIterator f, InputIterator l)))
  ConstRefVector<PeakArrayType> pl2(pl.begin()+1,pl.end()-1);
  TEST_EQUAL(pl2.size(), 3)
  TEST_REAL_SIMILAR(pl2[0].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl2[1].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl2[2].getIntensity(), 0.5)
END_SECTION

START_SECTION((bool operator==(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl.size(), pl2.size())
  TEST_EQUAL(pl == pl2 , true)
END_SECTION

START_SECTION((bool operator!=(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl.size(), pl2.size())
  TEST_EQUAL(pl != pl2 , false)
END_SECTION

START_SECTION((bool operator<(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl < pl2, false)
  pl2.push_back(Peak1D());
  TEST_EQUAL(pl < pl2 , true)
END_SECTION

START_SECTION((bool operator>(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl > pl2, false)
  pl2.erase(pl2.end()-1);
  TEST_EQUAL(pl > pl2 , true)
END_SECTION

START_SECTION((bool operator<=(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl <= pl2, true)
  pl2.push_back(Peak1D());
  TEST_EQUAL(pl <= pl2 , true)
  pl2.erase(pl2.begin()+1,pl2.end()-2);
  TEST_EQUAL(pl <= pl2 , false)
END_SECTION

START_SECTION((bool operator>=(const ConstRefVector &array) const))
  ConstRefVector<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl >= pl2, true)
  pl2.erase(pl2.end()-1);
  TEST_EQUAL(pl >= pl2 , true)
  pl2.insert(pl2.end(),2,pl2.front());
  TEST_EQUAL(pl >= pl2 , false)
END_SECTION

START_SECTION((void clear()))
  pl.clear();

  TEST_EQUAL(pl.size(), 0)
END_SECTION

Peak1D peak11;
peak11.setIntensity(4713.0);
START_SECTION((void resize(size_type new_size)))
  pl.resize(4,peak11);

  TEST_EQUAL(pl.size(), 4)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 4713.0)
  TEST_REAL_SIMILAR(pl[3].getIntensity(), 4713.0)
END_SECTION

START_SECTION((void resize(size_type new_size, const ValueType &t)))
  ConstRefVector<PeakArrayType> pl;
  Peak1D peak;
  peak.getPosition()[0] = 0.0;
  peak.setIntensity(2.5f);
  pl.resize(2,peak);

  TEST_EQUAL(pl.size(), 2)
  TEST_EQUAL(pl[0].getIntensity() == peak.getIntensity(),true)
  TEST_EQUAL(pl[0].getPosition() == peak.getPosition(),true)
  TEST_EQUAL(pl[1].getIntensity() == peak.getIntensity(),true)
  TEST_EQUAL(pl[1].getPosition() == peak.getPosition(),true)
END_SECTION

START_SECTION((ConstRefVector(ContainerType &p)))
  PeakArrayType pa(5);
  ConstRefVector<PeakArrayType> pl(pa);

 	for (Size i=0; i<pa.size(); ++i)
 	{
 		TEST_EQUAL(pa[i]== pl[i],true)
 	}
END_SECTION

START_SECTION((template <class InputIterator> void assign(InputIterator f , InputIterator l)))
  ConstRefVector<PeakArrayType> dpa2;
  dpa2.push_back(peak1);
  dpa2.push_back(peak2);
  dpa2.push_back(peak3);
  TEST_EQUAL(pl.size(), 4)
  pl.assign(dpa2.begin(),dpa2.end());
  TEST_EQUAL(pl.size(), 3)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 1.0)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.5)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
END_SECTION

START_SECTION((void assign(size_type n, const ValueType &x)))
  pl.assign(5,peak3);
  TEST_EQUAL(pl.size(), 5)
  TEST_REAL_SIMILAR(pl[0].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[1].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[2].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[3].getIntensity(), 0.01)
  TEST_REAL_SIMILAR(pl[4].getIntensity(), 0.01)
END_SECTION

START_SECTION((Iterator erase(Iterator first,Iterator last)))
  TEST_EQUAL(pl.size(), 5)
  pl.erase(pl.begin(),pl.end());

  TEST_EQUAL(pl.size(), 0)
END_SECTION

START_SECTION((void sortByPosition()))
	ConstRefVector<PeakArray2DType> dpa2;
	Peak2D p1(peak4);
	p1.setIntensity(1.0f);
	Peak2D p2(peak5);
	p2.setIntensity(2.0f);
	Peak2D p3(peak6);
	p3.setIntensity(3.0f);
	Peak2D p4;
	p4.getPosition()[0]=4.3;
	p4.getPosition()[1]=4711;
	p4.setIntensity(4.0f);
	Peak2D p5;
	p5.getPosition()[1]=4711;
	p5.setIntensity(5.0f);
	Peak2D p6;
	p6.getPosition()[1]=4711;
	p6.setIntensity(6.0f);
	dpa2.push_back(p1);
	dpa2.push_back(p2);
	dpa2.push_back(p3);
	dpa2.push_back(p4);
	dpa2.push_back(p5);
	dpa2.push_back(p6);
	dpa2.sortByPosition();
	TEST_REAL_SIMILAR(dpa2[0].getIntensity(), 2.0)
	TEST_REAL_SIMILAR(dpa2[1].getIntensity(), 5.0)
	TEST_REAL_SIMILAR(dpa2[2].getIntensity(), 6.0)
	TEST_REAL_SIMILAR(dpa2[3].getIntensity(), 1.0)
	TEST_REAL_SIMILAR(dpa2[4].getIntensity(), 4.0)
	TEST_REAL_SIMILAR(dpa2[5].getIntensity(), 3.0)
END_SECTION

START_SECTION((template <typename ComparatorType> void sortByComparator(ComparatorType const &comparator=ComparatorType())))
  pl2.sortByComparator<Peak2D::PositionLess>();
  TEST_EQUAL(pl2.size(), 3)

  TEST_REAL_SIMILAR(pl2[1].getIntensity(), peak4.getIntensity())
  TEST_REAL_SIMILAR(pl2[1].getPosition()[0], peak4.getPosition()[0])
  TEST_REAL_SIMILAR(pl2[1].getPosition()[1], peak4.getPosition()[1])

  TEST_REAL_SIMILAR(pl2[0].getIntensity(), peak5.getIntensity())
  TEST_REAL_SIMILAR(pl2[0].getPosition()[0], peak5.getPosition()[0])
  TEST_REAL_SIMILAR(pl2[0].getPosition()[1], peak5.getPosition()[1])

  TEST_REAL_SIMILAR(pl2[2].getIntensity(), peak6.getIntensity())
  TEST_REAL_SIMILAR(pl2[2].getPosition()[0], peak6.getPosition()[0])
  TEST_REAL_SIMILAR(pl2[2].getPosition()[1], peak6.getPosition()[1])

  // ----------------

  ConstRefVector<PeakArray2DType> dpa2;
	Peak2D p1(peak4);
	p1.setIntensity(1.0f);
	Peak2D p2(peak5);
	p2.setIntensity(2.0f);
	Peak2D p3(peak6);
	p3.setIntensity(3.0f);
	Peak2D p4;
	p4.getPosition()[0]=4.3;
	p4.getPosition()[1]=4711;
	p4.setIntensity(4.0f);
	Peak2D p5;
	p5.getPosition()[1]=4711;
	p5.setIntensity(5.0f);
	Peak2D p6;
	p6.getPosition()[1]=4711;
	p6.setIntensity(6.0f);
	dpa2.push_back(p1);
	dpa2.push_back(p2);
	dpa2.push_back(p3);
	dpa2.push_back(p4);
	dpa2.push_back(p5);
	dpa2.push_back(p6);


	dpa2.sortByComparator<Peak2D::MZLess >(Peak2D::MZLess());
	TEST_REAL_SIMILAR(dpa2[0].getIntensity(), 3.0)
	TEST_REAL_SIMILAR(dpa2[1].getIntensity(), 2.0)
	TEST_REAL_SIMILAR(dpa2[2].getIntensity(), 1.0)
	TEST_REAL_SIMILAR(dpa2[3].getIntensity(), 4.0)
	TEST_REAL_SIMILAR(dpa2[4].getIntensity(), 5.0)
	TEST_REAL_SIMILAR(dpa2[5].getIntensity(), 6.0)
END_SECTION

START_SECTION([EXTRA] Container without special members for sorting)
	vector<Int> vec(5);
	ConstRefVector<vector<Int> > ref_vec(vec);
	TEST_EQUAL(ref_vec.size(),5)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

