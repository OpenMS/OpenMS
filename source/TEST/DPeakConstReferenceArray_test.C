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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/Peak2D.h>

///////////////////////////
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPeakArray< Peak1D > PeakArrayType;
typedef DPeakArray< Peak2D > PeakArray2DType;

START_TEST(DPeakConstReferenceArray, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DPeakConstReferenceArray<PeakArrayType>* ptr = 0;
CHECK((DPeakConstReferenceArray()))
	ptr = new DPeakConstReferenceArray<PeakArrayType>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~DPeakConstReferenceArray()))
	delete ptr;
RESULT

CHECK((DPeakConstReferenceArray(const DPeakConstReferenceArray& p)))
  DPeakConstReferenceArray<PeakArrayType> pl;
  Peak1D peak1;
  Peak1D peak2;
  peak1.setIntensity(1.0);
  pl.push_back(peak1);
  peak2.setIntensity(2.0);
  pl.push_back(peak2);
  
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl2.size(), 2)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.0)
RESULT

CHECK((DPeakConstReferenceArray& operator=(const DPeakConstReferenceArray &rhs)))
  DPeakConstReferenceArray<PeakArrayType> pl;
  Peak1D peak1;
  Peak1D peak2;
  peak1.setIntensity(1.0);
  pl.push_back(peak1);
  peak2.setIntensity(2.0);
  pl.push_back(peak2);
  
  DPeakConstReferenceArray<PeakArrayType> pl2;
  pl2 = pl;
  TEST_EQUAL(pl2.size(), 2)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.0)
RESULT

DPeakConstReferenceArray<PeakArrayType> pl;

Peak1D peak1;
peak1.setPosition(2.0);
peak1.setIntensity(1.0);

Peak1D peak2;
peak2.setPosition(0.0);
peak2.setIntensity(0.5);

Peak1D peak3;
peak3.setPosition(10.5);
peak3.setIntensity(0.01);

CHECK((size_type size() const))
  TEST_EQUAL(pl.size(), 0)
  
  pl.push_back(peak1);
  TEST_EQUAL(pl.size(), 1)
RESULT

CHECK((void push_back(const PeakType& x)))
  pl.push_back(peak2);
  TEST_EQUAL(pl.size(), 2)
RESULT

CHECK((size_type max_size() const))
  DPeakConstReferenceArray<PeakArrayType>::size_type max = pl.max_size();
  pl.push_back(peak3);
  TEST_EQUAL(pl.max_size() == max, true)
RESULT

CHECK((bool empty() const))
  TEST_EQUAL(pl.empty(), false)
RESULT

CHECK([EXTRA] ConstIterator begin() const)
  const DPeakConstReferenceArray<PeakArrayType>& c_pl(pl);
  TEST_EQUAL(c_pl.size(), 3)
  ABORT_IF(c_pl.size() != 3)
  TEST_REAL_EQUAL(c_pl.begin()->getIntensity(), peak1.getIntensity())
  TEST_REAL_EQUAL(c_pl.begin()->getPosition()[0], peak1.getPosition()[0])
RESULT

CHECK([EXTRA] ConstIterator end() const)
  const DPeakConstReferenceArray<PeakArrayType>& c_pl(pl);
  TEST_EQUAL(c_pl.size(), 3)
  ABORT_IF(c_pl.size() != 3)
  bool result = (c_pl.begin() == c_pl.end());
  TEST_EQUAL(result, false)
  const DPeakConstReferenceArray<PeakArrayType> empty;
  result = (empty.begin() == empty.end());
  TEST_EQUAL(result, true)
  std::vector<Peak1D> v(c_pl.size());
  std::copy(c_pl.begin(), c_pl.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_EQUAL(v[0].getIntensity(), peak1.getIntensity())
  TEST_REAL_EQUAL(v[0].getPosition()[0], peak1.getPosition()[0])

  TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
  TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])

  TEST_REAL_EQUAL(v[2].getIntensity(), peak3.getIntensity())
  TEST_REAL_EQUAL(v[2].getPosition()[0], peak3.getPosition()[0])
RESULT

CHECK((void sortByIntensity()))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  pl2.sortByIntensity();
  TEST_EQUAL(pl2.size(), 3)
  
  std::vector<Peak1D> v(pl2.size());
  std::copy(pl2.begin(), pl2.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_EQUAL(v[2].getIntensity(), peak1.getIntensity())
  TEST_REAL_EQUAL(v[2].getPosition()[0], peak1.getPosition()[0])

  TEST_REAL_EQUAL(v[1].getIntensity(), peak2.getIntensity())
  TEST_REAL_EQUAL(v[1].getPosition()[0], peak2.getPosition()[0])

  TEST_REAL_EQUAL(v[0].getIntensity(), peak3.getIntensity())
  TEST_REAL_EQUAL(v[0].getPosition()[0], peak3.getPosition()[0])
RESULT


DPeakConstReferenceArray<PeakArray2DType> pl2;

Peak2D peak4;
peak4.getPosition()[0] = 2.0;
peak4.getPosition()[1] = 3.0;
peak4.setIntensity(1.0);
pl2.push_back(peak4);

Peak2D peak5;
peak5.getPosition()[0] = 0.0;
peak5.getPosition()[1] = 2.5;
peak5.setIntensity(0.5);
pl2.push_back(peak5);

Peak2D peak6;
peak6.getPosition()[0] = 10.5;
peak6.getPosition()[1] = 0.0;
peak6.setIntensity(0.01);
pl2.push_back(peak6);

CHECK((void sortByNthPosition(UInt i) throw (Exception::NotImplemented)))
  pl2.sortByNthPosition(0);
  TEST_EQUAL(pl2.size(), 3)
  
  std::vector<Peak2D> v(pl2.size());
  std::copy(pl2.begin(), pl2.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_EQUAL(v[1].getIntensity(), peak4.getIntensity())
  TEST_REAL_EQUAL(v[1].getPosition()[0], peak4.getPosition()[0])
  TEST_REAL_EQUAL(v[1].getPosition()[1], peak4.getPosition()[1])

  TEST_REAL_EQUAL(v[0].getIntensity(), peak5.getIntensity())
  TEST_REAL_EQUAL(v[0].getPosition()[0], peak5.getPosition()[0])
  TEST_REAL_EQUAL(v[0].getPosition()[1], peak5.getPosition()[1])

  TEST_REAL_EQUAL(v[2].getIntensity(), peak6.getIntensity())
  TEST_REAL_EQUAL(v[2].getPosition()[0], peak6.getPosition()[0])
  TEST_REAL_EQUAL(v[2].getPosition()[1], peak6.getPosition()[1])

  pl2.sortByNthPosition(1);
  TEST_EQUAL(pl2.size(), 3)
  
  std::copy(pl2.begin(), pl2.end(), v.begin());
  TEST_EQUAL(v.size(), 3)
  ABORT_IF(v.size() != 3)
  TEST_REAL_EQUAL(v[2].getIntensity(), peak4.getIntensity())
  TEST_REAL_EQUAL(v[2].getPosition()[0], peak4.getPosition()[0])
  TEST_REAL_EQUAL(v[2].getPosition()[1], peak4.getPosition()[1])

  TEST_REAL_EQUAL(v[1].getIntensity(), peak5.getIntensity())
  TEST_REAL_EQUAL(v[1].getPosition()[0], peak5.getPosition()[0])
  TEST_REAL_EQUAL(v[1].getPosition()[1], peak5.getPosition()[1])

  TEST_REAL_EQUAL(v[0].getIntensity(), peak6.getIntensity())
  TEST_REAL_EQUAL(v[0].getPosition()[0], peak6.getPosition()[0])
  TEST_REAL_EQUAL(v[0].getPosition()[1], peak6.getPosition()[1])
RESULT

CHECK((template < typename ComparatorType > void sortByComparator ()))
  pl2.sortByComparator<Peak2D::PositionLess>();
  TEST_EQUAL(pl2.size(), 3)
  
  TEST_REAL_EQUAL(pl2[1].getIntensity(), peak4.getIntensity())
  TEST_REAL_EQUAL(pl2[1].getPosition()[0], peak4.getPosition()[0])
  TEST_REAL_EQUAL(pl2[1].getPosition()[1], peak4.getPosition()[1])

  TEST_REAL_EQUAL(pl2[0].getIntensity(), peak5.getIntensity())
  TEST_REAL_EQUAL(pl2[0].getPosition()[0], peak5.getPosition()[0])
  TEST_REAL_EQUAL(pl2[0].getPosition()[1], peak5.getPosition()[1])

  TEST_REAL_EQUAL(pl2[2].getIntensity(), peak6.getIntensity())
  TEST_REAL_EQUAL(pl2[2].getPosition()[0], peak6.getPosition()[0])
  TEST_REAL_EQUAL(pl2[2].getPosition()[1], peak6.getPosition()[1])
RESULT

CHECK((Iterator begin()))
  DPeakConstReferenceArray<PeakArrayType>::Iterator it = pl.begin();
  TEST_REAL_EQUAL(it->getIntensity(), 1.0)
  TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
RESULT

CHECK((Iterator end()))
  DPeakConstReferenceArray<PeakArrayType>::Iterator it = pl.end()-1;
  TEST_REAL_EQUAL(it->getIntensity(), 0.01)
  TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
RESULT

CHECK((ConstIterator begin() const))
  DPeakConstReferenceArray<PeakArrayType>::ConstIterator it = pl.begin();
  TEST_REAL_EQUAL(it->getIntensity(), 1.0)
  TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
RESULT

CHECK((ConstIterator end() const))
  DPeakConstReferenceArray<PeakArrayType>::ConstIterator it = pl.end();
  --it;
  TEST_REAL_EQUAL(it->getIntensity(), 0.01)
  TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
RESULT

CHECK((ReverseIterator rbegin()))
  DPeakConstReferenceArray<PeakArrayType>::ReverseIterator it = pl.rbegin();
  TEST_REAL_EQUAL(it->getIntensity(), 0.01)
  TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
RESULT

CHECK((ReverseIterator rend()))
  DPeakConstReferenceArray<PeakArrayType>::ReverseIterator it = pl.rend()-1;
  TEST_REAL_EQUAL(it->getIntensity(), 1.0)
  TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
RESULT

CHECK((ConstReverseIterator rbegin() const))
  DPeakConstReferenceArray<PeakArrayType>::ConstReverseIterator it = pl.rbegin();
  TEST_REAL_EQUAL(it->getIntensity(), 0.01)
  TEST_REAL_EQUAL(it->getPosition()[0], 10.5)
RESULT

CHECK((ConstReverseIterator rend() const))
  DPeakConstReferenceArray<PeakArrayType>::ConstReverseIterator it = pl.rend()-1;
  TEST_REAL_EQUAL(it->getIntensity(), 1.0)
  TEST_REAL_EQUAL(it->getPosition()[0], 2.0)
RESULT

CHECK((size_type capacity() const))
  TEST_EQUAL(pl.capacity(), 3)
  TEST_EQUAL(pl.size(), 3)
RESULT

Peak1D peak7;
peak7.getPosition()[0] = 1.1;
peak7.setIntensity(1.1);

CHECK((void reserve(size_type n)))
  pl.reserve(4);
  TEST_EQUAL(pl.size(), 3)
  TEST_EQUAL(pl.capacity(), 4)

  pl.push_back(peak7);
  
  TEST_EQUAL(pl.size(), 4)
  TEST_EQUAL(pl.capacity(), 4)
RESULT

CHECK((const_reference operator [](size_type n) const))
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[2].getPosition()[0], 10.5)
    
  TEST_REAL_EQUAL(pl[3].getIntensity(), 1.1)
  TEST_REAL_EQUAL(pl[3].getPosition()[0], 1.1)
RESULT

CHECK((DPeakConstReferenceArray(size_type n)))
  DPeakConstReferenceArray<PeakArrayType> pl2(2);
  
  TEST_REAL_EQUAL(pl2.size(), 2)
RESULT

CHECK((DPeakConstReferenceArray(size_type n, const PeakType &element)))
  Peak2D peak;
  peak.getPosition()[0] = 1.1;
  peak.setIntensity(5.1);
  DPeakConstReferenceArray<PeakArray2DType> pl2(3, peak);
  TEST_REAL_EQUAL(pl2.size(), 3)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 5.1)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 5.1)
  TEST_REAL_EQUAL(pl2[2].getIntensity(), 5.1)
RESULT

CHECK((const_reference front() const))
  Peak1D peak;
  peak = pl.front();
 
  TEST_REAL_EQUAL(peak.getIntensity(), 1.0)
  TEST_REAL_EQUAL(peak.getPosition()[0], 2) 
RESULT

CHECK((const_reference back() const))
 	Peak1D peak;
	peak = pl.back();
    
  TEST_REAL_EQUAL(peak.getIntensity(), 1.1)
  TEST_REAL_EQUAL(peak.getPosition()[0], 1.1)
RESULT

CHECK((void pop_back()))
  TEST_REAL_EQUAL(pl.size(), 4)
  pl.pop_back();
  TEST_REAL_EQUAL(pl.size(), 3)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
RESULT

Peak1D peak8;
peak8.getPosition()[0] = 2.0;
peak8.setIntensity(1.0);

Peak1D peak9;
peak9.getPosition()[0] = 0.0;
peak9.setIntensity(2.5);

CHECK((void swap(DPeakConstReferenceArray &array)))
  DPeakConstReferenceArray<PeakArrayType> pl2;
  
  pl2.push_back(peak8);
  pl2.push_back(peak9);

  TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.5) 
  TEST_REAL_EQUAL(pl2.size(), 2)
  TEST_REAL_EQUAL(pl.size(), 3)
  
  pl.swap(pl2);
  
  TEST_REAL_EQUAL(pl2.size(), 3)
  TEST_REAL_EQUAL(pl.size(), 2)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl2[2].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 2.5)
  
  swap(pl,pl2);
  
  TEST_REAL_EQUAL(pl.size(), 3)
  TEST_REAL_EQUAL(pl2.size(), 2)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 2.5)
RESULT

Peak1D peak10;
peak10.setIntensity(4712.0);
CHECK((Iterator insert(Iterator pos, const PeakType &element)))
  TEST_REAL_EQUAL(pl.size(), 3)
  pl.insert(pl.end(),peak10);
  
  TEST_REAL_EQUAL(pl.size(), 4)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[3].getIntensity(), 4712.0)
RESULT

CHECK((Iterator erase(Iterator pos)))
  TEST_REAL_EQUAL(pl.size(), 4)
  pl.erase(pl.end()-1);
   
  TEST_REAL_EQUAL(pl.size(), 3)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
RESULT

CHECK((void insert(Iterator pos, size_type n, const PeakType &element)))
  peak10.setIntensity(4714.0);
  TEST_REAL_EQUAL(pl.size(), 3)
  pl.insert(pl.begin(),3,peak10);
  
  TEST_REAL_EQUAL(pl.size(), 6)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 4714.0)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 4714.0)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 4714.0)
  TEST_REAL_EQUAL(pl[3].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[4].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[5].getIntensity(), 0.01)
RESULT

CHECK((template <class InputIterator> void insert(Iterator pos, InputIterator f, InputIterator l)))
  pl.erase(pl.begin(),pl.begin()+3);
  TEST_REAL_EQUAL(pl.size(), 3)
  pl.insert(pl.begin(),pl.begin()+1,pl.end());
   
  TEST_REAL_EQUAL(pl.size(), 5)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl[3].getIntensity(), 0.5)
  TEST_REAL_EQUAL(pl[4].getIntensity(), 0.01)
RESULT

CHECK((template <class InputIterator> DPeakConstReferenceArray(InputIterator f, InputIterator l)))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl.begin()+1,pl.end()-1);
  TEST_REAL_EQUAL(pl2.size(), 3)
  TEST_REAL_EQUAL(pl2[0].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl2[1].getIntensity(), 1.0)
  TEST_REAL_EQUAL(pl2[2].getIntensity(), 0.5)
RESULT

CHECK((bool operator==(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl.size(), pl2.size())
  TEST_EQUAL(pl == pl2 , true)
RESULT

CHECK((bool operator!=(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl.size(), pl2.size())
  TEST_EQUAL(pl != pl2 , false)
RESULT

CHECK((bool operator<(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl < pl2, false)
  pl2.push_back(Peak1D());
  TEST_EQUAL(pl < pl2 , true)
RESULT

CHECK((bool operator>(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl > pl2, false)
  pl2.erase(pl2.end()-1);
  TEST_EQUAL(pl > pl2 , true)
RESULT

CHECK((bool operator<=(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl <= pl2, true)
  pl2.push_back(Peak1D());
  TEST_EQUAL(pl <= pl2 , true)
  pl2.erase(pl2.begin()+1,pl2.end()-2);
  TEST_EQUAL(pl <= pl2 , false)
RESULT

CHECK((bool operator>=(const DPeakConstReferenceArray &array) const))
  DPeakConstReferenceArray<PeakArrayType> pl2(pl);
  TEST_EQUAL(pl >= pl2, true)
  pl2.erase(pl2.end()-1);
  TEST_EQUAL(pl >= pl2 , true)
  pl2.insert(pl2.end(),2,pl2.front());
  TEST_EQUAL(pl >= pl2 , false)
RESULT

CHECK((void clear()))
  pl.clear();
  
  TEST_REAL_EQUAL(pl.size(), 0)
RESULT

Peak1D peak11;
peak11.setIntensity(4713.0); 
CHECK((void resize(size_type new_size)))
  pl.resize(4,peak11);
    
  TEST_EQUAL(pl.size(), 4)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 4713.0)
  TEST_REAL_EQUAL(pl[3].getIntensity(), 4713.0)
RESULT

CHECK((void resize(size_type new_size, const PeakType& t)))
  DPeakConstReferenceArray<PeakArrayType> pl;
  Peak1D peak;
  peak.getPosition()[0] = 0.0;
  peak.setIntensity(2.5);
  pl.resize(2,peak);
    
  TEST_EQUAL(pl.size(), 2)
  TEST_EQUAL(pl[0].getIntensity() == peak.getIntensity(),true)
  TEST_EQUAL(pl[0].getPosition() == peak.getPosition(),true)
  TEST_EQUAL(pl[1].getIntensity() == peak.getIntensity(),true)
  TEST_EQUAL(pl[1].getPosition() == peak.getPosition(),true)
RESULT

CHECK((void setBaseContainerPointer(const BaseMapType *base_map_pointer)))
  DPeakConstReferenceArray<PeakArrayType> pl;
  PeakArrayType pa;
  pl.setBaseContainerPointer(&pa);
  
  TEST_EQUAL(pl.getBaseContainerPointer() == &pa, true)  
RESULT

CHECK((void setBaseContainer(const BaseMapType &base_map)))
  DPeakConstReferenceArray<PeakArrayType> pl;
  PeakArrayType pa;
  pl.setBaseContainer(pa);
  
  TEST_EQUAL(pl.getBaseContainer() == pa , true)  
RESULT

CHECK((const BaseMapType* getBaseContainerPointer() const))
  DPeakConstReferenceArray<PeakArrayType> pl;
    
  TEST_EQUAL(pl.getBaseContainerPointer() == 0 , true)  
RESULT

CHECK((const BaseMapType& getBaseContainer() const))
  DPeakConstReferenceArray<PeakArrayType> pl;
  PeakArrayType pa;
  pl.setBaseContainer(pa);
  
  TEST_EQUAL(pl.getBaseContainer() == pa , true)  

RESULT

CHECK((DPeakConstReferenceArray(BaseMapType &p)))
  PeakArrayType pa;
  DPeakConstReferenceArray<PeakArrayType> pl(pa);
 
  TEST_EQUAL(*(pl.getBaseContainerPointer()) == pa , true)  
RESULT

CHECK((template <class InputIterator> void assign(InputIterator f , InputIterator l)))
  DPeakConstReferenceArray<PeakArrayType> dpa2;
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

CHECK((void assign(size_type n, const PeakType &x)))
  pl.assign(5,peak3);
  TEST_EQUAL(pl.size(), 5)
  TEST_REAL_EQUAL(pl[0].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[1].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[2].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[3].getIntensity(), 0.01)
  TEST_REAL_EQUAL(pl[4].getIntensity(), 0.01)
RESULT

CHECK((Iterator erase(Iterator first,Iterator last)))
  TEST_REAL_EQUAL(pl.size(), 5)
  pl.erase(pl.begin(),pl.end());
   
  TEST_REAL_EQUAL(pl.size(), 0)
RESULT

CHECK((void sortByPosition()))
	DPeakConstReferenceArray<PeakArray2DType> dpa2;
	Peak2D p1(peak4);
	p1.setIntensity(1);
	Peak2D p2(peak5);
	p2.setIntensity(2);
	Peak2D p3(peak6);
	p3.setIntensity(3);
	Peak2D p4;
	p4.getPosition()[0]=4.3;
	p4.getPosition()[1]=4711;
	p4.setIntensity(4);
	Peak2D p5;
	p5.getPosition()[1]=4711;
	p5.setIntensity(5);
	Peak2D p6;
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
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 6.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 1.0)
	TEST_REAL_EQUAL(dpa2[4].getIntensity(), 4.0)
	TEST_REAL_EQUAL(dpa2[5].getIntensity(), 3.0)
RESULT

CHECK( template < typename ComparatorType > void sortByComparator ( ComparatorType const & comparator ))
	DPeakConstReferenceArray<PeakArray2DType> dpa2;
	Peak2D p1(peak4);
	p1.setIntensity(1);
	Peak2D p2(peak5);
	p2.setIntensity(2);
	Peak2D p3(peak6);
	p3.setIntensity(3);
	Peak2D p4;
	p4.getPosition()[0]=4.3;
	p4.getPosition()[1]=4711;
	p4.setIntensity(4);
	Peak2D p5;
	p5.getPosition()[1]=4711;
	p5.setIntensity(5);
	Peak2D p6;
	p6.getPosition()[1]=4711;
	p6.setIntensity(6);
	dpa2.push_back(p1);
	dpa2.push_back(p2);
	dpa2.push_back(p3);
	dpa2.push_back(p4);
	dpa2.push_back(p5);
	dpa2.push_back(p6);
	

	dpa2.sortByComparator<Peak2D::NthPositionLess<1> >(Peak2D::NthPositionLess<1>());
	TEST_REAL_EQUAL(dpa2[0].getIntensity(), 3.0)
	TEST_REAL_EQUAL(dpa2[1].getIntensity(), 2.0)
	TEST_REAL_EQUAL(dpa2[2].getIntensity(), 1.0)
	TEST_REAL_EQUAL(dpa2[3].getIntensity(), 4.0)
	TEST_REAL_EQUAL(dpa2[4].getIntensity(), 5.0)
	TEST_REAL_EQUAL(dpa2[5].getIntensity(), 6.0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

