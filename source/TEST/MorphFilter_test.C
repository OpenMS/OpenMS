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
#include <OpenMS/KERNEL/DPeakArray.h>


///////////////////////////
#include <OpenMS/FILTERING/BASELINE/MorphFilter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MorphFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MorphFilter* ptr = 0;
CHECK(MorphFilter())
  ptr = new MorphFilter();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MorphFilter())
  delete ptr;
RESULT

CHECK(MorphFilter& operator=(const MorphFilter& m))
  Param p;
  p.setValue("struc_elem_length", 3);
  
  MorphFilter m;
  m.setParam(p);
  m.setStrucElemSize(4);
  
  MorphFilter m_copy;
  m_copy = m;
  
  TEST_EQUAL(p == m_copy.getParam(),true)
  TEST_EQUAL(m_copy.getStrucElemSize(),4)
RESULT

CHECK(MorphFilter(const MorphFilter& m))
  Param p;
  p.setValue("struc_elem_length", 3);
  
  MorphFilter m;
  m.setParam(p);
  m.setStrucElemSize(4);
  
  MorphFilter m_copy(m);
  
  TEST_EQUAL(p == m_copy.getParam(),true)
  TEST_REAL_EQUAL(m_copy.getStrucElemSize(),4)
RESULT

CHECK(MorphFilter(const Param& parameters))
  Param p;
  p.setValue("struc_elem_length", 3);
  MorphFilter m(p);
  
  TEST_EQUAL(p == m.getParam(), true)
RESULT

CHECK(const Param& getParam() const)
  Param p;
  p.setValue("struc_elem_length", 3);
  MorphFilter m(p);
  
  TEST_EQUAL(p == m.getParam(), true)
RESULT

CHECK(const float& getStrucElemSize() const)
  MorphFilter m;
  
  TEST_REAL_EQUAL(m.getStrucElemSize(),3)
RESULT

CHECK(float& getStrucElemSize())
  MorphFilter m;
  m.getStrucElemSize() = 3.5;
  
  TEST_REAL_EQUAL(m.getStrucElemSize(),3.5)
RESULT

CHECK((template< typename InputPeakIterator, typename OutputPeakContainer > void dilatation(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)))
  DPeakArray<1, DRawDataPoint<1> > raw(5);
  raw[0].getIntensity() = 0;
  raw[1].getIntensity() = 1;
  raw[2].getIntensity() = 1;
  raw[3].getIntensity() = 1;
  raw[4].getIntensity() = 0;

  raw[0].getPos() = 0;
  raw[1].getPos() = 1;
  raw[2].getPos() = 2;
  raw[3].getPos() = 3;
  raw[4].getPos() = 4;
  
  DPeakArray< 1,DRawDataPoint<1> > filtered;
  
  MorphFilter m;
  unsigned int struc_length = 3;
  
  m.dilatation(raw.begin(),raw.end(),filtered,struc_length);
  
  TEST_EQUAL(filtered[0].getIntensity(), 1)
  TEST_EQUAL(filtered[1].getIntensity(), 1)
  TEST_EQUAL(filtered[2].getIntensity(), 1)
  TEST_EQUAL(filtered[3].getIntensity(), 1)
  TEST_EQUAL(filtered[4].getIntensity(), 1)    
RESULT

CHECK((template< typename InputPeakIterator, typename OutputPeakContainer > void erosion(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)))
  DPeakArray<1, DRawDataPoint<1> > raw(5);
  raw[0].getIntensity() = 0;
  raw[1].getIntensity() = 1;
  raw[2].getIntensity() = 1;
  raw[3].getIntensity() = 1;
  raw[4].getIntensity() = 0;

  raw[0].getPos() = 0;
  raw[1].getPos() = 1;
  raw[2].getPos() = 2;
  raw[3].getPos() = 3;
  raw[4].getPos() = 4;
  
  DPeakArray< 1,DRawDataPoint<1> > filtered;
  
  MorphFilter m;
  unsigned int struc_length = 3;
  
  m.erosion(raw.begin(),raw.end(),filtered,struc_length);
  
  TEST_EQUAL(filtered[0].getIntensity(), 0)
  TEST_EQUAL(filtered[1].getIntensity(), 0)
  TEST_EQUAL(filtered[2].getIntensity(), 1)
  TEST_EQUAL(filtered[3].getIntensity(), 0)
  TEST_EQUAL(filtered[4].getIntensity(), 0)    
RESULT

CHECK(void setParam(const Param& param))
  Param p;
  p.setValue("struc_elem_length", 3);
  MorphFilter m;
  m.setParam(p);
  
  TEST_EQUAL(p == m.getParam(), true)
RESULT

CHECK(void setStrucElemSize(const float& struc_size))
  MorphFilter m;
  m.setStrucElemSize(3.5);
  
  TEST_REAL_EQUAL(m.getStrucElemSize(),3.5)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



