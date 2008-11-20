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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/Peak1D.h>


///////////////////////////
#include <OpenMS/FILTERING/BASELINE/MorphFilter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MorphFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MorphFilter* ptr = 0;
START_SECTION((MorphFilter()))
  ptr = new MorphFilter();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MorphFilter()))
  delete ptr;
END_SECTION

START_SECTION((template< typename InputPeakIterator, typename OutputPeakContainer > void dilatation(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)))
  std::vector<Peak1D > raw(5);
  raw[0].setIntensity(0);
  raw[1].setIntensity(1);
  raw[2].setIntensity(1);
  raw[3].setIntensity(1);
  raw[4].setIntensity(0);

  raw[0].setMZ(0);
  raw[1].setMZ(1);
  raw[2].setMZ(2);
  raw[3].setMZ(3);
  raw[4].setMZ(4);
  
  std::vector<Peak1D > filtered;
  
  MorphFilter m;
  UInt struc_length = 3;
  
  m.dilatation(raw.begin(),raw.end(),filtered,struc_length);
  
  TEST_EQUAL(filtered[0].getIntensity(), 1)
  TEST_EQUAL(filtered[1].getIntensity(), 1)
  TEST_EQUAL(filtered[2].getIntensity(), 1)
  TEST_EQUAL(filtered[3].getIntensity(), 1)
  TEST_EQUAL(filtered[4].getIntensity(), 1)    
END_SECTION

START_SECTION((template< typename InputPeakIterator, typename OutputPeakContainer > void erosion(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)))
  std::vector<Peak1D > raw(5);
  raw[0].setIntensity(0);
  raw[1].setIntensity(1);
  raw[2].setIntensity(1);
  raw[3].setIntensity(1);
  raw[4].setIntensity(0);

  raw[0].setMZ(0);
  raw[1].setMZ(1);
  raw[2].setMZ(2);
  raw[3].setMZ(3);
  raw[4].setMZ(4);
  
  std::vector<Peak1D > filtered;
  
  MorphFilter m;
  UInt struc_length = 3;
  
  m.erosion(raw.begin(),raw.end(),filtered,struc_length);
  
  TEST_EQUAL(filtered[0].getIntensity(), 0)
  TEST_EQUAL(filtered[1].getIntensity(), 0)
  TEST_EQUAL(filtered[2].getIntensity(), 1)
  TEST_EQUAL(filtered[3].getIntensity(), 0)
  TEST_EQUAL(filtered[4].getIntensity(), 0)    
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



