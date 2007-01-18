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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ExtractSignalRegions.h>
#include <OpenMS/KERNEL/DPeakArray.h>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST("ExtractSignalRegions<D,Container>", "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExtractSignalRegions* esr_ptr;
CHECK((ExtractSignalRegions()))
  esr_ptr = new ExtractSignalRegions;
  TEST_NOT_EQUAL(esr_ptr, 0)
RESULT

CHECK((~ExtractSignalRegions()))
    delete esr_ptr;
RESULT

CHECK((ExtractSignalRegions(const Param& parameters)))
  Param param;
  param.setValue("Split:DaltonPerSplit",4);
  ExtractSignalRegions esr(param);

  TEST_EQUAL(esr.getDaltonPerSplit(),4)
RESULT

CHECK((ExtractSignalRegions(const ExtractSignalRegions& e)))
  ExtractSignalRegions esr;
  esr.getDaltonPerSplit()=10;

  ExtractSignalRegions esr_copy(esr);
  TEST_EQUAL(esr_copy.getDaltonPerSplit(),10)
RESULT

CHECK((ExtractSignalRegions& operator=(const ExtractSignalRegions& e)))
  ExtractSignalRegions esr;
  esr.getDaltonPerSplit() = 10;

  ExtractSignalRegions esr_copy;
  esr_copy = esr;
  TEST_EQUAL(esr_copy.getDaltonPerSplit(),10)
RESULT

CHECK((Param& getParam()))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);

  ExtractSignalRegions esr(param);
  TEST_REAL_EQUAL((esr.getParam()) == param, true)
RESULT

CHECK((const Param& getParam() const))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  const ExtractSignalRegions esr(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK((const float& getDaltonPerSplit() const))
  const ExtractSignalRegions esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
RESULT

CHECK((Param& getParam()))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  const ExtractSignalRegions esr(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK((const float& getDaltonPerSplit() const))
  ExtractSignalRegions esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
  esr.getDaltonPerSplit() = 123;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 123)
  esr.getDaltonPerSplit() = 0;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 0)
RESULT

CHECK((void setParam(const Param& param)))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  ExtractSignalRegions esr;
  esr.setParam(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK((void setDaltonPerSplit(const float& dalton_per_split)))
  ExtractSignalRegions esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
  esr.setDaltonPerSplit(123);
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 123)
  esr.setDaltonPerSplit(0);
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 0)
RESULT

CHECK((template< typename InputPeakIterator > void splitScan(InputPeakIterator it_begin, InputPeakIterator it_end, double noise_level, std::vector<InputPeakIterator>& splitted_array)))
  typedef DPeakArray<1,DRawDataPoint<1> > RawData;
  ExtractSignalRegions esr;
  esr.setDaltonPerSplit(2);
  RawData raw;
  raw.resize(20);
  vector<RawData::iterator> split_vector;

  int i;
  for (i=0; i < 6; ++i)
  {
    DPosition<1> pos;
    pos=i*0.5;
    raw[i].setPosition(pos);
    raw[i].setIntensity(1);
  }

  for (; i < 14; ++i)
  {
    DPosition<1> pos;
    pos=i*0.5;
    raw[i].setPosition(pos);

    if ((i==6) || (i==8) || (i==13))
      {
        raw[i].setIntensity(50);
      }
    if (i==9)
      {
        raw[i].setIntensity(40);
      }
    if ((i==10) || (i==12))
      {
        raw[i].setIntensity(70);
      }
    if ((i==7) || (i==11))
      {
        raw[i].setIntensity(90);
      }
  }

  for (; i < 20; ++i)
  {
    DPosition<1> pos;
    pos=i*0.5;
    raw[i].setPosition(pos);
    raw[i].setIntensity(1);
  }

  esr.splitScan(raw.begin(),raw.end(),2,split_vector);

  TEST_EQUAL(split_vector.size(),2);
  TEST_REAL_EQUAL(split_vector[0]->getPos(),2.5);
  TEST_REAL_EQUAL((split_vector[1]-1)->getPos(),9.5);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
