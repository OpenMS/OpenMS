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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/DExtractSignalRegions.h>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST("DExtractSignalRegions<D,Container>", "$Id: DExtractSignalRegions_test.C,v 1.3 2006/04/19 08:38:31 elange Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DExtractSignalRegions<1>* esr_ptr;
CHECK(DExtractSignalRegions())
  esr_ptr = new DExtractSignalRegions<1>;
  TEST_NOT_EQUAL(esr_ptr, 0)
RESULT

CHECK(~DExtractSignalRegions())
    delete esr_ptr;
RESULT

CHECK(DExtractSignalRegions(const Param& parameters))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  DExtractSignalRegions<1,DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > > esr(param);

  TEST_EQUAL(esr.getDaltonPerSplit(),4)
RESULT

CHECK(DExtractSignalRegions(const DExtractSignalRegions& e))
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  esr.getMZdim()=1;
  esr.getRTdim()=0;
  esr.getDaltonPerSplit()=10;

  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr_copy(esr);
  TEST_EQUAL(esr_copy.getMZdim(),1)
  TEST_EQUAL(esr_copy.getRTdim(),0)
  TEST_EQUAL(esr_copy.getDaltonPerSplit(),10)
RESULT

CHECK(DExtractSignalRegions& operator=(const DExtractSignalRegions& e))
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  esr.getMZdim()=1;
  esr.getRTdim()=0;
  esr.getDaltonPerSplit() = 10;

  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr_copy;
  esr_copy = esr;
  TEST_EQUAL(esr_copy.getMZdim(), 1)
  TEST_EQUAL(esr_copy.getRTdim(),0)
  TEST_EQUAL(esr_copy.getDaltonPerSplit(),10)
RESULT

CHECK(Param& getParam())
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);

  DExtractSignalRegions<1,DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > > esr(param);
  TEST_REAL_EQUAL((esr.getParam()) == param, true)
RESULT

CHECK(const Param& getParam() const)
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  const DExtractSignalRegions<1,DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > > esr(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK(const int getDaltonPerSplit() const)
  const DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
RESULT

CHECK(const int getMZdim() const)
  const DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getMZdim(), 1)
RESULT

CHECK(const int getRTdim() const)
  const DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getRTdim(), 0)
RESULT

CHECK(Param& getParam())
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  const DExtractSignalRegions<1,DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > > esr(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK(int getDaltonPerSplit())
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
  esr.getDaltonPerSplit() = 123;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 123)
  esr.getDaltonPerSplit() = 0;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 0)
RESULT

CHECK(int getMZdim())
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getMZdim(), 1)
  esr.getMZdim() = 0;
  TEST_REAL_EQUAL(esr.getMZdim(), 0)
  esr.getMZdim() = 1;
  TEST_REAL_EQUAL(esr.getMZdim(), 1)
RESULT

CHECK(int getRTdim())
   DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getRTdim(), 0)
  esr.getRTdim() = 1;
  TEST_REAL_EQUAL(esr.getRTdim(), 1)
  esr.getRTdim() = 0;
  TEST_REAL_EQUAL(esr.getRTdim(), 0)
RESULT

CHECK(void setParam(const Param& param))
  Param param;
  param.setValue("PeakPickingParameter:Split:DaltonPerSplit",4);
  DExtractSignalRegions<1,DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > > esr;
  esr.setParam(param);

  TEST_REAL_EQUAL(esr.getParam() == param, true)
RESULT

CHECK(void setDaltonPerSplit(const int dalton_per_split))
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 10)
  esr.setDaltonPerSplit(123);
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 123)
  esr.setDaltonPerSplit(0);
  TEST_REAL_EQUAL(esr.getDaltonPerSplit(), 0)
RESULT

CHECK(void setMZdim(const int mz_dim))
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getMZdim(), 1)
  esr.setMZdim(0);
  TEST_REAL_EQUAL(esr.getMZdim(), 0)
  esr.setMZdim(1);
  TEST_REAL_EQUAL(esr.getMZdim(), 1)
RESULT

CHECK(void setRTdim(const int rt_dim))
  DExtractSignalRegions<2, DPeakArrayNonPolymorphic<2> > esr;
  TEST_REAL_EQUAL(esr.getRTdim(), 0)
  esr.setRTdim(1);
  TEST_REAL_EQUAL(esr.getRTdim(), 1)
  esr.setRTdim(0);
  TEST_REAL_EQUAL(esr.getRTdim(), 0)
RESULT

CHECK((void splitScan(ConstIterator it_begin, ConstIterator it_end, double noise_level, IteratorVector &splitted_array)))
  typedef  DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > RawData;
  DExtractSignalRegions<1, RawData > esr;
  esr.setDaltonPerSplit(2);
  RawData raw;
  raw.resize(20);
  vector<RawData::const_iterator> split_vector;

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

  esr.splitScan(raw.begin(), raw.end(),2,split_vector);

  TEST_EQUAL(split_vector.size(),4);
  TEST_REAL_EQUAL(split_vector[0]->getPosition()[0],1.5);
  TEST_REAL_EQUAL(split_vector[1]->getPosition()[0],5.);
  TEST_REAL_EQUAL(split_vector[2]->getPosition()[0],4.5);
  TEST_REAL_EQUAL(split_vector[3]->getPosition()[0],8.);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
