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
// $Id: DSignalToNoiseEstimatorWindowing_test.C,v 1.1 2006/04/24 13:44:14 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorWindowing.h>

///////////////////////////

START_TEST(DSignalToNoiseEstimator, "$Id: DSignalToNoiseEstimatorWindowing_test.C,v 1.1 2006/04/24 13:44:14 elange Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DSignalToNoiseEstimatorWindowing<1>* ptr;
CHECK(DSignalToNoiseEstimator())
  ptr = new DSignalToNoiseEstimatorWindowing<1>;
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DSignalToNoiseEstimatorWindowing())
  delete ptr;
RESULT

CHECK(DSignalToNoiseEstimatorWindowing(const Param& parameters))
  Param p;
  p.setValue("SignalToNoiseEstimationParameter:Bucket",6);
  p.setValue("SignalToNoiseEstimationParameter:Window",20);

  DSignalToNoiseEstimatorWindowing<1> sne(p);
  TEST_REAL_EQUAL(sne.getBucketSize(),6);
  TEST_REAL_EQUAL(sne.getWindowSize(),20);
RESULT

CHECK(DSignalToNoiseEstimatorWindowing(const DSignalToNoiseEstimatorWindowing& ne))
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.setBucketSize(6);
  sne.setWindowSize(20);

  DSignalToNoiseEstimatorWindowing<1> sne2(sne);
  TEST_EQUAL(sne2.getBucketSize(),6);
  TEST_REAL_EQUAL(sne2.getWindowSize(),20);
RESULT

CHECK((void init(ConstIterator it_begin, ConstIterator it_end)))
  DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > raw_data(2);

  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.init(raw_data.begin(),raw_data.end());

  TEST_EQUAL((sne.getFirstDataPoint() == raw_data.begin()),true);
  TEST_EQUAL((sne.getLastDataPoint() == raw_data.end()),true);
RESULT

CHECK(DSignalToNoiseEstimatorWindowing& operator=(const DSignalToNoiseEstimatorWindowing& ne))
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.setBucketSize(6);
  sne.setWindowSize(20);

  DSignalToNoiseEstimatorWindowing<1> sne2;
  sne2 = sne;
  TEST_EQUAL(sne2.getBucketSize(),6);
  TEST_REAL_EQUAL(sne2.getWindowSize(),20);
RESULT

CHECK(const int getBucketSize() const)
  const DSignalToNoiseEstimatorWindowing<1> sne;

  TEST_EQUAL(sne.getBucketSize(),10);
RESULT

CHECK(const int getWindowSize() const)
  const DSignalToNoiseEstimatorWindowing<1> sne;

  TEST_EQUAL(sne.getWindowSize(),700);
RESULT

CHECK(double getSignalToNoise(ConstIterator data_point) throw(Exception::OutOfRange))
  DPeakArrayNonPolymorphic<1,DRawDataPoint<1> > raw_data;
  int i;
  for (i=0; i < 6; ++i)
  {
    DRawDataPoint<1> p;
    DPosition<1> pos = i;
    if ((i == 2) || (i == 4))
    {
      p.setIntensity(100);
    }
    else
    {
      p.setIntensity(0);
    }
    p.setPosition(pos);
    raw_data.push_back(p);
  }

  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.setBucketSize(3);
  sne.setWindowSize(10);
  sne.init(raw_data.begin(),raw_data.end());

  for (i=0; i < (int)raw_data.size(); ++i)
  {
    DPeakArrayNonPolymorphic<1,DRawDataPoint<1> >::const_iterator first = raw_data.begin() + i;

    if ((i == 2) || (i == 4))
    {
      TEST_REAL_EQUAL(sne.getSignalToNoise(first),1);
    }
    else
    {
      TEST_REAL_EQUAL(sne.getSignalToNoise(first),0);
    }
  }
RESULT

CHECK(int getBucketSize())
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.getBucketSize() = 4;

  TEST_EQUAL(sne.getBucketSize(),4);
RESULT

CHECK(int getWindowSize())
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.getWindowSize() = 48;

  TEST_EQUAL(sne.getWindowSize(),48);
RESULT

CHECK(void setBucketSize(const int bucket_size))
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.setBucketSize(4);

  TEST_EQUAL(sne.getBucketSize(),4);
RESULT

CHECK(void setWindowSize(const int window_size))
  DSignalToNoiseEstimatorWindowing<1> sne;
  sne.setWindowSize(48);

  TEST_EQUAL(sne.getWindowSize(),48);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

