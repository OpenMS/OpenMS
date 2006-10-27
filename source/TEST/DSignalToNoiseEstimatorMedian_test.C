// // -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMedian.h>

///////////////////////////

START_TEST(DSignalToNoiseEstimatorMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DSignalToNoiseEstimatorMedian<1>* ptr;
CHECK(DSignalToNoiseEstimatorMedian())
  ptr = new DSignalToNoiseEstimatorMedian<1>;
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DSignalToNoiseEstimatorMedian())
  delete ptr;
RESULT

CHECK(DSignalToNoiseEstimatorMedian(const Param& parameters))
  Param p;
  p.setValue("SignalToNoiseEstimationParameter:Window",20);
  p.setValue("SignalToNoiseEstimationParameter:Median_perc",0.5);

  DSignalToNoiseEstimatorMedian<1> sne(p);
  TEST_REAL_EQUAL(sne.getWindowSize(),20);
  TEST_REAL_EQUAL(sne.getFactor(),0.5);
RESULT

CHECK(DSignalToNoiseEstimatorMedian(const DSignalToNoiseEstimatorMedian& ne))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setFactor(0.5);
  sne.setWindowSize(20);

  DSignalToNoiseEstimatorMedian<1> sne2(sne);
  TEST_EQUAL(sne2.getFactor(),0.5);
  TEST_REAL_EQUAL(sne2.getWindowSize(),20);
RESULT

CHECK((void init(ConstIterator it_begin, ConstIterator it_end)))

  DPeakArray<1,DRawDataPoint<1> > raw_data(2);

  DSignalToNoiseEstimatorMedian<1,  DPeakArray<1,DRawDataPoint<1> >::const_iterator > sne;
  sne.init(raw_data.begin(),raw_data.end());

  TEST_EQUAL((sne.getFirstDataPoint() == raw_data.begin()),true);
  TEST_EQUAL((sne.getLastDataPoint() == raw_data.end()),true);

RESULT

CHECK(DSignalToNoiseEstimatorMedian& operator=(const DSignalToNoiseEstimatorMedian& ne))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setFactor(0.5);
  sne.setWindowSize(20);

  DSignalToNoiseEstimatorMedian<1> sne2;
  sne2 = sne;
  TEST_EQUAL(sne2.getFactor(),0.5);
  TEST_REAL_EQUAL(sne2.getWindowSize(),20);
RESULT

CHECK(const int getFactor() const)
  const DSignalToNoiseEstimatorMedian<1> sne;

  TEST_EQUAL(sne.getFactor(),1);
RESULT

CHECK(const int getWindowSize() const)
  const DSignalToNoiseEstimatorMedian<1> sne;

  TEST_EQUAL(sne.getWindowSize(),100);
RESULT

CHECK(double getSignalToNoise(ConstIterator data_point) throw(Exception::OutOfRange))
  DPeakArray<1,DRawDataPoint<1> > raw_data;
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

  DSignalToNoiseEstimatorMedian<1,  DPeakArray<1,DRawDataPoint<1> >::const_iterator > sne;
  sne.setFactor(1);
  sne.setWindowSize(1);
  sne.init(raw_data.begin(),raw_data.end());

  for (i=0; i < (int)raw_data.size(); ++i)
  {
    DPeakArray<1,DRawDataPoint<1> >::const_iterator first = raw_data.begin() + i;

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

CHECK(int getFactor())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getFactor() = 0.5;

  TEST_REAL_EQUAL(sne.getFactor(),0.5);
RESULT

CHECK(int getWindowSize())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getWindowSize() = 48;

  TEST_EQUAL(sne.getWindowSize(),48);
RESULT

CHECK(void setWinodowSize(const int window_size))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setWindowSize(4);

  TEST_EQUAL(sne.getWindowSize(),4);
RESULT

CHECK(void setFactor(const float factor))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setFactor(0.5);

  TEST_EQUAL(sne.getFactor(),0.5);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

