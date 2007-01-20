// // -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimatorMedian.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

///////////////////////////

START_TEST(DSignalToNoiseEstimatorMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DSignalToNoiseEstimatorMedian<1>* ptr;

// CTORS

CHECK((DSignalToNoiseEstimatorMedian( double win_len = (double) DEFAULT_WINLEN, int bin_count = DEFAULT_BINCOUNT, int min_required_elements = DEFAULT_MIN_REQUIRED_ELEMENTS, double noise_for_empty_window= (double) DEFAULT_NOISE_ON_EMTPY_WINDOW)))
  // include std ctor here (checker.php does not recognize this as default ctor)
  ptr = new DSignalToNoiseEstimatorMedian<1>;
  TEST_NOT_EQUAL(ptr, 0)

  DSignalToNoiseEstimatorMedian<1> sne(31.1, 33);
  TEST_REAL_EQUAL(sne.getWinLen(), 31.1);
  TEST_EQUAL(sne.getBinCount(), 33);
  TEST_EQUAL(sne.getMinReqElements(), DSignalToNoiseEstimatorMedian<>::DEFAULT_MIN_REQUIRED_ELEMENTS);
  TEST_REAL_EQUAL(sne.getNoiseForEmtpyWindow(), DSignalToNoiseEstimatorMedian<>::DEFAULT_NOISE_ON_EMTPY_WINDOW);
  TEST_REAL_EQUAL(sne.getMaxIntensity(), -1);
  TEST_REAL_EQUAL(sne.getAutoMaxIntensity(), DSignalToNoiseEstimatorMedian<>::DEFAULT_MAXINTENSITY_BYSTDEV);
  TEST_REAL_EQUAL(sne.getAutoMode(), DSignalToNoiseEstimatorMedian<>::AUTOMAXBYSTDEV);
RESULT

CHECK(~DSignalToNoiseEstimatorMedian())
  delete ptr;
RESULT

CHECK((DSignalToNoiseEstimatorMedian( double win_len, int bin_count, int min_required_elements, double noise_for_empty_window, double auto_max_intensity, int auto_mode = AUTOMAXBYSTDEV)))
  DSignalToNoiseEstimatorMedian<1> sne2(31.1, 33, 11, 2.1, 190);
  TEST_REAL_EQUAL(sne2.getWinLen(), 31.1);
  TEST_EQUAL(sne2.getBinCount(), 33);
  TEST_EQUAL(sne2.getMinReqElements(), 11);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2.1);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 190);
  TEST_REAL_EQUAL(sne2.getAutoMaxIntensity(), DSignalToNoiseEstimatorMedian<>::DEFAULT_MAXINTENSITY_BYSTDEV);
  TEST_REAL_EQUAL(sne2.getAutoMode(), DSignalToNoiseEstimatorMedian<>::AUTOMAXBYSTDEV);
RESULT

CHECK((DSignalToNoiseEstimatorMedian( double win_len, int bin_count, int min_required_elements, double noise_for_empty_window, int max_intensity)))
  DSignalToNoiseEstimatorMedian<1> sne3(31.1, 33, 11, 2.1, 95.4, 1);
  TEST_REAL_EQUAL(sne3.getWinLen(), 31.1);
  TEST_EQUAL(sne3.getBinCount(), 33);
  TEST_EQUAL(sne3.getMinReqElements(), 11);
  TEST_REAL_EQUAL(sne3.getNoiseForEmtpyWindow(), 2.1);
  TEST_REAL_EQUAL(sne3.getMaxIntensity(), -1);
  TEST_REAL_EQUAL(sne3.getAutoMaxIntensity(), 95.4);
  TEST_REAL_EQUAL(sne3.getAutoMode(), 1);
RESULT


CHECK(DSignalToNoiseEstimatorMedian(const Param& parameters))
  Param p;
  p.setValue("WindowLength", 31);
  p.setValue("BinCount", 33);
  p.setValue("MinReqElementsInWindow", 4);
  p.setValue("NoiseEmptyWindow", 2);
  p.setValue("MaxIntensity", 45);
  p.setValue("AutoMode", 0);
  p.setValue("AutoMaxIntensity", 3);

  DSignalToNoiseEstimatorMedian<1> sne(p);
  TEST_REAL_EQUAL(sne.getWinLen(), 31);
  TEST_REAL_EQUAL(sne.getBinCount(), 33);
  TEST_REAL_EQUAL(sne.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne.getAutoMaxIntensity(), 3);
  TEST_REAL_EQUAL(sne.getAutoMode(), 0);
  TEST_REAL_EQUAL(sne.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(DSignalToNoiseEstimatorMedian(const DSignalToNoiseEstimatorMedian& ne))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setMaxIntensity(45);
  sne.setAutoMaxIntensity(3);
  sne.setAutoMode(0);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  DSignalToNoiseEstimatorMedian<1> sne2(sne);
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxIntensity(), 3);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 0);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(void init(PeakIterator it_begin, PeakIterator it_end))

  DPeakArray<1,DRawDataPoint<1> > raw_data(2);

  DSignalToNoiseEstimatorMedian<1,  DPeakArray<1,DRawDataPoint<1> >::const_iterator > sne;
  sne.init(raw_data.begin(),raw_data.end());

  TEST_EQUAL((sne.getFirstDataPoint() == raw_data.begin()),true);
  TEST_EQUAL((sne.getLastDataPoint()  == raw_data.end()),true);

RESULT

CHECK(DSignalToNoiseEstimatorMedian& operator=(const DSignalToNoiseEstimatorMedian& ne))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setMaxIntensity(45);
  sne.setAutoMaxIntensity(3);
  sne.setAutoMode(0);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  DSignalToNoiseEstimatorMedian<1> sne2;
  sne2 = sne;
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxIntensity(), 3);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 0);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(const double& getWinLen() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getWinLen(), DSignalToNoiseEstimatorMedian<>::DEFAULT_WINLEN);
RESULT

CHECK(const int& getBinCount() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getBinCount(), DSignalToNoiseEstimatorMedian<>::DEFAULT_BINCOUNT);
RESULT

CHECK(const double& getMaxIntensity() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getMaxIntensity(), -1);
RESULT

CHECK(const double& getAutoMaxIntensity() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getAutoMaxIntensity(), DSignalToNoiseEstimatorMedian<>::DEFAULT_MAXINTENSITY_BYSTDEV);
RESULT

CHECK(const int& getAutoMode() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getAutoMode(), DSignalToNoiseEstimatorMedian<>::AUTOMAXBYSTDEV);
RESULT

CHECK(const int& getMinReqElements() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getMinReqElements(), DSignalToNoiseEstimatorMedian<>::DEFAULT_MIN_REQUIRED_ELEMENTS);
RESULT

CHECK(const double& getNoiseForEmtpyWindow() const)
  const DSignalToNoiseEstimatorMedian<1> sne;
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), DSignalToNoiseEstimatorMedian<>::DEFAULT_NOISE_ON_EMTPY_WINDOW);
RESULT



CHECK(double getSignalToNoise(PeakIterator data_point))

  // A container for the raw data 
  MSSpectrum <DRawDataPoint <1> > raw_data;
  DTAFile dta_file;
  dta_file.load("./data/DSignalToNoise_test.dta", raw_data);
  
    
  DSignalToNoiseEstimatorMedian<1,  MSSpectrum <DRawDataPoint <1> >::const_iterator > sne(40);  //winLen of 40 Th
  sne.setMinReqElements(10);
  sne.setNoiseForEmtpyWindow(2);
  sne.init(raw_data.begin(),raw_data.end());

  DPeakArray<1,DRawDataPoint<1> >::const_iterator first = raw_data.begin() + 1;
  TEST_REAL_EQUAL(sne.getSignalToNoise(first), 1.07393);

  first = raw_data.begin() + 15;
  TEST_REAL_EQUAL(sne.getSignalToNoise(first), 2.14786);


RESULT

// MUTABLE GET-METHODS  

CHECK(double& getWinLen())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getWinLen() = 31;
  TEST_REAL_EQUAL(sne.getWinLen(), 31);
RESULT

CHECK(int& getBinCount())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getBinCount() = 50;
  TEST_EQUAL(sne.getBinCount(), 50);
RESULT

CHECK(double& getMaxIntensity())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getMaxIntensity() = 110.5;
  TEST_REAL_EQUAL(sne.getMaxIntensity(), 110.5);
RESULT

CHECK(double& getAutoMaxIntensity())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getAutoMaxIntensity() = 3.35;
  TEST_REAL_EQUAL(sne.getAutoMaxIntensity(), 3.35);
RESULT

CHECK(int& getAutoMode())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getAutoMode() = DSignalToNoiseEstimatorMedian<>::AUTOMAXBYPERCENT;
  TEST_EQUAL(sne.getAutoMode(), DSignalToNoiseEstimatorMedian<>::AUTOMAXBYPERCENT);
RESULT

CHECK(int& getMinReqElements())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getMinReqElements() = 11;
  TEST_EQUAL(sne.getMinReqElements(), 11);
RESULT


CHECK(double& getNoiseForEmtpyWindow())
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.getNoiseForEmtpyWindow() = 2.5;
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), 2.5);
RESULT

// SET METHODS

CHECK(void setWinLen(const double& win_len))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setWinLen(31);
  TEST_REAL_EQUAL(sne.getWinLen(), 31);
RESULT

CHECK(void setBinCount(const int& bin_count))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setBinCount(50);
  TEST_EQUAL(sne.getBinCount(), 50);
RESULT

CHECK(void setMaxIntensity(const double& max_intensity))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setMaxIntensity(110.5);
  TEST_REAL_EQUAL(sne.getMaxIntensity(), 110.5);
RESULT

CHECK(void setAutoMaxIntensity(const double& auto_max_intensity))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setAutoMaxIntensity(3.35);
  TEST_REAL_EQUAL(sne.getAutoMaxIntensity(), 3.35);
RESULT

CHECK(void setAutoMode(const int& auto_mode))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setAutoMode(DSignalToNoiseEstimatorMedian<>::AUTOMAXBYPERCENT);
  TEST_EQUAL(sne.getAutoMode(), DSignalToNoiseEstimatorMedian<>::AUTOMAXBYPERCENT);
RESULT

CHECK(void setMinReqElements(const int& min_required_elements))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setMinReqElements(11);
  TEST_EQUAL(sne.getMinReqElements(), 11);
RESULT


CHECK(void setNoiseForEmtpyWindow(const double& noise_for_empty_window))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setNoiseForEmtpyWindow(2.5);
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), 2.5);
RESULT

/// overwritten base-members

CHECK(void setFirstDataPoint(const PeakIterator& first))
  DPeakArray<1,DRawDataPoint<1> > raw_data;
  DRawDataPoint<1> p;
  DPosition<1> pos = 130;
  p.setIntensity(100);
  p.setPosition(pos);
  raw_data.push_back(p);

  DSignalToNoiseEstimatorMedian<1,  DPeakArray<1,DRawDataPoint<1> >::const_iterator > sne;
  sne.setFirstDataPoint(raw_data.begin());
  //TEST_EQUAL(sne.getFirstDataPoint(), raw_data.begin());
RESULT

CHECK(void setLastDataPoint(const PeakIterator& last))
  DPeakArray<1,DRawDataPoint<1> > raw_data;
  DRawDataPoint<1> p;
  DPosition<1> pos = 130;
  p.setIntensity(100);
  p.setPosition(pos);
  raw_data.push_back(p);

  DSignalToNoiseEstimatorMedian<1,  DPeakArray<1,DRawDataPoint<1> >::const_iterator > sne;
  sne.setLastDataPoint(raw_data.end());
  //TEST_EQUAL(sne.getLastDataPoint(), raw_data.end());
RESULT

CHECK(void setMZdim(const int& mz_dim))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setMZdim(1);
  TEST_EQUAL(sne.getMZdim(), 1);
RESULT

CHECK(void setRTdim(const int& rt_dim))
  DSignalToNoiseEstimatorMedian<1> sne;
  sne.setRTdim(-1);
  TEST_EQUAL(sne.getRTdim(), -1);
RESULT

CHECK(void setParam(const Param& param))
  DSignalToNoiseEstimatorMedian<1> sne;
  Param param_;
  param_.setValue("WindowLength", 32);
  param_.setValue("BinCount", (double) 21);
  param_.setValue("MinReqElementsInWindow", (double) 4);
  param_.setValue("NoiseEmptyWindow", 2.1);
  param_.setValue("MaxIntensity", 1400);
  param_.setValue("AutoMode", (double) 1);
  param_.setValue("AutoMaxIntensity", 97);
  sne.setParam(param_);    
      
  TEST_EQUAL(sne.getParam(), param_);
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

