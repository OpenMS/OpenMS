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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SignalToNoiseEstimatorMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMedian* ptr = 0;
CHECK(SignalToNoiseEstimatorMedian())
        ptr = new SignalToNoiseEstimatorMedian();
        TEST_NOT_EQUAL(ptr, 0)

        SignalToNoiseEstimatorMedian sne;
        TEST_REAL_EQUAL(sne.getWinLen(), 200);
        TEST_EQUAL(sne.getBinCount(), 30);
        TEST_EQUAL(sne.getMinReqElements(), 10);
        TEST_REAL_EQUAL(sne.getNoiseForEmtpyWindow(), 2);
        TEST_REAL_EQUAL(sne.getMaxIntensity(), -1);
        TEST_REAL_EQUAL(sne.getAutoMode(), 0);
        TEST_REAL_EQUAL(sne.getAutoMaxPercentile(), 95);
        TEST_REAL_EQUAL(sne.getAutoMaxStdevFactor(), 3);
RESULT

CHECK(~SignalToNoiseEstimatorMedian())
        delete ptr;
RESULT

CHECK(SignalToNoiseEstimatorMedian& operator=(const SignalToNoiseEstimatorMedian &source))
  SignalToNoiseEstimatorMedian sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setMaxIntensity(45);
  sne.setAutoMaxPercentile(90);
  sne.setAutoMaxStdevFactor(3.1);
  sne.setAutoMode(1);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  SignalToNoiseEstimatorMedian sne2 = sne;
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxPercentile(), 90);
  TEST_REAL_EQUAL(sne2.getAutoMaxStdevFactor(), 3.1);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 1);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(SignalToNoiseEstimatorMedian(const SignalToNoiseEstimatorMedian &source))
  SignalToNoiseEstimatorMedian sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setMaxIntensity(45);
  sne.setAutoMaxPercentile(90);
  sne.setAutoMaxStdevFactor(3.1);
  sne.setAutoMode(1);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  SignalToNoiseEstimatorMedian sne2(sne);
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxPercentile(), 90);
  TEST_REAL_EQUAL(sne2.getAutoMaxStdevFactor(), 3.1);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 1);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(virtual ~SignalToNoiseEstimatorMedian())
  // 
RESULT

CHECK(DoubleReal getMaxIntensity() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getMaxIntensity(), -1);
RESULT

CHECK(void setMaxIntensity(DoubleReal max_intensity))
  SignalToNoiseEstimatorMedian sne;
  sne.setMaxIntensity(100);
  TEST_EQUAL(sne.getMaxIntensity(), 100);
RESULT


CHECK(DoubleReal getAutoMaxStdevFactor() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getAutoMaxStdevFactor(), 3);
RESULT

CHECK(void setAutoMaxStdevFactor(DoubleReal value))
  SignalToNoiseEstimatorMedian sne;
  sne.setAutoMaxStdevFactor(100);
  TEST_EQUAL(sne.getAutoMaxStdevFactor(), 100);
RESULT


CHECK(DoubleReal getAutoMaxPercentile() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getAutoMaxPercentile(), 95);
RESULT

CHECK(void setAutoMaxPercentile(DoubleReal value))
  SignalToNoiseEstimatorMedian sne;
  sne.setAutoMaxPercentile(100);
  TEST_EQUAL(sne.getAutoMaxPercentile(), 100);
RESULT


CHECK(Int getAutoMode() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getAutoMode(), 0);
RESULT

CHECK(void setAutoMode(Int auto_mode))
  SignalToNoiseEstimatorMedian sne;
  sne.setAutoMode(100);
  TEST_EQUAL(sne.getAutoMode(), 100);
RESULT


CHECK(DoubleReal getWinLen() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getWinLen(), 200);
RESULT

CHECK(void setWinLen(DoubleReal win_len))
  SignalToNoiseEstimatorMedian sne;
  sne.setWinLen(100);
  TEST_EQUAL(sne.getWinLen(), 100);
RESULT


CHECK(Int getBinCount() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getBinCount(), 30);
RESULT

CHECK(void setBinCount(Int bin_count))
  SignalToNoiseEstimatorMedian sne;
  sne.setBinCount(100);
  TEST_EQUAL(sne.getBinCount(), 100);
RESULT


CHECK(Int getMinReqElements() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getMinReqElements(), 10);
RESULT

CHECK(void setMinReqElements(Int min_required_elements))
  SignalToNoiseEstimatorMedian sne;
  sne.setMinReqElements(100);
  TEST_EQUAL(sne.getMinReqElements(), 100);
RESULT


CHECK(DoubleReal getNoiseForEmtpyWindow() const )
  const SignalToNoiseEstimatorMedian sne;
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window))
  SignalToNoiseEstimatorMedian sne;
  sne.setNoiseForEmtpyWindow(100);
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), 100);
RESULT


CHECK(virtual double getSignalToNoise(PeakIterator data_point))
  // A container for the raw data 
  MSSpectrum < > raw_data;
  MSSpectrum< >::const_iterator it;
  DTAFile dta_file;
  dta_file.load("./data/SignalToNoiseEstimator_test.dta", raw_data);
  
    
  SignalToNoiseEstimatorMedian sne;  
  sne.setWinLen(40);
  sne.setMinReqElements(10);
  sne.setNoiseForEmtpyWindow(2);
  sne.init(raw_data.begin(),raw_data.end());

  MSSpectrum < > stn_data;
  dta_file.load("./data/SignalToNoiseEstimatorMedian_test.out", stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_EQUAL (stn_data[i].getIntensity(), sne.getSignalToNoise(it));
        
    
    //Peak1D peak = (*it);
    //peak.setIntensity(sne.getSignalToNoise(it));
    //stn_data.push_back(peak);
    ++i;
  }

  //dta_file.store("./data/SignalToNoiseEstimatorMedian_test.tmp", stn_data);
  
  //TEST_FILE("./data/SignalToNoiseEstimatorMedian_test.tmp", "./data/SignalToNoiseEstimatorMedian_test.out");


RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


