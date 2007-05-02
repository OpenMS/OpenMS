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
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SignalToNoiseEstimatorMeanIterative, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMeanIterative< >* ptr = 0;
CHECK(SignalToNoiseEstimatorMeanIterative())
        ptr = new SignalToNoiseEstimatorMeanIterative<>;
        TEST_NOT_EQUAL(ptr, 0)

        SignalToNoiseEstimatorMeanIterative<> sne;
        TEST_REAL_EQUAL(sne.getWinLen(), 200);
        TEST_EQUAL(sne.getBinCount(), 30);
        TEST_REAL_EQUAL(sne.getSTDEVMultiplier(), 3);        
        TEST_EQUAL(sne.getMinReqElements(), 10);
        TEST_REAL_EQUAL(sne.getNoiseForEmtpyWindow(), 2);
        TEST_REAL_EQUAL(sne.getMaxIntensity(), -1);
        TEST_REAL_EQUAL(sne.getAutoMode(), 0);
        TEST_REAL_EQUAL(sne.getAutoMaxPercentile(), 95);
        TEST_REAL_EQUAL(sne.getAutoMaxStdevFactor(), 3);
RESULT

CHECK(~SignalToNoiseEstimatorMeanIterative())
        delete ptr;
RESULT

CHECK(SignalToNoiseEstimatorMeanIterative& operator=(const SignalToNoiseEstimatorMeanIterative &source))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setSTDEVMultiplier(3.4);
  sne.setMaxIntensity(45);
  sne.setAutoMaxPercentile(90);
  sne.setAutoMaxStdevFactor(3.1);
  sne.setAutoMode(1);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  SignalToNoiseEstimatorMeanIterative<> sne2 = sne;
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getSTDEVMultiplier(), 3.4);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxPercentile(), 90);
  TEST_REAL_EQUAL(sne2.getAutoMaxStdevFactor(), 3.1);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 1);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative &source))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setWinLen(31);
  sne.setBinCount(33);
  sne.setSTDEVMultiplier(3.4);
  sne.setMaxIntensity(45);
  sne.setAutoMaxPercentile(90);
  sne.setAutoMaxStdevFactor(3.1);
  sne.setAutoMode(1);
  sne.setMinReqElements(4);
  sne.setNoiseForEmtpyWindow(2);  

  SignalToNoiseEstimatorMeanIterative<> sne2(sne);
  TEST_REAL_EQUAL(sne2.getWinLen(), 31);
  TEST_REAL_EQUAL(sne2.getBinCount(), 33);
  TEST_REAL_EQUAL(sne2.getSTDEVMultiplier(), 3.4);
  TEST_REAL_EQUAL(sne2.getMaxIntensity(), 45);
  TEST_REAL_EQUAL(sne2.getAutoMaxPercentile(), 90);
  TEST_REAL_EQUAL(sne2.getAutoMaxStdevFactor(), 3.1);
  TEST_REAL_EQUAL(sne2.getAutoMode(), 1);
  TEST_REAL_EQUAL(sne2.getMinReqElements(), 4);
  TEST_REAL_EQUAL(sne2.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(virtual ~SignalToNoiseEstimatorMeanIterative())
  // 
RESULT

CHECK(DoubleReal getMaxIntensity() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getMaxIntensity(), -1);
RESULT

CHECK(void setMaxIntensity(DoubleReal max_intensity))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setMaxIntensity(100);
  TEST_EQUAL(sne.getMaxIntensity(), 100);
RESULT


CHECK(DoubleReal getAutoMaxStdevFactor() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getAutoMaxStdevFactor(), 3);
RESULT

CHECK(void setAutoMaxStdevFactor(DoubleReal value))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setAutoMaxStdevFactor(100);
  TEST_EQUAL(sne.getAutoMaxStdevFactor(), 100);
RESULT


CHECK(DoubleReal getAutoMaxPercentile() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getAutoMaxPercentile(), 95);
RESULT

CHECK(void setAutoMaxPercentile(DoubleReal value))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setAutoMaxPercentile(100);
  TEST_EQUAL(sne.getAutoMaxPercentile(), 100);
RESULT


CHECK(Int getAutoMode() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getAutoMode(), 0);
RESULT

CHECK(void setAutoMode(Int auto_mode))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setAutoMode(100);
  TEST_EQUAL(sne.getAutoMode(), 100);
RESULT


CHECK(DoubleReal getWinLen() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getWinLen(), 200);
RESULT

CHECK(void setWinLen(DoubleReal win_len))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setWinLen(100);
  TEST_EQUAL(sne.getWinLen(), 100);
RESULT


CHECK(Int getBinCount() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getBinCount(), 30);
RESULT

CHECK(void setBinCount(Int bin_count))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setBinCount(100);
  TEST_EQUAL(sne.getBinCount(), 100);
RESULT


CHECK(DoubleReal getSTDEVMultiplier() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getSTDEVMultiplier(), 3);
RESULT

CHECK(void setSTDEVMultiplier(DoubleReal stdev))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setSTDEVMultiplier(100);
  TEST_EQUAL(sne.getSTDEVMultiplier(), 100);
RESULT


CHECK(Int getMinReqElements() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getMinReqElements(), 10);
RESULT

CHECK(void setMinReqElements(Int min_required_elements))
  SignalToNoiseEstimatorMeanIterative<> sne;
  sne.setMinReqElements(100);
  TEST_EQUAL(sne.getMinReqElements(), 100);
RESULT


CHECK(DoubleReal getNoiseForEmtpyWindow() const )
  const SignalToNoiseEstimatorMeanIterative<> sne;
  TEST_EQUAL(sne.getNoiseForEmtpyWindow(), 2);
RESULT

CHECK(void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window))
  SignalToNoiseEstimatorMeanIterative<> snee;
  snee.setNoiseForEmtpyWindow(100);
  TEST_EQUAL(snee.getNoiseForEmtpyWindow(), 100);

//RESULT

PRECISION(0.5)

// this is a protected method, but we need to test it via init(), which is defined in the base class
//CHECK(virtual void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) throw(Exception::InvalidValue))
  
  MSSpectrum < > raw_data;
  MSSpectrum< >::const_iterator it;
  DTAFile dta_file;
  dta_file.load("./data/SignalToNoiseEstimator_test.dta", raw_data);
  
    
  SignalToNoiseEstimatorMeanIterative<> sne;  
  //sne.setLogType(SignalToNoiseEstimatorMeanIterative<>::CMD);  
  sne.setWinLen(40.1);
  sne.setMinReqElements(10);
  sne.setNoiseForEmtpyWindow(2);
  sne.init(raw_data.begin(),raw_data.end());

  MSSpectrum < > stn_data;
  dta_file.load("./data/SignalToNoiseEstimatorMeanIterative_test.out", stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_EQUAL (stn_data[i].getIntensity(), sne.getSignalToNoise(it));
    
    //Peak1D peak = (*it);
    //peak.setIntensity(sne.getSignalToNoise(it));
    //stn_data.push_back(peak);
    ++i;
  }

  //dta_file.store("./data/SignalToNoiseEstimatorMeanIterative_test.tmp", stn_data);
  
  //TEST_FILE("./data/SignalToNoiseEstimatorMeanIterative_test.tmp", "./data/SignalToNoiseEstimatorMeanIterative_test.out");
  
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


