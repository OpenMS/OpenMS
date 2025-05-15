// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SignalToNoiseEstimatorMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMedian< >* ptr = nullptr;
SignalToNoiseEstimatorMedian< >* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimatorMedian()))
	ptr = new SignalToNoiseEstimatorMedian<>;
	TEST_NOT_EQUAL(ptr, nullPointer)
	SignalToNoiseEstimatorMedian<> sne;
END_SECTION

START_SECTION((SignalToNoiseEstimatorMedian& operator=(const SignalToNoiseEstimatorMedian &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMedian<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMedian<> sne2 = sne;
	NOT_TESTABLE
END_SECTION

START_SECTION((SignalToNoiseEstimatorMedian(const SignalToNoiseEstimatorMedian &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMedian<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMedian<> sne2(sne);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual ~SignalToNoiseEstimatorMedian()))
	delete ptr;
END_SECTION


START_SECTION([EXTRA](virtual void init(const Container& c)))

  MSSpectrum raw_data;
  MSSpectrum::const_iterator it;
  DTAFile dta_file;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimator_test.dta"), raw_data);
  
    
  SignalToNoiseEstimatorMedian< MSSpectrum > sne;
	Param p;
	p.setValue("win_len", 40.0);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(raw_data);

  MSSpectrum stn_data;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMedian_test.out"), stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_SIMILAR (stn_data[i].getIntensity(), sne.getSignalToNoise(i));
        
    
    //Peak1D peak = (*it);
    //peak.setIntensity(sne.getSignalToNoise(it));
    //stn_data.push_back(peak);
    ++i;
  }

  //dta_file.store("./data/SignalToNoiseEstimatorMedian_test.tmp", stn_data);
  
  //TEST_FILE_EQUAL("./data/SignalToNoiseEstimatorMedian_test.tmp", "./data/SignalToNoiseEstimatorMedian_test.out");


END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


