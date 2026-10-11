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


START_SECTION([EXTRA] auto_mode 1 (AUTOMAXBYPERCENT) with a zero intensity and on an empty spectrum)
{
  // 10 peaks: one at intensity 0, eight at 1, one at 100; the window covers all of them.
  // Percentile histogram: maximum 100 -> bin size 1; intensities 0 and 1 fall into bin 0 (9 peaks),
  // 9 = int(95% of 10) peaks are reached in bin 0 -> max_intensity = 0.5 -> histogram bin size 1 (minimum).
  // Median histogram (30 bins): bin 0 holds 1, bin 1 holds 8, bin 29 holds 1; the 5th element lies in bin 1:
  // noise = 1 + (5 - 1) / 8 = 1.5.
  MSSpectrum s;
  for (Size i = 0; i < 10; ++i)
  {
    s.push_back(Peak1D(100.0 + i, i == 0 ? 0.0f : (i == 9 ? 100.0f : 1.0f)));
  }
  SignalToNoiseEstimatorMedian<MSSpectrum> sne;
  Param p = sne.getParameters();
  p.setValue("auto_mode", 1);
  p.setValue("win_len", 1000.0);
  p.setValue("min_required_elements", 1);
  p.setValue("write_log_messages", "false");
  sne.setParameters(p);
  sne.init(s);
  TEST_REAL_SIMILAR(sne.getSignalToNoise(0), 0.0)
  TEST_REAL_SIMILAR(sne.getSignalToNoise(1), 1.0 / 1.5)
  TEST_REAL_SIMILAR(sne.getSignalToNoise(9), 100.0 / 1.5)

  MSSpectrum empty;
  sne.init(empty); // must not dereference end()
  TEST_EQUAL(empty.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] intensities far above max_intensity land in the last histogram bin)
{
  // max_intensity 1 -> bin size 1 (minimum); 3e9 / 1 exceeds INT_MAX and must be clamped to the last bin (29),
  // not converted to int first. Bin 0 holds the 0.5 peak, bin 29 the nine 3e9 peaks; the 5th element lies in
  // bin 29: noise = 29 + (5 - 1) / 9.
  MSSpectrum s;
  s.push_back(Peak1D(100.0, 0.5f));
  for (Size i = 1; i < 10; ++i)
  {
    s.push_back(Peak1D(100.0 + i, 3.0e9f));
  }
  SignalToNoiseEstimatorMedian<MSSpectrum> sne;
  Param p = sne.getParameters();
  p.setValue("auto_mode", -1);
  p.setValue("max_intensity", 1);
  p.setValue("win_len", 1000.0);
  p.setValue("min_required_elements", 1);
  p.setValue("write_log_messages", "false");
  sne.setParameters(p);
  sne.init(s);
  const double noise = 29.0 + 4.0 / 9.0;
  TEST_REAL_SIMILAR(sne.getSignalToNoise(0), 0.5 / noise)
  TEST_REAL_SIMILAR(sne.getSignalToNoise(1), 3.0e9 / noise)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


