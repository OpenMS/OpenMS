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
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestableSignalToNoiseEstimatorMeanIterative : public SignalToNoiseEstimatorMeanIterative<MSSpectrum>
{
public:
  double maxIntensity() const { return max_intensity_; }
};

// #define DEBUG_TEST
// #undef  DEBUG_TEST

START_TEST(SignalToNoiseEstimatorMeanIterative, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMeanIterative< >* ptr = nullptr;
SignalToNoiseEstimatorMeanIterative< >* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimatorMeanIterative()))
        ptr = new SignalToNoiseEstimatorMeanIterative<>;
	TEST_NOT_EQUAL(ptr, nullPointer)
        SignalToNoiseEstimatorMeanIterative<> sne;
END_SECTION

START_SECTION([EXTRA] auto_mode 1 bins small intensities and handles a zero maximum)
{
  MSSpectrum s;
  const float intensities[] = {0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  for (Size i = 0; i < 5; ++i)
  {
    s.push_back(Peak1D(100.0 + i, intensities[i]));
  }

  TestableSignalToNoiseEstimatorMeanIterative sne;
  Param p = sne.getParameters();
  p.setValue("auto_mode", 1);
  p.setValue("auto_max_percentile", 80);
  p.setValue("win_len", 1000.0);
  p.setValue("min_required_elements", 1);
  sne.setParameters(p);
  sne.init(s);

  // Bins 0, 25, 50 and 75 contain the first four peaks; the 80th percentile
  // therefore estimates the maximum as (75 + 0.5) * 0.01.
  TEST_REAL_SIMILAR(sne.maxIntensity(), 0.755)
  for (Size i = 0; i < 5; ++i)
  {
    TEST_REAL_SIMILAR(sne.getSignalToNoise(i), intensities[i])
  }

  MSSpectrum zero_spectrum;
  for (Size i = 0; i < 5; ++i)
  {
    zero_spectrum.push_back(Peak1D(100.0 + i, 0.0f));
  }
  sne.init(zero_spectrum);
  TEST_REAL_SIMILAR(sne.maxIntensity(), 0.5)
  for (Size i = 0; i < 5; ++i)
  {
    TEST_REAL_SIMILAR(sne.getSignalToNoise(i), 0.0)
  }

  MSSpectrum empty;
  sne.init(empty);
  TEST_EQUAL(empty.size(), 0)
}
END_SECTION



START_SECTION((SignalToNoiseEstimatorMeanIterative& operator=(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2 = sne;
	NOT_TESTABLE
END_SECTION

START_SECTION((SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2(sne);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual ~SignalToNoiseEstimatorMeanIterative()))
  delete ptr;
END_SECTION


START_SECTION([EXTRA](virtual void init(const Container& c)))

	TOLERANCE_ABSOLUTE(0.5)

  MSSpectrum raw_data;
  MSSpectrum::const_iterator it;
  DTAFile dta_file;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimator_test.dta"), raw_data);

  SignalToNoiseEstimatorMeanIterative<> sne;
	Param p;
	p.setValue("win_len", 40.1);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(raw_data);

  MSSpectrum stn_data;

#ifdef DEBUG_TEST
  MSSpectrum stn_data__;
#endif

  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.out"), stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_SIMILAR (stn_data[i].getIntensity(), sne.getSignalToNoise(i));
#ifdef DEBUG_TEST
    Peak1D peak = (*it);
    peak.setIntensity(stn_data[i].getIntensity() / sne.getSignalToNoise(it));
    stn_data__.push_back(peak);
#endif
    ++i;
  }

#ifdef DEBUG_TEST
  dta_file.store(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.debug"), stn_data__);
#endif
  
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


