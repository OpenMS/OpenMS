// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSignalToNoiseEstimator
  : public SignalToNoiseEstimator< >
{
  public:
  TestSignalToNoiseEstimator()
    : SignalToNoiseEstimator< >()
  {
  }

  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf)
  : SignalToNoiseEstimator< >(bpf)
  {
  }

  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
    if (&bpf==this) return *this;

    SignalToNoiseEstimator< >::operator=(bpf);

    return *this;
  }

  protected:

  void computeSTN_(const MSSpectrum& C)
      throw() override
  {
    if (C.begin() == C.end())
    {
      std::cout << "bla";
    }
    // do nothing here...
  }

};

START_TEST(SignalToNoiseEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = nullptr;
TestSignalToNoiseEstimator* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimator()))
	ptr = new TestSignalToNoiseEstimator();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((SignalToNoiseEstimator(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec);
  TestSignalToNoiseEstimator sne_copy(sne);
	NOT_TESTABLE
END_SECTION


START_SECTION((SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec);
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
	NOT_TESTABLE
END_SECTION


START_SECTION((virtual ~SignalToNoiseEstimator()))
	delete ptr;
END_SECTION


START_SECTION((virtual void init(const Container& c)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec);
	NOT_TESTABLE
END_SECTION


START_SECTION((virtual double getSignalToNoise(const Size index)))
  // hard to do without implementing computeSTN_ properly
	NOT_TESTABLE
END_SECTION


START_SECTION((float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans, const double percentile)))
{
  // Build an interleaved MS1/MS2 experiment where the two levels have clearly
  // different intensities. The picker must sample only non-empty spectra of the
  // requested MS level, so the result is deterministic despite the random draw.
  MSExperiment exp;
  for (int i = 0; i < 20; ++i)
  {
    MSSpectrum ms1;
    ms1.setMSLevel(1);
    for (int k = 0; k < 10; ++k)
    {
      Peak1D p;
      p.setMZ(100.0 + k);
      p.setIntensity(1000.0f);
      ms1.push_back(p);
    }
    exp.addSpectrum(ms1);

    MSSpectrum ms2;
    ms2.setMSLevel(2);
    for (int k = 0; k < 10; ++k)
    {
      Peak1D p;
      p.setMZ(100.0 + k);
      p.setIntensity(1.0f);
      ms2.push_back(p);
    }
    exp.addSpectrum(ms2);

    MSSpectrum empty_ms1;
    empty_ms1.setMSLevel(1);
    exp.addSpectrum(empty_ms1);
  }

  // all MS1 spectra carry intensity 1000.0 -> any sampled MS1 scan yields 1000.0
  TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(exp, 1, 10, 80), 1000.0)
  // all MS2 spectra carry intensity 1.0 -> any sampled MS2 scan yields 1.0
  TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(exp, 2, 10, 80), 1.0)
  // empty experiment must return 0.0 (no matching scans)
  TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(MSExperiment(), 1, 10, 80), 0.0)
  // requesting no scans must return 0.0 instead of dividing by zero
  TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(exp, 1, 0, 80), 0.0)

  // A sole eligible spectrum at the final experiment index must be sampled.
  MSExperiment final_candidate_exp;
  final_candidate_exp.addSpectrum(exp[1]);  // wrong MS level
  final_candidate_exp.addSpectrum(exp[2]);  // empty requested MS level
  final_candidate_exp.addSpectrum(exp[0]);  // sole eligible spectrum
  TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(final_candidate_exp, 1, 1, 100), 1000.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
