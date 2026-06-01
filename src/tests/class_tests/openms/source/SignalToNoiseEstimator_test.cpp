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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


