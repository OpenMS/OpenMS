// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSpectrumCompareFunctor, "$Id$")

/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

START_SECTION(BinnedSpectrumCompareFunctor())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(~BinnedSpectrumCompareFunctor())
{
  NOT_TESTABLE
}
END_SECTION

//interface class is not testable

START_SECTION((BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual double operator()(const BinnedSpectrum &spec) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void registerChildren()))
{
  BinnedSpectrumCompareFunctor* c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSharedPeakCount");
  TEST_EQUAL(c1->getName(), "BinnedSharedPeakCount")
  delete c1;
  c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSpectralContrastAngle");
  TEST_EQUAL(c1->getName(), "BinnedSpectralContrastAngle")
  delete c1;
  c1 = Factory<BinnedSpectrumCompareFunctor>::create("BinnedSumAgreeingIntensities");
  TEST_EQUAL(c1->getName(), "BinnedSumAgreeingIntensities")
  delete c1;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(BinnedSpectrumCompareFunctor::getProductName(), "BinnedSpectrumCompareFunctor")
}
END_SECTION

START_SECTION(([BinnedSpectrumCompareFunctor::IncompatibleBinning] IncompatibleBinning(const char *file, int line, const char *function, const char *message="compared spectra have different settings in binsize and/or binspread")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([BinnedSpectrumCompareFunctor::IncompatibleBinning] virtual ~IncompatibleBinning()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



