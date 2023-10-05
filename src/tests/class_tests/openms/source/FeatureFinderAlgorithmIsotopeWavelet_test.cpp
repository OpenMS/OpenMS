// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>

START_TEST(FeatureFinderAlgorithmIsotopeWavelet, "$Id$")

using namespace OpenMS;
using namespace std;

typedef FeatureFinderAlgorithmIsotopeWavelet FFASS;

FFASS* ptr = nullptr;
FFASS* nullPointer = nullptr;
FeatureFinderAlgorithm* ffA_nullPointer = nullptr;

START_SECTION((FeatureFinderAlgorithmIsotopeWavelet()))
	ptr = new FFASS;
  TEST_NOT_EQUAL(ptr,nullPointer)
END_SECTION

START_SECTION(IsotopeWaveletTransform<PeakType>::TransSpectrum* prepareHRDataCuda(const UInt i, IsotopeWaveletTransform< PeakType > *iwt))
	NOT_TESTABLE
END_SECTION

START_SECTION(MSSpectrum* createHRData(const UInt i))
	NOT_TESTABLE
END_SECTION

START_SECTION(virtual ~FeatureFinderAlgorithmIsotopeWavelet())
	delete ptr;
END_SECTION

START_SECTION(void run())
	NOT_TESTABLE
END_SECTION

START_SECTION((static FeatureFinderAlgorithm<PeakType>* create()))
  FeatureFinderAlgorithm* ptr2 = FFASS::create();
  TEST_NOT_EQUAL(ptr2,ffA_nullPointer)
	delete ptr2;
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(FFASS::getProductName(),"isotope_wavelet")
END_SECTION

END_TEST
