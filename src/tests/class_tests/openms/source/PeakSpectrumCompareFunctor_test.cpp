// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Factory.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakSpectrumCompareFunctor, "$Id$")

/////////////////////////////////////////////////////////////

// pure interface class cannot test this

START_SECTION(PeakSpectrumCompareFunctor())
  NOT_TESTABLE
END_SECTION

START_SECTION(PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source))
  NOT_TESTABLE
END_SECTION

START_SECTION(~PeakSpectrumCompareFunctor())
	NOT_TESTABLE
END_SECTION

START_SECTION(PeakSpectrumCompareFunctor& operator = (const PeakSpectrumCompareFunctor& source))
  NOT_TESTABLE
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
  NOT_TESTABLE
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
  NOT_TESTABLE
END_SECTION

START_SECTION(static void registerChildren())
  PeakSpectrumCompareFunctor* c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumCheapDPCorr");
  delete c1;
  c1 = Factory<PeakSpectrumCompareFunctor>::create("SpectrumPrecursorComparator");
	TEST_EQUAL(c1->getName(), "SpectrumPrecursorComparator")
	delete c1;
	c1 = Factory<PeakSpectrumCompareFunctor>::create("ZhangSimilarityScore");
	TEST_EQUAL(c1->getName(), "ZhangSimilarityScore")
	delete c1;
	c1 = Factory<PeakSpectrumCompareFunctor>::create("SteinScottImproveScore");
	TEST_EQUAL(c1->getName(), "SteinScottImproveScore");
	delete c1;
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(PeakSpectrumCompareFunctor::getProductName(), "PeakSpectrumCompareFunctor")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
