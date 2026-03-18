// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
