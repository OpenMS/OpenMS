// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaboliteSpectralMatching, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaboliteSpectralMatching* ptr = nullptr;
MetaboliteSpectralMatching* null_ptr = nullptr;
START_SECTION(MetaboliteSpectralMatching())
{
	ptr = new MetaboliteSpectralMatching();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MetaboliteSpectralMatching())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~MetaboliteSpectralMatching()))
{
  // TODO
}
END_SECTION

START_SECTION((double computeHyperScore(MSSpectrum, MSSpectrum, const double &, const double &)))
{
  // TODO
}
END_SECTION

START_SECTION((void run(PeakMap &, MzTab &)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



