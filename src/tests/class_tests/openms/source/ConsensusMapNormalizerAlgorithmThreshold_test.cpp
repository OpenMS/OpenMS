// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse$
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmThreshold.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMapNormalizerAlgorithmThreshold, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmThreshold* ptr = nullptr;
ConsensusMapNormalizerAlgorithmThreshold* nullPointer = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmThreshold())
{
  ptr = new ConsensusMapNormalizerAlgorithmThreshold();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmThreshold())
{
	delete ptr;
}
END_SECTION

START_SECTION((static std::vector<double> computeCorrelation(const ConsensusMap &map, const double &ratio_threshold)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap &map, const std::vector< double > &ratios)))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

