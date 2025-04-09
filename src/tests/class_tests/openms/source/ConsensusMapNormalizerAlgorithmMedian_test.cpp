// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMapNormalizerAlgorithmMedian, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmMedian* ptr = nullptr;
ConsensusMapNormalizerAlgorithmMedian* null_ptr = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmMedian())
{
	ptr = new ConsensusMapNormalizerAlgorithmMedian();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmMedian())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~ConsensusMapNormalizerAlgorithmMedian()))
{
  // TODO
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap &map)))
{
  // TODO
}
END_SECTION

START_SECTION((static std::vector<double> computeNormalizationFactors(const ConsensusMap &map)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



