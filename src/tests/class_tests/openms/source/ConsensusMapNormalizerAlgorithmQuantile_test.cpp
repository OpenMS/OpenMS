// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmQuantile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMapNormalizerAlgorithmQuantile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMapNormalizerAlgorithmQuantile* ptr = nullptr;
ConsensusMapNormalizerAlgorithmQuantile* null_ptr = nullptr;
START_SECTION(ConsensusMapNormalizerAlgorithmQuantile())
{
	ptr = new ConsensusMapNormalizerAlgorithmQuantile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ConsensusMapNormalizerAlgorithmQuantile())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~ConsensusMapNormalizerAlgorithmQuantile()))
{
  // TODO
}
END_SECTION

START_SECTION((static void normalizeMaps(ConsensusMap &map)))
{
  // TODO
}
END_SECTION

START_SECTION((static void resample(const std::vector< double > &data_in, std::vector< double > &data_out, UInt n_resampling_points)))
{
  // TODO
}
END_SECTION

START_SECTION((static void extractIntensityVectors(const ConsensusMap &map, std::vector< std::vector< double > > &out_intensities)))
{
  // TODO
}
END_SECTION

START_SECTION((static void setNormalizedIntensityValues(const std::vector< std::vector< double > > &feature_ints, ConsensusMap &map)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



