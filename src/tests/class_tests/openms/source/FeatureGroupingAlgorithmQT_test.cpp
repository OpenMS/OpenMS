// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmQT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmQT* ptr = nullptr;
FeatureGroupingAlgorithmQT* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmQT()))
	ptr = new FeatureGroupingAlgorithmQT();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmQT()))
	delete ptr;
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap >& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

START_SECTION(([EXTRA] group() with only empty FeatureMaps and use_identifications does not crash))
{
  // Regression test for #10310: ProteomicsLFQ links the runs of a fraction with
  // use_identifications=true. When no run of the fraction keeps an identification (e.g. none
  // passes the FDR filter), feature detection finds nothing and every map is empty. The
  // QTClusterFinder then read the first element of the empty mass range (segfault).
  std::vector<FeatureMap> maps(3); // three feature-empty maps
  // metadata carried by the (feature-empty) maps must NOT be silently dropped:
  // postprocess_() still transfers protein / unassigned IDs.
  for (Size i = 0; i < maps.size(); ++i)
  {
    ProteinIdentification prot;
    prot.setIdentifier("run" + std::to_string(i));
    maps[i].getProteinIdentifications().push_back(prot);
    PeptideIdentification upep;
    upep.setIdentifier("run" + std::to_string(i));
    maps[i].getUnassignedPeptideIdentifications().push_back(upep);
  }

  FeatureGroupingAlgorithmQT algo;
  Param param = algo.getParameters();
  param.setValue("use_identifications", "true");
  algo.setParameters(param);

  ConsensusMap out;
  algo.group(maps, out);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getProteinIdentifications().size(), 3)
  TEST_EQUAL(out.getUnassignedPeptideIdentifications().size(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



