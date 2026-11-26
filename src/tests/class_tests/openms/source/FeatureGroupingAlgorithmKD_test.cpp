// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureGroupingAlgorithmKD, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmKD* ptr = nullptr;
FeatureGroupingAlgorithmKD* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmKD()))
  ptr = new FeatureGroupingAlgorithmKD();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmKD()))
  delete ptr;
END_SECTION

START_SECTION((virtual void group(const std::vector<FeatureMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

START_SECTION((group() with two empty FeatureMaps does not crash))
{
  // Regression test: two input maps with no features yield an empty mass range.
  // Previously group_() dereferenced the empty mass range (front()/back()/size()-1)
  // -> crash/UB. After the fix an early guard runs postprocess_() and returns.
  FeatureGroupingAlgorithmKD algo;
  std::vector<FeatureMap> maps(2); // two feature-empty maps
  // metadata carried by the (feature-empty) maps must NOT be silently dropped by the
  // empty-input early return: postprocess_() still transfers protein / unassigned IDs.
  ProteinIdentification prot;
  prot.setIdentifier("run0");
  maps[0].getProteinIdentifications().push_back(prot);
  PeptideIdentification upep;
  upep.setIdentifier("run0");
  maps[0].getUnassignedPeptideIdentifications().push_back(upep);

  ConsensusMap out;
  algo.group(maps, out);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getProteinIdentifications().size(), 1)
  TEST_EQUAL(out.getUnassignedPeptideIdentifications().size(), 1)
}
END_SECTION

START_SECTION((group() ranks link partners by their warped RTs))
{
  // Map 1 is shifted by +40 s relative to map 0. 60 compounds give clean anchors
  // for the internal LOWESS warp. For 10 more compounds, map 1 also contains a decoy
  // at the unshifted RT (and almost the same m/z). Both the true partner and the
  // decoy fall into the link window, so the RT distance decides between them:
  // after warping, the true partner is the closer one; on the original RTs the
  // decoy would be.
  FeatureGroupingAlgorithmKD algo;
  Param p = algo.getParameters();
  p.setValue("mz_unit", "Da");
  p.setValue("link:mz_tol", 0.01);
  p.setValue("link:rt_tol", 60.0);
  p.setValue("warp:enabled", "true");
  p.setValue("warp:mz_tol", 0.01);
  p.setValue("warp:rt_tol", 100.0);
  algo.setParameters(p);

  const Size n_anchor = 60;
  const Size n_decoy = 10;
  std::vector<FeatureMap> maps(2);
  for (Size i = 0; i < n_anchor + n_decoy; ++i)
  {
    const double mz = 200.0 + 10.0 * i;
    const double rt = 100.0 + 30.0 * i;
    Feature f0;
    f0.setMZ(mz);
    f0.setRT(rt);
    f0.setIntensity(1e5);
    f0.setUniqueId(1000 + i);
    maps[0].push_back(f0);
    Feature f1; // true partner, shifted by +40 s
    f1.setMZ(mz);
    f1.setRT(rt + 40.0);
    f1.setIntensity(1e5);
    f1.setUniqueId(2000 + i);
    maps[1].push_back(f1);
    if (i >= n_anchor)
    {
      Feature decoy; // at the unshifted RT
      decoy.setMZ(mz + 0.001);
      decoy.setRT(rt);
      decoy.setIntensity(1e5);
      decoy.setUniqueId(3000 + i);
      maps[1].push_back(decoy);
    }
  }

  ConsensusMap out;
  algo.group(maps, out);

  Size n_true = 0;
  Size n_decoy_linked = 0;
  for (const ConsensusFeature& cf : out)
  {
    UInt64 id0 = 0;
    UInt64 id1 = 0;
    for (const FeatureHandle& fh : cf.getFeatures())
    {
      if (fh.getMapIndex() == 0)
      {
        id0 = fh.getUniqueId();
      }
      else
      {
        id1 = fh.getUniqueId();
      }
    }
    if (id0 >= 1000 + n_anchor && id0 < 2000) // compounds that have a decoy
    {
      if (id1 == id0 + 1000)
      {
        ++n_true;
      }
      else if (id1 == id0 + 2000)
      {
        ++n_decoy_linked;
      }
    }
  }
  TEST_EQUAL(n_true, n_decoy)
  TEST_EQUAL(n_decoy_linked, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
