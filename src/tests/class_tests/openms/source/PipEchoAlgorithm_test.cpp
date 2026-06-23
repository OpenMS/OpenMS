// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PipEchoAlgorithm.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  // Build a feature; a non-empty sequence attaches a peptide ID, turning the
  // feature into a PIP-ECHO "donor" (identified). Empty sequence -> "acceptor".
  Feature makeFeature(double rt, double mz, Int charge, double intensity,
                      const std::string& seq, UInt64 uid)
  {
    Feature f;
    f.setRT(rt);
    f.setMZ(mz);
    f.setCharge(charge);
    f.setIntensity(intensity);
    f.setUniqueId(uid);
    if (!seq.empty())
    {
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(seq));
      hit.setCharge(charge);
      PeptideIdentification pid;
      pid.setHits({hit});
      f.setPeptideIdentifications({pid});
    }
    return f;
  }

  FeatureMap makeRun(const std::string& path, const vector<Feature>& feats)
  {
    FeatureMap fm;
    for (const auto& f : feats) { fm.push_back(f); }
    fm.setPrimaryMSRunPath({path});
    return fm;
  }
}

START_TEST(PipEchoAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PipEchoAlgorithm* ptr = nullptr;
PipEchoAlgorithm* nullPointer = nullptr;
START_SECTION((PipEchoAlgorithm()))
  ptr = new PipEchoAlgorithm();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~PipEchoAlgorithm()))
  delete ptr;
END_SECTION

START_SECTION([EXTRA] parameter defaults)
{
  PipEchoAlgorithm algo;
  TEST_EQUAL(algo.getName(), "PipEcho")

  Param p = algo.getParameters();
  // m/z window defaults to ppm (must match the QT linker scale, not 10 Da).
  TEST_EQUAL(p.getValue("distance_MZ:unit").toString(), "ppm")
  TEST_REAL_SIMILAR(double(p.getValue("distance_MZ:max_difference")), 10.0)
  TEST_REAL_SIMILAR(double(p.getValue("distance_RT:max_difference")), 100.0)
  TEST_REAL_SIMILAR(double(p.getValue("fdr")), 0.05)
  // A seed parameter exists so decoy selection is reproducible.
  TEST_EQUAL(p.exists("random_seed"), true)
}
END_SECTION

START_SECTION((virtual void group(const std::vector<FeatureMap>& features, ConsensusMap& consensus)))
{
  // Two aligned runs sharing identified peptides, plus unidentified
  // (acceptor) features that match-between-runs should consider for transfer.
  vector<Feature> runA = {
    makeFeature(50.0, 500.00, 2, 1.0e6, "PEPTIDEAAK", 1),
    makeFeature(80.0, 600.00, 2, 2.0e6, "PEPTIDEBBR", 2),
    makeFeature(110.0, 700.00, 2, 3.0e6, "PEPTIDECCK", 3),
    makeFeature(80.0, 600.00, 2, 1.5e6, "", 4) // unidentified
  };
  vector<Feature> runB = {
    makeFeature(51.0, 500.00, 2, 1.1e6, "PEPTIDEAAK", 11),
    makeFeature(80.5, 600.00, 2, 1.9e6, "", 12), // acceptor near PEPTIDEBBR
    makeFeature(110.5, 700.00, 2, 2.9e6, "", 13) // acceptor near PEPTIDECCK
  };

  vector<FeatureMap> maps = {makeRun("runA.mzML", runA),
                             makeRun("runB.mzML", runB)};

  PipEchoAlgorithm algo;
  Param param = algo.getParameters();
  param.setValue("fdr", 1.0);
  algo.setParameters(param);

  ConsensusMap consensus;
  algo.group(maps, consensus);

  // Every identified donor peptide yields at least one consensus feature, so
  // the output is never empty.
  TEST_EQUAL(consensus.empty(), false)

  Size transfer_count = 0;
  for (const ConsensusFeature& cf : consensus)
  {
    if (cf.metaValueExists("mbr_transfer_map_indices"))
    {
      TEST_EQUAL(cf.metaValueExists("mbr_transfer_qvalues"), true)

      IntList map_indices = cf.getMetaValue("mbr_transfer_map_indices");
      DoubleList q_values = cf.getMetaValue("mbr_transfer_qvalues");
      TEST_EQUAL(map_indices.size(), q_values.size())

      for (Size i = 0; i < map_indices.size(); ++i)
      {
        TEST_EQUAL(map_indices[i], 1)
        TEST_REAL_SIMILAR(q_values[i], 0.5)
      }

      transfer_count += map_indices.size();
    }
  }
  TEST_EQUAL(transfer_count, 2)

  // Determinism: with the default fixed seed, a second run on identical input
  // must produce an identical consensus map (guards the seeded-RNG contract).
  PipEchoAlgorithm algo2;
  algo2.setParameters(param);
  ConsensusMap consensus2;
  algo2.group(maps, consensus2);

  TEST_EQUAL(consensus.size(), consensus2.size())
  bool identical = consensus.size() == consensus2.size();
  for (Size i = 0; identical && i < consensus.size(); ++i)
  {
    identical = (consensus[i].getRT() == consensus2[i].getRT())
                && (consensus[i].getMZ() == consensus2[i].getMZ())
                && (consensus[i].size() == consensus2[i].size());
  }
  TEST_EQUAL(identical, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
