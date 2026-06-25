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

#include <set>

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

  // Attach ion-mobility annotations (as ProteomicsLFQ keeps them) so the
  // feature participates in the optional IM scoring feature.
  Feature withIM(Feature f, double im_median, double im_min, double im_max)
  {
    f.setMetaValue("IM_median", im_median);
    f.setMetaValue("IM_min", im_min);
    f.setMetaValue("IM_max", im_max);
    return f;
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
  // FDR-estimability floor: transfers are dropped below this many decoys.
  TEST_EQUAL(int(p.getValue("min_decoys")), 20)
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
  // fdr=1.0 disables FDR control (the escape hatch): all transfers are kept
  // regardless of how few decoys this tiny input generates, so we can exercise
  // the transfer-provenance meta values without needing >= min_decoys decoys.
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

START_SECTION([EXTRA] conservative fallback drops transfers when the FDR is not estimable)
{
  // Same tiny input as above, which generates 0 decoy transfers.  At an FDR
  // cutoff < 1.0 the transfer FDR cannot be resolved (decoys < min_decoys),
  // so the conservative fallback must drop ALL transferred features while
  // keeping every direct identification (donor).
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
  // Use fdr=0.5 so this test is NOT vacuous: with 0 decoys the q-value of both
  // transfers collapses to 0.5, which the ordinary q-value filter would KEEP
  // (0.5 <= 0.5). So the only thing that can remove them is the conservative
  // FDR-estimability gate (0 decoys < min_decoys, and fdr < 1.0). A tighter
  // cutoff like 0.01 would pass even without the gate, hiding regressions.
  param.setValue("fdr", 0.5);
  algo.setParameters(param);

  ConsensusMap consensus;
  algo.group(maps, consensus);

  // Direct identifications (donors) are retained in full: the three donor
  // peptide groups (PEPTIDEAAK across both runs, PEPTIDEBBR, PEPTIDECCK) yield
  // three consensus features with four feature handles total (2 + 1 + 1) and
  // NOT one transferred handle.  Contrast the escape-hatch run above, where the
  // two transfers raise the handle count to six -- so this also proves the
  // transfers were actually dropped, not merely that donors survived.
  TEST_EQUAL(consensus.empty(), false)
  TEST_EQUAL(consensus.size(), 3)

  Size handle_count = 0;
  for (const ConsensusFeature& cf : consensus) { handle_count += cf.size(); }
  TEST_EQUAL(handle_count, 4)

  // No transferred features survive, so no consensus feature carries the MBR
  // transfer provenance meta values.
  for (const ConsensusFeature& cf : consensus)
  {
    TEST_EQUAL(cf.metaValueExists("mbr_transfer_map_indices"), false)
    TEST_EQUAL(cf.metaValueExists("mbr_transfer_qvalues"), false)
  }
}
END_SECTION

START_SECTION([EXTRA] ion mobility disambiguates competing acceptors when present)
{
  // When the data carries ion mobility, PIP-ECHO uses it as a 4th scoring
  // feature. Here a single donor (PEPTIDEBBR, in runA) has two equally good
  // candidate acceptors in runB -- identical RT, m/z and intensity -- differing
  // ONLY in ion mobility: one matches the donor's mobility, the other is far
  // off. The IM-matching acceptor must win the transfer.
  //
  // Note: RunStatistics (and thus the IM tolerance) is built from the ACCEPTOR
  // run's identified donors, so runB carries two IM-annotated donors to give a
  // data-driven IM tolerance.
  const UInt64 BBR_DONOR = 2, X_MATCH = 12, Y_MISMATCH = 13;

  // Ion mobility is enabled only when EVERY run has >= 2 IM-annotated donors;
  // runA gets extra IM donors (unrelated to the BBR transfer). Three, not two,
  // to also clear the pre-existing SummaryStatistics size-2 edge case (#9659).
  vector<Feature> runA = {
    withIM(makeFeature(80.0, 600.00, 2, 2.0e6, "PEPTIDEBBR", BBR_DONOR), 1.00, 0.99, 1.01),
    withIM(makeFeature(30.0, 450.00, 2, 1.0e6, "PEPTIDEZZK", 3), 0.85, 0.845, 0.855),
    withIM(makeFeature(130.0, 650.00, 2, 1.2e6, "PEPTIDEWWK", 5), 1.05, 1.045, 1.055)
  };
  vector<Feature> runB = {
    // >= 3 identified donors (with IM widths) so an IM tolerance can be built.
    // (Three, not two: Math::SummaryStatistics in the run statistics has a
    // pre-existing empty-range crash at exactly 2 values; real runs have many.)
    withIM(makeFeature(50.0, 500.00, 2, 1.0e6, "PEPTIDEAAK", 11), 0.90, 0.895, 0.905),
    withIM(makeFeature(110.0, 700.00, 2, 3.0e6, "PEPTIDECCK", 14), 1.10, 1.095, 1.105),
    withIM(makeFeature(65.0, 550.00, 2, 1.4e6, "PEPTIDEDDK", 15), 0.95, 0.945, 0.955),
    // Competing acceptors for PEPTIDEBBR: identical but for ion mobility.
    withIM(makeFeature(80.5, 600.00, 2, 2.0e6, "", X_MATCH), 1.00, 0.99, 1.01),  // matches donor IM
    withIM(makeFeature(80.5, 600.00, 2, 2.0e6, "", Y_MISMATCH), 1.30, 1.29, 1.31) // far-off IM
  };

  vector<FeatureMap> maps = {makeRun("runA.mzML", runA),
                             makeRun("runB.mzML", runB)};

  PipEchoAlgorithm algo;
  Param param = algo.getParameters();
  param.setValue("fdr", 1.0); // escape hatch: keep the transfer regardless of decoy count
  algo.setParameters(param);

  ConsensusMap consensus;
  algo.group(maps, consensus);

  // Locate the consensus group anchored by the PEPTIDEBBR donor and collect the
  // unique ids of its feature handles.
  bool found_bbr_group = false;
  bool has_match = false, has_mismatch = false;
  for (const ConsensusFeature& cf : consensus)
  {
    std::set<UInt64> uids;
    for (const FeatureHandle& h : cf) { uids.insert(h.getUniqueId()); }
    if (uids.count(BBR_DONOR) == 0) continue;

    found_bbr_group = true;
    has_match = uids.count(X_MATCH) > 0;
    has_mismatch = uids.count(Y_MISMATCH) > 0;
  }

  TEST_EQUAL(found_bbr_group, true)
  // The ion-mobility-matching acceptor is transferred; the mismatched one is not.
  TEST_EQUAL(has_match, true)
  TEST_EQUAL(has_mismatch, false)

  // The far-off-IM acceptor must not leak into ANY consensus group.
  bool mismatch_anywhere = false;
  for (const ConsensusFeature& cf : consensus)
  {
    for (const FeatureHandle& h : cf)
    {
      if (h.getUniqueId() == Y_MISMATCH) mismatch_anywhere = true;
    }
  }
  TEST_EQUAL(mismatch_anywhere, false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
