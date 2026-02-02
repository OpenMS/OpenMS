// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OpenSearchModificationAnalysis, "$Id$")

/////////////////////////////////////////////////////////////

OpenSearchModificationAnalysis* ptr = nullptr;
OpenSearchModificationAnalysis* nullPointer = nullptr;

START_SECTION((OpenSearchModificationAnalysis()))
  ptr = new OpenSearchModificationAnalysis();
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
END_SECTION

// Test analyzeDeltaMassPatterns: histogram and charge counts
START_SECTION(analyzeDeltaMassPatterns)
{
  OpenSearchModificationAnalysis osma;
  PeptideIdentificationList pep_ids;

  // Create synthetic peptide IDs with known DeltaMass values
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", 15.995);  // Oxidation
    hit.setCharge(2);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", 15.995);
    hit.setCharge(3);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", 42.011);  // Acetylation
    hit.setCharge(2);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }

  auto [histogram, charge_counts] = osma.analyzeDeltaMassPatterns(pep_ids, false);

  // Check histogram has entries near 15.995 and 42.011
  bool found_oxidation = false;
  bool found_acetylation = false;
  for (const auto& [mass, count] : histogram)
  {
    if (std::abs(mass - 15.995) < 0.001)
    {
      found_oxidation = true;
      TEST_REAL_SIMILAR(count, 2.0);
    }
    if (std::abs(mass - 42.011) < 0.001)
    {
      found_acetylation = true;
      TEST_REAL_SIMILAR(count, 1.0);
    }
  }
  TEST_TRUE(found_oxidation);
  TEST_TRUE(found_acetylation);

  // Check charge counts: oxidation should have 2 unique charges
  for (const auto& [mass, count] : charge_counts)
  {
    if (std::abs(mass - 15.995) < 0.001)
    {
      TEST_EQUAL(count, 2);
    }
  }
}
END_SECTION

// Test mapDeltaMassesToModifications: mapping oxidation mass to known modification
START_SECTION(mapDeltaMassesToModifications)
{
  OpenSearchModificationAnalysis osma;

  // Build histogram with oxidation mass
  OpenSearchModificationAnalysis::DeltaMassHistogram histogram(
    OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
  histogram[15.995] = 10.0;

  OpenSearchModificationAnalysis::DeltaMassToChargeCount charge_counts(
    OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
  charge_counts[15.995] = 2;

  // Create peptide IDs to annotate
  PeptideIdentificationList pep_ids;
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", 15.995);
    hit.setCharge(2);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }

  auto summaries = osma.mapDeltaMassesToModifications(
    histogram, charge_counts, pep_ids, 0.02, false, "");

  // Should find at least one modification summary
  TEST_TRUE(!summaries.empty());

  // Check that PTM metavalue was set on the hit
  const auto& hits = pep_ids[0].getHits();
  TEST_TRUE(hits[0].metaValueExists("PTM"));
  // The PTM should not be empty (it should map to a known modification)
  String ptm = hits[0].getMetaValue("PTM");
  TEST_TRUE(!ptm.empty());
}
END_SECTION

// Test smoothing: verify smoothed histogram with sigma=0.01 produces reasonable output
START_SECTION(analyzeDeltaMassPatterns_with_smoothing)
{
  OpenSearchModificationAnalysis osma;
  PeptideIdentificationList pep_ids;

  // Create a cluster of hits around oxidation mass and a separate noise floor
  double base_mass = 15.995;
  // Strong signal: many hits at the exact same binned mass (will accumulate in histogram)
  for (int i = 0; i < 50; ++i)
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", base_mass);
    hit.setCharge(2);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }
  // Noise: scattered hits at various masses (one per bin)
  for (double m = 1.0; m < 50.0; m += 0.5)
  {
    if (std::abs(m - base_mass) < 1.0) continue; // skip near signal
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setMetaValue("DeltaMass", m);
    hit.setCharge(2);
    pid.setHits({hit});
    pep_ids.push_back(pid);
  }

  auto [histogram, charge_counts] = osma.analyzeDeltaMassPatterns(pep_ids, true);

  // After smoothing + peak finding, we should have a peak near the oxidation mass
  bool found_peak = false;
  for (const auto& [mass, count] : histogram)
  {
    if (std::abs(mass - base_mass) < 0.01)
    {
      found_peak = true;
    }
  }
  TEST_TRUE(found_peak);
}
END_SECTION

END_TEST
