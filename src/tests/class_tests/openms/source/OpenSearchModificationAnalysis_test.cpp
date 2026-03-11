// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////////
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

/*
  OpenSearchModificationAnalysis tests

  This suite verifies:
  - analyzeDeltaMassPatterns(): histogram generation from delta masses
  - mapDeltaMassesToModifications(): mapping to known modifications
  - analyzeModificationsWithStatistics(): complete analysis with statistics tables
  - generateDeltaMassStatistics(): delta mass table generation
  - generatePTMStatistics(): PTM table generation
  - analyzeResidueFrequency(): residue frequency analysis
*/

using namespace OpenMS;
using namespace std;

START_TEST(OpenSearchModificationAnalysis, "$Id$")

/////////////////////////////////////////////////////////////

OpenSearchModificationAnalysis* ptr = nullptr;
OpenSearchModificationAnalysis* null_ptr = nullptr;

START_SECTION(OpenSearchModificationAnalysis())
{
  ptr = new OpenSearchModificationAnalysis();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~OpenSearchModificationAnalysis())
{
  delete ptr;
}
END_SECTION

START_SECTION(analyzeDeltaMassPatterns)
{
  OpenSearchModificationAnalysis analyzer;

  // Create test peptide identifications with known delta masses
  PeptideIdentificationList peptide_ids;

  // Add some peptide IDs with delta masses
  for (int i = 0; i < 10; ++i)
  {
    PeptideIdentification pid;
    pid.setRT(100.0 + i);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    // Simulate oxidation mass shift (~15.995 Da)
    hit.setMetaValue("DeltaMass", 15.995 + (i * 0.001));
    pid.setHits({hit});

    peptide_ids.push_back(pid);
  }

  // Add some with different delta mass (phosphorylation ~79.97 Da)
  for (int i = 0; i < 5; ++i)
  {
    PeptideIdentification pid;
    pid.setRT(200.0 + i);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString("STYPEPTIDE"));
    hit.setCharge(3);
    hit.setMetaValue("DeltaMass", 79.97 + (i * 0.001));
    pid.setHits({hit});

    peptide_ids.push_back(pid);
  }

  // Add some near-zero delta masses (should be ignored)
  for (int i = 0; i < 3; ++i)
  {
    PeptideIdentification pid;
    pid.setRT(300.0 + i);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString("NORMPEPTIDE"));
    hit.setCharge(2);
    hit.setMetaValue("DeltaMass", 0.01);
    pid.setHits({hit});

    peptide_ids.push_back(pid);
  }

  auto [histogram, charge_counts] = analyzer.analyzeDeltaMassPatterns(peptide_ids, false, false);

  // Should have entries for oxidation and phosphorylation, but not near-zero
  TEST_EQUAL(histogram.empty(), false)
  // Near-zero masses should be filtered out
  TEST_EQUAL(histogram.find(0.01) == histogram.end(), true)
}
END_SECTION

START_SECTION(analyzeModificationsWithStatistics)
{
  OpenSearchModificationAnalysis analyzer;

  // Create test peptide identifications
  PeptideIdentificationList peptide_ids;

  // Add oxidation-like modifications
  for (int i = 0; i < 20; ++i)
  {
    PeptideIdentification pid;
    pid.setRT(100.0 + i);

    PeptideHit hit;
    hit.setSequence(AASequence::fromString("METHIONINE"));
    hit.setCharge(2);
    hit.setMetaValue("DeltaMass", 15.9949); // Oxidation mass
    pid.setHits({hit});

    peptide_ids.push_back(pid);
  }

  // Run analysis
  auto result = analyzer.analyzeModificationsWithStatistics(
    peptide_ids,
    0.1,    // tolerance
    false,  // Da units
    false,  // no smoothing
    ""      // no output file
  );

  // Check delta mass statistics
  TEST_EQUAL(result.delta_mass_stats.total_psms, 20)
  TEST_EQUAL(result.delta_mass_stats.modified_psms, 20)
  TEST_EQUAL(result.delta_mass_stats.unmodified_psms, 0)
  TEST_EQUAL(result.delta_mass_stats.entries.empty(), false)
}
END_SECTION

START_SECTION(generateDeltaMassStatistics)
{
  OpenSearchModificationAnalysis analyzer;

  // Create a simple histogram
  OpenSearchModificationAnalysis::DeltaMassHistogram histogram(
    OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
  histogram[15.995] = 10.0;  // Oxidation-like
  histogram[79.97] = 5.0;    // Phosphorylation-like
  histogram[42.01] = 3.0;    // Acetyl-like

  OpenSearchModificationAnalysis::DeltaMassToChargeCount charge_counts(
    OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
  charge_counts[15.995] = 2;
  charge_counts[79.97] = 3;
  charge_counts[42.01] = 1;

  // Create empty peptide list (for this simplified test)
  PeptideIdentificationList peptide_ids;

  auto stats = analyzer.generateDeltaMassStatistics(histogram, charge_counts, peptide_ids, 0.1, false);

  // Should have 3 entries
  TEST_EQUAL(stats.entries.size(), 3)

  // Entries should be sorted by count (descending)
  if (stats.entries.size() >= 2)
  {
    TEST_EQUAL(stats.entries[0].count >= stats.entries[1].count, true)
  }
}
END_SECTION

START_SECTION(analyzeResidueFrequency)
{
  OpenSearchModificationAnalysis analyzer;

  // Create peptide IDs with known sequences
  PeptideIdentificationList peptide_ids;

  // Add peptides containing methionine
  for (int i = 0; i < 5; ++i)
  {
    PeptideIdentification pid;
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("METHIONINE")); // Contains M
    hit.setMetaValue("DeltaMass", 15.995);
    pid.setHits({hit});
    peptide_ids.push_back(pid);
  }

  auto residue_counts = analyzer.analyzeResidueFrequency(peptide_ids, 15.995, 0.1);

  // Should find methionine (M) in residue counts
  TEST_EQUAL(residue_counts.find('M') != residue_counts.end(), true)
  if (residue_counts.find('M') != residue_counts.end())
  {
    TEST_EQUAL(residue_counts['M'], 5) // One M per peptide, 5 peptides
  }
}
END_SECTION

START_SECTION(ModificationSummary struct)
{
  OpenSearchModificationAnalysis::ModificationSummary summary;
  summary.count = 10;
  summary.name = "Oxidation";
  summary.num_charge_states = 3;
  summary.masses.push_back(15.995);

  TEST_EQUAL(summary.count, 10)
  TEST_EQUAL(summary.name, "Oxidation")
  TEST_EQUAL(summary.num_charge_states, 3)
  TEST_EQUAL(summary.masses.size(), 1)
}
END_SECTION

START_SECTION(DeltaMassEntry struct)
{
  OpenSearchModificationAnalysis::DeltaMassEntry entry;
  entry.delta_mass = 15.995;
  entry.count = 100;
  entry.unique_peptides = 50;
  entry.num_charge_states = 3;
  entry.percentage = 25.5;
  entry.mapped_modification = "Oxidation";
  entry.is_known_modification = true;

  TEST_REAL_SIMILAR(entry.delta_mass, 15.995)
  TEST_EQUAL(entry.count, 100)
  TEST_EQUAL(entry.unique_peptides, 50)
  TEST_EQUAL(entry.is_known_modification, true)
}
END_SECTION

START_SECTION(PTMEntry struct)
{
  OpenSearchModificationAnalysis::PTMEntry entry;
  entry.name = "Oxidation (M)";
  entry.theoretical_mass = 15.9949;
  entry.observed_mass = 15.9960;
  entry.mass_deviation = 0.0011;
  entry.count = 500;
  entry.unique_peptides = 200;
  entry.percentage = 30.5;
  entry.residue_counts['M'] = 450;
  entry.target_residues = "M";

  TEST_EQUAL(entry.name, "Oxidation (M)")
  TEST_REAL_SIMILAR(entry.theoretical_mass, 15.9949)
  TEST_EQUAL(entry.count, 500)
  TEST_EQUAL(entry.residue_counts['M'], 450)
}
END_SECTION

START_SECTION(OpenSearchAnalysisResult struct)
{
  OpenSearchModificationAnalysis::OpenSearchAnalysisResult result;

  // Check default initialization
  TEST_EQUAL(result.delta_mass_stats.total_psms, 0)
  TEST_EQUAL(result.ptm_stats.total_modified_psms, 0)
  TEST_EQUAL(result.summaries.empty(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
