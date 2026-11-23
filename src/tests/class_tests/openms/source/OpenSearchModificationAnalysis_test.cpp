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
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

/*
  OpenSearchModificationAnalysis tests

  This suite verifies:
  - analyzeDeltaMassPatterns(): histogram generation from delta masses
  - mapDeltaMassesToModifications(): mapping to known modifications (the nearest one is named,
    "A//B" only for candidates that are equally close; PTM annotation, +1Da and combination lookups follow)
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

START_SECTION(mapDeltaMassesToModifications)
{
  OpenSearchModificationAnalysis analyzer;
  const ModificationsDB* mod_db = ModificationsDB::getInstance();
  using Summary = OpenSearchModificationAnalysis::ModificationSummary;

  // Two Unimod modifications lie within 0.02 Da of +15.9949 Da: oxidation (+15.994915 Da) and, 0.0178 Da lower,
  // thiocarboxylic acid (+15.977156 Da). The latter is the entry of lowest mass in the window.
  const std::string oxidation = mod_db->getModification("Oxidation (M)")->getFullName();
  const std::string thiocarboxylic = "thiocarboxylic acid";
  const double oxidation_mass = 15.994915;
  const double thiocarboxylic_mass = 15.977156;
  {
    // The premise of this section: both are in the table that the analysis maps to (Unimod search modifications)
    std::vector<const ResidueModification*> candidates;
    mod_db->searchModificationsByDiffMonoMass(candidates, 15.9949, 0.02);
    bool has_oxidation = false, has_thiocarboxylic = false;
    for (const ResidueModification* mod : candidates)
    {
      if (mod->getUniModRecordId() <= 0) continue;
      has_oxidation |= (mod->getFullName() == oxidation);
      has_thiocarboxylic |= (mod->getFullName() == thiocarboxylic);
    }
    TEST_EQUAL(has_oxidation, true)
    TEST_EQUAL(has_thiocarboxylic, true)
  }

  // PSMs with the given delta masses
  auto make_ids = [](const std::vector<double>& delta_masses)
  {
    PeptideIdentificationList ids;
    for (double delta_mass : delta_masses)
    {
      PeptideHit hit;
      hit.setSequence(AASequence::fromString("PEPTIDEK"));
      hit.setCharge(2);
      hit.setMetaValue("DeltaMass", delta_mass);
      PeptideIdentification pid;
      pid.setHits({hit});
      ids.push_back(pid);
    }
    return ids;
  };

  // Maps clusters of 10 PSMs (one charge state) at the given masses
  auto map_clusters = [&](const std::vector<double>& cluster_masses, PeptideIdentificationList& ids, double tolerance = 0.02, bool ppm = false)
  {
    OpenSearchModificationAnalysis::DeltaMassHistogram histogram(OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
    OpenSearchModificationAnalysis::DeltaMassToChargeCount charges(OpenSearchModificationAnalysis::FuzzyDoubleComparator(1e-9));
    for (double mass : cluster_masses)
    {
      histogram[mass] = 10.0;
      charges[mass] = 1;
    }
    return analyzer.mapDeltaMassesToModifications(histogram, charges, ids, tolerance, ppm);
  };

  auto names = [](const std::vector<Summary>& summaries)
  {
    std::set<std::string> result;
    for (const Summary& summary : summaries) result.insert(summary.name);
    return result;
  };

  // A cluster at the mass of oxidation is named oxidation, not the lowest-mass entry within the tolerance
  // (and not "thiocarboxylic acid//Oxidation or Hydroxylation": the second is 0.0178 Da away, not equally close)
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({15.9949}, ids);
    TEST_EQUAL(summaries.size(), 1)
    TEST_STRING_EQUAL(summaries[0].name, oxidation)
    TEST_EQUAL(summaries[0].count, 10)
    TEST_EQUAL(summaries[0].masses.size(), 1)
    TEST_REAL_SIMILAR(summaries[0].masses[0], oxidation_mass)
  }

  // The same holds when the nearest entry is the lowest one, with a higher one in the window (this was "A//B" before)
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({15.9773}, ids);
    TEST_EQUAL(summaries.size(), 1)
    TEST_STRING_EQUAL(summaries[0].name, thiocarboxylic)
    TEST_REAL_SIMILAR(summaries[0].masses[0], thiocarboxylic_mass)
  }

  // Candidates that are equally close are reported together, lower mass first. The distances may differ by up to
  // the width of a delta mass histogram bin (0.0005 Da), which a cluster mass cannot resolve
  {
    const double midpoint = 0.5 * (oxidation_mass + thiocarboxylic_mass);
    for (double offset : {0.0, 0.0001, -0.0001}) // distances differ by 0, 0.0002 and 0.0002 Da
    {
      PeptideIdentificationList ids;
      const auto summaries = map_clusters({midpoint + offset}, ids);
      TEST_EQUAL(summaries.size(), 1)
      TEST_STRING_EQUAL(summaries[0].name, thiocarboxylic + "//" + oxidation)
      TEST_REAL_SIMILAR(summaries[0].masses[0], midpoint + offset) // the cluster mass: none of the two is better
    }
    // distances differ by 0.002 Da: the nearer one wins on its own
    PeptideIdentificationList ids_up;
    const auto up = map_clusters({midpoint + 0.001}, ids_up);
    TEST_EQUAL(up.size(), 1)
    TEST_STRING_EQUAL(up[0].name, oxidation)
    PeptideIdentificationList ids_down;
    const auto down = map_clusters({midpoint - 0.001}, ids_down);
    TEST_EQUAL(down.size(), 1)
    TEST_STRING_EQUAL(down[0].name, thiocarboxylic)
  }

  // Only entries within the tolerance are candidates: with 0.005 Da, thiocarboxylic acid is out of reach
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({15.9949}, ids, 0.005);
    TEST_EQUAL(summaries.size(), 1)
    TEST_STRING_EQUAL(summaries[0].name, oxidation)
  }

  // A ppm tolerance maps with the cap on the tolerance (0.02 Da), as before
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({15.9949}, ids, 10.0, true);
    TEST_EQUAL(summaries.size(), 1)
    TEST_STRING_EQUAL(summaries[0].name, oxidation)
  }

  // The PSMs are annotated with the modification nearest to their delta mass. Both clusters have the other
  // modification within the tolerance, and the one of lower mass must not win for all PSMs
  {
    PeptideIdentificationList ids = make_ids({15.9949, 15.9949, 15.9949, 15.9773, 15.9773});
    const auto summaries = map_clusters({15.9949, 15.9773}, ids);
    TEST_EQUAL(summaries.size(), 2)
    TEST_EQUAL(names(summaries).count(oxidation), 1)
    TEST_EQUAL(names(summaries).count(thiocarboxylic), 1)
    const std::vector<std::string> expected = {oxidation, oxidation, oxidation, thiocarboxylic, thiocarboxylic};
    for (Size i = 0; i < ids.size(); ++i)
    {
      TEST_STRING_EQUAL(ids[i].getHits()[0].getMetaValue("PTM").toString(), expected[i])
    }
  }

  // A cluster 1 Da above a modification that was found is explained as that modification plus 1 Da: the nearest one,
  // not the one of lowest mass in the window around (cluster - 1 Da), which is thiocarboxylic acid here
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({15.9773, 15.9949, 16.9949}, ids);
    const auto found = names(summaries);
    TEST_EQUAL(found.size(), 3)
    TEST_EQUAL(found.count(oxidation), 1)
    TEST_EQUAL(found.count(thiocarboxylic), 1)
    TEST_EQUAL(found.count(oxidation + "+1Da"), 1)
    TEST_EQUAL(found.count(thiocarboxylic + "+1Da"), 0)
  }

  // A cluster without a modification in the window is matched to the combination of two modifications whose summed mass
  // is nearest, also if that is lower than the cluster mass (there are none at about -57 Da, so it is a combination)
  {
    std::vector<std::string> mod_names;
    mod_db->getAllSearchModifications(mod_names);
    std::set<double> single_masses;
    for (const std::string& mod_name : mod_names)
    {
      const ResidueModification* mod = mod_db->getModification(mod_name);
      if (!mod->getFullName().contains("substitution")) single_masses.insert(mod->getDiffMonoMass());
    }
    const std::vector<double> masses(single_masses.begin(), single_masses.end());

    const double cluster = -57.0095;
    for (double single : masses) TEST_EQUAL(std::abs(single - cluster) > 0.02, true)

    double nearest_distance = std::numeric_limits<double>::max();
    double nearest_sum = 0.0;
    for (Size i = 0; i < masses.size(); ++i)
    {
      for (Size j = i; j < masses.size(); ++j)
      {
        const double sum = masses[i] + masses[j];
        if (std::abs(sum - cluster) < nearest_distance)
        {
          nearest_distance = std::abs(sum - cluster);
          nearest_sum = sum;
        }
      }
    }
    TEST_EQUAL(nearest_distance < 0.02, true)
    TEST_EQUAL(nearest_sum < cluster, true) // the nearest combination is below the cluster mass, not the first one above it

    PeptideIdentificationList ids;
    const auto summaries = map_clusters({cluster}, ids);
    TEST_EQUAL(summaries.size(), 1)
    TEST_EQUAL(summaries[0].name.find("++") != std::string::npos, true)
    TEST_EQUAL(summaries[0].masses.size(), 1)
    TEST_REAL_SIMILAR(summaries[0].masses[0], nearest_sum)
  }

  // No modification, no combination within the tolerance
  {
    PeptideIdentificationList ids;
    const auto summaries = map_clusters({123.456}, ids, 0.0005);
    TEST_EQUAL(summaries.size(), 1)
    TEST_STRING_EQUAL(summaries[0].name, "Unknown123.0")
  }
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

  // The mapped modification is the nearest one within the tolerance, as in mapDeltaMassesToModifications():
  // 15.9949 Da is oxidation, not thiocarboxylic acid (+15.9772 Da, the lowest mass within 0.02 Da)
  const ModificationsDB* mod_db = ModificationsDB::getInstance();
  const std::string oxidation = mod_db->getModification("Oxidation (M)")->getFullName();
  const double oxidation_mass = 15.994915;
  const double thiocarboxylic_mass = 15.977156;
  const double midpoint = 0.5 * (oxidation_mass + thiocarboxylic_mass);
  histogram.clear();
  charge_counts.clear();
  for (double mass : {15.9949, 15.9773, midpoint})
  {
    histogram[mass] = 10.0;
    charge_counts[mass] = 1;
  }
  stats = analyzer.generateDeltaMassStatistics(histogram, charge_counts, peptide_ids, 0.02, false);
  TEST_EQUAL(stats.entries.size(), 3)
  for (const auto& entry : stats.entries)
  {
    TEST_EQUAL(entry.is_known_modification, true)
    if (std::abs(entry.delta_mass - 15.9949) < 1e-6)
    {
      TEST_STRING_EQUAL(entry.mapped_modification, oxidation)
    }
    else if (std::abs(entry.delta_mass - 15.9773) < 1e-6)
    {
      TEST_STRING_EQUAL(entry.mapped_modification, "thiocarboxylic acid")
    }
    else // equally close
    {
      TEST_STRING_EQUAL(entry.mapped_modification, "thiocarboxylic acid//" + oxidation)
    }
  }

  // The counts are of the target PSMs, and the unique peptides of an entry are those of its histogram bin (not of all
  // PSMs within the precursor tolerance, which in an open search is the whole window): an oxidation bin with two PSMs of
  // one peptide, an acetyl bin with one PSM each of two peptides, an unmodified PSM and a decoy PSM in the oxidation bin.
  peptide_ids.clear();
  auto add = [&peptide_ids](const std::string& sequence, double delta_mass, bool decoy)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sequence));
    hit.setCharge(2);
    hit.setMetaValue("DeltaMass", delta_mass);
    hit.setMetaValue("target_decoy", decoy ? "decoy" : "target");
    PeptideIdentification id;
    id.setHits({hit});
    peptide_ids.push_back(id);
  };
  add("PEPTIDEK", 15.9951, false);
  add("PEPTIDEK", 15.9949, false);
  add("PEPTIDER", 42.0106, false);
  add("ELVISK", 42.0104, false);
  add("SAMPLER", 0.001, false);
  add("KEDITPEP", 15.9950, true);
  auto [patterns, charges] = analyzer.analyzeDeltaMassPatterns(peptide_ids);
  stats = analyzer.generateDeltaMassStatistics(patterns, charges, peptide_ids, 500.0, false); // an open search window
  TEST_EQUAL(stats.total_psms, 5)
  TEST_EQUAL(stats.modified_psms, 4)
  TEST_EQUAL(stats.unmodified_psms, 1)
  TEST_EQUAL(stats.decoy_psms, 1)
  TEST_EQUAL(stats.entries.size(), 2)
  for (const auto& entry : stats.entries)
  {
    if (std::abs(entry.delta_mass - 15.995) < 1e-6)
    {
      TEST_EQUAL(entry.count, 2)
      TEST_EQUAL(entry.unique_peptides, 1)
    }
    else
    {
      TEST_REAL_SIMILAR(entry.delta_mass, 42.0105)
      TEST_EQUAL(entry.count, 2)
      TEST_EQUAL(entry.unique_peptides, 2)
    }
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
