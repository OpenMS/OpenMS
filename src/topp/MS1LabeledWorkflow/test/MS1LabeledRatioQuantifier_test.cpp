// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledRatioQuantifier.h"

#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

namespace
{
void addPeptide(ConsensusMap& map, const std::string& sequence, Size run, const std::vector<double>& intensities)
{
  ConsensusFeature feature;
  for (Size channel = 0; channel < intensities.size(); ++channel)
  {
    FeatureHandle handle;
    handle.setMapIndex(run * 3 + channel);
    handle.setIntensity(intensities[channel]);
    feature.insert(handle);
  }
  PeptideHit hit;
  hit.setSequence(AASequence::fromString(sequence));
  PeptideEvidence evidence;
  evidence.setProteinAccession("P1");
  hit.addPeptideEvidence(evidence);
  PeptideIdentification id;
  id.insertHit(hit);
  feature.getPeptideIdentifications().push_back(id);
  map.push_back(feature);
}

// a feature of @p sequence (protein @p accession) with the channel intensities of every run (three channels per run),
// measured by targeted quantification in @p targeted_runs
void addFeature(ConsensusMap& map, const std::string& sequence, const std::string& accession,
                const std::map<Size, std::vector<double>>& intensities, const StringList& targeted_runs = {})
{
  ConsensusFeature feature;
  for (const auto& [run, channels] : intensities)
  {
    for (Size channel = 0; channel < channels.size(); ++channel)
    {
      FeatureHandle handle;
      handle.setMapIndex(run * 3 + channel);
      handle.setIntensity(channels[channel]);
      feature.insert(handle);
    }
  }
  if (!targeted_runs.empty()) { feature.setMetaValue(MS1LabeledRatioQuantifier::TARGETED_RUNS, targeted_runs); }
  PeptideHit hit;
  hit.setSequence(AASequence::fromString(sequence));
  PeptideEvidence evidence;
  evidence.setProteinAccession(accession);
  hit.addPeptideEvidence(evidence);
  PeptideIdentification id;
  id.insertHit(hit);
  feature.getPeptideIdentifications().push_back(id);
  map.push_back(feature);
}
} // namespace

START_TEST(MS1LabeledRatioQuantifier, "$Id$")

START_SECTION((void run(ConsensusMap&, const ExperimentalDesign&, ProteinIdentification&)))
{
  ConsensusMap map;
  map.setExperimentType("labeled_MS1");
  ExperimentalDesign design;
  ExperimentalDesign::MSFileSection files;
  for (unsigned run = 0; run < 2; ++run)
  {
    for (unsigned channel = 0; channel < 3; ++channel)
    {
      auto& header = map.getColumnHeaders()[run * 3 + channel];
      header.filename = "run" + std::to_string(run) + ".mzML";
      header.setMetaValue("channel_id", channel);
      ExperimentalDesign::MSFileSectionEntry entry;
      entry.path = header.filename;
      entry.fraction_group = run + 1;
      entry.label = channel + 1;
      entry.sample = run * 3 + channel;
      files.push_back(entry);
    }
  }
  design.setMSFileSection(files);
  ProteinIdentification proteins;
  ProteinIdentification::ProteinGroup group;
  group.accessions = {"P1"};
  proteins.getIndistinguishableProteins() = {group};

  // Disjoint partial triplexes: neither M/L nor H/L has the required two peptides.
  addPeptide(map, "PEPTIDEK", 0, {100, 200, 0});
  addPeptide(map, "AAAAK", 0, {100, 0, 300});
  MS1LabeledRatioQuantifier quantifier;
  quantifier.run(map, design, proteins);
  TEST_EQUAL(quantifier.getProteinGroupRatios().size(), 0)
  TEST_EQUAL(proteins.getIndistinguishableProteins()[0].getFloatDataArrays().size(), 0)

  // A passing comparison in fraction group 2 must not admit a reference-only row in group 1.
  addPeptide(map, "PEPTIDEK", 1, {100, 200, 0});
  addPeptide(map, "AAAAK", 1, {100, 400, 0});
  quantifier.run(map, design, proteins);
  const auto& ratios = quantifier.getProteinGroupRatios().at("P1");
  TEST_EQUAL(ratios.size(), 2)
  TEST_EQUAL(ratios[0].fraction_group, 2)
  TEST_EQUAL(ratios[0].channel, 1)
  TEST_EQUAL(ratios[0].count, 2)
  TEST_REAL_SIMILAR(ratios[0].ratio, 1.0)
  TEST_EQUAL(ratios[1].fraction_group, 2)
  TEST_EQUAL(ratios[1].channel, 2)
  TEST_EQUAL(ratios[1].count, 2)
  TEST_REAL_SIMILAR(ratios[1].ratio, 3.0)

  // Channel 2 can be the reference; the same rule applies independently of channel ordering.
  Param param = quantifier.getParameters();
  param.setValue("reference_channel", 2);
  param.setValue("normalize", "false");
  quantifier.setParameters(param);
  ConsensusMap alternate = map;
  alternate.resize(0);
  addPeptide(alternate, "PEPTIDEK", 0, {100, 200, 0});
  addPeptide(alternate, "AAAAK", 0, {0, 200, 300});
  quantifier.run(alternate, design, proteins);
  TEST_EQUAL(quantifier.getProteinGroupRatios().size(), 0)
  // Re-running also clears earlier annotations from the protein group.
  TEST_EQUAL(proteins.getIndistinguishableProteins()[0].getFloatDataArrays().size(), 0)

  addPeptide(alternate, "VVVVK", 0, {100, 400, 0});
  quantifier.run(alternate, design, proteins);
  const auto& alternate_ratios = quantifier.getProteinGroupRatios().at("P1");
  TEST_EQUAL(alternate_ratios.size(), 2)
  TEST_EQUAL(alternate_ratios[0].channel, 1)
  TEST_EQUAL(alternate_ratios[0].count, 2)
  TEST_REAL_SIMILAR(alternate_ratios[0].ratio, 0.375)
  TEST_EQUAL(alternate_ratios[1].channel, 2)
  TEST_EQUAL(alternate_ratios[1].count, 3)
  TEST_REAL_SIMILAR(alternate_ratios[1].ratio, 1.0)
  TEST_EQUAL(proteins.getIndistinguishableProteins()[0].getFloatDataArrays().size(), 1)
}
END_SECTION

START_SECTION(([EXTRA] void run(ConsensusMap&, const ExperimentalDesign&, ProteinIdentification&) with evidences of targeted quantification))
{
  // two fractions (runs) of one fraction group, three channels each
  ConsensusMap map;
  map.setExperimentType("labeled_MS1");
  ExperimentalDesign design;
  ExperimentalDesign::MSFileSection files;
  for (unsigned run = 0; run < 2; ++run)
  {
    for (unsigned channel = 0; channel < 3; ++channel)
    {
      auto& header = map.getColumnHeaders()[run * 3 + channel];
      header.filename = "/data/run" + std::to_string(run) + ".mzML";
      header.setMetaValue("channel_id", channel);
      ExperimentalDesign::MSFileSectionEntry entry;
      entry.path = header.filename;
      entry.fraction_group = 1;
      entry.fraction = run + 1;
      entry.label = channel + 1;
      entry.sample = channel;
      files.push_back(entry);
    }
  }
  design.setMSFileSection(files);
  ProteinIdentification proteins;
  ProteinIdentification::ProteinGroup p1, p2;
  p1.accessions = {"P1"};
  p2.accessions = {"P2"};
  proteins.getIndistinguishableProteins() = {p1, p2};

  // P1: two detected peptides (M/L 2 and 4) and a targeted one (0.5). PEPAK's targeted evidence in run 1 sits on the same
  // feature as its detected one in run 0.
  addFeature(map, "PEPAK", "P1", {{0, {100, 200, 0}}, {1, {100, 100, 0}}}, {"run1.mzML"});
  addFeature(map, "PEPBK", "P1", {{0, {100, 400, 0}}});
  addFeature(map, "PEPCK", "P1", {{0, {100, 50, 0}}}, {"run0.mzML"});
  // P2: one detected peptide (4), one targeted (0.25)
  addFeature(map, "PEPDK", "P2", {{0, {100, 400, 0}}});
  addFeature(map, "PEPEK", "P2", {{1, {100, 25, 0}}}, {"run1.mzML"});

  MS1LabeledRatioQuantifier quantifier; // min_ratio_count 2, normalized
  quantifier.run(map, design, proteins);

  // peptide ratios: the targeted evidence only for a peptide without another one
  const auto& peptides = quantifier.getPeptideRatios();
  const auto& pep_a = peptides.at(AASequence::fromString("PEPAK"));
  TEST_EQUAL(pep_a.size(), 2) // reference and channel 2
  TEST_EQUAL(pep_a[1].channel, 2)
  TEST_EQUAL(pep_a[1].count, 1)
  TEST_EQUAL(pep_a[1].targeted_count, 0)
  TEST_REAL_SIMILAR(pep_a[1].ratio, 2.0) // not the median of 2 and 1
  const auto& pep_c = peptides.at(AASequence::fromString("PEPCK"));
  TEST_EQUAL(pep_c[1].count, 1)
  TEST_EQUAL(pep_c[1].targeted_count, 1)
  TEST_REAL_SIMILAR(pep_c[1].ratio, 0.5)
  // normalized by the median of the detected peptide ratios (2, 4, 4), not of all five (0.25, 0.5, 2, 4, 4)
  TEST_REAL_SIMILAR(pep_c[1].normalized_ratio, 0.125)

  // P1 has two peptide ratios without targeted quantification: they alone give its ratio
  const auto& r1 = quantifier.getProteinGroupRatios().at("P1");
  TEST_EQUAL(r1.size(), 2)
  TEST_EQUAL(r1[0].channel, 1)
  TEST_EQUAL(r1[0].count, 2)
  TEST_EQUAL(r1[0].targeted_count, 0)
  TEST_EQUAL(r1[1].channel, 2)
  TEST_EQUAL(r1[1].count, 2)
  TEST_EQUAL(r1[1].targeted_count, 0)
  TEST_REAL_SIMILAR(r1[1].ratio, 3.0)
  TEST_REAL_SIMILAR(r1[1].normalized_ratio, 0.75)

  // P2 has one: the targeted one is added, which reaches the minimum of two
  const auto& r2 = quantifier.getProteinGroupRatios().at("P2");
  TEST_EQUAL(r2.size(), 2)
  TEST_EQUAL(r2[0].channel, 1)
  TEST_EQUAL(r2[0].count, 2)
  TEST_EQUAL(r2[0].targeted_count, 1)
  TEST_EQUAL(r2[1].channel, 2)
  TEST_EQUAL(r2[1].count, 2)
  TEST_EQUAL(r2[1].targeted_count, 1)
  TEST_REAL_SIMILAR(r2[1].ratio, 2.125)
  TEST_REAL_SIMILAR(r2[1].normalized_ratio, 0.53125)

  // with a minimum of one, P2's detected peptide suffices, and its ratio is that peptide's
  Param param = quantifier.getParameters();
  param.setValue("min_ratio_count", 1);
  quantifier.setParameters(param);
  quantifier.run(map, design, proteins);
  const auto& r2_min1 = quantifier.getProteinGroupRatios().at("P2");
  TEST_EQUAL(r2_min1[1].count, 1)
  TEST_EQUAL(r2_min1[1].targeted_count, 0)
  TEST_REAL_SIMILAR(r2_min1[1].ratio, 4.0)

  // without the targeted features, P1's and P2's detected ratios are the same as above
  ConsensusMap detected_only = map;
  detected_only.clear(false);
  for (const ConsensusFeature& feature : map)
  {
    if (!feature.metaValueExists(MS1LabeledRatioQuantifier::TARGETED_RUNS)) { detected_only.push_back(feature); }
  }
  addFeature(detected_only, "PEPAK", "P1", {{0, {100, 200, 0}}});
  quantifier.run(detected_only, design, proteins);
  TEST_REAL_SIMILAR(quantifier.getProteinGroupRatios().at("P1")[1].ratio, 3.0)
  TEST_REAL_SIMILAR(quantifier.getProteinGroupRatios().at("P1")[1].normalized_ratio, 0.75)
  TEST_REAL_SIMILAR(quantifier.getProteinGroupRatios().at("P2")[1].ratio, 4.0)
}
END_SECTION

END_TEST
