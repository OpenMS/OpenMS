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

END_TEST
