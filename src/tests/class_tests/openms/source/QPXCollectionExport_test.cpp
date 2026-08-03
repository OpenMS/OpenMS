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
#include <OpenMS/FORMAT/QPXCollectionExport.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ConsensusMapArrowExport.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  /// A label-free map over @p paths, one fraction each, all in one sample, with one quantified
  /// protein group. This is the shape every section below perturbs.
  ConsensusMap makeMap(const std::vector<std::string>& paths, Size n_samples = 1)
  {
    ConsensusMap cmap;
    cmap.setExperimentType("label-free");
    for (Size m = 0; m < paths.size(); ++m)
    {
      ConsensusMap::ColumnHeader ch;
      ch.filename = paths[m];
      ch.label = "label-free";
      ch.setMetaValue("fraction", static_cast<int>(m) + 1);
      ch.setMetaValue("fraction_group", 1);
      ch.setMetaValue("sample_name", "S1");
      cmap.getColumnHeaders()[m] = ch;
    }

    ProteinIdentification prot;
    prot.setIdentifier("merged");
    prot.setScoreType("q-value");
    prot.setHigherScoreBetter(false);
    prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
    ProteinHit ph; ph.setAccession("PROT_A");
    prot.setHits({ph});
    ProteinIdentification::ProteinGroup grp;
    grp.accessions = {"PROT_A"};
    grp.probability = 0.01;
    grp.getFloatDataArrays().resize(1);
    grp.getFloatDataArrays()[0].setName("abundances");
    for (Size i = 0; i < n_samples; ++i) { grp.getFloatDataArrays()[0].push_back(100.0f * (i + 1)); }
    grp.getFloatDataArrays().resize(2);
    grp.getFloatDataArrays()[1].setName("fraction_group_level_abundance");
    grp.getFloatDataArrays()[1].push_back(100.0f);
    grp.getIntegerDataArrays().resize(2);
    grp.getIntegerDataArrays()[0].setName("fraction_group_level_fraction_group");
    grp.getIntegerDataArrays()[0].push_back(1);
    grp.getIntegerDataArrays()[1].setName("fraction_group_level_label");
    grp.getIntegerDataArrays()[1].push_back(1);
    prot.insertIndistinguishableProteins(grp);
    cmap.setProteinIdentifications({prot});
    return cmap;
  }

  /// One PSM naming file @p merge_index of the merged run; merge_index < 0 omits the metavalue.
  PeptideIdentification makePSM(int merge_index, const std::string& seq = "PEPTIDEK")
  {
    PeptideIdentification pid;
    pid.setIdentifier("merged");
    pid.setScoreType("q-value");
    pid.setHigherScoreBetter(false);
    if (merge_index >= 0) { pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, merge_index); }
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setCharge(2);
    hit.setScore(0.01);
    PeptideEvidence ev; ev.setProteinAccession("PROT_A");
    hit.setPeptideEvidences({ev});
    pid.setHits({hit});
    return pid;
  }

  void addFeature(ConsensusMap& cmap, Size map_index, const std::vector<PeptideIdentification>& pids)
  {
    ConsensusFeature cf;
    cf.setMZ(500.0); cf.setRT(100.0 + map_index); cf.setCharge(2);
    BaseFeature bf; bf.setIntensity(10.0f); bf.setMZ(500.0); bf.setRT(100.0 + map_index); bf.setCharge(2);
    cf.insert(map_index, bf);
    cf.setPeptideIdentifications(pids);
    cmap.push_back(cf);
  }

  /// A one-sample design over @p paths, optionally omitting one (path, label) row.
  ExperimentalDesign makeDesign(const std::vector<std::string>& paths, int omit_row = -1)
  {
    ExperimentalDesign::MSFileSection section;
    for (Size m = 0; m < paths.size(); ++m)
    {
      if (static_cast<int>(m) == omit_row) { continue; }
      ExperimentalDesign::MSFileSectionEntry e;
      e.path = paths[m];
      e.fraction = static_cast<unsigned>(m) + 1;
      e.fraction_group = 1;
      e.label = 1;
      e.sample = 0;
      e.sample_name = "S1";
      section.push_back(e);
    }
    ExperimentalDesign::SampleSection samples;
    samples.addSample("S1");
    ExperimentalDesign d;
    d.setMSFileSection(section);
    d.setSampleSection(samples);
    return d;
  }
}

START_TEST(QPXCollectionExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static bool requireExportable(const ConsensusMap&, const ExperimentalDesign&)))
{
  // A collection all three views can write passes, and each view individually agrees.
  const std::vector<std::string> paths = {"/d/F1.mzML", "/d/F2.mzML"};
  ConsensusMap cmap = makeMap(paths);
  addFeature(cmap, 0, {makePSM(0)});
  addFeature(cmap, 1, {makePSM(1)});
  const ExperimentalDesign design = makeDesign(paths);

  TEST_TRUE(QPXCollectionExport::requireExportable(cmap, design))
  TEST_NOT_EQUAL(ConsensusMapArrowExport::exportToArrow(cmap), nullptr)
  TEST_NOT_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, design), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses missing protein-identification prerequisites))
{
  const std::vector<std::string> paths = {"/d/F1.mzML"};

  ConsensusMap no_ids = makeMap(paths);
  no_ids.setProteinIdentifications({});
  TEST_FALSE(QPXCollectionExport::requireExportable(no_ids, makeDesign(paths)))

  ConsensusMap no_groups = makeMap(paths);
  no_groups.getProteinIdentifications()[0].getIndistinguishableProteins().clear();
  TEST_FALSE(QPXCollectionExport::requireExportable(no_groups, makeDesign(paths)))
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - rejects the same duplicate run identifier as the feature writer))
{
  const std::vector<std::string> paths = {"/d/F1.mzML"};
  ConsensusMap cmap = makeMap(paths);
  auto duplicate = cmap.getProteinIdentifications()[0];
  duplicate.setPrimaryMSRunPath({"/d/other.mzML"});
  cmap.getProteinIdentifications().push_back(duplicate);

  TEST_EXCEPTION(Exception::InvalidValue,
                 QPXCollectionExport::requireExportable(cmap, makeDesign(paths)))
  TEST_EXCEPTION(Exception::InvalidValue, ConsensusMapArrowExport::exportToArrow(cmap))
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses what the feature view refuses))
{
  // Divergent identifications: the feature view throws, so the collection must throw too --
  // otherwise the throw would land after quantms.feature.parquet is already on disk.
  const std::vector<std::string> paths = {"/d/F1.mzML"};
  ConsensusMap cmap = makeMap(paths);
  addFeature(cmap, 0, {makePSM(0, "PEPTIDEK"), makePSM(0, "OTHERPEPTIDEK")});
  const ExperimentalDesign design = makeDesign(paths);

  TEST_EXCEPTION(Exception::IllegalArgument, QPXCollectionExport::requireExportable(cmap, design))
  TEST_EXCEPTION(Exception::IllegalArgument, ConsensusMapArrowExport::exportToArrow(cmap))
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses what the psm view refuses))
{
  // A PSM of a merged run with no id_merge_index. The psm view refuses it, but only after the
  // feature file has been written -- and the feature view does NOT catch this one, because it
  // validates only each feature's winning identification and this PSM is unassigned.
  const std::vector<std::string> paths = {"/d/F1.mzML", "/d/F2.mzML"};
  ConsensusMap cmap = makeMap(paths);
  addFeature(cmap, 0, {makePSM(0)});
  cmap.getUnassignedPeptideIdentifications() = {makePSM(-1, "STRAYPEPTIDEK")};
  const ExperimentalDesign design = makeDesign(paths);

  // The feature view accepts it -- which is exactly why the collection-level check is needed.
  TEST_NOT_EQUAL(ConsensusMapArrowExport::exportToArrow(cmap), nullptr)
  TEST_EXCEPTION(Exception::MissingInformation, QPXCollectionExport::requireExportable(cmap, design))
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses an out-of-range design sample))
{
  const std::vector<std::string> paths = {"/d/F1.mzML"};
  ConsensusMap cmap = makeMap(paths);
  ExperimentalDesign design = makeDesign(paths);
  auto section = design.getMSFileSection();
  section[0].sample = 1; // sample section contains only index 0
  design.setMSFileSection(section);

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, design))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, design), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses incomplete unit-level annotation arrays))
{
  const std::vector<std::string> paths = {"/d/F1.mzML"};
  ConsensusMap cmap = makeMap(paths);
  cmap.getProteinIdentifications()[0].getIndistinguishableProteins()[0]
    .getIntegerDataArrays().pop_back();

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, makeDesign(paths)))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, makeDesign(paths)), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses conflicting label names within a fraction group))
{
  const std::vector<std::string> paths = {"/d/F1.mzML", "/d/F2.mzML"};
  ConsensusMap cmap = makeMap(paths);
  cmap.setExperimentType("labeled_MS2");
  auto& first = cmap.getColumnHeaders()[0];
  first.label = "tmt6plex_126";
  first.setMetaValue("channel_name", "126");
  first.setMetaValue("channel_id", 0);
  auto& second = cmap.getColumnHeaders()[1];
  second.label = "itraq4plex_114";
  second.setMetaValue("channel_name", "114");
  second.setMetaValue("channel_id", 0);

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, makeDesign(paths)))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, makeDesign(paths)), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses two channels with the same canonical label))
{
  const std::vector<std::string> paths = {"/d/plex.mzML"};
  ConsensusMap cmap = makeMap(paths, 2);
  cmap.setExperimentType("labeled_MS2");

  auto& first = cmap.getColumnHeaders()[0];
  first.label = "tmt6plex_126";
  first.setMetaValue("channel_name", "126");
  first.setMetaValue("channel_id", 0);
  ConsensusMap::ColumnHeader second = first;
  second.setMetaValue("channel_id", 1);
  cmap.getColumnHeaders()[1] = second;

  auto& group = cmap.getProteinIdentifications()[0].getIndistinguishableProteins()[0];
  group.getFloatDataArrays()[1].push_back(200.0f);
  group.getIntegerDataArrays()[0].push_back(1);
  group.getIntegerDataArrays()[1].push_back(2);

  ExperimentalDesign::MSFileSection section;
  for (UInt label = 1; label <= 2; ++label)
  {
    ExperimentalDesign::MSFileSectionEntry entry;
    entry.path = paths[0];
    entry.fraction = 1;
    entry.fraction_group = 1;
    entry.label = label;
    entry.sample = label - 1;
    entry.sample_name = "S" + StringUtils::toStr(label);
    section.push_back(entry);
  }
  ExperimentalDesign::SampleSection samples;
  samples.addSample("S1");
  samples.addSample("S2");
  ExperimentalDesign design(section, samples);

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, design))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, design), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses a ragged design))
{
  // The pg view refuses a unit whose runs do not all carry every label it publishes -- the
  // shape of the shipped ExperimentalDesign_input_2_wrong.tsv, where one fraction has no
  // label-1 row. Raggedness needs at least two labels, so this is a 2-plex over 2 fractions.
  // Without the preflight that refusal lands after two files are written.
  const std::vector<std::string> paths = {"/d/F1.mzML", "/d/F2.mzML"};
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  const char* names[2] = {"126", "127N"};
  Size mi = 0;
  for (Size f = 0; f < paths.size(); ++f)
  {
    for (Size c = 0; c < 2; ++c)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = paths[f];
      h.label = std::string("tmt6plex_") + names[c];
      h.setMetaValue("channel_name", names[c]);
      h.setMetaValue("channel_id", static_cast<int>(c));
      cmap.getColumnHeaders()[mi++] = h;
    }
  }
  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath(StringList(paths.begin(), paths.end()));
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  grp.getFloatDataArrays().resize(1);
  grp.getFloatDataArrays()[0].setName("abundances");
  grp.getFloatDataArrays()[0].push_back(100.0f);
  grp.getFloatDataArrays()[0].push_back(200.0f);
  grp.getFloatDataArrays().resize(2);
  grp.getFloatDataArrays()[1].setName("fraction_group_level_abundance");
  grp.getFloatDataArrays()[1].push_back(100.0f);
  grp.getFloatDataArrays()[1].push_back(200.0f);
  grp.getIntegerDataArrays().resize(2);
  grp.getIntegerDataArrays()[0].setName("fraction_group_level_fraction_group");
  grp.getIntegerDataArrays()[0].push_back(1);
  grp.getIntegerDataArrays()[0].push_back(1);
  grp.getIntegerDataArrays()[1].setName("fraction_group_level_label");
  grp.getIntegerDataArrays()[1].push_back(1);
  grp.getIntegerDataArrays()[1].push_back(2);
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  ExperimentalDesign::MSFileSection section;
  for (Size f = 0; f < paths.size(); ++f)
  {
    for (unsigned lab = 1; lab <= 2; ++lab)
    {
      if (f == 1 && lab == 1) { continue; }            // the missing row
      ExperimentalDesign::MSFileSectionEntry e;
      e.path = paths[f];
      e.fraction = static_cast<unsigned>(f) + 1;
      e.fraction_group = 1;
      e.label = lab;
      e.sample = lab - 1;
      e.sample_name = "S" + StringUtils::toStr(lab);
      section.push_back(e);
    }
  }
  ExperimentalDesign::SampleSection samples;
  samples.addSample("S1");
  samples.addSample("S2");
  ExperimentalDesign ragged;
  ragged.setMSFileSection(section);
  ragged.setSampleSection(samples);

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, ragged))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, ragged), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses a design that is not the one used to quantify))
{
  // Two samples in the design, one abundance in the group.
  const std::vector<std::string> paths = {"/d/F1.mzML"};
  ConsensusMap cmap = makeMap(paths, /*n_samples=*/1);
  addFeature(cmap, 0, {makePSM(0)});

  ExperimentalDesign::MSFileSection section;
  ExperimentalDesign::MSFileSectionEntry e;
  e.path = paths[0]; e.fraction = 1; e.fraction_group = 1; e.label = 1; e.sample = 0; e.sample_name = "S1";
  section.push_back(e);
  ExperimentalDesign::SampleSection samples;
  samples.addSample("S1");
  samples.addSample("S2");                       // a sample no MS file references
  ExperimentalDesign foreign;
  foreign.setMSFileSection(section);
  foreign.setSampleSection(samples);
  TEST_EQUAL(foreign.getNumberOfSamples(), 2)

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, foreign))
  TEST_EQUAL(ProteinGroupArrowExport::exportToArrow(cmap, foreign), nullptr)
}
END_SECTION

START_SECTION(([EXTRA] requireExportable - refuses a channel QPX cannot label))
{
  // intensities[].label is a documented join key. The feature view refuses an unnameable
  // channel; the pg view invents a token. The collection refuses, so neither can write.
  ConsensusMap cmap;
  cmap.setExperimentType("labeled_MS2");
  for (Size c = 0; c < 2; ++c)
  {
    ConsensusMap::ColumnHeader h;
    h.filename = "/d/run.mzML";
    h.label = std::string("not_a_method_") + (c == 0 ? "126" : "127N");
    h.setMetaValue("channel_name", c == 0 ? "126" : "127N");
    h.setMetaValue("channel_id", static_cast<int>(c));
    cmap.getColumnHeaders()[c] = h;
  }
  ProteinIdentification prot;
  prot.setIdentifier("merged");
  prot.setScoreType("q-value");
  prot.setPrimaryMSRunPath({"/d/run.mzML"});
  ProteinHit ph; ph.setAccession("PROT_A");
  prot.setHits({ph});
  ProteinIdentification::ProteinGroup grp;
  grp.accessions = {"PROT_A"};
  grp.probability = 0.01;
  prot.insertIndistinguishableProteins(grp);
  cmap.setProteinIdentifications({prot});

  TEST_FALSE(QPXCollectionExport::requireExportable(cmap, makeDesign({"/d/run.mzML"})))
  TEST_EQUAL(ConsensusMapArrowExport::exportToArrow(cmap), nullptr)
}
END_SECTION

END_TEST
