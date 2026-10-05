// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/test_config.h>
#include <algorithm>

using namespace OpenMS;
using ID = IdentificationData;
using Adapter = IdentificationDataAdapter;

namespace
{
ProteinIdentification protein()
{
  ProteinIdentification result;
  result.setIdentifier("search");
  result.setSearchEngine("test-search");
  result.setSearchEngineVersion("1");
  result.setScoreType("protein score");
  result.setPrimaryMSRunPath({"/exact/a.mzML", "/other/a.mzML"});
  auto parameters = result.getSearchParameters();
  parameters.db = "database.fasta";
  parameters.setMetaValue("empty list", StringList {});
  result.setSearchParameters(parameters);
  ProteinHit hit;
  hit.setAccession("P1");
  hit.setSequence("PEPTIDE");
  hit.setScore(0.95);
  hit.setRank(2);
  hit.setCoverage(50);
  hit.setDescription("description");
  hit.setMetaValue("user", "value");
  result.insertHit(hit);
  ProteinIdentification::ProteinGroup group;
  group.probability = 0.9;
  group.accessions = {"P1", "unlisted-alias"};
  result.getProteinGroups().push_back(group);
  return result;
}
PeptideIdentification peptide(const std::string& score_name = "PEP", double score = 0.1)
{
  PeptideIdentification result;
  result.setIdentifier("search");
  result.setScoreType(score_name);
  result.setHigherScoreBetter(false);
  result.setRT(42);
  result.setMZ(501);
  result.setSpectrumReference("controllerType=0 controllerNumber=1 scan=42");
  result.setMetaValue("empty", DataValue::EMPTY);
  PeptideHit hit(score, 2, 2, AASequence::fromString("PEPTIDE"));
  hit.setPeptideEvidences({PeptideEvidence("P1", 0, 6, '[', ']')});
  hit.setPeakAnnotations({{"y1", 1, 150.1, 500.0}});
  hit.setMetaValue("values", DoubleList {1.0, 2.0});
  hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
  result.insertHit(hit);
  return result;
}
} // namespace

START_TEST(IdentificationDataAdapter, "$Id$")

START_SECTION((static ImportResult importLegacy(const std::vector<ProteinIdentification>&, const PeptideIdentificationList&)))
{
  const auto original = protein();
  const auto first = peptide();
  auto second = peptide("q-value", 0.01);
  second.setMetaValue("id_merge_index", 1);
  auto empty = peptide();
  empty.getHits().clear();
  PeptideIdentificationList peptides {first, second, empty};
  auto imported = Adapter::importLegacy({original}, peptides);
  TEST_EQUAL(imported.data.getRuns().size(), 2)
  TEST_EQUAL(imported.queries.size(), 3)
  const auto& run = imported.data.getRuns()[0];
  TEST_EQUAL(run.getSourceBlocks()[0].source.path, "")
  TEST_EQUAL(run.getSourceBlocks()[0].source.primary_files.size(), 2)
  TEST_EQUAL(imported.data.getRuns()[1].getSourceBlocks()[0].source.path, "/other/a.mzML")
  TEST_EQUAL(run.getNumberOfIdentifications(), 2)
  TEST_EQUAL(run.getNumberOfMatches(), 1)
  TEST_EQUAL(run.getProcessingMetadata().getHits().size(), 0)
  TEST_EQUAL(run.getParents()->size(), 1)
  TEST_FALSE(imported.data.getInferenceResults()[0].inputs[0].membership_known)
  TEST_EQUAL(imported.data.getInferenceResults()[0].inputs.size(), 2)
  TEST_TRUE(imported.data.getInferenceResults()[0].proteins == original)
  const auto exported = Adapter::toLegacy(imported.data);
  TEST_TRUE(exported.losses.empty())
  TEST_EQUAL(exported.proteins.size(), 1)
  TEST_TRUE(exported.proteins[0] == original)
  // Source blocks preserve scientific order within a run; mappings preserve the
  // original cross-run order without requiring a dataset-wide sort.
  TEST_TRUE(exported.peptides[0] == first)
  TEST_TRUE(exported.peptides[1] == empty)
  TEST_TRUE(exported.peptides[2] == second)
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::fromLegacy({original, original}, peptides))
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::fromLegacy({}, peptides))
  auto malformed = first;
  malformed.setMetaValue("id_merge_index", 2);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::fromLegacy({original}, {malformed}))
}
END_SECTION

START_SECTION((static LegacyResult toLegacy(const IdentificationData&, const ExportOptions&)))
{
  auto data = Adapter::fromLegacy({protein()}, {peptide()});
  auto& run = data.getRun("search");
  const auto match_id = run.getSourceBlocks()[0].identifications[0].getMatches()[0].getId();
  auto second = data.getInferenceResults()[0];
  second.identifier = "second";
  data.addInferenceResult(second);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(data))
  Adapter::ExportOptions options;
  options.inference_result = "legacy:search";
  TEST_EQUAL(Adapter::toLegacy(data, options).proteins.size(), 1)
  // Filtering does not silently remove the selected protein inference output.
  run.eraseMatches([](const auto&) { return true; }, true);
  auto retained = Adapter::toLegacy(data, options);
  TEST_EQUAL(retained.peptides.size(), 1)
  TEST_EQUAL(retained.peptides[0].getHits().size(), 0)
  TEST_TRUE(retained.proteins[0] == protein())
  TEST_TRUE(run.findMatch(match_id) == nullptr)
  data.clearInferenceResults();
  options.inference_result.reset();
  auto reduced = Adapter::toLegacy(data, options);
  TEST_EQUAL(reduced.proteins[0].getProteinGroups().size(), 0)
  TEST_EQUAL(reduced.proteins[0].getHits().size(), 1)

  auto unsupported = Adapter::fromLegacy({protein()}, {peptide()});
  auto& edited = unsupported.getRun("search");
  const auto& query = edited.getSourceBlocks()[0].identifications[0];
  edited.setSelectedMatch(query.getId(), query.getMatches()[0].getId());
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(unsupported))
  options.loss_policy = Adapter::LossPolicy::ALLOW;
  TEST_EQUAL(Adapter::toLegacy(unsupported, options).losses.size(), 1)
}
END_SECTION

START_SECTION([EXTRA] idXML round trip preserves supported imported values)
{
  std::vector<ProteinIdentification> proteins;
  PeptideIdentificationList peptides;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), proteins, peptides);
  auto native = Adapter::fromLegacy(proteins, peptides);
  auto exported = Adapter::toLegacy(native);
  TEST_EQUAL(exported.proteins.size(), proteins.size())
  TEST_EQUAL(exported.peptides.size(), peptides.size())
  for (Size i = 0; i < proteins.size(); ++i)
    TEST_TRUE(exported.proteins[i] == proteins[i])
  for (const auto& item : peptides)
    TEST_TRUE(std::find(exported.peptides.begin(), exported.peptides.end(), item) != exported.peptides.end())
}
END_SECTION

START_SECTION([EXTRA] custom modifications retain definitions and reject conflicting chemistry)
{
  ResidueModification definition;
  definition.setId("OwningAdapterTest");
  definition.setOrigin('P');
  definition.setTermSpecificity(ResidueModification::ANYWHERE);
  definition.setFullId();
  definition.setDiffFormula(EmpiricalFormula("CH2"));
  definition.setDiffMonoMass(EmpiricalFormula("CH2").getMonoWeight());
  ModificationsDB::getInstance()->registerDefinition(definition);
  auto item = peptide();
  item.getHits()[0].setSequence(AASequence::fromString("P(OwningAdapterTest)EPTIDE"));
  auto data = Adapter::fromLegacy({protein()}, {item});
  auto& run = data.getRun("search");
  TEST_TRUE(run.getProcessingMetadata().getSearchParameters().metaValueExists(Constants::UserParam::MODIFICATION_DEFINITIONS))
  const auto exported = Adapter::toLegacy(data);
  TEST_TRUE(exported.peptides[0] == item)
  auto configuration = run.getProcessingMetadata();
  definition.setDiffMonoMass(definition.getDiffMonoMass() + 1.0);
  configuration.getSearchParameters().setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, definition.toDefinitionString());
  run.setProcessingMetadata(configuration);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(data))
}
END_SECTION

END_TEST
