// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/test_config.h>
#include <filesystem>
#include <fstream>
using namespace OpenMS;
using ID = IdentificationData;
START_TEST(IdentificationDataMigration, "$Id$")
START_SECTION((released OMS schemas load into owning records and roundtrip natively))
{
  const auto root = std::filesystem::path(OPENMS_GET_TEST_DATA_PATH("dummy")).parent_path().parent_path().parent_path().parent_path();
  for (const auto& name : {"JSONExporter_RNA.oms", "IDMerger_6_input1.oms", "IDMerger_6_input2.oms", "MapAlignerIdentification_8_input1.oms"})
  {
    const auto input = (root / "topp" / name).string();
    ID data;
    OMSFile().load(input, data);
    data.validate();
    TEST_EQUAL(data.empty(), false);
    TEST_EQUAL(data.getRuns().front().getNumberOfMatches() > 0, true);
    std::string output;
    NEW_TMP_FILE(output);
    OMSFile().store(output, data);
    ID restored;
    OMSFile().load(output, restored);
    TEST_EQUAL(restored == data, true);
  }
}
END_SECTION
START_SECTION((incompatible legacy PSM layouts fail while compatible inference retains both grouping types))
{
  const auto root = std::filesystem::path(OPENMS_GET_TEST_DATA_PATH("dummy")).parent_path().parent_path().parent_path().parent_path();
  ID incompatible;
  TEST_EXCEPTION(Exception::InvalidValue, OMSFile().load((root / "topp" / "JSONExporter_protein.oms").string(), incompatible));
  TEST_EQUAL(incompatible.empty(), true);
  // Derived from the released protein OMS fixture: reversed PSM scores are negated
  // and assigned its MOWSE definition, giving one explicit common orientation.
  ID data;
  OMSFile().load(OPENMS_GET_TEST_DATA_PATH("IdentificationDataMigration_protein.oms"), data);
  TEST_EQUAL(data.getRuns().front().getNumberOfMatches(), 5);
  TEST_EQUAL(data.getInferenceResults().size(), 1);
  const auto& proteins = data.getInferenceResults()[0].proteins;
  TEST_EQUAL(proteins.getProteinGroups().size(), 1);
  TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 1);
  TEST_REAL_SIMILAR(proteins.getProteinGroups()[0].probability, 0.88);
  std::string output;
  NEW_TMP_FILE(output);
  OMSFile().store(output, data);
  ID copy;
  OMSFile().load(output, copy);
  TEST_EQUAL(data == copy, true);
}
END_SECTION
START_SECTION((explicit sequence catalogs roundtrip without weakening the PSM score contract))
{
  ID data;
  auto& run = data.addRun("digestion", ID::MoleculeKind::OLIGONUCLEOTIDE);
  ProteinIdentification processing;
  processing.setMetaValue("identification:catalog", "true");
  run.setProcessingMetadata(processing);
  ID::Observation observation;
  observation.data_id = "digest=1";
  const auto query = run.addIdentification(run.addSource({}), observation);
  ID::MatchData match;
  match.encoding = ID::Encoding::NA_SEQUENCE;
  match.representation = "ACUGp";
  run.addMatch(query, match);
  std::string output;
  NEW_TMP_FILE(output);
  OMSFile().store(output, data);
  ID copy;
  OMSFile().load(output, copy);
  TEST_EQUAL(data == copy, true);
  processing.removeMetaValue("identification:catalog");
  copy.getRun("digestion").setProcessingMetadata(processing);
  TEST_EXCEPTION(Exception::InvalidValue, copy.validate());
  ID::ScoreDefinition score;
  score.name = "not a catalog score";
  run.addScore(score);
  TEST_EXCEPTION(Exception::InvalidValue, data.validate());
}
END_SECTION
START_SECTION((failed OMS writes and loads preserve existing data))
{
  const auto root = std::filesystem::path(OPENMS_GET_TEST_DATA_PATH("dummy")).parent_path().parent_path().parent_path().parent_path();
  ID data;
  OMSFile().load((root / "topp" / "JSONExporter_RNA.oms").string(), data);
  std::string output;
  NEW_TMP_FILE(output);
  OMSFile().store(output, data);
  OMSFile().store(output, data); // replacing an existing file remains supported
  const auto expected = data;
  const auto& original = data.getRuns().front();
  ID::Run second("second", original.getMoleculeKind());
  for (const auto& definition : original.getScoreDefinitions())
    second.addScore(definition);
  second.setPrimaryScore(second.getScoreId(original.getPrimaryScore()->value));
  data.addRun(std::move(second));
  ID::ScoreDefinition extra;
  extra.name = "inconsistent";
  data.getRun("second").addScore(extra);
  TEST_EXCEPTION(Exception::InvalidValue, OMSFile().store(output, data));
  ID restored;
  OMSFile().load(output, restored);
  TEST_EQUAL(restored == expected, true);
  std::string broken;
  NEW_TMP_FILE(broken);
  {
    std::ofstream file(broken);
    file << "not SQLite";
  }
  TEST_EXCEPTION(Exception::FailedAPICall, OMSFile().load(broken, restored));
  TEST_EQUAL(restored == expected, true);
}
END_SECTION
START_SECTION((legacy accurate - mass features retain molecular and match associations))
{
  FeatureMap map;
  OMSFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabMFile_input_1.oms"), map);
  TEST_EQUAL(map.empty(), false);
  TEST_EQUAL(map.getIdentificationData().empty(), false);
  Size matches = 0;
  for (const auto& feature : map)
    for (const auto& reference : feature.getIDMatches())
    {
      const auto* run = map.getIdentificationData().findRunByUuid(reference.run_uuid);
      TEST_EQUAL(run != nullptr, true);
      TEST_EQUAL(run->findMatch(reference.match) != nullptr, true);
      ++matches;
    }
  TEST_EQUAL(matches > 0, true);
  auto tab = MzTabM::exportFeatureMapToMzTabM(map);
  TEST_EQUAL(tab.getMSmallMoleculeEvidenceSectionRows().empty(), false);
  std::string output;
  NEW_TMP_FILE(output);
  OMSFile().store(output, map);
  FeatureMap copy;
  OMSFile().load(output, copy);
  TEST_EQUAL(copy == map, true);
}
END_SECTION
START_SECTION((pooled FDR, inclusive filtering and decoy evidence removal use stable values))
{
  ID data;
  ID::ScoreDefinition raw;
  raw.name = "raw";
  for (const auto& name : {"A", "B"})
  {
    auto& run = data.addRun(name);
    run.setPrimaryScore(run.addScore(raw));
    auto source = run.addSource({});
    ID::ParentRecord target, decoy;
    target.identity = {"db", "target"};
    target.target_decoy = ID::TargetDecoy::TARGET;
    decoy.identity = {"db", "decoy"};
    decoy.target_decoy = ID::TargetDecoy::DECOY;
    run.setParents(std::vector {target, decoy});
    ID::Observation observation;
    observation.data_id = "scan=1";
    observation.rt = 10;
    const auto query = run.addIdentification(source, observation);
    ID::MatchData match;
    match.representation = "PEPTIDE";
    match.target_decoy = ID::TargetDecoy::BOTH;
    match.parent_evidence = {{{"db", "target"}, {}, {}, "", ""}, {{"db", "decoy"}, {}, {}, "", ""}};
    run.addMatch(query, match, {100.0});
    match.target_decoy = ID::TargetDecoy::DECOY;
    match.parent_evidence.erase(match.parent_evidence.begin());
    run.addMatch(query, match, {50.0});
  }
  FalseDiscoveryRate fdr;
  auto parameters = fdr.getParameters();
  parameters.setValue("use_all_hits", "true");
  fdr.setParameters(parameters);
  const auto q = fdr.applyToObservationMatches(data, raw);
  data.validate();
  TEST_EQUAL(q.higher_better, false);
  TEST_EQUAL(data.getScoreDefinitions().size(), 2);
  IDFilter::removeDecoys(data);
  TEST_EQUAL(data.getRun("A").getNumberOfMatches(), 1);
  const auto& match = data.getRun("A").getSourceBlocks()[0].identifications[0].getMatches()[0];
  TEST_EQUAL(match.parent_evidence.size(), 1);
  TEST_EQUAL(match.parent_evidence[0].parent.accession, "target");
  IDFilter::filterObservationMatchesByScore(data, raw, 100);
  TEST_EQUAL(data.getRun("B").getNumberOfMatches(), 1);
  IDFilter::filterObservationMatchesByScore(data, raw, 101);
  TEST_EQUAL(data.getRun("B").getNumberOfMatches(), 0);
}
END_SECTION
END_TEST
