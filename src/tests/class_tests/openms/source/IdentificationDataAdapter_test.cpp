// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
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
  TEST_EXCEPTION(Exception::InvalidValue, Adapter::importLegacy({original}, {first, peptide("q-value", 0.01)}))
  auto second = peptide("PEP", 0.01);
  second.setMetaValue("id_merge_index", 1);
  auto empty = peptide();
  empty.getHits().clear();
  PeptideIdentificationList peptides {first, second, empty};
  auto imported = Adapter::importLegacy({original}, peptides);
  TEST_EQUAL(imported.data.getRuns().size(), 1)
  TEST_EQUAL(imported.queries.size(), 3)
  const auto& run = imported.data.getRuns()[0];
  // The sources are the files of the legacy run, then an unknown file for the identifications without an index.
  ABORT_IF(run.getSources().size() != 3)
  TEST_EQUAL(run.getSources()[0].file.path, "/exact/a.mzML")
  TEST_EQUAL(run.getSources()[0].identifications.size(), 0)
  TEST_EQUAL(run.getSources()[1].file.path, "/other/a.mzML")
  TEST_EQUAL(run.getSources()[1].identifications.size(), 1)
  TEST_EQUAL(run.getSources()[1].identifications[0].metaValueExists("id_merge_index"), false)
  TEST_EQUAL(run.getSources()[2].file.path, "")
  TEST_EQUAL(run.getSources()[2].identifications.size(), 2)
  TEST_EQUAL(run.getSettings().metaValueExists("spectra_data"), false)
  TEST_EQUAL(run.getNumberOfIdentifications(), 3)
  TEST_EQUAL(run.getNumberOfMatches(), 2)
  // The settings keep the search engine; the proteins of the legacy run are an inference result.
  TEST_EQUAL(run.getSettings().software, "test-search")
  // The database of the legacy search is a database of the run; the search settings no longer name it.
  TEST_EQUAL(run.getSettings().search.db, "")
  ABORT_IF(run.getDatabases().size() != 1)
  TEST_EQUAL(run.getDatabases()[0].path, "database.fasta")
  TEST_EQUAL(run.getDatabaseSequences()->size(), 1)
  TEST_EQUAL(run.getDatabaseSequences()->at(0).accession, "P1")
  TEST_EQUAL(imported.data.getInferenceResults()[0].inputs[0].selection, "Imported legacy run-level provenance")
  TEST_EQUAL(imported.data.getInferenceResults()[0].inputs.size(), 1)
  auto without_files = original;
  without_files.removeMetaValue("spectra_data");
  TEST_TRUE(imported.data.getInferenceResults()[0].proteins == without_files)
  const auto exported = Adapter::toLegacy(imported.data);
  TEST_TRUE(exported.losses.empty())
  TEST_EQUAL(exported.proteins.size(), 1)
  TEST_TRUE(exported.proteins[0] == original)
  // Identifications are exported by source, in the order of the run's files; mappings preserve the
  // original cross-run order without requiring a dataset-wide sort.
  TEST_TRUE(exported.peptides[0] == second)
  TEST_TRUE(exported.peptides[1] == first)
  TEST_TRUE(exported.peptides[2] == empty)
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
  // The identification has no index into the two files, so it is in the unknown source after them.
  const auto match_id = run.getSources().back().identifications[0].getMatches()[0].getId();
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
  const auto& query = edited.getSources().back().identifications[0];
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
  TEST_EXCEPTION(Exception::InvalidValue, Adapter::fromLegacy(proteins, peptides))
  // The fixture deliberately mixes higher/lower-is-better for MOWSE.
  // Roundtrip a compatible subset without relabelling its scientific values.
  const auto first = peptides.front();
  peptides.erase(std::remove_if(peptides.begin(), peptides.end(), [&](const auto& item) {
    return item.getScoreType() != first.getScoreType() || item.isHigherScoreBetter() != first.isHigherScoreBetter();
  }), peptides.end());
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

START_SECTION([EXTRA] meta values that fields hold are stored once and restored on export)
{
  // The spectrum reference and the target/decoy states are fields of the owning model; their legacy meta values
  // are dropped on import and written back on export. A spelling or value type that export would not reproduce stays.
  auto peptides = PeptideIdentificationList {peptide()};
  auto other = peptides.front().getHits().front();
  other.setSequence(AASequence::fromString("PEPTIDER"));
  other.setMetaValue("target_decoy", "Target");
  peptides.front().insertHit(other);
  auto integer_reference = peptide();
  integer_reference.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, 42);
  peptides.push_back(integer_reference);
  auto proteins = std::vector<ProteinIdentification> {protein()};
  proteins.front().getHits().front().setMetaValue("target_decoy", "decoy");

  const auto imported = Adapter::importLegacy(proteins, peptides);
  const auto& run = imported.data.getRuns().front();
  const auto& query = run.getIdentification(imported.queries[0].query);
  TEST_EQUAL(query.data_id, "controllerType=0 controllerNumber=1 scan=42")
  TEST_FALSE(query.metaValueExists(Constants::UserParam::SPECTRUM_REFERENCE))
  TEST_TRUE(query.getMatches()[0].target_decoy == ID::TargetDecoy::TARGET)
  TEST_FALSE(query.getMatches()[0].metaValueExists("target_decoy"))
  TEST_TRUE(query.getMatches()[1].target_decoy == ID::TargetDecoy::TARGET)
  TEST_EQUAL(query.getMatches()[1].getMetaValue("target_decoy").toString(), "Target")
  const auto& integer_query = run.getIdentification(imported.queries[1].query);
  TEST_EQUAL(integer_query.data_id, "42")
  TEST_EQUAL(integer_query.getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE).valueType() == DataValue::INT_VALUE, true)
  const auto& sequence = run.getDatabaseSequences()->front();
  TEST_TRUE(sequence.target_decoy == ID::TargetDecoy::DECOY)
  TEST_FALSE(sequence.metaValueExists("target_decoy"))

  auto exported = Adapter::toLegacy(imported.data);
  TEST_EQUAL(exported.peptides.size(), peptides.size())
  TEST_TRUE(exported.peptides[0] == peptides[0])
  TEST_EQUAL(exported.proteins.front().getHits().front().getMetaValue("target_decoy").toString(), "decoy")
}
END_SECTION

START_SECTION([EXTRA] derived scores belong to their recorded producer and not to the search engine)
{
  // Two engines whose PSMs were rescored with posterior error probabilities.
  auto comet = protein();
  comet.setIdentifier("comet");
  comet.setSearchEngine("Comet");
  comet.setSearchEngineVersion("2024.01");
  auto msgf = protein();
  msgf.setIdentifier("msgf");
  msgf.setSearchEngine("MSGFPlus");
  msgf.setSearchEngineVersion("2024.07");
  auto from_comet = peptide("Posterior Error Probability", 0.01);
  from_comet.setIdentifier("comet");
  auto from_msgf = peptide("Posterior Error Probability", 0.02);
  from_msgf.setIdentifier("msgf");
  PeptideIdentificationList peptides {from_comet, from_msgf};

  // Without a recorded producer the PEPs are attributed to different engines and cannot share a schema.
  TEST_EXCEPTION(Exception::InvalidValue, Adapter::fromLegacy({comet, msgf}, peptides))

  comet.setScoreSoftware("Posterior Error Probability", "IDPosteriorErrorProbability", "3.7.0");
  msgf.setScoreSoftware("Posterior Error Probability", "IDPosteriorErrorProbability", "3.7.0");
  const auto data = Adapter::fromLegacy({comet, msgf}, peptides);
  TEST_EQUAL(data.getRuns().size(), 2)
  const auto& definition = data.getScoreDefinitions().at(0);
  TEST_EQUAL(definition.name, "Posterior Error Probability")
  TEST_EQUAL(definition.software, "IDPosteriorErrorProbability")
  TEST_EQUAL(definition.software_version, "3.7.0")
  TEST_EQUAL(data.getRuns()[0].getScoreDefinitions() == data.getRuns()[1].getScoreDefinitions(), true)
  // The engines themselves remain in the run configuration.
  TEST_EQUAL(data.getRun("comet").getSettings().software, "Comet")

  // Strict export keeps both the engine and the producer.
  const auto exported = Adapter::toLegacy(data);
  TEST_EQUAL(exported.losses.empty(), true)
  TEST_EQUAL(exported.proteins.size(), 2)
  for (const auto& run : exported.proteins)
  {
    TEST_EQUAL(run.getSearchEngine() == "Comet" || run.getSearchEngine() == "MSGFPlus", true)
    TEST_EQUAL(run.getScoreSoftware("Posterior Error Probability").first, "IDPosteriorErrorProbability")
  }

  // Different producer versions are different definitions and are still rejected.
  msgf.setScoreSoftware("Posterior Error Probability", "IDPosteriorErrorProbability", "3.6.0");
  TEST_EXCEPTION(Exception::InvalidValue, Adapter::fromLegacy({comet, msgf}, peptides))
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
  TEST_TRUE(run.getSettings().search.metaValueExists(Constants::UserParam::MODIFICATION_DEFINITIONS))
  const auto exported = Adapter::toLegacy(data);
  TEST_TRUE(exported.peptides[0] == item)
  auto settings = run.getSettings();
  definition.setDiffMonoMass(definition.getDiffMonoMass() + 1.0);
  settings.search.setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, definition.toDefinitionString());
  run.setSettings(settings);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::toLegacy(data))
}
END_SECTION

START_SECTION([EXTRA] the sources of a run are the files of its legacy run)
{
  // A file listed twice keeps both positions, and a file without identifications keeps its source.
  auto legacy_run = protein();
  legacy_run.setPrimaryMSRunPath({"a.mzML", "b.mzML", "a.mzML"});
  legacy_run.setPrimaryMSRunPath({"a.raw", "b.raw"}, true);
  auto third = peptide();
  third.setMetaValue("id_merge_index", 2);
  auto first = peptide("PEP", 0.2);
  first.setMetaValue("id_merge_index", 0);
  auto imported = Adapter::fromLegacy({legacy_run}, {third, first});
  const auto& run = imported.getRuns()[0];
  ABORT_IF(run.getSources().size() != 3)
  TEST_EQUAL(run.getSources()[0].file.path, "a.mzML")
  TEST_EQUAL(run.getSources()[1].file.path, "b.mzML")
  TEST_EQUAL(run.getSources()[2].file.path, "a.mzML")
  TEST_EQUAL(run.getSources()[0].identifications.size(), 1)
  TEST_EQUAL(run.getSources()[1].identifications.size(), 0)
  TEST_EQUAL(run.getSources()[2].identifications.size(), 1)
  TEST_EQUAL(run.getSources()[2].identifications[0].metaValueExists("id_merge_index"), false)
  TEST_EQUAL(run.getSettings().metaValueExists("spectra_data"), false)
  StringList raw;
  raw = run.getSettings().getMetaValue("spectra_data_raw").toStringList();
  TEST_EQUAL(ListUtils::concatenate(raw, ","), "a.raw,b.raw")
  TEST_EQUAL(ListUtils::concatenate(Adapter::legacyFiles(run), ","), "a.mzML,b.mzML,a.mzML")

  // Export rebuilds the file list and the indices from the sources.
  auto exported = Adapter::toLegacy(imported);
  TEST_TRUE(exported.losses.empty())
  StringList files;
  exported.proteins[0].getPrimaryMSRunPath(files);
  TEST_EQUAL(ListUtils::concatenate(files, ","), "a.mzML,b.mzML,a.mzML")
  ABORT_IF(exported.peptides.size() != 2)
  TEST_TRUE(exported.peptides[0] == first)
  TEST_TRUE(exported.peptides[1] == third)

  // A stale index in the metadata of an identification does not override its source.
  ID stale;
  auto& stale_run = stale.addRun("stale");
  stale_run.addScore({"PEP", "", false});
  stale_run.setPrimaryScore(stale_run.getScoreId(0));
  ID::SourceFile file;
  file.path = "a.mzML";
  stale_run.addSource(file);
  file.path = "b.mzML";
  const auto b = stale_run.addSource(file);
  ID::Observation observation;
  observation.setMetaValue("id_merge_index", 0);
  ID::MatchData match;
  match.representation = "PEPTIDE";
  stale_run.addMatch(stale_run.addIdentification(b, observation), match, {0.1});
  exported = Adapter::toLegacy(stale);
  TEST_EQUAL(static_cast<Int>(exported.peptides[0].getMetaValue("id_merge_index")), 1)

  // A single-file run needs no index: its identifications belong to its file either way.
  legacy_run.setPrimaryMSRunPath({"a.mzML"});
  imported = Adapter::fromLegacy({legacy_run}, {peptide(), first});
  TEST_EQUAL(imported.getRuns()[0].getSources().size(), 1)
  TEST_EQUAL(imported.getRuns()[0].getSources()[0].identifications.size(), 2)
  exported = Adapter::toLegacy(imported);
  TEST_EQUAL(exported.peptides[1].metaValueExists("id_merge_index"), false)

  // Without files, the identifications have an unknown source, and the export lists no files.
  legacy_run.removeMetaValue("spectra_data");
  imported = Adapter::fromLegacy({legacy_run}, {peptide()});
  TEST_EQUAL(imported.getRuns()[0].getSources().size(), 1)
  TEST_EQUAL(imported.getRuns()[0].getSources()[0].file.path, "")
  exported = Adapter::toLegacy(imported);
  TEST_EQUAL(exported.proteins[0].metaValueExists("spectra_data"), false)
  first.setMetaValue("id_merge_index", 0);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::fromLegacy({legacy_run}, {first}))

  // The files of a run are its sources, so its settings cannot list them as well.
  ID::Run files_in_metadata("search");
  auto settings = Adapter::settingsFromLegacy(protein());
  TEST_EQUAL(settings.metaValueExists("spectra_data"), false)
  settings.setMetaValue("spectra_data", StringList {"a.mzML"});
  TEST_EXCEPTION(Exception::InvalidValue, files_in_metadata.setSettings(settings))
}
END_SECTION

END_TEST
