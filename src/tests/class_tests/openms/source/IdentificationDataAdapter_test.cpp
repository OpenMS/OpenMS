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
  // The protein's sequence, description and metadata are stored once, in the database sequence; the inference
  // result keeps the inference values, and proteinHits() completes them.
  const auto& inference = imported.data.getInferenceResults()[0];
  ABORT_IF(inference.proteins.getHits().size() != 1)
  TEST_EQUAL(inference.proteins.getHits()[0].getSequence(), "")
  TEST_EQUAL(inference.proteins.getHits()[0].getDescription(), "")
  TEST_FALSE(inference.proteins.getHits()[0].metaValueExists("user"))
  TEST_REAL_SIMILAR(inference.proteins.getHits()[0].getScore(), 0.95)
  TEST_EQUAL(run.getDatabaseSequences()->at(0).sequence, "PEPTIDE")
  TEST_EQUAL(run.getDatabaseSequences()->at(0).description, "description")
  TEST_TRUE(Adapter::proteinHits(imported.data, inference) == without_files.getHits())
  auto stripped = inference.proteins;
  stripped.setHits(without_files.getHits());
  TEST_TRUE(stripped == without_files)
  const auto exported = Adapter::toLegacy(imported.data);
  TEST_TRUE(exported.losses.empty())
  TEST_EQUAL(exported.proteins.size(), 1)
  TEST_TRUE(exported.proteins[0] == original)
  // Identifications are exported in their legacy order, across the sources of the run (query IDs follow that order).
  TEST_TRUE(exported.peptides[0] == first)
  TEST_TRUE(exported.peptides[1] == second)
  TEST_TRUE(exported.peptides[2] == empty)
  TEST_TRUE(exported.peptides == peptides)
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

START_SECTION([EXTRA] export keeps the order of protein runs and of identifications across runs)
{
  // As in IDMerger outputs: the identifications of the second protein run come first, and a protein run without
  // identifications precedes both.
  auto empty_run = protein();
  empty_run.setIdentifier("empty");
  auto a = protein();
  a.setIdentifier("a");
  auto b = protein();
  b.setIdentifier("b");
  PeptideIdentificationList peptides;
  for (const auto* run : {&b, &b, &a, &b})
  {
    auto item = peptide("PEP", 0.01 * (peptides.size() + 1));
    item.setIdentifier(run->getIdentifier());
    item.setRT(peptides.size());
    peptides.push_back(item);
  }
  const auto imported = Adapter::importLegacy({empty_run, a, b}, peptides);
  TEST_EQUAL(imported.data.getRuns().size(), 3)
  TEST_EQUAL(imported.data.getRuns()[0].getIdentifier(), "empty")
  TEST_EQUAL(imported.data.getRuns()[1].getIdentifier(), "a")
  TEST_EQUAL(imported.data.getRuns()[2].getIdentifier(), "b")
  const auto exported = Adapter::toLegacy(imported.data);
  ABORT_IF(exported.proteins.size() != 3)
  TEST_EQUAL(exported.proteins[0].getIdentifier(), "empty")
  TEST_EQUAL(exported.proteins[1].getIdentifier(), "a")
  TEST_EQUAL(exported.proteins[2].getIdentifier(), "b")
  TEST_TRUE(exported.peptides == peptides)
  for (Size i = 0; i < peptides.size(); ++i)
    TEST_TRUE(exported.queries[i] == imported.queries[i])

  // Runs that were created independently share query IDs; such runs are exported one after the other.
  ID merged = Adapter::fromLegacy({a}, {peptides[2]});
  auto other = Adapter::fromLegacy({b}, {peptides[0], peptides[1]});
  merged.merge(other);
  const auto sequential = Adapter::toLegacy(merged);
  ABORT_IF(sequential.peptides.size() != 3)
  TEST_TRUE(sequential.peptides[0] == peptides[2])
  TEST_TRUE(sequential.peptides[1] == peptides[0])
  TEST_TRUE(sequential.peptides[2] == peptides[1])
}
END_SECTION

START_SECTION([EXTRA] legacy values that native records cannot hold survive a round trip)
{
  // Repeated and empty protein accessions: no database sequences, the inference result keeps the complete hits.
  auto repeated = protein();
  repeated.insertHit(repeated.getHits().front());
  ProteinHit unnamed;
  unnamed.setSequence("SEQ");
  repeated.insertHit(unnamed);
  // An evidence without an accession (flanking residues of an unknown protein) and target_decoy values that the
  // legacy getters do not know.
  auto item = peptide();
  auto hit = item.getHits().front();
  hit.setPeptideEvidences({PeptideEvidence("P1", 0, 6, '[', ']'), PeptideEvidence("", 3, 9, 'K', 'R')});
  hit.setMetaValue("target_decoy", "");
  item.setHits({hit});
  repeated.getHits().front().setMetaValue("target_decoy", "");
  const auto imported = Adapter::fromLegacy({repeated}, {item});
  const auto& run = imported.getRuns().front();
  TEST_FALSE(run.getDatabaseSequences().has_value())
  TEST_EQUAL(imported.getInferenceResults().front().proteins.getHits().size(), 3)
  const auto& match = run.getSources().back().identifications.front().getMatches().front();
  TEST_TRUE(match.target_decoy == ID::TargetDecoy::UNKNOWN)
  TEST_EQUAL(match.getMetaValue("target_decoy").toString(), "")
  TEST_EQUAL(match.sequence_evidence.size(), 2)
  TEST_EQUAL(match.sequence_evidence[1].accession, "")
  const auto exported = Adapter::toLegacy(imported);
  TEST_TRUE(exported.proteins.front() == repeated)
  TEST_TRUE(exported.peptides.front() == item)

  // Two bundles split from one run keep its identifier; merged, each run is exported under its own name.
  auto first = protein();
  first.setIdentifier("split");
  auto second = first;
  second.getHits().front().setAccession("P2");
  auto in_first = peptide();
  in_first.setIdentifier("split");
  auto in_second = peptide("PEP", 0.2);
  in_second.setIdentifier("split");
  auto merged = Adapter::fromLegacy({first}, {in_first});
  merged.merge(Adapter::fromLegacy({second}, {in_second}));
  const auto legacy = Adapter::toLegacy(merged);
  ABORT_IF(legacy.proteins.size() != 2)
  TEST_EQUAL(legacy.proteins[0].getIdentifier(), "split")
  TEST_EQUAL(legacy.proteins[1].getIdentifier(), merged.getRuns()[1].getIdentifier())
  TEST_NOT_EQUAL(legacy.proteins[1].getIdentifier(), "split")
  TEST_EQUAL(legacy.peptides[1].getIdentifier(), legacy.proteins[1].getIdentifier())
  TEST_EQUAL(legacy.proteins[1].getHits().front().getAccession(), "P2")
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
  // ProteinHit keeps its description as the meta value "Description".
  TEST_EQUAL(sequence.description, "description")
  TEST_FALSE(sequence.metaValueExists("Description"))

  auto exported = Adapter::toLegacy(imported.data);
  TEST_EQUAL(exported.peptides.size(), peptides.size())
  TEST_TRUE(exported.peptides[0] == peptides[0])
  TEST_EQUAL(exported.proteins.front().getHits().front().getMetaValue("target_decoy").toString(), "decoy")
  TEST_EQUAL(exported.proteins.front().getHits().front().getDescription(), "description")
}
END_SECTION

START_SECTION([EXTRA] a protein run that export rebuilds from its run is not kept as an inference result)
{
  // As a search engine writes it: the proteins of the matches, without scores, and the PSM score as score type.
  ProteinIdentification search;
  search.setIdentifier("search");
  search.setSearchEngine("test-search");
  search.setScoreType("PEP");
  search.setHigherScoreBetter(false);
  search.setPrimaryMSRunPath({"/exact/a.mzML"});
  auto parameters = search.getSearchParameters();
  parameters.db = "database.fasta";
  search.setSearchParameters(parameters);
  ProteinHit described;
  described.setAccession("P1");
  described.setDescription("first protein");
  described.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  search.insertHit(described);
  ProteinHit decoy;
  decoy.setAccession("DECOY_P1");
  decoy.setTargetDecoyType(ProteinHit::TargetDecoyType::DECOY);
  search.insertHit(decoy);
  ProteinHit empty_description;
  empty_description.setAccession("P2");
  empty_description.setDescription("");
  search.insertHit(empty_description);
  const PeptideIdentificationList peptides {peptide()};

  auto data = Adapter::fromLegacy({search}, peptides);
  TEST_EQUAL(data.getInferenceResults().size(), 0)
  const auto& sequences = *data.getRuns().front().getDatabaseSequences();
  ABORT_IF(sequences.size() != 3)
  TEST_FALSE(sequences[0].metaValueExists("Description"))
  TEST_FALSE(sequences[1].metaValueExists("Description"))
  // An explicitly empty description is metadata that export restores as such.
  TEST_TRUE(sequences[2].metaValueExists("Description"))
  auto exported = Adapter::toLegacy(data);
  ABORT_IF(exported.proteins.size() != 1)
  TEST_TRUE(exported.proteins[0] == search)
  ABORT_IF(exported.peptides.size() != 1)
  TEST_TRUE(exported.peptides[0] == peptides[0])

  // The exported protein list follows edits of the run.
  auto& run = data.getRun("search");
  auto targets = *run.getDatabaseSequences();
  std::erase_if(targets, [](const auto& sequence) { return sequence.target_decoy == ID::TargetDecoy::DECOY; });
  run.setDatabaseSequences(targets);
  TEST_EQUAL(Adapter::toLegacy(data).proteins[0].getHits().size(), 2)

  // A score type other than the PSM score (empty, or the search engine score after rescoring) is kept in the settings.
  for (const auto& [type, higher_better] : {std::pair {std::string(), true}, std::pair {std::string("hyperscore"), true}})
  {
    auto rescored = search;
    rescored.setScoreType(type);
    rescored.setHigherScoreBetter(higher_better);
    const auto imported = Adapter::fromLegacy({rescored}, peptides);
    TEST_EQUAL(imported.getInferenceResults().size(), 0)
    const auto& settings = imported.getRuns().front().getSettings();
    TEST_EQUAL(settings.getMetaValue("identification:legacy_protein_score_type").toString(), type)
    TEST_EQUAL(settings.getMetaValue("identification:legacy_protein_higher_score_better").toString(), "true")
    TEST_TRUE(Adapter::toLegacy(imported).proteins[0] == rescored)
  }
  TEST_FALSE(data.getRuns().front().getSettings().metaValueExists("identification:legacy_protein_score_type"))

  // Protein scores, or a protein run without PSMs, are inference results.
  auto scored = search;
  scored.getHits()[0].setScore(0.9);
  auto with_scores = Adapter::fromLegacy({scored}, peptides);
  ABORT_IF(with_scores.getInferenceResults().size() != 1)
  TEST_EQUAL(with_scores.getInferenceResults()[0].identifier, "legacy:search")
  TEST_FALSE(with_scores.getRuns().front().getSettings().metaValueExists("identification:legacy_protein_score_type"))
  TEST_TRUE(Adapter::toLegacy(with_scores).proteins[0] == scored)
  TEST_EQUAL(Adapter::fromLegacy({search}, {}).getInferenceResults().size(), 1)
}
END_SECTION

START_SECTION((static void keepLegacyProteinScoreType(IdentificationData::Run& run) and static std::string legacyIdentifier(const IdentificationData::Run& run)))
{
  // As a search engine writes it: the protein run takes the PSM score type.
  ProteinIdentification search;
  search.setIdentifier("search");
  search.setSearchEngine("test-search");
  search.setScoreType("PEP");
  search.setHigherScoreBetter(false);
  const PeptideIdentificationList peptides {peptide()};
  auto data = Adapter::fromLegacy({search}, peptides);
  auto& run = data.getRun("search");
  TEST_EQUAL(Adapter::legacyIdentifier(run), "search")
  // A new primary score would become the score type of the legacy protein run, unless it is kept.
  const auto previous = run.getScoreDefinition(*run.getPrimaryScore());
  ID::ScoreDefinition q_value;
  q_value.name = "q-value";
  q_value.higher_better = false;
  const auto score = run.addScore(q_value);
  for (const auto& source : run.getSources())
    for (const auto& query : source.identifications)
      for (const auto& match : query.getMatches())
        run.setScore(match.getId(), score, 0.01);
  Adapter::keepLegacyProteinScoreType(run);
  run.setPrimaryScore(score);
  data.removeScore(previous);
  const auto exported = Adapter::toLegacy(data);
  ABORT_IF(exported.proteins.size() != 1)
  TEST_EQUAL(exported.proteins[0].getScoreType(), "PEP")
  TEST_EQUAL(exported.proteins[0].isHigherScoreBetter(), false)
  TEST_EQUAL(exported.peptides[0].getScoreType(), "q-value")
  // A second call keeps what the first recorded.
  Adapter::keepLegacyProteinScoreType(run);
  TEST_EQUAL(Adapter::toLegacy(data).proteins[0].getScoreType(), "PEP")
}
END_SECTION

START_SECTION((static void replacePrimaryScore(IdentificationData& data, const IdentificationData::ScoreDefinition& definition, ...)))
{
  ProteinIdentification search;
  search.setIdentifier("search");
  search.setSearchEngine("test-search");
  search.setScoreType("PEP");
  search.setHigherScoreBetter(false);
  const PeptideIdentificationList peptides {peptide("PEP", 0.1)};
  auto data = Adapter::fromLegacy({search}, peptides);
  ID::ScoreDefinition q_value;
  q_value.name = "q-value";
  q_value.higher_better = false;
  // As legacy rescoring: the new score is the main score, the previous one the meta value "<name>_score".
  Adapter::replacePrimaryScore(data, q_value, [](const ID::Run&, const ID::Match&, double previous) { return previous / 10; }, "_score");
  TEST_EQUAL(data.getScoreDefinitions().size(), 1)
  auto exported = Adapter::toLegacy(data);
  TEST_EQUAL(exported.proteins[0].getScoreType(), "PEP")
  ABORT_IF(exported.peptides.size() != 1)
  TEST_EQUAL(exported.peptides[0].getScoreType(), "q-value")
  TEST_REAL_SIMILAR(exported.peptides[0].getHits()[0].getScore(), 0.01)
  TEST_REAL_SIMILAR(exported.peptides[0].getHits()[0].getMetaValue("PEP_score"), 0.1)
  // An existing meta value of the previous score's name with another value stays (keep_different).
  ID::ScoreDefinition pep;
  pep.name = "PEP";
  pep.higher_better = false;
  Adapter::replacePrimaryScore(data, pep, [](const ID::Run&, const ID::Match&, double) { return 0.5; }, "", true, "PEP_score");
  exported = Adapter::toLegacy(data);
  TEST_EQUAL(exported.peptides[0].getScoreType(), "PEP")
  TEST_REAL_SIMILAR(exported.peptides[0].getHits()[0].getMetaValue("PEP_score"), 0.1)
  TEST_REAL_SIMILAR(exported.peptides[0].getHits()[0].getMetaValue("PEP_score~"), 0.01)
}
END_SECTION

START_SECTION([EXTRA] a protein list with empty or repeated accessions has no database sequences)
{
  // Legacy runs allow such lists (e.g. compound identifications); database sequences need distinct accessions.
  for (const auto& accessions : {std::vector<std::string> {"P1", "P1"}, std::vector<std::string> {"P1", ""}})
  {
    ProteinIdentification search;
    search.setIdentifier("search");
    search.setSearchEngine("test-search");
    search.setScoreType("PEP");
    search.setHigherScoreBetter(false);
    for (const auto& accession : accessions)
    {
      ProteinHit hit;
      hit.setAccession(accession);
      search.insertHit(hit);
    }
    const PeptideIdentificationList peptides {peptide()};
    const auto data = Adapter::fromLegacy({search}, peptides);
    TEST_FALSE(data.getRuns().front().getDatabaseSequences().has_value())
    // The inference result keeps the protein hits.
    TEST_EQUAL(data.getInferenceResults().size(), 1)
    const auto exported = Adapter::toLegacy(data);
    ABORT_IF(exported.proteins.size() != 1)
    TEST_TRUE(exported.proteins[0].getHits() == search.getHits())
    ABORT_IF(exported.peptides.size() != 1)
    TEST_TRUE(exported.peptides[0] == peptides[0])
  }
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
  // The identifications keep their legacy order, not that of their files.
  ABORT_IF(exported.peptides.size() != 2)
  TEST_TRUE(exported.peptides[0] == third)
  TEST_TRUE(exported.peptides[1] == first)

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

  // A single-file run needs no index: its identifications belong to its file either way. Export writes none,
  // so an index that names the file stays metadata, and the identifications come back as they were.
  legacy_run.setPrimaryMSRunPath({"a.mzML"});
  imported = Adapter::fromLegacy({legacy_run}, {peptide(), first});
  TEST_EQUAL(imported.getRuns()[0].getSources().size(), 1)
  TEST_EQUAL(imported.getRuns()[0].getSources()[0].identifications.size(), 2)
  exported = Adapter::toLegacy(imported);
  TEST_EQUAL(exported.peptides[0].metaValueExists("id_merge_index"), false)
  TEST_EQUAL(static_cast<Int>(exported.peptides[1].getMetaValue("id_merge_index")), 0)
  TEST_TRUE(exported.peptides[1] == first)
  // An index that does not name the file of a single-file run is not exported.
  auto& single = imported.getRun(imported.getRuns()[0].getIdentifier());
  const auto& query = single.getSources()[0].identifications[1];
  auto wrong = query.getObservation();
  wrong.setMetaValue("id_merge_index", 3);
  single.replaceObservation(query.getId(), wrong);
  exported = Adapter::toLegacy(imported);
  TEST_EQUAL(exported.peptides[1].metaValueExists("id_merge_index"), false)

  // A legacy run with an empty file list keeps it; one without a file list gets none.
  auto no_files = protein();
  no_files.setMetaValue("spectra_data", DataValue(StringList()));
  auto empty_list = Adapter::toLegacy(Adapter::fromLegacy({no_files}, {peptide()}));
  TEST_TRUE(empty_list.proteins[0] == no_files)
  no_files.removeMetaValue("spectra_data");
  auto no_list = Adapter::toLegacy(Adapter::fromLegacy({no_files}, {peptide()}));
  TEST_TRUE(no_list.proteins[0] == no_files)
  // Settings may keep an empty file list, but must not list files.
  ID::RunSettings listed;
  listed.setMetaValue("spectra_data", DataValue(StringList()));
  ID with_list;
  with_list.addRun("listed").setSettings(listed);
  listed.setMetaValue("spectra_data", DataValue(StringList {"a.mzML"}));
  TEST_EXCEPTION(Exception::InvalidValue, with_list.getRun("listed").setSettings(listed))

  // Without files, the identifications have an unknown source, and the export lists no files.
  legacy_run.removeMetaValue("spectra_data");
  imported = Adapter::fromLegacy({legacy_run}, {peptide()});
  TEST_EQUAL(imported.getRuns()[0].getSources().size(), 1)
  TEST_EQUAL(imported.getRuns()[0].getSources()[0].file.path, "")
  exported = Adapter::toLegacy(imported);
  TEST_EQUAL(exported.proteins[0].metaValueExists("spectra_data"), false)
  first.setMetaValue("id_merge_index", 0);
  TEST_EXCEPTION(Exception::InvalidParameter, Adapter::fromLegacy({legacy_run}, {first}))

  // Identifications without hits and score type get a run of their own, which import splits off the legacy run
  // and which shares its files: export lists them once (e.g. in QualityControl and ProteomicsLFQ consensus maps).
  auto one_file = protein();
  one_file.setPrimaryMSRunPath({"a.mzML"});
  PeptideIdentification unscored;
  unscored.setIdentifier("search");
  unscored.setRT(43);
  const auto split = Adapter::fromLegacy({one_file}, {peptide(), unscored});
  TEST_EQUAL(split.getRuns().size(), 2)
  const auto rejoined = Adapter::toLegacy(split);
  TEST_TRUE(rejoined.losses.empty())
  ABORT_IF(rejoined.proteins.size() != 1 || rejoined.peptides.size() != 2)
  TEST_TRUE(rejoined.proteins[0] == one_file)
  TEST_TRUE(rejoined.peptides[0] == peptide())
  TEST_TRUE(rejoined.peptides[1] == unscored)

  // The files of a run are its sources, so its settings cannot list them as well.
  ID::Run files_in_metadata("search");
  auto settings = Adapter::settingsFromLegacy(protein());
  TEST_EQUAL(settings.metaValueExists("spectra_data"), false)
  settings.setMetaValue("spectra_data", StringList {"a.mzML"});
  TEST_EXCEPTION(Exception::InvalidValue, files_in_metadata.setSettings(settings))
}
END_SECTION

END_TEST
