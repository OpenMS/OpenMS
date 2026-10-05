// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include "all_casters.h"

#include <OpenMS/ANALYSIS/ID/IdentificationDataInference.h>
#include <OpenMS/FORMAT/IdentificationDataFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <nanobind/operators.h>

namespace pyopenms_identification
{
namespace nb = nanobind;
using ID = OpenMS::IdentificationData;
using File = OpenMS::IdentificationDataFile;
using Adapter = OpenMS::IdentificationDataAdapter;
using Inference = OpenMS::IdentificationDataInference;

// All exposed records are owned Python values. Never retain a pointer into a
// run's vector: appending/filtering its owner may invalidate such a pointer.
template<typename T, typename... Bases>
auto valueClass(nb::handle scope, const char* name)
{
  return nb::class_<T, Bases...>(scope, name, "Owning value. Nested records and containers are returned as copies; assign them back after editing.")
    .def(nb::init<>())
    .def(nb::init<const T&>())
    .def("__copy__", [](const T& self) { return T(self); })
    .def("__deepcopy__", [](const T& self, nb::dict) { return T(self); }, nb::arg("memo"));
}

template<typename Class, typename T, typename Value>
void field(Class& cls, const char* name, Value T::* member)
{
  // Explicit value return also protects nested structs and vector elements
  // when Python subsequently replaces the containing field.
  cls.def_prop_rw(
    name, [member](const T& self) -> Value { return self.*member; }, [member](T& self, const Value& value) { self.*member = value; },
    "Owned value; assign an edited nested record or container back to this property.");
}

inline void bind(nb::module_& m)
{
  auto data = valueClass<ID>(m, "IdentificationData");
  nb::enum_<ID::MoleculeKind>(data, "MoleculeKind")
    .value("PEPTIDE", ID::MoleculeKind::PEPTIDE)
    .value("OLIGONUCLEOTIDE", ID::MoleculeKind::OLIGONUCLEOTIDE)
    .value("COMPOUND", ID::MoleculeKind::COMPOUND);
  nb::enum_<ID::Encoding>(data, "Encoding")
    .value("AA_SEQUENCE", ID::Encoding::AA_SEQUENCE)
    .value("NA_SEQUENCE", ID::Encoding::NA_SEQUENCE)
    .value("SMILES", ID::Encoding::SMILES)
    .value("INCHI", ID::Encoding::INCHI)
    .value("DATABASE_ID", ID::Encoding::DATABASE_ID);
  nb::enum_<ID::TargetDecoy>(data, "TargetDecoy")
    .value("UNKNOWN", ID::TargetDecoy::UNKNOWN)
    .value("TARGET", ID::TargetDecoy::TARGET)
    .value("DECOY", ID::TargetDecoy::DECOY)
    .value("BOTH", ID::TargetDecoy::BOTH);
  nb::enum_<ID::ScoreScope>(data, "ScoreScope")
    .value("MATCH", ID::ScoreScope::MATCH)
    .value("PEPTIDE", ID::ScoreScope::PEPTIDE)
    .value("PROTEIN", ID::ScoreScope::PROTEIN)
    .value("PROTEIN_GROUP", ID::ScoreScope::PROTEIN_GROUP)
    .value("OTHER", ID::ScoreScope::OTHER);
  nb::enum_<ID::InferencePolicy>(data, "InferencePolicy")
    .value("PRESERVE", ID::InferencePolicy::PRESERVE)
    .value("DISCARD", ID::InferencePolicy::DISCARD);
  auto queryid = valueClass<ID::QueryId>(data, "QueryId");
  queryid.def_ro("value", &ID::QueryId::value).def(nb::self == nb::self).def(nb::self != nb::self);
  queryid.def(nb::init<OpenMS::UInt64>(), nb::arg("value"));
  auto matchid = valueClass<ID::MatchId>(data, "MatchId");
  matchid.def_ro("value", &ID::MatchId::value).def(nb::self == nb::self).def(nb::self != nb::self);
  matchid.def(nb::init<OpenMS::UInt64>(), nb::arg("value"));
  auto queryreference = valueClass<ID::QueryReference>(data, "QueryReference");
  field(queryreference, "run_uuid", &ID::QueryReference::run_uuid);
  field(queryreference, "query", &ID::QueryReference::query);
  auto matchreference = valueClass<ID::MatchReference>(data, "MatchReference");
  field(matchreference, "run_uuid", &ID::MatchReference::run_uuid);
  field(matchreference, "match", &ID::MatchReference::match);
  auto moleculeidentity = valueClass<ID::MoleculeIdentity>(data, "MoleculeIdentity");
  field(moleculeidentity, "encoding", &ID::MoleculeIdentity::encoding);
  field(moleculeidentity, "representation", &ID::MoleculeIdentity::representation);
  auto scoreid = valueClass<ID::ScoreId>(data, "ScoreId");
  scoreid.def_ro("value", &ID::ScoreId::value).def(nb::self == nb::self).def(nb::self != nb::self);
  auto sourceid = valueClass<ID::SourceId>(data, "SourceId");
  sourceid.def_ro("value", &ID::SourceId::value).def(nb::self == nb::self).def(nb::self != nb::self);
  auto qualifiedaccession = valueClass<ID::QualifiedAccession>(data, "QualifiedAccession");
  field(qualifiedaccession, "database", &ID::QualifiedAccession::database);
  field(qualifiedaccession, "accession", &ID::QualifiedAccession::accession);
  auto scoredefinition = valueClass<ID::ScoreDefinition>(data, "ScoreDefinition");
  field(scoredefinition, "name", &ID::ScoreDefinition::name);
  field(scoredefinition, "accession", &ID::ScoreDefinition::accession);
  field(scoredefinition, "higher_better", &ID::ScoreDefinition::higher_better);
  field(scoredefinition, "scope", &ID::ScoreDefinition::scope);
  field(scoredefinition, "software", &ID::ScoreDefinition::software);
  field(scoredefinition, "software_version", &ID::ScoreDefinition::software_version);
  field(scoredefinition, "parameters", &ID::ScoreDefinition::parameters);
  field(scoredefinition, "calibration", &ID::ScoreDefinition::calibration);
  field(scoredefinition, "aggregation", &ID::ScoreDefinition::aggregation);
  auto sourcefile = valueClass<ID::SourceFile, OpenMS::MetaInfoInterface>(data, "SourceFile");
  field(sourcefile, "identifier", &ID::SourceFile::identifier);
  field(sourcefile, "path", &ID::SourceFile::path);
  field(sourcefile, "primary_files", &ID::SourceFile::primary_files);
  auto parentevidence = valueClass<ID::ParentEvidence>(data, "ParentEvidence");
  field(parentevidence, "parent", &ID::ParentEvidence::parent);
  field(parentevidence, "start", &ID::ParentEvidence::start);
  field(parentevidence, "end", &ID::ParentEvidence::end);
  field(parentevidence, "before", &ID::ParentEvidence::before);
  field(parentevidence, "after", &ID::ParentEvidence::after);
  auto parentrecord = valueClass<ID::ParentRecord, OpenMS::MetaInfoInterface>(data, "ParentRecord");
  field(parentrecord, "identity", &ID::ParentRecord::identity);
  field(parentrecord, "sequence", &ID::ParentRecord::sequence);
  field(parentrecord, "description", &ID::ParentRecord::description);
  field(parentrecord, "target_decoy", &ID::ParentRecord::target_decoy);
  auto observation = valueClass<ID::Observation, OpenMS::MetaInfoInterface>(data, "Observation");
  field(observation, "data_id", &ID::Observation::data_id);
  field(observation, "rt", &ID::Observation::rt);
  field(observation, "mz", &ID::Observation::mz);
  auto matchdata = valueClass<ID::MatchData, OpenMS::MetaInfoInterface>(data, "MatchData");
  field(matchdata, "representation", &ID::MatchData::representation);
  field(matchdata, "encoding", &ID::MatchData::encoding);
  field(matchdata, "charge", &ID::MatchData::charge);
  field(matchdata, "calculated_mz", &ID::MatchData::calculated_mz);
  field(matchdata, "target_decoy", &ID::MatchData::target_decoy);
  field(matchdata, "name", &ID::MatchData::name);
  field(matchdata, "formula", &ID::MatchData::formula);
  field(matchdata, "identifiers", &ID::MatchData::identifiers);
  field(matchdata, "adduct", &ID::MatchData::adduct);
  field(matchdata, "parent_evidence", &ID::MatchData::parent_evidence);
  field(matchdata, "peak_annotations", &ID::MatchData::peak_annotations);
  auto inferenceinput = valueClass<ID::InferenceInput>(data, "InferenceInput");
  field(inferenceinput, "run_identifier", &ID::InferenceInput::run_identifier);
  field(inferenceinput, "run_uuid", &ID::InferenceInput::run_uuid);
  field(inferenceinput, "score", &ID::InferenceInput::score);
  field(inferenceinput, "selection", &ID::InferenceInput::selection);
  auto inferenceresult = valueClass<ID::InferenceResult>(data, "InferenceResult");
  field(inferenceresult, "identifier", &ID::InferenceResult::identifier);
  field(inferenceresult, "proteins", &ID::InferenceResult::proteins);
  field(inferenceresult, "parent_score", &ID::InferenceResult::parent_score);
  field(inferenceresult, "group_score", &ID::InferenceResult::group_score);
  field(inferenceresult, "qualified_accessions", &ID::InferenceResult::qualified_accessions);
  field(inferenceresult, "inputs", &ID::InferenceResult::inputs);

  valueClass<ID::Match, ID::MatchData>(data, "Match")
    .def("get_id", &ID::Match::getId)
    .def("get_data", [](const ID::Match& self) { return ID::MatchData(self.getData()); })
    .def("get_scores", &ID::Match::getScores);
  valueClass<ID::Identification, ID::Observation>(data, "Identification")
    .def("get_id", &ID::Identification::getId)
    .def("get_observation", [](const ID::Identification& self) { return ID::Observation(self.getObservation()); })
    .def("get_matches", [](const ID::Identification& self) { return self.getMatches(); })
    .def("get_selected_match", &ID::Identification::getSelectedMatch);
  auto block = valueClass<ID::SourceBlock>(data, "SourceBlock");
  field(block, "id", &ID::SourceBlock::id);
  field(block, "source", &ID::SourceBlock::source);
  field(block, "identifications", &ID::SourceBlock::identifications);
  valueClass<ID::ScoreView>(data, "ScoreView")
    .def("__call__", &ID::ScoreView::operator(), nb::arg("match"))
    .def("get_definition", [](const ID::ScoreView& self) { return self.getDefinition(); });

  auto run = valueClass<ID::Run>(data, "Run");
  run.def(nb::init<std::string, ID::MoleculeKind>(), nb::arg("identifier"), nb::arg("kind") = ID::MoleculeKind::PEPTIDE)
    .def("get_identifier", &ID::Run::getIdentifier)
    .def("get_uuid", &ID::Run::getUuid)
    .def("get_molecule_kind", &ID::Run::getMoleculeKind)
    .def("get_processing_metadata", [](const ID::Run& self) { return self.getProcessingMetadata(); })
    .def("set_processing_metadata", &ID::Run::setProcessingMetadata, nb::arg("metadata"))
    .def("get_parents", [](const ID::Run& self) { return self.getParents(); })
    .def("set_parents", &ID::Run::setParents, nb::arg("parents"))
    .def("get_source_blocks", [](const ID::Run& self) { return self.getSourceBlocks(); })
    .def("get_score_definitions", [](const ID::Run& self) { return self.getScoreDefinitions(); })
    .def("add_source", &ID::Run::addSource, nb::arg("source"))
    .def("get_source_id", &ID::Run::getSourceId, nb::arg("index"))
    .def("add_score", &ID::Run::addScore, nb::arg("definition"))
    .def("get_score_id", &ID::Run::getScoreId, nb::arg("index"))
    .def("find_score", &ID::Run::findScore, nb::arg("definition"))
    .def(
      "get_score_definition", [](const ID::Run& self, ID::ScoreId score) { return self.getScoreDefinition(score); }, nb::arg("score"))
    .def("bind_score", &ID::Run::bindScore, nb::arg("score"))
    .def("get_primary_score", &ID::Run::getPrimaryScore)
    .def("set_primary_score", &ID::Run::setPrimaryScore, nb::arg("score"))
    .def("add_identification", &ID::Run::addIdentification, nb::arg("source"), nb::arg("observation"))
    .def("add_match", &ID::Run::addMatch, nb::arg("query"), nb::arg("data"), nb::arg("scores") = std::vector<std::optional<double>> {})
    .def(
      "find_identification",
      [](const ID::Run& self, ID::QueryId id) -> std::optional<ID::Identification> {
        const auto* found = self.findIdentification(id);
        return found ? std::optional<ID::Identification>(*found) : std::nullopt;
      },
      nb::arg("id"))
    .def(
      "find_match",
      [](const ID::Run& self, ID::MatchId id) -> std::optional<ID::Match> {
        const auto* found = self.findMatch(id);
        return found ? std::optional<ID::Match>(*found) : std::nullopt;
      },
      nb::arg("id"))
    .def(
      "get_identification", [](const ID::Run& self, ID::QueryId id) { return self.getIdentification(id); }, nb::arg("id"))
    .def(
      "get_match", [](const ID::Run& self, ID::MatchId id) { return self.getMatch(id); }, nb::arg("id"))
    .def("get_score", &ID::Run::getScore, nb::arg("match"), nb::arg("score"))
    .def("set_score", &ID::Run::setScore, nb::arg("match"), nb::arg("score"), nb::arg("value"))
    .def("set_selected_match", &ID::Run::setSelectedMatch, nb::arg("query"), nb::arg("selected"))
    .def("replace_observation", &ID::Run::replaceObservation, nb::arg("query"), nb::arg("observation"))
    .def(
      "replace_match", [](ID::Run& self, ID::MatchId id, const ID::MatchData& value) { self.replaceMatch(id, value); }, nb::arg("match"),
      nb::arg("data"))
    .def(
      "replace_match",
      [](ID::Run& self, ID::MatchId id, const ID::MatchData& value, const std::vector<std::optional<double>>& scores) {
        self.replaceMatch(id, value, scores);
      },
      nb::arg("match"), nb::arg("data"), nb::arg("scores"))
    .def(
      "filter_matches",
      [](ID::Run& self, nb::callable keep, bool keep_empty) {
        return self.filterMatches([&](const ID::Match& match) { return nb::cast<bool>(keep(nb::cast(match, nb::rv_policy::copy))); }, keep_empty);
      },
      nb::arg("keep"), nb::arg("keep_empty_queries") = false)
    .def(
      "erase_matches",
      [](ID::Run& self, nb::callable remove, bool keep_empty) {
        return self.eraseMatches([&](const ID::Match& match) { return nb::cast<bool>(remove(nb::cast(match, nb::rv_policy::copy))); }, keep_empty);
      },
      nb::arg("remove"), nb::arg("keep_empty_queries") = false)
    .def("retain_best", &ID::Run::retainBest, nb::arg("score"), nb::arg("keep_ties") = true, nb::arg("keep_empty_queries") = false)
    .def(
      "transform_matches",
      [](ID::Run& self, nb::callable transform) {
        self.transformMatches([&](ID::MatchData& payload) {
          auto owned = nb::cast(payload, nb::rv_policy::copy);
          transform(owned);
          payload = nb::cast<ID::MatchData>(owned);
        });
      },
      nb::arg("transform"))
    .def("get_number_of_identifications", &ID::Run::getNumberOfIdentifications)
    .def("prepare_lookup_indexes", &ID::Run::prepareLookupIndexes)
    .def("get_number_of_matches", &ID::Run::getNumberOfMatches)
    .def(
      "get_identification_for_match",
      [](const ID::Run& self, ID::MatchId match) { return ID::Identification(self.getIdentificationForMatch(match)); }, nb::arg("match"))
    .def("get_next_query_id", &ID::Run::getNextQueryId)
    .def("get_next_match_id", &ID::Run::getNextMatchId)
    .def("import_identification", &ID::Run::importIdentification, nb::arg("source"), nb::arg("id"), nb::arg("observation"))
    .def("import_match", &ID::Run::importMatch, nb::arg("query"), nb::arg("id"), nb::arg("data"),
         nb::arg("scores") = std::vector<std::optional<double>> {})
    .def("restore_identity", &ID::Run::restoreIdentity, nb::arg("uuid"), nb::arg("next_query"), nb::arg("next_match"))
    .def("reserve_match_id", &ID::Run::reserveMatchId, nb::arg("id"))
    .def("validate", &ID::Run::validate);

  data
    .def(
      "add_run", [](ID& self, const std::string& identifier, ID::MoleculeKind kind) { return ID::Run(self.addRun(identifier, kind)); },
      nb::arg("identifier"), nb::arg("kind") = ID::MoleculeKind::PEPTIDE, "Add a run and return an owned copy. Commit later edits with replace_run.")
    .def(
      "add_run", [](ID& self, ID::Run value) { return ID::Run(self.addRun(std::move(value))); }, nb::arg("run"))
    .def(
      "get_run", [](const ID& self, const std::string& identifier) { return ID::Run(self.getRun(identifier)); }, nb::arg("identifier"),
      "Return an owned run copy. Commit edits with replace_run.")
    .def("replace_run", &ID::replaceRun, nb::arg("run"))
    .def(
      "find_run_by_uuid",
      [](const ID& self, const std::string& uuid) -> std::optional<ID::Run> {
        const auto* found = self.findRunByUuid(uuid);
        return found ? std::optional<ID::Run>(*found) : std::nullopt;
      },
      nb::arg("uuid"))
    .def("get_runs", [](const ID& self) { return std::vector<ID::Run>(self.getRuns().begin(), self.getRuns().end()); })
    .def("get_inference_results", [](const ID& self) { return self.getInferenceResults(); })
    .def("add_inference_result", &ID::addInferenceResult, nb::arg("result"))
    .def("clear_inference_results", &ID::clearInferenceResults)
    .def("empty", &ID::empty)
    .def("clear", &ID::clear)
    .def("merge", &ID::merge, nb::arg("other"))
    .def(nb::self == nb::self)
    .def(nb::self != nb::self)
    .def(
      "filter_matches",
      [](ID& self, nb::callable keep, ID::InferencePolicy policy, bool keep_empty) {
        return self.filterMatches([&](const ID::Match& match) { return nb::cast<bool>(keep(nb::cast(match, nb::rv_policy::copy))); }, policy,
                                  keep_empty);
      },
      nb::arg("keep"), nb::arg("inference_policy"), nb::arg("keep_empty_queries") = false)
    .def("get_score_definitions", &ID::getScoreDefinitions, nb::rv_policy::copy)
    .def("get_primary_score_definition", &ID::getPrimaryScoreDefinition)
    .def("set_primary_score", &ID::setPrimaryScore, nb::arg("definition"))
    .def("validate", &ID::validate)
    .def("swap", &ID::swap, nb::arg("other"));

  // --- IdentificationDataFile ---
  auto file = valueClass<File>(m, "IdentificationDataFile");
  auto options = valueClass<File::Options>(file, "Options");
  field(options, "batch_rows", &File::Options::batch_rows);
  field(options, "row_group_rows", &File::Options::row_group_rows);
  field(options, "batch_bytes", &File::Options::batch_bytes);
  field(options, "row_group_bytes", &File::Options::row_group_bytes);
  field(options, "max_record_bytes", &File::Options::max_record_bytes);
  field(options, "threads", &File::Options::threads);
  auto projection = valueClass<File::Projection>(file, "Projection");
  field(projection, "molecule", &File::Projection::molecule);
  field(projection, "evidence", &File::Projection::evidence);
  field(projection, "annotations", &File::Projection::annotations);
  field(projection, "metadata", &File::Projection::metadata);
  field(projection, "all_scores", &File::Projection::all_scores);
  field(projection, "score_ids", &File::Projection::score_ids);
  auto scanoptions = valueClass<File::ScanOptions>(file, "ScanOptions");
  field(scanoptions, "buffering", &File::ScanOptions::buffering);
  field(scanoptions, "projection", &File::ScanOptions::projection);
  field(scanoptions, "runs", &File::ScanOptions::runs);
  field(scanoptions, "validate_unique_ids", &File::ScanOptions::validate_unique_ids);
  auto rundescriptor = valueClass<File::RunDescriptor>(file, "RunDescriptor");
  field(rundescriptor, "identifier", &File::RunDescriptor::identifier);
  field(rundescriptor, "uuid", &File::RunDescriptor::uuid);
  field(rundescriptor, "molecule_kind", &File::RunDescriptor::molecule_kind);
  field(rundescriptor, "scores", &File::RunDescriptor::scores);
  field(rundescriptor, "sources", &File::RunDescriptor::sources);
  field(rundescriptor, "primary_score", &File::RunDescriptor::primary_score);
  field(rundescriptor, "query_count", &File::RunDescriptor::query_count);
  field(rundescriptor, "match_count", &File::RunDescriptor::match_count);
  field(rundescriptor, "next_query_id", &File::RunDescriptor::next_query_id);
  field(rundescriptor, "next_match_id", &File::RunDescriptor::next_match_id);
  auto queryrecord = valueClass<File::QueryRecord>(file, "QueryRecord");
  field(queryrecord, "query_id", &File::QueryRecord::query_id);
  field(queryrecord, "source_id", &File::QueryRecord::source_id);
  field(queryrecord, "data", &File::QueryRecord::data);
  field(queryrecord, "selected_match_id", &File::QueryRecord::selected_match_id);
  auto matchrecord = valueClass<File::MatchRecord>(file, "MatchRecord");
  field(matchrecord, "match_id", &File::MatchRecord::match_id);
  field(matchrecord, "query_id", &File::MatchRecord::query_id);
  field(matchrecord, "data", &File::MatchRecord::data);
  field(matchrecord, "scores", &File::MatchRecord::scores);
  auto scanstatistics = valueClass<File::ScanStatistics>(file, "ScanStatistics");
  field(scanstatistics, "queries", &File::ScanStatistics::queries);
  field(scanstatistics, "matches", &File::ScanStatistics::matches);
  field(scanstatistics, "descriptor_bytes", &File::ScanStatistics::descriptor_bytes);

  file
    .def_static(
      "store", [](const std::string& path, const ID& values, const File::Options& options) { File::store(path, values, options); }, nb::arg("path"),
      nb::arg("data"), nb::arg("options") = File::Options {})
    .def_static(
      "load",
      [](const std::string& path, const File::Options& options) {
        ID values;
        File::load(path, values, options);
        return values;
      },
      nb::arg("path"), nb::arg("options") = File::Options {})
    .def_static(
      "load", [](const std::string& path, ID& values, const File::Options& options) { File::load(path, values, options); }, nb::arg("path"),
      nb::arg("data"), nb::arg("options") = File::Options {})
    .def_static(
      "load_run",
      [](const std::string& path, const std::string& selected, const File::Options& options) { return File::loadRun(path, selected, options); },
      nb::arg("path"), nb::arg("run"), nb::arg("options") = File::Options {})
    .def_static("inspect", &File::inspect, nb::arg("path"))
    .def_static(
      "scan",
      [](const std::string& path, const File::ScanOptions& options, nb::object queries, nb::object matches) {
        File::QueryCallback query_callback;
        File::MatchCallback match_callback;
        if (! queries.is_none())
        {
          if (! PyCallable_Check(queries.ptr())) { throw nb::type_error("queries must be callable or None"); }
          query_callback
            = [&](const std::string& uuid, const std::vector<File::QueryRecord>& batch) { queries(uuid, nb::cast(batch, nb::rv_policy::copy)); };
        }
        if (! matches.is_none())
        {
          if (! PyCallable_Check(matches.ptr())) { throw nb::type_error("matches must be callable or None"); }
          match_callback
            = [&](const std::string& uuid, const std::vector<File::MatchRecord>& batch) { matches(uuid, nb::cast(batch, nb::rv_policy::copy)); };
        }
        return File::scan(path, options, query_callback, match_callback);
      },
      nb::arg("path"), nb::arg("options"), nb::arg("queries") = nb::none(), nb::arg("matches") = nb::none(),
      "Deliver owned batches; callbacks may retain records after the scan returns.")
    .def_static(
      "filter",
      [](const std::string& input, const std::string& output, nb::callable keep, ID::InferencePolicy policy, bool keep_empty,
         const File::Options& options) {
        File::filter(
          input, output,
          [&](const std::string& uuid, const File::MatchRecord& record) { return nb::cast<bool>(keep(uuid, nb::cast(record, nb::rv_policy::copy))); },
          policy, keep_empty, options);
      },
      nb::arg("input"), nb::arg("output"), nb::arg("keep"), nb::arg("inference_policy"), nb::arg("keep_empty_queries") = false,
      nb::arg("options") = File::Options {});

  auto oms = valueClass<OpenMS::OMSFile>(m, "OMSFile");
  oms.def(
       "store", [](OpenMS::OMSFile& self, const std::string& path, const ID& data) { self.store(path, data); }, nb::arg("path"), nb::arg("data"))
    .def(
      "store", [](OpenMS::OMSFile& self, const std::string& path, const OpenMS::FeatureMap& map) { self.store(path, map); }, nb::arg("path"),
      nb::arg("map"))
    .def(
      "store", [](OpenMS::OMSFile& self, const std::string& path, const OpenMS::ConsensusMap& map) { self.store(path, map); }, nb::arg("path"),
      nb::arg("map"))
    .def(
      "load",
      [](OpenMS::OMSFile& self, const std::string& path) {
        ID result;
        self.load(path, result);
        return result;
      },
      nb::arg("path"))
    .def(
      "load_feature_map",
      [](OpenMS::OMSFile& self, const std::string& path) {
        OpenMS::FeatureMap result;
        self.load(path, result);
        return result;
      },
      nb::arg("path"))
    .def(
      "load_consensus_map",
      [](OpenMS::OMSFile& self, const std::string& path) {
        OpenMS::ConsensusMap result;
        self.load(path, result);
        return result;
      },
      nb::arg("path"));

  // --- IdentificationDataAdapter ---
  auto adapter = valueClass<Adapter>(m, "IdentificationDataAdapter");
  nb::enum_<Adapter::LossPolicy>(adapter, "LossPolicy").value("STRICT", Adapter::LossPolicy::STRICT).value("ALLOW", Adapter::LossPolicy::ALLOW);
  nb::enum_<Adapter::MissingLinkPolicy>(adapter, "MissingLinkPolicy")
    .value("REJECT", Adapter::MissingLinkPolicy::REJECT)
    .value("PRUNE", Adapter::MissingLinkPolicy::PRUNE);
  auto adapter_exportoptions = valueClass<Adapter::ExportOptions>(adapter, "ExportOptions");
  field(adapter_exportoptions, "loss_policy", &Adapter::ExportOptions::loss_policy);
  field(adapter_exportoptions, "include_inference", &Adapter::ExportOptions::include_inference);
  field(adapter_exportoptions, "inference_result", &Adapter::ExportOptions::inference_result);
  adapter.attr("QueryReference") = data.attr("QueryReference");
  auto adapter_importresult = valueClass<Adapter::ImportResult>(adapter, "ImportResult");
  field(adapter_importresult, "data", &Adapter::ImportResult::data);
  field(adapter_importresult, "queries", &Adapter::ImportResult::queries);
  auto adapter_legacyresult = valueClass<Adapter::LegacyResult>(adapter, "LegacyResult");
  field(adapter_legacyresult, "proteins", &Adapter::LegacyResult::proteins);
  field(adapter_legacyresult, "peptides", &Adapter::LegacyResult::peptides);
  field(adapter_legacyresult, "queries", &Adapter::LegacyResult::queries);
  field(adapter_legacyresult, "losses", &Adapter::LegacyResult::losses);
  auto adapter_featureassociation = valueClass<Adapter::FeatureAssociation>(adapter, "FeatureAssociation");
  field(adapter_featureassociation, "unassigned", &Adapter::FeatureAssociation::unassigned);
  field(adapter_featureassociation, "feature_path", &Adapter::FeatureAssociation::feature_path);
  field(adapter_featureassociation, "query", &Adapter::FeatureAssociation::query);
  field(adapter_featureassociation, "matches", &Adapter::FeatureAssociation::matches);
  auto adapter_featureimportresult = valueClass<Adapter::FeatureImportResult>(adapter, "FeatureImportResult");
  field(adapter_featureimportresult, "data", &Adapter::FeatureImportResult::data);
  field(adapter_featureimportresult, "associations", &Adapter::FeatureImportResult::associations);

  adapter.def_static("import_legacy", &Adapter::importLegacy, nb::arg("proteins"), nb::arg("peptides"))
    .def_static("from_legacy", &Adapter::fromLegacy, nb::arg("proteins"), nb::arg("peptides"))
    .def_static(
      "to_legacy", [](const ID& values, const Adapter::ExportOptions& options) { return Adapter::toLegacy(values, options); }, nb::arg("data"),
      nb::arg("options") = Adapter::ExportOptions {})
    .def_static("materialize_peptide", &Adapter::materializePeptide, nb::arg("run"), nb::arg("match"), nb::arg("score"))
    .def_static("from_feature_map", &Adapter::fromFeatureMap, nb::arg("map"))
    .def_static("from_consensus_map", &Adapter::fromConsensusMap, nb::arg("map"))
    .def_static(
      "reconcile_associations",
      [](const ID& values, std::vector<Adapter::FeatureAssociation> associations, Adapter::MissingLinkPolicy policy) {
        const auto removed = Adapter::reconcileAssociations(values, associations, policy);
        return std::make_pair(std::move(associations), removed);
      },
      nb::arg("data"), nb::arg("associations"), nb::arg("policy"), "Return (updated associations, removed link count).")
    .def_static("apply_to_feature_map", &Adapter::applyToFeatureMap, nb::arg("data"), nb::arg("associations"), nb::arg("map"), nb::arg("options"),
                nb::arg("policy"))
    .def_static("apply_to_consensus_map", &Adapter::applyToConsensusMap, nb::arg("data"), nb::arg("associations"), nb::arg("map"), nb::arg("options"),
                nb::arg("policy"));

  // --- IdentificationDataInference ---
  auto inference = valueClass<Inference>(m, "IdentificationDataInference");
  nb::enum_<Inference::ProbabilityType>(inference, "ProbabilityType")
    .value("POSTERIOR_ERROR_PROBABILITY", Inference::ProbabilityType::POSTERIOR_ERROR_PROBABILITY)
    .value("POSTERIOR_PROBABILITY", Inference::ProbabilityType::POSTERIOR_PROBABILITY);
  auto inference_input = valueClass<Inference::Input>(inference, "Input");
  field(inference_input, "run_uuid", &Inference::Input::run_uuid);
  field(inference_input, "score", &Inference::Input::score);
  field(inference_input, "probability", &Inference::Input::probability);
  inference
    .def_static(
      "infer",
      [](const ID& values, const std::vector<Inference::Input>& inputs, const std::string& identifier) {
        return Inference::infer(values, inputs, identifier);
      },
      nb::arg("data"), nb::arg("inputs"), nb::arg("identifier"))
    .def_static(
      "infer",
      [](const ID& values, const std::vector<Inference::Input>& inputs, const std::string& identifier, const OpenMS::Param& parameters) {
        return Inference::infer(values, inputs, identifier, parameters);
      },
      nb::arg("data"), nb::arg("inputs"), nb::arg("identifier"), nb::arg("parameters"))
    .def_static(
      "retain_proteins",
      [](ID::InferenceResult& result, const std::vector<ID::QualifiedAccession>& retained) {
        Inference::retainProteins(result, {retained.begin(), retained.end()});
      },
      nb::arg("result"), nb::arg("retained"));
}
} // namespace pyopenms_identification
