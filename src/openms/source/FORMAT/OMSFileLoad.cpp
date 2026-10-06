// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include "OMSIdentificationData.h"

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/OMSFileLoad.h>
#include <OpenMS/FORMAT/OMSFileStore.h> // for "raiseDBError_"
#include <SQLiteCpp/Database.h>
#include <charconv>
#include <nlohmann/json.hpp> // for JSON export
#include <sqlite3.h>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
namespace
{
  UInt64 parseRecordId(const std::string& text)
  {
    UInt64 result = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), result);
    if (error != std::errc {} || end != text.data() + text.size() || ! result)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid record ID", text);
    return result;
  }
} // namespace
// initialize lookup table:
map<std::string, std::string> OMSFileLoad::export_order_by_
  = {{"version", ""},
     {"ID_IdentifiedCompound", "molecule_id"},
     {"ID_ParentMatch", "molecule_id, parent_id, start_pos, end_pos"},
     {"ID_ParentGroup_ParentSequence", "group_id, parent_id"},
     {"ID_ProcessingStep_InputFile", "processing_step_id, input_file_id"},
     {"ID_ProcessingSoftware_AssignedScore", "software_id, score_type_order"},
     {"ID_ObservationMatch_PeakAnnotation", "parent_id, processing_step_id, peak_mz, peak_annotation"},
     {"FEAT_ConvexHull", "feature_id, hull_index, point_index"},
     {"FEAT_ObservationMatch", "feature_id"},
     {"FEAT_Query", "feature_id, run_uuid, query_id"},
     {"FEAT_MapMetaData", "unique_id"}};


OMSFileLoad::OMSFileLoad(const std::string& filename, LogType log_type): db_(make_unique<SQLite::Database>(filename))
{
  setLogType(log_type);

  // read version number:
  try
  {
    auto version = db_->execAndGet("SELECT OMSFile FROM version");
    version_number_ = version.getInt();
    if (version_number_ < 1 || version_number_ > 6)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unsupported OMS version", std::to_string(version_number_));
  }
  catch (...)
  {
    raiseDBError_(db_->getErrorMsg(), __LINE__, OPENMS_PRETTY_FUNCTION, "error reading file format version number");
  }
}


OMSFileLoad::~OMSFileLoad()
{
}


bool OMSFileLoad::isEmpty_(const SQLite::Statement& query)
{ return query.getQuery().empty(); }


DataValue OMSFileLoad::makeDataValue_(const SQLite::Statement& query)
{
  DataValue::DataType type = DataValue::EMPTY_VALUE;
  int type_index = query.getColumn("data_type_id").getInt();
  if (type_index > 0) type = DataValue::DataType(type_index - 1);
  if (type == DataValue::STRING_LIST && query.getColumn("value").getType() == SQLITE_BLOB)
  {
    const auto column = query.getColumn("value");
    const auto* bytes = static_cast<const unsigned char*>(column.getBlob());
    Size position = 0, size = column.getBytes();
    const auto integer = [&]() {
      if (size - position < 8) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Truncated string list", "");
      UInt64 result = 0;
      for (Size i = 0; i < 8; ++i)
        result |= static_cast<UInt64>(bytes[position++]) << (i * 8);
      return result;
    };
    const auto count = integer();
    if (count > (size - position) / 8) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid string list count", "");
    StringList values;
    values.reserve(count);
    for (UInt64 i = 0; i < count; ++i)
    {
      const auto length = integer();
      if (length > size - position) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid string list length", "");
      values.emplace_back(reinterpret_cast<const char*>(bytes + position), length);
      position += length;
    }
    if (position != size) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Trailing string list bytes", "");
    return DataValue(values);
  }
  std::string value = query.getColumn("value").getString();
  switch (type)
  {
    case DataValue::STRING_VALUE:
      return DataValue(value);
    case DataValue::INT_VALUE:
      return DataValue(StringUtils::toInt64(value));
    case DataValue::DOUBLE_VALUE:
      return DataValue(StringUtils::toDouble(value));
    // converting lists to std::string adds square brackets - remove them:
    case DataValue::STRING_LIST:
      value = StringUtils::substr(value, 1, value.size() - 2);
      return DataValue(ListUtils::create<std::string>(value));
    case DataValue::INT_LIST:
      value = StringUtils::substr(value, 1, value.size() - 2);
      return DataValue(ListUtils::create<int>(value));
    case DataValue::DOUBLE_LIST:
      value = StringUtils::substr(value, 1, value.size() - 2);
      return DataValue(ListUtils::create<double>(value));
    default: // DataValue::EMPTY_VALUE (avoid warning about missing return)
      return DataValue();
    }
}


  bool OMSFileLoad::prepareQueryMetaInfo_(SQLite::Statement& query,
                                          const std::string& parent_table)
  {
    std::string table_name = parent_table + "_MetaInfo";
    if (!db_->tableExists(table_name)) return false;

    std::string sql_select =
    "SELECT * FROM " + table_name + " AS MI " \
    "WHERE MI.parent_id = :id";

    if (version_number_ < 4)
    {
      sql_select =
      "SELECT * FROM " + table_name + " AS MI " \
      "JOIN DataValue AS DV ON MI.data_value_id = DV.id "   \
      "WHERE MI.parent_id = :id";
    }
    query = SQLite::Statement(*db_, sql_select);
    return true;
  }


  void OMSFileLoad::handleQueryMetaInfo_(SQLite::Statement& query,
                                         MetaInfoInterface& info,
                                         Key parent_id)
  {
    query.bind(":id", parent_id);
    while (query.executeStep())
    {
      DataValue value = makeDataValue_(query);
      if (version_number_ >= 6)
      {
        const auto unit_type = query.getColumn("unit_type").getInt();
        if (unit_type < 0 || unit_type > static_cast<int>(DataValue::OTHER))
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid metadata unit type", "");
        value.setUnitType(static_cast<DataValue::UnitType>(unit_type));
        value.setUnit(query.getColumn("unit").getInt());
      }
      info.setMetaValue(query.getColumn("name").getString(), value);
    }
    query.reset(); // get ready for new executeStep()
  }


  void OMSFileLoad::loadLegacyIdentifications_(IdentificationData& data)
  {
    // Read released SQLite schemas directly into owning values. These maps exist only
    // during import; no iterator-reference graph survives the file boundary.
    const auto rows = [&](const std::string& table, const auto& consume) {
      if (! db_->tableExists(table)) return;
      SQLite::Statement query(*db_, "SELECT * FROM " + table);
      while (query.executeStep())
        consume(query);
    };
    const auto metadata = [&](const std::string& table, Key id, MetaInfoInterface& value) {
      SQLite::Statement query(*db_, "");
      if (prepareQueryMetaInfo_(query, table)) handleQueryMetaInfo_(query, value, id);
    };
    const auto kind = [](int value) {
      if (value == 1) return ID::MoleculeKind::PEPTIDE;
      if (value == 2) return ID::MoleculeKind::COMPOUND;
      if (value == 3) return ID::MoleculeKind::OLIGONUCLEOTIDE;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid legacy molecule type", std::to_string(value));
    };
    std::map<Key, ID::ScoreDefinition> definitions;
    if (db_->tableExists("ID_ScoreType"))
    {
      SQLite::Statement query(*db_, "SELECT S.*, C.accession, C.name FROM ID_ScoreType S JOIN CVTerm C ON S.cv_term_id=C.id ORDER BY S.id");
      while (query.executeStep())
      {
        ID::ScoreDefinition definition;
        definition.name = query.getColumn("name").getString();
        definition.accession = query.getColumn("accession").getString();
        definition.higher_better = query.getColumn("higher_better").getInt();
        definitions.emplace(query.getColumn("id").getInt64(), std::move(definition));
      }
    }
    std::map<Key, ID::SourceFile> sources;
    rows("ID_InputFile", [&](auto& query) {
      ID::SourceFile source;
      source.path = query.getColumn("name").getString();
      source.identifier = query.getColumn("experimental_design_id").getString();
      source.primary_files = ListUtils::create<std::string>(query.getColumn("primary_files").getString());
      sources.emplace(query.getColumn("id").getInt64(), std::move(source));
    });
    std::map<Key, ProteinIdentification> software, parameters, steps;
    rows("ID_ProcessingSoftware", [&](auto& query) {
      auto& value = software[query.getColumn("id").getInt64()];
      value.setSearchEngine(query.getColumn("name").getString());
      value.setSearchEngineVersion(query.getColumn("version").getString());
    });
    rows("ID_DBSearchParam", [&](auto& query) {
      auto& value = parameters[query.getColumn("id").getInt64()];
      auto& param = value.getSearchParameters();
      param.db = query.getColumn("database").getString();
      param.db_version = query.getColumn("database_version").getString();
      param.taxonomy = query.getColumn("taxonomy").getString();
      param.charges = query.getColumn("charges").getString();
      param.mass_type = query.getColumn("mass_type_average").getInt() ? ProteinIdentification::PeakMassType::AVERAGE
                                                                      : ProteinIdentification::PeakMassType::MONOISOTOPIC;
      param.fixed_modifications = ListUtils::create<std::string>(query.getColumn("fixed_mods").getString());
      param.variable_modifications = ListUtils::create<std::string>(query.getColumn("variable_mods").getString());
      param.precursor_mass_tolerance = query.getColumn("precursor_mass_tolerance").getDouble();
      param.fragment_mass_tolerance = query.getColumn("fragment_mass_tolerance").getDouble();
      param.precursor_mass_tolerance_ppm = query.getColumn("precursor_tolerance_ppm").getInt();
      param.fragment_mass_tolerance_ppm = query.getColumn("fragment_tolerance_ppm").getInt();
      param.missed_cleavages = query.getColumn("missed_cleavages").getUInt();
      auto enzyme = query.getColumn("digestion_enzyme").getString();
      if (! enzyme.empty() && kind(query.getColumn("molecule_type_id").getInt()) == ID::MoleculeKind::PEPTIDE)
        param.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);
      else if (! enzyme.empty())
        param.setMetaValue("rna_enzyme", enzyme);
      param.setMetaValue("legacy:min_length", query.getColumn("min_length").getInt());
      param.setMetaValue("legacy:max_length", query.getColumn("max_length").getInt());
      if (version_number_ > 1) param.setMetaValue("legacy:enzyme_term_specificity", query.getColumn("enzyme_term_specificity").getString());
    });
    rows("ID_ProcessingStep", [&](auto& query) {
      auto id = query.getColumn("id").getInt64();
      auto& value = steps[id];
      value = software.at(query.getColumn("software_id").getInt64());
      const auto param = query.getColumn("search_param_id");
      if (! param.isNull()) value.setSearchParameters(parameters.at(param.getInt64()).getSearchParameters());
      const auto time = query.getColumn("date_time").getString();
      if (! time.empty())
      {
        DateTime dt;
        dt.set(time);
        value.setDateTime(dt);
      }
      metadata("ID_ProcessingStep", id, value);
    });
    // Preserve every processing step as typed run metadata, including input paths.
    rows("ID_ProcessingStep_InputFile", [&](auto& query) {
      auto& step = steps.at(query.getColumn("processing_step_id").getInt64());
      std::vector<std::string> paths;
      step.getPrimaryMSRunPath(paths);
      paths.push_back(sources.at(query.getColumn("input_file_id").getInt64()).path);
      step.setPrimaryMSRunPath(paths);
    });
    const auto applied = [&](const std::string& table, Key id, MetaInfoInterface& value) {
      std::map<Key, double> scores;
      if (! db_->tableExists(table + "_AppliedProcessingStep")) return scores;
      SQLite::Statement query(*db_, "SELECT * FROM " + table + "_AppliedProcessingStep WHERE parent_id=:id ORDER BY processing_step_order");
      query.bind(":id", id);
      while (query.executeStep())
      {
        const auto score = query.getColumn("score_type_id");
        if (score.isNull()) continue;
        const auto score_id = score.getInt64();
        const auto number = query.getColumn("score").getDouble();
        const auto step = query.getColumn("processing_step_id");
        value.setMetaValue("legacy:score:" + std::to_string(step.isNull() ? 0 : step.getInt64()) + ":" + definitions.at(score_id).name, number);
        scores[score_id] = number; // latest application is the materialized score
      }
      return scores;
    };
    std::set<std::string> databases;
    for (const auto& [id, step] : steps)
      if (! step.getSearchParameters().db.empty()) databases.insert(step.getSearchParameters().db);
    const std::string database = databases.size() == 1 ? *databases.begin() : std::string {};
    std::map<Key, std::pair<ID::MoleculeKind, ID::ParentRecord>> parents;
    std::map<Key, std::map<Key, double>> parent_scores;
    rows("ID_ParentSequence", [&](auto& query) {
      auto id = query.getColumn("id").getInt64();
      ID::ParentRecord parent;
      parent.identity = {database, query.getColumn("accession").getString()};
      parent.sequence = query.getColumn("sequence").getString();
      parent.description = query.getColumn("description").getString();
      parent.target_decoy = query.getColumn("is_decoy").getInt() ? ID::TargetDecoy::DECOY : ID::TargetDecoy::TARGET;
      metadata("ID_ParentSequence", id, parent);
      parent.setMetaValue("coverage", query.getColumn("coverage").getDouble());
      parent_scores[id] = applied("ID_ParentSequence", id, parent);
      parents.emplace(id, std::make_pair(kind(query.getColumn("molecule_type_id").getInt()), std::move(parent)));
    });
    std::map<Key, std::pair<ID::MoleculeKind, ID::MatchData>> molecules;
    rows("ID_IdentifiedMolecule", [&](auto& query) {
      auto id = query.getColumn("id").getInt64();
      auto type = kind(query.getColumn("molecule_type_id").getInt());
      ID::MatchData match;
      match.encoding = type == ID::MoleculeKind::PEPTIDE           ? ID::Encoding::AA_SEQUENCE
                       : type == ID::MoleculeKind::OLIGONUCLEOTIDE ? ID::Encoding::NA_SEQUENCE
                                                                   : ID::Encoding::DATABASE_ID;
      match.representation = query.getColumn("identifier").getString();
      metadata("ID_IdentifiedMolecule", id, match);
      applied("ID_IdentifiedMolecule", id, match);
      compatibility_molecules_[id] = {match.encoding, match.representation};
      molecules.emplace(id, std::make_pair(type, std::move(match)));
    });
    rows("ID_IdentifiedCompound", [&](auto& query) {
      auto& match = molecules.at(query.getColumn("molecule_id").getInt64()).second;
      match.identifiers.push_back({database, match.representation});
      match.name = query.getColumn("name").getString();
      auto formula = query.getColumn("formula").getString();
      if (! formula.empty()) match.formula = formula;
      auto smiles = query.getColumn("smile").getString(), inchi = query.getColumn("inchi").getString();
      if (! smiles.empty() && smiles != "null")
      {
        match.encoding = ID::Encoding::SMILES;
        match.representation = smiles;
      }
      else if (! inchi.empty() && inchi != "null")
      {
        match.encoding = ID::Encoding::INCHI;
        match.representation = inchi;
      }
      compatibility_molecules_[query.getColumn("molecule_id").getInt64()] = {match.encoding, match.representation};
    });
    rows("ID_ParentMatch", [&](auto& query) {
      auto& match = molecules.at(query.getColumn("molecule_id").getInt64()).second;
      const auto& parent = parents.at(query.getColumn("parent_id").getInt64()).second;
      ID::ParentEvidence evidence;
      evidence.parent = parent.identity;
      auto start = query.getColumn("start_pos"), end = query.getColumn("end_pos");
      if (! start.isNull() && start.getInt64() >= 0) evidence.start = start.getInt64();
      if (! end.isNull() && end.getInt64() >= 0) evidence.end = end.getInt64();
      evidence.before = query.getColumn("left_neighbor").getString();
      evidence.after = query.getColumn("right_neighbor").getString();
      match.parent_evidence.push_back(std::move(evidence));
      if (match.target_decoy == ID::TargetDecoy::UNKNOWN) match.target_decoy = parent.target_decoy;
      else if (match.target_decoy != parent.target_decoy)
        match.target_decoy = ID::TargetDecoy::BOTH;
    });
    std::map<Key, std::pair<Key, ID::Observation>> observations;
    rows("ID_Observation", [&](auto& query) {
      auto id = query.getColumn("id").getInt64();
      ID::Observation value;
      value.data_id = query.getColumn("data_id").getString();
      auto rt = query.getColumn("rt"), mz = query.getColumn("mz");
      if (! rt.isNull()) value.rt = rt.getDouble();
      if (! mz.isNull()) value.mz = mz.getDouble();
      metadata("ID_Observation", id, value);
      observations.emplace(id, std::make_pair(query.getColumn("input_file_id").getInt64(), std::move(value)));
    });
    std::map<Key, AdductInfo> adducts;
    rows("AdductInfo", [&](auto& query) {
      adducts.emplace(query.getColumn("id").getInt64(),
                      AdductInfo(query.getColumn("name").getString(), EmpiricalFormula(query.getColumn("formula").getString()),
                                 query.getColumn("charge").getInt(), query.getColumn("mol_multiplier").getInt()));
    });
    // Protein/group-only scores do not belong in the PSM column contract.
    std::map<Key, ID::ScoreDefinition> match_definitions;
    rows("ID_ObservationMatch_AppliedProcessingStep", [&](auto& query) {
      const auto score = query.getColumn("score_type_id");
      if (! score.isNull()) match_definitions.emplace(score.getInt64(), definitions.at(score.getInt64()));
    });
    ID loaded;
    std::map<ID::MoleculeKind, ID::Run*> runs;
    std::map<std::pair<ID::MoleculeKind, Key>, ID::QueryId> queries;
    const auto get_run = [&](ID::MoleculeKind type) -> ID::Run& {
      auto found = runs.find(type);
      if (found != runs.end()) return *found->second;
      auto& run = loaded.addRun("OMS:" + std::to_string(static_cast<int>(type)), type);
      ProteinIdentification processing;
      if (! steps.empty()) processing = steps.rbegin()->second;
      processing.setIdentifier(run.getIdentifier());
      for (const auto& [id, step] : steps)
      {
        const auto prefix = "legacy:processing:" + std::to_string(id) + ":";
        processing.setMetaValue(prefix + "software", step.getSearchEngine());
        processing.setMetaValue(prefix + "version", step.getSearchEngineVersion());
        processing.setMetaValue(prefix + "date", step.getDateTime().get());
        std::vector<std::string> paths;
        step.getPrimaryMSRunPath(paths);
        processing.setMetaValue(prefix + "inputs", paths);
        std::vector<std::string> keys;
        step.getKeys(keys);
        for (const auto& key : keys)
          processing.setMetaValue(prefix + key, step.getMetaValue(key));
      }
      run.setProcessingMetadata(processing);
      for (const auto& [id, definition] : match_definitions)
        run.addScore(definition);
      std::vector<ID::ParentRecord> selected;
      for (const auto& [id, parent] : parents)
        if (parent.first == type) selected.push_back(parent.second);
      if (! selected.empty()) run.setParents(std::move(selected));
      runs.emplace(type, &run);
      return run;
    };
    std::map<std::pair<ID::MoleculeKind, Key>, ID::SourceId> source_ids;
    const auto get_query = [&](ID::MoleculeKind type, Key observation) {
      const auto key = std::make_pair(type, observation);
      auto found = queries.find(key);
      if (found != queries.end()) return found->second;
      auto& run = get_run(type);
      const auto& item = observations.at(observation);
      const auto skey = std::make_pair(type, item.first);
      auto source = source_ids.find(skey);
      if (source == source_ids.end()) source = source_ids.emplace(skey, run.addSource(sources.at(item.first))).first;
      auto query = run.addIdentification(source->second, item.second);
      queries.emplace(key, query);
      return query;
    };
    rows("ID_ObservationMatch", [&](auto& query) {
      const auto id = query.getColumn("id").getInt64();
      auto [type, match] = molecules.at(query.getColumn("identified_molecule_id").getInt64());
      auto& run = get_run(type);
      auto observation = get_query(type, query.getColumn("observation_id").getInt64());
      match.charge = query.getColumn("charge").getInt();
      metadata("ID_ObservationMatch", id, match);
      const auto values = applied("ID_ObservationMatch", id, match);
      const auto adduct = query.getColumn("adduct_id");
      if (! adduct.isNull())
      {
        auto value = adducts.at(adduct.getInt64());
        if (value.getCharge() != match.charge && type == ID::MoleculeKind::OLIGONUCLEOTIDE)
          value = AdductInfo(value.getName(), value.getEmpiricalFormula() + EmpiricalFormula("H") * (match.charge - value.getCharge()), match.charge,
                             value.getMolMultiplier());
        match.adduct = value;
      }
      if (type == ID::MoleculeKind::COMPOUND && match.formula && match.adduct)
        match.calculated_mz = match.adduct->getMZ(EmpiricalFormula(*match.formula).getMonoWeight());
      if (db_->tableExists("ID_ObservationMatch_PeakAnnotation"))
      {
        SQLite::Statement annotations(*db_, "SELECT * FROM ID_ObservationMatch_PeakAnnotation WHERE parent_id=:id");
        annotations.bind(":id", id);
        while (annotations.executeStep())
        {
          PeptideHit::PeakAnnotation annotation;
          annotation.annotation = annotations.getColumn("peak_annotation").getString();
          annotation.charge = annotations.getColumn("peak_charge").getInt();
          annotation.mz = annotations.getColumn("peak_mz").getDouble();
          annotation.intensity = annotations.getColumn("peak_intensity").getDouble();
          match.peak_annotations.push_back(std::move(annotation));
        }
      }
      std::vector<std::optional<double>> scores;
      for (const auto& [key, definition] : match_definitions)
      {
        auto found = values.find(key);
        scores.push_back(found == values.end() ? std::nullopt : std::optional<double>(found->second));
      }
      const auto match_id = run.addMatch(observation, match, scores);
      compatibility_matches_[id] = {run.getUuid(), match_id};
    });
    // Scoreless sequence catalogs (e.g. saved RNA digestion) remain owning values.
    if (compatibility_matches_.empty() && ! molecules.empty())
      for (const auto& [id, molecule] : molecules)
      {
        auto& run = get_run(molecule.first);
        auto processing = run.getProcessingMetadata();
        processing.setMetaValue("identification:catalog", "true");
        run.setProcessingMetadata(processing);
        if (run.getSourceBlocks().empty()) run.addSource(ID::SourceFile {});
        ID::Observation observation;
        observation.data_id = "catalog=" + std::to_string(id);
        const auto query = run.addIdentification(run.getSourceId(0), observation);
        run.addMatch(query, molecule.second, std::vector<std::optional<double>>(match_definitions.size()));
      }
    if (runs.empty() && ! parents.empty())
      for (const auto& [id, parent] : parents)
        get_run(parent.first);
    // Keep empty observations too; they are not tied to a particular candidate type.
    if (! observations.empty() && runs.empty()) get_run(ID::MoleculeKind::PEPTIDE);
    for (const auto& [id, observation] : observations)
      if (std::none_of(queries.begin(), queries.end(), [id = id](const auto& item) { return item.first.second == id; }))
        get_query(runs.begin()->first, id);
    // Honor the preferred score of the most recent software when it is complete.
    std::vector<Size> preferred;
    if (db_->tableExists("ID_ProcessingSoftware_AssignedScore") && db_->tableExists("ID_ProcessingStep"))
    {
      SQLite::Statement query(*db_, "SELECT A.score_type_id FROM ID_ProcessingSoftware_AssignedScore A JOIN ID_ProcessingStep P ON "
                                    "A.software_id=P.software_id ORDER BY P.date_time DESC, A.score_type_order");
      while (query.executeStep())
      {
        auto found = match_definitions.find(query.getColumn(0).getInt64());
        if (found != match_definitions.end()) preferred.push_back(std::distance(match_definitions.begin(), found));
      }
    }
    for (Size column = match_definitions.size(); column > 0; --column)
      preferred.push_back(column - 1);
    for (Size column : preferred)
    {
      bool complete = true;
      for (const auto& run : loaded.getRuns())
        for (const auto& source : run.getSourceBlocks())
          for (const auto& query : source.identifications)
            for (const auto& match : query.getMatches())
              complete = complete && run.getScore(match.getId(), run.getScoreId(column)).has_value();
      if (complete)
      {
        for (const auto& run : loaded.getRuns())
          loaded.getRun(run.getIdentifier()).setPrimaryScore(run.getScoreId(column));
        break;
      }
    }
    for (const auto& run : loaded.getRuns())
      if (run.getNumberOfMatches() && ! run.getScoreDefinitions().empty() && ! run.getPrimaryScore())
        throw Exception::InvalidValue(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Legacy OMS has no common, fully populated primary PSM score; split score layouts or select a common score before conversion",
          run.getIdentifier());
    for (const auto& [type, run] : runs)
    {
      ID::InferenceResult inference;
      inference.identifier = "OMS:parents:" + run->getIdentifier();
      inference.proteins = run->getProcessingMetadata();
      std::optional<Key> chosen;
      for (const auto& [id, parent] : parents)
        if (parent.first == type && ! parent_scores.at(id).empty()) chosen = parent_scores.at(id).rbegin()->first;
      if (! chosen || db_->tableExists("ID_ParentGroupSet")) continue;
      auto definition = definitions.at(*chosen);
      definition.scope = ID::ScoreScope::PROTEIN;
      inference.parent_score = definition;
      inference.proteins.setScoreType(definition.name);
      inference.proteins.setHigherScoreBetter(definition.higher_better);
      std::vector<ProteinHit> hits;
      for (const auto& [id, parent] : parents)
        if (parent.first == type)
        {
          ProteinHit hit;
          static_cast<MetaInfoInterface&>(hit) = parent.second;
          hit.setAccession(parent.second.identity.accession);
          hit.setSequence(parent.second.sequence);
          hit.setDescription(parent.second.description);
          hit.setCoverage(static_cast<double>(parent.second.getMetaValue("coverage")) * 100.0);
          hit.setTargetDecoyType(parent.second.target_decoy == ID::TargetDecoy::DECOY ? ProteinHit::TargetDecoyType::DECOY
                                                                                      : ProteinHit::TargetDecoyType::TARGET);
          const auto score = parent_scores.at(id).find(*chosen);
          if (score != parent_scores.at(id).end()) hit.setScore(score->second);
          else
            hit.setMetaValue("legacy:missing_parent_score", "true");
          inference.qualified_accessions[hit.getAccession()] = parent.second.identity;
          hits.push_back(std::move(hit));
        }
      inference.proteins.setHits(hits);
      inference.inputs.push_back({run->getIdentifier(), run->getUuid(), std::nullopt, "Legacy OMS parent scores"});
      loaded.addInferenceResult(std::move(inference));
    }
    std::vector<ID::InferenceResult> group_inferences;
    std::map<std::vector<Key>, Size> protein_groupings;
    rows("ID_ParentGroupSet", [&](auto& query) {
      const auto id = query.getColumn("id").getInt64();
      ID::InferenceResult inference;
      inference.identifier = "OMS:groups:" + std::to_string(id);
      inference.proteins.setIdentifier(query.getColumn("label").getString());
      metadata("ID_ParentGroupSet", id, inference.proteins);
      applied("ID_ParentGroupSet", id, inference.proteins);
      std::optional<Key> parent_primary;
      for (const auto& [parent_id, values] : parent_scores)
        if (! values.empty()) parent_primary = values.rbegin()->first;
      if (parent_primary)
      {
        auto definition = definitions.at(*parent_primary);
        definition.scope = ID::ScoreScope::PROTEIN;
        inference.parent_score = definition;
        inference.proteins.setScoreType(definition.name);
        inference.proteins.setHigherScoreBetter(definition.higher_better);
      }
      std::vector<ProteinHit> hits;
      for (const auto& [parent_id, parent] : parents)
      {
        ProteinHit hit;
        static_cast<MetaInfoInterface&>(hit) = parent.second;
        hit.setAccession(parent.second.identity.accession);
        hit.setSequence(parent.second.sequence);
        hit.setDescription(parent.second.description);
        hit.setCoverage(static_cast<double>(parent.second.getMetaValue("coverage")) * 100.0);
        if (parent_primary)
        {
          auto score = parent_scores.at(parent_id).find(*parent_primary);
          if (score != parent_scores.at(parent_id).end()) hit.setScore(score->second);
          else
            hit.setMetaValue("legacy:missing_parent_score", "true");
        }
        hit.setTargetDecoyType(parent.second.target_decoy == ID::TargetDecoy::DECOY ? ProteinHit::TargetDecoyType::DECOY
                                                                                    : ProteinHit::TargetDecoyType::TARGET);
        hits.push_back(hit);
        inference.qualified_accessions[hit.getAccession()] = parent.second.identity;
      }
      inference.proteins.setHits(hits);
      SQLite::Statement groups(*db_, "SELECT * FROM ID_ParentGroup WHERE grouping_id=:id ORDER BY id");
      groups.bind(":id", id);
      std::map<Key, ProteinIdentification::ProteinGroup> selected;
      while (groups.executeStep())
      {
        auto group_id = groups.getColumn("id").getInt64();
        auto& group = selected[group_id];
        const auto score = groups.getColumn("score_type_id");
        if (! score.isNull())
        {
          auto definition = definitions.at(score.getInt64());
          definition.scope = ID::ScoreScope::PROTEIN_GROUP;
          if (inference.group_score && *inference.group_score != definition)
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "Legacy group contains multiple score types; export separate inference results first");
          inference.group_score = definition;
          group.probability = groups.getColumn("score").getDouble();
        }
      }
      for (auto& [group_id, group] : selected)
      {
        SQLite::Statement members(*db_, "SELECT parent_id FROM ID_ParentGroup_ParentSequence WHERE group_id=:id");
        members.bind(":id", group_id);
        while (members.executeStep())
          group.accessions.push_back(parents.at(members.getColumn(0).getInt64()).second.identity.accession);
        if (query.getColumn("label").getString() == "indistinguishable proteins")
          inference.proteins.getIndistinguishableProteins().push_back(std::move(group));
        else
          inference.proteins.getProteinGroups().push_back(std::move(group));
      }
      for (const auto& run : loaded.getRuns())
        inference.inputs.push_back({run.getIdentifier(), run.getUuid(), std::nullopt, "Legacy OMS grouping"});
      const auto label = query.getColumn("label").getString();
      std::vector<Key> step_ids;
      if (db_->tableExists("ID_ParentGroupSet_AppliedProcessingStep"))
      {
        SQLite::Statement steps(
          *db_, "SELECT DISTINCT processing_step_id FROM ID_ParentGroupSet_AppliedProcessingStep WHERE parent_id=:id ORDER BY processing_step_id");
        steps.bind(":id", id);
        while (steps.executeStep())
          step_ids.push_back(steps.getColumn(0).isNull() ? 0 : steps.getColumn(0).getInt64());
      }
      inference.proteins.setMetaValue("legacy:grouping:" + std::to_string(id) + ":label", label);
      const bool conventional = ! step_ids.empty() && (label == "protein groups" || label == "indistinguishable proteins");
      const auto previous = protein_groupings.find(step_ids);
      if (conventional && previous != protein_groupings.end())
      {
        auto& combined = group_inferences[previous->second];
        if (combined.parent_score != inference.parent_score || combined.group_score != inference.group_score)
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Incompatible scores in paired legacy protein groupings");
        auto& groups = combined.proteins.getProteinGroups();
        const auto& added = inference.proteins.getProteinGroups();
        groups.insert(groups.end(), added.begin(), added.end());
        auto& indistinguishable = combined.proteins.getIndistinguishableProteins();
        const auto& added_indistinguishable = inference.proteins.getIndistinguishableProteins();
        indistinguishable.insert(indistinguishable.end(), added_indistinguishable.begin(), added_indistinguishable.end());
        combined.proteins.setMetaValue("legacy:grouping:" + std::to_string(id) + ":label", label);
      }
      else
      {
        if (conventional) protein_groupings.emplace(step_ids, group_inferences.size());
        group_inferences.push_back(std::move(inference));
      }
    });
    for (auto& inference : group_inferences)
      loaded.addInferenceResult(std::move(inference));
    loaded.validate();
    data = std::move(loaded);
    identification_data_ = &data;
  }

  void OMSFileLoad::load(IdentificationData& data)
  {
    if (version_number_ < 6)
    {
      loadLegacyIdentifications_(data);
      return;
    }
    loadOMSIdentifications(*db_, data);
    identification_data_ = &data;
  }

  template <class MapType>
  std::string OMSFileLoad::loadMapMetaDataTemplate_(MapType& features)
  {
    if (!db_->tableExists("FEAT_MapMetaData")) return "";

    SQLite::Statement query(*db_, "SELECT * FROM FEAT_MapMetaData");
    query.executeStep(); // there should be only one row
    Key id = query.getColumn("unique_id").getInt64();
    features.setUniqueId(id);
    features.setIdentifier(query.getColumn("identifier").getString());
    features.setLoadedFilePath(query.getColumn("file_path").getString());
    // The "file_type" column stores a FileTypes type *name* (e.g. "featureXML").
    // setLoadedFileType() takes a file *path* and detects the type via its
    // content, so passing the type name made it try to open a file called
    // "featureXML" (FileNotFound). There is no public setter for the
    // FileTypes::Type enum; the loaded file path is already restored above, so
    // we leave the (rarely used) loaded file type at its default here, matching
    // develop's effective behaviour.
    SQLite::Statement query_meta(*db_, "");
    if (prepareQueryMetaInfo_(query_meta, "FEAT_MapMetaData"))
    {
      handleQueryMetaInfo_(query_meta, features, id);
    }
    if (version_number_ < 5) return ""; // "experiment_type" column doesn't exist yet
    return query.getColumn("experiment_type").getString(); // for consensus map only
  }

  // template specializations:
  template std::string OMSFileLoad::loadMapMetaDataTemplate_<FeatureMap>(FeatureMap&);
  template std::string OMSFileLoad::loadMapMetaDataTemplate_<ConsensusMap>(ConsensusMap&);


  void OMSFileLoad::loadMapMetaData_(FeatureMap& features)
  {
    loadMapMetaDataTemplate_(features);
  }

  void OMSFileLoad::loadMapMetaData_(ConsensusMap& consensus)
  {
    std::string experiment_type = loadMapMetaDataTemplate_(consensus);
    consensus.setExperimentType(experiment_type);
  }


  void OMSFileLoad::loadDataProcessing_(vector<DataProcessing>& data_processing)
  {
    if (!db_->tableExists("FEAT_DataProcessing")) return;

    // "position" column was removed in schema version 3:
    std::string order_by = version_number_ > 2 ? "id" : "position";
    SQLite::Statement query(*db_, "SELECT * FROM FEAT_DataProcessing ORDER BY " + order_by + " ASC");

    SQLite::Statement subquery_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info, "FEAT_DataProcessing");

    while (query.executeStep())
    {
      DataProcessing proc;
      Software sw(query.getColumn("software_name").getString(),
                  query.getColumn("software_version").getString());
      proc.setSoftware(sw);
      vector<std::string> actions =
        ListUtils::create<std::string>(query.getColumn("processing_actions").getString());
      for (const std::string& action : actions)
      {
        auto pos = find(begin(DataProcessing::NamesOfProcessingAction),
                          end(DataProcessing::NamesOfProcessingAction), action);
        if (pos != end(DataProcessing::NamesOfProcessingAction))
        {
          Size index = pos - begin(DataProcessing::NamesOfProcessingAction);
          proc.getProcessingActions().insert(DataProcessing::ProcessingAction(index));
        }
        else // @TODO: throw an exception here?
        {
          OPENMS_LOG_ERROR << "Error: unknown data processing action '" << action << "' - skipping";
        }
      }
      DateTime time;
      time.set(query.getColumn("completion_time").getString());
      proc.setCompletionTime(time);
      if (have_meta_info)
      {
        Key id = query.getColumn("id").getInt64();
        handleQueryMetaInfo_(subquery_info, proc, id);
      }
      data_processing.push_back(proc);
    }
  }


  BaseFeature OMSFileLoad::makeBaseFeature_(int id, SQLite::Statement& query_feat,
                                            SQLite::Statement& query_meta,
                                            SQLite::Statement& query_match)
  {
    BaseFeature feature;
    feature.setRT(query_feat.getColumn("rt").getDouble());
    feature.setMZ(query_feat.getColumn("mz").getDouble());
    feature.setIntensity(query_feat.getColumn("intensity").getDouble());
    feature.setCharge(query_feat.getColumn("charge").getInt());
    feature.setWidth(query_feat.getColumn("width").getDouble());
    // setWidth adds the featureXML compatibility key; only persisted metadata
    // should be restored when loading OMS.
    feature.removeMetaValue("FWHM");
    string quality_column = (version_number_ < 5) ? "overall_quality" : "quality";
    feature.setQuality(query_feat.getColumn(quality_column.c_str()).getDouble());
    feature.setUniqueId(query_feat.getColumn("unique_id").getInt64());
    if (id == -1) return feature; // stop here for feature handles (in consensus maps)

    if (version_number_ >= 6)
    {
      auto encoding = query_feat.getColumn("primary_encoding");
      if (! encoding.isNull())
      {
        auto value = encoding.getInt();
        if (value < 0 || value > static_cast<int>(ID::Encoding::DATABASE_ID))
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid molecular encoding", "");
        feature.setPrimaryID({static_cast<ID::Encoding>(value), query_feat.getColumn("primary_representation").getString()});
      }
    }
    if (version_number_ < 6)
    {
      auto primary = query_feat.getColumn("primary_molecule_id");
      if (! primary.isNull()) feature.setPrimaryID(compatibility_molecules_.at(primary.getInt64()));
    }
    // meta data:
    if (!isEmpty_(query_meta))
    {
      handleQueryMetaInfo_(query_meta, feature, id);
    }
    // ID matches:
    if (!isEmpty_(query_match))
    {
      query_match.bind(":id", id);
      while (query_match.executeStep())
      {
        ID::MatchReference reference = version_number_ < 6
                                         ? compatibility_matches_.at(query_match.getColumn("observation_match_id").getInt64())
                                         : ID::MatchReference {query_match.getColumn("run_uuid").getString(),
                                                               ID::MatchId {parseRecordId(query_match.getColumn("match_id").getString())}};
        const auto* run = identification_data_ ? identification_data_->findRunByUuid(reference.run_uuid) : nullptr;
        if (! run || ! run->findMatch(reference.match))
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing feature match", "");
        feature.addIDMatch(reference);
      }
      query_match.reset(); // get ready for new executeStep()
    }
    if (db_->tableExists("FEAT_Query"))
    {
      SQLite::Statement links(*db_, "SELECT run_uuid, query_id FROM FEAT_Query WHERE feature_id = :id");
      links.bind(":id", id);
      while (links.executeStep())
      {
        ID::QueryReference reference {links.getColumn("run_uuid").getString(), ID::QueryId {parseRecordId(links.getColumn("query_id").getString())}};
        const auto* run = identification_data_ ? identification_data_->findRunByUuid(reference.run_uuid) : nullptr;
        if (! run || ! run->findIdentification(reference.query))
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing feature query", "");
        feature.addIDQuery(reference);
      }
    }
    return feature;
  }


  void OMSFileLoad::prepareQueriesBaseFeature_(SQLite::Statement& query_meta,
                                               SQLite::Statement& query_match)
  {
    // the "main" query is different for Feature/ConsensusFeature, so don't include it here
    string main_table = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature";
    prepareQueryMetaInfo_(query_meta, main_table);
    if (db_->tableExists("FEAT_ObservationMatch"))
    {
      query_match = SQLite::Statement(*db_, "SELECT * FROM FEAT_ObservationMatch WHERE feature_id = :id");
    }
  }


  Feature OMSFileLoad::loadFeatureAndSubordinates_(
    SQLite::Statement& query_feat, SQLite::Statement& query_meta,
    SQLite::Statement& query_match, SQLite::Statement& query_hull)
  {
    int id = query_feat.getColumn("id").getInt();
    Feature feature(makeBaseFeature_(id, query_feat, query_meta, query_match));
    // Feature-specific attributes:
    feature.setQuality(0, query_feat.getColumn("rt_quality").getDouble());
    feature.setQuality(1, query_feat.getColumn("mz_quality").getDouble());
    // convex hulls:
    if (!isEmpty_(query_hull))
    {
      query_hull.bind(":id", id);
      while (query_hull.executeStep())
      {
        Size hull_index = query_hull.getColumn("hull_index").getUInt();
        // first row should have max. hull index (sorted descending):
        if (feature.getConvexHulls().size() <= hull_index)
        {
          feature.getConvexHulls().resize(hull_index + 1);
        }
        ConvexHull2D::PointType point(query_hull.getColumn("point_x").getDouble(),
                                      query_hull.getColumn("point_y").getDouble());
        // @TODO: this may be inefficient (see implementation of "addPoint"):
        feature.getConvexHulls()[hull_index].addPoint(point);
      }
      query_hull.reset(); // get ready for new executeStep()
    }
    // subordinates:
    string from = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature JOIN FEAT_Feature ON id = feature_id";
    SQLite::Statement query_sub(*db_, "SELECT * FROM " + from + " WHERE subordinate_of = " + StringUtils::toStr(id) + " ORDER BY id ASC");
    while (query_sub.executeStep())
    {
      Feature sub = loadFeatureAndSubordinates_(query_sub, query_meta,
                                                query_match, query_hull);
      feature.getSubordinates().push_back(sub);
    }
    return feature;
  }


  void OMSFileLoad::loadFeatures_(FeatureMap& features)
  {
    if (!db_->tableExists("FEAT_Feature")) return;

    // start with top-level features only:
    string from = (version_number_ < 5) ? "FEAT_Feature" : "FEAT_BaseFeature JOIN FEAT_Feature ON id = feature_id";
    SQLite::Statement query_feat(*db_, "SELECT * FROM " + from + " WHERE subordinate_of IS NULL ORDER BY id ASC");
    // prepare sub-queries (optional - corresponding tables may not be present):
    SQLite::Statement query_meta(*db_, "");
    SQLite::Statement query_match(*db_, "");
    prepareQueriesBaseFeature_(query_meta, query_match);
    SQLite::Statement query_hull(*db_, "");
    if (db_->tableExists("FEAT_ConvexHull"))
    {
      query_hull = SQLite::Statement(*db_, "SELECT * FROM FEAT_ConvexHull WHERE feature_id = :id " \
                                     "ORDER BY hull_index DESC, point_index ASC");
    }

    while (query_feat.executeStep())
    {
      Feature feature = loadFeatureAndSubordinates_(query_feat, query_meta,
                                                    query_match, query_hull);
      features.push_back(feature);
    }
  }


  void OMSFileLoad::load(FeatureMap& features)
  {
    load(features.getIdentificationData()); // load IDs, if any
    startProgress(0, 3, "Reading feature data from file");
    loadMapMetaData_(features);
    nextProgress();
    loadDataProcessing_(features.getDataProcessing());
    nextProgress();
    loadFeatures_(features);
    features.updateRanges();
    endProgress();
  }


  void OMSFileLoad::loadConsensusFeatures_(ConsensusMap& consensus)
  {
    if (!db_->tableExists("FEAT_FeatureHandle")) return;

    // start with top-level features only:
    SQLite::Statement query_feat(*db_, "SELECT * FROM FEAT_BaseFeature LEFT JOIN FEAT_FeatureHandle ON id = feature_id ORDER BY id ASC");
    // prepare sub-queries (optional - corresponding tables may not be present):
    SQLite::Statement query_meta(*db_, "");
    SQLite::Statement query_match(*db_, "");
    prepareQueriesBaseFeature_(query_meta, query_match);
    SQLite::Statement query_ratio(*db_, "");
    if (db_->tableExists("FEAT_ConsensusRatio"))
    {
      query_ratio = SQLite::Statement(*db_, "SELECT * FROM FEAT_ConsensusRatio WHERE feature_id = :id " \
                                      "ORDER BY ratio_index DESC");
    }

    while (query_feat.executeStep())
    {
      if (query_feat.getColumn("subordinate_of").isNull()) // ConsensusFeature
      {
        int id = query_feat.getColumn("id").getInt();
        ConsensusFeature feature(makeBaseFeature_(id, query_feat, query_meta, query_match));
        consensus.push_back(feature);
        if (!isEmpty_(query_ratio))
        {
          query_ratio.bind(":id", id);
          while (query_ratio.executeStep())
          {
            Size ratio_index = query_ratio.getColumn("ratio_index").getUInt();
            // first row should have max. hull index (sorted descending):
            if (feature.getRatios().size() <= ratio_index)
            {
              feature.getRatios().resize(ratio_index + 1);
            }
            ConsensusFeature::Ratio& ratio = feature.getRatios()[ratio_index];
            ratio.ratio_value_ = query_ratio.getColumn("ratio_value").getDouble();
            ratio.denominator_ref_ = query_ratio.getColumn("denominator_ref").getString();
            ratio.numerator_ref_ = query_ratio.getColumn("numerator_ref").getString();
            ratio.description_ = ListUtils::create<std::string>(query_ratio.getColumn("description").getString());
          }
          query_ratio.reset(); // get ready for new executeStep()
        }
      }
      else // FeatureHandle
      {
        BaseFeature feature(makeBaseFeature_(-1, query_feat, query_meta, query_match));
        UInt64 map_index = query_feat.getColumn("map_index").getInt64();
        FeatureHandle handle(map_index, feature);
        consensus.back().insert(handle);
      }
    }
  }


  void OMSFileLoad::loadConsensusColumnHeaders_(ConsensusMap& consensus)
  {
    consensus.getColumnHeaders().clear();
    if (!db_->tableExists("FEAT_ConsensusColumnHeader")) return;

    SQLite::Statement query(*db_, "SELECT * FROM FEAT_ConsensusColumnHeader");
    SQLite::Statement query_info(*db_, "");
    bool have_meta_info = prepareQueryMetaInfo_(query_info, "FEAT_ConsensusColumnHeader");
    while (query.executeStep())
    {
      UInt64 id = query.getColumn("id").getInt64();
      ConsensusMap::ColumnHeader header;
      header.filename = query.getColumn("filename").getString();
      header.label = query.getColumn("label").getString();
      header.size = query.getColumn("size").getInt64();
      header.unique_id = query.getColumn("unique_id").getInt64();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(query_info, header, id);
      }
      consensus.getColumnHeaders()[id] = header;
    }
  }


  void OMSFileLoad::load(ConsensusMap& consensus)
  {
    load(consensus.getIdentificationData()); // load IDs, if any
    startProgress(0, 4, "Reading feature data from file");
    loadMapMetaData_(consensus);
    nextProgress();
    loadConsensusColumnHeaders_(consensus);
    nextProgress();
    loadDataProcessing_(consensus.getDataProcessing());
    nextProgress();
    loadConsensusFeatures_(consensus);
    consensus.updateRanges();
    endProgress();
  }


  // file-local helper (was private method, moved here to avoid nlohmann include in header)
  static nlohmann::json exportTableToJSON_(SQLite::Database& db, const std::string& table, const std::string& order_by)
  {
    using json = nlohmann::json;
    // code based on: https://stackoverflow.com/a/18067555
    std::string sql = "SELECT * FROM " + table;
    if (!order_by.empty())
    {
      sql += " ORDER BY " + order_by;
    }

    SQLite::Statement query(db, sql);

    json array = json::array();
    while (query.executeStep())
    {
      json record = json::object();
      for (int i = 0; i < query.getColumnCount(); ++i)
      {
        // @TODO: this will repeat field names for every row -
        // avoid this with separate "header" and "rows" (array)?

        // sqlite stores each cell based on the actual value, not the declared column type;
        // thus, we could use query.getColumnDeclaredType(i), but it would incur conversion
        switch (query.getColumn(i).getType())
        {
          case SQLITE_INTEGER: record[query.getColumnName(i)] = query.getColumn(i).getInt64(); break;
          case SQLITE_FLOAT: record[query.getColumnName(i)] = query.getColumn(i).getDouble(); break;
          case SQLITE_BLOB: {
            const auto column = query.getColumn(i);
            const auto* bytes = static_cast<const unsigned char*>(column.getBlob());
            std::string hex;
            hex.reserve(column.getBytes() * 2);
            constexpr char digits[] = "0123456789abcdef";
            for (int b = 0; b < column.getBytes(); ++b)
            {
              hex += digits[bytes[b] >> 4];
              hex += digits[bytes[b] & 15];
            }
            record[query.getColumnName(i)] = hex;
            break;
          }
          case SQLITE_NULL: record[query.getColumnName(i)] = ""; break;
          case SQLITE3_TEXT: record[query.getColumnName(i)] = query.getColumn(i).getText(); break;
          default:
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
      }
      array.push_back(record);
    }
    return array;
  }


  void OMSFileLoad::exportToJSON(ostream& output)
  {
    using json = nlohmann::json;
    // @TODO: this constructs the whole JSON file in memory - write directly to stream instead?
    // (more code, but would use less memory)
    json json_data = json::object();
    // get names of all tables (except SQLite-internal ones) in the database:
    SQLite::Statement query(*db_, "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name");
    while (query.executeStep())
    {
      std::string table = query.getColumn("name").getString();
      std::string order_by = "id"; // row order for most tables
      // special cases regarding ordering, e.g. tables without "id" column:
      if (StringUtils::hasSuffix(table, "_MetaInfo"))
      {
        order_by = "parent_id, name";
      }
      else if (StringUtils::hasSuffix(table, "_AppliedProcessingStep"))
      {
        order_by = "parent_id, processing_step_order, score_type_id";
      }
      else if (auto pos = export_order_by_.find(table); pos != export_order_by_.end())
      {
        order_by = pos->second;
      }
      json_data[table] = exportTableToJSON_(*db_, table, order_by);
    }

    output << json_data.dump(4) << '\n';
  }
}
