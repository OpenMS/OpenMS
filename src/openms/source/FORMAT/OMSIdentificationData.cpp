// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
#include "OMSIdentificationData.h"

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <SQLiteCpp/SQLiteCpp.h>
#include <bit>
#include <charconv>
#include <cmath>
#include <limits>
#include <utility>
#include <variant>

namespace OpenMS::Internal
{
namespace
{
  using ID = IdentificationData;
  using Key = int64_t;
  using Row = SQLite::Statement;

  [[noreturn]] void invalid(const std::string& message)
  { throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid OMS identification data", message); }

  // SQLite normalizes negative zero and converts NaNs to NULL. Ordinary numbers
  // remain SQL REALs; an auxiliary INTEGER is populated only for exceptional bits.
  struct Real
  { std::optional<double> value; };
  using Value = std::variant<std::nullptr_t, Key, std::string, Real>;
  template<class T>
  Key integer(T value)
  { return static_cast<Key>(value); }
  Value optionalText(const std::optional<std::string>& value)
  { return value ? Value(*value) : Value(nullptr); }
  Value unsignedText(const std::optional<UInt64>& value)
  { return value ? Value(std::to_string(*value)) : Value(nullptr); }
  UInt64 unsignedValue(const std::string& text)
  {
    UInt64 value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc {} || end != text.data() + text.size()) invalid("Invalid unsigned integer: " + text);
    return value;
  }
  std::string str(const Row& row, const char* name)
  { return row.getColumn(name).getString(); }
  Key number(const Row& row, const char* name)
  {
    const auto value = row.getColumn(name);
    if (! value.isInteger()) invalid(std::string("Expected an integer: ") + name);
    return value.getInt64();
  }
  template<class T>
  T checkedNumber(const Row& row, const char* name)
  {
    const auto value = number(row, name);
    if (! std::in_range<T>(value)) invalid(std::string("Integer out of range: ") + name);
    return static_cast<T>(value);
  }
  bool null(const Row& row, const char* name)
  { return row.getColumn(name).isNull(); }
  std::optional<std::string> optionalString(const Row& row, const char* name)
  { return null(row, name) ? std::nullopt : std::optional(str(row, name)); }
  std::optional<UInt64> optionalUnsigned(const Row& row, const char* name)
  { return null(row, name) ? std::nullopt : std::optional(unsignedValue(str(row, name))); }
  std::optional<double> optionalReal(const Row& row, const std::string& name)
  {
    const auto bits = row.getColumn((name + "_bits").c_str());
    if (! bits.isNull()) return std::bit_cast<double>(bits.getInt64());
    const auto value = row.getColumn(name.c_str());
    if (value.isNull()) return std::nullopt;
    return value.getDouble();
  }
  double real(const Row& row, const std::string& name)
  {
    const auto value = optionalReal(row, name);
    if (! value) invalid("Missing numeric value: " + name);
    return *value;
  }
  std::string dateText(const DateTime& date)
  { return date.isNull() ? "" : date.toString("yyyy-MM-ddThh:mm:ss.zzz"); }
  DateTime dateValue(const std::string& text)
  {
    DateTime date;
    if (! text.empty()) date.set(text);
    return date;
  }
  EmpiricalFormula formula(const Row& row, const char* composition, const char* charge)
  {
    EmpiricalFormula value(str(row, composition));
    value.setCharge(checkedNumber<Int>(row, charge));
    return value;
  }

  class Writer
  {
  public:
    explicit Writer(SQLite::Database& database): db(database)
    {
    }
    SQLite::Database& db;
    std::map<std::string, std::unique_ptr<Row>> inserts;
    std::map<std::string, Key> counters;
    std::map<std::string, Key> metadata_names;
    Key nameId(const std::string& name)
    {
      const auto found = metadata_names.find(name);
      if (found != metadata_names.end()) return found->second;
      const auto id = insert("ID_MetadataName", {name});
      metadata_names.emplace(name, id);
      return id;
    }
    void table(const std::string& name, const std::string& columns)
    { db.exec("CREATE TABLE " + name + " (id INTEGER PRIMARY KEY, " + columns + ")"); }
    Key insert(const std::string& table, const std::vector<Value>& values)
    {
      auto& statement = inserts[table];
      if (! statement)
      {
        std::string sql = "INSERT INTO " + table + " VALUES (?";
        for (const auto& value : values)
          sql += std::holds_alternative<Real>(value) ? ",?,?" : ",?";
        statement = std::make_unique<Row>(db, sql + ")");
      }
      const Key id = ++counters[table];
      statement->bind(1, id);
      int column = 2;
      for (const auto& value : values)
      {
        std::visit(
          [&](const auto& item) {
            using T = std::decay_t<decltype(item)>;
            if constexpr (std::is_same_v<T, Real>)
            {
              if (item.value)
              {
                statement->bind(column++, *item.value);
                if (! std::isfinite(*item.value) || (*item.value == 0 && std::signbit(*item.value)))
                  statement->bind(column++, std::bit_cast<Key>(*item.value));
                else
                  statement->bind(column++);
              }
              else
              {
                statement->bind(column++);
                statement->bind(column++);
              }
            }
            else if constexpr (std::is_same_v<T, std::nullptr_t>)
              statement->bind(column++);
            else
              statement->bind(column++, item);
          },
          value);
      }
      statement->exec();
      statement->reset();
      return id;
    }
    void meta(const std::string& owner, Key parent, const MetaInfoInterface& info)
    {
      std::vector<std::string> keys;
      info.getKeys(keys);
      for (const auto& key : keys)
      {
        const auto& value = info.getMetaValue(key);
        Value text = nullptr, count = nullptr;
        Real floating {};
        switch (value.valueType())
        {
          case DataValue::STRING_VALUE:
            text = value.toString();
            break;
          case DataValue::INT_VALUE:
            count = static_cast<Key>(value);
            break;
          case DataValue::DOUBLE_VALUE:
            floating.value = static_cast<double>(value);
            break;
          default:
            break;
        }
        const auto id = insert("ID_Metadata", {nameId(owner), parent, nameId(key), integer(value.valueType()), integer(value.getUnitType()),
                                               integer(value.getUnit()), text, count, floating});
        switch (value.valueType())
        {
          case DataValue::STRING_LIST:
            for (const auto& item : value.toStringList())
              insert("ID_MetadataItem", {id, item, nullptr, Real {}});
            break;
          case DataValue::INT_LIST:
            for (const auto item : value.toIntList())
              insert("ID_MetadataItem", {id, nullptr, integer(item), Real {}});
            break;
          case DataValue::DOUBLE_LIST:
            for (const auto item : value.toDoubleList())
              insert("ID_MetadataItem", {id, nullptr, nullptr, Real {item}});
            break;
          default:
            break;
        }
      }
    }
  };

  // Ordered cursors merge PSMs and their child rows with their owners. Metadata
  // scalars share these scans; variable-length metadata lists use an indexed lookup.
  class Reader
  {
  public:
    explicit Reader(SQLite::Database& database): db(database)
    {
      Row names(db, "SELECT id,name FROM ID_MetadataName ORDER BY id");
      while (names.executeStep())
      {
        const auto id = number(names, "id");
        auto name = str(names, "name");
        metadata_name_ids.emplace(name, id);
        metadata_names.emplace(id, std::move(name));
      }
    }
    SQLite::Database& db;
    struct Cursor
    {
      Row row;
      bool available;
      Cursor(SQLite::Database& db, const std::string& sql): row(db, sql), available(row.executeStep())
      {
      }
    };
    std::map<std::string, std::unique_ptr<Cursor>> cursors;
    std::map<Key, std::string> metadata_names;
    std::map<std::string, Key> metadata_name_ids;
    Key metadata_read = 0;
    std::unique_ptr<Row> metadata_items;
    template<class F>
    void children(const std::string& table, Key parent, F consume)
    {
      auto& cursor = cursors[table];
      if (! cursor) cursor = std::make_unique<Cursor>(db, "SELECT * FROM " + table + " ORDER BY parent_id,id");
      while (cursor->available && number(cursor->row, "parent_id") <= parent)
      {
        if (number(cursor->row, "parent_id") != parent) invalid("Orphan or unordered rows in " + table);
        consume(cursor->row);
        cursor->available = cursor->row.executeStep();
      }
    }
    void meta(const std::string& owner, Key parent, MetaInfoInterface& info)
    {
      auto& cursor = cursors["metadata:" + owner];
      if (! cursor)
      {
        const auto found = metadata_name_ids.find(owner);
        const auto owner_id = found == metadata_name_ids.end() ? 0 : found->second;
        cursor = std::make_unique<Cursor>(db, "SELECT * FROM ID_Metadata WHERE owner_id=" + std::to_string(owner_id) + " ORDER BY parent_id,name_id");
      }
      while (cursor->available && number(cursor->row, "parent_id") <= parent)
      {
        auto& row = cursor->row;
        if (number(row, "parent_id") != parent) invalid("Orphan metadata for " + owner);
        const auto type = number(row, "type");
        DataValue value;
        switch (type)
        {
          case DataValue::STRING_VALUE:
            value = str(row, "text_value");
            break;
          case DataValue::INT_VALUE:
            value = number(row, "integer_value");
            break;
          case DataValue::DOUBLE_VALUE:
            value = real(row, "real_value");
            break;
          case DataValue::EMPTY_VALUE:
            break;
          case DataValue::STRING_LIST:
          case DataValue::INT_LIST:
          case DataValue::DOUBLE_LIST: {
            // Different owner types are consumed independently. Reuse one indexed
            // lookup for their lists; no statement is prepared per list value.
            if (! metadata_items) metadata_items = std::make_unique<Row>(db, "SELECT * FROM ID_MetadataItem WHERE parent_id=? ORDER BY id");
            auto& items = *metadata_items;
            items.bind(1, number(row, "id"));
            StringList strings;
            IntList integers;
            DoubleList reals;
            while (items.executeStep())
            {
              if (type == DataValue::STRING_LIST) strings.push_back(str(items, "text_value"));
              else if (type == DataValue::INT_LIST)
                integers.push_back(checkedNumber<Int>(items, "integer_value"));
              else
                reals.push_back(real(items, "real_value"));
            }
            items.reset();
            if (type == DataValue::STRING_LIST) value = strings;
            else if (type == DataValue::INT_LIST)
              value = integers;
            else
              value = reals;
            break;
          }
          default:
            invalid("Invalid metadata type");
        }
        const auto unit = number(row, "unit_type");
        if (unit < 0 || unit > DataValue::OTHER) invalid("Invalid metadata unit type");
        value.setUnitType(static_cast<DataValue::UnitType>(unit));
        value.setUnit(checkedNumber<Int>(row, "unit"));
        info.setMetaValue(metadata_names.at(number(row, "name_id")), value);
        ++metadata_read;
        cursor->available = row.executeStep();
      }
    }
    void finish()
    {
      if (db.execAndGet("SELECT COUNT(*) FROM ID_Metadata").getInt64() != metadata_read) invalid("Unowned metadata rows");
      for (const auto& [table, cursor] : cursors)
        if (cursor->available) invalid("Unconsumed rows in " + table);
    }
  };

  void schema(Writer& out, Size scores)
  {
    // Intern names once, including owner names. Repeating both strings in every
    // metadata row and index otherwise dominates a million-PSM SQLite database.
    out.table("ID_MetadataName", "name TEXT NOT NULL UNIQUE");
    out.table("ID_Metadata", "owner_id INTEGER NOT NULL REFERENCES ID_MetadataName(id),parent_id INTEGER NOT NULL,name_id INTEGER NOT NULL "
                             "REFERENCES ID_MetadataName(id),type INTEGER NOT NULL,unit_type INTEGER NOT NULL,unit INTEGER NOT NULL,text_value "
                             "TEXT,integer_value INTEGER,real_value REAL,real_value_bits INTEGER,UNIQUE(owner_id,parent_id,name_id)");
    out.table("ID_MetadataItem",
              "parent_id INTEGER NOT NULL REFERENCES ID_Metadata(id),text_value TEXT,integer_value INTEGER,real_value REAL,real_value_bits INTEGER");
    out.table("ID_Score", "name TEXT NOT NULL,accession TEXT NOT NULL,higher_better INTEGER NOT NULL,scope INTEGER NOT NULL,software TEXT NOT "
                          "NULL,software_version TEXT NOT NULL,calibration TEXT NOT NULL,aggregation TEXT NOT NULL");
    out.table("ID_PSMScore", "definition_id INTEGER NOT NULL REFERENCES ID_Score(id)");
    out.table(
      "ID_Processing",
      "identifier TEXT NOT NULL,search_engine TEXT NOT NULL,search_engine_version TEXT NOT NULL,date TEXT NOT NULL,score_type TEXT NOT "
      "NULL,higher_better INTEGER NOT NULL,significance_threshold REAL,significance_threshold_bits INTEGER,db TEXT NOT NULL,db_version TEXT NOT "
      "NULL,taxonomy TEXT NOT NULL,charges TEXT NOT NULL,mass_type INTEGER NOT NULL,missed_cleavages INTEGER NOT NULL,fragment_tolerance "
      "REAL,fragment_tolerance_bits INTEGER,fragment_ppm INTEGER NOT NULL,precursor_tolerance REAL,precursor_tolerance_bits INTEGER,precursor_ppm "
      "INTEGER NOT NULL,specificity INTEGER NOT NULL,enzyme_name TEXT NOT NULL,enzyme_regex TEXT NOT NULL,enzyme_description TEXT NOT "
      "NULL,n_term_gain TEXT NOT NULL,n_term_charge INTEGER NOT NULL,c_term_gain TEXT NOT NULL,c_term_charge INTEGER NOT NULL,psi_id TEXT NOT "
      "NULL,xtandem_id TEXT NOT NULL,comet_id INTEGER NOT NULL,msgf_id INTEGER NOT NULL,omssa_id INTEGER NOT NULL");
    out.table("ID_ProcessingString", "parent_id INTEGER NOT NULL REFERENCES ID_Processing(id),kind INTEGER NOT NULL,value TEXT NOT NULL");
    out.table("ID_Run", "identifier TEXT NOT NULL UNIQUE,uuid TEXT NOT NULL UNIQUE,kind INTEGER NOT NULL,processing_id INTEGER NOT NULL REFERENCES "
                        "ID_Processing(id),has_scores INTEGER NOT NULL,primary_score INTEGER REFERENCES ID_PSMScore(id),has_parents INTEGER NOT "
                        "NULL,next_query TEXT NOT NULL,next_match TEXT NOT NULL");
    out.table("ID_Source", "parent_id INTEGER NOT NULL REFERENCES ID_Run(id),identifier TEXT NOT NULL,path TEXT NOT NULL");
    out.table("ID_SourceFile", "parent_id INTEGER NOT NULL REFERENCES ID_Source(id),path TEXT NOT NULL");
    out.table("ID_Query", "parent_id INTEGER NOT NULL REFERENCES ID_Source(id),query_id TEXT NOT NULL,data_id TEXT NOT NULL,rt REAL,rt_bits "
                          "INTEGER,mz REAL,mz_bits INTEGER,selected_match TEXT");
    std::string match
      = "parent_id INTEGER NOT NULL REFERENCES ID_Query(id),match_id TEXT NOT NULL,encoding INTEGER NOT NULL,representation TEXT NOT NULL,charge "
        "INTEGER NOT NULL,calculated_mz REAL,calculated_mz_bits INTEGER,target_decoy INTEGER NOT NULL,name TEXT NOT NULL,formula TEXT,adduct_name "
        "TEXT,adduct_formula TEXT,adduct_formula_charge INTEGER,adduct_charge INTEGER,adduct_multiplier INTEGER";
    for (Size i = 0; i < scores; ++i)
      match += ",score_" + std::to_string(i) + " REAL,score_" + std::to_string(i) + "_bits INTEGER";
    out.table("ID_Match", match);
    out.table("ID_MatchIdentifier", "parent_id INTEGER NOT NULL REFERENCES ID_Match(id),database TEXT NOT NULL,accession TEXT NOT NULL");
    out.table("ID_Evidence", "parent_id INTEGER NOT NULL REFERENCES ID_Match(id),database TEXT NOT NULL,accession TEXT NOT NULL,start TEXT,end "
                             "TEXT,before TEXT NOT NULL,after TEXT NOT NULL");
    out.table("ID_PeakAnnotation", "parent_id INTEGER NOT NULL REFERENCES ID_Match(id),annotation TEXT NOT NULL,charge INTEGER NOT NULL,mz "
                                   "REAL,mz_bits INTEGER,intensity REAL,intensity_bits INTEGER");
    out.table("ID_Parent", "parent_id INTEGER NOT NULL REFERENCES ID_Run(id),database TEXT NOT NULL,accession TEXT NOT NULL,target_decoy INTEGER NOT "
                           "NULL,sequence TEXT NOT NULL,description TEXT NOT NULL,UNIQUE(parent_id,database,accession)");
    out.table("ID_Inference", "identifier TEXT NOT NULL UNIQUE,processing_id INTEGER NOT NULL REFERENCES ID_Processing(id),parent_score_id INTEGER "
                              "REFERENCES ID_Score(id),group_score_id INTEGER REFERENCES ID_Score(id)");
    out.table("ID_Input", "parent_id INTEGER NOT NULL REFERENCES ID_Inference(id),run_identifier TEXT NOT NULL,run_uuid TEXT NOT NULL,score_id "
                          "INTEGER REFERENCES ID_Score(id),selection TEXT NOT NULL");
    out.table("ID_Accession", "parent_id INTEGER NOT NULL REFERENCES ID_Inference(id),alias TEXT NOT NULL,database TEXT NOT NULL,accession TEXT NOT "
                              "NULL,UNIQUE(parent_id,alias)");
    out.table("ID_Protein", "parent_id INTEGER NOT NULL REFERENCES ID_Processing(id),accession TEXT NOT NULL,score REAL,score_bits INTEGER,rank "
                            "INTEGER NOT NULL,sequence TEXT NOT NULL,coverage REAL,coverage_bits INTEGER");
    out.table("ID_ProteinGroup",
              "parent_id INTEGER NOT NULL REFERENCES ID_Processing(id),kind INTEGER NOT NULL,probability REAL,probability_bits INTEGER");
    out.table("ID_GroupMember", "parent_id INTEGER NOT NULL REFERENCES ID_ProteinGroup(id),accession TEXT NOT NULL");
    out.table(
      "ID_Modification",
      "parent_id INTEGER NOT NULL REFERENCES ID_Protein(id),position TEXT NOT NULL,mod_id TEXT NOT NULL,full_id TEXT NOT NULL,psi_mod_accession TEXT "
      "NOT NULL,unimod_record_id INTEGER NOT NULL,full_name TEXT NOT NULL,name TEXT NOT NULL,term_specificity INTEGER NOT NULL,origin INTEGER NOT "
      "NULL,classification INTEGER NOT NULL,provenance INTEGER NOT NULL,average_mass REAL,average_mass_bits INTEGER,mono_mass REAL,mono_mass_bits "
      "INTEGER,diff_average_mass REAL,diff_average_mass_bits INTEGER,diff_mono_mass REAL,diff_mono_mass_bits INTEGER,formula TEXT NOT "
      "NULL,diff_formula TEXT NOT NULL,diff_formula_charge INTEGER NOT NULL");
    out.table("ID_ModificationItem", "parent_id INTEGER NOT NULL REFERENCES ID_Modification(id),kind INTEGER NOT NULL,text_value TEXT,charge "
                                     "INTEGER,real_value REAL,real_value_bits INTEGER");
    out.table("ID_GroupArray", "parent_id INTEGER NOT NULL REFERENCES ID_ProteinGroup(id),kind INTEGER NOT NULL,name TEXT NOT NULL");
    out.table(
      "ID_ArrayItem",
      "parent_id INTEGER NOT NULL REFERENCES ID_GroupArray(id),text_value TEXT,integer_value INTEGER,real_value REAL,real_value_bits INTEGER");
    out.table("ID_ArrayProcessing",
              "parent_id INTEGER NOT NULL REFERENCES ID_GroupArray(id),is_null INTEGER NOT NULL,software TEXT,version TEXT,date TEXT");
    out.table("ID_ArrayAction", "parent_id INTEGER NOT NULL REFERENCES ID_ArrayProcessing(id),action INTEGER NOT NULL");
  }

  Key writeScore(Writer& out, const ID::ScoreDefinition& score)
  {
    auto id = out.insert("ID_Score", {score.name, score.accession, integer(score.higher_better), integer(score.scope), score.software,
                                      score.software_version, score.calibration, score.aggregation});
    out.meta("ID_Score", id, score.parameters);
    return id;
  }
  ID::ScoreDefinition readScore(Reader& in, const Row& row)
  {
    ID::ScoreDefinition score;
    score.name = str(row, "name");
    score.accession = str(row, "accession");
    score.higher_better = number(row, "higher_better");
    score.scope = static_cast<ID::ScoreScope>(number(row, "scope"));
    score.software = str(row, "software");
    score.software_version = str(row, "software_version");
    score.calibration = str(row, "calibration");
    score.aggregation = str(row, "aggregation");
    in.meta("ID_Score", number(row, "id"), score.parameters);
    return score;
  }

  void writeModifications(Writer& out, Key protein, const ProteinHit& hit)
  {
    for (const auto& [position, mod] : hit.getModifications())
    {
      const auto id
        = out.insert("ID_Modification",
                     {protein, std::to_string(position), mod.getId(), mod.getFullId(), mod.getPSIMODAccession(), integer(mod.getUniModRecordId()),
                      mod.getFullName(), mod.getName(), integer(mod.getTermSpecificity()), integer(static_cast<unsigned char>(mod.getOrigin())),
                      integer(mod.getSourceClassification()), integer(mod.getProvenance()), Real {mod.getAverageMass()}, Real {mod.getMonoMass()},
                      Real {mod.getDiffAverageMass()}, Real {mod.getDiffMonoMass()}, mod.getFormula(), mod.getDiffFormula().toString(),
                      integer(mod.getDiffFormula().getCharge())});
      for (const auto& item : mod.getSynonyms())
        out.insert("ID_ModificationItem", {id, Key(0), item, nullptr, Real {}});
      for (const auto& item : mod.getNeutralLossDiffFormulas())
        out.insert("ID_ModificationItem", {id, Key(1), item.toString(), integer(item.getCharge()), Real {}});
      for (double item : mod.getNeutralLossMonoMasses())
        out.insert("ID_ModificationItem", {id, Key(2), nullptr, nullptr, Real {item}});
      for (double item : mod.getNeutralLossAverageMasses())
        out.insert("ID_ModificationItem", {id, Key(3), nullptr, nullptr, Real {item}});
    }
  }
  void readModifications(Reader& in, Key protein, ProteinHit& hit)
  {
    std::set<std::pair<Size, ResidueModification>> mods;
    in.children("ID_Modification", protein, [&](const Row& row) {
      ResidueModification mod;
      mod.setId(str(row, "mod_id"));
      mod.setPSIMODAccession(str(row, "psi_mod_accession"));
      mod.setUniModRecordId(checkedNumber<Int>(row, "unimod_record_id"));
      mod.setFullName(str(row, "full_name"));
      mod.setName(str(row, "name"));
      const auto specificity = number(row, "term_specificity"), classification = number(row, "classification"),
                 provenance = number(row, "provenance");
      if (specificity < 0 || specificity >= ResidueModification::NUMBER_OF_TERM_SPECIFICITY || classification < 0
          || classification >= ResidueModification::NUMBER_OF_SOURCE_CLASSIFICATIONS || provenance < 0
          || provenance >= ResidueModification::NUMBER_OF_PROVENANCE)
        invalid("Invalid modification enum");
      mod.setTermSpecificity(static_cast<ResidueModification::TermSpecificity>(specificity));
      mod.setOrigin(static_cast<char>(number(row, "origin")));
      mod.setSourceClassification(static_cast<ResidueModification::SourceClassification>(classification));
      mod.setProvenance(static_cast<ResidueModification::Provenance>(provenance));
      mod.setAverageMass(real(row, "average_mass"));
      mod.setMonoMass(real(row, "mono_mass"));
      mod.setDiffAverageMass(real(row, "diff_average_mass"));
      mod.setDiffMonoMass(real(row, "diff_mono_mass"));
      mod.setFormula(str(row, "formula"));
      mod.setDiffFormula(formula(row, "diff_formula", "diff_formula_charge"));
      std::set<std::string> synonyms;
      std::vector<EmpiricalFormula> formulas;
      std::vector<double> mono, average;
      in.children("ID_ModificationItem", number(row, "id"), [&](const Row& item) {
        switch (number(item, "kind"))
        {
          case 0:
            if (! synonyms.insert(str(item, "text_value")).second) invalid("Duplicate modification synonym");
            break;
          case 1:
            formulas.push_back(formula(item, "text_value", "charge"));
            break;
          case 2:
            mono.push_back(real(item, "real_value"));
            break;
          case 3:
            average.push_back(real(item, "real_value"));
            break;
          default:
            invalid("Invalid modification item kind");
        }
      });
      mod.setSynonyms(synonyms);
      mod.setNeutralLossDiffFormulas(formulas);
      mod.setNeutralLossMonoMasses(mono);
      mod.setNeutralLossAverageMasses(average);
      if (! str(row, "full_id").empty()) mod.setFullId(str(row, "full_id"));
      if (! mods.emplace(unsignedValue(str(row, "position")), std::move(mod)).second) invalid("Duplicate protein modification");
    });
    hit.setModifications(mods);
  }

  template<class Arrays>
  void writeArrays(Writer& out, Key group, Key kind, const Arrays& arrays)
  {
    for (const auto& array : arrays)
    {
      const auto id = out.insert("ID_GroupArray", {group, kind, array.getName()});
      out.meta("ID_GroupArray", id, array);
      for (const auto& value : array)
      {
        if constexpr (std::is_same_v<typename Arrays::value_type::value_type, std::string>) out.insert("ID_ArrayItem", {id, value, nullptr, Real {}});
        else if constexpr (std::is_floating_point_v<typename Arrays::value_type::value_type>)
        {
          // Float NaN payloads must not be changed by widening to double.
          out.insert("ID_ArrayItem", {id, nullptr, integer(std::bit_cast<UInt32>(value)), Real {static_cast<double>(value)}});
        }
        else
          out.insert("ID_ArrayItem", {id, nullptr, integer(value), Real {}});
      }
      for (const auto& processing : array.getDataProcessing())
      {
        if (! processing)
        {
          out.insert("ID_ArrayProcessing", {id, Key(1), nullptr, nullptr, nullptr});
          continue;
        }
        const auto& software = processing->getSoftware();
        const auto pid
          = out.insert("ID_ArrayProcessing", {id, Key(0), software.getName(), software.getVersion(), dateText(processing->getCompletionTime())});
        out.meta("ID_ArrayProcessing", pid, *processing);
        out.meta("ID_ArraySoftware", pid, software);
        for (auto action : processing->getProcessingActions())
          out.insert("ID_ArrayAction", {pid, integer(action)});
      }
    }
  }
  template<class Arrays>
  void readArray(Reader& in, const Row& row, Arrays& arrays)
  {
    typename Arrays::value_type array;
    const auto id = number(row, "id");
    array.setName(str(row, "name"));
    in.meta("ID_GroupArray", id, array);
    in.children("ID_ArrayItem", id, [&](const Row& item) {
      if constexpr (std::is_same_v<typename Arrays::value_type::value_type, std::string>) array.push_back(str(item, "text_value"));
      else if constexpr (std::is_floating_point_v<typename Arrays::value_type::value_type>)
        array.push_back(std::bit_cast<float>(checkedNumber<UInt32>(item, "integer_value")));
      else
        array.push_back(checkedNumber<Int>(item, "integer_value"));
    });
    std::vector<DataProcessingPtr> processing;
    in.children("ID_ArrayProcessing", id, [&](const Row& item) {
      if (number(item, "is_null"))
      {
        processing.push_back(nullptr);
        return;
      }
      auto value = std::make_shared<DataProcessing>();
      auto& software = value->getSoftware();
      software.setName(str(item, "software"));
      software.setVersion(str(item, "version"));
      value->setCompletionTime(dateValue(str(item, "date")));
      in.meta("ID_ArrayProcessing", number(item, "id"), *value);
      in.meta("ID_ArraySoftware", number(item, "id"), software);
      std::set<DataProcessing::ProcessingAction> actions;
      in.children("ID_ArrayAction", number(item, "id"), [&](const Row& action) {
        const auto kind = number(action, "action");
        if (kind < 0 || kind >= DataProcessing::SIZE_OF_PROCESSINGACTION
            || ! actions.insert(static_cast<DataProcessing::ProcessingAction>(kind)).second)
          invalid("Invalid array processing action");
      });
      value->setProcessingActions(actions);
      processing.push_back(std::move(value));
    });
    array.setDataProcessing(processing);
    arrays.push_back(std::move(array));
  }

  Key writeProcessing(Writer& out, const ProteinIdentification& processing)
  {
    const auto& sp = processing.getSearchParameters();
    const auto& enzyme = sp.digestion_enzyme;
    const auto id = out.insert("ID_Processing", {processing.getIdentifier(),
                                                 processing.getSearchEngine(),
                                                 processing.getSearchEngineVersion(),
                                                 dateText(processing.getDateTime()),
                                                 processing.getScoreType(),
                                                 integer(processing.isHigherScoreBetter()),
                                                 Real {processing.getSignificanceThreshold()},
                                                 sp.db,
                                                 sp.db_version,
                                                 sp.taxonomy,
                                                 sp.charges,
                                                 integer(sp.mass_type),
                                                 integer(sp.missed_cleavages),
                                                 Real {sp.fragment_mass_tolerance},
                                                 integer(sp.fragment_mass_tolerance_ppm),
                                                 Real {sp.precursor_mass_tolerance},
                                                 integer(sp.precursor_mass_tolerance_ppm),
                                                 integer(sp.enzyme_term_specificity),
                                                 enzyme.getName(),
                                                 enzyme.getRegEx(),
                                                 enzyme.getRegExDescription(),
                                                 enzyme.getNTermGain().toString(),
                                                 integer(enzyme.getNTermGain().getCharge()),
                                                 enzyme.getCTermGain().toString(),
                                                 integer(enzyme.getCTermGain().getCharge()),
                                                 enzyme.getPSIID(),
                                                 enzyme.getXTandemID(),
                                                 integer(enzyme.getCometID()),
                                                 integer(enzyme.getMSGFID()),
                                                 integer(enzyme.getOMSSAID())});
    out.meta("ID_Processing", id, processing);
    out.meta("ID_SearchParameters", id, sp);
    for (const auto& value : sp.fixed_modifications)
      out.insert("ID_ProcessingString", {id, Key(0), value});
    for (const auto& value : sp.variable_modifications)
      out.insert("ID_ProcessingString", {id, Key(1), value});
    for (const auto& value : enzyme.getSynonyms())
      out.insert("ID_ProcessingString", {id, Key(2), value});
    for (const auto& protein : processing.getHits())
    {
      const auto pid = out.insert("ID_Protein", {id, protein.getAccession(), Real {protein.getScore()}, integer(protein.getRank()),
                                                 protein.getSequence(), Real {protein.getCoverage()}});
      out.meta("ID_Protein", pid, protein);
      writeModifications(out, pid, protein);
    }
    const auto groups = [&](const auto& values, Key kind) {
      for (const auto& group : values)
      {
        const auto gid = out.insert("ID_ProteinGroup", {id, kind, Real {group.probability}});
        for (const auto& accession : group.accessions)
          out.insert("ID_GroupMember", {gid, accession});
        writeArrays(out, gid, 0, group.getFloatDataArrays());
        writeArrays(out, gid, 1, group.getIntegerDataArrays());
        writeArrays(out, gid, 2, group.getStringDataArrays());
      }
    };
    groups(processing.getProteinGroups(), 0);
    groups(processing.getIndistinguishableProteins(), 1);
    return id;
  }

  ProteinIdentification readProcessing(Reader& in, const Row& row)
  {
    ProteinIdentification processing;
    const auto id = number(row, "id");
    processing.setIdentifier(str(row, "identifier"));
    processing.setSearchEngine(str(row, "search_engine"));
    processing.setSearchEngineVersion(str(row, "search_engine_version"));
    processing.setDateTime(dateValue(str(row, "date")));
    processing.setScoreType(str(row, "score_type"));
    processing.setHigherScoreBetter(number(row, "higher_better"));
    processing.setSignificanceThreshold(real(row, "significance_threshold"));
    auto& sp = processing.getSearchParameters();
    sp.db = str(row, "db");
    sp.db_version = str(row, "db_version");
    sp.taxonomy = str(row, "taxonomy");
    sp.charges = str(row, "charges");
    const auto mass_type = number(row, "mass_type"), specificity = number(row, "specificity");
    if (mass_type < 0 || mass_type >= static_cast<Key>(ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE)
        || (specificity != 8 && specificity != 9 && (specificity < 0 || specificity > 3)))
      invalid("Invalid search parameter enum");
    sp.mass_type = static_cast<ProteinIdentification::PeakMassType>(mass_type);
    sp.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(specificity);
    sp.missed_cleavages = checkedNumber<UInt>(row, "missed_cleavages");
    sp.fragment_mass_tolerance = real(row, "fragment_tolerance");
    sp.fragment_mass_tolerance_ppm = number(row, "fragment_ppm");
    sp.precursor_mass_tolerance = real(row, "precursor_tolerance");
    sp.precursor_mass_tolerance_ppm = number(row, "precursor_ppm");
    std::set<std::string> synonyms;
    in.children("ID_ProcessingString", id, [&](const Row& item) {
      switch (number(item, "kind"))
      {
        case 0:
          sp.fixed_modifications.push_back(str(item, "value"));
          break;
        case 1:
          sp.variable_modifications.push_back(str(item, "value"));
          break;
        case 2:
          if (! synonyms.insert(str(item, "value")).second) invalid("Duplicate enzyme synonym");
          break;
        default:
          invalid("Invalid search parameter list kind");
      }
    });
    sp.digestion_enzyme = Protease(str(row, "enzyme_name"), str(row, "enzyme_regex"), synonyms, str(row, "enzyme_description"),
                                   formula(row, "n_term_gain", "n_term_charge"), formula(row, "c_term_gain", "c_term_charge"), str(row, "psi_id"),
                                   str(row, "xtandem_id"), checkedNumber<Int>(row, "comet_id"), checkedNumber<Int>(row, "msgf_id"),
                                   checkedNumber<Int>(row, "omssa_id"));
    in.meta("ID_Processing", id, processing);
    in.meta("ID_SearchParameters", id, sp);
    in.children("ID_Protein", id, [&](const Row& protein) {
      ProteinHit hit(real(protein, "score"), checkedNumber<UInt>(protein, "rank"), str(protein, "accession"), str(protein, "sequence"));
      hit.setCoverage(real(protein, "coverage"));
      in.meta("ID_Protein", number(protein, "id"), hit);
      readModifications(in, number(protein, "id"), hit);
      processing.insertHit(std::move(hit));
    });
    in.children("ID_ProteinGroup", id, [&](const Row& row) {
      ProteinIdentification::ProteinGroup group;
      const auto gid = number(row, "id"), kind = number(row, "kind");
      group.probability = real(row, "probability");
      in.children("ID_GroupMember", gid, [&](const Row& member) { group.accessions.push_back(str(member, "accession")); });
      in.children("ID_GroupArray", gid, [&](const Row& array) {
        switch (number(array, "kind"))
        {
          case 0:
            readArray(in, array, group.getFloatDataArrays());
            break;
          case 1:
            readArray(in, array, group.getIntegerDataArrays());
            break;
          case 2:
            readArray(in, array, group.getStringDataArrays());
            break;
          default:
            invalid("Invalid group array kind");
        }
      });
      if (kind == 0) processing.insertProteinGroup(group);
      else if (kind == 1)
        processing.insertIndistinguishableProteins(group);
      else
        invalid("Invalid protein group kind");
    });
    return processing;
  }

  void writeRun(Writer& out, const ID::Run& run)
  {
    const auto processing = writeProcessing(out, run.getProcessingMetadata());
    const auto id = out.insert(
      "ID_Run", {run.getIdentifier(), run.getUuid(), integer(run.getMoleculeKind()), processing, integer(! run.getScoreDefinitions().empty()),
                 run.getPrimaryScore() ? Value(integer(run.getPrimaryScore()->value) + 1) : Value(nullptr), integer(run.getParents().has_value()),
                 std::to_string(run.getNextQueryId()), std::to_string(run.getNextMatchId())});
    if (run.getParents())
      for (const auto& parent : *run.getParents())
      {
        const auto pid = out.insert(
          "ID_Parent", {id, parent.identity.database, parent.identity.accession, integer(parent.target_decoy), parent.sequence, parent.description});
        out.meta("ID_Parent", pid, parent);
      }
    for (const auto& block : run.getSourceBlocks())
    {
      const auto& source = block.source;
      const auto sid = out.insert("ID_Source", {id, source.identifier, source.path});
      out.meta("ID_Source", sid, source);
      for (const auto& path : source.primary_files)
        out.insert("ID_SourceFile", {sid, path});
      for (const auto& observation : block.identifications)
      {
        const auto qid
          = out.insert("ID_Query", {sid, std::to_string(observation.getId().value), observation.data_id, Real {observation.rt}, Real {observation.mz},
                                    observation.getSelectedMatch() ? Value(std::to_string(observation.getSelectedMatch()->value)) : Value(nullptr)});
        out.meta("ID_Query", qid, observation);
        for (const auto& match : observation.getMatches())
        {
          std::vector<Value> values {qid,
                                     std::to_string(match.getId().value),
                                     integer(match.encoding),
                                     match.representation,
                                     integer(match.charge),
                                     Real {match.calculated_mz},
                                     integer(match.target_decoy),
                                     match.name,
                                     optionalText(match.formula)};
          if (match.adduct)
          {
            const auto& a = *match.adduct;
            values.insert(values.end(), {a.getName(), a.getEmpiricalFormula().toString(), integer(a.getEmpiricalFormula().getCharge()),
                                         integer(a.getCharge()), integer(a.getMolMultiplier())});
          }
          else
            values.insert(values.end(), 5, nullptr);
          for (const auto score : match.getScoreValues())
            values.push_back(Real {std::isnan(score) ? std::nullopt : std::optional(score)});
          // Scoreless catalogs share the table but leave all score columns NULL.
          for (Size i = match.getScoreValues().size(); i < static_cast<Size>(out.counters["ID_PSMScore"]); ++i)
            values.push_back(Real {});
          const auto mid = out.insert("ID_Match", values);
          out.meta("ID_Match", mid, match);
          for (const auto& identity : match.identifiers)
            out.insert("ID_MatchIdentifier", {mid, identity.database, identity.accession});
          for (const auto& evidence : match.parent_evidence)
            out.insert("ID_Evidence", {mid, evidence.parent.database, evidence.parent.accession, unsignedText(evidence.start),
                                       unsignedText(evidence.end), evidence.before, evidence.after});
          for (const auto& peak : match.peak_annotations)
            out.insert("ID_PeakAnnotation", {mid, peak.annotation, integer(peak.charge), Real {peak.mz}, Real {peak.intensity}});
        }
      }
    }
  }

  ID::Run readRun(Reader& in, const Row& row, const std::vector<ID::ScoreDefinition>& scores, std::map<Key, ProteinIdentification>& processing)
  {
    ID::Run run(str(row, "identifier"), static_cast<ID::MoleculeKind>(number(row, "kind")));
    run.setProcessingMetadata(processing.at(number(row, "processing_id")));
    const auto id = number(row, "id");
    const bool scored = number(row, "has_scores");
    if (scored)
      for (const auto& definition : scores)
        run.addScore(definition);
    if (! null(row, "primary_score"))
    {
      const auto primary = number(row, "primary_score");
      if (! scored || primary < 1 || primary > static_cast<Key>(scores.size())) invalid("Invalid primary score");
      run.setPrimaryScore(run.getScoreId(static_cast<UInt32>(primary - 1)));
    }
    std::vector<ID::ParentRecord> parents;
    in.children("ID_Parent", id, [&](const Row& item) {
      ID::ParentRecord parent;
      parent.identity = {str(item, "database"), str(item, "accession")};
      parent.target_decoy = static_cast<ID::TargetDecoy>(number(item, "target_decoy"));
      parent.sequence = str(item, "sequence");
      parent.description = str(item, "description");
      in.meta("ID_Parent", number(item, "id"), parent);
      parents.push_back(std::move(parent));
    });
    if (number(row, "has_parents")) run.setParents(std::move(parents));
    else if (! parents.empty())
      invalid("Unexpected parent catalog");
    in.children("ID_Source", id, [&](const Row& item) {
      const auto sid = number(item, "id");
      ID::SourceFile source;
      source.identifier = str(item, "identifier");
      source.path = str(item, "path");
      in.meta("ID_Source", sid, source);
      in.children("ID_SourceFile", sid, [&](const Row& path) { source.primary_files.push_back(str(path, "path")); });
      const auto source_id = run.addSource(source);
      in.children("ID_Query", sid, [&](const Row& query) {
        const auto qid = number(query, "id");
        ID::Observation observation;
        observation.data_id = str(query, "data_id");
        observation.rt = optionalReal(query, "rt");
        observation.mz = optionalReal(query, "mz");
        in.meta("ID_Query", qid, observation);
        const auto query_id = run.importIdentification(source_id, ID::QueryId {unsignedValue(str(query, "query_id"))}, std::move(observation));
        in.children("ID_Match", qid, [&](const Row& row) {
          const auto mid = number(row, "id");
          ID::MatchData match;
          match.encoding = static_cast<ID::Encoding>(number(row, "encoding"));
          match.representation = str(row, "representation");
          match.charge = checkedNumber<Int>(row, "charge");
          match.calculated_mz = optionalReal(row, "calculated_mz");
          match.target_decoy = static_cast<ID::TargetDecoy>(number(row, "target_decoy"));
          match.name = str(row, "name");
          match.formula = optionalString(row, "formula");
          if (! null(row, "adduct_name"))
            match.adduct.emplace(str(row, "adduct_name"), formula(row, "adduct_formula", "adduct_formula_charge"),
                                 checkedNumber<Int>(row, "adduct_charge"), checkedNumber<UInt>(row, "adduct_multiplier"));
          in.meta("ID_Match", mid, match);
          in.children("ID_MatchIdentifier", mid,
                      [&](const Row& identity) { match.identifiers.push_back({str(identity, "database"), str(identity, "accession")}); });
          in.children("ID_Evidence", mid, [&](const Row& evidence) {
            match.parent_evidence.push_back({{str(evidence, "database"), str(evidence, "accession")},
                                             optionalUnsigned(evidence, "start"),
                                             optionalUnsigned(evidence, "end"),
                                             str(evidence, "before"),
                                             str(evidence, "after")});
          });
          in.children("ID_PeakAnnotation", mid, [&](const Row& row) {
            PeptideHit::PeakAnnotation peak;
            peak.annotation = str(row, "annotation");
            peak.charge = checkedNumber<Int>(row, "charge");
            peak.mz = real(row, "mz");
            peak.intensity = real(row, "intensity");
            match.peak_annotations.push_back(std::move(peak));
          });
          std::vector<std::optional<double>> values;
          for (Size i = 0; i < scores.size(); ++i)
          {
            const auto score = optionalReal(row, "score_" + std::to_string(i));
            if (scored) values.push_back(score);
            else if (score)
              invalid("Score in scoreless catalog");
          }
          run.importMatch(query_id, ID::MatchId {unsignedValue(str(row, "match_id"))}, std::move(match), values);
        });
        if (! null(query, "selected_match")) run.setSelectedMatch(query_id, ID::MatchId {unsignedValue(str(query, "selected_match"))});
      });
    });
    run.restoreIdentity(str(row, "uuid"), unsignedValue(str(row, "next_query")), unsignedValue(str(row, "next_match")));
    return run;
  }
} // namespace

void storeOMSIdentifications(SQLite::Database& db, const IdentificationData& data)
{
  data.validate();
  const std::vector<ID::ScoreDefinition>* scores = nullptr;
  for (const auto& run : data.getRuns())
    if (! run.getScoreDefinitions().empty())
    {
      scores = &run.getScoreDefinitions();
      break;
    }
  Writer out(db);
  schema(out, scores ? scores->size() : 0);
  if (scores)
    for (const auto& score : *scores)
      out.insert("ID_PSMScore", {writeScore(out, score)});
  for (const auto& run : data.getRuns())
    writeRun(out, run);
  for (const auto& inference : data.getInferenceResults())
  {
    const auto processing = writeProcessing(out, inference.proteins);
    const auto parent_score = inference.parent_score ? Value(writeScore(out, *inference.parent_score)) : Value(nullptr);
    const auto group_score = inference.group_score ? Value(writeScore(out, *inference.group_score)) : Value(nullptr);
    const auto id = out.insert("ID_Inference", {inference.identifier, processing, parent_score, group_score});
    for (const auto& input : inference.inputs)
      out.insert("ID_Input",
                 {id, input.run_identifier, input.run_uuid, input.score ? Value(writeScore(out, *input.score)) : Value(nullptr), input.selection});
    for (const auto& [alias, identity] : inference.qualified_accessions)
      out.insert("ID_Accession", {id, alias, identity.database, identity.accession});
  }
  // Build child indexes after bulk insertion. Insertion order supplies vector order;
  // persistent query/match IDs are separate strings covering the full UInt64 range.
  Row tables(db, "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'ID_%' ORDER BY name");
  std::vector<std::string> children;
  while (tables.executeStep())
  {
    const auto table = str(tables, "name");
    if (table == "ID_Metadata") continue;
    Row columns(db, "PRAGMA table_info(" + table + ")");
    while (columns.executeStep())
      if (str(columns, "name") == "parent_id") children.push_back(table);
  }
  for (const auto& table : children)
    db.exec("CREATE INDEX " + table + "_parent ON " + table + "(parent_id,id)");
  // A view resolves the small name dictionary for ad hoc SQL inspection.
  db.exec("CREATE VIEW ID_MetadataValue AS SELECT m.id,o.name AS owner,m.parent_id,n.name AS "
          "name,m.type,m.unit_type,m.unit,m.text_value,m.integer_value,m.real_value,m.real_value_bits FROM ID_Metadata m JOIN ID_MetadataName o ON "
          "o.id=m.owner_id JOIN ID_MetadataName n ON n.id=m.name_id");
  // The metadata uniqueness index also supplies owner/record/key scan order.
}

void loadOMSIdentifications(SQLite::Database& db, IdentificationData& data)
{
  // SQLite does not enforce constraints retroactively on an externally edited file.
  Row integrity(db, "PRAGMA foreign_key_check");
  if (integrity.executeStep()) invalid("Broken foreign key in " + integrity.getColumn(0).getString());
  Reader in(db);
  std::map<Key, ID::ScoreDefinition> definitions;
  Row score_rows(db, "SELECT * FROM ID_Score ORDER BY id");
  while (score_rows.executeStep())
    definitions.emplace(number(score_rows, "id"), readScore(in, score_rows));
  std::vector<ID::ScoreDefinition> scores;
  Row schema_rows(db, "SELECT * FROM ID_PSMScore ORDER BY id");
  while (schema_rows.executeStep())
  {
    if (number(schema_rows, "id") != static_cast<Key>(scores.size()) + 1) invalid("Noncontiguous score schema");
    scores.push_back(definitions.at(number(schema_rows, "definition_id")));
  }
  std::map<Key, ProteinIdentification> processing;
  Row processing_rows(db, "SELECT * FROM ID_Processing ORDER BY id");
  while (processing_rows.executeStep())
    processing.emplace(number(processing_rows, "id"), readProcessing(in, processing_rows));
  Row runs(db, "SELECT * FROM ID_Run ORDER BY id");
  while (runs.executeStep())
    data.addRun(readRun(in, runs, scores, processing));
  Row results(db, "SELECT * FROM ID_Inference ORDER BY id");
  while (results.executeStep())
  {
    ID::InferenceResult result;
    const auto id = number(results, "id");
    result.identifier = str(results, "identifier");
    result.proteins = std::move(processing.at(number(results, "processing_id")));
    if (! null(results, "parent_score_id")) result.parent_score = definitions.at(number(results, "parent_score_id"));
    if (! null(results, "group_score_id")) result.group_score = definitions.at(number(results, "group_score_id"));
    in.children("ID_Input", id, [&](const Row& row) {
      ID::InferenceInput input;
      input.run_identifier = str(row, "run_identifier");
      input.run_uuid = str(row, "run_uuid");
      input.selection = str(row, "selection");
      if (! null(row, "score_id")) input.score = definitions.at(number(row, "score_id"));
      result.inputs.push_back(std::move(input));
    });
    in.children("ID_Accession", id, [&](const Row& row) {
      if (! result.qualified_accessions.emplace(str(row, "alias"), ID::QualifiedAccession {str(row, "database"), str(row, "accession")}).second)
        invalid("Duplicate qualified accession");
    });
    data.addInferenceResult(std::move(result));
  }
  in.finish();
  data.validate();
}
} // namespace OpenMS::Internal
