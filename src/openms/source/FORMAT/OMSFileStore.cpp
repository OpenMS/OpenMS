// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFileStore.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <SQLiteCpp/Database.h>
#include <SQLiteCpp/Transaction.h>

#include <sqlite3.h>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
  constexpr int version_number = 5; // increase this whenever the DB schema changes!

  void raiseDBError_(const String& error, int line, const char* function,
                     const String& context, const String& query)
  {
    String msg = context + ": " + error;
    if (!query.empty())
    {
      msg += String("\nQuery was: ") + query;
    }
    throw Exception::FailedAPICall(__FILE__, line, function, msg);
  }

  bool execAndReset(SQLite::Statement& query, int expected_modifications)
  {
    auto ret = query.exec();
    query.reset();
    return ret == expected_modifications;
  }

  void execWithExceptionAndReset(SQLite::Statement& query, int expected_modifications, int line, const char* function, const char* context)
  {
    if (!execAndReset(query, expected_modifications))
    {
      raiseDBError_(query.getErrorMsg(), line, function, context);
    }
  }

  OMSFileStore::OMSFileStore(const String& filename, LogType log_type)
  {
    setLogType(log_type);
    File::remove(filename); // nuke the file (SQLite cannot overwrite it)
    db_ = make_unique<SQLite::Database>(filename, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE); // throws on error
    // foreign key constraints are disabled by default - turn them on:
    // @TODO: performance impact? (seems negligible, but should be tested more)
    db_->exec("PRAGMA foreign_keys = ON");
    // disable synchronous filesystem access and the rollback journal to greatly
    // increase write performance - since we write a new output file every time,
    // we don't have to worry about database consistency:
    db_->exec("PRAGMA synchronous = OFF");
    db_->exec("PRAGMA journal_mode = OFF");
    db_->exec("PRAGMA foreign_keys = ON");
    db_->exec("PRAGMA foreign_keys = ON");
  }

  OMSFileStore::~OMSFileStore() = default;

  void OMSFileStore::createTable_(const String& name, const String& definition, bool may_exist)
  {
    String sql_create = "CREATE TABLE ";
    if (may_exist) sql_create += "IF NOT EXISTS ";
    sql_create += name + " (" + definition + ")";
    db_->exec(sql_create);
  }


  void OMSFileStore::storeVersionAndDate_()
  {
    createTable_("version",
                 "OMSFile INT NOT NULL, "       \
                 "date TEXT NOT NULL, "         \
                 "OpenMS TEXT, "                \
                 "build_date TEXT");

    SQLite::Statement query(*db_, "INSERT INTO version VALUES ("  \
                                 ":format_version, "             \
                                 "datetime('now'), "             \
                                 ":openms_version, "             \
                                 ":build_date)");
    query.bind(":format_version", version_number);
    query.bind(":openms_version", VersionInfo::getVersion());
    query.bind(":build_date", VersionInfo::getTime());
    query.exec();
  }


  void OMSFileStore::createTableMoleculeType_()
  {
    createTable_("ID_MoleculeType",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "molecule_type TEXT UNIQUE NOT NULL");
    auto sql_insert =
      "INSERT INTO ID_MoleculeType VALUES "     \
      "(1, 'PROTEIN'), "                        \
      "(2, 'COMPOUND'), "                       \
      "(3, 'RNA')";
    db_->exec(sql_insert);
  }

  void OMSFileStore::createTableDataValue_DataType_()
  {
    createTable_("DataValue_DataType",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "data_type TEXT UNIQUE NOT NULL");
    auto sql_insert =
      "INSERT INTO DataValue_DataType VALUES " \
      "(1, 'STRING_VALUE'), "                  \
      "(2, 'INT_VALUE'), "                     \
      "(3, 'DOUBLE_VALUE'), "                  \
      "(4, 'STRING_LIST'), "                   \
      "(5, 'INT_LIST'), "                      \
      "(6, 'DOUBLE_LIST')";
    db_->exec(sql_insert);
  }

  void OMSFileStore::createTableCVTerm_()
  {
    createTable_("CVTerm",
                 "id INTEGER PRIMARY KEY NOT NULL, "        \
                 "accession TEXT UNIQUE, "                  \
                 "name TEXT NOT NULL, "                     \
                 "cv_identifier_ref TEXT, "                 \
                 // does this constrain "name" if "accession" is NULL?
                 "UNIQUE (accession, name)");
    // @TODO: add support for unit and value
    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT OR IGNORE INTO CVTerm VALUES ("   \
                                 "NULL, "                                  \
                                 ":accession, "                            \
                                 ":name, "                                 \
                                 ":cv_identifier_ref)");
    prepared_queries_.emplace("CVTerm", std::move(query));
    // alternative query if CVTerm already exists:
    auto query2 = make_unique<SQLite::Statement>(*db_, "SELECT id FROM CVTerm "                 \
                                  "WHERE accession = :accession AND name = :name");
    prepared_queries_.emplace("CVTerm_2", std::move(query2));
  }


  OMSFileStore::Key OMSFileStore::storeCVTerm_(const CVTerm& cv_term)
  {
    // this assumes the "CVTerm" table exists already!
    auto& query = *prepared_queries_["CVTerm"];
    if (cv_term.getAccession().empty()) // use NULL for empty accessions
    {
      query.bind(":accession");
    }
    else
    {
      query.bind(":accession", cv_term.getAccession());
    }
    query.bind(":name", cv_term.getName());
    query.bind(":cv_identifier_ref", cv_term.getCVIdentifierRef());
    if (execAndReset(query, 1)) // one row was inserted
    {
      return db_->getLastInsertRowid();
    }

    // else: insert has failed, record must already exist - get the key:
    auto& alt_query = *prepared_queries_["CVTerm_2"];
    alt_query.reset();                  // get ready for a new execution
    if (cv_term.getAccession().empty()) // use NULL for empty accessions
    {
      alt_query.bind(":accession");
    }
    else
    {
      alt_query.bind(":accession", cv_term.getAccession());
    }
    alt_query.bind(":name", cv_term.getName());
    if (!alt_query.executeStep())
    {
      raiseDBError_(alt_query.getErrorMsg(), __LINE__, OPENMS_PRETTY_FUNCTION, "error querying database");
    }

    return Key(alt_query.getColumn(0).getInt64());
  }


  void OMSFileStore::createTableMetaInfo_(const String& parent_table, const String& key_column)
  {
    if (!db_->tableExists("DataValue_DataType")) createTableDataValue_DataType_();

    String parent_ref = parent_table + " (" + key_column + ")";
    String table = parent_table + "_MetaInfo";
    // for the data_type_id, empty values are represented using NULL
    createTable_(
      table,
      "parent_id INTEGER NOT NULL, "                            \
      "name TEXT NOT NULL, "                                    \
      "data_type_id INTEGER, "                                  \
      "value TEXT, "                                            \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_ref + ", " \
      "FOREIGN KEY (data_type_id) REFERENCES DataValue_DataType (id), " \
      "PRIMARY KEY (parent_id, name)");

    // @TODO: add support for units
    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO " + table +
                                 " VALUES ("                    \
                                 ":parent_id, "                 \
                                 ":name, "                      \
                                 ":data_type_id, "              \
                                 ":value)");
    prepared_queries_.emplace(table, std::move(query));
  }


  void OMSFileStore::storeMetaInfo_(const MetaInfoInterface& info, const String& parent_table, Key parent_id)
  {
    if (info.isMetaEmpty()) return;

    // this assumes the "..._MetaInfo" and "DataValue_DataType" tables exist already!
    auto& query = *prepared_queries_[parent_table + "_MetaInfo"];
    query.bind(":parent_id", parent_id);
    // this is inefficient, but MetaInfoInterface doesn't support iteration:
    vector<String> info_keys;
    info.getKeys(info_keys);
    for (const String& info_key : info_keys)
    {
      query.bind(":name", info_key);

      const DataValue& value = info.getMetaValue(info_key);
      if (value.isEmpty()) // use NULL as the type for empty values
      {
        query.bind(":data_type_id");
      }
      else
      {
        query.bind(":data_type_id", int(value.valueType()) + 1);
      }
      query.bind(":value", value.toString());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
    }
  }


  void OMSFileStore::createTableAppliedProcessingStep_(const String& parent_table)
  {
    String table = parent_table + "_AppliedProcessingStep";
    createTable_(
      table,
      "parent_id INTEGER NOT NULL, "                                    \
      "processing_step_id INTEGER, "                                    \
      "processing_step_order INTEGER NOT NULL, "                        \
      "score_type_id INTEGER, "                                         \
      "score REAL, "                                                    \
      "UNIQUE (parent_id, processing_step_id, score_type_id), "         \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_table + " (id), "  \
      "FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id), "      \
      "FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)");
    // @TODO: add constraint that "processing_step_id" and "score_type_id" can't both be NULL
    // @TODO: add constraint that "processing_step_order" must match "..._id"?
    // @TODO: normalize table? (splitting into multiple tables is awkward here)
    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO " + table +
                                 " VALUES ("  \
                                 ":parent_id, "
                                 ":processing_step_id, "
                                 ":processing_step_order, "
                                 ":score_type_id, "
                                 ":score)");
    prepared_queries_.emplace(table, std::move(query));
  }


  void OMSFileStore::storeAppliedProcessingStep_(const ID::AppliedProcessingStep& step, Size step_order,
                                                 const String& parent_table, Key parent_id)
  {
    // this assumes the "..._AppliedProcessingStep" table exists already!
    auto& query = *prepared_queries_[parent_table + "_AppliedProcessingStep"];
    query.bind(":parent_id", parent_id);
    query.bind(":processing_step_order", int(step_order));
    if (step.processing_step_opt)
    {
      query.bind(":processing_step_id", processing_step_keys_[&(**step.processing_step_opt)]);
      if (step.scores.empty()) // insert processing step information only
      {
        query.bind(":score_type_id"); // NULL
        query.bind(":score"); // NULL
        execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
    }
    else // use NULL for missing processing step reference
    {
      query.bind(":processing_step_id");
    }
    for (const auto& score_pair : step.scores)
    {
      query.bind(":score_type_id", score_type_keys_[&(*score_pair.first)]);
      query.bind(":score", score_pair.second);
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
    }
  }


  void OMSFileStore::storeScoreTypes_(const IdentificationData& id_data)
  {
    if (id_data.getScoreTypes().empty()) return;

    createTableCVTerm_();
    createTable_(
      "ID_ScoreType",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "cv_term_id INTEGER NOT NULL, "                                   \
      "higher_better NUMERIC NOT NULL CHECK (higher_better in (0, 1)), " \
      "FOREIGN KEY (cv_term_id) REFERENCES CVTerm (id)");

    SQLite::Statement query(*db_, "INSERT INTO ID_ScoreType VALUES ("       \
                  ":id, "                                   \
                  ":cv_term_id, "                           \
                  ":higher_better)");
    Key id = 1;
    for (const ID::ScoreType& score_type : id_data.getScoreTypes())
    {
      Key cv_id = storeCVTerm_(score_type.cv_term);
      query.bind(":id", id);
      query.bind(":cv_term_id", cv_id);
      query.bind(":higher_better", int(score_type.higher_better));
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      score_type_keys_[&score_type] = id;
      ++id;
    }
  }


  void OMSFileStore::storeInputFiles_(const IdentificationData& id_data)
  {
    if (id_data.getInputFiles().empty()) return;

    createTable_("ID_InputFile",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "name TEXT UNIQUE NOT NULL, "        \
                 "experimental_design_id TEXT, "      \
                 "primary_files TEXT");

    SQLite::Statement query(*db_, "INSERT INTO ID_InputFile VALUES ("  \
                  ":id, "                              \
                  ":name, "                            \
                  ":experimental_design_id, "          \
                  ":primary_files)");
    Key id = 1;
    for (const ID::InputFile& input : id_data.getInputFiles())
    {
      query.bind(":id", id);
      query.bind(":name", input.name);
      query.bind(":experimental_design_id",
                      input.experimental_design_id);
      // @TODO: what if a primary file name contains ","?
      String primary_files = ListUtils::concatenate(input.primary_files, ",");
      query.bind(":primary_files", primary_files);
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      input_file_keys_[&input] = id;
      ++id;
    }
  }


  void OMSFileStore::storeProcessingSoftwares_(const IdentificationData& id_data)
  {
    if (id_data.getProcessingSoftwares().empty()) return;

    createTable_("ID_ProcessingSoftware",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "name TEXT NOT NULL, "               \
                 "version TEXT, "                     \
                 "UNIQUE (name, version)");

    SQLite::Statement query(*db_, "INSERT INTO ID_ProcessingSoftware VALUES ("  \
                  ":id, "                                           \
                  ":name, "                                         \
                  ":version)");
    bool any_scores = false; // does any software have assigned scores stored?
    Key id = 1;
    for (const ID::ProcessingSoftware& software : id_data.getProcessingSoftwares())
    {
      if (!software.assigned_scores.empty()) any_scores = true;
      query.bind(":id", id);
      query.bind(":name", software.getName());
      query.bind(":version", software.getVersion());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      processing_software_keys_[&software] = id;
      ++id;
    }
    if (any_scores)
    {
      createTable_(
        "ID_ProcessingSoftware_AssignedScore",
        "software_id INTEGER NOT NULL, "                                \
        "score_type_id INTEGER NOT NULL, "                              \
        "score_type_order INTEGER NOT NULL, "                           \
        "UNIQUE (software_id, score_type_id), "                         \
        "UNIQUE (software_id, score_type_order), "                      \
        "FOREIGN KEY (software_id) REFERENCES ID_ProcessingSoftware (id), " \
        "FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id)");

      SQLite::Statement query2(*db_,
        "INSERT INTO ID_ProcessingSoftware_AssignedScore VALUES ("      \
        ":software_id, "                                                \
        ":score_type_id, "                                              \
        ":score_type_order)");
      for (const ID::ProcessingSoftware& software : id_data.getProcessingSoftwares())
      {
        query2.bind(":software_id", processing_software_keys_[&software]);
        Size counter = 0;
        for (ID::ScoreTypeRef score_type_ref : software.assigned_scores)
        {
          query2.bind(":score_type_id", score_type_keys_[&(*score_type_ref)]);
          query2.bind(":score_type_order", int(++counter));
          execWithExceptionAndReset(query2, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }
      }
    }
  }


  void OMSFileStore::storeDBSearchParams_(const IdentificationData& id_data)
  {
    if (id_data.getDBSearchParams().empty()) return;

    if (!db_->tableExists("ID_MoleculeType")) createTableMoleculeType_();

    createTable_(
      "ID_DBSearchParam",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "molecule_type_id INTEGER NOT NULL, "                             \
      "mass_type_average NUMERIC NOT NULL CHECK (mass_type_average in (0, 1)) DEFAULT 0, " \
      "database TEXT, "                                                 \
      "database_version TEXT, "                                         \
      "taxonomy TEXT, "                                                 \
      "charges TEXT, "                                                  \
      "fixed_mods TEXT, "                                               \
      "variable_mods TEXT, "                                            \
      "precursor_mass_tolerance REAL, "                                 \
      "fragment_mass_tolerance REAL, "                                  \
      "precursor_tolerance_ppm NUMERIC NOT NULL CHECK (precursor_tolerance_ppm in (0, 1)) DEFAULT 0, " \
      "fragment_tolerance_ppm NUMERIC NOT NULL CHECK (fragment_tolerance_ppm in (0, 1)) DEFAULT 0, " \
      "digestion_enzyme TEXT, "                                         \
      "enzyme_term_specificity TEXT, " // new in version 2!
      "missed_cleavages NUMERIC, "                                      \
      "min_length NUMERIC, "                                            \
      "max_length NUMERIC, "                                            \
      "FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)");

    SQLite::Statement query(*db_, "INSERT INTO ID_DBSearchParam VALUES (" \
                  ":id, "                                 \
                  ":molecule_type_id, "                   \
                  ":mass_type_average, "                  \
                  ":database, "                           \
                  ":database_version, "                   \
                  ":taxonomy, "                           \
                  ":charges, "                            \
                  ":fixed_mods, "                         \
                  ":variable_mods, "                      \
                  ":precursor_mass_tolerance, "           \
                  ":fragment_mass_tolerance, "            \
                  ":precursor_tolerance_ppm, "            \
                  ":fragment_tolerance_ppm, "             \
                  ":digestion_enzyme, "                   \
                  ":enzyme_term_specificity, "            \
                  ":missed_cleavages, "                   \
                  ":min_length, "                         \
                  ":max_length)");
    Key id = 1;
    for (const ID::DBSearchParam& param : id_data.getDBSearchParams())
    {
      query.bind(":id", id);
      query.bind(":molecule_type_id", int(param.molecule_type) + 1);
      query.bind(":mass_type_average", int(param.mass_type));
      query.bind(":database", param.database);
      query.bind(":database_version", param.database_version);
      query.bind(":taxonomy", param.taxonomy);
      String charges = ListUtils::concatenate(param.charges, ",");
      query.bind(":charges", charges);
      String fixed_mods = ListUtils::concatenate(param.fixed_mods, ",");
      query.bind(":fixed_mods", fixed_mods);
      String variable_mods = ListUtils::concatenate(param.variable_mods, ",");
      query.bind(":variable_mods", variable_mods);
      query.bind(":precursor_mass_tolerance", param.precursor_mass_tolerance);
      query.bind(":fragment_mass_tolerance", param.fragment_mass_tolerance);
      query.bind(":precursor_tolerance_ppm", int(param.precursor_tolerance_ppm));
      query.bind(":fragment_tolerance_ppm", int(param.fragment_tolerance_ppm));
      if (param.digestion_enzyme != nullptr)
      {
        query.bind(":digestion_enzyme", param.digestion_enzyme->getName());
      }
      else // bind NULL value
      {
        query.bind(":digestion_enzyme");
      }
      query.bind(":enzyme_term_specificity",
                      EnzymaticDigestion::NamesOfSpecificity[param.enzyme_term_specificity]);
      query.bind(":missed_cleavages", uint32_t(param.missed_cleavages));
      query.bind(":min_length", uint32_t(param.min_length));
      query.bind(":max_length", uint32_t(param.max_length));
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      search_param_keys_[&param] = id;
      ++id;
    }
  }


  void OMSFileStore::storeProcessingSteps_(const IdentificationData& id_data)
  {
    if (id_data.getProcessingSteps().empty()) return;

    createTable_(
      "ID_ProcessingStep",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "software_id INTEGER NOT NULL, "                                  \
      "date_time TEXT, "                                                \
      "search_param_id INTEGER, "                                       \
      "FOREIGN KEY (search_param_id) REFERENCES ID_DBSearchParam (id)");
    // @TODO: add support for processing actions
    // @TODO: store primary files in a separate table (like input files)?
    // @TODO: store (optional) search param reference in a separate table?

    SQLite::Statement query(*db_, "INSERT INTO ID_ProcessingStep VALUES ("      \
                  ":id, "                                       \
                  ":software_id, "                              \
                  ":date_time, "                                \
                  ":search_param_id)");
    bool any_input_files = false;
    Key id = 1;
    // use iterator here because we need one to look up the DB search params:
    for (ID::ProcessingStepRef step_ref = id_data.getProcessingSteps().begin();
         step_ref != id_data.getProcessingSteps().end(); ++step_ref)
    {
      const ID::ProcessingStep& step = *step_ref;
      if (!step.input_file_refs.empty()) any_input_files = true;
      query.bind(":id", id);
      query.bind(":software_id", processing_software_keys_[&(*step.software_ref)]);
      query.bind(":date_time", step.date_time.get());
      auto pos = id_data.getDBSearchSteps().find(step_ref);
      if (pos != id_data.getDBSearchSteps().end())
      {
        query.bind(":search_param_id", search_param_keys_[&(*pos->second)]);
      }
      else
      {
        query.bind(":search_param_id"); // NULL
      }
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      processing_step_keys_[&step] = id;
      ++id;
    }
    if (any_input_files)
    {
      createTable_(
        "ID_ProcessingStep_InputFile",
        "processing_step_id INTEGER NOT NULL, "                         \
        "input_file_id INTEGER NOT NULL, "                              \
        "FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id), " \
        "FOREIGN KEY (input_file_id) REFERENCES ID_InputFile (id), "      \
        "UNIQUE (processing_step_id, input_file_id)");

      SQLite::Statement query2(*db_, "INSERT INTO ID_ProcessingStep_InputFile VALUES (" \
                    ":processing_step_id, "                             \
                    ":input_file_id)");

      for (const ID::ProcessingStep& step : id_data.getProcessingSteps())
      {
        query2.bind(":processing_step_id", processing_step_keys_[&step]);
        for (ID::InputFileRef input_file_ref : step.input_file_refs)
        {
          query2.bind(":input_file_id", input_file_keys_[&(*input_file_ref)]);
          execWithExceptionAndReset(query2, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }
      }
    }
    storeMetaInfos_(id_data.getProcessingSteps(), "ID_ProcessingStep", processing_step_keys_);
  }


  void OMSFileStore::storeObservations_(const IdentificationData& id_data)
  {
    if (id_data.getObservations().empty()) return;

    createTable_("ID_Observation",
                 "id INTEGER PRIMARY KEY NOT NULL, "                    \
                 "data_id TEXT NOT NULL, "                              \
                 "input_file_id INTEGER NOT NULL, "                     \
                 "rt REAL, "                                            \
                 "mz REAL, "                                            \
                 "UNIQUE (data_id, input_file_id), "                    \
                 "FOREIGN KEY (input_file_id) REFERENCES ID_InputFile (id)");

    SQLite::Statement query(*db_, "INSERT INTO ID_Observation VALUES (" \
                  ":id, "                             \
                  ":data_id, "                        \
                  ":input_file_id, "                  \
                  ":rt, "                             \
                  ":mz)");
    Key id = 1;
    for (const ID::Observation& obs : id_data.getObservations())
    {
      query.bind(":id", id);
      query.bind(":data_id", obs.data_id);
      query.bind(":input_file_id", input_file_keys_[&(*obs.input_file)]);
      if (obs.rt == obs.rt)
      {
        query.bind(":rt", obs.rt);
      }
      else // NaN
      {
        query.bind(":rt"); // NULL
      }
      if (obs.mz == obs.mz)
      {
        query.bind(":mz", obs.mz);
      }
      else // NaN
      {
        query.bind(":mz"); // NULL
      }
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      observation_keys_[&obs] = id;
      ++id;
    }
    storeMetaInfos_(id_data.getObservations(), "ID_Observation", observation_keys_);
  }


  void OMSFileStore::storeParentSequences_(const IdentificationData& id_data)
  {
    if (id_data.getParentSequences().empty()) return;

    if (!db_->tableExists("ID_MoleculeType")) createTableMoleculeType_();

    createTable_(
      "ID_ParentSequence",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "accession TEXT UNIQUE NOT NULL, "                                \
      "molecule_type_id INTEGER NOT NULL, "                             \
      "sequence TEXT, "                                                 \
      "description TEXT, "                                              \
      "coverage REAL, "                                                 \
      "is_decoy NUMERIC NOT NULL CHECK (is_decoy in (0, 1)) DEFAULT 0, " \
      "FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)");

    SQLite::Statement query(*db_, "INSERT INTO ID_ParentSequence VALUES ("  \
                  ":id, "                                   \
                  ":accession, "                            \
                  ":molecule_type_id, "                     \
                  ":sequence, "                             \
                  ":description, "                          \
                  ":coverage, "                             \
                  ":is_decoy)");
    Key id = 1;
    for (const ID::ParentSequence& parent : id_data.getParentSequences())
    {
      query.bind(":id", id);
      query.bind(":accession", parent.accession);
      query.bind(":molecule_type_id", int(parent.molecule_type) + 1);
      query.bind(":sequence", parent.sequence);
      query.bind(":description", parent.description);
      query.bind(":coverage", parent.coverage);
      query.bind(":is_decoy", int(parent.is_decoy));
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      parent_sequence_keys_[&parent] = id;
      ++id;
    }
    storeScoredProcessingResults_(id_data.getParentSequences(), "ID_ParentSequence", parent_sequence_keys_);
  }


  void OMSFileStore::storeParentGroupSets_(const IdentificationData& id_data)
  {
    if (id_data.getParentGroupSets().empty()) return;

    createTable_("ID_ParentGroupSet",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "label TEXT UNIQUE");

    createTable_(
      "ID_ParentGroup",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "grouping_id INTEGER NOT NULL, "                                  \
      "score_type_id INTEGER, "                                         \
      "score REAL, "                                                    \
      "UNIQUE (id, score_type_id), "                                    \
      "FOREIGN KEY (grouping_id) REFERENCES ID_ParentGroupSet (id)");

    createTable_(
      "ID_ParentGroup_ParentSequence",
      "group_id INTEGER NOT NULL, "                                     \
      "parent_id INTEGER NOT NULL, "                                    \
      "UNIQUE (group_id, parent_id), "                                  \
      "FOREIGN KEY (group_id) REFERENCES ID_ParentGroup (id), "         \
      "FOREIGN KEY (parent_id) REFERENCES ID_ParentSequence (id)");

    SQLite::Statement query_grouping(*db_, "INSERT INTO ID_ParentGroupSet VALUES ("     \
                           ":id, "                                      \
                           ":label)");

    SQLite::Statement query_group(*db_, "INSERT INTO ID_ParentGroup VALUES ("          \
                        ":id, "                                        \
                        ":grouping_id, "                               \
                        ":score_type_id, "                             \
                        ":score)");

   SQLite::Statement query_parent(*db_, "INSERT INTO ID_ParentGroup_ParentSequence VALUES (" \
                         ":group_id, "                                  \
                         ":parent_id)");

    Key grouping_id = 1;
    Key group_id = 1;
    for (const ID::ParentGroupSet& grouping : id_data.getParentGroupSets())
    {
      query_grouping.bind(":id", grouping_id);
      query_grouping.bind(":label", grouping.label);
      execWithExceptionAndReset(query_grouping, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      for (const ID::ParentGroup& group : grouping.groups)
      {
        query_group.bind(":id", group_id);
        query_group.bind(":grouping_id", grouping_id);
        if (group.scores.empty()) // store group with an empty score
        {
          query_group.bind(":score_type_id");
          query_group.bind(":score");
          execWithExceptionAndReset(query_group, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }
        else // store group multiple times with different scores
        {
          for (const auto& score_pair : group.scores)
          {
            query_group.bind(":score_type_id", score_type_keys_[&(*score_pair.first)]);
            query_group.bind(":score", score_pair.second);
            execWithExceptionAndReset(query_group, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
          }
        }

        query_parent.bind(":group_id", group_id);
        for (ID::ParentSequenceRef parent_ref : group.parent_refs)
        {
          query_parent.bind(":parent_id", parent_sequence_keys_[&(*parent_ref)]);
          execWithExceptionAndReset(query_parent, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }
        ++group_id;
      }
      parent_grouping_keys_[&grouping] = grouping_id;
      ++grouping_id;
    }

    storeScoredProcessingResults_(id_data.getParentGroupSets(), "ID_ParentGroupSet", parent_grouping_keys_);
  }


  void OMSFileStore::createTableIdentifiedMolecule_()
  {
    if (!db_->tableExists("ID_MoleculeType")) createTableMoleculeType_();

    // use one table for all types of identified molecules to allow foreign key
    // references from the input match table:
    createTable_(
      "ID_IdentifiedMolecule",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "molecule_type_id INTEGER NOT NULL, "                             \
      "identifier TEXT NOT NULL, "                                      \
      "UNIQUE (molecule_type_id, identifier), "                         \
      "FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)");
    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO ID_IdentifiedMolecule VALUES ("          \
                  ":id, "                                               \
                  ":molecule_type_id, "                                 \
                  ":identifier)");
    prepared_queries_.emplace("ID_IdentifiedMolecule", std::move(query));
  }


  void OMSFileStore::storeIdentifiedCompounds_(const IdentificationData& id_data)
  {
    if (id_data.getIdentifiedCompounds().empty()) return;

    if (!db_->tableExists("ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }
    auto& query_molecule = *prepared_queries_["ID_IdentifiedMolecule"];
    query_molecule.bind(":molecule_type_id", int(ID::MoleculeType::COMPOUND) + 1);

    createTable_(
      "ID_IdentifiedCompound",
      "molecule_id INTEGER UNIQUE NOT NULL , "                          \
      "formula TEXT, "                                                  \
      "name TEXT, "                                                     \
      "smile TEXT, "                                                    \
      "inchi TEXT, "                                                    \
      "FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)");
    SQLite::Statement query_compound(*db_, "INSERT INTO ID_IdentifiedCompound VALUES (" \
                           ":molecule_id, "                             \
                           ":formula, "                                 \
                           ":name, "                                    \
                           ":smile, "                                   \
                           ":inchi)");
    Key id = 1;
    for (const ID::IdentifiedCompound& compound : id_data.getIdentifiedCompounds())
    {
      query_molecule.bind(":id", id);
      query_molecule.bind(":identifier", compound.identifier);
      execWithExceptionAndReset(query_molecule, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      query_compound.bind(":molecule_id", id);
      query_compound.bind(":formula", compound.formula.toString());
      query_compound.bind(":name", compound.name);
      query_compound.bind(":smile", compound.name);
      query_compound.bind(":inchi", compound.inchi);
      execWithExceptionAndReset(query_compound, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");


      identified_compound_keys_[&compound] = id;
      ++id;
    }
    storeScoredProcessingResults_(id_data.getIdentifiedCompounds(), "ID_IdentifiedMolecule",
                                  identified_compound_keys_);
  }


  void OMSFileStore::storeIdentifiedSequences_(const IdentificationData& id_data)
  {
    if (id_data.getIdentifiedPeptides().empty() &&
        id_data.getIdentifiedOligos().empty()) return;

    if (!db_->tableExists("ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }
    auto& query = *prepared_queries_["ID_IdentifiedMolecule"];

    bool any_parent_matches = false;
    // identified compounds get stored earlier and use the same key range as peptides/oligos:
    Key id = id_data.getIdentifiedCompounds().size() + 1;
    // store peptides:
    query.bind(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    for (const ID::IdentifiedPeptide& peptide : id_data.getIdentifiedPeptides())
    {
      if (!peptide.parent_matches.empty()) any_parent_matches = true;
      query.bind(":id", id);
      query.bind(":identifier", peptide.sequence.toString());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      identified_peptide_keys_[&peptide] = id;
      ++id;
    }
    storeScoredProcessingResults_(id_data.getIdentifiedPeptides(), "ID_IdentifiedMolecule",
                                  identified_peptide_keys_);
    // store RNA oligos:
    query.bind(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    for (const ID::IdentifiedOligo& oligo : id_data.getIdentifiedOligos())
    {
      if (!oligo.parent_matches.empty()) any_parent_matches = true;
      query.bind(":id", id); // use address as primary key
      query.bind(":identifier", oligo.sequence.toString());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      identified_oligo_keys_[&oligo] = id;
      ++id;
    }
    storeScoredProcessingResults_(id_data.getIdentifiedOligos(), "ID_IdentifiedMolecule",
                                  identified_oligo_keys_);

    if (any_parent_matches)
    {
      createTableParentMatches_();
      for (const ID::IdentifiedPeptide& peptide : id_data.getIdentifiedPeptides())
      {
        if (!peptide.parent_matches.empty())
        {
          storeParentMatches_(peptide.parent_matches, identified_peptide_keys_[&peptide]);
        }
      }
      for (const ID::IdentifiedOligo& oligo : id_data.getIdentifiedOligos())
      {
        if (!oligo.parent_matches.empty())
        {
          storeParentMatches_(oligo.parent_matches, identified_oligo_keys_[&oligo]);
        }
      }
    }
  }


  void OMSFileStore::createTableParentMatches_()
  {
    createTable_(
      "ID_ParentMatch",
      "molecule_id INTEGER NOT NULL, "                                  \
      "parent_id INTEGER NOT NULL, "                                    \
      "start_pos NUMERIC, "                                             \
      "end_pos NUMERIC, "                                               \
      "left_neighbor TEXT, "                                            \
      "right_neighbor TEXT, "                                           \
      "UNIQUE (molecule_id, parent_id, start_pos, end_pos), "           \
      "FOREIGN KEY (parent_id) REFERENCES ID_ParentSequence (id), "     \
      "FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)");
    // prepare query for inserting data:
    auto query = make_unique<SQLite::Statement>(*db_, "INSERT INTO ID_ParentMatch VALUES (" \
                                                ":molecule_id, "        \
                                                ":parent_id, "          \
                                                ":start_pos, "          \
                                                ":end_pos, "            \
                                                ":left_neighbor, "      \
                                                ":right_neighbor)");
    prepared_queries_.emplace("ID_ParentMatch", std::move(query));
  }


  void OMSFileStore::storeParentMatches_(const ID::ParentMatches& matches, Key molecule_id)
  {
    // this assumes the "ID_ParentMatch" table exists already!
    auto& query = *prepared_queries_["ID_ParentMatch"];
    // @TODO: cache the prepared query between function calls somehow?
    query.bind(":molecule_id", molecule_id);
    for (const auto& pair : matches)
    {
      query.bind(":parent_id", parent_sequence_keys_[&(*pair.first)]);
      for (const auto& match : pair.second)
      {
        if (match.start_pos != ID::ParentMatch::UNKNOWN_POSITION)
        {
          query.bind(":start_pos", uint32_t(match.start_pos));
        }
        else // use NULL value
        {
          query.bind(":start_pos");
        }
        if (match.end_pos != ID::ParentMatch::UNKNOWN_POSITION)
        {
          query.bind(":end_pos", uint32_t(match.end_pos));
        }
        else // use NULL value
        {
          query.bind(":end_pos");
        }
        query.bind(":left_neighbor", match.left_neighbor);
        query.bind(":right_neighbor", match.right_neighbor);
        execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
    }
  }


  void OMSFileStore::storeAdducts_(const IdentificationData& id_data)
  {
    if (id_data.getAdducts().empty()) return;

    createTable_("AdductInfo",
                 "id INTEGER PRIMARY KEY NOT NULL, "                    \
                 "name TEXT, "                                          \
                 "formula TEXT NOT NULL, "                              \
                 "charge INTEGER NOT NULL, "                            \
                 "mol_multiplier INTEGER NOT NULL CHECK (mol_multiplier > 0) DEFAULT 1, " \
                 "UNIQUE (formula, charge)");

    SQLite::Statement query(*db_, "INSERT INTO AdductInfo VALUES (" \
                            ":id, "                                 \
                            ":name, "                               \
                            ":formula, "                            \
                            ":charge, "                             \
                            ":mol_multiplier)");
    Key id = 1;
    for (const AdductInfo& adduct : id_data.getAdducts())
    {
      query.bind(":id", id);
      query.bind(":name", adduct.getName());
      query.bind(":formula", adduct.getEmpiricalFormula().toString());
      query.bind(":charge", adduct.getCharge());
      query.bind(":mol_multiplier", adduct.getMolMultiplier());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      adduct_keys_[&adduct] = id;
      ++id;
    }
  }


  OMSFileStore::Key OMSFileStore::getDatabaseKey_(const ID::IdentifiedMolecule& molecule_var)
  {
    switch (molecule_var.getMoleculeType())
    {
      case ID::MoleculeType::PROTEIN:
        return identified_peptide_keys_[&(*molecule_var.getIdentifiedPeptideRef())];
      case ID::MoleculeType::COMPOUND:
        return identified_compound_keys_[&(*molecule_var.getIdentifiedCompoundRef())];
      case ID::MoleculeType::RNA:
        return identified_oligo_keys_[&(*molecule_var.getIdentifiedOligoRef())];
      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
  }


  void OMSFileStore::storeObservationMatches_(const IdentificationData& id_data)
  {
    if (id_data.getObservationMatches().empty()) return;

    String table_def =
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "identified_molecule_id INTEGER NOT NULL, "                       \
      "observation_id INTEGER NOT NULL, "                               \
      "adduct_id INTEGER, "                                             \
      "charge INTEGER, "                                                \
      "FOREIGN KEY (identified_molecule_id) REFERENCES ID_IdentifiedMolecule (id), " \
      "FOREIGN KEY (observation_id) REFERENCES ID_Observation (id)";
    // add foreign key constraint if the adduct table exists (having the
    // constraint without the table would cause an error on data insertion):
    if (db_->tableExists("AdductInfo"))
    {
      table_def += ", FOREIGN KEY (adduct_id) REFERENCES AdductInfo (id)";
    }
    createTable_("ID_ObservationMatch", table_def);

    SQLite::Statement query(*db_, "INSERT INTO ID_ObservationMatch VALUES ("    \
                  ":id, "                                       \
                  ":identified_molecule_id, "                   \
                  ":observation_id, "                           \
                  ":adduct_id, "                                \
                  ":charge)");
    bool any_peak_annotations = false;
    Key id = 1;
    for (const ID::ObservationMatch& match : id_data.getObservationMatches())
    {
      if (!match.peak_annotations.empty()) any_peak_annotations = true;
      query.bind(":id", id);
      query.bind(":identified_molecule_id", getDatabaseKey_(match.identified_molecule_var));
      query.bind(":observation_id", observation_keys_[&(*match.observation_ref)]);
      if (match.adduct_opt)
      {
        query.bind(":adduct_id", adduct_keys_[&(**match.adduct_opt)]);
      }
      else // bind NULL value
      {
        query.bind(":adduct_id");
      }
      query.bind(":charge", match.charge);
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      observation_match_keys_[&match] = id;
      ++id;
    }
    storeScoredProcessingResults_(id_data.getObservationMatches(), "ID_ObservationMatch", observation_match_keys_);

    if (any_peak_annotations)
    {
      createTable_(
        "ID_ObservationMatch_PeakAnnotation",
        "parent_id INTEGER NOT NULL, "                                  \
        "processing_step_id INTEGER, "                                  \
        "peak_annotation TEXT, "                                        \
        "peak_charge INTEGER, "                                         \
        "peak_mz REAL, "                                                \
        "peak_intensity REAL, "                                         \
        "FOREIGN KEY (parent_id) REFERENCES ID_ObservationMatch (id), " \
        "FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)");

      SQLite::Statement query2(*db_,
        "INSERT INTO ID_ObservationMatch_PeakAnnotation VALUES ("   \
        ":parent_id, "                                              \
        ":processing_step_id, "                                     \
        ":peak_annotation, "                                        \
        ":peak_charge, "                                            \
        ":peak_mz, "                                                \
        ":peak_intensity)");

      for (const ID::ObservationMatch& match : id_data.getObservationMatches())
      {
        if (match.peak_annotations.empty()) continue;
        query2.bind(":parent_id", observation_match_keys_[&match]);
        for (const auto& pair : match.peak_annotations)
        {
          if (pair.first) // processing step given
          {
            query2.bind(":processing_step_id", processing_step_keys_[&(**pair.first)]);
          }
          else // use NULL value
          {
            query2.bind(":processing_step_id");
          }
          for (const auto& peak_ann : pair.second)
          {
            query2.bind(":peak_annotation", peak_ann.annotation);
            query2.bind(":peak_charge", peak_ann.charge);
            query2.bind(":peak_mz", peak_ann.mz);
            query2.bind(":peak_intensity", peak_ann.intensity);
            execWithExceptionAndReset(query2, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
          }
        }
      }
      // create index on parent_id column
      db_->exec("CREATE INDEX PeakAnnotation_parent_id ON ID_ObservationMatch_PeakAnnotation (parent_id)");
    }
  }


  void OMSFileStore::store(const IdentificationData& id_data)
  {
    startProgress(0, 13, "Writing identification data to file");
    // generally, create tables only if we have data to write - no empty ones!
    auto body = [&]() {
      storeVersionAndDate_();
      nextProgress(); // 1
      storeInputFiles_(id_data);
      nextProgress(); // 2
      storeScoreTypes_(id_data);
      nextProgress(); // 3
      storeProcessingSoftwares_(id_data);
      nextProgress(); // 4
      storeDBSearchParams_(id_data);
      nextProgress(); // 5
      storeProcessingSteps_(id_data);
      nextProgress(); // 6
      storeObservations_(id_data);
      nextProgress(); // 7
      storeParentSequences_(id_data);
      nextProgress(); // 8
      storeParentGroupSets_(id_data);
      nextProgress(); // 9
      storeIdentifiedCompounds_(id_data);
      nextProgress(); // 10
      storeIdentifiedSequences_(id_data);
      nextProgress(); // 11
      storeAdducts_(id_data);
      nextProgress(); // 12
      storeObservationMatches_(id_data);
    };

    if (sqlite3_get_autocommit(db_->getHandle()) == 1)
    { // allow a transaction, otherwise another on is already in flight
      SQLite::Transaction transaction(*db_); // avoid SQLite's "implicit transactions", improve runtime
      body();
      transaction.commit();
    }
    else
    {
      body();
    }
    endProgress();
    // @TODO: store input match groups
  }


  void OMSFileStore::createTableBaseFeature_(bool with_metainfo, bool with_idmatches)
  {
    createTable_("FEAT_BaseFeature",
                 "id INTEGER PRIMARY KEY NOT NULL, "                    \
                 "rt REAL, "                                            \
                 "mz REAL, "                                            \
                 "intensity REAL, "                                     \
                 "charge INTEGER, "                                     \
                 "width REAL, "                                         \
                 "quality REAL, "                                       \
                 "unique_id INTEGER, "                                  \
                 "primary_molecule_id INTEGER, "                        \
                 "subordinate_of INTEGER, "                             \
                 "FOREIGN KEY (primary_molecule_id) REFERENCES ID_IdentifiedMolecule (id), " \
                 "FOREIGN KEY (subordinate_of) REFERENCES FEAT_BaseFeature (id), " \
                 "CHECK (id > subordinate_of)"); // check to prevent cycles

    auto query = make_unique<SQLite::Statement>(*db_,
                                                "INSERT INTO FEAT_BaseFeature VALUES (" \
                                                ":id, "                 \
                                                ":rt, "                 \
                                                ":mz, "                 \
                                                ":intensity, "          \
                                                ":charge, "             \
                                                ":width, "              \
                                                ":quality, "            \
                                                ":unique_id, "          \
                                                ":primary_molecule_id, " \
                                                ":subordinate_of)");
    prepared_queries_.emplace("FEAT_BaseFeature", std::move(query));

    if (with_metainfo)
    {
      createTableMetaInfo_("FEAT_BaseFeature");
    }
    if (with_idmatches)
    {
      createTable_("FEAT_ObservationMatch",
                   "feature_id INTEGER NOT NULL, "                      \
                   "observation_match_id INTEGER NOT NULL, "            \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id), " \
                   "FOREIGN KEY (observation_match_id) REFERENCES ID_ObservationMatch (id)");
      query = make_unique<SQLite::Statement>(*db_,
                                             "INSERT INTO FEAT_ObservationMatch VALUES (" \
                                             ":feature_id, "            \
                                             ":observation_match_id)");
      prepared_queries_.emplace("FEAT_ObservationMatch", std::move(query));
    }
  }


  void OMSFileStore::storeBaseFeature_(const BaseFeature& feature, int feature_id, int parent_id)
  {
    auto& query_feat = *prepared_queries_["FEAT_BaseFeature"];
    query_feat.bind(":id", feature_id);
    query_feat.bind(":rt", feature.getRT());
    query_feat.bind(":mz", feature.getMZ());
    query_feat.bind(":intensity", feature.getIntensity());
    query_feat.bind(":charge", feature.getCharge());
    query_feat.bind(":width", feature.getWidth());
    query_feat.bind(":quality", feature.getQuality());
    query_feat.bind(":unique_id", int64_t(feature.getUniqueId()));
    if (feature.hasPrimaryID())
    {
      query_feat.bind(":primary_molecule_id", getDatabaseKey_(feature.getPrimaryID()));
    }
    else // use NULL value
    {
      query_feat.bind(":primary_molecule_id");
    }
    if (parent_id >= 0) // feature is a subordinate
    {
      query_feat.bind(":subordinate_of", parent_id);
    }
    else // use NULL value
    {
      query_feat.bind(":subordinate_of");
    }
    execWithExceptionAndReset(query_feat, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    // store ID observation matches:
    if (!feature.getIDMatches().empty())
    {
      auto& query_match = *prepared_queries_["FEAT_ObservationMatch"];
      query_match.bind(":feature_id", feature_id);
      for (ID::ObservationMatchRef ref : feature.getIDMatches())
      {
        query_match.bind(":observation_match_id", observation_match_keys_[&(*ref)]);
        execWithExceptionAndReset(query_match, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
    }

    storeMetaInfo_(feature, "FEAT_BaseFeature", feature_id);
  }


  void OMSFileStore::storeFeatureAndSubordinates_(
    const Feature& feature, int& feature_id, int parent_id)
  {
    storeBaseFeature_(feature, feature_id, parent_id);

    auto& query_feat = *prepared_queries_["FEAT_Feature"];
    query_feat.bind(":feature_id", feature_id);
    query_feat.bind(":rt_quality", feature.getQuality(0));
    query_feat.bind(":mz_quality", feature.getQuality(1));
    execWithExceptionAndReset(query_feat, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    // store convex hulls:
    const vector<ConvexHull2D>& hulls = feature.getConvexHulls();
    if (!hulls.empty())
    {
      auto& query_hull = *prepared_queries_["FEAT_ConvexHull"];
      query_hull.bind(":feature_id", feature_id);
      for (size_t i = 0; i < hulls.size(); ++i)
      {
        query_hull.bind(":hull_index", (int64_t)i);
        for (size_t j = 0; j < hulls[i].getHullPoints().size(); ++j)
        {
          const ConvexHull2D::PointType& point = hulls[i].getHullPoints()[j];
          query_hull.bind(":point_index", (int64_t)j);
          query_hull.bind(":point_x", point.getX());
          query_hull.bind(":point_y", point.getY());
          execWithExceptionAndReset(query_hull, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
        }

      }
    }
    // recurse into subordinates:
    parent_id = feature_id;
    ++feature_id; // variable is passed by reference, so effect is global
    for (const Feature& sub : feature.getSubordinates())
    {
      storeFeatureAndSubordinates_(sub, feature_id, parent_id);
    }
  }


  void OMSFileStore::storeFeatures_(const FeatureMap& features)
  {
    if (features.empty()) return;

    // create table(s) for BaseFeature parent class:
    // any meta infos on features?
    bool any_metainfo = anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.isMetaEmpty();
    });
    // any ID observations on features?
    bool any_idmatches = anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.getIDMatches().empty();
    });
    createTableBaseFeature_(any_metainfo, any_idmatches);

    createTable_("FEAT_Feature",
                 "feature_id INTEGER NOT NULL, "                        \
                 "rt_quality REAL, "                                    \
                 "mz_quality REAL, "                                    \
                 "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
    auto query = make_unique<SQLite::Statement>(*db_,
                                                "INSERT INTO FEAT_Feature VALUES (" \
                                                ":feature_id, "         \
                                                ":rt_quality, "         \
                                                ":mz_quality)");
    prepared_queries_.emplace("FEAT_Feature", std::move(query));

    // any convex hulls on features?
    if (anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.getConvexHulls().empty();
    }))
    {
      createTable_("FEAT_ConvexHull",
                   "feature_id INTEGER NOT NULL, "                      \
                   "hull_index INTEGER NOT NULL CHECK (hull_index >= 0), " \
                   "point_index INTEGER NOT NULL CHECK (point_index >= 0), " \
                   "point_x REAL, "                                     \
                   "point_y REAL, "                                     \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
      auto query2 = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_ConvexHull VALUES (" \
                                                   ":feature_id, "      \
                                                   ":hull_index, "      \
                                                   ":point_index, "     \
                                                   ":point_x, "         \
                                                   ":point_y)");
      prepared_queries_.emplace("FEAT_ConvexHull", std::move(query2));
    }

    // features and their subordinates are stored in DFS-like order:
    int feature_id = 0;
    for (const Feature& feat : features)
    {
      storeFeatureAndSubordinates_(feat, feature_id, -1);
      nextProgress();
    }
  }


  template <class MapType>
  void OMSFileStore::storeMapMetaData_(const MapType& features,
                                       const String& experiment_type)
  {
    createTable_("FEAT_MapMetaData",
                 "unique_id INTEGER PRIMARY KEY, "  \
                 "identifier TEXT, "                \
                 "file_path TEXT, "                 \
                 "file_type TEXT, "
                 "experiment_type TEXT"); // ConsensusMap only
    // @TODO: worth using a prepared query for just one insert?
    SQLite::Statement query(*db_,
                            "INSERT INTO FEAT_MapMetaData VALUES (" \
                            ":unique_id, "                          \
                            ":identifier, "                         \
                            ":file_path, "                          \
                            ":file_type, "                          \
                            ":experiment_type)");
    query.bind(":unique_id", int64_t(features.getUniqueId()));
    query.bind(":identifier", features.getIdentifier());
    query.bind(":file_path", features.getLoadedFilePath());
    String file_type = FileTypes::typeToName(features.getLoadedFileType());
    query.bind(":file_type", file_type);
    if (!experiment_type.empty())
    {
      query.bind(":experiment_type", experiment_type);
    }

    execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

    if (!features.isMetaEmpty())
    {
      createTableMetaInfo_("FEAT_MapMetaData", "unique_id");
      storeMetaInfo_(features, "FEAT_MapMetaData", int64_t(features.getUniqueId()));
    }
  }

  // template specializations:
  template void OMSFileStore::storeMapMetaData_<FeatureMap>(const FeatureMap&, const String&);
  template void OMSFileStore::storeMapMetaData_<ConsensusMap>(const ConsensusMap&, const String&);


  void OMSFileStore::storeDataProcessing_(const vector<DataProcessing>& data_processing)
  {
    if (data_processing.empty()) return;

    createTable_("FEAT_DataProcessing",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "software_name TEXT, "                 \
                 "software_version TEXT, "              \
                 "processing_actions TEXT, "            \
                 "completion_time TEXT");
    // "id" is needed to connect to meta info table (see "storeMetaInfos_");
    // "position" is position in the vector ("index" is a reserved word in SQL)
    SQLite::Statement query(*db_, "INSERT INTO FEAT_DataProcessing VALUES (" \
                            ":id, "                                     \
                            ":software_name, "                          \
                            ":software_version, "                       \
                            ":processing_actions, "                     \
                            ":completion_time)");

    Key id = 1;
    for (const DataProcessing& proc : data_processing)
    {
      query.bind(":id", id);
      query.bind(":software_name", proc.getSoftware().getName());
      query.bind(":software_version", proc.getSoftware().getVersion());
      String actions;
      for (DataProcessing::ProcessingAction action : proc.getProcessingActions())
      {
        if (!actions.empty()) actions += ","; // @TODO: use different separator?
        actions += DataProcessing::NamesOfProcessingAction[action];
      }
      query.bind(":processing_actions", actions);
      query.bind(":completion_time", proc.getCompletionTime().get());
      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      feat_processing_keys_[&proc] = id;
      ++id;
    }
    storeMetaInfos_(data_processing, "FEAT_DataProcessing", feat_processing_keys_);
  }


  void OMSFileStore::store(const FeatureMap& features)
  {
    SQLite::Transaction transaction(*db_); // avoid SQLite's "implicit transactions", improve runtime
    if (features.getIdentificationData().empty())
    {
      storeVersionAndDate_();
    }
    else
    {
      store(features.getIdentificationData());
    }
    startProgress(0, features.size() + 2, "Writing feature data to file");
    storeMapMetaData_(features);
    nextProgress();
    storeDataProcessing_(features.getDataProcessing());
    nextProgress();
    storeFeatures_(features);
    transaction.commit();
    endProgress();
  }


  void OMSFileStore::storeConsensusColumnHeaders_(const ConsensusMap& consensus)
  {
    if (consensus.getColumnHeaders().empty()) return; // shouldn't be empty in practice

    createTable_("FEAT_ConsensusColumnHeader",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "filename TEXT, "                      \
                 "label TEXT, "                         \
                 "size INTEGER, "                       \
                 "unique_id INTEGER");
    if (any_of(consensus.getColumnHeaders().begin(), consensus.getColumnHeaders().end(),
               [](const auto& pair){
                 return !pair.second.isMetaEmpty();
               }))
    {
      createTableMetaInfo_("FEAT_ConsensusColumnHeader");
    }

    SQLite::Statement query(*db_,
                            "INSERT INTO FEAT_ConsensusColumnHeader VALUES (" \
                            ":id, "                                     \
                            ":filename, "                               \
                            ":label, "                                  \
                            ":size, "                                   \
                            ":unique_id)");
    for (const auto& pair : consensus.getColumnHeaders())
    {
      Key id = int64_t(pair.first);
      query.bind(":id", id);
      query.bind(":filename", pair.second.filename);
      query.bind(":label", pair.second.label);
      query.bind(":size", int64_t(pair.second.size));
      query.bind(":unique_id", int64_t(pair.second.unique_id));

      execWithExceptionAndReset(query, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");

      storeMetaInfo_(pair.second, "FEAT_ConsensusColumnHeader", id);
    }
  }


  void OMSFileStore::storeConsensusFeatures_(const ConsensusMap& consensus)
  {
    if (consensus.empty()) return;

    // create table(s) for BaseFeature parent class:
    // any meta infos on features?
    bool any_metainfo = any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.isMetaEmpty();
    });
    // any ID observations on features?
    bool any_idmatches = any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.getIDMatches().empty();
    });
    createTableBaseFeature_(any_metainfo, any_idmatches);

    createTable_("FEAT_FeatureHandle",
                 "feature_id INTEGER NOT NULL, "                        \
                 "map_index INTEGER NOT NULL, "                                    \
                 "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
    SQLite::Statement query_handle(*db_,
                                   "INSERT INTO FEAT_FeatureHandle VALUES (" \
                                   ":feature_id, "                      \
                                   ":map_index)");

    // any ratios on consensus features?
    unique_ptr<SQLite::Statement> query_ratio; // only assign if needed below
    if (any_of(consensus.begin(), consensus.end(), [](const ConsensusFeature& feature) {
      return !feature.getRatios().empty();
    }))
    {
      createTable_("FEAT_ConsensusRatio",
                   "feature_id INTEGER NOT NULL, "                      \
                   "ratio_index INTEGER NOT NULL CHECK (ratio_index >= 0), " \
                   "ratio_value REAL, "                                 \
                   "denominator_ref TEXT, "                             \
                   "numerator_ref TEXT, "                               \
                   "description TEXT, "                                 \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_BaseFeature (id)");
      query_ratio = make_unique<SQLite::Statement>(*db_, "INSERT INTO FEAT_ConsensusRatio VALUES (" \
                                                   ":feature_id, "      \
                                                   ":ratio_index, "     \
                                                   ":ratio_value, "     \
                                                   ":denominator_ref, " \
                                                   ":numerator_ref, "   \
                                                   ":description)");
    }

    // consensus features and their subfeatures are stored in DFS-like order:
    int feature_id = 0;
    for (const ConsensusFeature& feat : consensus)
    {
      storeBaseFeature_(feat, feature_id, -1);
      int parent_id = feature_id;
      for (const FeatureHandle& handle : feat.getFeatures())
      {
        storeBaseFeature_(BaseFeature(handle), ++feature_id, parent_id);
        query_handle.bind(":feature_id", feature_id);
        query_handle.bind(":map_index", int64_t(handle.getMapIndex()));
        execWithExceptionAndReset(query_handle, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
      for (uint32_t i = 0; i < feat.getRatios().size(); ++i)
      {
        const ConsensusFeature::Ratio& ratio = feat.getRatios()[i];
        query_ratio->bind(":feature_id", feature_id);
        query_ratio->bind(":ratio_index", i);
        query_ratio->bind(":ratio_value", ratio.ratio_value_);
        query_ratio->bind(":denominator_ref", ratio.denominator_ref_);
        query_ratio->bind(":numerator_ref", ratio.numerator_ref_);
        query_ratio->bind(":description", ListUtils::concatenate(ratio.description_, ","));
        execWithExceptionAndReset(*query_ratio, 1, __LINE__, OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
      nextProgress();
      ++feature_id;
    }
  }


  void OMSFileStore::store(const ConsensusMap& consensus)
  {
    SQLite::Transaction transaction(*db_); // avoid SQLite's "implicit transactions", improve runtime
    if (consensus.getIdentificationData().empty())
    {
      storeVersionAndDate_();
    }
    else
    {
      store(consensus.getIdentificationData());
    }
    startProgress(0, consensus.size() + 3, "Writing consensus feature data to file");
    storeMapMetaData_(consensus, consensus.getExperimentType());
    nextProgress();
    storeConsensusColumnHeaders_(consensus);
    nextProgress();
    storeDataProcessing_(consensus.getDataProcessing());
    nextProgress();
    storeConsensusFeatures_(consensus);
    transaction.commit();
    endProgress();
  }
}
