// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFileStore.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlError>
// strangely, this is needed for type conversions in "QSqlQuery::bindValue":
#include <QtSql/QSqlQueryModel>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS::Internal
{
  int version_number = 1;

  void raiseDBError_(const QSqlError& error, int line,
                     const char* function, const String& context)
  {
    String msg = context + ": " + error.text();
    throw Exception::FailedAPICall(__FILE__, line, function, msg);
  }


  bool tableExists_(const String& db_name, const String& name)
  {
    QSqlDatabase db = QSqlDatabase::database(db_name.toQString());
    return db.tables(QSql::Tables).contains(name.toQString());
  }


  OMSFileStore::OMSFileStore(const String& filename, LogType log_type):
    db_name_("store_" + filename.toQString() + "_" +
             QString::number(UniqueIdGenerator::getUniqueId()))
  {
    setLogType(log_type);

    // delete output file if present:
    File::remove(filename);

    // open database:
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", db_name_);
    db.setDatabaseName(filename.toQString());
    if (!db.open())
    {
      raiseDBError_(db.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error opening SQLite database");
      // if d'tor doesn't get called, DB connection (db_name_) doesn't get
      // removed, but that shouldn't be a big problem
    }

    // configure database settings:
    QSqlQuery query(db);
    // foreign key constraints are disabled by default - turn them on:
    // @TODO: performance impact? (seems negligible, but should be tested more)
    if (!query.exec("PRAGMA foreign_keys = ON"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error configuring database");
    }
    // disable synchronous filesystem access and the rollback journal to greatly
    // increase write performance - since we write a new output file every time,
    // we don't have to worry about database consistency:
    if (!query.exec("PRAGMA synchronous = OFF"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error configuring database");
    }
    if (!query.exec("PRAGMA journal_mode = OFF"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error configuring database");
    }
  }


  OMSFileStore::~OMSFileStore()
  {
    QSqlDatabase::database(db_name_).close();
    QSqlDatabase::removeDatabase(db_name_);
  }


  void OMSFileStore::createTable_(const String& name,
                                  const String& definition,
                                  bool may_exist)
  {
    QString sql_create = "CREATE TABLE ";
    if (may_exist) sql_create += "IF NOT EXISTS ";
    sql_create += name.toQString() + " (" + definition.toQString() + ")";
    QSqlQuery query(QSqlDatabase::database(db_name_));
    if (!query.exec(sql_create))
    {
      String msg = "error creating database table " + name;
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION, msg);
    }
  }


  void OMSFileStore::storeVersionAndDate_()
  {
    createTable_("version",
                 "OMSFile INT NOT NULL, "       \
                 "date TEXT NOT NULL, "         \
                 "OpenMS TEXT, "                \
                 "build_date TEXT");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO version VALUES ("  \
                  ":format_version, "             \
                  "datetime('now'), "             \
                  ":openms_version, "             \
                  ":build_date)");
    query.bindValue(":format_version", version_number);
    query.bindValue(":openms_version", VersionInfo::getVersion().toQString());
    query.bindValue(":build_date", VersionInfo::getTime().toQString());
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
  }


  void OMSFileStore::createTableMoleculeType_()
  {
    createTable_("ID_MoleculeType",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "molecule_type TEXT UNIQUE NOT NULL");
    QString sql_insert =
      "INSERT INTO ID_MoleculeType VALUES "     \
      "(1, 'PROTEIN'), "                        \
      "(2, 'COMPOUND'), "                       \
      "(3, 'RNA')";
    QSqlQuery query(QSqlDatabase::database(db_name_));
    if (!query.exec(sql_insert))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
  }


  void OMSFileStore::createTableDataValue_()
  {
    createTable_("DataValue_DataType",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "data_type TEXT UNIQUE NOT NULL");
    QString sql_insert =
      "INSERT INTO DataValue_DataType VALUES " \
      "(1, 'STRING_VALUE'), "                  \
      "(2, 'INT_VALUE'), "                     \
      "(3, 'DOUBLE_VALUE'), "                  \
      "(4, 'STRING_LIST'), "                   \
      "(5, 'INT_LIST'), "                      \
      "(6, 'DOUBLE_LIST')";
    QSqlQuery query(QSqlDatabase::database(db_name_));
    if (!query.exec(sql_insert))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
    createTable_(
      "DataValue",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "data_type_id INTEGER, "                                          \
      "value TEXT, "                                                    \
      "FOREIGN KEY (data_type_id) REFERENCES DataValue_DataType (id)");
    // @TODO: add support for units
    // prepare query for inserting data:
    query.prepare("INSERT INTO DataValue VALUES ("           \
                  "NULL, "                                   \
                  ":data_type, "                             \
                  ":value)");
    prepared_queries_["DataValue"] = query;
  }


  OMSFileStore::Key OMSFileStore::storeDataValue_(const DataValue& value)
  {
    // this assumes the "DataValue" table exists already!
    // @TODO: split this up and make several tables for different types?
    QSqlQuery& query = prepared_queries_["DataValue"];
    if (value.isEmpty()) // use NULL as the type for empty values
    {
      query.bindValue(":data_type", QVariant(QVariant::Int));
    }
    else
    {
      query.bindValue(":data_type", int(value.valueType()) + 1);
    }
    query.bindValue(":value", value.toQString());
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
    return query.lastInsertId().toLongLong();
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
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT OR IGNORE INTO CVTerm VALUES ("   \
                  "NULL, "                                  \
                  ":accession, "                            \
                  ":name, "                                 \
                  ":cv_identifier_ref)");
    prepared_queries_["CVTerm"] = query;
    // alternative query if CVTerm already exists:
    query.prepare("SELECT id FROM CVTerm "                          \
                  "WHERE accession = :accession AND name = :name");
    prepared_queries_["CVTerm_2"] = query;
  }


  OMSFileStore::Key OMSFileStore::storeCVTerm_(const CVTerm& cv_term)
  {
    // this assumes the "CVTerm" table exists already!
    QSqlQuery& query = prepared_queries_["CVTerm"];
    if (cv_term.getAccession().empty()) // use NULL for empty accessions
    {
      query.bindValue(":accession", QVariant(QVariant::String));
    }
    else
    {
      query.bindValue(":accession", cv_term.getAccession().toQString());
    }
    query.bindValue(":name", cv_term.getName().toQString());
    query.bindValue(":cv_identifier_ref",
                    cv_term.getCVIdentifierRef().toQString());
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error updating database");
    }
    if (query.lastInsertId().isValid())
    {
      return query.lastInsertId().toLongLong();
    }
    // else: insert has failed, record must already exist - get the key:
    QSqlQuery& alt_query = prepared_queries_["CVTerm_2"];
    if (cv_term.getAccession().empty()) // use NULL for empty accessions
    {
      alt_query.bindValue(":accession", QVariant(QVariant::String));
    }
    else
    {
      alt_query.bindValue(":accession", cv_term.getAccession().toQString());
    }
    alt_query.bindValue(":name", cv_term.getName().toQString());
    if (!alt_query.exec() || !alt_query.next())
    {
      raiseDBError_(alt_query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error querying database");
    }
    return alt_query.value(0).toLongLong();
  }


  void OMSFileStore::createTableMetaInfo_(const String& parent_table,
                                          const String& key_column)
  {
    if (!tableExists_(db_name_, "DataValue")) createTableDataValue_();

    String parent_ref = parent_table + " (" + key_column + ")";
    String table = parent_table + "_MetaInfo";
    createTable_(
      table,
      "parent_id INTEGER NOT NULL, "                            \
      "name TEXT NOT NULL, "                                    \
      "data_value_id INTEGER NOT NULL, "                        \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_ref + ", " \
      "FOREIGN KEY (data_value_id) REFERENCES DataValue (id), " \
      "PRIMARY KEY (parent_id, name)");

    // prepare query for inserting data:
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO " + table.toQString() + " VALUES ("  \
                  ":parent_id, "                                    \
                  ":name, "                                         \
                  ":data_value_id)");
    prepared_queries_[table] = query;
  }


  void OMSFileStore::storeMetaInfo_(const MetaInfoInterface& info,
                                    const String& parent_table,
                                    Key parent_id)
  {
    if (info.isMetaEmpty()) return;

    // this assumes the "..._MetaInfo" and "DataValue" tables exist already!
    QSqlQuery& query = prepared_queries_[parent_table + "_MetaInfo"];
    query.bindValue(":parent_id", parent_id);
    // this is inefficient, but MetaInfoInterface doesn't support iteration:
    vector<String> info_keys;
    info.getKeys(info_keys);
    for (const String& info_key : info_keys)
    {
      query.bindValue(":name", info_key.toQString());
      Key value_id = storeDataValue_(info.getMetaValue(info_key));
      query.bindValue(":data_value_id", value_id);
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO " + table.toQString() + " VALUES ("  \
                  ":parent_id, "                                    \
                  ":processing_step_id, "                           \
                  ":processing_step_order, "                        \
                  ":score_type_id, "                                \
                  ":score)");
    prepared_queries_[table] = query;
  }


  void OMSFileStore::storeAppliedProcessingStep_(
    const ID::AppliedProcessingStep& step, Size step_order,
    const String& parent_table, Key parent_id)
  {
    // this assumes the "..._AppliedProcessingStep" table exists already!
    QSqlQuery& query = prepared_queries_[parent_table + "_AppliedProcessingStep"];
    query.bindValue(":parent_id", parent_id);
    query.bindValue(":processing_step_order", int(step_order));
    if (step.processing_step_opt)
    {
      query.bindValue(":processing_step_id",
                      Key(&(**step.processing_step_opt)));
      if (step.scores.empty()) // insert processing step information only
      {
        query.bindValue(":score_type_id", QVariant(QVariant::Int)); // NULL
        query.bindValue(":score", QVariant(QVariant::Double)); // NULL
        if (!query.exec())
        {
          raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                        "error inserting data");
        }
      }
    }
    else // use NULL for missing processing step reference
    {
      query.bindValue(":processing_step_id", QVariant(QVariant::Int));
    }
    for (const auto& score_pair : step.scores)
    {
      query.bindValue(":score_type_id", Key(&(*score_pair.first)));
      query.bindValue(":score", score_pair.second);
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ScoreType VALUES ("       \
                  ":id, "                                   \
                  ":cv_term_id, "                           \
                  ":higher_better)");
    for (const ID::ScoreType& score_type : id_data.getScoreTypes())
    {
      Key cv_id = storeCVTerm_(score_type.cv_term);
      query.bindValue(":id", Key(&score_type)); // use address as primary key
      query.bindValue(":cv_term_id", cv_id);
      query.bindValue(":higher_better", int(score_type.higher_better));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_InputFile VALUES ("  \
                  ":id, "                              \
                  ":name, "                            \
                  ":experimental_design_id, "          \
                  ":primary_files)");
    for (const ID::InputFile& input : id_data.getInputFiles())
    {
      query.bindValue(":id", Key(&input));
      query.bindValue(":name", input.name.toQString());
      query.bindValue(":experimental_design_id",
                      input.experimental_design_id.toQString());
      // @TODO: what if a primary file name contains ","?
      String primary_files = ListUtils::concatenate(input.primary_files);
      query.bindValue(":primary_files", primary_files.toQString());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ProcessingSoftware VALUES ("  \
                  ":id, "                                           \
                  ":name, "                                         \
                  ":version)");
    bool any_scores = false; // does any software have assigned scores stored?
    for (const ID::ProcessingSoftware& software : id_data.getProcessingSoftwares())
    {
      if (!software.assigned_scores.empty()) any_scores = true;
      query.bindValue(":id", Key(&software));
      query.bindValue(":name", software.getName().toQString());
      query.bindValue(":version", software.getVersion().toQString());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

      query.prepare(
        "INSERT INTO ID_ProcessingSoftware_AssignedScore VALUES ("      \
        ":software_id, "                                                \
        ":score_type_id, "                                              \
        ":score_type_order)");
      for (const ID::ProcessingSoftware& software : id_data.getProcessingSoftwares())
      {
        query.bindValue(":software_id", Key(&software));
        Size counter = 0;
        for (ID::ScoreTypeRef score_type_ref : software.assigned_scores)
        {
          query.bindValue(":score_type_id", Key(&(*score_type_ref)));
          query.bindValue(":score_type_order", int(++counter));
          if (!query.exec())
          {
            raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                          "error inserting data");
          }
        }
      }
    }
  }


  void OMSFileStore::storeDBSearchParams_(const IdentificationData& id_data)
  {
    if (id_data.getDBSearchParams().empty()) return;

    if (!tableExists_(db_name_, "ID_MoleculeType")) createTableMoleculeType_();

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
      "missed_cleavages NUMERIC, "                                      \
      "min_length NUMERIC, "                                            \
      "max_length NUMERIC, "                                            \
      "FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_DBSearchParam VALUES (" \
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
                  ":missed_cleavages, "                   \
                  ":min_length, "                         \
                  ":max_length)");
    for (const ID::DBSearchParam& param : id_data.getDBSearchParams())
    {
      query.bindValue(":id", Key(&param));
      query.bindValue(":molecule_type_id", int(param.molecule_type) + 1);
      query.bindValue(":mass_type_average", int(param.mass_type));
      query.bindValue(":database", param.database.toQString());
      query.bindValue(":database_version", param.database_version.toQString());
      query.bindValue(":taxonomy", param.taxonomy.toQString());
      String charges = ListUtils::concatenate(param.charges, ",");
      query.bindValue(":charges", charges.toQString());
      String fixed_mods = ListUtils::concatenate(param.fixed_mods, ",");
      query.bindValue(":fixed_mods", fixed_mods.toQString());
      String variable_mods = ListUtils::concatenate(param.variable_mods, ",");
      query.bindValue(":variable_mods", variable_mods.toQString());
      query.bindValue(":precursor_mass_tolerance",
                      param.precursor_mass_tolerance);
      query.bindValue(":fragment_mass_tolerance",
                      param.fragment_mass_tolerance);
      query.bindValue(":precursor_tolerance_ppm",
                      int(param.precursor_tolerance_ppm));
      query.bindValue(":fragment_tolerance_ppm",
                      int(param.fragment_tolerance_ppm));
      if (param.digestion_enzyme != nullptr)
      {
        query.bindValue(":digestion_enzyme",
                        param.digestion_enzyme->getName().toQString());
      }
      else // bind NULL value
      {
        query.bindValue(":digestion_enzyme", QVariant(QVariant::String));
      }
      query.bindValue(":missed_cleavages", uint(param.missed_cleavages));
      query.bindValue(":min_length", uint(param.min_length));
      query.bindValue(":max_length", uint(param.max_length));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ProcessingStep VALUES ("  \
                  ":id, "                                       \
                  ":software_id, "                              \
                  ":date_time, "                                \
                  ":search_param_id)");
    bool any_input_files = false;
    // use iterator here because we need one to look up the DB search params:
    for (ID::ProcessingStepRef step_ref = id_data.getProcessingSteps().begin();
         step_ref != id_data.getProcessingSteps().end(); ++step_ref)
    {
      const ID::ProcessingStep& step = *step_ref;
      if (!step.input_file_refs.empty()) any_input_files = true;
      query.bindValue(":id", Key(&step));
      query.bindValue(":software_id", Key(&(*step.software_ref)));
      query.bindValue(":date_time", step.date_time.get().toQString());
      auto pos = id_data.getDBSearchSteps().find(step_ref);
      if (pos != id_data.getDBSearchSteps().end())
      {
        query.bindValue(":search_param_id", Key(&(*pos->second)));
      }
      else
      {
        query.bindValue(":search_param_id", QVariant(QVariant::Int)); // NULL
      }
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
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

      query.prepare("INSERT INTO ID_ProcessingStep_InputFile VALUES (" \
                    ":processing_step_id, "                             \
                    ":input_file_id)");

      for (const ID::ProcessingStep& step : id_data.getProcessingSteps())
      {
        query.bindValue(":processing_step_id", Key(&step));
        for (ID::InputFileRef input_file_ref : step.input_file_refs)
        {
          query.bindValue(":input_file_id", Key(&(*input_file_ref)));
          if (!query.exec())
          {
            raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                          "error inserting data");
          }
        }
      }
    }
    storeMetaInfos_(id_data.getProcessingSteps(), "ID_ProcessingStep");
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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_Observation VALUES (" \
                  ":id, "                             \
                  ":data_id, "                        \
                  ":input_file_id, "                  \
                  ":rt, "                             \
                  ":mz)");
    for (const ID::Observation& obs : id_data.getObservations())
    {
      query.bindValue(":id", Key(&obs)); // use address as primary key
      query.bindValue(":data_id", obs.data_id.toQString());
      query.bindValue(":input_file_id", Key(&(*obs.input_file)));

      if (obs.rt == obs.rt)
      {
        query.bindValue(":rt", obs.rt);
      }
      else // NaN
      {
        query.bindValue(":rt", QVariant(QVariant::Double)); // NULL
      }
      if (obs.mz == obs.mz)
      {
        query.bindValue(":mz", obs.mz);
      }
      else // NaN
      {
        query.bindValue(":mz", QVariant(QVariant::Double)); // NULL
      }
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
    storeMetaInfos_(id_data.getObservations(), "ID_Observation");
  }


  void OMSFileStore::storeParentSequences_(const IdentificationData& id_data)
  {
    if (id_data.getParentSequences().empty()) return;

    if (!tableExists_(db_name_, "ID_MoleculeType")) createTableMoleculeType_();

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

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ParentSequence VALUES ("  \
                  ":id, "                                   \
                  ":accession, "                            \
                  ":molecule_type_id, "                     \
                  ":sequence, "                             \
                  ":description, "                          \
                  ":coverage, "                             \
                  ":is_decoy)");
    for (const ID::ParentSequence& parent : id_data.getParentSequences())
    {
      query.bindValue(":id", Key(&parent)); // use address as primary key
      query.bindValue(":accession", parent.accession.toQString());
      query.bindValue(":molecule_type_id", int(parent.molecule_type) + 1);
      query.bindValue(":sequence", parent.sequence.toQString());
      query.bindValue(":description", parent.description.toQString());
      query.bindValue(":coverage", parent.coverage);
      query.bindValue(":is_decoy", int(parent.is_decoy));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
    storeScoredProcessingResults_(id_data.getParentSequences(), "ID_ParentSequence");
  }


  void OMSFileStore::storeParentGroupSets_(const IdentificationData& id_data)
  {
    if (id_data.getParentGroupSets().empty()) return;

    createTable_("ID_ParentGroupSet",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "label TEXT, "                       \
                 "grouping_order INTEGER NOT NULL");

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
      "FOREIGN KEY (group_id) REFERENCES ID_ParentGroup (id), " \
      "FOREIGN KEY (parent_id) REFERENCES ID_ParentSequence (id)");

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query_grouping(db);
    query_grouping.prepare("INSERT INTO ID_ParentGroupSet VALUES (" \
                           ":id, "                                      \
                           ":label, "                                   \
                           ":grouping_order)");

    QSqlQuery query_group(db);
    query_group.prepare("INSERT INTO ID_ParentGroup VALUES ("  \
                        ":id, "                                        \
                        ":grouping_id, "                               \
                        ":score_type_id, "                             \
                        ":score)");

    QSqlQuery query_parent(db);
    query_parent.prepare(
      "INSERT INTO ID_ParentGroup_ParentSequence VALUES ("   \
      ":group_id, "                                                  \
      ":parent_id)");

    Size counter = 0;
    for (const ID::ParentGroupSet& grouping : id_data.getParentGroupSets())
    {
      Key grouping_id = Key(&grouping);
      query_grouping.bindValue(":id", grouping_id);
      query_grouping.bindValue(":label", grouping.label.toQString());
      query_grouping.bindValue(":grouping_order", int(++counter));
      if (!query_grouping.exec())
      {
        raiseDBError_(query_grouping.lastError(), __LINE__,
                      OPENMS_PRETTY_FUNCTION, "error inserting data");
      }

      for (const ID::ParentGroup& group : grouping.groups)
      {
        Key group_id = Key(&group);
        query_group.bindValue(":id", group_id);
        query_group.bindValue(":grouping_id", grouping_id);
        if (group.scores.empty()) // store group with an empty score
        {
          query_group.bindValue(":score_type_id", QVariant(QVariant::Int));
          query_group.bindValue(":score", QVariant(QVariant::Double));
          if (!query_group.exec())
          {
            raiseDBError_(query_group.lastError(), __LINE__,
                          OPENMS_PRETTY_FUNCTION, "error inserting data");
          }
        }
        else // store group multiple times with different scores
        {
          for (const auto& score_pair : group.scores)
          {
            query_group.bindValue(":score_type_id", Key(&(*score_pair.first)));
            query_group.bindValue(":score", score_pair.second);
            if (!query_group.exec())
            {
              raiseDBError_(query_group.lastError(), __LINE__,
                            OPENMS_PRETTY_FUNCTION, "error inserting data");
            }
          }
        }

        query_parent.bindValue(":group_id", group_id);
        for (ID::ParentSequenceRef parent_ref : group.parent_refs)
        {
          query_parent.bindValue(":parent_id", Key(&(*parent_ref)));
          if (!query_parent.exec())
          {
            raiseDBError_(query_parent.lastError(), __LINE__,
                          OPENMS_PRETTY_FUNCTION, "error inserting data");
          }
        }
      }
    }

    storeScoredProcessingResults_(id_data.getParentGroupSets(), "ID_ParentGroupSet");
  }


  void OMSFileStore::createTableIdentifiedMolecule_()
  {
    if (!tableExists_(db_name_, "ID_MoleculeType")) createTableMoleculeType_();

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
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_IdentifiedMolecule VALUES ("          \
                  ":id, "                                               \
                  ":molecule_type_id, "                                 \
                  ":identifier)");
    prepared_queries_["ID_IdentifiedMolecule"] = query;
  }


  void OMSFileStore::storeIdentifiedCompounds_(const IdentificationData& id_data)
  {
    if (id_data.getIdentifiedCompounds().empty()) return;

    if (!tableExists_(db_name_, "ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }
    QSqlQuery& query_molecule = prepared_queries_["ID_IdentifiedMolecule"];
    query_molecule.bindValue(":molecule_type_id",
                             int(ID::MoleculeType::COMPOUND) + 1);

    createTable_(
      "ID_IdentifiedCompound",
      "molecule_id INTEGER UNIQUE NOT NULL , "                          \
      "formula TEXT, "                                                  \
      "name TEXT, "                                                     \
      "smile TEXT, "                                                    \
      "inchi TEXT, "                                                    \
      "FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)");
    QSqlQuery query_compound(QSqlDatabase::database(db_name_));
    query_compound.prepare("INSERT INTO ID_IdentifiedCompound VALUES (" \
                           ":molecule_id, "                             \
                           ":formula, "                                 \
                           ":name, "                                    \
                           ":smile, "                                   \
                           ":inchi)");
    for (const ID::IdentifiedCompound& compound : id_data.getIdentifiedCompounds())
    {
      // use address as primary key:
      query_molecule.bindValue(":id", Key(&compound));
      query_molecule.bindValue(":identifier", compound.identifier.toQString());
      if (!query_molecule.exec())
      {
        raiseDBError_(query_molecule.lastError(), __LINE__,
                      OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
      query_compound.bindValue(":molecule_id", Key(&compound));
      query_compound.bindValue(":formula",
                               compound.formula.toString().toQString());
      query_compound.bindValue(":name", compound.name.toQString());
      query_compound.bindValue(":smile", compound.name.toQString());
      query_compound.bindValue(":inchi", compound.inchi.toQString());
      if (!query_compound.exec())
      {
        raiseDBError_(query_compound.lastError(), __LINE__,
                      OPENMS_PRETTY_FUNCTION, "error inserting data");
      }
    }
    storeScoredProcessingResults_(id_data.getIdentifiedCompounds(),
                                  "ID_IdentifiedMolecule");
  }


  void OMSFileStore::storeIdentifiedSequences_(const IdentificationData& id_data)
  {
    if (id_data.getIdentifiedPeptides().empty() &&
        id_data.getIdentifiedOligos().empty()) return;

    if (!tableExists_(db_name_, "ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }
    QSqlQuery& query = prepared_queries_["ID_IdentifiedMolecule"];

    bool any_parent_matches = false;
    // store peptides:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    for (const ID::IdentifiedPeptide& peptide : id_data.getIdentifiedPeptides())
    {
      if (!peptide.parent_matches.empty()) any_parent_matches = true;
      query.bindValue(":id", Key(&peptide)); // use address as primary key
      query.bindValue(":identifier", peptide.sequence.toString().toQString());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
    storeScoredProcessingResults_(id_data.getIdentifiedPeptides(),
                                  "ID_IdentifiedMolecule");
    // store RNA oligos:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    for (const ID::IdentifiedOligo& oligo : id_data.getIdentifiedOligos())
    {
      if (!oligo.parent_matches.empty()) any_parent_matches = true;
      query.bindValue(":id", Key(&oligo)); // use address as primary key
      query.bindValue(":identifier",
                      QString::fromStdString(oligo.sequence.toString()));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
    storeScoredProcessingResults_(id_data.getIdentifiedOligos(),
                                  "ID_IdentifiedMolecule");

    if (any_parent_matches)
    {
      createTableParentMatches_();
      for (const ID::IdentifiedPeptide& peptide : id_data.getIdentifiedPeptides())
      {
        if (peptide.parent_matches.empty()) continue;
        storeParentMatches_(peptide.parent_matches, Key(&peptide));
      }
      for (const ID::IdentifiedOligo& oligo : id_data.getIdentifiedOligos())
      {
        if (oligo.parent_matches.empty()) continue;
        storeParentMatches_(oligo.parent_matches, Key(&oligo));
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
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ParentMatch VALUES ("         \
                  ":molecule_id, "                              \
                  ":parent_id, "                                \
                  ":start_pos, "                                \
                  ":end_pos, "                                  \
                  ":left_neighbor, "                            \
                  ":right_neighbor)");
    prepared_queries_["ID_ParentMatch"] = query;
  }


  void OMSFileStore::storeParentMatches_(const ID::ParentMatches& matches,
                                         Key molecule_id)
  {
    // this assumes the "ID_ParentMatch" table exists already!
    QSqlQuery& query = prepared_queries_["ID_ParentMatch"];
    // @TODO: cache the prepared query between function calls somehow?
    query.bindValue(":molecule_id", molecule_id);
    for (const auto& pair : matches)
    {
      query.bindValue(":parent_id", Key(&(*pair.first)));
      for (const auto& match : pair.second)
      {
        if (match.start_pos != ID::ParentMatch::UNKNOWN_POSITION)
        {
          query.bindValue(":start_pos", uint(match.start_pos));
        }
        else // use NULL value
        {
          query.bindValue(":start_pos", QVariant(QVariant::Int));
        }
        if (match.end_pos != ID::ParentMatch::UNKNOWN_POSITION)
        {
          query.bindValue(":end_pos", uint(match.end_pos));
        }
        else // use NULL value
        {
          query.bindValue(":end_pos", QVariant(QVariant::Int));
        }
        query.bindValue(":left_neighbor", match.left_neighbor.toQString());
        query.bindValue(":right_neighbor", match.right_neighbor.toQString());
        if (!query.exec())
        {
          raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                        "error inserting data");
        }
      }
    }
  }


  void OMSFileStore::storeAdducts_(const IdentificationData& id_data)
  {
    if (id_data.getAdducts().empty()) return;

    createTable_(
      "AdductInfo",
      "id INTEGER PRIMARY KEY NOT NULL, " \
      "name TEXT, " \
      "formula TEXT NOT NULL, " \
      "charge INTEGER NOT NULL, " \
      "mol_multiplier INTEGER NOT NULL CHECK (mol_multiplier > 0) DEFAULT 1, " \
      "UNIQUE (formula, charge)");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO AdductInfo VALUES (" \
                  ":id, "                           \
                  ":name, "                         \
                  ":formula, "                      \
                  ":charge, "                       \
                  ":mol_multiplier)");
    for (const AdductInfo& adduct : id_data.getAdducts())
    {
      query.bindValue(":id", Key(&adduct));
      query.bindValue(":name", adduct.getName().toQString());
      query.bindValue(":formula", adduct.getEmpiricalFormula().toString().toQString());
      query.bindValue(":charge", adduct.getCharge());
      query.bindValue(":mol_multiplier", adduct.getMolMultiplier());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
  }


  OMSFileStore::Key OMSFileStore::getAddress_(const ID::IdentifiedMolecule& molecule_var)
  {
    switch (molecule_var.getMoleculeType())
    {
      case ID::MoleculeType::PROTEIN:
        return Key(&(*molecule_var.getIdentifiedPeptideRef()));
      case ID::MoleculeType::COMPOUND:
        return Key(&(*molecule_var.getIdentifiedCompoundRef()));
      case ID::MoleculeType::RNA:
        return Key(&(*molecule_var.getIdentifiedOligoRef()));
      default:
        throw Exception::NotImplemented(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION);
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
    if (tableExists_(db_name_, "AdductInfo"))
    {
      table_def += ", FOREIGN KEY (adduct_id) REFERENCES AdductInfo (id)";
    }
    createTable_("ID_ObservationMatch", table_def);

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ObservationMatch VALUES ("    \
                  ":id, "                                       \
                  ":identified_molecule_id, "                   \
                  ":observation_id, "                           \
                  ":adduct_id, "                                \
                  ":charge)");
    bool any_peak_annotations = false;
    for (const ID::ObservationMatch& match : id_data.getObservationMatches())
    {
      if (!match.peak_annotations.empty()) any_peak_annotations = true;
      query.bindValue(":id", Key(&match)); // use address as primary key
      query.bindValue(":identified_molecule_id", getAddress_(match.identified_molecule_var));
      query.bindValue(":observation_id", Key(&(*match.observation_ref)));
      if (match.adduct_opt)
      {
        query.bindValue(":adduct_id", Key(&(**match.adduct_opt)));
      }
      else // bind NULL value
      {
        query.bindValue(":adduct_id", QVariant(QVariant::Int));
      }
      query.bindValue(":charge", match.charge);
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
    storeScoredProcessingResults_(id_data.getObservationMatches(), "ID_ObservationMatch");

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

      query.prepare(
        "INSERT INTO ID_ObservationMatch_PeakAnnotation VALUES (" \
        ":parent_id, "                                              \
        ":processing_step_id, "                                     \
        ":peak_annotation, "                                        \
        ":peak_charge, "                                            \
        ":peak_mz, "                                                \
        ":peak_intensity)");

      for (const ID::ObservationMatch& match : id_data.getObservationMatches())
      {
        if (match.peak_annotations.empty()) continue;
        query.bindValue(":parent_id", Key(&match));
        for (const auto& pair : match.peak_annotations)
        {
          if (pair.first) // processing step given
          {
            query.bindValue(":processing_step_id", Key(&(**pair.first)));
          }
          else // use NULL value
          {
            query.bindValue(":processing_step_id", QVariant(QVariant::Int));
          }
          for (const auto& peak_ann : pair.second)
          {
            query.bindValue(":peak_annotation",
                            peak_ann.annotation.toQString());
            query.bindValue(":peak_charge", peak_ann.charge);
            query.bindValue(":peak_mz", peak_ann.mz);
            query.bindValue(":peak_intensity", peak_ann.intensity);
            if (!query.exec())
            {
              raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                            "error inserting data");
            }
          }
        }
      }
    }
    // create index on parent_id column
    query.exec("CREATE INDEX PeakAnnotation_parent_id ON ID_ObservationMatch_PeakAnnotation (parent_id)");
  }


  void OMSFileStore::store(const IdentificationData& id_data)
  {
    QSqlDatabase db = QSqlDatabase::database(db_name_);
    startProgress(0, 13, "Writing identification data to file");
    // generally, create tables only if we have data to write - no empty ones!
    db.transaction(); // avoid SQLite's "implicit transactions", improve runtime
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
    db.commit();
    endProgress();
    // @TODO: store input match groups
  }


  void OMSFileStore::storeFeatureAndSubordinates_(
    const Feature& feature, int& feature_id, int parent_id)
  {
    QSqlQuery& query_feat = prepared_queries_["FEAT_Feature"];
    query_feat.bindValue(":id", feature_id);
    query_feat.bindValue(":rt", feature.getRT());
    query_feat.bindValue(":mz", feature.getMZ());
    query_feat.bindValue(":intensity", feature.getIntensity());
    query_feat.bindValue(":charge", feature.getCharge());
    query_feat.bindValue(":width", feature.getWidth());
    query_feat.bindValue(":overall_quality", feature.getOverallQuality());
    query_feat.bindValue(":rt_quality", feature.getQuality(0));
    query_feat.bindValue(":mz_quality", feature.getQuality(1));
    query_feat.bindValue(":unique_id", qint64(feature.getUniqueId()));
    if (feature.hasPrimaryID())
    {
      query_feat.bindValue(":primary_molecule_id", getAddress_(feature.getPrimaryID()));
    }
    else // use NULL value
    {
      query_feat.bindValue(":primary_molecule_id", QVariant(QVariant::Int));
    }
    if (parent_id >= 0) // feature is a subordinate
    {
      query_feat.bindValue(":subordinate_of", parent_id);
    }
    else // use NULL value
    {
      query_feat.bindValue(":subordinate_of", QVariant(QVariant::Int));
    }
    if (!query_feat.exec())
    {
      raiseDBError_(query_feat.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
    storeMetaInfo_(feature, "FEAT_Feature", feature_id);
    // store convex hulls:
    const vector<ConvexHull2D>& hulls = feature.getConvexHulls();
    if (!hulls.empty())
    {
      QSqlQuery& query_hull = prepared_queries_["FEAT_ConvexHull"];
      query_hull.bindValue(":feature_id", feature_id);
      for (uint i = 0; i < hulls.size(); ++i)
      {
        query_hull.bindValue(":hull_index", i);
        for (uint j = 0; j < hulls[i].getHullPoints().size(); ++j)
        {
          const ConvexHull2D::PointType& point = hulls[i].getHullPoints()[j];
          query_hull.bindValue(":point_index", j);
          query_hull.bindValue(":point_x", point.getX());
          query_hull.bindValue(":point_y", point.getY());
          if (!query_hull.exec())
          {
            raiseDBError_(query_hull.lastError(), __LINE__,
                          OPENMS_PRETTY_FUNCTION, "error inserting data");
          }
        }

      }
    }
    // store ID input items:
    if (!feature.getIDMatches().empty())
    {
      QSqlQuery& query_match = prepared_queries_["FEAT_ObservationMatch"];
      query_match.bindValue(":feature_id", feature_id);
      for (ID::ObservationMatchRef ref : feature.getIDMatches())
      {
        query_match.bindValue(":observation_match_id", Key(&(*ref)));
        if (!query_match.exec())
        {
          raiseDBError_(query_match.lastError(), __LINE__,
                        OPENMS_PRETTY_FUNCTION, "error inserting data");
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

    createTable_("FEAT_Feature",
                 "id INTEGER PRIMARY KEY NOT NULL, "                    \
                 "rt REAL, "                                            \
                 "mz REAL, "                                            \
                 "intensity REAL, "                                     \
                 "charge INTEGER, "                                     \
                 "width REAL, "                                         \
                 "overall_quality REAL, "                               \
                 "rt_quality REAL, "                                    \
                 "mz_quality REAL, "                                    \
                 "unique_id INTEGER, "                                  \
                 "primary_molecule_id INTEGER, "                        \
                 "subordinate_of INTEGER, "                             \
                 "FOREIGN KEY (primary_molecule_id) REFERENCES ID_IdentifiedMolecule (id), " \
                 "FOREIGN KEY (subordinate_of) REFERENCES FEAT_Feature (id), " \
                 "CHECK (id > subordinate_of)"); // check to prevent cycles

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO FEAT_Feature VALUES ("      \
                  ":id, "                                  \
                  ":rt, "                                  \
                  ":mz, "                                  \
                  ":intensity, "                           \
                  ":charge, "                              \
                  ":width, "                               \
                  ":overall_quality, "                     \
                  ":rt_quality, "                          \
                  ":mz_quality, "                          \
                  ":unique_id, "                           \
                  ":primary_molecule_id, "                 \
                  ":subordinate_of)");
    prepared_queries_["FEAT_Feature"] = query;
    // any meta infos on features?
    if (anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.isMetaEmpty();
    }))
    {
      createTableMetaInfo_("FEAT_Feature");
    }
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
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_Feature (id)");
      query.prepare("INSERT INTO FEAT_ConvexHull VALUES ("      \
                    ":feature_id, "                             \
                    ":hull_index, "                             \
                    ":point_index, "                            \
                    ":point_x, "                                \
                    ":point_y)");
      prepared_queries_["FEAT_ConvexHull"] = query;
    }
    // any ID observations on features?
    if (anyFeaturePredicate_(features, [](const Feature& feature) {
      return !feature.getIDMatches().empty();
    }))
    {
      createTable_("FEAT_ObservationMatch",
                   "feature_id INTEGER NOT NULL, "                      \
                   "observation_match_id INTEGER NOT NULL, "            \
                   "FOREIGN KEY (feature_id) REFERENCES FEAT_Feature (id), " \
                   "FOREIGN KEY (observation_match_id) REFERENCES ID_ObservationMatch (id)");
      query.prepare("INSERT INTO FEAT_ObservationMatch VALUES (" \
                    ":feature_id, "                              \
                    ":observation_match_id)");
      prepared_queries_["FEAT_ObservationMatch"] = query;
    }

    // features and their subordinates are stored in DFS-like order:
    int feature_id = 0;
    for (const Feature& feat : features)
    {
      storeFeatureAndSubordinates_(feat, feature_id, -1);
      nextProgress();
    }
  }


  void OMSFileStore::storeMapMetaData_(const FeatureMap& features)
  {
    createTable_("FEAT_MapMetaData",
                 "unique_id INTEGER PRIMARY KEY, "  \
                 "identifier TEXT, "                \
                 "file_path TEXT, "                 \
                 "file_type TEXT");
    QSqlQuery query(QSqlDatabase::database(db_name_));
    // @TODO: worth using a prepared query for just one insert?
    query.prepare("INSERT INTO FEAT_MapMetaData VALUES (" \
                  ":unique_id, "                          \
                  ":identifier, "                         \
                  ":file_path, "                          \
                  ":file_type)");
    query.bindValue(":unique_id", qint64(features.getUniqueId()));
    query.bindValue(":identifier", features.getIdentifier().toQString());
    query.bindValue(":file_path", features.getLoadedFilePath().toQString());
    String file_type = FileTypes::typeToName(features.getLoadedFileType());
    query.bindValue(":file_type", file_type.toQString());

    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error inserting data");
    }
    if (!features.isMetaEmpty())
    {
      createTableMetaInfo_("FEAT_MapMetaData", "unique_id");
      storeMetaInfo_(features, "FEAT_MapMetaData", qint64(features.getUniqueId()));
    }
  }


  void OMSFileStore::storeDataProcessing_(const FeatureMap& features)
  {
    if (features.getDataProcessing().empty()) return;

    createTable_("FEAT_DataProcessing",
                 "id INTEGER PRIMARY KEY NOT NULL, "    \
                 "position INTEGER NOT NULL, "          \
                 "software_name TEXT, "                 \
                 "software_version TEXT, "              \
                 "processing_actions TEXT, "            \
                 "completion_time TEXT");
    // "id" is needed to connect to meta info table (see "storeMetaInfos_");
    // "position" is position in the vector ("index" is a reserved word in SQL)
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO FEAT_DataProcessing VALUES (" \
                  ":id, "                                    \
                  ":position, "                              \
                  ":software_name, "                         \
                  ":software_version, "                      \
                  ":processing_actions, "                    \
                  ":completion_time)");

    int index = 0;
    for (const DataProcessing& proc : features.getDataProcessing())
    {
      query.bindValue(":id", Key(&proc));
      query.bindValue(":position", index);
      query.bindValue(":software_name", proc.getSoftware().getName().toQString());
      query.bindValue(":software_version", proc.getSoftware().getVersion().toQString());
      String actions;
      for (DataProcessing::ProcessingAction action : proc.getProcessingActions())
      {
        if (!actions.empty()) actions += ","; // @TODO: use different separator?
        actions += DataProcessing::NamesOfProcessingAction[action];
      }
      query.bindValue(":processing_actions", actions.toQString());
      query.bindValue(":completion_time", proc.getCompletionTime().get().toQString());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
      index++;
    }
    storeMetaInfos_(features.getDataProcessing(), "FEAT_DataProcessing");
  }


  void OMSFileStore::store(const FeatureMap& features)
  {
    QSqlDatabase db = QSqlDatabase::database(db_name_);
    db.transaction(); // avoid SQLite's "implicit transactions", improve runtime
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
    storeDataProcessing_(features);
    nextProgress();
    storeFeatures_(features);
    db.commit();
    endProgress();
  }
}
