// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/FORMAT/OMSFile.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/SYSTEM/File.h>

// strangely, this is needed for type conversions in "QSqlQuery::bindValue":
#include <QtSql/QSqlQueryModel>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS
{
  int version_number = 1;

  void OMSFile::raiseDBError_(const QSqlError& error, int line,
                              const char* function, const String& context)
  {
    String msg = context + ": " + error.text();
    throw Exception::FailedAPICall(__FILE__, line, function, msg);
  }


  bool OMSFile::tableExists_(const String& db_name, const String& name)
  {
    QSqlDatabase db = QSqlDatabase::database(db_name.toQString());
    return db.tables(QSql::Tables).contains(name.toQString());
  }


  OMSFile::OMSFileStore::OMSFileStore(const String& filename,
                                      const IdentificationData& id_data):
    db_name_("store_" + filename.toQString() + "_" +
             QString::number(UniqueIdGenerator::getUniqueId())),
    id_data_(id_data)
  {
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


  OMSFile::OMSFileStore::~OMSFileStore()
  {
    QSqlDatabase::database(db_name_).close();
    QSqlDatabase::removeDatabase(db_name_);
  }


  void OMSFile::OMSFileStore::createTable_(const String& name,
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


  void OMSFile::OMSFileStore::storeVersionAndDate()
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


  void OMSFile::OMSFileStore::createTableMoleculeType_()
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


  void OMSFile::OMSFileStore::createTableDataValue_()
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
  }


  OMSFile::Key OMSFile::OMSFileStore::storeDataValue_(const DataValue& value)
  {
    // this assumes the "DataValue" table exists already!
    // @TODO: split this up and make several tables for different types?
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO DataValue VALUES (" \
                  "NULL, "                         \
                  ":data_type, "                   \
                  ":value)");
    // @TODO: cache the prepared query between function calls somehow?
    if (!value.isEmpty()) // use NULL as the type for empty values
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


  void OMSFile::OMSFileStore::createTableCVTerm_()
  {
    createTable_("CVTerm",
                 "id INTEGER PRIMARY KEY NOT NULL, "        \
                 "accession TEXT UNIQUE, "                  \
                 "name TEXT NOT NULL, "                     \
                 "cv_identifier_ref TEXT, "                 \
                 // does this constrain "name" if "accession" is NULL?
                 "UNIQUE (accession, name)");
    // @TODO: add support for unit and value
  }


  OMSFile::Key OMSFile::OMSFileStore::storeCVTerm_(const CVTerm& cv_term)
  {
    // this assumes the "CVTerm" table exists already!
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT OR IGNORE INTO CVTerm VALUES ("   \
                  "NULL, "                                  \
                  ":accession, "                            \
                  ":name, "                                 \
                  ":cv_identifier_ref)");
    if (!cv_term.getAccession().empty()) // use NULL for empty accessions
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
    query.prepare("SELECT id FROM CVTerm "                          \
                  "WHERE accession = :accession AND name = :name");
    if (!cv_term.getAccession().empty()) // use NULL for empty accessions
    {
      query.bindValue(":accession", cv_term.getAccession().toQString());
    }
    query.bindValue(":name", cv_term.getName().toQString());
    if (!query.exec() || !query.next())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error querying database");
    }
    return query.value(0).toLongLong();
  }


  void OMSFile::OMSFileStore::createTableMetaInfo_(const String& parent_table)
  {
    if (!tableExists_(db_name_, "DataValue")) createTableDataValue_();

    createTable_(
      parent_table + "_MetaInfo",
      "parent_id INTEGER NOT NULL, "                                    \
      "name TEXT NOT NULL, "                                            \
      "data_value_id INTEGER NOT NULL, "                                \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_table + " (id), "  \
      "FOREIGN KEY (data_value_id) REFERENCES DataValue (id), "         \
      "UNIQUE (parent_id, name)");
  }


  void OMSFile::OMSFileStore::storeMetaInfo_(const MetaInfoInterface& info,
                                             const String& parent_table,
                                             Key parent_id)
  {
    if (info.isMetaEmpty()) return;

    // this assumes the "..._MetaInfo" and "DataValue" tables exist already!
    String table = parent_table + "_MetaInfo";

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO " + table.toQString() + " VALUES ("  \
                  ":parent_id, "                                    \
                  ":name, "                                         \
                  ":data_value_id)");
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


  void OMSFile::OMSFileStore::createTableAppliedProcessingStep_(
    const String& parent_table)
  {
    createTable_(
      parent_table + "_AppliedProcessingStep",
      "parent_id INTEGER NOT NULL, "                                    \
      "processing_step_id INTEGER, "                                    \
      "processing_step_order INTEGER NOT NULL, "                        \
      "score_type_id INTEGER, "                                         \
      "score REAL, "                                                    \
      "UNIQUE (parent_id, processing_step_id, score_type_id), "         \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_table + " (id), "  \
      "FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id), "      \
      "FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)");
    // @TODO: add constraint that "processing_step_id" and "score_type_id"
    // can't both be NULL
    // @TODO: add constraint that "processing_step_order" must match "..._id"?
    // @TODO: normalize table? (splitting into multiple tables is awkward here)
  }


  void OMSFile::OMSFileStore::storeAppliedProcessingStep_(
    const ID::AppliedProcessingStep& step, Size step_order,
    const String& parent_table, Key parent_id)
  {
    // this assumes the "..._AppliedProcessingStep" table exists already!
    String table = parent_table + "_AppliedProcessingStep";

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO " + table.toQString() + " VALUES ("  \
                  ":parent_id, "                                    \
                  ":processing_step_id, "                           \
                  ":processing_step_order, "                        \
                  ":score_type_id, "                                \
                  ":score)");
    query.bindValue(":parent_id", parent_id);
    query.bindValue(":processing_step_order", int(step_order));
    if (step.processing_step_opt)
    {
      query.bindValue(":processing_step_id",
                      Key(&(**step.processing_step_opt)));
      if (step.scores.empty()) // insert processing step information only
      {
        if (!query.exec())
        {
          raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                        "error inserting data");
        }
      }
    } // else: use NULL for missing processing step reference
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


  void OMSFile::OMSFileStore::storeScoreTypes()
  {
    if (id_data_.getScoreTypes().empty()) return;

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
    for (const ID::ScoreType& score_type : id_data_.getScoreTypes())
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


  void OMSFile::OMSFileStore::storeInputFiles()
  {
    if (id_data_.getInputFiles().empty()) return;

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
    for (const ID::InputFile& input : id_data_.getInputFiles())
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


  void OMSFile::OMSFileStore::storeProcessingSoftwares()
  {
    if (id_data_.getProcessingSoftwares().empty()) return;

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
    for (const ID::ProcessingSoftware& software :
           id_data_.getProcessingSoftwares())
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
        "INSERT INTO ID_ProcessingSoftware_AssignedScore VALUES ("  \
        ":software_id, "                                                \
        ":score_type_id, "                                              \
        ":score_type_order)");
      for (const ID::ProcessingSoftware& software :
             id_data_.getProcessingSoftwares())
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


  void OMSFile::OMSFileStore::storeDBSearchParams()
  {
    if (id_data_.getDBSearchParams().empty()) return;

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
    for (const ID::DBSearchParam& param : id_data_.getDBSearchParams())
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


  void OMSFile::OMSFileStore::storeProcessingSteps()
  {
    if (id_data_.getProcessingSteps().empty()) return;

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
    for (ID::ProcessingStepRef step_ref =
           id_data_.getProcessingSteps().begin(); step_ref !=
           id_data_.getProcessingSteps().end(); ++step_ref)
    {
      const ID::ProcessingStep& step = *step_ref;
      if (!step.input_file_refs.empty()) any_input_files = true;
      query.bindValue(":id", Key(&step));
      query.bindValue(":software_id", Key(&(*step.software_ref)));
      query.bindValue(":date_time", step.date_time.get().toQString());
      auto pos = id_data_.getDBSearchSteps().find(step_ref);
      if (pos != id_data_.getDBSearchSteps().end())
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

      for (const ID::ProcessingStep& step :
             id_data_.getProcessingSteps())
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
    storeMetaInfos_(id_data_.getProcessingSteps(),
                    "ID_ProcessingStep");
  }


  void OMSFile::OMSFileStore::storeInputItems()
  {
    if (id_data_.getInputItems().empty()) return;

    createTable_("ID_InputItem",
                 "id INTEGER PRIMARY KEY NOT NULL, "                    \
                 "data_id TEXT NOT NULL, "                              \
                 "input_file_id INTEGER, "                              \
                 "rt REAL, "                                            \
                 "mz REAL, "                                            \
                 "UNIQUE (data_id, input_file_id), "                    \
                 "FOREIGN KEY (input_file_id) REFERENCES ID_InputFile (id)");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_InputItem VALUES (" \
                  ":id, "                             \
                  ":data_id, "                        \
                  ":input_file_id, "                  \
                  ":rt, "                             \
                  ":mz)");
    for (const ID::InputItem& input_item : id_data_.getInputItems())
    {
      query.bindValue(":id", Key(&input_item)); // use address as primary key
      query.bindValue(":data_id", input_item.data_id.toQString());
      if (input_item.input_file_opt)
      {
        query.bindValue(":input_file_id", Key(&(**input_item.input_file_opt)));
      }
      else
      {
        query.bindValue(":input_file_id", QVariant(QVariant::String)); // NULL
      }
      if (input_item.rt == input_item.rt)
      {
        query.bindValue(":rt", input_item.rt);
      }
      else // NaN
      {
        query.bindValue(":rt", QVariant(QVariant::Double)); // NULL
      }
      if (input_item.mz == input_item.mz)
      {
        query.bindValue(":mz", input_item.mz);
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
    storeMetaInfos_(id_data_.getInputItems(), "ID_InputItem");
  }


  void OMSFile::OMSFileStore::storeParentSequences()
  {
    if (id_data_.getParentSequences().empty()) return;

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
    for (const ID::ParentSequence& parent : id_data_.getParentSequences())
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
    storeScoredProcessingResults_(id_data_.getParentSequences(),
                                  "ID_ParentSequence");
  }


  void OMSFile::OMSFileStore::storeParentGroupSets()
  {
    if (id_data_.getParentGroupSets().empty()) return;

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
    for (const ID::ParentGroupSet& grouping :
           id_data_.getParentGroupSets())
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

    storeScoredProcessingResults_(id_data_.getParentGroupSets(),
                                  "ID_ParentGroupSet");
  }


  void OMSFile::OMSFileStore::createTableIdentifiedMolecule_()
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
  }


  void OMSFile::OMSFileStore::storeIdentifiedCompounds()
  {
    if (id_data_.getIdentifiedCompounds().empty()) return;

    if (!tableExists_(db_name_, "ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }

    createTable_(
      "ID_IdentifiedCompound",
      "molecule_id UNIQUE INTEGER NOT NULL, "                           \
      "formula TEXT, "                                                  \
      "name TEXT, "                                                     \
      "smile TEXT, "                                                    \
      "inchi TEXT, "                                                    \
      "FOREIGN KEY (molecule_id) REFERENCES ID_IdentifiedMolecule (id)");

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query_molecule(db);
    query_molecule.prepare("INSERT INTO ID_IdentifiedMolecule VALUES (" \
                           ":id, "                                      \
                           ":molecule_type_id, "                        \
                           ":identifier)");
    query_molecule.bindValue(":molecule_type_id",
                             int(ID::MoleculeType::COMPOUND) + 1);
    QSqlQuery query_compound(db);
    query_compound.prepare("INSERT INTO ID_IdentifiedCompound VALUES (" \
                           ":molecule_id, "                             \
                           ":formula, "                                 \
                           ":name, "                                    \
                           ":smile, "                                   \
                           ":inchi)");
    for (const ID::IdentifiedCompound& compound :
           id_data_.getIdentifiedCompounds())
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
    storeScoredProcessingResults_(id_data_.getIdentifiedCompounds(),
                                  "ID_IdentifiedMolecule");
  }


  void OMSFile::OMSFileStore::storeIdentifiedSequences()
  {
    if (id_data_.getIdentifiedPeptides().empty() &&
        id_data_.getIdentifiedOligos().empty()) return;

    if (!tableExists_(db_name_, "ID_IdentifiedMolecule"))
    {
      createTableIdentifiedMolecule_();
    }

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_IdentifiedMolecule VALUES (" \
                  ":id, "                                      \
                  ":molecule_type_id, "                        \
                  ":identifier)");
    bool any_parent_matches = false;
    // store peptides:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    for (const ID::IdentifiedPeptide& peptide :
           id_data_.getIdentifiedPeptides())
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
    storeScoredProcessingResults_(id_data_.getIdentifiedPeptides(),
                                  "ID_IdentifiedMolecule");
    // store RNA oligos:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    for (const ID::IdentifiedOligo& oligo : id_data_.getIdentifiedOligos())
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
    storeScoredProcessingResults_(id_data_.getIdentifiedOligos(),
                                  "ID_IdentifiedMolecule");

    if (any_parent_matches)
    {
      createTableParentMatches_();
      for (const ID::IdentifiedPeptide& peptide :
             id_data_.getIdentifiedPeptides())
      {
        if (peptide.parent_matches.empty()) continue;
        storeParentMatches_(peptide.parent_matches, Key(&peptide));
      }
      for (const ID::IdentifiedOligo& oligo : id_data_.getIdentifiedOligos())
      {
        if (oligo.parent_matches.empty()) continue;
        storeParentMatches_(oligo.parent_matches, Key(&oligo));
      }
    }
  }


  void OMSFile::OMSFileStore::createTableParentMatches_()
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
  }


  void OMSFile::OMSFileStore::storeParentMatches_(
    const ID::ParentMatches& matches, Key molecule_id)
  {
    // this assumes the "ID_ParentMatch" table exists already!
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ParentMatch VALUES (" \
                  ":molecule_id, "                              \
                  ":parent_id, "                                \
                  ":start_pos, "                                \
                  ":end_pos, "                                  \
                  ":left_neighbor, "                            \
                  ":right_neighbor)");
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


  void OMSFile::OMSFileStore::storeAdducts()
  {
    if (id_data_.getAdducts().empty()) return;

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
    for (const AdductInfo& adduct : id_data_.getAdducts())
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


  void OMSFile::OMSFileStore::storeInputMatches()
  {
    if (id_data_.getInputMatches().empty()) return;

    String table_def =
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "identified_molecule_id INTEGER NOT NULL, "                       \
      "input_item_id INTEGER NOT NULL, "                                \
      "adduct_id INTEGER, "                                             \
      "charge INTEGER, "                                                \
      "FOREIGN KEY (identified_molecule_id) REFERENCES ID_IdentifiedMolecule (id), " \
      "FOREIGN KEY (input_item_id) REFERENCES ID_InputItem (id)";
    // add foreign key constraint if the adduct table exists (having the
    // constraint without the table would cause an error on data insertion):
    if (tableExists_(db_name_, "AdductInfo"))
    {
      table_def += ", FOREIGN KEY (adduct_id) REFERENCES AdductInfo (id)";
    }
    createTable_("ID_InputMatch", table_def);

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_InputMatch VALUES ("  \
                  ":id, "                                       \
                  ":identified_molecule_id, "                   \
                  ":input_item_id, "                            \
                  ":adduct_id, "                                \
                  ":charge)");
    bool any_peak_annotations = false;
    for (const ID::InputMatch& match :
           id_data_.getInputMatches())
    {
      if (!match.peak_annotations.empty()) any_peak_annotations = true;
      query.bindValue(":id", Key(&match)); // use address as primary key
      Key molecule_id;
      const ID::IdentifiedMolecule& molecule_var =
        match.identified_molecule_var;
      switch (molecule_var.getMoleculeType())
      {
      case ID::MoleculeType::PROTEIN:
        molecule_id = Key(&(*molecule_var.getIdentifiedPeptideRef()));
        break;
      case ID::MoleculeType::COMPOUND:
        molecule_id = Key(&(*molecule_var.getIdentifiedCompoundRef()));
        break;
      case ID::MoleculeType::RNA:
        molecule_id = Key(&(*molecule_var.getIdentifiedOligoRef()));
        break;
      default:
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION,
                                         "invalid molecule type");
      }
      query.bindValue(":identified_molecule_id", molecule_id);
      query.bindValue(":input_item_id", Key(&(*match.input_item_ref)));
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
    storeScoredProcessingResults_(id_data_.getInputMatches(),
                                  "ID_InputMatch");

    if (any_peak_annotations)
    {
      createTable_(
        "ID_InputMatch_PeakAnnotation",
        "parent_id INTEGER NOT NULL, "                                  \
        "processing_step_id INTEGER, "                                  \
        "peak_annotation TEXT, "                                        \
        "peak_charge INTEGER, "                                         \
        "peak_mz REAL, "                                                \
        "peak_intensity REAL, "                                         \
        "FOREIGN KEY (parent_id) REFERENCES ID_InputMatch (id), " \
        "FOREIGN KEY (processing_step_id) REFERENCES ID_ProcessingStep (id)");

      query.prepare(
        "INSERT INTO ID_InputMatch_PeakAnnotation VALUES (" \
        ":parent_id, "                                              \
        ":processing_step_id, "                                     \
        ":peak_annotation, "                                        \
        ":peak_charge, "                                            \
        ":peak_mz, "                                                \
        ":peak_intensity)");

      for (const ID::InputMatch& match :
             id_data_.getInputMatches())
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
  }


  void OMSFile::store(const String& filename, const IdentificationData& id_data)
  {
    OMSFileStore helper(filename, id_data);
    startProgress(0, 13, "Writing data to file");
    // generally, create tables only if we have data to write - no empty ones!
    helper.storeVersionAndDate();
    nextProgress();
    helper.storeInputFiles();
    nextProgress();
    helper.storeScoreTypes();
    nextProgress();
    helper.storeProcessingSoftwares();
    nextProgress();
    helper.storeDBSearchParams();
    nextProgress();
    helper.storeProcessingSteps();
    nextProgress();
    helper.storeInputItems();
    nextProgress();
    helper.storeParentSequences();
    nextProgress();
    helper.storeParentGroupSets();
    nextProgress();
    helper.storeIdentifiedCompounds();
    nextProgress();
    helper.storeIdentifiedSequences();
    nextProgress();
    helper.storeAdducts();
    nextProgress();
    helper.storeInputMatches();
    endProgress();
    // @TODO: store input match groups
  }


  OMSFile::OMSFileLoad::OMSFileLoad(const String& filename,
                                    IdentificationData& id_data):
    db_name_("load_" + filename.toQString() + "_" +
             QString::number(UniqueIdGenerator::getUniqueId())),
    id_data_(id_data)
  {
    // open database:
    QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", db_name_);
    db.setDatabaseName(filename.toQString());
    db.setConnectOptions("QSQLITE_OPEN_READONLY");
    if (!db.open())
    {
      raiseDBError_(db.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error opening SQLite database");
      // if d'tor doesn't get called, DB connection (db_name_) doesn't get
      // removed, but that shouldn't be a big problem
    }
  }


  OMSFile::OMSFileLoad::~OMSFileLoad()
  {
    QSqlDatabase::database(db_name_).close();
    QSqlDatabase::removeDatabase(db_name_);
  }


  // currently not needed:
  // CVTerm OMSFile::OMSFileLoad::loadCVTerm_(int id)
  // {
  //   // this assumes that the "CVTerm" table exists!
  //   QSqlQuery query(db_);
  //   query.setForwardOnly(true);
  //   QString sql_select = "SELECT * FROM CVTerm WHERE id = " + QString(id);
  //   if (!query.exec(sql_select) || !query.next())
  //   {
  //     raiseDBError_(model.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
  //                   "error reading from database");
  //   }
  //   return CVTerm(query.value("accession").toString(),
  //                 query.value("name").toString(),
  //                 query.value("cv_identifier_ref").toString());
  // }


  void OMSFile::OMSFileLoad::loadScoreTypes()
  {
    if (!tableExists_(db_name_, "ID_ScoreType")) return;
    if (!tableExists_(db_name_, "CVTerm")) // every score type is a CV term
    {
      String msg = "required database table 'CVTerm' not found";
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    // careful - both joined tables have an "id" field, need to exclude one:
    if (!query.exec("SELECT S.*, C.accession, C.name, C.cv_identifier_ref " \
                    "FROM ID_ScoreType AS S JOIN CVTerm AS C "          \
                    "ON S.cv_term_id = C.id"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      CVTerm cv_term(query.value("accession").toString(),
                     query.value("name").toString(),
                     query.value("cv_identifier_ref").toString());
      bool higher_better = query.value("higher_better").toInt();
      ID::ScoreType score_type(cv_term, higher_better);
      ID::ScoreTypeRef ref = id_data_.registerScoreType(score_type);
      score_type_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadInputFiles()
  {
    if (!tableExists_(db_name_, "ID_InputFile")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_InputFile"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::InputFile input(query.value("name").toString(),
                          query.value("experimental_design_id").toString());
      String primary_files = query.value("primary_files").toString();
      vector<String> pf_list = ListUtils::create<String>(primary_files);
      input.primary_files.insert(pf_list.begin(), pf_list.end());
      ID::InputFileRef ref = id_data_.registerInputFile(input);
      input_file_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadProcessingSoftwares()
  {
    if (!tableExists_(db_name_, "ID_ProcessingSoftware")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ProcessingSoftware"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    bool have_scores = tableExists_(db_name_,
                                    "ID_ProcessingSoftware_AssignedScore");
    QSqlQuery subquery(db);
    if (have_scores)
    {
      subquery.setForwardOnly(true);
      subquery.prepare("SELECT score_type_id "                         \
                       "FROM ID_ProcessingSoftware_AssignedScore " \
                       "WHERE software_id = :id ORDER BY score_type_order ASC");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      ID::ProcessingSoftware software(query.value("name").toString(),
                                          query.value("version").toString());
      if (have_scores)
      {
        subquery.bindValue(":id", id);
        if (!subquery.exec())
        {
          raiseDBError_(subquery.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                        "error reading from database");
        }
        while (subquery.next())
        {
          Key score_type_id = subquery.value(0).toLongLong();
          software.assigned_scores.push_back(score_type_refs_[score_type_id]);
        }
      }
      ID::ProcessingSoftwareRef ref =
        id_data_.registerProcessingSoftware(software);
      processing_software_refs_[id] = ref;
    }
  }


  DataValue OMSFile::OMSFileLoad::makeDataValue_(const QSqlQuery& query)
  {
    DataValue::DataType type = DataValue::EMPTY_VALUE;
    int type_index = query.value("data_type_id").toInt();
    if (type_index > 0) type = DataValue::DataType(type_index - 1);
    String value = query.value("value").toString();
    switch (type)
    {
    case DataValue::STRING_VALUE:
      return DataValue(value);
    case DataValue::INT_VALUE:
      return DataValue(value.toInt());
    case DataValue::DOUBLE_VALUE:
      return DataValue(value.toDouble());
    case DataValue::STRING_LIST:
      return DataValue(ListUtils::create<String>(value));
    case DataValue::INT_LIST:
      return DataValue(ListUtils::create<int>(value));
    case DataValue::DOUBLE_LIST:
      return DataValue(ListUtils::create<double>(value));
    default: // DataValue::EMPTY_VALUE (avoid warning about missing return)
      return DataValue();
    }
  }


  bool OMSFile::OMSFileLoad::prepareQueryMetaInfo_(QSqlQuery& query,
                                                   const String& parent_table)
  {
    String table_name = parent_table + "_MetaInfo";
    if (!tableExists_(db_name_, table_name)) return false;

    query.setForwardOnly(true);
    QString sql_select =
      "SELECT * FROM " + table_name.toQString() + " AS MI " \
      "JOIN DataValue AS DV ON MI.data_value_id = DV.id "   \
      "WHERE MI.parent_id = :id";
    query.prepare(sql_select);
    return true;
  }


  bool OMSFile::OMSFileLoad::prepareQueryAppliedProcessingStep_(
    QSqlQuery& query, const String& parent_table)
  {
    String table_name = parent_table + "_AppliedProcessingStep";
    if (!tableExists_(db_name_, table_name)) return false;

    // query.setForwardOnly(true);
    QString sql_select = "SELECT * FROM " + table_name.toQString() +
      " WHERE parent_id = :id ORDER BY processing_step_order ASC";
    query.prepare(sql_select);
    return true;
  }


  void OMSFile::OMSFileLoad::handleQueryMetaInfo_(QSqlQuery& query,
                                                  MetaInfoInterface& info,
                                                  Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      DataValue value = makeDataValue_(query);
      info.setMetaValue(query.value("name").toString(), value);
    }
  }


  void OMSFile::OMSFileLoad::handleQueryAppliedProcessingStep_(
    QSqlQuery& query,
    IdentificationDataInternal::ScoredProcessingResult& result,
    Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::AppliedProcessingStep step;
      QVariant step_id_opt = query.value("processing_step_id");
      if (!step_id_opt.isNull())
      {
        step.processing_step_opt =
          processing_step_refs_[step_id_opt.toLongLong()];
      }
      QVariant score_type_opt = query.value("score_type_id");
      if (!score_type_opt.isNull())
      {
        step.scores[score_type_refs_[score_type_opt.toLongLong()]] =
          query.value("score").toDouble();
      }
      result.addProcessingStep(step); // this takes care of merging the steps
    }
  }


  void OMSFile::OMSFileLoad::loadDBSearchParams()
  {
    if (!tableExists_(db_name_, "ID_DBSearchParam")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_DBSearchParam"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      ID::DBSearchParam param;
      int molecule_type_index = query.value("molecule_type_id").toInt() - 1;
      param.molecule_type = ID::MoleculeType(molecule_type_index);
      int mass_type_index = query.value("mass_type_average").toInt();
      param.mass_type = ID::MassType(mass_type_index);
      param.database = query.value("database").toString();
      param.database_version = query.value("database_version").toString();
      param.taxonomy = query.value("taxonomy").toString();
      vector<Int> charges =
        ListUtils::create<Int>(query.value("charges").toString());
      param.charges.insert(charges.begin(), charges.end());
      vector<String> fixed_mods =
        ListUtils::create<String>(query.value("fixed_mods").toString());
      param.fixed_mods.insert(fixed_mods.begin(), fixed_mods.end());
      vector<String> variable_mods =
        ListUtils::create<String>(query.value("variable_mods").toString());
      param.variable_mods.insert(variable_mods.begin(), variable_mods.end());
      param.precursor_mass_tolerance =
        query.value("precursor_mass_tolerance").toDouble();
      param.fragment_mass_tolerance =
        query.value("fragment_mass_tolerance").toDouble();
      param.precursor_tolerance_ppm =
        query.value("precursor_tolerance_ppm").toInt();
      param.fragment_tolerance_ppm =
        query.value("fragment_tolerance_ppm").toInt();
      String enzyme = query.value("digestion_enzyme").toString();
      if (!enzyme.empty())
      {
        if (param.molecule_type == ID::MoleculeType::PROTEIN)
        {
          param.digestion_enzyme = ProteaseDB::getInstance()->getEnzyme(enzyme);
        }
        else if (param.molecule_type == ID::MoleculeType::RNA)
        {
          param.digestion_enzyme = RNaseDB::getInstance()->getEnzyme(enzyme);
        }
      }
      param.missed_cleavages = query.value("missed_cleavages").toUInt();
      param.min_length = query.value("min_length").toUInt();
      param.max_length = query.value("max_length").toUInt();
      ID::SearchParamRef ref = id_data_.registerDBSearchParam(param);
      search_param_refs_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadProcessingSteps()
  {
    if (!tableExists_(db_name_, "ID_ProcessingStep")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ProcessingStep"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_file(db);
    bool have_input_files = tableExists_(db_name_,
                                         "ID_ProcessingStep_InputFile");
    if (have_input_files)
    {
      subquery_file.setForwardOnly(true);
      subquery_file.prepare("SELECT input_file_id "                 \
                            "FROM ID_ProcessingStep_InputFile " \
                            "WHERE processing_step_id = :id");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                              "ID_ProcessingStep");
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      Key software_id = query.value("software_id").toLongLong();
      ID::ProcessingStep step(processing_software_refs_[software_id]);
      String date_time = query.value("date_time").toString();
      if (!date_time.empty()) step.date_time.set(date_time);
      if (have_input_files)
      {
        subquery_file.bindValue(":id", id);
        if (!subquery_file.exec())
        {
          raiseDBError_(subquery_file.lastError(), __LINE__,
                        OPENMS_PRETTY_FUNCTION, "error reading from database");
        }
        while (subquery_file.next())
        {
          Key input_file_id = subquery_file.value(0).toLongLong();
          // the foreign key constraint should ensure that look-up succeeds:
          step.input_file_refs.push_back(input_file_refs_[input_file_id]);
        }
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, step, id);
      }
      ID::ProcessingStepRef ref;
      QVariant opt_search_param_id = query.value("search_param_id");
      if (opt_search_param_id.isNull()) // no DB search params available
      {
        ref = id_data_.registerProcessingStep(step);
      }
      else
      {
        ID::SearchParamRef search_param_ref =
          search_param_refs_[opt_search_param_id.toLongLong()];
        ref = id_data_.registerProcessingStep(step, search_param_ref);
      }
      processing_step_refs_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadInputItems()
  {
    if (!tableExists_(db_name_, "ID_InputItem")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_InputItem"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_InputItem");

    while (query.next())
    {
      ID::InputItem input_item(query.value("data_id").toString());
      QVariant input_file_id = query.value("input_file_id");
      if (!input_file_id.isNull())
      {
        input_item.input_file_opt =
          input_file_refs_[input_file_id.toLongLong()];
      }
      QVariant rt = query.value("rt");
      if (!rt.isNull()) input_item.rt = rt.toDouble();
      QVariant mz = query.value("mz");
      if (!mz.isNull()) input_item.mz = mz.toDouble();
      Key id = query.value("id").toLongLong();
      if (have_meta_info) handleQueryMetaInfo_(subquery_info, input_item, id);
      ID::InputItemRef ref = id_data_.registerInputItem(input_item);
      input_item_refs_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadParentSequences()
  {
    if (!tableExists_(db_name_, "ID_ParentSequence")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentSequence"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentSequence");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step, "ID_ParentSequence");

    while (query.next())
    {
      String accession = query.value("accession").toString();
      ID::ParentSequence parent(accession);
      int molecule_type_index = query.value("molecule_type_id").toInt() - 1;
      parent.molecule_type = ID::MoleculeType(molecule_type_index);
      parent.sequence = query.value("sequence").toString();
      parent.description = query.value("description").toString();
      parent.coverage = query.value("coverage").toDouble();
      parent.is_decoy = query.value("is_decoy").toInt();
      Key id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, parent, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, parent, id);
      }
      ID::ParentSequenceRef ref = id_data_.registerParentSequence(parent);
      parent_refs_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadParentGroupSets()
  {
    if (!tableExists_(db_name_, "ID_ParentGroupSet")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentGroupSet "  \
                    "ORDER BY grouping_order ASC"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentGroupSet");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_ParentGroupSet");

    QSqlQuery subquery_group(db);
    subquery_group.setForwardOnly(true);
    subquery_group.prepare("SELECT * FROM ID_ParentGroup "  \
                           "WHERE grouping_id = :id");

    QSqlQuery subquery_parent(db);
    subquery_parent.setForwardOnly(true);
    subquery_parent.prepare(
      "SELECT parent_id FROM ID_ParentGroup_ParentSequence " \
      "WHERE group_id = :id");

    while (query.next())
    {
      ID::ParentGroupSet grouping(query.value("label").toString());
      Key grouping_id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, grouping, grouping_id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, grouping, grouping_id);
      }

      subquery_group.bindValue(":id", grouping_id);
      if (!subquery_group.exec())
      {
        raiseDBError_(subquery_group.lastError(), __LINE__,
                      OPENMS_PRETTY_FUNCTION, "error reading from database");
      }

      // get all groups in this grouping:
      map<Key, ID::ParentGroup> groups_map;
      while (subquery_group.next())
      {
        Key group_id = subquery_group.value("id").toLongLong();
        QVariant score_type_id = subquery_group.value("score_type_id");
        if (score_type_id.isNull()) // no scores
        {
          groups_map[group_id]; // insert empty group
        }
        else
        {
          ID::ScoreTypeRef ref = score_type_refs_[score_type_id.toLongLong()];
          groups_map[group_id].scores[ref] =
            subquery_group.value("score").toDouble();
        }
      }
      // get parent sequences in each group:
      for (auto& pair : groups_map)
      {
        subquery_parent.bindValue(":id", pair.first);
        if (!subquery_parent.exec())
        {
          raiseDBError_(subquery_parent.lastError(), __LINE__,
                        OPENMS_PRETTY_FUNCTION, "error reading from database");
        }
        while (subquery_parent.next())
        {
          Key parent_id = subquery_parent.value(0).toLongLong();
          pair.second.parent_refs.insert(
            parent_refs_[parent_id]);
        }
        grouping.groups.insert(pair.second);
      }

      id_data_.registerParentGroupSet(grouping);
    }
  }


  void OMSFile::OMSFileLoad::loadIdentifiedCompounds()
  {
    if (!tableExists_(db_name_, "ID_IdentifiedCompound")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    QString sql_select =
      "SELECT * FROM ID_IdentifiedMolecule JOIN ID_IdentifiedCompound " \
      "ON ID_IdentifiedMolecule.id = ID_IdentifiedCompound.molecule_id";
    if (!query.exec(sql_select))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_IdentifiedMolecule");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_IdentifiedMolecule");

    while (query.next())
    {
      ID::IdentifiedCompound compound(
        query.value("identifier").toString(),
        EmpiricalFormula(query.value("formula").toString()),
        query.value("name").toString(),
        query.value("smile").toString(),
        query.value("inchi").toString());
      Key id = query.value("id").toLongLong();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, compound, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, compound, id);
      }
      ID::IdentifiedCompoundRef ref =
        id_data_.registerIdentifiedCompound(compound);
      identified_molecule_vars_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::handleQueryParentMatch_(
    QSqlQuery& query, IdentificationData::ParentMatches& parent_matches,
    Key molecule_id)
  {
    query.bindValue(":id", molecule_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      ID::ParentSequenceRef ref =
        parent_refs_[query.value("parent_id").toLongLong()];
      ID::ParentMatch match;
      QVariant start_pos = query.value("start_pos");
      QVariant end_pos = query.value("end_pos");
      if (!start_pos.isNull()) match.start_pos = start_pos.toInt();
      if (!end_pos.isNull()) match.end_pos = end_pos.toInt();
      match.left_neighbor = query.value("left_neighbor").toString();
      match.right_neighbor = query.value("right_neighbor").toString();
      parent_matches[ref].insert(match);
    }
  }


  void OMSFile::OMSFileLoad::loadIdentifiedSequences()
  {
    if (!tableExists_(db_name_, "ID_IdentifiedMolecule")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    query.prepare("SELECT * FROM ID_IdentifiedMolecule "          \
                  "WHERE molecule_type_id = :molecule_type_id");
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_IdentifiedMolecule");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_IdentifiedMolecule");
    QSqlQuery subquery_parent(db);
    bool have_parent_matches = tableExists_(db_name_,
                                            "ID_ParentMatch");
    if (have_parent_matches)
    {
      subquery_parent.setForwardOnly(true);
      subquery_parent.prepare("SELECT * FROM ID_ParentMatch " \
                              "WHERE molecule_id = :id");
    }

    // load peptides:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::PROTEIN) + 1);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      String sequence = query.value("identifier").toString();
      ID::IdentifiedPeptide peptide(AASequence::fromString(sequence));
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, peptide, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, peptide, id);
      }
      if (have_parent_matches)
      {
        handleQueryParentMatch_(subquery_parent, peptide.parent_matches, id);
      }
      ID::IdentifiedPeptideRef ref =
        id_data_.registerIdentifiedPeptide(peptide);
      identified_molecule_vars_[id] = ref;
    }

    // load RNA oligos:
    query.bindValue(":molecule_type_id", int(ID::MoleculeType::RNA) + 1);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      String sequence = query.value("identifier").toString();
      ID::IdentifiedOligo oligo(NASequence::fromString(sequence));
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, oligo, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, oligo, id);
      }
      if (have_parent_matches)
      {
        handleQueryParentMatch_(subquery_parent, oligo.parent_matches, id);
      }
      ID::IdentifiedOligoRef ref = id_data_.registerIdentifiedOligo(oligo);
      identified_molecule_vars_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::handleQueryPeakAnnotation_(
    QSqlQuery& query, ID::InputMatch& match, Key parent_id)
  {
    query.bindValue(":id", parent_id);
    if (!query.exec())
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      QVariant processing_step_id = query.value("processing_step_id");
      boost::optional<ID::ProcessingStepRef> processing_step_opt = boost::none;
      if (!processing_step_id.isNull())
      {
        processing_step_opt =
          processing_step_refs_[processing_step_id.toLongLong()];
      }
      PeptideHit::PeakAnnotation ann;
      ann.annotation = query.value("peak_annotation").toString();
      ann.charge = query.value("peak_charge").toInt();
      ann.mz = query.value("peak_mz").toDouble();
      ann.intensity = query.value("peak_intensity").toDouble();
      match.peak_annotations[processing_step_opt].push_back(ann);
    }
  }


  void OMSFile::OMSFileLoad::loadAdducts()
  {
    if (!tableExists_(db_name_, "AdductInfo")) return;

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM AdductInfo"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }

    while (query.next())
    {
      EmpiricalFormula formula(query.value("formula").toString());
      AdductInfo adduct(query.value("name").toString(), formula,
                        query.value("charge").toInt(),
                        query.value("mol_multiplier").toInt());
      ID::AdductRef ref = id_data_.registerAdduct(adduct);
      adduct_refs_[query.value("id").toLongLong()] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadInputMatches()
  {
    if (!tableExists_(db_name_, "ID_InputMatch")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_InputMatch"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    // @TODO: can we combine handling of meta info and applied processing steps?
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_InputMatch");
    QSqlQuery subquery_step(db);
    bool have_applied_steps =
      prepareQueryAppliedProcessingStep_(subquery_step,
                                         "ID_InputMatch");
    QSqlQuery subquery_ann(db);
    bool have_peak_annotations =
      tableExists_(db_name_, "ID_InputMatch_PeakAnnotation");
    if (have_peak_annotations)
    {
      subquery_ann.setForwardOnly(true);
      subquery_ann.prepare(
        "SELECT * FROM ID_InputMatch_PeakAnnotation " \
        "WHERE parent_id = :id");
    }

    while (query.next())
    {
      Key id = query.value("id").toLongLong();
      Key molecule_id = query.value("identified_molecule_id").toLongLong();
      Key query_id = query.value("input_item_id").toLongLong();
      ID::InputMatch match(identified_molecule_vars_[molecule_id],
                                   input_item_refs_[query_id],
                                   query.value("charge").toInt());
      QVariant adduct_id = query.value("adduct_id"); // adduct is optional
      if (!adduct_id.isNull())
      {
        match.adduct_opt = adduct_refs_[adduct_id.toLongLong()];
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, match, id);
      }
      if (have_applied_steps)
      {
        handleQueryAppliedProcessingStep_(subquery_step, match, id);
      }
      if (have_peak_annotations)
      {
        handleQueryPeakAnnotation_(subquery_ann, match, id);
      }
      ID::InputMatchRef ref = id_data_.registerInputMatch(match);
      input_match_refs_[id] = ref;
    }
  }


  void OMSFile::load(const String& filename, IdentificationData& id_data)
  {
    OMSFileLoad helper(filename, id_data);
    startProgress(0, 12, "Reading data from file");
    helper.loadInputFiles();
    nextProgress();
    helper.loadScoreTypes();
    nextProgress();
    helper.loadProcessingSoftwares();
    nextProgress();
    helper.loadDBSearchParams();
    nextProgress();
    helper.loadProcessingSteps();
    nextProgress();
    helper.loadInputItems();
    nextProgress();
    helper.loadParentSequences();
    nextProgress();
    helper.loadParentGroupSets();
    nextProgress();
    helper.loadIdentifiedCompounds();
    nextProgress();
    helper.loadIdentifiedSequences();
    nextProgress();
    helper.loadAdducts();
    nextProgress();
    helper.loadInputMatches();
    endProgress();
    // @TODO: load parent sequence groups and input match groups
  }

}
