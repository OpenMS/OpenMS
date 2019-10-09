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
// $Authors: Julianus Pfeuffer, Oliver Alka, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFile.h>

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
    return query.lastInsertId().toInt();
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
    if (query.lastInsertId().isValid()) return query.lastInsertId().toInt();
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
    return query.value(0).toInt();
  }


  void OMSFile::OMSFileStore::createTableMetaInfo_(const String& parent_table)
  {
    if (!tableExists_(db_name_, "DataValue")) createTableDataValue_();

    createTable_(
      parent_table + "_MetaInfo",
      "parent_id INTEGER PRIMARY KEY NOT NULL, "                        \
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
      "score_type_id INTEGER, "                                         \
      "score REAL, "                                                    \
      "UNIQUE (parent_id, processing_step_id, score_type_id), "         \
      "FOREIGN KEY (parent_id) REFERENCES " + parent_table + " (id), "  \
      "FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id), "      \
      "FOREIGN KEY (processing_step_id) REFERENCES ID_DataProcessingStep (id)");
    // @TODO: add constraint that "processing_step_id" and "score_type_id"
    // can't both be NULL
    // @TODO: normalize table? (splitting into multiple tables is awkward here)
  }


  void OMSFile::OMSFileStore::storeAppliedProcessingStep_(
    const ID::AppliedProcessingStep& step, const String& parent_table,
    Key parent_id)
  {
    // this assumes the "..._AppliedProcessingStep" table exists already!
    String table = parent_table + "_AppliedProcessingStep";

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO " + table.toQString() + " VALUES ("  \
                  ":parent_id, "                                    \
                  ":processing_step_id, "                           \
                  ":score_type_id, "
                  ":score)");
    query.bindValue(":parent_id", parent_id);
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
    }
    for (auto score_pair : step.scores)
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
                 "file TEXT UNIQUE NOT NULL");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_InputFile VALUES ("  \
                  ":id, "                              \
                  ":file)");
    for (const String& input_file : id_data_.getInputFiles())
    {
      query.bindValue(":id", Key(&input_file));
      query.bindValue(":file", input_file.toQString());
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error inserting data");
      }
    }
  }


  void OMSFile::OMSFileStore::storeDataProcessingSoftwares()
  {
    if (id_data_.getDataProcessingSoftwares().empty()) return;

    createTable_("ID_DataProcessingSoftware",
                 "id INTEGER PRIMARY KEY NOT NULL, "  \
                 "name TEXT NOT NULL, "               \
                 "version TEXT, "                     \
                 "UNIQUE (name, version)");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_DataProcessingSoftware VALUES ("  \
                  ":id, "                                           \
                  ":name, "                                         \
                  ":version)");
    bool any_scores = false; // does any software have assigned scores stored?
    for (const ID::DataProcessingSoftware& software :
           id_data_.getDataProcessingSoftwares())
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
        "ID_DataProcessingSoftware_AssignedScore",
        "software_id INTEGER NOT NULL, "                                \
        "score_type_id INTEGER NOT NULL, "                              \
        "UNIQUE (software_id, score_type_id), "                         \
        "FOREIGN KEY (software_id) REFERENCES ID_DataProcessingSoftware (id), " \
        "FOREIGN KEY (score_type_id) REFERENCES ID_ScoreType (id)");

      query.prepare(
        "INSERT INTO ID_DataProcessingSoftware_AssignedScore VALUES ("  \
        ":software_id, "                                                \
        ":score_type_id)");
      for (const ID::DataProcessingSoftware& software :
             id_data_.getDataProcessingSoftwares())
      {
        query.bindValue(":software_id", Key(&software));
        for (ID::ScoreTypeRef score_type_ref : software.assigned_scores)
        {
          query.bindValue(":score_type_id", Key(&(*score_type_ref)));
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


  void OMSFile::OMSFileStore::storeDataProcessingSteps()
  {
    if (id_data_.getDataProcessingSteps().empty()) return;

    createTable_(
      "ID_DataProcessingStep",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "software_id INTEGER NOT NULL, "                                  \
      "primary_files TEXT, "                                            \
      "date_time TEXT, "                                                \
      "search_param_id INTEGER, "                                       \
      "FOREIGN KEY (search_param_id) REFERENCES ID_DBSearchParam (id)");
    // @TODO: add support for processing actions
    // @TODO: store primary files in a separate table (like input files)?
    // @TODO: store (optional) search param reference in a separate table?

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_DataProcessingStep VALUES ("  \
                  ":id, "                                       \
                  ":software_id, "                              \
                  ":primary_files, "                            \
                  ":data_time, "                                \
                  ":search_param_id)");
    bool any_input_files = false;
    // use iterator here because we need one to look up the DB search params:
    for (ID::ProcessingStepRef step_ref =
           id_data_.getDataProcessingSteps().begin(); step_ref !=
           id_data_.getDataProcessingSteps().end(); ++step_ref)
    {
      const ID::DataProcessingStep& step = *step_ref;
      if (!step.input_file_refs.empty()) any_input_files = true;
      query.bindValue(":id", Key(&step));
      query.bindValue(":software_id", Key(&(*step.software_ref)));
      // @TODO: what if a primary file name contains ","?
      String primary_files = ListUtils::concatenate(step.primary_files, ",");
      query.bindValue(":primary_files", primary_files.toQString());
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
        "ID_DataProcessingStep_InputFile",
        "processing_step_id INTEGER NOT NULL, "                         \
        "input_file_id INTEGER NOT NULL, "                              \
        "FOREIGN KEY processing_step_id REFERENCES ID_DataProcessingStep (id), " \
        "FOREIGN KEY input_file_id REFERENCES ID_InputFile (id), "      \
        "UNIQUE (processing_step_id, input_file_id)");

      query.prepare("INSERT INTO ID_DataProcessingStep_InputFile VALUES (" \
                    ":processing_step_id, "                             \
                    ":input_file_id)");

      for (const ID::DataProcessingStep& step :
             id_data_.getDataProcessingSteps())
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
    storeMetaInfos_(id_data_.getDataProcessingSteps(),
                    "ID_DataProcessingStep");
  }


  void OMSFile::OMSFileStore::storeParentMolecules()
  {
    if (id_data_.getParentMolecules().empty()) return;

    if (!tableExists_(db_name_, "ID_MoleculeType")) createTableMoleculeType_();

    createTable_(
      "ID_ParentMolecule",
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "accession TEXT UNIQUE NOT NULL, "                                \
      "molecule_type_id INTEGER NOT NULL, "                             \
      "sequence TEXT, "                                                 \
      "description TEXT, "                                              \
      "coverage REAL, "                                                 \
      "is_decoy NUMERIC NOT NULL CHECK (is_decoy in (0, 1)) DEFAULT 0, " \
      "FOREIGN KEY (molecule_type_id) REFERENCES ID_MoleculeType (id)");

    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.prepare("INSERT INTO ID_ParentMolecule VALUES ("  \
                  ":id, "                                   \
                  ":accession, "                            \
                  ":molecule_type_id, "                     \
                  ":sequence, "                             \
                  ":description, "                          \
                  ":coverage, "                             \
                  ":is_decoy)");
    for (const ID::ParentMolecule& parent : id_data_.getParentMolecules())
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
    storeScoredProcessingResults_(id_data_.getParentMolecules(),
                                  "ID_ParentMolecule");
  }


  void OMSFile::store(const String& filename, const IdentificationData& id_data)
  {
    OMSFileStore helper(filename, id_data);
    // generally, create tables only if we have data to write - no empty ones!
    helper.storeVersionAndDate();
    helper.storeInputFiles();
    helper.storeScoreTypes();
    helper.storeDataProcessingSoftwares();
    helper.storeDBSearchParams();
    helper.storeDataProcessingSteps();
    helper.storeParentMolecules();
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
    score_type_refs_.clear();
    if (!tableExists_(db_name_, "ID_ScoreType")) return;
    if (!tableExists_(db_name_, "CVTerm")) // every score type is a CV term
    {
      String msg = "required database table 'CVTerm' not found";
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    QSqlQuery query(QSqlDatabase::database(db_name_));
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ScoreType JOIN CVTerm "   \
                    "ON ID_ScoreType.cv_term_id = CVTerm.id"))
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
      score_type_refs_[query.value("id").toInt()] = ref;
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
      ID::InputFileRef ref =
        id_data_.registerInputFile(query.value("file").toString());
      input_file_refs_[query.value("id").toInt()] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadDataProcessingSoftwares()
  {
    if (!tableExists_(db_name_, "ID_DataProcessingSoftware")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_DataProcessingSoftware"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    bool have_scores = tableExists_(db_name_,
                                    "ID_DataProcessingSoftware_AssignedScore");
    QSqlQuery subquery(db);
    if (have_scores)
    {
      subquery.setForwardOnly(true);
      subquery.prepare("SELECT score_type_id "                         \
                       "FROM ID_DataProcessingSoftware_AssignedScore " \
                       "WHERE software_id = :id");
    }
    while (query.next())
    {
      Key id = query.value("id").toInt();
      ID::DataProcessingSoftware software(query.value("name").toString(),
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
          Key score_type_id = query.value(0).toInt();
          // the foreign key constraint should ensure that look-up succeeds:
          software.assigned_scores.push_back(score_type_refs_[score_type_id]);
        }
        // order in the vector should be the same as in the table:
        reverse(software.assigned_scores.begin(),
                software.assigned_scores.end());
      }
      ID::ProcessingSoftwareRef ref =
        id_data_.registerDataProcessingSoftware(software);
      processing_software_refs_[query.value("id").toInt()] = ref;
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
      Key id = query.value("id").toInt();
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
      param.variable_mods.insert(fixed_mods.begin(), fixed_mods.end());
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


  void OMSFile::OMSFileLoad::loadDataProcessingSteps()
  {
    if (!tableExists_(db_name_, "ID_DataProcessingStep")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_DataProcessingStep"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_file(db);
    bool have_input_files = tableExists_(db_name_,
                                         "ID_DataProcessingStep_InputFile");
    if (have_input_files)
    {
      subquery_file.setForwardOnly(true);
      subquery_file.prepare("SELECT input_file_id "                 \
                            "FROM ID_DataProcessingStep_InputFile " \
                            "WHERE processing_step_id = :id");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                              "ID_DataProcessingStep");
    while (query.next())
    {
      Key id = query.value("id").toInt();
      Key software_id = query.value("software_id").toInt();
      ID::DataProcessingStep step(processing_software_refs_[software_id]);
      String primary_files = query.value("primary_files").toString();
      step.primary_files = ListUtils::create<String>(primary_files);
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
          Key input_file_id = subquery_file.value(0).toInt();
          // the foreign key constraint should ensure that look-up succeeds:
          step.input_file_refs.push_back(input_file_refs_[input_file_id]);
        }
        // order in the vector should be the same as in the table:
        reverse(step.input_file_refs.begin(), step.input_file_refs.end());
      }
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, step, id);
      }
      ID::ProcessingStepRef ref;
      QVariant opt_search_param_id = query.value("search_param_id");
      if (opt_search_param_id.isNull()) // no DB search params available
      {
        ref = id_data_.registerDataProcessingStep(step);
      }
      else
      {
        ID::SearchParamRef search_param_ref =
          search_param_refs_[opt_search_param_id.toInt()];
        ref = id_data_.registerDataProcessingStep(step, search_param_ref);
      }
      processing_step_refs_[id] = ref;
    }
  }


  void OMSFile::OMSFileLoad::loadParentMolecules()
  {
    if (!tableExists_(db_name_, "ID_ParentMolecule")) return;

    QSqlDatabase db = QSqlDatabase::database(db_name_);
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentMolecule"))
    {
      raiseDBError_(query.lastError(), __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    QSqlQuery subquery_info(db);
    bool have_meta_info = prepareQueryMetaInfo_(subquery_info,
                                                "ID_ParentMolecule");
    while (query.next())
    {
      String accession = query.value("accession").toString();
      ID::ParentMolecule parent(accession);
      int molecule_type_index = query.value("molecule_type_id").toInt() - 1;
      parent.molecule_type = ID::MoleculeType(molecule_type_index);
      parent.sequence = query.value("sequence").toString();
      parent.description = query.value("description").toString();
      parent.coverage = query.value("coverage").toDouble();
      parent.is_decoy = query.value("is_decoy").toInt();
      if (have_meta_info)
      {
        handleQueryMetaInfo_(subquery_info, parent, query.value("id").toInt());
      }
      id_data_.registerParentMolecule(parent);
    }
  }


  void OMSFile::load(const String& filename, IdentificationData& id_data)
  {
    OMSFileLoad helper(filename, id_data);
    helper.loadInputFiles();
    helper.loadScoreTypes();
    helper.loadDataProcessingSoftwares();
    helper.loadDBSearchParams();
    helper.loadDataProcessingSteps();
    helper.loadParentMolecules();
  }

}
