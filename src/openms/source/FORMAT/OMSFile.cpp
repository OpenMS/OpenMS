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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer, Oliver Alka, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <QtSql/QSqlQuery>
#include <QtSql/QSqlQueryModel>
#include <QtSql/QSqlRecord>

using namespace std;

namespace OpenMS
{
  int version_number = 1;

  void OMSFile::raiseDBError_(const QSqlError& error, QSqlDatabase& db,
                               int line, const char* function,
                               const String& context)
  {
    // clean-up:
    db.close();
    QSqlDatabase::removeDatabase(db.connectionName());

    String msg = context + ": " + error.text();
    throw Exception::FailedAPICall(__FILE__, line, function, msg);
  }


  bool OMSFile::tableExists_(QSqlDatabase& db, const String& name)
  {
    return db.tables(QSql::Tables).contains(name.toQString());
  }


  void OMSFile::storeVersionAndDate_(QSqlDatabase& db)
  {
    QSqlQuery query(db);
    QString sql_create =
      "CREATE TABLE version ("                           \
      "OMSFile INT NOT NULL, "                           \
      "date TEXT NOT NULL, "                             \
      "OpenMS TEXT, "                                    \
      "build_date TEXT)";
    if (!query.exec(sql_create))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error creating database table");
    }
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
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error updating database");
    }
  }


  String OMSFile::getMoleculeTypeAbbrev_(IdentificationData::MoleculeType molecule_type)
  {
    switch (molecule_type)
    {
    case IdentificationData::MoleculeType::PROTEIN:
      return "PRO";
    case IdentificationData::MoleculeType::COMPOUND:
      return "COM";
    case IdentificationData::MoleculeType::RNA:
      return "RNA";
    default:
      return "";
    }
  }


  int OMSFile::storeCVTerm_(const CVTerm& cv_term, QSqlDatabase& db)
  {
    // this assumes the "CVTerm" table exists already!
    QSqlQuery query(db);
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
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
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
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error querying database");
    }
    return query.value(0).toInt();
  }


  void OMSFile::storeScoreTypes_(const IdentificationData& id_data,
                                 QSqlDatabase& db)
  {
    if (id_data.getScoreTypes().empty()) return;

    QSqlQuery query(db);
    QString sql_create =
      "CREATE TABLE CVTerm ("                               \
      "id INTEGER PRIMARY KEY NOT NULL, "                   \
      "accession TEXT UNIQUE, "                             \
      "name TEXT NOT NULL, "                                \
      "cv_identifier_ref TEXT, "                            \
      // does this constrain "name" if "accession" is NULL?
      "UNIQUE (accession, name))";
    // @TODO: add support for unit and value
    if (!query.exec(sql_create))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error creating database table");
    }
    sql_create =
      "CREATE TABLE ID_ScoreType ("                                     \
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "cv_term_id INTEGER NOT NULL, "                                   \
      "higher_better NUMERIC NOT NULL CHECK (higher_better in (0, 1)), " \
      "FOREIGN KEY (cv_term_id) REFERENCES CVTerm (id))";
    if (!query.exec(sql_create))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error creating database table");
    }
    query.prepare("INSERT INTO ID_ScoreType VALUES ("  \
                  ":id, "                                   \
                  ":cv_term_id, "                           \
                  ":higher_better)");
    for (const IdentificationData::ScoreType& score_type :
           id_data.getScoreTypes())
    {
      int cv_id = storeCVTerm_(score_type.cv_term, db);
      query.bindValue(":id", qint64(&score_type)); // use address as primary key
      query.bindValue(":cv_term_id", cv_id);
      query.bindValue(":higher_better", int(score_type.higher_better));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error updating database");
      }
    }
  }


  void OMSFile::storeParentMolecules_(const IdentificationData& id_data,
                                      QSqlDatabase& db)
  {
    if (id_data.getParentMolecules().empty()) return;

    QSqlQuery query(db);
    QString sql_create =
      "CREATE TABLE ID_ParentMolecule ("                                \
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "accession TEXT UNIQUE NOT NULL, "                                \
      "molecule_type TEXT NOT NULL CHECK (molecule_type IN ('PRO', 'COM', 'RNA')), " \
      "sequence TEXT, "                                                 \
      "description TEXT, "                                              \
      "coverage REAL, "                                                 \
      "is_decoy NUMERIC NOT NULL CHECK (is_decoy in (0, 1)) DEFAULT 0)";
    if (!query.exec(sql_create))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error creating database table");
    }
    query.prepare("INSERT INTO ID_ParentMolecule VALUES ("  \
                  ":id, "                                   \
                  ":accession, "                            \
                  ":molecule_type, "                        \
                  ":sequence, "                             \
                  ":description, "                          \
                  ":coverage, "                             \
                  ":is_decoy)");
    for (const IdentificationData::ParentMolecule& parent :
           id_data.getParentMolecules())
    {
      query.bindValue(":id", qint64(&parent)); // use address as primary key
      query.bindValue(":accession", parent.accession.toQString());
      query.bindValue(":molecule_type",
                      getMoleculeTypeAbbrev_(parent.molecule_type).toQString());
      query.bindValue(":sequence", parent.sequence.toQString());
      query.bindValue(":description", parent.description.toQString());
      query.bindValue(":coverage", parent.coverage);
      query.bindValue(":is_decoy", int(parent.is_decoy));
      if (!query.exec())
      {
        raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error updating database");
      }
    }
  }


  void OMSFile::store(const String& filename, const IdentificationData& id_data)
  {
    // delete output file if present:
    File::remove(filename);

    // open database:
    QString connection = "store_" + filename.toQString();
    { // extra scope to avoid warning from "QSqlDatabase::removeDatabase"
      QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connection);
      db.setDatabaseName(filename.toQString());
      if (!db.open())
      {
        raiseDBError_(db.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error opening SQLite database");
      }
      // generally, create tables only if we have data to write - no empty ones!

      storeVersionAndDate_(db);

      storeScoreTypes_(id_data, db);

      storeParentMolecules_(id_data, db);
    }

    QSqlDatabase::removeDatabase(connection);
  }

/*

      // IdentifiedCompound
      "CREATE TABLE IDENTIFIEDCOMPOUND(" \
      "id INT PRIMARY KEY NOT NULL," \
      "identifier TEXT UNIQUE NOT NULL," \
      "formula TEXT DEFAULT NULL," \
      "name TEXT DEFAULT NULL," \
      "smile TEXT DEFAULT NULL," \
      "inchi TEXT DEFAULT NULL);" \

      // IdentifiedSequence
      "CREATE TABLE IDENTIFIEDSEQUENCE(" \
      "id INT PRIMARY KEY NOT NULL," \
      "sequence TEXT NOT NULL," \
      "seq_type TEXT NOT NULL CHECK (seq_type IN ('PRO', 'RNA'))," \
      "UNIQUE(sequence, seq_type));" \

      // ParentMoleculeGroup
      "CREATE TABLE PARENTMOLECULEGROUP(" \
      "id INT PRIMARY KEY NOT NULL," \
      "scoretype INT NOT NULL," \
      "scorevalue REAL NOT NULL);" \

      // ParentMoleculeGrouping
      "CREATE TABLE PARENTMOLECULEGROUPING(" \
      "id INT PRIMARY KEY NOT NULL," \
      "label TEXT NOT NULL);" \

      // Scores
      "CREATE TABLE SCORES(" \
      "id INT PRIMARY KEY NOT NULL," \
      "scorename TEXT NOT NULL," \
      "cvterm TEXT NOT NULL);" \

      // AppliedProcessingSteps
      "CREATE TABLE APPLIEDPROCESSINGSTEPS(" \
      "id INT PRIMARY KEY NOT NULL," \
      "dateprocessingstep INT NOT NULL);" \

      // DataProcessingStep
      "CREATE TABLE DATAPROCESSINGSTEP(" \
      "id INT PRIMARY KEY NOT NULL," \
      "software TEXT NOT NULL," \
      "inputfile TEXT NOT NULL," \
      "datetime TEXT NOT NULL," \
      "actions TEXT NOT NULL," \
      "dbsearchparam INT NULL);" \

      // DBSearchParam
      "CREATE TABLE DBSEARCHPARAM(" \
      "id INT PRIMARY KEY NOT NULL," \
      "moleculetype TEXT NOT NULL," \
      "masstype TEXT NOT NULL," \
      "database TEXT NOT NULL," \
      "dbversion TEXT NOT NULL," \
      "taxonomy TEXT," \
      "charges TEXT NOT NULL," \
      "fixedmods TEXT NOT NULL," \
      "variablemods TEXT NOT NULL," \
      "precursormasstolerance INT NOT NULL," \
      "fragmentmasstolerance INT NOT NULL," \
      "precursorrmasstoleranceppm BIT NOT NULL," \
      "fragmentmasstoleranceppm BIT NOT NULL," \
      "digestionenyzme TEXT NOT NULL," \
      "missedcleavages INT NOT NULL," \
      "minlength INT NOT NULL," \
      "maxlength INT NOT NULL);" \

      // ProcessingSoftware
      "CREATE TABLE PROCESSINGSOFTWARE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "NAME TEXT NOT NULL," \
      "VERSION TEXT NOT NULL);" \

      // DataQuery
      "CREATE TABLE DATAQUERY(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "DATA_ID TEXT NOT NULL," \
      "INPUT_FILE_REF TEXT," \
      "RT REAL NOT NULL," \
      "MZ REAL NOT NULL);" \

      // MoleculeQueryMatch
      "CREATE TABLE MOLECULEQUERYMATCH(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "IDENTIFIED_MOLECULE_REF TEXT NOT NULL," \
      "MOLECULETYPE TEXT NOT NULL," \
      "DATAQUERYREF TEXT NOT NULL," \
      "CHARGE INT NOT NULL);" \

      // PeakAnnotations
      "CREATE TABLE PEAKANNOTATIONS(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "ANNOTATIONS TEXT NOT NULL," \
      "CHARGE INT NOT NULL," \
      "MZ REAL NOT NULL," \
      "INTENSITY REAL NOT NULL);" \

      // ParentMoleculeMatch
      "CREATE TABLE PARENTMOLECULEMATCH(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "START_POS TEXT NOT NULL," \
      "END_POS TEXT NOT NULL," \
      "LEFT_NEIGHBOR TEXT NOT NULL," \
      "RIGHT_NEIGHBOR TEXT NOT NULL);" \

      // Mapping:  Scores AppliedProcessingSteps
      "CREATE TABLE SCORES_APPLIEDPROCESSINGSTEPS_MAPPING(" \
      "APSSID INT NOT NULL," \
      "SCORETYPEID INT NOT NULL," \
      "SCORE REAL NOT NULL);" \

      // Mapping: ParentMolecule AppliedProcessingStep
      "CREATE TABLE PARENTMOLECULE_APPLIEDPROCESSINGSTEP_MAPPING(" \
      "PMID INT NOT NULL," \
      "APSID INT NOT NULL);" \

      // Mapping: ParentMolecule PrarentMoleculeGrouping
      "CREATE TABLE PARENTMOLECULE_PARENTMOLECULEGROUP_MAPPING(" \
      "PMLID INT PRIMARY KEY NOT NULL," \
      "PMGID INT NOT NULL);" \

      // Mapping: AppliedProcessingStep MoleculeQueryMatch
      "CREATE TABLE APPLIEDPROCESSINGSTEP_MOLECULEQUERYMATCH_MAPPING(" \
      "APSID INT PRIMARY KEY NOT NULL," \
      "MQMID INT NOT NULL," \
      "SCORE REAL NOT NULL," \
      "SCORETYPE TEXT NOT NULL);" \

      // Mapping: PeakAnnotation MoleculeQueryMatch
        // PeakAnnotationID
        // DataProcessingID
        // MoleculeQueryMatchID
      "CREATE TABLE PEAKANNOTATION_MOLECULEQUERYMATCH_MAPPING(" \
      "PAID INT PRIMARY KEY NOT NULL," \
      "DPID INT NOT NULL," \
      "MQMID REAL NOT NULL);" \

      // Mapping: ParentMolecule IdentifiedSequence
        // ParentMoleculeID
        // IdentifiedSequenceID
      "CREATE TABLE PARENTMOLECULE_IDENTIFIEDSEQUENCE_MAPPING(" \
      "PMID INT PRIMARY KEY NOT NULL," \
      "ISID INT NOT NULL);" \

      // Mapping: IdentfiedSquence AppliedProcessingStep Score ScoreType
        // IdentifiedSequenceID
        // AppliedProcessingStep
        // Score
        // ScoreType
      "CREATE TABLE IDS_APS_SCORE_SCORETYPE_MAPPING(" \
      "ISID INT PRIMARY KEY NOT NULL," \
      "APSID INT NOT NULL," \
      "SCORE REAL NOT NULL," \
      "SCORETYPE TEXT NOT NULL);" \

      // Mapping: QueryMatchGroup MoleculeQueryMatch
      // QueryMatchGroupID
      // MoleculeQueryMatchID
      "CREATE TABLE QUERYMATCHGROUP_MOLECULEQUERYMATH(" \
      "QMDID INT PRIMARY KEY NOT NULL," \
      "MQMID INT NOT NULL);"

      // MetaboInterface STRING
      "CREATE TABLE METABOINTERFACE_STRING(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE TEXT NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);"

      // MetaboInterface INT
      "CREATE TABLE METABOINTERFACE_INT(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE INT NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);"

      // MetaboInterface REAL
      "CREATE TABLE METABOINTERFACE_REAL(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE REAL NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);"

      // MetaboInterface SRINGLIST
      "CREATE TABLE METABOINTERFACE_STRINGLIST(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE REAL NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);"

      // MetaboInterface INTLIST
      "CREATE TABLE METABOINTERFACE_INTLIST(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE REAL NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);"

      // MetaboInterface DOUBLELIST
      "CREATE TABLE METABOINTERFACE_DOUBLELIST(" \
      "ID INT NOT NULL," \
      "METAVALUE TEXT NOT NULL," \
      "VALUE REAL NOT NULL," \
      "OBJECT TEXT NOT NULL," \
      "OBJECTREFERENCE INT NOT NULL);" ;

      // TODO: add MetaInfo Tables
      // TODO: how would the ref work best

    // Execute SQL create statement
    conn.executeStatement(create_sql);

    // Prepare insert statements

    // Index maps (memory map)
    */


  IdentificationData::MoleculeType OMSFile::getMoleculeTypeFromAbbrev_(const String& abbrev)
  {
    if (abbrev == "PRO") return IdentificationData::MoleculeType::PROTEIN;
    if (abbrev == "COM") return IdentificationData::MoleculeType::COMPOUND;
    if (abbrev == "RNA") return IdentificationData::MoleculeType::RNA;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "invalid abbreviation for a molecule type",
                                  abbrev);
  }


  void OMSFile::loadParentMolecules_(IdentificationData& id_data,
                                     QSqlDatabase& db)
  {
    if (!tableExists_(db, "ID_ParentMolecule")) return;

    QSqlQueryModel model;
    model.setQuery("SELECT * FROM ID_ParentMolecule", db);
    if (model.lastError().isValid())
    {
      raiseDBError_(model.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                     "error reading from database");
    }
    for (int i = 0; i < model.rowCount(); ++i)
    {
      const QSqlRecord& row = model.record(i);
      String accession = row.value("accession").toString();
      IdentificationData::ParentMolecule parent(accession);
      parent.molecule_type =
        getMoleculeTypeFromAbbrev_(row.value("molecule_type").toString());
      parent.sequence = row.value("sequence").toString();
      parent.description = row.value("description").toString();
      parent.coverage = row.value("coverage").toDouble();
      parent.is_decoy = row.value("is_decoy").toInt();
      id_data.registerParentMolecule(parent);
    }
  }


  void OMSFile::load(const String& filename, IdentificationData& id_data)
  {
    // open database:
    QString connection = "load_" + filename.toQString();
    { // extra scope to avoid warning from "QSqlDatabase::removeDatabase"
      QSqlDatabase db = QSqlDatabase::addDatabase("QSQLITE", connection);
      db.setDatabaseName(filename.toQString());
      if (!db.open())
      {
        raiseDBError_(db.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "error opening SQLite database");
      }

      loadParentMolecules_(id_data, db);
    }

    QSqlDatabase::removeDatabase(connection);
  }

  /*
        // ParentMolecule table
        "CREATE TABLE PARENTMOLECULE(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "ACCESSION TEXT NOT NULL," \
        "SEQUENCE TEXT NOT NULL," \
        "DESCRIPTION TEXT NOT NULL," \
        "COVERAGE REAL NOT NULL," \
        "DECOY BIT NOT NULL);" \

        // ParentMoleculeGroup table
        // TODO Use Mapping table (Molecule to Group)?
        "CREATE TABLE PARENTMOLECULEGROUP(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "SCORETYPEREF INT NOT NULL,"
        "SCORE REAL NOT NULL);" \

        // ParentMolecule to ParentMoleculeGroup
        "CREATE TABLE PARENTMOLECULE_PARENTMOLECULEGROUP_MAPPING(" \
        "PMOLID INT NOT NULL," \
        "PMOLGRPID INT NOT NULL);" \

        // ParentMoleculeGrouping table
        "CREATE TABLE PARENTMOLECULEGROUPING(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "LABEL TEXT NOT NULL);" \

        // AppliedProcessingStep
        "CREATE TABLE APPLIEDPROCESSINGSTEP(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "PROCESSING_STEP_OPT TEXT);"

        // ScoresForAppliedProcessingStep
        "CREATE TABLE SCORESFORAPPLIEDPROCESSINGSTEP(" \
        "APSID INT PRIMARY KEY NOT NULL," \
        "SCORETYPEREF INT NOT NULL," \
        "SCORE REAL NOT NULL);";

        // DataProcessingStep
       "CREATE TABLE DATAPROCESSINGSTEP(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "SOFTWARE_REF TEXT NOT NULL," \
        "INPUT_FILE_REFS TEXT NOT NULL," \ // TODO: vector <InputFileRef>
        "PRIMARY_FILES REAL NOT NULL," \ // TODO: vector <String>
        "DATETIME TEXT NOT NULL," \ // TODO: DateTime
        "ACTIONS REAL NOT NULL);" \ //TODO: std::set< DataProcessing::ProcessingAction >

      // DBSearchParam
      "CREATE TABLE DBSEARCHPARAM(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "MOLECULETYPE TEXT NOT NULL," \
        "MASSTYPE TEXT NOT NULL," \
        "DATABASE TEXT NOT NULL," \
        "DATABASE_VERSION TEXT NOT NULL," \
        "TAXONOMY TEXT NOT NULL," \
        "CHARGES TEXT NOT NULL," \ //TODO: set<int>
      "FIXED_MODS TEXT NOT NULL," \ //TODO: set<String>
      "VARIABLE_MODS TEXT NOT NULL," \ //TODO: set<String>
      "PRECURSOR_MASS_TOLERANCE TEXT NOT NULL," \
        "FRAGMENT_MASS_TOLERANCE TEXT NOT NULL," \
        "PRECURSOR_TOLERANCE_PPM INT NOT NULL," \ //TODO: bool
      "FRAGMENT_TOLERANCE_PPM INT NOT NULL," \ //TODO: bool
      "DIGESTION_ENZYME TEXT NOT NULL," \ //TODO: class DigestionEnzyme
      "MISSED_CLEAVAGES INT NOT NULL," \
        "MIN_LENGTH INT NOT NULL, \
        "MAX_LENGTH INT NOT NULL") \

      // ScoreType
      "CREATE TABLE SCORETYPE(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "CV_TERM TEXT NOT NULL," \ //TODO: CVTerm
      "NAME TEXT NOT NULL," \
        "HIGHER_BETTER INT NOT NULL"); //TODO: bool
      */
}

