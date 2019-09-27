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
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <QtSql/QSqlQuery>
// strangely, this is needed for type conversions in "QSqlQuery::bindValue":
#include <QtSql/QSqlQueryModel>

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
    query.prepare("INSERT INTO ID_ScoreType VALUES ("       \
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


  IdentificationData::MoleculeType OMSFile::getMoleculeTypeFromAbbrev_(const String& abbrev)
  {
    if (abbrev == "PRO") return IdentificationData::MoleculeType::PROTEIN;
    if (abbrev == "COM") return IdentificationData::MoleculeType::COMPOUND;
    if (abbrev == "RNA") return IdentificationData::MoleculeType::RNA;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "invalid abbreviation for a molecule type",
                                  abbrev);
  }


  // currently not needed:
  // CVTerm OMSFile::loadCVTerm_(int id, QSqlDatabase& db)
  // {
  //   // this assumes that the "CVTerm" table exists!
  //   QSqlQuery query(db);
  //   query.setForwardOnly(true);
  //   QString sql_select = "SELECT * FROM CVTerm WHERE id = " + QString(id);
  //   if (!query.exec(sql_select) || !query.next())
  //   {
  //     raiseDBError_(model.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
  //                   "error reading from database");
  //   }
  //   return CVTerm(query.value("accession").toString(),
  //                 query.value("name").toString(),
  //                 query.value("cv_identifier_ref").toString());
  // }


  void OMSFile::loadScoreTypes_(IdentificationData& id_data, QSqlDatabase& db)
  {
    if (!tableExists_(db, "ID_ScoreType")) return;
    if (!tableExists_(db, "CVTerm")) // every score type is a CV term
    {
      String msg = "required database table 'CVTerm' not found";
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ScoreType JOIN CVTerm "   \
                    "ON ID_ScoreType.cv_term_id = CVTerm.id"))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      CVTerm cv_term(query.value("accession").toString(),
                     query.value("name").toString(),
                     query.value("cv_identifier_ref").toString());
      bool higher_better = query.value("higher_better").toInt();
      IdentificationData::ScoreType score_type(cv_term, higher_better);
      id_data.registerScoreType(score_type);
    }
  }


  void OMSFile::loadParentMolecules_(IdentificationData& id_data,
                                     QSqlDatabase& db)
  {
    if (!tableExists_(db, "ID_ParentMolecule")) return;

    QSqlQuery query(db);
    query.setForwardOnly(true);
    if (!query.exec("SELECT * FROM ID_ParentMolecule"))
    {
      raiseDBError_(query.lastError(), db, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "error reading from database");
    }
    while (query.next())
    {
      String accession = query.value("accession").toString();
      IdentificationData::ParentMolecule parent(accession);
      parent.molecule_type =
        getMoleculeTypeFromAbbrev_(query.value("molecule_type").toString());
      parent.sequence = query.value("sequence").toString();
      parent.description = query.value("description").toString();
      parent.coverage = query.value("coverage").toDouble();
      parent.is_decoy = query.value("is_decoy").toInt();
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

      loadScoreTypes_(id_data, db);

      loadParentMolecules_(id_data, db);
    }

    QSqlDatabase::removeDatabase(connection);
  }

}
