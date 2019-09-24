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

#include <sqlite3.h>

using namespace std;

namespace OpenMS
{
  namespace Sql = Internal::SqliteHelper;

  int version_number = 1;

  static int callback(void* /* NotUsed */, int argc, char** argv, char** azColName)
  {
    int i;
    for (i = 0; i < argc; i++)
    {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
  }


  void OMSFile::storeVersionAndDate_(SqliteConnector& con)
  {
    String sql_create = "CREATE TABLE version ("    \
      "OMSFile INT NOT NULL, "                      \
      "date TEXT NOT NULL, "                        \
      "OpenMS TEXT, "                               \
      "build_date TEXT);";
    con.executeStatement(sql_create);
    stringstream sql_insert;
    sql_insert << "INSERT INTO version VALUES ("
               << version_number << ", "
               << "datetime('now'), '"
               << VersionInfo::getVersion() << "', '"
               // this uses a non-standard date/time format:
               << VersionInfo::getTime() << "');";
    con.executeStatement(sql_insert);
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


  void OMSFile::storeParentMolecules_(const IdentificationData& id_data,
                                      SqliteConnector& con)
  {
    if (id_data.getParentMolecules().empty()) return;
    String sql_create = "CREATE TABLE ID_ParentMolecule ("              \
      "id INTEGER PRIMARY KEY NOT NULL, "                               \
      "accession TEXT UNIQUE NOT NULL, "                                \
      "molecule_type TEXT NOT NULL CHECK (molecule_type IN ('PRO', 'COM', 'RNA')), " \
      "sequence TEXT, "                                                 \
      "description TEXT, "                                              \
      "coverage REAL, "                                    \
      "is_decoy NUMERIC NOT NULL CHECK (is_decoy in (0, 1)) DEFAULT 0);";
    con.executeStatement(sql_create);
    for (const IdentificationData::ParentMolecule& parent :
           id_data.getParentMolecules())
    {
      stringstream sql_insert;
      sql_insert << "INSERT INTO ID_ParentMolecule VALUES ("
                 // use address as primary key:
                 << Int64(&parent) << ", '"
                 << parent.accession << "', '"
                 << getMoleculeTypeAbbrev_(parent.molecule_type) << "', '"
                 << parent.sequence << "', '"
                 << parent.description << "', "
                 << parent.coverage << ", "
                 << int(parent.is_decoy) << ");";
      con.executeStatement(sql_insert);
    }
  }


  void OMSFile::store(const String& filename, const IdentificationData& id_data)
  {
    // delete output file if present:
    File::remove(filename);

    // open database:
    SqliteConnector con(filename.c_str());
    // generally, create tables only if we have data to write - no empty ones!

    storeVersionAndDate_(con);

    storeParentMolecules_(id_data, con);

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


  void OMSFile::loadParentMolecules_(SqliteConnector& con,
                                     IdentificationData& id_data)
  {
    if (!con.tableExists("ID_ParentMolecule")) return;

    sqlite3_stmt* result;
    con.prepareStatement(&result, "SELECT * FROM ID_ParentMolecule;");
    sqlite3_step(result);
    while (sqlite3_column_type(result, 0) != SQLITE_NULL)
    {
      String accession = sqlite3_column_text(result, 1);
      String mt_abbrev = sqlite3_column_text(result, 2);
      IdentificationData::MoleculeType molecule_type =
        getMoleculeTypeFromAbbrev_(mt_abbrev);
      String sequence = sqlite3_column_text(result, 3);
      String description = sqlite3_column_text(result, 4);
      double coverage = sqlite3_column_double(result, 5);
      bool is_decoy = sqlite3_column_int(result, 6);
      IdentificationData::ParentMolecule parent(accession, molecule_type,
                                                sequence, description, coverage,
                                                is_decoy);
      id_data.registerParentMolecule(parent);
      sqlite3_step(result);
    }
    sqlite3_finalize(result);
  }


  void OMSFile::load(const String& filename, IdentificationData& id_data)
  {
    SqliteConnector con(filename.c_str());

    loadParentMolecules_(con, id_data);
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

