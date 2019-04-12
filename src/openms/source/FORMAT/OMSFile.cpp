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
// $Authors: Julianus Pfeuffer, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <sqlite3.h>
#include <OpenMS/FORMAT/SqliteConnector.h>

namespace OpenMS
{

  namespace Sql = Internal::SqliteHelper;

  static int callback(void * /* NotUsed */, int argc, char **argv, char **azColName)
  {
    int i;
    for (i = 0; i < argc; i++)
    {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
  }

  void OMSFile::store(const char* filename, IdentificationData& id_data)
  {

    // delete file if present
    remove(filename);

    // Open database
    SqliteConnector conn(filename);
    //db = conn.getDB();

    // Create SQL structure
    const char* create_sql =
      "CREATE TABLE VERSION(" \
      "ID INT NOT NULL);" \

      // New datastructure

      // IdentifiedCompound
      "CREATE TABLE IDENTIFIEDCOMPOUND(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "IDENTIFIER TEXT NOT NULL," \
      "FORMULA TEXT NOT NULL," \
      "NAME TEXT NOT NULL," \
      "SMILE TEXT NOT NULL," \
      "INCHI TEXT NOT NULL);" \

      // IdentifiedSequence
      "CREATE TABLE IDENTIFIEDSEQUENCE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "SEQUENCE TEXT NOT NULL," \
      "SEQTYPE TEXT NOT NULL);" \

      // ParentMolecule
      "CREATE TABLE PARENTMOLECULE(" \
      "ID INT NOT NULL," \
      "ACCESSION TEXT NOT NULL," \
      "SEQUENCE TEXT NOT NULL," \
      "DECOY BIT NOT NULL," \
      "GROUPID INT NOT NULL," \
      "MOLECULETYPE TEXT NOT NULL," \
      "COVERAGE TEXT NOT NULL);" \

      // ParentMoleculeGroup
      "CREATE TABLE PARENTMOLECULEGROUPS(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "SCORETYPE INT NOT NULL," \
      "SCOREVALUE REAL NOT NULL);" \

      // ParentMoleculeGrouping
      "CREATE TABLE PARENTMOLECULEGROUPING(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "LABEL TEXT NOT NULL);" \

      // Scores
      "CREATE TABLE SCORES(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "SCORENAME TEXT NOT NULL," \
      "CVTERM TEXT NOT NULL);" \

      // AppliedProcessingSteps
      "CREATE TABLE APPLIEDPROCESSINGSTEPS(" \
      "ID INT NOT NULL," \
      "DATEPROCESSINGSTEP INT NOT NULL);" \

      // DataProcessingStep
      "CREATE TABLE DATAPROCESSINGSTEP(" \
      "ID INT NOT NULL," \
      "SOFTWARE TEXT NOT NULL," \
      "INPUTFILE TEXT NOT NULL," \
      "DATATIME TEXT NOT NULL," \
      "ACTIONS TEXT NOT NULL," \
      "DBSEARCHPARAM INT NULL);" \

      // DBSearchParam
      "CREATE TABLE DBSEARCHPARAM(" \
      "ID INT NOT NULL," \
      "MOLECULETYPE TEXT NOT NULL," \
      "MASSTYPE TEXT NOT NULL," \
      "DATABASE TEXT NOT NULL," \
      "DBVERSION TEXT NOT NULL," \
      "TAXONOMY TEXT," \
      "CHARGES TEXT NOT NULL," \
      "FIXED_MODIFICATIONS TEXT NOT NULL," \
      "VARIABLE_MODIFICATIONS TEYT NOT NULL," \
      "PRECURSORMASSTOLERANCE INT NOT NULL," \
      "FRAGMENTMASSTOLERANCE INT NOT NULL," \
      "PRECURSORRMASSTOLERANCEPPM BIT NOT NULL," \
      "FRAGMENTMASSTOLERANCEPPM BIT NOT NULL," \
      "DIGESTION_ENYZME TEXT NOT NULL," \
      "MISSED_CLEAVAGES INT NOT NULL," \
      "MIN_LENGTH INT NOT NULL," \
      "MAX_LENGTH INT NOT NULL);" \

      // ProcessingSoftware
      "CREATE TABLE PROCESSINGSOFTWARE(" \
      "ID INT NOT NULL," \
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

  }

  void OMSFile::load(const char* filename, IdentificationData& id_data)
  {

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

