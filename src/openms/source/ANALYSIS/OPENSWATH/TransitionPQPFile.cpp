// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

namespace OpenMS
{

  TransitionPQPFile::TransitionPQPFile() :
    TransitionTSVFile()
  {
  }

  TransitionPQPFile::~TransitionPQPFile()
  {
  }

  static int callback(void * /* NotUsed */, int argc, char **argv, char **azColName){
    int i;
    for (i=0; i<argc; i++)
    {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return(0);
  }

  void TransitionPQPFile::readPQPInput_(const char* filename, std::vector<TSVTransition>& transition_list, bool legacy_traml_id)
  {
    sqlite3 *db;
    sqlite3_stmt * cntstmt;
    sqlite3_stmt * stmt;
    int rc;
    std::string select_sql;

    // Use legacy TraML identifiers for precursors (transition_group_id) and transitions (transition_name)?
    std::string traml_id = "ID";
    if (legacy_traml_id)
    {
      traml_id = "TRAML_ID";
    }

    // Open database
    rc = sqlite3_open(filename, &db);
    if ( rc )
    {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
    }

    // Count transitions
    sqlite3_prepare_v2(db, "SELECT COUNT(*) FROM TRANSITION;", -1, &cntstmt, nullptr);
    sqlite3_step( cntstmt );
    int num_transitions = sqlite3_column_int( cntstmt, 0 );
    sqlite3_finalize(cntstmt);

    // Get peptides
    select_sql = "SELECT " \
                  "PRECURSOR.PRECURSOR_MZ AS precursor, " \
                  "TRANSITION.PRODUCT_MZ AS product, " \
                  "PRECURSOR.LIBRARY_RT AS rt_calibrated, " \
                  "TRANSITION." + traml_id + " AS transition_name, " \
                  "-1 AS CE, " \
                  "TRANSITION.LIBRARY_INTENSITY AS library_intensity, " \
                  "PRECURSOR." + traml_id + " AS group_id, " \
                  "TRANSITION.DECOY AS decoy, " \
                  "PEPTIDE.UNMODIFIED_SEQUENCE AS PeptideSequence, " \
                  "PROTEIN_AGGREGATED.PROTEIN_ACCESSION AS ProteinName, " \
                  "NULL AS Annotation, " \
                  "PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName, " \
                  "NULL AS CompoundName, " \
                  "NULL AS SMILES, " \
                  "NULL AS SumFormula, " \
                  "PRECURSOR.CHARGE AS precursor_charge, " \
                  "PRECURSOR.GROUP_LABEL AS peptide_group_label, " \
                  "NULL AS label_type, " \
                  "TRANSITION.CHARGE AS fragment_charge, " \
                  "TRANSITION.ORDINAL AS fragment_nr, " \
                  "NULL AS fragment_mzdelta, " \
                  "NULL AS fragment_modification, " \
                  "TRANSITION.TYPE AS fragment_type, " \
                  "NULL AS uniprot_id, " \
                  "TRANSITION.DETECTING AS detecting_transition, " \
                  "TRANSITION.IDENTIFYING AS identifying_transition, " \
                  "TRANSITION.QUANTIFYING AS quantifying_transition, " \
                  "NULL AS peptidoforms " \
                  "FROM PRECURSOR " \
                  "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID " \
                  "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID " \
                  "INNER JOIN " \
                  "(SELECT PEPTIDE_ID, GROUP_CONCAT(PROTEIN_ACCESSION,';') AS PROTEIN_ACCESSION FROM PROTEIN INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID GROUP BY PEPTIDE_ID) AS PROTEIN_AGGREGATED ON PEPTIDE.ID = PROTEIN_AGGREGATED.PEPTIDE_ID ";

    // Get compounds
    select_sql += "UNION SELECT " \
                  "PRECURSOR.PRECURSOR_MZ AS precursor, " \
                  "TRANSITION.PRODUCT_MZ AS product, " \
                  "PRECURSOR.LIBRARY_RT AS rt_calibrated, " \
                  "TRANSITION." + traml_id + " AS transition_name, " \
                  "-1 AS CE, " \
                  "TRANSITION.LIBRARY_INTENSITY AS library_intensity, " \
                  "PRECURSOR." + traml_id + " AS group_id, " \
                  "TRANSITION.DECOY AS decoy, " \
                  "NULL AS PeptideSequence, " \
                  "NULL AS ProteinName, " \
                  "NULL AS Annotation, " \
                  "NULL AS FullPeptideName, " \
                  "COMPOUND.COMPOUND_NAME AS CompoundName, " \
                  "COMPOUND.SMILES AS SMILES, " \
                  "COMPOUND.SUM_FORMULA AS SumFormula, " \
                  "PRECURSOR.CHARGE AS precursor_charge, " \
                  "PRECURSOR.GROUP_LABEL AS peptide_group_label, " \
                  "NULL AS label_type, " \
                  "TRANSITION.CHARGE AS fragment_charge, " \
                  "TRANSITION.ORDINAL AS fragment_nr, " \
                  "NULL AS fragment_mzdelta, " \
                  "NULL AS fragment_modification, " \
                  "TRANSITION.TYPE AS fragment_type, " \
                  "NULL AS uniprot_id, " \
                  "TRANSITION.DETECTING AS detecting_transition, " \
                  "TRANSITION.IDENTIFYING AS identifying_transition, " \
                  "TRANSITION.QUANTIFYING AS quantifying_transition, " \
                  "NULL AS peptidoforms " \
                  "FROM PRECURSOR " \
                  "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID " \
                  "INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID; ";

    // Execute SQL select statement
    sqlite3_prepare_v2(db, select_sql.c_str(), -1, &stmt, nullptr);
    sqlite3_step( stmt );

    Size progress = 0;
    startProgress(0, num_transitions, "reading PQP file");
    // Convert SQLite data to TSVTransition data structure
    while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
    {
      setProgress(progress++);
      TSVTransition mytransition;

      if (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        mytransition.precursor = sqlite3_column_double( stmt, 0 );
      }
      if (sqlite3_column_type( stmt, 1 ) != SQLITE_NULL)
      {
        mytransition.product = sqlite3_column_double( stmt, 1 );
      }
      if (sqlite3_column_type( stmt, 2 ) != SQLITE_NULL)
      {
        mytransition.rt_calibrated = sqlite3_column_double( stmt, 2 );
      }
      if (sqlite3_column_type( stmt, 3 ) != SQLITE_NULL)
      {
        mytransition.transition_name = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 3 )));
      }
      if (sqlite3_column_type( stmt, 4 ) != SQLITE_NULL)
      {
        mytransition.CE = sqlite3_column_double( stmt, 4 );
      }
      if (sqlite3_column_type( stmt, 5 ) != SQLITE_NULL)
      {
        mytransition.library_intensity = sqlite3_column_double( stmt, 5 );
      }
      if (sqlite3_column_type( stmt, 6 ) != SQLITE_NULL)
      {
        mytransition.group_id = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 6 )));
      }
      if (sqlite3_column_type( stmt, 7 ) != SQLITE_NULL)
      {
        mytransition.decoy = sqlite3_column_int( stmt, 7 );
      }
      if (sqlite3_column_type( stmt, 8 ) != SQLITE_NULL)
      {
        mytransition.PeptideSequence = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 8 )));
      }
      if (sqlite3_column_type( stmt, 9 ) != SQLITE_NULL)
      {
        mytransition.ProteinName = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 9 )));
      }
      if (sqlite3_column_type( stmt, 10 ) != SQLITE_NULL)
      {
        mytransition.Annotation = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 10 )));
      }
      if (sqlite3_column_type( stmt, 11 ) != SQLITE_NULL)
      {
        mytransition.FullPeptideName = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 11 )));
      }
      if (sqlite3_column_type( stmt, 12 ) != SQLITE_NULL)
      {
        mytransition.CompoundName = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 12 )));
      }
      if (sqlite3_column_type( stmt, 13 ) != SQLITE_NULL)
      {
        mytransition.SMILES = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 13 )));
      }
      if (sqlite3_column_type( stmt, 14 ) != SQLITE_NULL)
      {
        mytransition.SumFormula = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 14 )));
      }
      if (sqlite3_column_type( stmt, 15 ) != SQLITE_NULL)
      {
        mytransition.precursor_charge = sqlite3_column_int( stmt, 15 );
      }
      if (sqlite3_column_type( stmt, 16 ) != SQLITE_NULL)
      {
        mytransition.peptide_group_label = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 16 )));
      }
      if (sqlite3_column_type( stmt, 17 ) != SQLITE_NULL)
      {
        mytransition.label_type = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 17 )));
      }
      if (sqlite3_column_type( stmt, 18 ) != SQLITE_NULL)
      {
        mytransition.fragment_charge = sqlite3_column_int( stmt, 18 );
      }
      if (sqlite3_column_type( stmt, 19 ) != SQLITE_NULL)
      {
        mytransition.fragment_nr = sqlite3_column_int( stmt, 19 );
      }
      if (sqlite3_column_type( stmt, 20 ) != SQLITE_NULL)
      {
        mytransition.fragment_mzdelta = sqlite3_column_double( stmt, 20 );
      }
      if (sqlite3_column_type( stmt, 21 ) != SQLITE_NULL)
      {
        mytransition.fragment_modification = sqlite3_column_int( stmt, 21 );
      }
      if (sqlite3_column_type( stmt, 22 ) != SQLITE_NULL)
      {
        mytransition.fragment_type = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 22 )));
      }
      if (sqlite3_column_type( stmt, 23 ) != SQLITE_NULL)
      {
        mytransition.uniprot_id = std::string(reinterpret_cast<const char*>(sqlite3_column_text( stmt, 23 )));
      }
      if (sqlite3_column_type( stmt, 24 ) != SQLITE_NULL)
      {
        mytransition.detecting_transition = sqlite3_column_int( stmt, 24 );
      }
      if (sqlite3_column_type( stmt, 25 ) != SQLITE_NULL)
      {
        mytransition.identifying_transition = sqlite3_column_int( stmt, 25 );
      }
      if (sqlite3_column_type( stmt, 26 ) != SQLITE_NULL)
      {
        mytransition.quantifying_transition = sqlite3_column_int( stmt, 26 );
      }

      transition_list.push_back(mytransition);
      sqlite3_step( stmt );
    }
    endProgress();

    sqlite3_finalize(stmt);
    sqlite3_close(db);

  }

  void TransitionPQPFile::writePQPOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    sqlite3 *db;
    char *zErrMsg = nullptr;
    int  rc;

    // delete file if present
    remove(filename);

    // Open database
    rc = sqlite3_open(filename, &db);
    if ( rc )
    {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
    }

    // Create SQL structure
    const char* create_sql =
      // protein table
      // OpenSWATH proteomics workflows
      "CREATE TABLE PROTEIN(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "PROTEIN_ACCESSION TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // peptide_protein_mapping table
      // OpenSWATH proteomics workflows
      "CREATE TABLE PEPTIDE_PROTEIN_MAPPING(" \
      "PEPTIDE_ID INT NOT NULL," \
      "PROTEIN_ID INT NOT NULL);" \

      // peptide table
      // OpenSWATH proteomics workflows
      "CREATE TABLE PEPTIDE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "UNMODIFIED_SEQUENCE TEXT NOT NULL," \
      "MODIFIED_SEQUENCE TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // precursor_peptide_mapping table
      // OpenSWATH proteomics workflows
      "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(" \
      "PRECURSOR_ID INT NOT NULL," \
      "PEPTIDE_ID INT NOT NULL);" \

      // compound table
      // OpenSWATH metabolomics workflows
      "CREATE TABLE COMPOUND(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "COMPOUND_NAME TEXT NOT NULL," \
      "SUM_FORMULA TEXT NOT NULL," \
      "SMILES TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // precursor_compound_mapping table
      // OpenSWATH metabolomics workflows
      "CREATE TABLE PRECURSOR_COMPOUND_MAPPING(" \
      "PRECURSOR_ID INT NOT NULL," \
      "COMPOUND_ID INT NOT NULL);" \

      // precursor table
      "CREATE TABLE PRECURSOR(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "TRAML_ID TEXT NULL," \
      "GROUP_LABEL TEXT NULL," \
      "PRECURSOR_MZ REAL NOT NULL," \
      "CHARGE INT NULL," \
      "LIBRARY_INTENSITY REAL NULL," \
      "LIBRARY_RT REAL NULL," \
      "DECOY INT NOT NULL);" \

      // transition_precursor_mapping table
      "CREATE TABLE TRANSITION_PRECURSOR_MAPPING(" \
      "TRANSITION_ID INT NOT NULL," \
      "PRECURSOR_ID INT NOT NULL);" \

      // transition_peptide_mapping table
      // IPF proteomics workflows
      "CREATE TABLE TRANSITION_PEPTIDE_MAPPING(" \
      "TRANSITION_ID INT NOT NULL," \
      "PEPTIDE_ID INT NOT NULL);" \

      // transition table
      "CREATE TABLE TRANSITION(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "TRAML_ID TEXT NULL," \
      "PRODUCT_MZ REAL NOT NULL," \
      "CHARGE INT NULL," \
      "TYPE CHAR(1) NULL," \
      "ORDINAL INT NULL," \
      "DETECTING INT NOT NULL," \
      "IDENTIFYING INT NOT NULL," \
      "QUANTIFYING INT NOT NULL," \
      "LIBRARY_INTENSITY REAL NULL," \
      "DECOY INT NOT NULL);";

    // Execute SQL create statement
    rc = sqlite3_exec(db, create_sql, callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Prepare insert statements

    // Index maps
    std::vector<std::string> group_vec, peptide_vec, compound_vec, protein_vec;
    std::map<std::string, int > group_map, peptide_map, compound_map, protein_map;
    std::map<int,double> precursor_mz_map;
    std::map<int,bool> precursor_decoy_map;

    std::stringstream insert_transition_sql, insert_transition_peptide_mapping_sql, insert_transition_precursor_mapping_sql;
    insert_transition_sql.precision(11);

    // OpenSWATH: Loop through TargetedExperiment to generate index maps for peptides
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      OpenMS::TargetedExperiment::Peptide peptide = targeted_exp.getPeptides()[i];
      std::string peptide_sequence = TargetedExperimentHelper::getAASequence(peptide).toUniModString();
      peptide_vec.push_back(peptide_sequence);
      group_vec.push_back(peptide.id);
    }

    // OpenSWATH: Loop through TargetedExperiment to generate index maps for compounds
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      OpenMS::TargetedExperiment::Compound compound = targeted_exp.getCompounds()[i];
      compound_vec.push_back(compound.id);
      group_vec.push_back(compound.id);
    }

    // OpenSWATH: Group set must be unique
    boost::erase(group_vec, boost::unique<boost::return_found_end>(boost::sort(group_vec)));
    int group_map_idx = 0;
    for (auto const & x : group_vec) { group_map[x] = group_map_idx; group_map_idx++; }

    // IPF: Loop through all transitions and generate peptidoform data structures
    std::vector<TransitionPQPFile::TSVTransition > transitions;
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      TransitionPQPFile::TSVTransition transition = convertTransition_(&targeted_exp.getTransitions()[i], targeted_exp);
      transitions.push_back(transition);

      std::copy( transition.peptidoforms.begin(), transition.peptidoforms.end(), std::inserter( peptide_vec, peptide_vec.end() ) );

      int group_set_index = group_map[transition.group_id];

      if (precursor_mz_map.find(group_set_index) == precursor_mz_map.end())
      {
        precursor_mz_map[group_set_index] = transition.precursor;
      }
      if (precursor_decoy_map.find(group_set_index) == precursor_decoy_map.end())
      {
        if (transition.detecting_transition == 1)
        {
          precursor_decoy_map[group_set_index] = transition.decoy;
        }
      }
    }

    // OpenSWATH: Peptide and compound sets must be unique
    boost::erase(peptide_vec, boost::unique<boost::return_found_end>(boost::sort(peptide_vec)));
    int peptide_map_idx = 0;
    for (auto const & x : peptide_vec) { peptide_map[x] = peptide_map_idx; peptide_map_idx++; }

    boost::erase(compound_vec, boost::unique<boost::return_found_end>(boost::sort(compound_vec)));
    int compound_map_idx = 0;
    for (auto const & x : compound_vec) { compound_map[x] = compound_map_idx; compound_map_idx++; }

    // OpenSWATH: Loop through TargetedExperiment to generate index maps for proteins
    for (Size i = 0; i < targeted_exp.getProteins().size(); i++)
    {
      OpenMS::TargetedExperiment::Protein protein = targeted_exp.getProteins()[i];
      protein_vec.push_back(protein.id);
    }

    // OpenSWATH: Protein set must be unique
    boost::erase(protein_vec, boost::unique<boost::return_found_end>(boost::sort(protein_vec)));
    int protein_map_idx = 0;
    for (auto const & x : protein_vec) { protein_map[x] = protein_map_idx; protein_map_idx++; }

    // OpenSWATH: Prepare transition inserts
    for (Size i = 0; i < transitions.size(); i++)
    {
      TransitionPQPFile::TSVTransition transition = transitions[i];

      // IPF: Generate transition-peptide mapping tables (one identification transition can map to multiple peptidoforms)
      for (Size j = 0; j < transition.peptidoforms.size(); j++)
      {
        insert_transition_peptide_mapping_sql << "INSERT INTO TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID, PEPTIDE_ID) VALUES (" << i << "," << peptide_map[transition.peptidoforms[j]] << "); ";
      }

      // OpenSWATH: Associate transitions with their precursors
      insert_transition_precursor_mapping_sql << "INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (" << i << "," << group_map[transition.group_id] << "); ";

      std::string transition_charge = "NULL"; // workaround for compounds with missing charge
      if (transition.fragment_charge != "NA")
      {
        transition_charge = transition.fragment_charge;
      }

      // OpenSWATH: Insert transition data
      insert_transition_sql << "INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (" << i << ",'" << transition.transition_name << "'," << transition.product << "," << transition_charge << ",'" << transition.fragment_type<< "'," << transition.fragment_nr << "," << transition.detecting_transition << "," << transition.identifying_transition << "," << transition.quantifying_transition << "," << transition.library_intensity << "," << transition.decoy << "); ";
    }

    std::stringstream insert_precursor_sql, insert_precursor_peptide_mapping, insert_precursor_compound_mapping;
    insert_precursor_sql.precision(11);
    std::vector<std::pair<int, int> > peptide_protein_map;

    // OpenSWATH: Prepare peptide precursor inserts
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      OpenMS::TargetedExperiment::Peptide peptide = targeted_exp.getPeptides()[i];
      std::string peptide_sequence = TargetedExperimentHelper::getAASequence(peptide).toUniModString();
      int group_set_index = group_map[peptide.id];
      int peptide_set_index = peptide_map[peptide_sequence];

      for (std::vector<String>::iterator it = peptide.protein_refs.begin(); it != peptide.protein_refs.end(); ++it)
      {
        peptide_protein_map.push_back(std::make_pair(peptide_set_index,protein_map[*it]));
      }

      insert_precursor_sql << "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, DECOY) VALUES (" << group_set_index << ",'" << peptide.id << "','" << peptide.getPeptideGroupLabel() << "'," << precursor_mz_map[group_set_index] << "," << peptide.getChargeState() << ",NULL," << peptide.getRetentionTime() << "," << precursor_decoy_map[group_set_index] << "); ";

      insert_precursor_peptide_mapping << "INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (" << group_set_index << "," << peptide_set_index << "); ";

    }

    // OpenSWATH: Prepare compound precursor inserts
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      OpenMS::TargetedExperiment::Compound compound = targeted_exp.getCompounds()[i];
      int group_set_index = group_map[compound.id];
      int compound_set_index = compound_map[compound.id];

      std::string compound_charge = "NULL"; // workaround for compounds with missing charge
      if (compound.hasCharge())
      {
        compound_charge = String(compound.getChargeState());
      }

      insert_precursor_sql << "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, DECOY) VALUES (" << group_set_index << ",'" << compound.id << "',NULL," << precursor_mz_map[group_set_index] << "," << compound_charge << ",NULL,NULL" << "," << precursor_decoy_map[group_set_index] << "); ";

      insert_precursor_compound_mapping << "INSERT INTO PRECURSOR_COMPOUND_MAPPING (PRECURSOR_ID, COMPOUND_ID) VALUES (" << group_set_index << "," << compound_set_index << "); ";

    }

    boost::erase(peptide_protein_map, boost::unique<boost::return_found_end>(boost::sort(peptide_protein_map)));
    // OpenSWATH: Prepare peptide-protein mapping inserts
    std::stringstream insert_peptide_protein_mapping;
    for (std::vector<std::pair<int, int> >::iterator it = peptide_protein_map.begin(); it != peptide_protein_map.end(); ++it)
    {
      insert_peptide_protein_mapping << "INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (" << it->first << "," << it->second << "); ";
    }

    // OpenSWATH: Prepare protein inserts
    std::stringstream insert_protein_sql;
    for (std::map<std::string, int >::iterator it = protein_map.begin(); it != protein_map.end(); ++it)
    {
      insert_protein_sql << "INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (" << it->second << ",'" << it->first << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare peptide inserts
    std::stringstream insert_peptide_sql;
    for (std::map<std::string, int >::iterator it = peptide_map.begin(); it != peptide_map.end(); ++it)
    {
      insert_peptide_sql << "INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (" << it->second << ",'" << AASequence::fromString(it->first).toUnmodifiedString() << "','" << it->first << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare compound inserts
    std::stringstream insert_compound_sql;
    for (std::map<std::string, int >::iterator it = compound_map.begin(); it != compound_map.end(); ++it)
    {
      OpenMS::TargetedExperiment::Compound compound = targeted_exp.getCompoundByRef(it->first);
      insert_compound_sql << "INSERT INTO COMPOUND (ID, COMPOUND_NAME, SUM_FORMULA, SMILES, DECOY) VALUES (" << it->second << ",'" << compound.id << "','" << compound.molecular_formula << "','" << compound.smiles_string << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare decoy updates
    std::stringstream update_decoys_sql;
    // Peptides
    update_decoys_sql << "UPDATE PEPTIDE SET DECOY = 1 WHERE ID IN (SELECT PEPTIDE.ID FROM PRECURSOR JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID =  PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID WHERE PRECURSOR.DECOY = 1); ";
    // Compounds
    update_decoys_sql << "UPDATE COMPOUND SET DECOY = 1 WHERE ID IN (SELECT COMPOUND.ID FROM PRECURSOR JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID =  PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID WHERE PRECURSOR.DECOY = 1); ";
    // Proteins
    update_decoys_sql << "UPDATE PROTEIN SET DECOY = 1 WHERE ID IN (SELECT PROTEIN.ID FROM PEPTIDE JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE.ID =  PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID JOIN PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID WHERE PEPTIDE.DECOY = 1); ";

    sqlite3_exec(db, "BEGIN TRANSACTION", nullptr, nullptr, &zErrMsg);

    // Execute SQL insert statement
    std::string insert_protein_sql_str = insert_protein_sql.str();
    rc = sqlite3_exec(db, insert_protein_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_peptide_protein_mapping_str = insert_peptide_protein_mapping.str();
    rc = sqlite3_exec(db, insert_peptide_protein_mapping_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_peptide_sql_str = insert_peptide_sql.str();
    rc = sqlite3_exec(db, insert_peptide_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_compound_sql_str = insert_compound_sql.str();
    rc = sqlite3_exec(db, insert_compound_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_precursor_peptide_mapping_str = insert_precursor_peptide_mapping.str();
    rc = sqlite3_exec(db, insert_precursor_peptide_mapping_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_precursor_compound_mapping_str = insert_precursor_compound_mapping.str();
    rc = sqlite3_exec(db, insert_precursor_compound_mapping_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_precursor_sql_str = insert_precursor_sql.str();
    rc = sqlite3_exec(db, insert_precursor_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_transition_sql_str = insert_transition_sql.str();
    rc = sqlite3_exec(db, insert_transition_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_transition_peptide_mapping_sql_str = insert_transition_peptide_mapping_sql.str();
    rc = sqlite3_exec(db, insert_transition_peptide_mapping_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL insert statement
    std::string insert_transition_precursor_mapping_sql_str = insert_transition_precursor_mapping_sql.str();
    rc = sqlite3_exec(db, insert_transition_precursor_mapping_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    // Execute SQL update statement
    std::string update_decoys_sql_str = update_decoys_sql.str();
    rc = sqlite3_exec(db, update_decoys_sql_str.c_str(), callback, nullptr, &zErrMsg);
    if ( rc != SQLITE_OK )
    {
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          zErrMsg);
    }

    sqlite3_exec(db, "END TRANSACTION", nullptr, nullptr, &zErrMsg);

    sqlite3_close(db);

  }

  // public methods
  void TransitionPQPFile::convertTargetedExperimentToPQP(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    if (targeted_exp.containsInvalidReferences())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Your input file contains invalid references, cannot process file.");
    }
    writePQPOutput_(filename, targeted_exp);
  }

  void TransitionPQPFile::convertPQPToTargetedExperiment(const char* filename, OpenMS::TargetedExperiment& targeted_exp, bool legacy_traml_id)
  {
    std::vector<TSVTransition> transition_list;
    readPQPInput_(filename, transition_list, legacy_traml_id);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionPQPFile::convertPQPToTargetedExperiment(const char* filename, OpenSwath::LightTargetedExperiment& targeted_exp, bool legacy_traml_id)
  {
    std::vector<TSVTransition> transition_list;
    readPQPInput_(filename, transition_list, legacy_traml_id);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

}

