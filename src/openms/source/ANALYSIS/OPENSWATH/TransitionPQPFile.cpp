// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

#include <sqlite3.h>
#include <OpenMS/FORMAT/SqliteConnector.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include <sstream>
#include <unordered_map>
#include <iostream>

namespace OpenMS
{

  namespace Sql = Internal::SqliteHelper;

  TransitionPQPFile::TransitionPQPFile() :
    TransitionTSVFile()
  {
  }

  TransitionPQPFile::~TransitionPQPFile() = default;

  void TransitionPQPFile::readPQPInput_(const char* filename, std::vector<TSVTransition>& transition_list, bool legacy_traml_id)
  {
    sqlite3 *db;
    sqlite3_stmt * cntstmt;
    sqlite3_stmt * stmt;
    std::string select_sql;

    // Use legacy TraML identifiers for precursors (transition_group_id) and transitions (transition_name)?
    std::string traml_id = "ID";
    if (legacy_traml_id)
    {
      traml_id = "TRAML_ID";
    }

    startProgress(0, 1, "reading PQP file (SQL warmup)");

    // Open database
    SqliteConnector conn(filename);
    db = conn.getDB();

    // Count transitions
    SqliteConnector::prepareStatement(db, &cntstmt, "SELECT COUNT(*) FROM TRANSITION;");
    sqlite3_step( cntstmt );
    int num_transitions = sqlite3_column_int(cntstmt, 0);
    sqlite3_finalize(cntstmt);

    String select_drift_time = "";
    bool drift_time_exists = SqliteConnector::columnExists(db, "PRECURSOR", "LIBRARY_DRIFT_TIME");
    if (drift_time_exists)
    {
      select_drift_time = ", PRECURSOR.LIBRARY_DRIFT_TIME AS drift_time ";
    }

    String select_gene = "";
    String select_gene_null = "";
    String join_gene = "";
    bool gene_exists = SqliteConnector::tableExists(db, "GENE");
    if (gene_exists)
    {
      select_gene = ", GENE_AGGREGATED.GENE_NAME AS gene_name ";
      select_gene_null = ", 'NA' AS gene_name ";
      join_gene = "INNER JOIN PEPTIDE_GENE_MAPPING ON PEPTIDE.ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID " \
                  "INNER JOIN " \
                  "(SELECT PEPTIDE_ID, GROUP_CONCAT(GENE_NAME,';') AS GENE_NAME " \
                  "FROM GENE " \
                  "INNER JOIN PEPTIDE_GENE_MAPPING ON GENE.ID = PEPTIDE_GENE_MAPPING.GENE_ID "\
                  "GROUP BY PEPTIDE_ID) " \
                  "AS GENE_AGGREGATED ON PEPTIDE.ID = GENE_AGGREGATED.PEPTIDE_ID ";
    }

    String select_annotation = "'' AS Annotation, ";
    bool annotation_exists = SqliteConnector::columnExists(db, "TRANSITION", "ANNOTATION");
    if (annotation_exists) select_annotation = "TRANSITION.ANNOTATION AS Annotation, ";

    String select_adducts = "'' AS Adducts, ";
    bool adducts_exists = SqliteConnector::columnExists(db, "COMPOUND", "ADDUCTS");
    if (adducts_exists) select_adducts = "COMPOUND.ADDUCTS AS Adducts, ";

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
                  + select_annotation + \
                  "PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName, " \
                  "NULL AS CompoundName, " \
                  "NULL AS SMILES, " \
                  "NULL AS SumFormula, " \
                  "NULL AS Adducts, " \
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
                  "PEPTIDE_AGGREGATED.PEPTIDOFORMS AS peptidoforms " + \
                  select_drift_time + \
                  select_gene + \
                  "FROM PRECURSOR " + \
                  join_gene + \
                  "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID " \
                  "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID " \
                  "INNER JOIN " \
                    "(SELECT PEPTIDE_ID, GROUP_CONCAT(PROTEIN_ACCESSION,';') AS PROTEIN_ACCESSION " \
                    "FROM PROTEIN " \
                    "INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID "\
                    "GROUP BY PEPTIDE_ID) " \
                    "AS PROTEIN_AGGREGATED ON PEPTIDE.ID = PROTEIN_AGGREGATED.PEPTIDE_ID " \
                  "LEFT OUTER JOIN " \
                    "(SELECT TRANSITION_ID, GROUP_CONCAT(MODIFIED_SEQUENCE,'|') AS PEPTIDOFORMS " \
                    "FROM TRANSITION_PEPTIDE_MAPPING "\
                    "INNER JOIN PEPTIDE ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID "\
                    "GROUP BY TRANSITION_ID) "\
                    "AS PEPTIDE_AGGREGATED ON TRANSITION.ID = PEPTIDE_AGGREGATED.TRANSITION_ID ";

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
                  + select_annotation + \
                  "NULL AS FullPeptideName, " \
                  "COMPOUND.COMPOUND_NAME AS CompoundName, " \
                  "COMPOUND.SMILES AS SMILES, " \
                  "COMPOUND.SUM_FORMULA AS SumFormula, " \
                  + select_adducts + \
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
                  "NULL AS peptidoforms " +
                  select_drift_time +
                  select_gene_null +
                  "FROM PRECURSOR " \
                  "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON PRECURSOR.ID = TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN TRANSITION ON TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID = TRANSITION.ID " \
                  "INNER JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID " \
                  "INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID; ";


    // Execute SQL select statement
    SqliteConnector::prepareStatement(db, &stmt, select_sql);
    sqlite3_step(stmt);
    endProgress();

    Size progress = 0;
    startProgress(0, num_transitions, "reading PQP file");
    // Convert SQLite data to TSVTransition data structure
    while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
    {
      setProgress(progress++);
      TSVTransition mytransition;

      Sql::extractValue<double>(&mytransition.precursor, stmt, 0);
      Sql::extractValue<double>(&mytransition.product, stmt, 1);
      Sql::extractValue<double>(&mytransition.rt_calibrated, stmt, 2);
      Sql::extractValue<std::string>(&mytransition.transition_name, stmt, 3);
      Sql::extractValue<double>(&mytransition.CE, stmt, 4);
      Sql::extractValue<double>(&mytransition.library_intensity, stmt, 5);
      Sql::extractValue<std::string>(&mytransition.group_id, stmt, 6);
      Sql::extractValue<int>((int*)&mytransition.decoy, stmt, 7);
      Sql::extractValue<std::string>(&mytransition.PeptideSequence, stmt, 8);
      String tmp_field;
      if (Sql::extractValue<std::string>(&tmp_field, stmt, 9)) tmp_field.split(';', mytransition.ProteinName);
      Sql::extractValue<std::string>(&mytransition.Annotation, stmt, 10);
      Sql::extractValue<std::string>(&mytransition.FullPeptideName, stmt, 11);
      Sql::extractValue<std::string>(&mytransition.CompoundName, stmt, 12);
      Sql::extractValue<std::string>(&mytransition.SMILES, stmt, 13);
      Sql::extractValue<std::string>(&mytransition.SumFormula, stmt, 14);
      Sql::extractValue<std::string>(&mytransition.Adducts, stmt, 15);
      Sql::extractValueIntStr(&mytransition.precursor_charge, stmt, 16);
      Sql::extractValue<std::string>(&mytransition.peptide_group_label, stmt, 17);
      Sql::extractValue<std::string>(&mytransition.label_type, stmt, 18);
      Sql::extractValueIntStr(&mytransition.fragment_charge, stmt, 19);
      Sql::extractValue<int>(&mytransition.fragment_nr, stmt, 20);
      Sql::extractValue<double>(&mytransition.fragment_mzdelta, stmt, 21);
      Sql::extractValue<int>(&mytransition.fragment_modification, stmt, 22);
      Sql::extractValue<std::string>(&mytransition.fragment_type, stmt, 23);
      if (Sql::extractValue<std::string>(&tmp_field, stmt, 24)) tmp_field.split(';', mytransition.uniprot_id);
      Sql::extractValue<int>((int*)&mytransition.detecting_transition, stmt, 25);
      Sql::extractValue<int>((int*)&mytransition.identifying_transition, stmt, 26);
      Sql::extractValue<int>((int*)&mytransition.quantifying_transition, stmt, 27);
      if (Sql::extractValue<std::string>(&tmp_field, stmt, 28)) tmp_field.split('|', mytransition.peptidoforms);
      // optional attributes only present in newer file versions
      if (drift_time_exists) Sql::extractValue<double>(&mytransition.drift_time, stmt, 29);
      if (gene_exists) Sql::extractValue<std::string>(&mytransition.GeneName, stmt, 30);

      if (mytransition.GeneName == "NA") mytransition.GeneName = "";

      transition_list.push_back(mytransition);
      sqlite3_step( stmt );
    }
    endProgress();

    sqlite3_finalize(stmt);
  }

  void TransitionPQPFile::writePQPOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    // delete file if present
    remove(filename);

    // Open database
    SqliteConnector conn(filename);

    // Create SQL structure
    const char* create_sql =
      "CREATE TABLE VERSION(" \
      "ID INT NOT NULL);" \

      // gene table
      // OpenSWATH proteomics workflows
      "CREATE TABLE GENE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "GENE_NAME TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // peptide_gene_mapping table
      // OpenSWATH proteomics workflows
      "CREATE TABLE PEPTIDE_GENE_MAPPING(" \
      "PEPTIDE_ID INT NOT NULL," \
      "GENE_ID INT NOT NULL);" \

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
      "ADDUCTS TEXT NOT NULL," \
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
      "LIBRARY_DRIFT_TIME REAL NULL," \
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
      "ANNOTATION TEXT NULL," \
      "ORDINAL INT NULL," \
      "DETECTING INT NOT NULL," \
      "IDENTIFYING INT NOT NULL," \
      "QUANTIFYING INT NOT NULL," \
      "LIBRARY_INTENSITY REAL NULL," \
      "DECOY INT NOT NULL);";

    // Execute SQL create statement
    conn.executeStatement(create_sql);

    // Prepare insert statements

    // Index maps
    std::vector<std::string> group_vec, peptide_vec, compound_vec, protein_vec;
    std::unordered_map<std::string, int > group_map, peptide_map, compound_map, protein_map, gene_map;
    std::unordered_map<int,double> precursor_mz_map;
    std::unordered_map<int,bool> precursor_decoy_map;


    // OpenSWATH: Loop through TargetedExperiment to generate index maps for peptides
    peptide_vec.reserve(targeted_exp.getPeptides().size());
    group_vec.reserve(targeted_exp.getPeptides().size()  + targeted_exp.getCompounds().size());
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      OpenMS::TargetedExperiment::Peptide peptide = targeted_exp.getPeptides()[i];
      std::string peptide_sequence = TargetedExperimentHelper::getAASequence(peptide).toUniModString();
      peptide_vec.push_back(peptide_sequence);
      group_vec.push_back(peptide.id);
    }

    // OpenSWATH: Loop through TargetedExperiment to generate index maps for compounds
    compound_vec.reserve(targeted_exp.getCompounds().size()); 
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      OpenMS::TargetedExperiment::Compound compound = targeted_exp.getCompounds()[i];
      compound_vec.push_back(compound.id);
      group_vec.push_back(compound.id);
    }

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

    // OpenSWATH: Group set must be unique
    boost::erase(group_vec, boost::unique<boost::return_found_end>(boost::sort(group_vec)));
    int group_map_idx = 0;
    for (auto const & x : group_vec) { group_map[x] = group_map_idx; group_map_idx++; }

    // IPF: Loop through all transitions and generate peptidoform data structures
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      std::vector<String> peptidoforms;
      String(targeted_exp.getTransitions()[i].getMetaValue("Peptidoforms")).split('|', peptidoforms);
      std::copy( peptidoforms.begin(), peptidoforms.end(),
          std::inserter( peptide_vec, peptide_vec.end() ) );
    }
    // OpenSWATH: Peptide and compound sets must be unique
    boost::erase(peptide_vec, boost::unique<boost::return_found_end>(boost::sort(peptide_vec)));
    int peptide_map_idx = 0;
    for (auto const & x : peptide_vec) { peptide_map[x] = peptide_map_idx; peptide_map_idx++; }

    {
      std::stringstream insert_transition_sql, insert_transition_peptide_mapping_sql, insert_transition_precursor_mapping_sql;
      insert_transition_sql.precision(11);
      for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
      {
        TransitionPQPFile::TSVTransition transition = convertTransition_(&targeted_exp.getTransitions()[i], targeted_exp);

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

        // IPF: Generate transition-peptide mapping tables (one identification transition can map to multiple peptidoforms)
        for (Size j = 0; j < transition.peptidoforms.size(); j++)
        {
          insert_transition_peptide_mapping_sql << "INSERT INTO TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID, PEPTIDE_ID) VALUES (" <<
            i << "," << peptide_map[transition.peptidoforms[j]] << "); ";
        }

        // OpenSWATH: Associate transitions with their precursors
        insert_transition_precursor_mapping_sql << "INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (" <<
          i << "," << group_map[transition.group_id] << "); ";

        std::string transition_charge = "NULL"; // workaround for compounds with missing charge
        if (transition.fragment_charge != "NA")
        {
          transition_charge = transition.fragment_charge;
        }

        // OpenSWATH: Insert transition data
        insert_transition_sql << "INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, " <<
          "DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (" << i << ",'" <<
          transition.transition_name << "'," <<
          transition.product << "," <<
          transition_charge << ",'" <<
          transition.fragment_type << "','" <<
          transition.Annotation <<"'," <<
          transition.fragment_nr << "," <<
          transition.detecting_transition << "," <<
          transition.identifying_transition << "," <<
          transition.quantifying_transition << "," <<
          transition.library_intensity << "," << transition.decoy << "); ";

        if (i % 50000 == 0)
        // if (i % 2 == 0) // for testing
        {
          conn.executeStatement("BEGIN TRANSACTION");
          conn.executeStatement(insert_transition_sql.str());
          conn.executeStatement(insert_transition_peptide_mapping_sql.str());
          conn.executeStatement(insert_transition_precursor_mapping_sql.str());
          conn.executeStatement("END TRANSACTION");
          insert_transition_sql.str("");
          insert_transition_sql.clear();
          insert_transition_peptide_mapping_sql.str("");
          insert_transition_peptide_mapping_sql.clear();
          insert_transition_precursor_mapping_sql.str("");
          insert_transition_precursor_mapping_sql.clear();
        }
      }
      conn.executeStatement("BEGIN TRANSACTION");
      conn.executeStatement(insert_transition_sql.str());
      conn.executeStatement(insert_transition_peptide_mapping_sql.str());
      conn.executeStatement(insert_transition_precursor_mapping_sql.str());
      conn.executeStatement("END TRANSACTION");
    }

    std::stringstream insert_precursor_sql, insert_precursor_peptide_mapping, insert_precursor_compound_mapping;
    insert_precursor_sql.precision(11);
    std::vector<std::pair<int, int> > peptide_protein_map;
    std::vector<std::pair<int, int> > peptide_gene_map;

    // OpenSWATH: Prepare peptide precursor inserts
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      OpenMS::TargetedExperiment::Peptide peptide = targeted_exp.getPeptides()[i];
      std::string peptide_sequence = TargetedExperimentHelper::getAASequence(peptide).toUniModString();
      int group_set_index = group_map[peptide.id];
      int peptide_set_index = peptide_map[peptide_sequence];

      for (const auto& it : peptide.protein_refs)
      {
        peptide_protein_map.emplace_back(peptide_set_index, protein_map[it]);
      }

      String gene_name = "NA";
      if (peptide.metaValueExists("GeneName"))
      {
        gene_name = peptide.getMetaValue("GeneName");
      }

      if (gene_map.find(gene_name) == gene_map.end()) gene_map[gene_name] = (int)gene_map.size();
      peptide_gene_map.emplace_back(peptide_set_index, gene_map[gene_name]);

      insert_precursor_sql <<
        "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, " <<
        "LIBRARY_DRIFT_TIME, LIBRARY_RT, DECOY) VALUES (" <<
        group_set_index << ",'" << peptide.id << "','" <<
        peptide.getPeptideGroupLabel() << "'," <<
        precursor_mz_map[group_set_index] << "," <<
        peptide.getChargeState() <<
        ",NULL," <<
        peptide.getDriftTime() << "," <<
        peptide.getRetentionTime() << "," <<
        precursor_decoy_map[group_set_index] << "); ";

      insert_precursor_peptide_mapping << "INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (" <<
        group_set_index << "," << peptide_set_index << "); ";

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

      insert_precursor_sql << "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, " <<
        "LIBRARY_DRIFT_TIME, LIBRARY_RT, DECOY) VALUES (" << group_set_index
        << ",'" << compound.id << "',NULL," <<
        precursor_mz_map[group_set_index] << "," <<
        compound_charge <<
        ",NULL," <<
        compound.getDriftTime() << "," <<
        compound.getRetentionTime() << "," <<
        precursor_decoy_map[group_set_index] << "); ";

      insert_precursor_compound_mapping << "INSERT INTO PRECURSOR_COMPOUND_MAPPING (PRECURSOR_ID, COMPOUND_ID) VALUES (" <<
        group_set_index << "," << compound_set_index << "); ";
    }

    boost::erase(peptide_protein_map, boost::unique<boost::return_found_end>(boost::sort(peptide_protein_map)));
    boost::erase(peptide_gene_map, boost::unique<boost::return_found_end>(boost::sort(peptide_gene_map)));

    // OpenSWATH: Prepare peptide-gene mapping inserts
    std::stringstream insert_peptide_gene_mapping;
    for (const auto& it : peptide_gene_map)
    {
      insert_peptide_gene_mapping << "INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (" <<
        it.first << "," << it.second << "); ";
    }
    // OpenSWATH: Prepare gene inserts
    std::stringstream insert_gene_sql;
    for (const auto& it : gene_map)
    {
      insert_gene_sql << "INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (" <<
        it.second << ",'" << it.first << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare peptide-protein mapping inserts
    std::stringstream insert_peptide_protein_mapping;
    for (const auto& it : peptide_protein_map)
    {
      insert_peptide_protein_mapping << "INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (" <<
        it.first << "," << it.second << "); ";
    }

    // OpenSWATH: Prepare protein inserts
    std::stringstream insert_protein_sql;
    for (const auto& it : protein_map)
    {
      insert_protein_sql << "INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (" <<
        it.second << ",'" << it.first << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare peptide inserts
    std::stringstream insert_peptide_sql;
    for (const auto& it : peptide_map)
    {
      insert_peptide_sql << "INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (" <<
        it.second << ",'" <<
        AASequence::fromString(it.first).toUnmodifiedString() << "','" <<
        it.first << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare compound inserts
    std::stringstream insert_compound_sql;
    for (const auto& it : compound_map)
    {
      String adducts;
      String compound_name;
      const auto& compound = targeted_exp.getCompoundByRef(it.first);
      if (compound.metaValueExists("Adducts"))
      {
        adducts = compound.getMetaValue("Adducts");
      }
      if (compound.metaValueExists("CompoundName"))
      {
        compound_name = compound.getMetaValue("CompoundName");
      }
      else
      {
        compound_name = compound.id;
      }
      insert_compound_sql << "INSERT INTO COMPOUND (ID, COMPOUND_NAME, SUM_FORMULA, SMILES, ADDUCTS, DECOY) VALUES (" <<
        it.second << ",'" <<
        compound_name << "','" <<
        compound.molecular_formula << "','" <<
        compound.smiles_string << "','" <<
        adducts << "'," <<
        0 << "); ";
    }

    // OpenSWATH: Prepare decoy updates
    std::stringstream update_decoys_sql;
    // Peptides
    update_decoys_sql << "UPDATE PEPTIDE SET DECOY = 1 WHERE ID IN " <<
      "(SELECT PEPTIDE.ID FROM PRECURSOR " <<
      "JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID " <<
      "JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID WHERE PRECURSOR.DECOY = 1); ";
    // Compounds
    update_decoys_sql << "UPDATE COMPOUND SET DECOY = 1 WHERE ID IN " <<
      "(SELECT COMPOUND.ID FROM PRECURSOR " <<
      "JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID " <<
      "JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID WHERE PRECURSOR.DECOY = 1); ";
    // Proteins
    update_decoys_sql << "UPDATE PROTEIN SET DECOY = 1 WHERE ID IN " <<
      "(SELECT PROTEIN.ID FROM PEPTIDE " <<
      "JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE.ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID " <<
      "JOIN PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID WHERE PEPTIDE.DECOY = 1); ";
    // Genes
    update_decoys_sql << "UPDATE GENE SET DECOY = 1 WHERE ID IN " <<
      "(SELECT GENE.ID FROM PEPTIDE " <<
      "JOIN PEPTIDE_GENE_MAPPING ON PEPTIDE.ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID " <<
      "JOIN GENE ON PEPTIDE_GENE_MAPPING.GENE_ID = GENE.ID WHERE PEPTIDE.DECOY = 1); ";

    conn.executeStatement("BEGIN TRANSACTION");

    // Execute SQL insert statement
    String insert_version = "INSERT INTO VERSION (ID) VALUES (3);";
    conn.executeStatement(insert_version);
    conn.executeStatement(insert_protein_sql.str());
    conn.executeStatement(insert_peptide_protein_mapping.str());
    conn.executeStatement(insert_gene_sql.str());
    conn.executeStatement(insert_peptide_gene_mapping.str());
    conn.executeStatement(insert_peptide_sql.str());
    conn.executeStatement(insert_compound_sql.str());
    conn.executeStatement(insert_precursor_peptide_mapping.str());
    conn.executeStatement(insert_precursor_compound_mapping.str());
    conn.executeStatement(insert_precursor_sql.str());
    conn.executeStatement(update_decoys_sql.str());
    conn.executeStatement("END TRANSACTION");
  }

  // public methods
  std::unordered_map<std::string, std::string> TransitionPQPFile::getPQPIDToTraMLIDMap(const char* filename, std::string tableName)
  {
    sqlite3 *db;
    sqlite3_stmt * cntstmt;
    sqlite3_stmt * stmt;
    std::string select_sql;
    std::unordered_map<std::string, std::string> out;

    // Open database
    SqliteConnector conn(filename);
    db = conn.getDB();

    // Count Precursors 
    SqliteConnector::prepareStatement(db, &cntstmt, "SELECT COUNT(*) FROM " + tableName + ";");
    sqlite3_step( cntstmt );
    sqlite3_finalize(cntstmt);

    std::string query = "SELECT ID, TRAML_ID FROM " + tableName + ";"; 

    // Execute SQL select statement
    SqliteConnector::prepareStatement(db, &stmt, query);
    sqlite3_step(stmt);

    while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
    {
      std::string traml_id, prec_id;

      Sql::extractValue<std::string>(&prec_id, stmt, 0);
      Sql::extractValue<std::string>(&traml_id, stmt, 1);

      out[traml_id] = prec_id;
      sqlite3_step( stmt );
    }
    return out;
  }

  void TransitionPQPFile::convertTargetedExperimentToPQP(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    if (targeted_exp.containsInvalidReferences())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Your input file contains invalid references, cannot process file.");
    }
    writePQPOutput_(filename, targeted_exp);
  }

  void TransitionPQPFile::convertPQPToTargetedExperiment(const char* filename,
                                                         OpenMS::TargetedExperiment& targeted_exp,
                                                         bool legacy_traml_id)
  {
    std::vector<TSVTransition> transition_list;
    readPQPInput_(filename, transition_list, legacy_traml_id);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionPQPFile::convertPQPToTargetedExperiment(const char* filename,
                                                         OpenSwath::LightTargetedExperiment& targeted_exp,
                                                         bool legacy_traml_id)
  {
    std::vector<TSVTransition> transition_list;
    readPQPInput_(filename, transition_list, legacy_traml_id);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

}
