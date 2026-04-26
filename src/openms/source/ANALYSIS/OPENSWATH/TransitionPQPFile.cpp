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
#include <OpenMS/CONCEPT/LogStream.h>

#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{

  namespace Sql = Internal::SqliteHelper;

  namespace
  {
    constexpr Size default_transitions_per_precursor = 6;
    constexpr Size max_protein_reserve = 100000;
    constexpr Size stable_stream_order_transition_limit = 100000;

    Size estimatePrecursorCountFromTransitions_(Size num_transitions)
    {
      return num_transitions / default_transitions_per_precursor + 1;
    }

    bool tableHasRows_(sqlite3* db, const String& table_name)
    {
      sqlite3_stmt* stmt = nullptr;
      SqliteConnector::prepareStatement(db, &stmt, "SELECT 1 FROM " + table_name + " LIMIT 1;");
      const bool has_rows = sqlite3_step(stmt) == SQLITE_ROW;
      sqlite3_finalize(stmt);
      return has_rows;
    }

    void configurePQPReadConnection_(sqlite3* db)
    {
      // Large PQP libraries are streamed once from a read-only connection. A bounded
      // SQLite page cache and mmap window reduce cold-cache wall time without
      // materializing the database or changing file contents.
      SqliteConnector::executeStatement(db, "PRAGMA cache_size = -524288;");
      SqliteConnector::executeStatement(db, "PRAGMA mmap_size = 8589934592;");
      SqliteConnector::executeStatement(db, "PRAGMA query_only = ON;");
#ifdef _OPENMP
      int sqlite_threads = omp_get_max_threads();
      if (sqlite_threads < 1)
      {
        sqlite_threads = 1;
      }
      SqliteConnector::executeStatement(db, String("PRAGMA threads = ") + String(sqlite_threads) + ";");
#endif
    }
  }

  TransitionPQPFile::TransitionPQPFile() :
    TransitionTSVFile()
  {
  }

  TransitionPQPFile::~TransitionPQPFile() = default;

  TransitionPQPFile::PQPSqlQueryInfo TransitionPQPFile::buildPQPSelectQuery_(sqlite3* db, bool legacy_traml_id, bool stable_order) const
  {
    PQPSqlQueryInfo info;

    // Use legacy TraML identifiers for precursors (transition_group_id) and transitions (transition_name)?
    // When legacy_traml_id=true, use TRAML_ID column (string identifiers from TraML)
    // When legacy_traml_id=false, use numeric ID column
    std::string traml_id = legacy_traml_id ? "TRAML_ID" : "ID";

    // Check for optional columns/tables
    String select_drift_time = "";
    info.drift_time_exists = SqliteConnector::columnExists(db, "PRECURSOR", "LIBRARY_DRIFT_TIME");
    if (info.drift_time_exists)
    {
      select_drift_time = ", PRECURSOR.LIBRARY_DRIFT_TIME AS drift_time ";
    }

    String select_gene = "";
    String select_gene_null = "";
    String join_gene = "";
    info.gene_exists = SqliteConnector::tableExists(db, "GENE");
    if (info.gene_exists)
    {
      // Check if the GENE table actually has any data (issue #8687).
      // An empty GENE table with INNER JOIN would return zero rows.
      info.gene_exists = tableHasRows_(db, "GENE");
    }
    if (info.gene_exists)
    {
      select_gene = ", GENE_AGGREGATED.GENE_NAME AS gene_name ";
      select_gene_null = ", 'NA' AS gene_name ";
      // Join only the per-peptide aggregate. A raw PEPTIDE_GENE_MAPPING join
      // duplicates transition rows for multi-gene peptides and later breaks
      // transition-group construction when native IDs collide.
      join_gene = "INNER JOIN " \
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

    info.compound_exists = SqliteConnector::tableExists(db, "PRECURSOR_COMPOUND_MAPPING") &&
                           tableHasRows_(db, "PRECURSOR_COMPOUND_MAPPING");

    info.peptidoforms_exists = SqliteConnector::tableExists(db, "TRANSITION_PEPTIDE_MAPPING") &&
                               tableHasRows_(db, "TRANSITION_PEPTIDE_MAPPING");

    String select_peptidoforms = "NULL AS peptidoforms ";
    String join_peptidoforms = "";
    if (info.peptidoforms_exists)
    {
      select_peptidoforms = "PEPTIDE_AGGREGATED.PEPTIDOFORMS AS peptidoforms ";
      join_peptidoforms = "LEFT OUTER JOIN " \
                    "(SELECT TRANSITION_ID, GROUP_CONCAT(MODIFIED_SEQUENCE,'|') AS PEPTIDOFORMS " \
                    "FROM TRANSITION_PEPTIDE_MAPPING "\
                    "INNER JOIN PEPTIDE ON TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID "\
                    "GROUP BY TRANSITION_ID) "\
                    "AS PEPTIDE_AGGREGATED ON TRANSITION.ID = PEPTIDE_AGGREGATED.TRANSITION_ID ";
    }

    // Build peptides query
    info.select_sql = "SELECT " \
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
                  + select_peptidoforms + \
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
                    "AS PROTEIN_AGGREGATED ON PEPTIDE.ID = PROTEIN_AGGREGATED.PEPTIDE_ID " +
                  join_peptidoforms;

    // Append compounds query
    if (info.compound_exists)
    {
      info.select_sql += "UNION ALL SELECT " \
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
                  "INNER JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID ";
    }

    if (stable_order)
    {
      // TargetedFileConverter and small TOPP round-trip tests historically saw the
      // implicit UNION ordering. Keep that deterministic order explicit while the
      // large-library streaming path can still avoid a global sort.
      info.select_sql += "ORDER BY precursor, product, rt_calibrated, transition_name, group_id ";
    }
    info.select_sql += "; ";

    return info;
  }

  void TransitionPQPFile::readPQPInput_(const char* filename, std::vector<TSVTransition>& transition_list, bool legacy_traml_id)
  {
    sqlite3 *db;
    sqlite3_stmt * cntstmt;
    sqlite3_stmt * stmt;

    startProgress(0, 1, "reading PQP file (SQL warmup)");

    // Open database
    SqliteConnector conn(filename, SqliteConnector::SqlOpenMode::READ_ONLY);
    db = conn.getDB();
    configurePQPReadConnection_(db);

    // Count transitions
    SqliteConnector::prepareStatement(db, &cntstmt, "SELECT COUNT(*) FROM TRANSITION;");
    sqlite3_step( cntstmt );
    const Size num_transitions = static_cast<Size>(sqlite3_column_int64(cntstmt, 0));
    sqlite3_finalize(cntstmt);
    transition_list.reserve(num_transitions);

    // Build SQL query using shared helper
    PQPSqlQueryInfo query_info = buildPQPSelectQuery_(db, legacy_traml_id, true);

    // Execute SQL select statement
    SqliteConnector::prepareStatement(db, &stmt, query_info.select_sql);
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
      int optional_col = 29;
      if (query_info.drift_time_exists)
      {
        Sql::extractValue<double>(&mytransition.drift_time, stmt, optional_col++);
      }
      if (query_info.gene_exists)
      {
        Sql::extractValue<std::string>(&mytransition.GeneName, stmt, optional_col++);
      }

      if (mytransition.GeneName == "NA") mytransition.GeneName = "";

      transition_list.push_back(mytransition);
      sqlite3_step( stmt );
    }
    endProgress();

    sqlite3_finalize(stmt);
  }

  void TransitionPQPFile::streamPQPToLightTargetedExperiment_(const char* filename, OpenSwath::LightTargetedExperiment& exp, bool legacy_traml_id)
  {
    // Sets for deduplication
    std::unordered_set<std::string> compound_seen;
    std::unordered_set<std::string> protein_seen;

    sqlite3 *db;
    sqlite3_stmt * cntstmt;
    sqlite3_stmt * stmt;

    startProgress(0, 1, "reading PQP file (SQL warmup)");

    // Open database
    SqliteConnector conn(filename, SqliteConnector::SqlOpenMode::READ_ONLY);
    db = conn.getDB();
    configurePQPReadConnection_(db);

    // Count transitions
    SqliteConnector::prepareStatement(db, &cntstmt, "SELECT COUNT(*) FROM TRANSITION;");
    sqlite3_step( cntstmt );
    const Size num_transitions = static_cast<Size>(sqlite3_column_int64(cntstmt, 0));
    sqlite3_finalize(cntstmt);
    exp.transitions.reserve(num_transitions);

    // Avoid full PRECURSOR/PROTEIN scans just for reserve sizes. The main
    // streaming query will read those tables once; reserve from the transition
    // count keeps allocation churn low without extra cold-cache I/O.
    const Size estimated_precursor_count = estimatePrecursorCountFromTransitions_(num_transitions);
    if (SqliteConnector::tableExists(db, "PRECURSOR"))
    {
      exp.compounds.reserve(estimated_precursor_count);
      compound_seen.reserve(estimated_precursor_count);
    }
    if (SqliteConnector::tableExists(db, "PROTEIN"))
    {
      const Size estimated_protein_count = std::min(estimated_precursor_count, max_protein_reserve);
      exp.proteins.reserve(estimated_protein_count);
      protein_seen.reserve(estimated_protein_count);
    }

    // Build SQL query using shared helper
    const bool stable_order = num_transitions <= stable_stream_order_transition_limit;
    PQPSqlQueryInfo query_info = buildPQPSelectQuery_(db, legacy_traml_id, stable_order);

    // Execute SQL select statement
    SqliteConnector::prepareStatement(db, &stmt, query_info.select_sql);
    sqlite3_step(stmt);
    endProgress();

    Size progress = 0;
    startProgress(0, num_transitions, "streaming PQP to LightTargetedExperiment");

    // Stream SQL results directly to LightTargetedExperiment
    while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
    {
      setProgress(progress++);

      // Extract values directly into variables
      double precursor_mz = 0, product_mz = 0, rt_calibrated = 0, library_intensity = 0, drift_time = -1;
      std::string transition_name, group_id, peptide_sequence, full_peptide_name;
      std::string compound_name, smiles, sum_formula, adducts_str;
      std::string peptide_group_label, fragment_type_str, gene_name;
      int decoy = 0, precursor_charge = 0, fragment_charge = 0, fragment_nr = -1;
      int detecting = 1, identifying = 0, quantifying = 1;
      String protein_names_str, peptidoforms_str;

      Sql::extractValue<double>(&precursor_mz, stmt, 0);
      Sql::extractValue<double>(&product_mz, stmt, 1);
      Sql::extractValue<double>(&rt_calibrated, stmt, 2);
      Sql::extractValue<std::string>(&transition_name, stmt, 3);
      // Skip CE (column 4) - not used in light path
      Sql::extractValue<double>(&library_intensity, stmt, 5);
      Sql::extractValue<std::string>(&group_id, stmt, 6);
      Sql::extractValue<int>(&decoy, stmt, 7);
      Sql::extractValue<std::string>(&peptide_sequence, stmt, 8);
      Sql::extractValue<String>(&protein_names_str, stmt, 9);
      // Skip Annotation (column 10) - reconstructed from fragment info
      Sql::extractValue<std::string>(&full_peptide_name, stmt, 11);
      Sql::extractValue<std::string>(&compound_name, stmt, 12);
      Sql::extractValue<std::string>(&smiles, stmt, 13);
      Sql::extractValue<std::string>(&sum_formula, stmt, 14);
      Sql::extractValue<std::string>(&adducts_str, stmt, 15);
      Sql::extractValue<int>(&precursor_charge, stmt, 16);
      Sql::extractValue<std::string>(&peptide_group_label, stmt, 17);
      // Skip label_type (column 18) - not in PQP
      Sql::extractValue<int>(&fragment_charge, stmt, 19);
      Sql::extractValue<int>(&fragment_nr, stmt, 20);
      // Skip fragment_mzdelta (column 21)
      // Skip fragment_modification (column 22)
      Sql::extractValue<std::string>(&fragment_type_str, stmt, 23);
      // Skip uniprot_id (column 24) - not in PQP
      Sql::extractValue<int>(&detecting, stmt, 25);
      Sql::extractValue<int>(&identifying, stmt, 26);
      Sql::extractValue<int>(&quantifying, stmt, 27);
      Sql::extractValue<String>(&peptidoforms_str, stmt, 28);
      int optional_col = 29;
      if (query_info.drift_time_exists)
      {
        Sql::extractValue<double>(&drift_time, stmt, optional_col++);
      }
      if (query_info.gene_exists)
      {
        Sql::extractValue<std::string>(&gene_name, stmt, optional_col++);
        if (gene_name == "NA") gene_name = "";
      }

      // Create LightTransition directly
      OpenSwath::LightTransition transition;
      transition.transition_name = transition_name;
      transition.peptide_ref = group_id;
      transition.library_intensity = library_intensity;
      transition.precursor_mz = precursor_mz;
      transition.product_mz = product_mz;
      transition.precursor_im = drift_time;
      transition.fragment_charge = static_cast<int8_t>(fragment_charge);
      transition.setDecoy(decoy != 0);
      transition.setDetectingTransition(detecting != 0);
      transition.setIdentifyingTransition(identifying != 0);
      transition.setQuantifyingTransition(quantifying != 0);
      transition.fragment_nr = static_cast<int16_t>(fragment_nr);
      transition.setFragmentType(fragment_type_str);
      if (!peptidoforms_str.empty())
      {
        std::vector<String> peptidoforms_tmp;
        peptidoforms_str.split('|', peptidoforms_tmp);
        transition.peptidoforms.assign(peptidoforms_tmp.begin(), peptidoforms_tmp.end());
      }

      exp.transitions.push_back(std::move(transition));

      // Create compound if needed
      if (compound_seen.insert(group_id).second)
      {
        OpenSwath::LightCompound compound;
        compound.id = group_id;
        compound.drift_time = drift_time;
        compound.rt = rt_calibrated;
        compound.charge = precursor_charge;
        compound.peptide_group_label = peptide_group_label;
        compound.gene_name = gene_name;

        bool is_peptide = compound_name.empty();
        std::vector<String> protein_names;
        if (is_peptide)
        {
          compound.sequence = full_peptide_name.empty() ? peptide_sequence : full_peptide_name;
          if (!protein_names_str.empty())
          {
            protein_names_str.split(';', protein_names);
            compound.protein_refs.assign(protein_names.begin(), protein_names.end());
          }
          // Parse modifications from sequence
          String sequence = full_peptide_name.empty() ? peptide_sequence : full_peptide_name;
          if (!sequence.empty())
          {
            try
            {
              AASequence aa_sequence = AASequence::fromString(sequence);
              if (aa_sequence.hasNTerminalModification())
              {
                OpenSwath::LightModification mod;
                mod.location = -1;
                mod.unimod_id = aa_sequence.getNTerminalModification()->getUniModRecordId();
                compound.modifications.push_back(mod);
              }
              if (aa_sequence.hasCTerminalModification())
              {
                OpenSwath::LightModification mod;
                mod.location = static_cast<int>(aa_sequence.size());
                mod.unimod_id = aa_sequence.getCTerminalModification()->getUniModRecordId();
                compound.modifications.push_back(mod);
              }
              for (Size i = 0; i != aa_sequence.size(); i++)
              {
                if (aa_sequence[i].isModified())
                {
                  OpenSwath::LightModification mod;
                  mod.location = static_cast<int>(i);
                  mod.unimod_id = aa_sequence.getResidue(i).getModification()->getUniModRecordId();
                  compound.modifications.push_back(mod);
                }
              }
            }
            catch (Exception::InvalidValue&)
            {
              OPENMS_LOG_DEBUG << "Could not parse modifications from sequence: " << sequence << std::endl;
            }
          }
        }
        else
        {
          compound.compound_name = compound_name;
          compound.sum_formula = sum_formula;
          compound.smiles = smiles;
          compound.adducts = adducts_str;
        }

        exp.compounds.push_back(std::move(compound));

        // Create proteins once per precursor instead of re-splitting the same
        // protein list for every transition row.
        if (is_peptide)
        {
          for (const auto& pname : protein_names)
          {
            const std::string protein_id = pname;
            if (protein_seen.insert(protein_id).second)
            {
              OpenSwath::LightProtein protein;
              protein.id = protein_id;
              protein.sequence = "";
              exp.proteins.push_back(std::move(protein));
            }
          }
        }
      }

      sqlite3_step(stmt);
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
    std::vector<std::string> gene_vec;
    gene_vec.reserve(gene_map.size());
    for (const auto& it : gene_map)
    {
      gene_vec.push_back(it.first);
    }
    boost::sort(gene_vec);
    for (const auto& gene_name : gene_vec)
    {
      insert_gene_sql << "INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (" <<
        gene_map[gene_name] << ",'" << gene_name << "'," << 0 << "); ";
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
    for (const auto& protein_id : protein_vec)
    {
      insert_protein_sql << "INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (" <<
        protein_map[protein_id] << ",'" << protein_id << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare peptide inserts
    std::stringstream insert_peptide_sql;
    for (const auto& peptide_sequence : peptide_vec)
    {
      insert_peptide_sql << "INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (" <<
        peptide_map[peptide_sequence] << ",'" <<
        AASequence::fromString(peptide_sequence).toUnmodifiedString() << "','" <<
        peptide_sequence << "'," << 0 << "); ";
    }

    // OpenSWATH: Prepare compound inserts
    std::stringstream insert_compound_sql;
    for (const auto& compound_id : compound_vec)
    {
      String adducts;
      String compound_name;
      const auto& compound = targeted_exp.getCompoundByRef(compound_id);
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
        compound_map[compound_id] << ",'" <<
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
    // Use streaming parser for memory efficiency (~5x reduction in peak memory)
    streamPQPToLightTargetedExperiment_(filename, targeted_exp, legacy_traml_id);
  }

  void TransitionPQPFile::convertLightTargetedExperimentToPQP(const char* filename, const OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    // delete file if present
    remove(filename);

    // Open database
    SqliteConnector conn(filename);

    // Create SQL structure (same as heavy version)
    const char* create_sql =
      "CREATE TABLE VERSION(" \
      "ID INT NOT NULL);" \

      // gene table
      "CREATE TABLE GENE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "GENE_NAME TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // peptide_gene_mapping table
      "CREATE TABLE PEPTIDE_GENE_MAPPING(" \
      "PEPTIDE_ID INT NOT NULL," \
      "GENE_ID INT NOT NULL);" \

      // protein table
      "CREATE TABLE PROTEIN(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "PROTEIN_ACCESSION TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // peptide_protein_mapping table
      "CREATE TABLE PEPTIDE_PROTEIN_MAPPING(" \
      "PEPTIDE_ID INT NOT NULL," \
      "PROTEIN_ID INT NOT NULL);" \

      // peptide table
      "CREATE TABLE PEPTIDE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "UNMODIFIED_SEQUENCE TEXT NOT NULL," \
      "MODIFIED_SEQUENCE TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // precursor_peptide_mapping table
      "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(" \
      "PRECURSOR_ID INT NOT NULL," \
      "PEPTIDE_ID INT NOT NULL);" \

      // compound table
      "CREATE TABLE COMPOUND(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "COMPOUND_NAME TEXT NOT NULL," \
      "SUM_FORMULA TEXT NOT NULL," \
      "SMILES TEXT NOT NULL," \
      "ADDUCTS TEXT NOT NULL," \
      "DECOY INT NOT NULL);" \

      // precursor_compound_mapping table
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

    // Build index maps
    std::vector<std::string> group_vec, peptide_vec, compound_vec, protein_vec;
    std::unordered_map<std::string, int> group_map, peptide_map, compound_map, protein_map, gene_map;
    std::unordered_map<int, double> precursor_mz_map;
    std::unordered_map<int, bool> precursor_decoy_map;

    // Loop through compounds to generate index maps
    for (const auto& compound : targeted_exp.compounds)
    {
      group_vec.push_back(compound.id);
      if (compound.isPeptide())
      {
        peptide_vec.push_back(compound.sequence);
      }
      else
      {
        compound_vec.push_back(compound.id);
      }
    }

    // Loop through proteins
    for (const auto& protein : targeted_exp.proteins)
    {
      protein_vec.push_back(protein.id);
    }

    // Loop through transitions and add peptidoforms
    for (const auto& tr : targeted_exp.transitions)
    {
      for (const auto& peptidoform : tr.peptidoforms)
      {
        peptide_vec.push_back(peptidoform);
      }
    }

    // Create unique sorted sets and maps
    boost::erase(compound_vec, boost::unique<boost::return_found_end>(boost::sort(compound_vec)));
    int compound_map_idx = 0;
    for (const auto& x : compound_vec) { compound_map[x] = compound_map_idx++; }

    boost::erase(protein_vec, boost::unique<boost::return_found_end>(boost::sort(protein_vec)));
    int protein_map_idx = 0;
    for (const auto& x : protein_vec) { protein_map[x] = protein_map_idx++; }

    boost::erase(group_vec, boost::unique<boost::return_found_end>(boost::sort(group_vec)));
    int group_map_idx = 0;
    for (const auto& x : group_vec) { group_map[x] = group_map_idx++; }

    boost::erase(peptide_vec, boost::unique<boost::return_found_end>(boost::sort(peptide_vec)));
    int peptide_map_idx = 0;
    for (const auto& x : peptide_vec) { peptide_map[x] = peptide_map_idx++; }

    // Build compound lookup
    std::map<std::string, const OpenSwath::LightCompound*> compound_lookup;
    for (const auto& compound : targeted_exp.compounds)
    {
      compound_lookup[compound.id] = &compound;
    }

    // Insert transitions
    {
      std::stringstream insert_transition_sql, insert_transition_peptide_mapping_sql, insert_transition_precursor_mapping_sql;
      insert_transition_sql.precision(11);

      for (Size i = 0; i < targeted_exp.transitions.size(); i++)
      {
        const auto& tr = targeted_exp.transitions[i];
        int group_set_index = group_map[tr.peptide_ref];

        if (precursor_mz_map.find(group_set_index) == precursor_mz_map.end())
        {
          precursor_mz_map[group_set_index] = tr.precursor_mz;
        }
        if (precursor_decoy_map.find(group_set_index) == precursor_decoy_map.end())
        {
          if (tr.isDetectingTransition())
          {
            precursor_decoy_map[group_set_index] = tr.getDecoy();
          }
        }

        // IPF: Generate transition-peptide mapping tables
        for (const auto& peptidoform : tr.peptidoforms)
        {
          insert_transition_peptide_mapping_sql << "INSERT INTO TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID, PEPTIDE_ID) VALUES (" <<
            i << "," << peptide_map[peptidoform] << "); ";
        }

        // Associate transitions with their precursors
        insert_transition_precursor_mapping_sql << "INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (" <<
          i << "," << group_set_index << "); ";

        std::string transition_charge = "NULL";
        if (tr.fragment_charge != 0)
        {
          transition_charge = String(static_cast<int>(tr.fragment_charge));
        }

        std::string fragment_type_str = tr.getFragmentType();
        std::string fragment_type_char = fragment_type_str.empty() ? "" : fragment_type_str.substr(0, 1);

        // Insert transition data
        insert_transition_sql << "INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ANNOTATION, ORDINAL, " <<
          "DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (" << i << ",'" <<
          tr.transition_name << "'," <<
          tr.product_mz << "," <<
          transition_charge << ",'" <<
          fragment_type_char << "','" <<
          tr.getAnnotation() << "'," <<
          tr.fragment_nr << "," <<
          tr.isDetectingTransition() << "," <<
          tr.isIdentifyingTransition() << "," <<
          tr.isQuantifyingTransition() << "," <<
          tr.library_intensity << "," << tr.getDecoy() << "); ";

        if (i % 50000 == 0 && i > 0)
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
    std::vector<std::pair<int, int>> peptide_protein_map_vec;
    std::vector<std::pair<int, int>> peptide_gene_map_vec;

    // Insert precursors (compounds)
    for (const auto& compound : targeted_exp.compounds)
    {
      int group_set_index = group_map[compound.id];

      if (compound.isPeptide())
      {
        int peptide_set_index = peptide_map[compound.sequence];

        for (const auto& prot_ref : compound.protein_refs)
        {
          if (protein_map.find(prot_ref) != protein_map.end())
          {
            peptide_protein_map_vec.emplace_back(peptide_set_index, protein_map[prot_ref]);
          }
        }

        std::string gene_name = compound.gene_name.empty() ? "NA" : compound.gene_name;
        if (gene_map.find(gene_name) == gene_map.end()) gene_map[gene_name] = (int)gene_map.size();
        peptide_gene_map_vec.emplace_back(peptide_set_index, gene_map[gene_name]);

        insert_precursor_sql <<
          "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, " <<
          "LIBRARY_DRIFT_TIME, LIBRARY_RT, DECOY) VALUES (" <<
          group_set_index << ",'" << compound.id << "','" <<
          compound.peptide_group_label << "'," <<
          precursor_mz_map[group_set_index] << "," <<
          compound.charge <<
          ",NULL," <<
          compound.drift_time << "," <<
          compound.rt << "," <<
          precursor_decoy_map[group_set_index] << "); ";

        insert_precursor_peptide_mapping << "INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (" <<
          group_set_index << "," << peptide_set_index << "); ";
      }
      else
      {
        int compound_set_index = compound_map[compound.id];

        std::string compound_charge = "NULL";
        if (compound.charge != 0)
        {
          compound_charge = String(compound.charge);
        }

        insert_precursor_sql << "INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, " <<
          "LIBRARY_DRIFT_TIME, LIBRARY_RT, DECOY) VALUES (" << group_set_index
          << ",'" << compound.id << "',NULL," <<
          precursor_mz_map[group_set_index] << "," <<
          compound_charge <<
          ",NULL," <<
          compound.drift_time << "," <<
          compound.rt << "," <<
          precursor_decoy_map[group_set_index] << "); ";

        insert_precursor_compound_mapping << "INSERT INTO PRECURSOR_COMPOUND_MAPPING (PRECURSOR_ID, COMPOUND_ID) VALUES (" <<
          group_set_index << "," << compound_set_index << "); ";
      }
    }

    boost::erase(peptide_protein_map_vec, boost::unique<boost::return_found_end>(boost::sort(peptide_protein_map_vec)));
    boost::erase(peptide_gene_map_vec, boost::unique<boost::return_found_end>(boost::sort(peptide_gene_map_vec)));

    // Prepare peptide-gene mapping inserts
    std::stringstream insert_peptide_gene_mapping;
    for (const auto& it : peptide_gene_map_vec)
    {
      insert_peptide_gene_mapping << "INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (" <<
        it.first << "," << it.second << "); ";
    }

    // Prepare gene inserts
    std::stringstream insert_gene_sql;
    std::vector<std::string> gene_vec;
    gene_vec.reserve(gene_map.size());
    for (const auto& it : gene_map)
    {
      gene_vec.push_back(it.first);
    }
    boost::sort(gene_vec);
    for (const auto& gene_name : gene_vec)
    {
      insert_gene_sql << "INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (" <<
        gene_map[gene_name] << ",'" << gene_name << "'," << 0 << "); ";
    }

    // Prepare peptide-protein mapping inserts
    std::stringstream insert_peptide_protein_mapping;
    for (const auto& it : peptide_protein_map_vec)
    {
      insert_peptide_protein_mapping << "INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (" <<
        it.first << "," << it.second << "); ";
    }

    // Prepare protein inserts
    std::stringstream insert_protein_sql;
    for (const auto& protein_id : protein_vec)
    {
      insert_protein_sql << "INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (" <<
        protein_map[protein_id] << ",'" << protein_id << "'," << 0 << "); ";
    }

    // Prepare peptide inserts
    std::stringstream insert_peptide_sql;
    for (const auto& peptide_sequence : peptide_vec)
    {
      std::string unmodified_seq;
      try
      {
        unmodified_seq = AASequence::fromString(peptide_sequence).toUnmodifiedString();
      }
      catch (Exception::InvalidValue&)
      {
        unmodified_seq = peptide_sequence;
      }
      insert_peptide_sql << "INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (" <<
        peptide_map[peptide_sequence] << ",'" <<
        unmodified_seq << "','" <<
        peptide_sequence << "'," << 0 << "); ";
    }

    // Prepare compound inserts
    std::stringstream insert_compound_sql;
    for (const auto& compound_id : compound_vec)
    {
      auto comp_it = compound_lookup.find(compound_id);
      std::string compound_name = compound_id;
      std::string sum_formula;
      std::string smiles;
      std::string adducts;
      if (comp_it != compound_lookup.end())
      {
        compound_name = comp_it->second->compound_name.empty() ? compound_id : comp_it->second->compound_name;
        sum_formula = comp_it->second->sum_formula;
        smiles = comp_it->second->smiles;
        adducts = comp_it->second->adducts;
      }
      insert_compound_sql << "INSERT INTO COMPOUND (ID, COMPOUND_NAME, SUM_FORMULA, SMILES, ADDUCTS, DECOY) VALUES (" <<
        compound_map[compound_id] << ",'" <<
        compound_name << "','" <<
        sum_formula << "','" <<
        smiles << "','" <<
        adducts << "'," <<
        0 << "); ";
    }

    // Prepare decoy updates
    std::stringstream update_decoys_sql;
    update_decoys_sql << "UPDATE PEPTIDE SET DECOY = 1 WHERE ID IN " <<
      "(SELECT PEPTIDE.ID FROM PRECURSOR " <<
      "JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID " <<
      "JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID WHERE PRECURSOR.DECOY = 1); ";
    update_decoys_sql << "UPDATE COMPOUND SET DECOY = 1 WHERE ID IN " <<
      "(SELECT COMPOUND.ID FROM PRECURSOR " <<
      "JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_COMPOUND_MAPPING.PRECURSOR_ID " <<
      "JOIN COMPOUND ON PRECURSOR_COMPOUND_MAPPING.COMPOUND_ID = COMPOUND.ID WHERE PRECURSOR.DECOY = 1); ";
    update_decoys_sql << "UPDATE PROTEIN SET DECOY = 1 WHERE ID IN " <<
      "(SELECT PROTEIN.ID FROM PEPTIDE " <<
      "JOIN PEPTIDE_PROTEIN_MAPPING ON PEPTIDE.ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID " <<
      "JOIN PROTEIN ON PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID = PROTEIN.ID WHERE PEPTIDE.DECOY = 1); ";
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

}
