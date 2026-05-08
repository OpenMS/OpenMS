// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger, Justin Sing $
// $Authors: George Rosenberger, Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OSWFile.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <sqlite3.h>

#include <algorithm>
#include <cstring> // for strcmp
#include <optional>
#include <sstream>
#include <unordered_map>
#include <utility> // for std::move

namespace OpenMS
{
  namespace Sql = Internal::SqliteHelper;
  using namespace std;

  namespace
  {
    void checkSqliteReturnCode_(sqlite3* db, const int rc, const String& action);
    void tryCreateIndex_(SqliteConnector& conn, const String& sql);

    void tryCreateIndexIfTableExists_(SqliteConnector& conn, const String& table_name, const String& sql)
    {
      if (!conn.tableExists(table_name))
      {
        return;
      }
      tryCreateIndex_(conn, sql);
    }

    void tryCreateIndex_(SqliteConnector& conn, const String& sql)
    {
      try
      {
        conn.executeStatement(sql);
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Failed to create SQLite index with statement '" << sql << "': " << e.getMessage() << '\n';
      }
    }

    void requireTable_(SqliteConnector& conn, const String& table_name, const String& message)
    {
      if (!conn.tableExists(table_name))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message + " Missing required table '" + table_name + "'.");
      }
    }

    void createIPFIndexes_(SqliteConnector& conn)
    {
      tryCreateIndexIfTableExists_(conn, "TRANSITION", "CREATE INDEX IF NOT EXISTS idx_transition_id ON TRANSITION (ID);");
      tryCreateIndexIfTableExists_(conn, "PRECURSOR", "CREATE INDEX IF NOT EXISTS idx_precursor_id ON PRECURSOR (ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE", "CREATE INDEX IF NOT EXISTS idx_feature_precursor_id ON FEATURE (PRECURSOR_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE", "CREATE INDEX IF NOT EXISTS idx_feature_id ON FEATURE (ID);");
      tryCreateIndexIfTableExists_(conn, "SCORE_MS2", "CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "SCORE_MS1", "CREATE INDEX IF NOT EXISTS idx_score_ms1_feature_id ON SCORE_MS1 (FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "SCORE_TRANSITION", "CREATE INDEX IF NOT EXISTS idx_score_transition_feature_id ON SCORE_TRANSITION (FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "SCORE_TRANSITION", "CREATE INDEX IF NOT EXISTS idx_score_transition_transition_id ON SCORE_TRANSITION (TRANSITION_ID);");
      tryCreateIndexIfTableExists_(conn, "TRANSITION_PEPTIDE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_transition_peptide_mapping_transition_id ON TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID);");
      tryCreateIndexIfTableExists_(conn, "TRANSITION_PEPTIDE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_transition_peptide_mapping_peptide_id ON TRANSITION_PEPTIDE_MAPPING (PEPTIDE_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_alignment_id ON FEATURE_MS2_ALIGNMENT (ALIGNMENT_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_reference_feature_id ON FEATURE_MS2_ALIGNMENT (REFERENCE_FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_aligned_feature_id ON FEATURE_MS2_ALIGNMENT (ALIGNED_FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_candidate_alignment_id ON FEATURE_MS2_ALIGNMENT_CANDIDATE (ALIGNMENT_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_candidate_reference_feature_id ON FEATURE_MS2_ALIGNMENT_CANDIDATE (REFERENCE_FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_candidate_aligned_feature_id ON FEATURE_MS2_ALIGNMENT_CANDIDATE (ALIGNED_FEATURE_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_candidate_selected ON FEATURE_MS2_ALIGNMENT_CANDIDATE (SELECTED);");
      tryCreateIndexIfTableExists_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "CREATE INDEX IF NOT EXISTS idx_feature_ms2_alignment_candidate_confidence ON FEATURE_MS2_ALIGNMENT_CANDIDATE (MAPPING_CONFIDENCE);");
      tryCreateIndexIfTableExists_(conn, "SCORE_ALIGNMENT", "CREATE INDEX IF NOT EXISTS idx_score_alignment_feature_id ON SCORE_ALIGNMENT (FEATURE_ID);");
    }

    void createLevelContextIndexes_(SqliteConnector& conn, const InferenceLevel level)
    {
      tryCreateIndexIfTableExists_(conn, "PEPTIDE", "CREATE INDEX IF NOT EXISTS idx_peptide_id ON PEPTIDE (ID);");
      tryCreateIndexIfTableExists_(conn, "PRECURSOR_PEPTIDE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_peptide_id ON PRECURSOR_PEPTIDE_MAPPING (PEPTIDE_ID);");
      tryCreateIndexIfTableExists_(conn, "PRECURSOR_PEPTIDE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_precursor_peptide_mapping_precursor_id ON PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID);");
      tryCreateIndexIfTableExists_(conn, "PRECURSOR", "CREATE INDEX IF NOT EXISTS idx_precursor_id_level_ctx ON PRECURSOR (ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE", "CREATE INDEX IF NOT EXISTS idx_feature_precursor_id_level_ctx ON FEATURE (PRECURSOR_ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE", "CREATE INDEX IF NOT EXISTS idx_feature_id_level_ctx ON FEATURE (ID);");
      tryCreateIndexIfTableExists_(conn, "FEATURE", "CREATE INDEX IF NOT EXISTS idx_feature_run_id_level_ctx ON FEATURE (RUN_ID);");
      tryCreateIndexIfTableExists_(conn, "SCORE_MS2", "CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id_level_ctx ON SCORE_MS2 (FEATURE_ID);");
      if (level == InferenceLevel::Protein)
      {
        tryCreateIndexIfTableExists_(conn, "PEPTIDE_PROTEIN_MAPPING", "CREATE INDEX IF NOT EXISTS idx_peptide_protein_mapping_peptide_id ON PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID);");
        tryCreateIndexIfTableExists_(conn, "PEPTIDE_PROTEIN_MAPPING", "CREATE INDEX IF NOT EXISTS idx_peptide_protein_mapping_protein_id ON PEPTIDE_PROTEIN_MAPPING (PROTEIN_ID);");
      }
      else if (level == InferenceLevel::Gene)
      {
        tryCreateIndexIfTableExists_(conn, "PEPTIDE_GENE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_peptide_gene_mapping_peptide_id ON PEPTIDE_GENE_MAPPING (PEPTIDE_ID);");
        tryCreateIndexIfTableExists_(conn, "PEPTIDE_GENE_MAPPING", "CREATE INDEX IF NOT EXISTS idx_peptide_gene_mapping_gene_id ON PEPTIDE_GENE_MAPPING (GENE_ID);");
      }
    }

    String prepareOutputFile_(const String& input_filename, const String& output_filename)
    {
      const String target_filename = output_filename.empty() ? input_filename : output_filename;
      if (target_filename != input_filename)
      {
        if (File::exists(target_filename) && !File::remove(target_filename))
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, target_filename);
        }
        if (!File::copy(input_filename, target_filename))
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, target_filename);
        }
      }
      return target_filename;
    }

    void bindNullableInt64_(sqlite3* db, sqlite3_stmt* stmt, const int column, const std::optional<Int64>& value)
    {
      int rc = SQLITE_OK;
      if (value.has_value())
      {
        rc = sqlite3_bind_int64(stmt, column, *value);
      }
      else
      {
        rc = sqlite3_bind_null(stmt, column);
      }
      checkSqliteReturnCode_(db, rc, "Failed to bind nullable INTEGER value");
    }

    String scoreTableName_(const InferenceLevel level)
    {
      switch (level)
      {
        case InferenceLevel::Peptide: return "SCORE_PEPTIDE";
        case InferenceLevel::Protein: return "SCORE_PROTEIN";
        case InferenceLevel::Gene: return "SCORE_GENE";
        default:
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level-context score tables only exist for peptide, protein, and gene inference.");
      }
    }

    String entityIdColumnName_(const InferenceLevel level)
    {
      switch (level)
      {
        case InferenceLevel::Peptide: return "PEPTIDE_ID";
        case InferenceLevel::Protein: return "PROTEIN_ID";
        case InferenceLevel::Gene: return "GENE_ID";
        default:
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level-context score tables only exist for peptide, protein, and gene inference.");
      }
    }

    void checkSqliteReturnCode_(sqlite3* db, const int rc, const String& action)
    {
      if (rc != SQLITE_OK)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, action + ": " + sqlite3_errmsg(db));
      }
    }

    void executePreparedStatement_(sqlite3* db, sqlite3_stmt* stmt, const String& action)
    {
      const int rc = sqlite3_step(stmt);
      if (rc != SQLITE_DONE)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, action + ": " + sqlite3_errmsg(db));
      }
      sqlite3_reset(stmt);
      sqlite3_clear_bindings(stmt);
    }
  } // namespace

  const std::array<std::string, (Size)OSWFile::OSWLevel::SIZE_OF_OSWLEVEL> OSWFile::names_of_oswlevel = { "ms1", "ms2", "transition" };

  void OSWFile::readToPIN(const std::string& in_osw,
                     const OSWLevel osw_level,
                     std::ostream& pin_output,
                     const double ipf_max_peakgroup_pep,
                     const double ipf_max_transition_isotope_overlap,
                     const double ipf_min_transition_sn)
  {
      sqlite3_stmt * stmt;
      std::string select_sql;

      // Open database
      SqliteConnector conn(in_osw);

      if (osw_level == OSWLevel::MS1)
      {
        select_sql = "SELECT *, RUN_ID || '_' || PRECURSOR.ID AS GROUP_ID " \
                      "FROM FEATURE_MS1 "\
                      "INNER JOIN (SELECT ID, PRECURSOR_ID, RUN_ID FROM FEATURE) AS FEATURE ON FEATURE_ID = FEATURE.ID "\
                      "INNER JOIN (SELECT ID, DECOY FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID "\
                      "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "\
                      "INNER JOIN (SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE) AS PEPTIDE ON "\
                        "PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID;";
      }
      else if (osw_level == OSWLevel::TRANSITION)
      {
        select_sql = "SELECT TRANSITION.DECOY AS DECOY, FEATURE_TRANSITION.*, "\
                        "RUN_ID || '_' || FEATURE_TRANSITION.FEATURE_ID || '_' || PRECURSOR_ID || '_' || TRANSITION_ID AS GROUP_ID, "\
                        "FEATURE_TRANSITION.FEATURE_ID || '_' || FEATURE_TRANSITION.TRANSITION_ID AS FEATURE_ID, "\
                        "'PEPTIDE' AS MODIFIED_SEQUENCE FROM FEATURE_TRANSITION "\
                        "INNER JOIN (SELECT RUN_ID, ID, PRECURSOR_ID FROM FEATURE) AS FEATURE ON FEATURE_TRANSITION.FEATURE_ID = FEATURE.ID " \
                        "INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID "\
                        "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "\
                        "INNER JOIN (SELECT ID, DECOY FROM TRANSITION) AS TRANSITION ON FEATURE_TRANSITION.TRANSITION_ID = TRANSITION.ID "\
                        "WHERE PEP <= " + OpenMS::String(ipf_max_peakgroup_pep) +
                          " AND VAR_ISOTOPE_OVERLAP_SCORE <= " + OpenMS::String(ipf_max_transition_isotope_overlap) +
                          " AND VAR_LOG_SN_SCORE > " + OpenMS::String(ipf_min_transition_sn) +
                          " AND PRECURSOR.DECOY == 0 ORDER BY FEATURE_ID, PRECURSOR_ID, TRANSITION_ID;";
      }
      else
      {
        // Peak group-level query including peptide sequence
        select_sql = "SELECT *, RUN_ID || '_' || PRECURSOR.ID AS GROUP_ID "\
                      "FROM FEATURE_MS2 "\
                      "INNER JOIN (SELECT ID, PRECURSOR_ID, RUN_ID FROM FEATURE) AS FEATURE ON FEATURE_ID = FEATURE.ID "\
                      "INNER JOIN (SELECT ID, DECOY FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID "\
                      "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "\
                      "INNER JOIN (SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE) AS PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID;";
      }

      // Execute SQL select statement
      conn.prepareStatement(&stmt, select_sql);
      sqlite3_step(stmt);

      int cols = sqlite3_column_count(stmt);

      // Generate features
      int k = 0;
      std::vector<std::string> group_id_index;
      OpenMS::String tmp;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        std::string psm_id;
        size_t scan_id = 0;
        int label = 0;
        std::string peptide;
        std::map<std::string, double> features;

        for (int i = 0; i < cols; i++)
        {
          if (strcmp(sqlite3_column_name(stmt, i), "FEATURE_ID") == 0)
          {
            Sql::extractValue<string>(&psm_id, stmt, i);
          }
          if (strcmp(sqlite3_column_name(stmt, i), "GROUP_ID") == 0)
          {
            std::string group_id(reinterpret_cast<const char*>(sqlite3_column_text(stmt, i)));
            auto it = std::find(group_id_index.begin(), group_id_index.end(), group_id);
            if (it != group_id_index.end())
            {
              scan_id = it - group_id_index.begin();
            }
            else
            {
              scan_id = group_id_index.size();
              group_id_index.emplace_back(group_id);
            }
          }
          if (strcmp(sqlite3_column_name(stmt, i), "DECOY") == 0)
          {
            if (sqlite3_column_int( stmt, i ) == 1)
            {
              label = -1;
            }
            else
            {
              label = 1;
            }
          }
          if (strcmp(sqlite3_column_name( stmt, i ), "MODIFIED_SEQUENCE") == 0)
          {
            Sql::extractValue<string>(&peptide, stmt, i);
          }
          if (strncmp(sqlite3_column_name( stmt, i ), "VAR_", 4) == 0)
          {
            features[OpenMS::String(sqlite3_column_name( stmt, i ))] = sqlite3_column_double( stmt, i );
          }
        }

        // Write output
        if (k == 0)
        {
          pin_output << "PSMId\tLabel\tScanNr";
          for (auto const &feat : features)
          {
            pin_output << "\t" << feat.first;
          }
          pin_output << "\tPeptide\tProteins\n";
        }
        pin_output << psm_id << "\t" << label << "\t" << scan_id;
        for (auto const &feat : features)
        {
          pin_output << "\t" << feat.second;
        }
        pin_output << "\t." << peptide << ".\tProt1" << "\n";

        sqlite3_step( stmt );
        k++;
      }

      sqlite3_finalize(stmt);

      if (k==0)
      {
        if (osw_level == OSWLevel::TRANSITION)
        {
          throw Exception::Precondition(__FILE__, __LINE__, __FUNCTION__,
              OpenMS::String("PercolatorAdapter needs to be applied on MS1 & MS2 levels before conducting transition-level scoring."));
        }
        else
        {
          throw Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, in_osw);
        }
      }

    }

    void OSWFile::writeFromPercolator(const std::string& in_osw,
                        const OSWFile::OSWLevel osw_level,
                        const std::map< std::string, PercolatorFeature >& features)
    {
      std::string table;
      std::string create_sql;

      if (osw_level == OSWLevel::MS1)
      {
        table = "SCORE_MS1";
        create_sql =  "DROP TABLE IF EXISTS " + table + "; " \
                      "CREATE TABLE " + table + "(" \
                      "FEATURE_ID INT NOT NULL," \
                      "SCORE DOUBLE NOT NULL," \
                      "QVALUE DOUBLE NOT NULL," \
                      "PEP DOUBLE NOT NULL);";

      }
      else if (osw_level == OSWLevel::TRANSITION)
      {
        table = "SCORE_TRANSITION";
        create_sql =  "DROP TABLE IF EXISTS " + table + "; " \
                      "CREATE TABLE " + table + "(" \
                      "FEATURE_ID INT NOT NULL," \
                      "TRANSITION_ID INT NOT NULL," \
                      "SCORE DOUBLE NOT NULL," \
                      "QVALUE DOUBLE NOT NULL," \
                      "PEP DOUBLE NOT NULL);";

      }
      else
      {
        table = "SCORE_MS2";
        create_sql =  "DROP TABLE IF EXISTS " + table + "; " \
                      "CREATE TABLE " + table + "(" \
                      "FEATURE_ID INT NOT NULL," \
                      "SCORE DOUBLE NOT NULL," \
                      "QVALUE DOUBLE NOT NULL," \
                      "PEP DOUBLE NOT NULL);";
      }

      std::vector<std::string> insert_sqls;
      for (auto const &feat : features)
      {
        std::stringstream insert_sql;
        insert_sql << "INSERT INTO " << table;
        if (osw_level == OSWLevel::TRANSITION)
        {
          std::vector<String> ids;
          String(feat.first).split("_", ids);
          insert_sql << " (FEATURE_ID, TRANSITION_ID, SCORE, QVALUE, PEP) VALUES (";
          insert_sql <<  ids[0] << ",";
          insert_sql <<  ids[1] << ",";
        }
        else
        {
          insert_sql << " (FEATURE_ID, SCORE, QVALUE, PEP) VALUES (";
          insert_sql <<  feat.first << ",";
        }
        insert_sql << feat.second.score << ",";
        insert_sql << feat.second.qvalue << ",";
        insert_sql << feat.second.posterior_error_prob << "); ";

        insert_sqls.push_back(insert_sql.str());
      }

      // Write to Sqlite database
      SqliteConnector conn(in_osw);
      conn.executeStatement(create_sql);
      conn.executeStatement("BEGIN TRANSACTION");
      for (size_t i = 0; i < insert_sqls.size(); i++)
      {
        conn.executeStatement(insert_sqls[i]);
      }
      conn.executeStatement("END TRANSACTION");
    }

    OSWFile::OSWFile(const String& filename)
      : filename_(filename),
        conn_(filename, SqliteConnector::SqlOpenMode::READ_ONLY)
    {
      has_SCOREMS2_ = conn_.tableExists("SCORE_MS2");
    }

    void OSWFile::readMinimal(OSWData& swath_result)
    {
      readMeta_(swath_result);

      readTransitions_(swath_result);

      String select_sql = "select PROTEIN.ID as prot_id, PROTEIN_ACCESSION as prot_accession from PROTEIN order by prot_id";
      sqlite3_stmt* stmt;
      conn_.prepareStatement(&stmt, select_sql);
      enum CBIG
      { // indices of respective columns in the query above
        I_PROTID,
        I_ACCESSION,
        SIZE_OF_CBIG
      };
      Sql::SqlState rc = Sql::nextRow(stmt);
      if (sqlite3_column_count(stmt) != SIZE_OF_CBIG)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Query was changed! Please report this bug!");
      }
      String accession;
      // protein loop
      while (rc == Sql::SqlState::SQL_ROW)
      {
        int id = Sql::extractInt(stmt, I_PROTID);
        accession = Sql::extractString(stmt, I_ACCESSION);
        swath_result.addProtein(OSWProtein(std::move(accession), id, {}));
        rc = Sql::nextRow(stmt, rc); // next row
      }

      sqlite3_finalize(stmt);
    }


    /**
    @brief populates the @p index'th protein with Peptides, unless the protein already has peptides

    Internally uses the proteins ID to search for cross referencing peptides and transitions in the OSW file.
    @throws Exception::InvalidValue if the ID is unknown
    */

    void OSWFile::readProtein(OSWData& swath_result, const Size index)
    {
      if (!swath_result.getProteins()[index].getPeptidePrecursors().empty())
      { // already populated
        return;
      }
      getFullProteins_(swath_result, index);
      if (swath_result.getProteins()[index].getPeptidePrecursors().empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ID is not known in OSWFile " + filename_, String(swath_result.getProteins()[index].getID()));
      }
    }

    void OSWFile::read(OSWData& swath_result)
    {
      readMeta_(swath_result);
      readTransitions_(swath_result);
      getFullProteins_(swath_result);      
    }


    enum ColProteinSelect
    { // indices of respective columns in the query below
      I_PROTID,
      I_ACCESSION,
      I_DECOY,
      I_MODSEQ,
      I_PRECID,
      I_PRECMZ,
      I_PRECQ,
      I_FEATID,
      I_EXPRT,
      I_DELTART,
      I_RTLEFT,
      I_RTRIGHT,
      I_TRID,
      I_QVALUE,
      SIZE_OF_ColProteinSelect
    };

    /// represents the state of an SQL row, which is
    /// updated partially whenever nested structures change
    struct LineState
    {
      // Layers of information. Whenever the id changes, we know a new item has begun
      // ... PROTEIN
      int prot_id;
      String accession;
      bool decoy;

      void setProt(sqlite3_stmt* stmt)
      {
        prot_id = Sql::extractInt(stmt, I_PROTID);
        accession = Sql::extractString(stmt, I_ACCESSION);
        decoy = Sql::extractBool(stmt, I_DECOY);
      }
      void updateProt(LineState& new_line)
      {
        prot_id = new_line.prot_id;
        accession = std::move(new_line.accession);
        decoy = new_line.decoy;
      }

      // ... PRECURSOR
      int prec_id;
      String seq;
      short chargePC;
      float precmz;

      void setPC(sqlite3_stmt* stmt)
      {
        prec_id = Sql::extractInt(stmt, I_PRECID);
        seq = Sql::extractString(stmt, I_MODSEQ);
        chargePC = (short)Sql::extractInt(stmt, I_PRECQ);
        precmz = Sql::extractFloat(stmt, I_PRECMZ);
      }
      void updatePC(LineState& new_line)
      {
        prec_id = new_line.prec_id;
        seq = std::move(new_line.seq);
        chargePC = new_line.chargePC;
        precmz = new_line.precmz;
      }

      // ... FEATURE
      Int64 feat_id; // in SQL, feature_id is a 63-bit integer...
      float rt_exp;
      float rt_lw;
      float rt_rw;
      float rt_delta;
      float qvalue;
      void setFeature(sqlite3_stmt* stmt)
      {
        feat_id = Sql::extractInt64(stmt, I_FEATID);
        rt_exp = Sql::extractFloat(stmt, I_EXPRT);
        rt_lw = Sql::extractFloat(stmt, I_RTLEFT);
        rt_rw = Sql::extractFloat(stmt, I_RTRIGHT);
        rt_delta = Sql::extractFloat(stmt, I_DELTART);
        qvalue = Sql::extractFloat(stmt, I_QVALUE);
      }
      void updateFeat(const LineState& new_line)
      {
        feat_id = new_line.feat_id;
        rt_exp = new_line.rt_exp;
        rt_lw = new_line.rt_lw;
        rt_rw = new_line.rt_rw;
        rt_delta = new_line.rt_delta;
        qvalue = new_line.qvalue;
      }
    };

    void initLine(LineState& current, sqlite3_stmt* stmt)
    {
      current.setProt(stmt);
      current.setPC(stmt);
      current.setFeature(stmt);
    }


    bool nextProtein(OSWProtein& prot, sqlite3_stmt* stmt, Sql::SqlState& rc, LineState& old_line)
    {
      LineState new_line;
      // PROTEIN
      std::vector<OSWPeptidePrecursor> precursors;
      OSWPeptidePrecursor new_pc;
      auto check_add_protein = [&](bool add_force = false)
      {
        precursors.push_back(new_pc); // the last PC already belonged to the old protein
        if (old_line.prot_id != new_line.prot_id || add_force)
        {
          prot = OSWProtein(old_line.accession, old_line.prot_id, std::move(precursors));
          old_line.updateProt(new_line);
          precursors.clear();
          return true;
        }
        return false;
      };
      // ... PRECURSOR
      std::vector<OSWPeakGroup> features;
      OSWPeakGroup new_feature;
      auto check_add_pc = [&](bool add_force = false)
      {
        features.push_back(std::move(new_feature)); // the last feature belonged to the old PC
        if (old_line.prec_id != new_line.prec_id || add_force)
        {
          new_pc = OSWPeptidePrecursor(old_line.seq, old_line.chargePC, old_line.decoy, old_line.precmz, std::move(features));
          old_line.updatePC(new_line);
          features.clear();
          return true;
        }
        return false;
      };

      // ... FEATURE
      std::vector<UInt32> transition_ids;
      UInt32 new_transition;
      auto check_add_feat = [&](bool add_force = false)
      {
        if (old_line.feat_id != new_line.feat_id || add_force)
        {
          new_feature = OSWPeakGroup(old_line.rt_exp, old_line.rt_lw, old_line.rt_rw, old_line.rt_delta, std::move(transition_ids), old_line.qvalue);
          old_line.updateFeat(new_line);
          transition_ids.clear();
          return true;
        }
        else
        { // if we enter the above block, we will parse the same sql row in the next iteration, so only add the tr-ID if its not a new block
          transition_ids.push_back(new_transition); // the current transition belongs to the current feature...
        }
        return false;
      };

      // protein loop
      while (rc == Sql::SqlState::SQL_ROW)
      {
        // precursor loop (peptide with charge)
        while (rc == Sql::SqlState::SQL_ROW)
        {
          // feature loop
          while (rc == Sql::SqlState::SQL_ROW)
          {
            new_transition = Sql::extractInt(stmt, I_TRID);
            new_line.setFeature(stmt);
            if (check_add_feat())
            {
              break; // new feature just started?--> check if new PC started as well.
            }
            rc = Sql::nextRow(stmt, rc); // next row
          }
          if (rc != Sql::SqlState::SQL_ROW) {
            // we are beyond last row; new feature is not yet made; so we forcibly do it now
            check_add_feat(true);    // add last feature
            check_add_pc(true);      // add last precursor
            check_add_protein(true); // add last protein
            return false; // this was the last protein
          }
          new_line.setPC(stmt);
          if (check_add_pc())
          {
            break; // new PC just started?--> check if if new protein started as well.
          }
        }
        new_line.setProt(stmt);
        if (check_add_protein())
        {
          return true; // current protein ended... but there are more..
        }
      }

      // we did not even enter the while-loops... so no data was there (but should have been)
      throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No rows available. Please report this as a bug!");
    }

    void OSWFile::getFullProteins_(OSWData& swath_result, Size index)
    {
      String protein_subselect;
      if (index == ALL_PROTEINS)
      {
        swath_result.clearProteins();
        protein_subselect = "PROTEIN";
      }
      else
      { //  do not use accession to filter -- its as slow as full query
        protein_subselect = "(select * from PROTEIN  where ID = " + String(swath_result.getProteins().at(index).getID()) + ") as PROTEIN";
      }
     

      // check of SCORE_MS2 table is available (for OSW files which underwent pyProphet)
      // set q_value to -1 if missing
      String MS2_select = (has_SCOREMS2_ ? "SCORE_MS2.QVALUE as qvalue" : "-1 as qvalue");
      String MS2_join = (has_SCOREMS2_ ? "inner join(select * from SCORE_MS2) as SCORE_MS2 on SCORE_MS2.FEATURE_ID = FEATURE.ID" : "");

      // assemble the protein-PeptidePrecursor-Feature hierarchy
      // note: when changing the query, make sure to keep the indices in ColProteinSelect in sync!!!
      String select_sql = "select PROTEIN.ID as prot_id, PROTEIN_ACCESSION as prot_accession, PROTEIN.DECOY as decoy, "
                          "       PEPTIDE.MODIFIED_SEQUENCE as modified_sequence,"
                          "       PRECURSOR.ID as prec_id, PRECURSOR.PRECURSOR_MZ as pc_mz, PRECURSOR.CHARGE as pc_charge,"
                          "       FEATURE.ID as feat_id, FEATURE.EXP_RT as rt_experimental, FEATURE.DELTA_RT as rt_delta, FEATURE.LEFT_WIDTH as rt_left_width, FEATURE.RIGHT_WIDTH as rt_right_width,"
                          "       FeatTrMap.TRANSITION_ID as tr_id, " +
        MS2_select +
        " FROM " + protein_subselect +
        " inner join(select* FROM PEPTIDE_PROTEIN_MAPPING) as PepProtMap on PepProtMap.PROTEIN_ID = PROTEIN.ID "
        " inner join(select ID, MODIFIED_SEQUENCE FROM PEPTIDE) as PEPTIDE on PEPTIDE.ID = PepProtMap.PEPTIDE_ID "
        " inner join(select * FROM PRECURSOR_PEPTIDE_MAPPING) as PrePepMap on PrePepMap.PEPTIDE_ID = PEPTIDE.ID "
        " inner join(select * from PRECURSOR) as PRECURSOR on PRECURSOR.ID = PrePepMap.PRECURSOR_ID "
        " inner join(select * from FEATURE) as FEATURE on FEATURE.PRECURSOR_ID = PRECURSOR.ID "
        " inner join(select * from FEATURE_TRANSITION) as FeatTrMap on FeatTrMap.FEATURE_ID = FEATURE.ID " +
        MS2_join +
        " order by prot_id, prec_id, feat_id, qvalue, tr_id ";


      sqlite3_stmt* stmt;
      conn_.prepareStatement(&stmt, select_sql);

      Sql::SqlState rc = Sql::nextRow(stmt);
      if (sqlite3_column_count(stmt) != SIZE_OF_ColProteinSelect)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Query was changed! Please report this bug!");
      }

      if (rc == Sql::SqlState::SQL_DONE)
      { // no data
        return;
      }

      LineState current_line;
      initLine(current_line, stmt);
      OSWProtein prot;

      if (index == ALL_PROTEINS)
      {
        bool has_more;
        do
        {
          has_more = nextProtein(prot, stmt, rc, current_line);
          swath_result.addProtein(std::move(prot));
        } while (has_more);
      }
      else // single protein
      {
        nextProtein(prot, stmt, rc, current_line);
        swath_result.setProtein(index, std::move(prot));
      }

      sqlite3_finalize(stmt);
    }

    void OSWFile::readMeta_(OSWData& data)
    {
      data.setSqlSourceFile(filename_);
      data.setRunID(getRunID());
    }

    std::map<Int64, String> OSWFile::readRunBasenames() const
    {
      SqliteConnector conn(filename_);
      requireTable_(conn, "RUN", "Run-name lookup requires run metadata.");

      const String query = "SELECT ID, FILENAME FROM RUN ORDER BY ID;";
      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, query);

      std::map<Int64, String> run_names;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        const Int64 run_id = Sql::extractInt64(stmt, 0);
        const String filename = Sql::extractString(stmt, 1);
        String run_name = File::stemName(filename);
        if (run_name.empty())
        {
          run_name = File::basename(filename);
        }
        if (run_name.empty())
        {
          run_name = "RUN_ID " + String(run_id);
        }
        run_names[run_id] = std::move(run_name);
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);
      return run_names;
    }

    UInt64 OSWFile::getRunID() const
    {
      SqliteConnector conn(filename_);
      Size nr_results = 0;

      std::string select_sql = "SELECT RUN.ID FROM RUN;";

      sqlite3_stmt* stmt;
      conn.prepareStatement(&stmt, select_sql);
      Sql::SqlState state = Sql::SqlState::SQL_ROW;
      UInt64 id;
      while ((state = Sql::nextRow(stmt, state)) == Sql::SqlState::SQL_ROW)
      {
        ++nr_results;
        id = Sql::extractInt64(stmt, 0);
      }
      // free memory
      sqlite3_finalize(stmt);

      if (nr_results != 1)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "File '" + filename_ + "' contains more than one run. This is currently not supported!");
      }
      return id;
    }

    std::vector<IPFPrecursorRow> OSWFile::readIPFPrecursorData(const PeptidoformInferenceConfig& config) const
    {
      SqliteConnector conn(filename_, SqliteConnector::SqlOpenMode::READWRITE);
      createIPFIndexes_(conn);

      requireTable_(conn, "PRECURSOR", "Peptidoform inference requires precursor information.");
      requireTable_(conn, "FEATURE", "Peptidoform inference requires candidate feature annotations.");
      requireTable_(conn, "SCORE_MS2", "Peptidoform inference requires MS2 peakgroup scores.");
      if (config.ipf_ms1_scoring)
      {
        requireTable_(conn, "SCORE_MS1", "Peptidoform inference with MS1 precursor scoring requires MS1 scores.");
      }
      if (config.ipf_ms1_scoring || config.ipf_ms2_scoring)
      {
        requireTable_(conn, "SCORE_TRANSITION", "Peptidoform inference requires transition-level scores.");
      }

      String query;
      if (!config.ipf_ms1_scoring && config.ipf_ms2_scoring)
      {
        query =
          "SELECT FEATURE.ID AS FEATURE_ID, "
          "       SCORE_MS2.PEP AS MS2_PEAKGROUP_PEP, "
          "       NULL AS MS1_PRECURSOR_PEP, "
          "       SCORE_TRANSITION.PEP AS MS2_PRECURSOR_PEP "
          "FROM PRECURSOR "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "INNER JOIN ("
          "  SELECT FEATURE_ID, PEP "
          "  FROM SCORE_TRANSITION "
          "  INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
          "  WHERE TRANSITION.TYPE = '' AND TRANSITION.DECOY = 0"
          ") AS SCORE_TRANSITION ON FEATURE.ID = SCORE_TRANSITION.FEATURE_ID "
          "WHERE PRECURSOR.DECOY = 0 AND SCORE_MS2.PEP < ? "
          "ORDER BY FEATURE.ID;";
      }
      else if (config.ipf_ms1_scoring && !config.ipf_ms2_scoring)
      {
        query =
          "SELECT FEATURE.ID AS FEATURE_ID, "
          "       SCORE_MS2.PEP AS MS2_PEAKGROUP_PEP, "
          "       SCORE_MS1.PEP AS MS1_PRECURSOR_PEP, "
          "       NULL AS MS2_PRECURSOR_PEP "
          "FROM PRECURSOR "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS1 ON FEATURE.ID = SCORE_MS1.FEATURE_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "WHERE PRECURSOR.DECOY = 0 AND SCORE_MS2.PEP < ? "
          "ORDER BY FEATURE.ID;";
      }
      else if (config.ipf_ms1_scoring && config.ipf_ms2_scoring)
      {
        query =
          "SELECT FEATURE.ID AS FEATURE_ID, "
          "       SCORE_MS2.PEP AS MS2_PEAKGROUP_PEP, "
          "       SCORE_MS1.PEP AS MS1_PRECURSOR_PEP, "
          "       SCORE_TRANSITION.PEP AS MS2_PRECURSOR_PEP "
          "FROM PRECURSOR "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS1 ON FEATURE.ID = SCORE_MS1.FEATURE_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "INNER JOIN ("
          "  SELECT FEATURE_ID, PEP "
          "  FROM SCORE_TRANSITION "
          "  INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
          "  WHERE TRANSITION.TYPE = '' AND TRANSITION.DECOY = 0"
          ") AS SCORE_TRANSITION ON FEATURE.ID = SCORE_TRANSITION.FEATURE_ID "
          "WHERE PRECURSOR.DECOY = 0 AND SCORE_MS2.PEP < ? "
          "ORDER BY FEATURE.ID;";
      }
      else
      {
        query =
          "SELECT FEATURE.ID AS FEATURE_ID, "
          "       SCORE_MS2.PEP AS MS2_PEAKGROUP_PEP, "
          "       NULL AS MS1_PRECURSOR_PEP, "
          "       NULL AS MS2_PRECURSOR_PEP "
          "FROM PRECURSOR "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "WHERE PRECURSOR.DECOY = 0 AND SCORE_MS2.PEP < ? "
          "ORDER BY FEATURE.ID;";
      }

      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, query);
      checkSqliteReturnCode_(conn.getDB(), sqlite3_bind_double(stmt, 1, config.ipf_max_peakgroup_pep),
                             "Failed to bind ipf_max_peakgroup_pep for precursor query");

      std::vector<IPFPrecursorRow> rows;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        IPFPrecursorRow row;
        row.feature_id = Sql::extractInt64(stmt, 0);
        row.ms2_peakgroup_pep = Sql::extractDouble(stmt, 1);
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL)
        {
          row.ms1_precursor_pep = Sql::extractDouble(stmt, 2);
        }
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL)
        {
          row.ms2_precursor_pep = Sql::extractDouble(stmt, 3);
        }
        rows.push_back(std::move(row));
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      OPENMS_LOG_INFO << "Read " << rows.size() << " precursor/IPF evidence rows." << std::endl;
      return rows;
    }

    std::vector<IPFTransitionRow> OSWFile::readIPFTransitionData(const PeptidoformInferenceConfig& config) const
    {
      SqliteConnector conn(filename_, SqliteConnector::SqlOpenMode::READWRITE);
      createIPFIndexes_(conn);

      requireTable_(conn, "SCORE_TRANSITION", "Peptidoform inference requires transition-level scores.");
      requireTable_(conn, "TRANSITION", "Peptidoform inference requires transition information.");
      requireTable_(conn, "TRANSITION_PEPTIDE_MAPPING", "Peptidoform inference requires transition-to-peptidoform mappings.");

      std::unordered_map<Int64, Int32> num_peptidoforms;
      {
        const String num_query =
          "SELECT SCORE_TRANSITION.FEATURE_ID, COUNT(DISTINCT TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID) AS NUM_PEPTIDOFORMS "
          "FROM SCORE_TRANSITION "
          "INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
          "INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION.ID = TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID "
          "WHERE TRANSITION.TYPE != '' AND TRANSITION.DECOY = 0 "
          "GROUP BY SCORE_TRANSITION.FEATURE_ID "
          "ORDER BY SCORE_TRANSITION.FEATURE_ID;";
        sqlite3_stmt* stmt = nullptr;
        conn.prepareStatement(&stmt, num_query);
        Sql::SqlState state = Sql::nextRow(stmt);
        while (state == Sql::SqlState::SQL_ROW)
        {
          num_peptidoforms[Sql::extractInt64(stmt, 0)] = static_cast<Int32>(Sql::extractInt(stmt, 1));
          state = Sql::nextRow(stmt, state);
        }
        sqlite3_finalize(stmt);
      }

      std::vector<std::pair<Int64, Int64>> bitmask_pairs;
      {
        const String bitmask_query =
          "SELECT DISTINCT TRANSITION.ID AS TRANSITION_ID, TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID "
          "FROM SCORE_TRANSITION "
          "INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
          "INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION.ID = TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID "
          "WHERE TRANSITION.TYPE != '' AND TRANSITION.DECOY = 0 "
          "ORDER BY TRANSITION.ID, TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID;";
        sqlite3_stmt* stmt = nullptr;
        conn.prepareStatement(&stmt, bitmask_query);
        Sql::SqlState state = Sql::nextRow(stmt);
        while (state == Sql::SqlState::SQL_ROW)
        {
          bitmask_pairs.emplace_back(Sql::extractInt64(stmt, 0), Sql::extractInt64(stmt, 1));
          state = Sql::nextRow(stmt, state);
        }
        sqlite3_finalize(stmt);
      }

      const String evidence_query =
        "SELECT SCORE_TRANSITION.FEATURE_ID, SCORE_TRANSITION.TRANSITION_ID, SCORE_TRANSITION.PEP "
        "FROM SCORE_TRANSITION "
        "INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
        "WHERE TRANSITION.TYPE != '' AND TRANSITION.DECOY = 0 AND SCORE_TRANSITION.PEP < ? "
        "ORDER BY SCORE_TRANSITION.FEATURE_ID, SCORE_TRANSITION.TRANSITION_ID;";
      const String candidate_query =
        "SELECT DISTINCT SCORE_TRANSITION.FEATURE_ID, TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID "
        "FROM SCORE_TRANSITION "
        "INNER JOIN TRANSITION ON SCORE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
        "INNER JOIN TRANSITION_PEPTIDE_MAPPING ON TRANSITION.ID = TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID "
        "WHERE TRANSITION.TYPE != '' AND TRANSITION.DECOY = 0 "
        "ORDER BY SCORE_TRANSITION.FEATURE_ID, TRANSITION_PEPTIDE_MAPPING.PEPTIDE_ID;";

      sqlite3_stmt* evidence_stmt = nullptr;
      sqlite3_stmt* candidate_stmt = nullptr;
      conn.prepareStatement(&evidence_stmt, evidence_query);
      conn.prepareStatement(&candidate_stmt, candidate_query);
      checkSqliteReturnCode_(conn.getDB(), sqlite3_bind_double(evidence_stmt, 1, config.ipf_max_transition_pep),
                             "Failed to bind ipf_max_transition_pep for transition query");

      struct CompactEvidenceRow
      {
        Int64 transition_id = -1;
        double pep = 1.0;
      };

      auto loadEvidenceGroup = [&](Int64& feature_id, std::vector<CompactEvidenceRow>& rows, Sql::SqlState& state) -> bool
      {
        rows.clear();
        if (state != Sql::SqlState::SQL_ROW)
        {
          return false;
        }
        feature_id = Sql::extractInt64(evidence_stmt, 0);
        while (state == Sql::SqlState::SQL_ROW && Sql::extractInt64(evidence_stmt, 0) == feature_id)
        {
          rows.push_back({Sql::extractInt64(evidence_stmt, 1), Sql::extractDouble(evidence_stmt, 2)});
          state = Sql::nextRow(evidence_stmt, state);
        }
        return true;
      };

      auto loadCandidateGroup = [&](Int64& feature_id, std::vector<Int64>& rows, Sql::SqlState& state) -> bool
      {
        rows.clear();
        if (state != Sql::SqlState::SQL_ROW)
        {
          return false;
        }
        feature_id = Sql::extractInt64(candidate_stmt, 0);
        while (state == Sql::SqlState::SQL_ROW && Sql::extractInt64(candidate_stmt, 0) == feature_id)
        {
          rows.push_back(Sql::extractInt64(candidate_stmt, 1));
          state = Sql::nextRow(candidate_stmt, state);
        }
        return true;
      };

      std::vector<IPFTransitionRow> rows;
      Size evidence_row_count = 0;
      Size output_feature_count = 0;
      Sql::SqlState evidence_state = Sql::nextRow(evidence_stmt);
      Sql::SqlState candidate_state = Sql::nextRow(candidate_stmt);
      Int64 evidence_feature_id = -1;
      Int64 candidate_feature_id = -1;
      std::vector<CompactEvidenceRow> evidence_rows;
      std::vector<Int64> candidate_rows;
      bool has_evidence = loadEvidenceGroup(evidence_feature_id, evidence_rows, evidence_state);
      bool has_candidates = loadCandidateGroup(candidate_feature_id, candidate_rows, candidate_state);

      while (has_evidence && has_candidates)
      {
        if (evidence_feature_id < candidate_feature_id)
        {
          has_evidence = loadEvidenceGroup(evidence_feature_id, evidence_rows, evidence_state);
          continue;
        }
        if (candidate_feature_id < evidence_feature_id)
        {
          has_candidates = loadCandidateGroup(candidate_feature_id, candidate_rows, candidate_state);
          continue;
        }

        const auto num_it = num_peptidoforms.find(evidence_feature_id);
        if (num_it != num_peptidoforms.end() && num_it->second > 0)
        {
          ++output_feature_count;
          evidence_row_count += evidence_rows.size();
          for (const auto& evidence_row : evidence_rows)
          {
            for (const Int64 peptide_id : candidate_rows)
            {
              const auto key = std::make_pair(evidence_row.transition_id, peptide_id);
              const bool supports_peptidoform = std::binary_search(bitmask_pairs.begin(), bitmask_pairs.end(), key);
              rows.push_back({evidence_feature_id, evidence_row.transition_id, peptide_id, evidence_row.pep,
                              supports_peptidoform ? 1 : 0, num_it->second, std::nullopt});
            }
            if (config.ipf_h0)
            {
              rows.push_back({evidence_feature_id, evidence_row.transition_id, -1, evidence_row.pep, 0, num_it->second, std::nullopt});
            }
          }
        }

        has_evidence = loadEvidenceGroup(evidence_feature_id, evidence_rows, evidence_state);
        has_candidates = loadCandidateGroup(candidate_feature_id, candidate_rows, candidate_state);
      }

      sqlite3_finalize(evidence_stmt);
      sqlite3_finalize(candidate_stmt);

      OPENMS_LOG_INFO << "Read " << evidence_row_count << " transition evidence rows across "
                      << output_feature_count << " eligible features and expanded them to "
                      << rows.size() << " IPF transition hypotheses." << std::endl;
      return rows;
    }

    std::vector<IPFAlignmentRow> OSWFile::readIPFAlignmentData(const PeptidoformInferenceConfig& config) const
    {
      if (!config.propagate_signal_across_runs)
      {
        return {};
      }

      SqliteConnector conn(filename_, SqliteConnector::SqlOpenMode::READWRITE);
      createIPFIndexes_(conn);
      requireTable_(conn, "FEATURE_MS2_ALIGNMENT_CANDIDATE", "Across-run peptidoform propagation requires candidate-based feature alignment results.");

      const String query =
        "SELECT DENSE_RANK() OVER (ORDER BY FEATURE_LIST.PRECURSOR_ID, FEATURE_LIST.ALIGNMENT_ID) AS ALIGNMENT_GROUP_ID, "
        "       FEATURE_LIST.FEATURE_ID "
        "FROM ("
        "  SELECT DISTINCT ALIGNMENT_ID, PRECURSOR_ID, REFERENCE_FEATURE_ID AS FEATURE_ID "
        "  FROM FEATURE_MS2_ALIGNMENT_CANDIDATE "
        "  WHERE SELECTED = 1 AND MAPPING_CONFIDENCE >= ? AND REFERENCE_FEATURE_ID != ALIGNED_FEATURE_ID AND ALIGNED_FEATURE_ID != -1 "
        "  UNION "
        "  SELECT DISTINCT ALIGNMENT_ID, PRECURSOR_ID, ALIGNED_FEATURE_ID AS FEATURE_ID "
        "  FROM FEATURE_MS2_ALIGNMENT_CANDIDATE "
        "  WHERE SELECTED = 1 AND MAPPING_CONFIDENCE >= ? AND REFERENCE_FEATURE_ID != ALIGNED_FEATURE_ID AND ALIGNED_FEATURE_ID != -1"
        ") AS FEATURE_LIST "
        "ORDER BY ALIGNMENT_GROUP_ID, FEATURE_LIST.FEATURE_ID;";

      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, query);
      checkSqliteReturnCode_(conn.getDB(), sqlite3_bind_double(stmt, 1, config.ipf_min_alignment_mapping_confidence),
                             "Failed to bind ipf_min_alignment_mapping_confidence for candidate alignment query");
      checkSqliteReturnCode_(conn.getDB(), sqlite3_bind_double(stmt, 2, config.ipf_min_alignment_mapping_confidence),
                             "Failed to bind ipf_min_alignment_mapping_confidence for candidate alignment query");

      std::vector<IPFAlignmentRow> rows;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        rows.push_back({Sql::extractInt64(stmt, 0), Sql::extractInt64(stmt, 1)});
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      OPENMS_LOG_INFO << "Read " << rows.size() << " candidate-derived alignment-group memberships for IPF propagation." << std::endl;
      return rows;
    }

    std::vector<IPFAlignmentRow> OSWFile::readIPFAlignmentData(double ipf_max_alignment_pep) const
    {
      SqliteConnector conn(filename_, SqliteConnector::SqlOpenMode::READWRITE);
      createIPFIndexes_(conn);

      requireTable_(conn, "FEATURE_MS2_ALIGNMENT", "Across-run peptidoform propagation requires legacy feature alignment results.");
      requireTable_(conn, "SCORE_ALIGNMENT", "Across-run peptidoform propagation requires legacy alignment scores.");

      const String query =
        "SELECT DENSE_RANK() OVER (ORDER BY FEATURE_LIST.PRECURSOR_ID, FEATURE_LIST.ALIGNMENT_ID) AS ALIGNMENT_GROUP_ID, "
        "       FEATURE_LIST.FEATURE_ID "
        "FROM ("
        "  SELECT DISTINCT ALIGNMENT_ID, PRECURSOR_ID, REFERENCE_FEATURE_ID AS FEATURE_ID "
        "  FROM FEATURE_MS2_ALIGNMENT "
        "  WHERE LABEL = 1 AND REFERENCE_FEATURE_ID != ALIGNED_FEATURE_ID "
        "  UNION "
        "  SELECT DISTINCT ALIGNMENT_ID, PRECURSOR_ID, ALIGNED_FEATURE_ID AS FEATURE_ID "
        "  FROM FEATURE_MS2_ALIGNMENT "
        "  WHERE LABEL = 1 AND REFERENCE_FEATURE_ID != ALIGNED_FEATURE_ID"
        ") AS FEATURE_LIST "
        "INNER JOIN ("
        "  SELECT DISTINCT FEATURE_ID "
        "  FROM SCORE_ALIGNMENT "
        "  WHERE PEP < ?"
        ") AS GOOD_ALIGNMENTS ON GOOD_ALIGNMENTS.FEATURE_ID = FEATURE_LIST.FEATURE_ID "
        "ORDER BY ALIGNMENT_GROUP_ID, FEATURE_LIST.FEATURE_ID;";

      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, query);
      checkSqliteReturnCode_(conn.getDB(), sqlite3_bind_double(stmt, 1, ipf_max_alignment_pep),
                             "Failed to bind ipf_max_alignment_pep for historical alignment query");

      std::vector<IPFAlignmentRow> rows;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        rows.push_back({Sql::extractInt64(stmt, 0), Sql::extractInt64(stmt, 1)});
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      OPENMS_LOG_INFO << "Read " << rows.size() << " historical alignment-group memberships for IPF propagation." << std::endl;
      return rows;
    }

    void OSWFile::writeIPFResults(const String& output_filename, const std::vector<IPFResultRow>& results) const
    {
      const String target_filename = prepareOutputFile_(filename_, output_filename);
      SqliteConnector conn(target_filename, SqliteConnector::SqlOpenMode::READWRITE);
      sqlite3* db = conn.getDB();

      conn.executeStatement("DROP TABLE IF EXISTS SCORE_IPF;");
      conn.executeStatement(
        "CREATE TABLE SCORE_IPF ("
        "FEATURE_ID INTEGER NOT NULL, "
        "PEPTIDE_ID INTEGER NOT NULL, "
        "PRECURSOR_PEAKGROUP_PEP DOUBLE NOT NULL, "
        "QVALUE DOUBLE NOT NULL, "
        "PEP DOUBLE NOT NULL);"
      );

      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt,
        "INSERT INTO SCORE_IPF (FEATURE_ID, PEPTIDE_ID, PRECURSOR_PEAKGROUP_PEP, QVALUE, PEP) "
        "VALUES (?, ?, ?, ?, ?);");

      conn.executeStatement("BEGIN TRANSACTION");
      try
      {
        for (const auto& row : results)
        {
          checkSqliteReturnCode_(db, sqlite3_bind_int64(stmt, 1, row.feature_id), "Failed to bind SCORE_IPF.FEATURE_ID");
          checkSqliteReturnCode_(db, sqlite3_bind_int64(stmt, 2, row.peptide_id), "Failed to bind SCORE_IPF.PEPTIDE_ID");
          checkSqliteReturnCode_(db, sqlite3_bind_double(stmt, 3, row.precursor_peakgroup_pep), "Failed to bind SCORE_IPF.PRECURSOR_PEAKGROUP_PEP");
          checkSqliteReturnCode_(db, sqlite3_bind_double(stmt, 4, row.qvalue), "Failed to bind SCORE_IPF.QVALUE");
          checkSqliteReturnCode_(db, sqlite3_bind_double(stmt, 5, row.pep), "Failed to bind SCORE_IPF.PEP");
          executePreparedStatement_(db, stmt, "Failed to insert SCORE_IPF row");
        }
        conn.executeStatement("COMMIT");
      }
      catch (...)
      {
        sqlite3_finalize(stmt);
        try
        {
          conn.executeStatement("ROLLBACK");
        }
        catch (...)
        {
        }
        throw;
      }
      sqlite3_finalize(stmt);

      OPENMS_LOG_DEBUG << "Wrote " << results.size() << " SCORE_IPF rows to '" << target_filename << "'." << std::endl;
    }

    std::vector<LevelContextInputRow> OSWFile::readLevelContextData(InferenceLevel level, InferenceContext context) const
    {
      if (level == InferenceLevel::Peptidoform)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level-context reading only supports peptide, protein, and gene inference.");
      }

      SqliteConnector conn(filename_, SqliteConnector::SqlOpenMode::READWRITE);
      createLevelContextIndexes_(conn, level);

      requireTable_(conn, "PEPTIDE", "Level-context inference requires peptide information.");
      requireTable_(conn, "PRECURSOR_PEPTIDE_MAPPING", "Level-context inference requires precursor-to-peptide mappings.");
      requireTable_(conn, "PRECURSOR", "Level-context inference requires precursor information.");
      requireTable_(conn, "FEATURE", "Level-context inference requires candidate feature information.");
      requireTable_(conn, "SCORE_MS2", "Level-context inference requires MS2 peakgroup scores.");
      if (level == InferenceLevel::Protein)
      {
        requireTable_(conn, "PROTEIN", "Protein inference requires protein information.");
        requireTable_(conn, "PEPTIDE_PROTEIN_MAPPING", "Protein inference requires peptide-to-protein mappings.");
      }
      else if (level == InferenceLevel::Gene)
      {
        requireTable_(conn, "GENE", "Gene inference requires gene information.");
        requireTable_(conn, "PEPTIDE_GENE_MAPPING", "Gene inference requires peptide-to-gene mappings.");
      }

      const bool global_context = context == InferenceContext::Global;
      const String run_select = global_context ? "NULL AS RUN_ID" : "FEATURE.RUN_ID AS RUN_ID";
      const String partition_by = global_context ? "ENTITY_ID" : "RUN_ID, ENTITY_ID";

      String inner_select;
      if (level == InferenceLevel::Peptide)
      {
        inner_select =
          "SELECT " + run_select + ", "
          "       PEPTIDE.ID AS ENTITY_ID, "
          "       PRECURSOR.DECOY AS DECOY, "
          "       SCORE_MS2.SCORE AS SCORE "
          "FROM PEPTIDE "
          "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID "
          "INNER JOIN PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "WHERE SCORE_MS2.SCORE IS NOT NULL";
      }
      else if (level == InferenceLevel::Protein)
      {
        inner_select =
          "SELECT " + run_select + ", "
          "       PROTEIN.ID AS ENTITY_ID, "
          "       PRECURSOR.DECOY AS DECOY, "
          "       SCORE_MS2.SCORE AS SCORE "
          "FROM PROTEIN "
          "INNER JOIN ("
          "  SELECT PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID AS PEPTIDE_ID, PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID "
          "  FROM ("
          "    SELECT PEPTIDE_ID, COUNT(*) AS NUM_PROTEINS "
          "    FROM PEPTIDE_PROTEIN_MAPPING "
          "    GROUP BY PEPTIDE_ID"
          "  ) AS PROTEINS_PER_PEPTIDE "
          "  INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEINS_PER_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID "
          "  WHERE NUM_PROTEINS = 1"
          ") AS UNIQUE_PEPTIDE_PROTEIN_MAPPING ON PROTEIN.ID = UNIQUE_PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID "
          "INNER JOIN PEPTIDE ON UNIQUE_PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID = PEPTIDE.ID "
          "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID "
          "INNER JOIN PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "WHERE SCORE_MS2.SCORE IS NOT NULL";
      }
      else
      {
        inner_select =
          "SELECT " + run_select + ", "
          "       GENE.ID AS ENTITY_ID, "
          "       PRECURSOR.DECOY AS DECOY, "
          "       SCORE_MS2.SCORE AS SCORE "
          "FROM GENE "
          "INNER JOIN ("
          "  SELECT PEPTIDE_GENE_MAPPING.PEPTIDE_ID AS PEPTIDE_ID, PEPTIDE_GENE_MAPPING.GENE_ID "
          "  FROM ("
          "    SELECT PEPTIDE_ID, COUNT(*) AS NUM_GENES "
          "    FROM PEPTIDE_GENE_MAPPING "
          "    GROUP BY PEPTIDE_ID"
          "  ) AS GENES_PER_PEPTIDE "
          "  INNER JOIN PEPTIDE_GENE_MAPPING ON GENES_PER_PEPTIDE.PEPTIDE_ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID "
          "  WHERE NUM_GENES = 1"
          ") AS UNIQUE_PEPTIDE_GENE_MAPPING ON GENE.ID = UNIQUE_PEPTIDE_GENE_MAPPING.GENE_ID "
          "INNER JOIN PEPTIDE ON UNIQUE_PEPTIDE_GENE_MAPPING.PEPTIDE_ID = PEPTIDE.ID "
          "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID "
          "INNER JOIN PRECURSOR ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = PRECURSOR.ID "
          "INNER JOIN FEATURE ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
          "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
          "WHERE SCORE_MS2.SCORE IS NOT NULL";
      }

      const String query =
        "WITH SCORED AS (" + inner_select + "), "
        "RANKED AS ("
        "  SELECT RUN_ID, ENTITY_ID, DECOY, SCORE, "
        "         ROW_NUMBER() OVER (PARTITION BY " + partition_by + " ORDER BY SCORE DESC) AS RN "
        "  FROM SCORED"
        ") "
        "SELECT RUN_ID, ENTITY_ID, DECOY, SCORE "
        "FROM RANKED "
        "WHERE RN = 1 "
        "ORDER BY RUN_ID, ENTITY_ID;";

      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, query);

      std::vector<LevelContextInputRow> rows;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        LevelContextInputRow row;
        if (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
        {
          row.run_id = Sql::extractInt64(stmt, 0);
        }
        row.entity_id = Sql::extractInt64(stmt, 1);
        row.decoy = Sql::extractInt(stmt, 2) != 0;
        row.score = Sql::extractDouble(stmt, 3);
        row.context = context;
        row.group_id = row.run_id.has_value() ? String(*row.run_id) + "_" + String(row.entity_id) : String(row.entity_id);
        rows.push_back(std::move(row));
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      OPENMS_LOG_INFO << "Read " << rows.size() << " best-score rows for "
                      << toString(level) << " inference in '" << toString(context)
                      << "' context." << std::endl;
      return rows;
    }

    void OSWFile::writeLevelContextResults(const String& output_filename,
                                           InferenceLevel level,
                                           InferenceContext context,
                                           const std::vector<LevelContextResultRow>& results) const
    {
      if (level == InferenceLevel::Peptidoform)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level-context writing only supports peptide, protein, and gene inference.");
      }

      const String target_filename = prepareOutputFile_(filename_, output_filename);
      SqliteConnector conn(target_filename, SqliteConnector::SqlOpenMode::READWRITE);
      sqlite3* db = conn.getDB();

      const String table_name = scoreTableName_(level);
      const String entity_column = entityIdColumnName_(level);
      const String context_value = toString(context);
      conn.executeStatement(
        "CREATE TABLE IF NOT EXISTS " + table_name + " ("
        "CONTEXT TEXT NOT NULL, "
        "RUN_ID INTEGER, "
        + entity_column + " INTEGER NOT NULL, "
        "SCORE DOUBLE NOT NULL, "
        "PVALUE DOUBLE NOT NULL, "
        "QVALUE DOUBLE NOT NULL, "
        "PEP DOUBLE NOT NULL);"
      );

      sqlite3_stmt* delete_stmt = nullptr;
      sqlite3_stmt* insert_stmt = nullptr;
      conn.prepareStatement(&delete_stmt, "DELETE FROM " + table_name + " WHERE CONTEXT = ?;");
      conn.prepareStatement(&insert_stmt,
        "INSERT INTO " + table_name + " (CONTEXT, RUN_ID, " + entity_column + ", SCORE, PVALUE, QVALUE, PEP) "
        "VALUES (?, ?, ?, ?, ?, ?, ?);");

      conn.executeStatement("BEGIN TRANSACTION");
      try
      {
        checkSqliteReturnCode_(db, sqlite3_bind_text(delete_stmt, 1, context_value.c_str(), -1, SQLITE_TRANSIENT),
                               "Failed to bind context for SCORE_* delete");
        executePreparedStatement_(db, delete_stmt, "Failed to delete existing level-context rows");

        for (const auto& row : results)
        {
          checkSqliteReturnCode_(db, sqlite3_bind_text(insert_stmt, 1, context_value.c_str(), -1, SQLITE_TRANSIENT),
                                 "Failed to bind SCORE_* CONTEXT");
          bindNullableInt64_(db, insert_stmt, 2, row.run_id);
          checkSqliteReturnCode_(db, sqlite3_bind_int64(insert_stmt, 3, row.entity_id), "Failed to bind SCORE_* entity id");
          checkSqliteReturnCode_(db, sqlite3_bind_double(insert_stmt, 4, row.score), "Failed to bind SCORE_* SCORE");
          checkSqliteReturnCode_(db, sqlite3_bind_double(insert_stmt, 5, row.pvalue), "Failed to bind SCORE_* PVALUE");
          checkSqliteReturnCode_(db, sqlite3_bind_double(insert_stmt, 6, row.qvalue), "Failed to bind SCORE_* QVALUE");
          checkSqliteReturnCode_(db, sqlite3_bind_double(insert_stmt, 7, row.pep), "Failed to bind SCORE_* PEP");
          executePreparedStatement_(db, insert_stmt, "Failed to insert level-context score row");
        }

        conn.executeStatement("COMMIT");
      }
      catch (...)
      {
        sqlite3_finalize(delete_stmt);
        sqlite3_finalize(insert_stmt);
        try
        {
          conn.executeStatement("ROLLBACK");
        }
        catch (...)
        {
        }
        throw;
      }
      sqlite3_finalize(delete_stmt);
      sqlite3_finalize(insert_stmt);

      OPENMS_LOG_DEBUG << "Wrote " << results.size() << " rows to " << table_name
                       << " for context '" << context_value << "' in '" << target_filename << "'." << std::endl;
    }

    void OSWFile::readTransitions_(OSWData& swath_result)
    {
      swath_result.clear();

      Size count = conn_.countTableRows("RUN");
      if (count != 1)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Database '" + filename_ + "' contains more than one RUN. This is currently not supported!");
      }

      // Grab transitions first
      // We do this separately, because the full sql query below will show transitions in duplicates, because many features might use the same XIC at different positions
      const StringList colnames_tr = { "ID", "PRODUCT_MZ", "TYPE", "DECOY", "ANNOTATION" };
      enum COLIDs
      {
        ID,
        PRODUCT_MZ,
        TYPE,
        DECOY,
        ANNOTATION
      };

      // does not make the query below any faster...
      //conn.executeStatement("ANALYZE");

      String select_transitions = "SELECT " + ListUtils::concatenate(colnames_tr, ",") + " FROM TRANSITION ORDER BY ID;";
      sqlite3_stmt* stmt;
      conn_.prepareStatement(&stmt, select_transitions);
      Sql::SqlState rc = Sql::nextRow(stmt);
      while (rc == Sql::SqlState::SQL_ROW)
      {
        OSWTransition tr(Sql::extractString(stmt, COLIDs::ANNOTATION),
          Sql::extractInt(stmt, COLIDs::ID),
          Sql::extractFloat(stmt, COLIDs::PRODUCT_MZ),
          Sql::extractChar(stmt, COLIDs::TYPE),
          Sql::extractInt(stmt, COLIDs::DECOY));
        swath_result.addTransition(std::move(tr));
        rc = Sql::nextRow(stmt);
      }
      sqlite3_finalize(stmt);
    }

} // namespace OpenMS
