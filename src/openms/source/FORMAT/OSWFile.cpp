// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OSWFile.h>

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <sqlite3.h>

#include <cstring> // for strcmp
#include <sstream>
#include <utility> // for std::move

namespace OpenMS
{
  namespace Sql = Internal::SqliteHelper;
  using namespace std;

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
      std::vector<std::basic_string<unsigned char>> group_id_index;
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
            auto it = std::find(group_id_index.begin(), group_id_index.end(), sqlite3_column_text(stmt, i));
            if (it != group_id_index.end())
            {
              scan_id = it - group_id_index.begin();
            }
            else
            {
              scan_id = group_id_index.size();
              group_id_index.emplace_back(sqlite3_column_text(stmt, i));
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
        conn_(filename, SqliteConnector::SqlOpenMode::READONLY)
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
