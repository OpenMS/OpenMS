// --------------------------------------------------------------------------
//           OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
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
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OSWFile.h>

#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <sqlite3.h>

#include <cstring> // for strcmp
#include <sstream>

#include <iostream> // remove

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
              group_id_index.push_back(sqlite3_column_text(stmt, i));
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

    void OSWFile::read(const String& filename, OSWData& swath_result)
    {
      swath_result.clear();

      // Open database
      SqliteConnector conn(filename);
      Sql::SqlState rc;
      Size count = Sql::countTableRows(conn, "RUN");
      if (count != 1)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Database '" + filename + "' contains more than one RUN. This is currently not supported!");
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

      String select_transitions = "SELECT " + ListUtils::concatenate(colnames_tr, ",") +" FROM TRANSITION ORDER BY ID;";
      sqlite3_stmt* stmt;
      conn.prepareStatement(&stmt, select_transitions);
      rc = Sql::nextRow(stmt);
      while (rc == Sql::SqlState::ROW)
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

      // check of SCORE_MS2 table is available (for OSW files which underwent pyProphet)
      // set q_value to -1 if missing
      bool has_SCORE_MS2 = Sql::countTableRows(conn, "SCORE_MS2") > 0;
      String MS2_select = (has_SCORE_MS2 ? "SCORE_MS2.QVALUE as qvalue" : "-1 as qvalue");
      String MS2_join = (has_SCORE_MS2 ? "inner join(select * from SCORE_MS2) as SCORE_MS2 on SCORE_MS2.FEATURE_ID = FEATURE.ID" : "");

      // assemble the protein-PeptidePrecursor-Feature hierachy
      String select_sql = "select PROTEIN.ID as prot_id, PROTEIN_ACCESSION as prot_accession, PROTEIN.DECOY as decoy, \
                                  PEPTIDE.MODIFIED_SEQUENCE as modified_sequence,\
                                  PRECURSOR.ID as prec_id, PRECURSOR.PRECURSOR_MZ as pc_mz, PRECURSOR.CHARGE as pc_charge,\
                                  FEATURE.ID as feat_id, FEATURE.EXP_RT as rt_expected, FEATURE.DELTA_RT as rt_delta, FEATURE.LEFT_WIDTH as rt_left_width, FEATURE.RIGHT_WIDTH as rt_right_width,\
                                  FeatTrMap.TRANSITION_ID as tr_id," +
                                  MS2_select + "\
        from PROTEIN\
        inner join(select* FROM PEPTIDE_PROTEIN_MAPPING) as PepProtMap on PepProtMap.PROTEIN_ID = PROTEIN.ID\
        inner join(select ID, MODIFIED_SEQUENCE FROM PEPTIDE) as PEPTIDE on PEPTIDE.ID = PepProtMap.PEPTIDE_ID\
        inner join(select * FROM PRECURSOR_PEPTIDE_MAPPING) as PrePepMap on PrePepMap.PEPTIDE_ID = PEPTIDE.ID\
        inner join(select * from PRECURSOR) as PRECURSOR on PRECURSOR.ID = PrePepMap.PRECURSOR_ID\
        inner join(select * from FEATURE) as FEATURE on FEATURE.PRECURSOR_ID = PRECURSOR.ID\
        inner join(select * from FEATURE_TRANSITION) as FeatTrMap on FeatTrMap.FEATURE_ID = FEATURE.ID " +
        MS2_join + "\
        order by prot_id, prec_id, feat_id, qvalue, tr_id";

      std::cout << select_sql << "\n\n";

      conn.prepareStatement(&stmt, select_sql);
      enum CBIG
      { // indices of respective columns in the query above
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
        SIZE_OF_CBIG
      };
      rc = Sql::nextRow(stmt);
      if (sqlite3_column_count(stmt) != SIZE_OF_CBIG)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Query was changed! Please report this bug!");
      }
      
      if (rc == Sql::SqlState::DONE)
      { // no data
        return;
      }

      // Layers of information. Whenever the id changes, we know a new item has begun
      // ... PROTEIN
      int prot_id_old{ Sql::extractInt(stmt, I_PROTID) }, prot_id_new;
      String accession_old = Sql::extractString(stmt, I_ACCESSION), accession_new;
      std::vector<OSWPeptidePrecursor> precursors;
      auto check_add_protein = [&]()
      {
        if (prot_id_old != prot_id_new)
        {
          swath_result.addProtein(OSWProtein(accession_old, precursors));
          prot_id_old = prot_id_new;
          accession_old = accession_new;
          precursors.clear();
        }
      };
      // ... PRECURSOR
      int prec_id_old{ Sql::extractInt(stmt, I_PRECID) }, prec_id_new;
      String seq_old{ Sql::extractString(stmt, I_MODSEQ) }, seq_new;
      short chargePC_old{ (short)Sql::extractInt(stmt, I_PRECQ) }, chargePC_new;
      bool decoy_old{ Sql::extractBool(stmt, I_DECOY) }, decoy_new;
      float precmz_old{ Sql::extractFloat(stmt, I_PRECMZ) }, precmz_new;
      std::vector<OSWPeakGroup> features;
      auto check_add_pc = [&]()
      {
        if (prec_id_old != prec_id_new)
        {
          precursors.push_back(OSWPeptidePrecursor(seq_old, chargePC_old, decoy_old, precmz_old, features));
          prec_id_old = prec_id_new;
          seq_old = seq_new;
          chargePC_old = chargePC_new;
          decoy_old = decoy_new;
          precmz_old = precmz_new;
          features.clear();
          return true;
        }
        return false;
      };
      // ... FEATURE
      Int64 feat_id_old{ Sql::extractInt64(stmt, I_FEATID) }, feat_id_new; // in SQL, feature_id is a 63-bit integer... parsing it would be expensive. So we just use the string.
      float rt_exp_old{ Sql::extractFloat(stmt, I_EXPRT) }, rt_exp_new;
      float rt_lw_old{ Sql::extractFloat(stmt, I_RTLEFT) }, rt_lw_new;
      float rt_rw_old{ Sql::extractFloat(stmt, I_RTRIGHT) }, rt_rw_new;
      float rt_delta_old{ Sql::extractFloat(stmt, I_DELTART) }, rt_delta_new;
      float qvalue_old{ Sql::extractFloat(stmt, I_QVALUE) }, qvalue_new;
      std::vector<UInt32> transition_ids;
      auto check_add_feat = [&]()
      {
        if (feat_id_old != feat_id_new)
        {
          features.push_back(OSWPeakGroup(rt_exp_old, rt_lw_old, rt_rw_old, rt_delta_old, transition_ids, qvalue_old));
          feat_id_old = feat_id_new;
          rt_exp_old = rt_exp_new;
          rt_lw_old = rt_lw_new;
          rt_rw_old = rt_rw_new;
          rt_delta_old = rt_delta_new;
          qvalue_old = qvalue_new;
          transition_ids.clear();
          return true;
        }
        return false;
      };

      // protein loop
      while (rc == Sql::SqlState::ROW)
      {
        prot_id_new = Sql::extractInt(stmt, I_PROTID);
        accession_new = Sql::extractString(stmt, I_ACCESSION);
        decoy_new = Sql::extractBool(stmt, I_DECOY);
        // precursor loop (peptide with charge)
        while (rc == Sql::SqlState::ROW)
        {
          prec_id_new = Sql::extractInt(stmt, I_PRECID);
          seq_new = Sql::extractString(stmt, I_MODSEQ);
          chargePC_new = Sql::extractInt(stmt, I_PRECQ);
          precmz_new = Sql::extractFloat(stmt, I_PRECMZ);
          // feature loop
          while (rc == Sql::SqlState::ROW)
          {
            feat_id_new = Sql::extractInt64(stmt, I_FEATID);
            transition_ids.push_back(Sql::extractInt(stmt, I_TRID));
            rt_exp_new = Sql::extractFloat(stmt, I_EXPRT);
            rt_lw_new = Sql::extractFloat(stmt, I_RTLEFT);
            rt_rw_new = Sql::extractFloat(stmt, I_RTRIGHT);
            rt_delta_new = Sql::extractFloat(stmt, I_DELTART);
            qvalue_new = Sql::extractFloat(stmt, I_QVALUE);
            if (check_add_feat()) break; // feature ended --> check if precursor ended as well.
            rc = Sql::nextRow(stmt); // next row
          }
          if (check_add_pc()) break; // PC ended --> check if protein ended as well.
          rc = Sql::nextRow(stmt); // next row
        }
        check_add_protein();
        rc = Sql::nextRow(stmt); // next row
      }
      check_add_feat();    // add last feature
      check_add_pc();      // add last precursor
      check_add_protein(); // add last protein
      sqlite3_finalize(stmt);
    }

} // namespace OpenMS
