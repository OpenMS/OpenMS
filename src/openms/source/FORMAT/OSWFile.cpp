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

#include <sqlite3.h>

#include <cstring> // for strcmp
#include <sstream>

namespace OpenMS
{
  namespace Sql = Internal::SqliteHelper;
  using namespace std;

  OSWFile::OSWFile()
  {
  }

  OSWFile::~OSWFile()
  {
  }

  void OSWFile::read(const std::string& in_osw,
                     const std::string& osw_level,
                     std::stringstream& pin_output,
                     const double& ipf_max_peakgroup_pep,
                     const double& ipf_max_transition_isotope_overlap,
                     const double& ipf_min_transition_sn)
  {
      sqlite3_stmt * stmt;
      std::string select_sql;

      // Open database
      SqliteConnector conn(in_osw);

      if (osw_level == "ms1")
      {
        select_sql = "SELECT *, RUN_ID || '_' || PRECURSOR.ID AS GROUP_ID " \
                      "FROM FEATURE_MS1 "\
                      "INNER JOIN (SELECT ID, PRECURSOR_ID, RUN_ID FROM FEATURE) AS FEATURE ON FEATURE_ID = FEATURE.ID "\
                      "INNER JOIN (SELECT ID, DECOY FROM PRECURSOR) AS PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID "\
                      "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "\
                      "INNER JOIN (SELECT ID, MODIFIED_SEQUENCE FROM PEPTIDE) AS PEPTIDE ON "\
                        "PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID;";
      }
      else if (osw_level == "transition")
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
        if (osw_level == "transition")
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

    void OSWFile::write(const std::string& in_osw,
                        const std::string& osw_level,
                        const std::map< std::string, std::vector<double> >& features)
    {
      std::string table;
      std::string create_sql;

      if (osw_level == "ms1")
      {
        table = "SCORE_MS1";
        create_sql =  "DROP TABLE IF EXISTS " + table + "; " \
                      "CREATE TABLE " + table + "(" \
                      "FEATURE_ID INT NOT NULL," \
                      "SCORE DOUBLE NOT NULL," \
                      "QVALUE DOUBLE NOT NULL," \
                      "PEP DOUBLE NOT NULL);";

      }
      else if (osw_level == "transition")
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
        if (osw_level == "transition") {
          std::vector<OpenMS::String> ids;
          OpenMS::String(feat.first).split("_", ids);
          insert_sql << "INSERT INTO " << table;
          insert_sql << " (FEATURE_ID, TRANSITION_ID, SCORE, QVALUE, PEP) VALUES (";
          insert_sql <<  ids[0] << ",";
          insert_sql <<  ids[1] << ",";
          insert_sql <<  feat.second[0] << ",";
          insert_sql <<  feat.second[1] << ",";
          insert_sql <<  feat.second[2] << "); ";
        }
        else
        {
          insert_sql << "INSERT INTO " << table;
          insert_sql << " (FEATURE_ID, SCORE, QVALUE, PEP) VALUES (";
          insert_sql <<  feat.first << ",";
          insert_sql <<  feat.second[0] << ",";
          insert_sql <<  feat.second[1] << ",";
          insert_sql <<  feat.second[2] << "); ";
        }

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

} // namespace OpenMS
