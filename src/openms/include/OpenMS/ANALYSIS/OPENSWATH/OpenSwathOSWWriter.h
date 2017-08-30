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
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHOSWWRITER_H
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHOSWWRITER_H

// Interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <OpenMS/KERNEL/FeatureMap.h>

#include <sqlite3.h>

#include <fstream>

namespace OpenMS
{

  /**
   * @brief Class to write out an OpenSwath OSW SQLite output (PyProphet input)
   *
   * The class can take a FeatureMap and create a set of string from it
   * suitable for output to OSW using the prepareLine function.
   *
   */
  class OPENMS_DLLAPI OpenSwathOSWWriter
  {
    String output_filename_;
    String input_filename_;
    OpenMS::UInt64 run_id_;
    bool doWrite_;
    bool use_ms1_traces_;
    bool sonar_;
    bool enable_uis_scoring_;

  public:

    OpenSwathOSWWriter(String output_filename, String input_filename = "inputfile", bool ms1_scores = false, bool sonar = false, bool uis_scores = false) :
      output_filename_(output_filename),
      input_filename_(input_filename),
      run_id_(OpenMS::UniqueIdGenerator::getUniqueId()),
      doWrite_(!output_filename.empty()),
      use_ms1_traces_(ms1_scores),
      sonar_(sonar),
      enable_uis_scoring_(uis_scores)
      {}

    static int callback(void * /* NotUsed */, int argc, char **argv, char **azColName){
      int i;
      for(i=0; i<argc; i++)
      {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }
      printf("\n");
      return(0);
    }

    bool isActive() 
    {
      return doWrite_;
    }

    /**
     * @brief Initializes file by generating SQLite tables
     *
     */
    void writeHeader()
    {
      sqlite3 *db;
      char *zErrMsg = 0;
      int  rc;

      // Open database
      rc = sqlite3_open(output_filename_.c_str(), &db);
      if( rc )
      {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      }

      // Create SQL structure
      const char * create_sql =
        "CREATE TABLE RUN(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "FILENAME TEXT NOT NULL); " \

        "CREATE TABLE FEATURE(" \
        "ID INT PRIMARY KEY NOT NULL," \
        "RUN_ID INT NOT NULL," \
        "PRECURSOR_ID INT NOT NULL," \
        "EXP_RT REAL NOT NULL," \
        "NORM_RT REAL NOT NULL," \
        "DELTA_RT REAL NOT NULL," \
        "LEFT_WIDTH REAL NOT NULL," \
        "RIGHT_WIDTH REAL NOT NULL); " \

        "CREATE TABLE FEATURE_MS1(" \
        "FEATURE_ID INT NOT NULL," \
        "AREA_INTENSITY REAL NOT NULL," \
        "APEX_INTENSITY REAL NOT NULL," \
        "VAR_MASSDEV_SCORE REAL NOT NULL," \
        "VAR_ISOTOPE_CORRELATION_SCORE REAL NOT NULL," \
        "VAR_ISOTOPE_OVERLAP_SCORE REAL NOT NULL," \
        "VAR_XCORR_COELUTION REAL NOT NULL," \
        "VAR_XCORR_SHAPE REAL NOT NULL); " \

        "CREATE TABLE FEATURE_MS2(" \
        "FEATURE_ID INT NOT NULL," \
        "AREA_INTENSITY REAL NOT NULL," \
        "APEX_INTENSITY REAL NOT NULL," \
        "VAR_BSERIES_SCORE REAL NOT NULL," \
        "VAR_DOTPROD_SCORE REAL NOT NULL," \
        "VAR_INTENSITY_SCORE REAL NOT NULL," \
        "VAR_ISOTOPE_CORRELATION_SCORE REAL NOT NULL," \
        "VAR_ISOTOPE_OVERLAP_SCORE REAL NOT NULL," \
        "VAR_LIBRARY_CORR REAL NOT NULL," \
        "VAR_LIBRARY_DOTPROD REAL NOT NULL," \
        "VAR_LIBRARY_MANHATTAN REAL NOT NULL," \
        "VAR_LIBRARY_RMSD REAL NOT NULL," \
        "VAR_LIBRARY_ROOTMEANSQUARE REAL NOT NULL," \
        "VAR_LIBRARY_SANGLE REAL NOT NULL," \
        "VAR_LOG_SN_SCORE REAL NOT NULL," \
        "VAR_MANHATTAN_SCORE REAL NOT NULL," \
        "VAR_MASSDEV_SCORE REAL NOT NULL," \
        "VAR_MASSDEV_SCORE_WEIGHTED REAL NOT NULL," \
        "VAR_NORM_RT_SCORE REAL NOT NULL," \
        "VAR_XCORR_COELUTION REAL NOT NULL," \
        "VAR_XCORR_COELUTION_WEIGHTED REAL NOT NULL," \
        "VAR_XCORR_SHAPE REAL NOT NULL," \
        "VAR_XCORR_SHAPE_WEIGHTED REAL NOT NULL," \
        "VAR_YSERIES_SCORE REAL NOT NULL," \
        "VAR_ELUTION_MODEL_FIT_SCORE REAL NULL," \
        "VAR_SONAR_LAG REAL NULL," \
        "VAR_SONAR_SHAPE REAL NULL," \
        "VAR_SONAR_LOG_SN REAL NULL," \
        "VAR_SONAR_LOG_DIFF REAL NULL," \
        "VAR_SONAR_LOG_TREND REAL NULL," \
        "VAR_SONAR_RSQ REAL NULL); " \

        "CREATE TABLE FEATURE_TRANSITION(" \
        "FEATURE_ID INT NOT NULL," \
        "TRANSITION_ID INT NOT NULL," \
        "AREA_INTENSITY REAL NOT NULL," \
        "APEX_INTENSITY REAL NOT NULL," \
        "VAR_LOG_INTENSITY REAL NULL," \
        "VAR_XCORR_COELUTION REAL NULL," \
        "VAR_XCORR_SHAPE REAL NULL," \
        "VAR_LOG_SN_SCORE REAL NULL," \
        "VAR_MASSDEV_SCORE REAL NULL," \
        "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
        "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL); " ;


      // Execute SQL create statement
      rc = sqlite3_exec(db, create_sql, callback, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        std::string error_message = zErrMsg;
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            error_message);
      }

      // Insert run_id information
      std::stringstream sql_run;
      sql_run << "INSERT INTO RUN (ID, FILENAME) VALUES ("
              << *(int64_t*)&run_id_ << ", '" // Conversion from UInt64 to int64_t to support SQLite
              << input_filename_ << "'); ";

      // Execute SQL insert statement
      rc = sqlite3_exec(db, sql_run.str().c_str(), callback, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        std::string error_message = zErrMsg;
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            error_message);
      }

      sqlite3_close(db);
    }

    /**
     * @brief Prepare a single line (feature) for output
     *
     * The result can be flushed to disk using writeLines (either line by line
     * or after collecting several lines).
     *
     * @param pep The compound (peptide/metabolite) used for extraction
     * @param transition The transition used for extraction 
     * @param output The feature map containing all features (each feature will generate one entry in the output)
     * @param id The transition group identifier (peptide/metabolite id)
     *
     * @returns A string to be written using writeLines
     *
     */
    String prepareLine(const OpenSwath::LightCompound& /* pep */,
        const OpenSwath::LightTransition* /* transition */,
        FeatureMap& output, String id)
    {
      std::stringstream sql, sql_feature, sql_feature_ms1, sql_feature_ms2, sql_feature_ms2_transition, sql_feature_uis_transition;

      for (FeatureMap::iterator feature_it = output.begin(); feature_it != output.end(); ++feature_it)
      {
        UInt64 uint64_feature_id = feature_it->getUniqueId();
        int64_t feature_id = *(int64_t*)&uint64_feature_id; // Conversion from UInt64 to int64_t to support SQLite

        for (std::vector<Feature>::iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
        {
          if (sub_it->metaValueExists("FeatureLevel") && sub_it->getMetaValue("FeatureLevel") == "MS2")
          {
            sql_feature_ms2_transition  << "INSERT INTO FEATURE_TRANSITION (FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, APEX_INTENSITY) VALUES (" 
                                        << feature_id << ", " 
                                        << sub_it->getMetaValue("native_id") << ", " 
                                        << sub_it->getIntensity() << ", " 
                                        << sub_it->getMetaValue("peak_apex_int") << "); ";
          }
          else if (sub_it->metaValueExists("FeatureLevel") && sub_it->getMetaValue("FeatureLevel") == "MS1")
          {
            sql_feature_ms1 << "INSERT INTO FEATURE_MS1 (FEATURE_ID, AREA_INTENSITY, APEX_INTENSITY, VAR_MASSDEV_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE) VALUES (" 
                            << feature_id << ", " 
                            << sub_it->getIntensity() << ", " 
                            << sub_it->getMetaValue("peak_apex_int") << ", " 
                            << feature_it->getMetaValue("var_ms1_ppm_diff") << ", " 
                            << feature_it->getMetaValue("var_ms1_isotope_correlation") << ", " 
                            << feature_it->getMetaValue("var_ms1_isotope_overlap") << ", " 
                            << feature_it->getMetaValue("var_ms1_xcorr_coelution") << ", " 
                            << feature_it->getMetaValue("var_ms1_xcorr_shape") << "); ";
          }
        }

        sql_feature << "INSERT INTO FEATURE (ID, RUN_ID, PRECURSOR_ID, EXP_RT, NORM_RT, DELTA_RT, LEFT_WIDTH, RIGHT_WIDTH) VALUES (" 
                    << feature_id << ", '" 
                    << *(int64_t*)&run_id_ << "', " 
                    << id << ", " 
                    << feature_it->getRT() << ", " 
                    << feature_it->getMetaValue("norm_RT") << ", " 
                    << feature_it->getMetaValue("delta_rt") << ", " 
                    << feature_it->getMetaValue("leftWidth") << ", " 
                    << feature_it->getMetaValue("rightWidth") << "); ";

        std::string var_elution_model_fit_score = "NULL", var_sonar_lag = "NULL", var_sonar_shape = "NULL", var_sonar_log_sn = "NULL", var_sonar_log_diff = "NULL", var_sonar_log_trend = "NULL", var_sonar_rsq = "NULL";

        if (!feature_it->getMetaValue("var_elution_model_fit_score").isEmpty())
        {
          var_elution_model_fit_score = feature_it->getMetaValue("var_elution_model_fit_score").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_lag").isEmpty())
        {
          var_sonar_lag = feature_it->getMetaValue("var_sonar_lag").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_shape").isEmpty())
        {
          var_sonar_shape = feature_it->getMetaValue("var_sonar_shape").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_log_sn").isEmpty())
        {
          var_sonar_log_sn = feature_it->getMetaValue("var_sonar_log_sn").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_log_diff").isEmpty())
        {
          var_sonar_log_diff = feature_it->getMetaValue("var_sonar_log_diff").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_log_trend").isEmpty())
        {
          var_sonar_log_trend = feature_it->getMetaValue("var_sonar_log_trend").toString();
        }
        if (!feature_it->getMetaValue("var_sonar_rsq").isEmpty())
        {
          var_sonar_rsq = feature_it->getMetaValue("var_sonar_rsq").toString();
        }

        sql_feature_ms2 << "INSERT INTO FEATURE_MS2 (FEATURE_ID, AREA_INTENSITY, APEX_INTENSITY, VAR_BSERIES_SCORE, VAR_DOTPROD_SCORE, VAR_INTENSITY_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE, VAR_LIBRARY_CORR, VAR_LIBRARY_DOTPROD, VAR_LIBRARY_MANHATTAN, VAR_LIBRARY_RMSD, VAR_LIBRARY_ROOTMEANSQUARE, VAR_LIBRARY_SANGLE, VAR_LOG_SN_SCORE, VAR_MANHATTAN_SCORE, VAR_MASSDEV_SCORE, VAR_MASSDEV_SCORE_WEIGHTED, VAR_NORM_RT_SCORE, VAR_XCORR_COELUTION,VAR_XCORR_COELUTION_WEIGHTED, VAR_XCORR_SHAPE, VAR_XCORR_SHAPE_WEIGHTED, VAR_YSERIES_SCORE, VAR_ELUTION_MODEL_FIT_SCORE, VAR_SONAR_LAG, VAR_SONAR_SHAPE, VAR_SONAR_LOG_SN, VAR_SONAR_LOG_DIFF, VAR_SONAR_LOG_TREND, VAR_SONAR_RSQ) VALUES (" 
                        << feature_id << ", " 
                        << feature_it->getIntensity() << ", " 
                        << feature_it->getMetaValue("peak_apices_sum") << ", " 
                        << feature_it->getMetaValue("var_bseries_score") << ", " 
                        << feature_it->getMetaValue("var_dotprod_score") << ", " 
                        << feature_it->getMetaValue("var_intensity_score") << ", " 
                        << feature_it->getMetaValue("var_isotope_correlation_score") << ", " 
                        << feature_it->getMetaValue("var_isotope_overlap_score") << ", " 
                        << feature_it->getMetaValue("var_library_corr") << ", " 
                        << feature_it->getMetaValue("var_library_dotprod") << ", " 
                        << feature_it->getMetaValue("var_library_manhattan") << ", " 
                        << feature_it->getMetaValue("var_library_rmsd") << ", " 
                        << feature_it->getMetaValue("var_library_rootmeansquare") << ", " 
                        << feature_it->getMetaValue("var_library_sangle") << ", " 
                        << feature_it->getMetaValue("var_log_sn_score") << ", " 
                        << feature_it->getMetaValue("var_manhatt_score") << ", " 
                        << feature_it->getMetaValue("var_massdev_score") << ", " 
                        << feature_it->getMetaValue("var_massdev_score_weighted") << ", " 
                        << feature_it->getMetaValue("var_norm_rt_score") << ", " 
                        << feature_it->getMetaValue("var_xcorr_coelution") << ", " 
                        << feature_it->getMetaValue("var_xcorr_coelution_weighted") << ", " 
                        << feature_it->getMetaValue("var_xcorr_shape") << ", " 
                        << feature_it->getMetaValue("var_xcorr_shape_weighted") << ", " 
                        << feature_it->getMetaValue("var_yseries_score") << ", " 
                        << var_elution_model_fit_score << ", " 
                        << var_sonar_lag << ", "
                        << var_sonar_shape << ", " 
                        << var_sonar_log_sn << ", " 
                        << var_sonar_log_diff << ", " 
                        << var_sonar_log_trend << ", " 
                        << var_sonar_rsq << "); ";

        if (enable_uis_scoring_)
        {
          std::vector<String> id_target_transition_names = ListUtils::create<String>((String)feature_it->getMetaValue("id_target_transition_names"),';');
          std::vector<double> id_target_area_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_area_intensity"),';');
          std::vector<double> id_target_apex_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_apex_intensity"),';');
          std::vector<double> id_target_log_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_log_intensity"),';');
          std::vector<double> id_target_ind_xcorr_coelution = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_xcorr_coelution"),';');
          std::vector<double> id_target_ind_xcorr_shape = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_xcorr_shape"),';');
          std::vector<double> id_target_ind_log_sn_score = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_log_sn_score"),';');
          std::vector<double> id_target_ind_massdev_score = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_massdev_score"),';');
          std::vector<double> id_target_ind_isotope_correlation = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_isotope_correlation"),';');
          std::vector<double> id_target_ind_isotope_overlap = ListUtils::create<double>((String)feature_it->getMetaValue("id_target_ind_isotope_overlap"),';');

          if ((String)feature_it->getMetaValue("id_target_num_transitions") != "")
          {
            for (int i = 0; i < feature_it->getMetaValue("id_target_num_transitions").toString().toInt(); ++i)
            {
              sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION (FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, APEX_INTENSITY, VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, VAR_MASSDEV_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE) VALUES (" 
                                          << feature_id << ", " 
                                          << id_target_transition_names[i] << ", " 
                                          << id_target_area_intensity[i] << ", " 
                                          << id_target_apex_intensity[i] << ", " 
                                          << id_target_log_intensity[i] << ", " 
                                          << id_target_ind_xcorr_coelution[i] << ", " 
                                          << id_target_ind_xcorr_shape[i] << ", " 
                                          << id_target_ind_log_sn_score[i] << ", " 
                                          << id_target_ind_massdev_score[i] << ", " 
                                          << id_target_ind_isotope_correlation[i] << ", " 
                                          << id_target_ind_isotope_overlap[i] << "); ";
            }
          }

          std::vector<String> id_decoy_transition_names = ListUtils::create<String>((String)feature_it->getMetaValue("id_decoy_transition_names"),';');
          std::vector<double> id_decoy_area_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_area_intensity"),';');
          std::vector<double> id_decoy_apex_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_apex_intensity"),';');
          std::vector<double> id_decoy_log_intensity = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_log_intensity"),';');
          std::vector<double> id_decoy_ind_xcorr_coelution = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_xcorr_coelution"),';');
          std::vector<double> id_decoy_ind_xcorr_shape = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_xcorr_shape"),';');
          std::vector<double> id_decoy_ind_log_sn_score = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_log_sn_score"),';');
          std::vector<double> id_decoy_ind_massdev_score = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_massdev_score"),';');
          std::vector<double> id_decoy_ind_isotope_correlation = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_isotope_correlation"),';');
          std::vector<double> id_decoy_ind_isotope_overlap = ListUtils::create<double>((String)feature_it->getMetaValue("id_decoy_ind_isotope_overlap"),';');

          if ((String)feature_it->getMetaValue("id_decoy_num_transitions") != "")
          {
            for (int i = 0; i < feature_it->getMetaValue("id_decoy_num_transitions").toString().toInt(); ++i)
            {
              sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION (FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, APEX_INTENSITY, VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, VAR_MASSDEV_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE) VALUES (" 
                                          << feature_id << ", " 
                                          << id_decoy_transition_names[i] << ", " 
                                          << id_decoy_area_intensity[i] << ", " 
                                          << id_decoy_apex_intensity[i] << ", " 
                                          << id_decoy_log_intensity[i] << ", " 
                                          << id_decoy_ind_xcorr_coelution[i] << ", " 
                                          << id_decoy_ind_xcorr_shape[i] << ", " 
                                          << id_decoy_ind_log_sn_score[i] << ", " 
                                          << id_decoy_ind_massdev_score[i] << ", " 
                                          << id_decoy_ind_isotope_correlation[i] << ", " 
                                          << id_decoy_ind_isotope_overlap[i] << "); ";
            }
          }
        }
      }

      if (enable_uis_scoring_)
      {
        sql << sql_feature.str() << sql_feature_ms1.str() << sql_feature_ms2.str() << sql_feature_uis_transition.str();
      }
      else
      {
        sql << sql_feature.str() << sql_feature_ms1.str() << sql_feature_ms2.str() << sql_feature_ms2_transition.str();
      }

      return(sql.str());
    }

    /**
     * @brief Write data to disk
     *
     * Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
     * 
     * @param to_osw_output Statements generated by prepareLine
     *
     * @note Try to call this function as little as possible (it opens a new
     * database connection each time)
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeLines(std::vector<String> to_osw_output)
    {
      sqlite3 *db;
      char *zErrMsg = 0;
      int  rc;
      // char *create_sql;

      // Open database
      rc = sqlite3_open(output_filename_.c_str(), &db);
      if( rc )
      {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      }

      sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);

      for (Size i = 0; i < to_osw_output.size(); i++)
      {
        rc = sqlite3_exec(db, to_osw_output[i].c_str(), callback, 0, &zErrMsg);
        if( rc != SQLITE_OK )
        {
          std::string error_message = zErrMsg;
          sqlite3_free(zErrMsg);
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              error_message);
        }
      }

      sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);

      sqlite3_close(db);
    }

  };

}

#endif // OPENMS_ANALYSIS_OPENSWATH_OPENSWATHOSWWRITER_H

