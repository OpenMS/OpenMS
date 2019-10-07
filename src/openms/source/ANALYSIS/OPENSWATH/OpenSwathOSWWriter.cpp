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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

#include <OpenMS/FORMAT/SqliteConnector.h>

#include <sqlite3.h>

namespace OpenMS
{

  bool OpenSwathOSWWriter::isActive() const
  {
    return doWrite_;
  }

  void OpenSwathOSWWriter::writeHeader()
  {
    // Open database
    SqliteConnector conn(output_filename_);

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
      "EXP_IM REAL, " \
      "NORM_RT REAL NOT NULL," \
      "DELTA_RT REAL NOT NULL," \
      "LEFT_WIDTH REAL NOT NULL," \
      "RIGHT_WIDTH REAL NOT NULL); " \

      // MS1-level scores
      "CREATE TABLE FEATURE_MS1(" \
      "FEATURE_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_CONTRAST_SCORE REAL NULL," \
      "VAR_MI_COMBINED_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL," \
      "VAR_IM_MS1_DELTA_SCORE REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_COELUTION_CONTRAST REAL NULL," \
      "VAR_XCORR_COELUTION_COMBINED REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_XCORR_SHAPE_CONTRAST REAL NULL," \
      "VAR_XCORR_SHAPE_COMBINED REAL NULL); " \

      // MS2-level scores
      "CREATE TABLE FEATURE_MS2(" \
      "FEATURE_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "TOTAL_AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "TOTAL_MI REAL NULL," \
      "VAR_BSERIES_SCORE REAL NULL," \
      "VAR_DOTPROD_SCORE REAL NULL," \
      "VAR_INTENSITY_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL," \
      "VAR_LIBRARY_CORR REAL NULL," \
      "VAR_LIBRARY_DOTPROD REAL NULL," \
      "VAR_LIBRARY_MANHATTAN REAL NULL," \
      "VAR_LIBRARY_RMSD REAL NULL," \
      "VAR_LIBRARY_ROOTMEANSQUARE REAL NULL," \
      "VAR_LIBRARY_SANGLE REAL NULL," \
      "VAR_LOG_SN_SCORE REAL NULL," \
      "VAR_MANHATTAN_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE_WEIGHTED REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_WEIGHTED_SCORE REAL NULL," \
      "VAR_MI_RATIO_SCORE REAL NULL," \
      "VAR_NORM_RT_SCORE REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_COELUTION_WEIGHTED REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_XCORR_SHAPE_WEIGHTED REAL NULL," \
      "VAR_YSERIES_SCORE REAL NULL," \
      "VAR_ELUTION_MODEL_FIT_SCORE REAL NULL," \

      "VAR_IM_XCORR_SHAPE REAL NULL," \
      "VAR_IM_XCORR_COELUTION REAL NULL," \
      "VAR_IM_DELTA_SCORE REAL NULL," \

      "VAR_SONAR_LAG REAL NULL," \
      "VAR_SONAR_SHAPE REAL NULL," \
      "VAR_SONAR_LOG_SN REAL NULL," \
      "VAR_SONAR_LOG_DIFF REAL NULL," \
      "VAR_SONAR_LOG_TREND REAL NULL," \
      "VAR_SONAR_RSQ REAL NULL); " \

      "CREATE TABLE FEATURE_PRECURSOR(" \
      "FEATURE_ID INT NOT NULL," \
      "ISOTOPE INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL);" \

      // Transition-level scores
      "CREATE TABLE FEATURE_TRANSITION(" \
      "FEATURE_ID INT NOT NULL," \
      "TRANSITION_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "TOTAL_AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "TOTAL_MI REAL NULL," \
      "VAR_INTENSITY_SCORE REAL NULL," \
      "VAR_INTENSITY_RATIO_SCORE REAL NULL," \
      "VAR_LOG_INTENSITY REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_LOG_SN_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_RATIO_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL);";

    // Execute SQL create statement
    conn.executeStatement(create_sql);

    // Insert run_id information
    std::stringstream sql_run;
    sql_run << "INSERT INTO RUN (ID, FILENAME) VALUES ("
            // Conversion from UInt64 to int64_t to support SQLite (and conversion to 63 bits)
            <<  static_cast<int64_t >(run_id_ & ~(1ULL << 63)) << ", '"
            << input_filename_ << "'); ";

    // Execute SQL insert statement
    conn.executeStatement(sql_run);
  }

  String OpenSwathOSWWriter::getScore(const Feature& feature, std::string score_name) const
  {
    String score = "NULL";
    if (!feature.getMetaValue(score_name).isEmpty())
    {
      score = feature.getMetaValue(score_name).toString();
    }
    if (score.toLower() == "nan") score = "NULL";
    if (score.toLower() == "-nan") score = "NULL";

    return score;
  }

  std::vector<String> OpenSwathOSWWriter::getSeparateScore(const Feature& feature, std::string score_name) const
  {
    std::vector<String> separated_scores;

    if (!feature.getMetaValue(score_name).isEmpty())
    {
      if (feature.getMetaValue(score_name).valueType() == DataValue::STRING_LIST)
      {
        separated_scores = feature.getMetaValue(score_name).toStringList();
      }
      else if (feature.getMetaValue(score_name).valueType() == DataValue::INT_LIST)
      {
        std::vector<int> int_separated_scores = feature.getMetaValue(score_name).toIntList();
        std::transform(int_separated_scores.begin(), int_separated_scores.end(), std::back_inserter(separated_scores), [](const int& num) { return String(num); });

      }
      else if (feature.getMetaValue(score_name).valueType() == DataValue::DOUBLE_LIST)
      {
        std::vector<double> double_separated_scores = feature.getMetaValue(score_name).toDoubleList();
        std::transform(double_separated_scores.begin(), double_separated_scores.end(), std::back_inserter(separated_scores), [](const double& num) { return String(num); });
      }
      else
      {
        separated_scores.push_back(feature.getMetaValue(score_name).toString());
      }
    }

    return separated_scores;
  }

  String OpenSwathOSWWriter::prepareLine(const OpenSwath::LightCompound& /* pep */,
                                         const OpenSwath::LightTransition* /* transition */,
                                         FeatureMap& output,
                                         String id) const
  {
    std::stringstream sql, sql_feature, sql_feature_ms1, sql_feature_ms1_precursor, sql_feature_ms2, sql_feature_ms2_transition, sql_feature_uis_transition;

    for (const auto& feature_it : output)
    {
      UInt64 uint64_feature_id = feature_it.getUniqueId();
      int64_t feature_id = static_cast<int64_t >(uint64_feature_id & ~(1ULL << 63)); // clear sign bit

      for (const auto& sub_it : feature_it.getSubordinates())
      {
        if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2")
        {
          std::string total_mi = "NULL"; // total_mi is not guaranteed to be set
          if (!sub_it.getMetaValue("total_mi").isEmpty())
          {
            total_mi = sub_it.getMetaValue("total_mi").toString();
          }
          sql_feature_ms2_transition  << "INSERT INTO FEATURE_TRANSITION "\
            "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, APEX_INTENSITY, TOTAL_MI) VALUES ("
                                      << feature_id << ", "
                                      << sub_it.getMetaValue("native_id") << ", "
                                      << sub_it.getIntensity() << ", "
                                      << sub_it.getMetaValue("total_xic") << ", "
                                      << sub_it.getMetaValue("peak_apex_int") << ", "
                                      << total_mi << "); ";
        }
        else if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS1" && sub_it.getIntensity() > 0.0)
        {
          std::vector<String> precursor_id;
          OpenMS::String(sub_it.getMetaValue("native_id")).split(OpenMS::String("Precursor_i"), precursor_id);
          sql_feature_ms1_precursor  << "INSERT INTO FEATURE_PRECURSOR (FEATURE_ID, ISOTOPE, AREA_INTENSITY, APEX_INTENSITY) VALUES ("
                                      << feature_id << ", "
                                      << precursor_id[1] << ", "
                                      << sub_it.getIntensity() << ", "
                                      << sub_it.getMetaValue("peak_apex_int") << "); ";
        }
      }

      // these will be missing if RT scoring is disabled
      double norm_rt = -1, delta_rt = -1;
      if (feature_it.metaValueExists("norm_RT") ) norm_rt = feature_it.getMetaValue("norm_RT");
      if (feature_it.metaValueExists("delta_rt") ) delta_rt = feature_it.getMetaValue("delta_rt");

      sql_feature << "INSERT INTO FEATURE (ID, RUN_ID, PRECURSOR_ID, EXP_RT, EXP_IM, NORM_RT, DELTA_RT, LEFT_WIDTH, RIGHT_WIDTH) VALUES ("
                  << feature_id << ", "
                  // Conversion from UInt64 to int64_t to support SQLite (and conversion to 63 bits)
                  <<  static_cast<int64_t >(run_id_ & ~(1ULL << 63)) << ", "
                  << id << ", "
                  << feature_it.getRT() << ", "
                  << getScore(feature_it, "im_drift") << ", "
                  << norm_rt << ", "
                  << delta_rt << ", "
                  << feature_it.getMetaValue("leftWidth") << ", "
                  << feature_it.getMetaValue("rightWidth") << "); ";

      sql_feature_ms2 << "INSERT INTO FEATURE_MS2 " \
        "(FEATURE_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, APEX_INTENSITY, TOTAL_MI, "\
        "VAR_BSERIES_SCORE, VAR_DOTPROD_SCORE, VAR_INTENSITY_SCORE, " \
        "VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE, VAR_LIBRARY_CORR,  "\
        "VAR_LIBRARY_DOTPROD, VAR_LIBRARY_MANHATTAN, VAR_LIBRARY_RMSD, VAR_LIBRARY_ROOTMEANSQUARE, "\
        "VAR_LIBRARY_SANGLE, VAR_LOG_SN_SCORE, VAR_MANHATTAN_SCORE, VAR_MASSDEV_SCORE, VAR_MASSDEV_SCORE_WEIGHTED, "\
        "VAR_MI_SCORE, VAR_MI_WEIGHTED_SCORE, VAR_MI_RATIO_SCORE, VAR_NORM_RT_SCORE, "\
        "VAR_XCORR_COELUTION,VAR_XCORR_COELUTION_WEIGHTED, VAR_XCORR_SHAPE, "\
        "VAR_XCORR_SHAPE_WEIGHTED, VAR_YSERIES_SCORE, VAR_ELUTION_MODEL_FIT_SCORE, "\
        "VAR_IM_XCORR_SHAPE, VAR_IM_XCORR_COELUTION, VAR_IM_DELTA_SCORE, " \
        "VAR_SONAR_LAG, VAR_SONAR_SHAPE, VAR_SONAR_LOG_SN, VAR_SONAR_LOG_DIFF, VAR_SONAR_LOG_TREND, VAR_SONAR_RSQ "\
        ") VALUES ("
                      << feature_id << ", "
                      << feature_it.getIntensity() << ", "
                      << getScore(feature_it, "total_xic") << ", "
                      << getScore(feature_it, "peak_apices_sum") << ", "
                      << getScore(feature_it, "total_mi") << ", "
                      << getScore(feature_it, "var_bseries_score") << ", "
                      << getScore(feature_it, "var_dotprod_score") << ", "
                      << getScore(feature_it, "var_intensity_score") << ", "
                      << getScore(feature_it, "var_isotope_correlation_score") << ", "
                      << getScore(feature_it, "var_isotope_overlap_score") << ", "
                      << getScore(feature_it, "var_library_corr") << ", "
                      << getScore(feature_it, "var_library_dotprod") << ", "
                      << getScore(feature_it, "var_library_manhattan") << ", "
                      << getScore(feature_it, "var_library_rmsd") << ", "
                      << getScore(feature_it, "var_library_rootmeansquare") << ", "
                      << getScore(feature_it, "var_library_sangle") << ", "
                      << getScore(feature_it, "var_log_sn_score") << ", "
                      << getScore(feature_it, "var_manhatt_score") << ", "
                      << getScore(feature_it, "var_massdev_score") << ", "
                      << getScore(feature_it, "var_massdev_score_weighted") << ", "
                      << getScore(feature_it, "var_mi_score") << ", "
                      << getScore(feature_it, "var_mi_weighted_score") << ", "
                      << getScore(feature_it, "var_mi_ratio_score") << ", "
                      << getScore(feature_it, "var_norm_rt_score") << ", "
                      << getScore(feature_it, "var_xcorr_coelution") << ", "
                      << getScore(feature_it, "var_xcorr_coelution_weighted") << ", "
                      << getScore(feature_it, "var_xcorr_shape") << ", "
                      << getScore(feature_it, "var_xcorr_shape_weighted") << ", "
                      << getScore(feature_it, "var_yseries_score") << ", "
                      << getScore(feature_it, "var_elution_model_fit_score") << ", "
                      << getScore(feature_it, "var_im_xcorr_shape") << ", "
                      << getScore(feature_it, "var_im_xcorr_coelution") << ", "
                      << getScore(feature_it, "var_im_delta_score") << ", "
                      << getScore(feature_it, "var_sonar_lag") << ", "
                      << getScore(feature_it, "var_sonar_shape") << ", "
                      << getScore(feature_it, "var_sonar_log_sn") << ", "
                      << getScore(feature_it, "var_sonar_log_diff") << ", "
                      << getScore(feature_it, "var_sonar_log_trend") << ", "
                      << getScore(feature_it, "var_sonar_rsq") << "); ";

      if (use_ms1_traces_)
      {
        sql_feature_ms1 << "INSERT INTO FEATURE_MS1 "\
          "(FEATURE_ID, AREA_INTENSITY, APEX_INTENSITY, "\
          " VAR_MASSDEV_SCORE, VAR_IM_MS1_DELTA_SCORE, "\
          " VAR_MI_SCORE, VAR_MI_CONTRAST_SCORE, VAR_MI_COMBINED_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, "\
          " VAR_ISOTOPE_OVERLAP_SCORE, VAR_XCORR_COELUTION, VAR_XCORR_COELUTION_CONTRAST, "\
          " VAR_XCORR_COELUTION_COMBINED, VAR_XCORR_SHAPE, VAR_XCORR_SHAPE_CONTRAST, VAR_XCORR_SHAPE_COMBINED "\
          ") VALUES ("
                        << feature_id << ", "
                        << getScore(feature_it, "ms1_area_intensity") << ", "
                        << getScore(feature_it, "ms1_apex_intensity") << ", "
                        << getScore(feature_it, "var_ms1_ppm_diff") << ", "
                        << getScore(feature_it, "var_im_ms1_delta_score") << ", "
                        << getScore(feature_it, "var_ms1_mi_score") << ", "
                        << getScore(feature_it, "var_ms1_mi_contrast_score") << ", "
                        << getScore(feature_it, "var_ms1_mi_combined_score") << ", "
                        << getScore(feature_it, "var_ms1_isotope_correlation") << ", "
                        << getScore(feature_it, "var_ms1_isotope_overlap") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_coelution") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_coelution_contrast") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_coelution_combined") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_shape") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_shape_contrast") << ", "
                        << getScore(feature_it, "var_ms1_xcorr_shape_combined") << "); ";
      }

      if (enable_uis_scoring_)
      {
        auto id_target_transition_names = getSeparateScore(feature_it, "id_target_transition_names");
        auto id_target_area_intensity = getSeparateScore(feature_it, "id_target_area_intensity");
        auto id_target_total_area_intensity = getSeparateScore(feature_it, "id_target_total_area_intensity");
        auto id_target_apex_intensity = getSeparateScore(feature_it, "id_target_apex_intensity");
        auto id_target_total_mi = getSeparateScore(feature_it, "id_target_apex_intensity");
        auto id_target_intensity_score = getSeparateScore(feature_it, "id_target_intensity_score");
        auto id_target_intensity_ratio_score = getSeparateScore(feature_it, "id_target_intensity_ratio_score");
        auto id_target_log_intensity = getSeparateScore(feature_it, "id_target_ind_log_intensity");
        auto id_target_ind_xcorr_coelution = getSeparateScore(feature_it, "id_target_ind_xcorr_coelution");
        auto id_target_ind_xcorr_shape = getSeparateScore(feature_it, "id_target_ind_xcorr_shape");
        auto id_target_ind_log_sn_score = getSeparateScore(feature_it, "id_target_ind_log_sn_score");
        auto id_target_ind_massdev_score = getSeparateScore(feature_it, "id_target_ind_massdev_score");
        auto id_target_ind_mi_score = getSeparateScore(feature_it, "id_target_ind_mi_score");
        auto id_target_ind_mi_ratio_score = getSeparateScore(feature_it, "id_target_ind_mi_ratio_score");
        auto id_target_ind_isotope_correlation = getSeparateScore(feature_it, "id_target_ind_isotope_correlation");
        auto id_target_ind_isotope_overlap = getSeparateScore(feature_it, "id_target_ind_isotope_overlap");

        if (feature_it.metaValueExists("id_target_num_transitions"))
        {
          int id_target_num_transitions = feature_it.getMetaValue("id_target_num_transitions");

          for (int i = 0; i < id_target_num_transitions; ++i)
          {
            sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION "\
              "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, "\
              " APEX_INTENSITY, TOTAL_MI, VAR_INTENSITY_SCORE, VAR_INTENSITY_RATIO_SCORE, "\
              " VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, "\
              " VAR_MASSDEV_SCORE, VAR_MI_SCORE, VAR_MI_RATIO_SCORE, "\
              " VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE "\
              ") VALUES ("
                                        << feature_id << ", "
                                        << id_target_transition_names[i] << ", "
                                        << id_target_area_intensity[i] << ", "
                                        << id_target_total_area_intensity[i] << ", "
                                        << id_target_apex_intensity[i] << ", "
                                        << id_target_total_mi[i] << ", "
                                        << id_target_intensity_score[i] << ", "
                                        << id_target_intensity_ratio_score[i] << ", "
                                        << id_target_log_intensity[i] << ", "
                                        << id_target_ind_xcorr_coelution[i] << ", "
                                        << id_target_ind_xcorr_shape[i] << ", "
                                        << id_target_ind_log_sn_score[i] << ", "
                                        << id_target_ind_massdev_score[i] << ", "
                                        << id_target_ind_mi_score[i] << ", "
                                        << id_target_ind_mi_ratio_score[i] << ", "
                                        << id_target_ind_isotope_correlation[i] << ", "
                                        << id_target_ind_isotope_overlap[i] << "); ";
          }
        }

        auto id_decoy_transition_names = getSeparateScore(feature_it, "id_decoy_transition_names");
        auto id_decoy_area_intensity = getSeparateScore(feature_it, "id_decoy_area_intensity");
        auto id_decoy_total_area_intensity = getSeparateScore(feature_it, "id_decoy_total_area_intensity");
        auto id_decoy_apex_intensity = getSeparateScore(feature_it, "id_decoy_apex_intensity");
        auto id_decoy_total_mi = getSeparateScore(feature_it, "id_decoy_total_mi");
        auto id_decoy_intensity_score = getSeparateScore(feature_it, "id_decoy_intensity_score");
        auto id_decoy_intensity_ratio_score = getSeparateScore(feature_it, "id_decoy_intensity_ratio_score");
        auto id_decoy_log_intensity = getSeparateScore(feature_it, "id_decoy_ind_log_intensity");
        auto id_decoy_ind_xcorr_coelution = getSeparateScore(feature_it, "id_decoy_ind_xcorr_coelution");
        auto id_decoy_ind_xcorr_shape = getSeparateScore(feature_it, "id_decoy_ind_xcorr_shape");
        auto id_decoy_ind_log_sn_score = getSeparateScore(feature_it, "id_decoy_ind_log_sn_score");
        auto id_decoy_ind_massdev_score = getSeparateScore(feature_it, "id_decoy_ind_massdev_score");
        auto id_decoy_ind_mi_score = getSeparateScore(feature_it, "id_decoy_ind_mi_score");
        auto id_decoy_ind_mi_ratio_score = getSeparateScore(feature_it, "id_decoy_ind_mi_ratio_score");
        auto id_decoy_ind_isotope_correlation = getSeparateScore(feature_it, "id_decoy_ind_isotope_correlation");
        auto id_decoy_ind_isotope_overlap = getSeparateScore(feature_it, "id_decoy_ind_isotope_overlap");

        if (feature_it.metaValueExists("id_decoy_num_transitions"))
        {
          int id_decoy_num_transitions = feature_it.getMetaValue("id_decoy_num_transitions");

          for (int i = 0; i < id_decoy_num_transitions; ++i)
          {
             sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION "\
                "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, "\
                " APEX_INTENSITY, TOTAL_MI, VAR_INTENSITY_SCORE, VAR_INTENSITY_RATIO_SCORE, "\
                " VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, "\
                " VAR_MASSDEV_SCORE, VAR_MI_SCORE, VAR_MI_RATIO_SCORE, "\
                " VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE) "\
                "VALUES ("
                                        << feature_id << ", "
                                        << id_decoy_transition_names[i] << ", "
                                        << id_decoy_area_intensity[i] << ", "
                                        << id_decoy_total_area_intensity[i] << ", "
                                        << id_decoy_apex_intensity[i] << ", "
                                        << id_decoy_total_mi[i] << ", "
                                        << id_decoy_intensity_score[i] << ", "
                                        << id_decoy_intensity_ratio_score[i] << ", "
                                        << id_decoy_log_intensity[i] << ", "
                                        << id_decoy_ind_xcorr_coelution[i] << ", "
                                        << id_decoy_ind_xcorr_shape[i] << ", "
                                        << id_decoy_ind_log_sn_score[i] << ", "
                                        << id_decoy_ind_massdev_score[i] << ", "
                                        << id_decoy_ind_mi_score[i] << ", "
                                        << id_decoy_ind_mi_ratio_score[i] << ", "
                                        << id_decoy_ind_isotope_correlation[i] << ", "
                                        << id_decoy_ind_isotope_overlap[i] << "); ";
          }
        }
      }
    }

    if (enable_uis_scoring_ && !sql_feature_uis_transition.str().empty() )
    {
      sql << sql_feature.str() << sql_feature_ms1.str() << sql_feature_ms1_precursor.str() << sql_feature_ms2.str() << sql_feature_uis_transition.str();
    }
    else
    {
      sql << sql_feature.str() << sql_feature_ms1.str() << sql_feature_ms1_precursor.str() << sql_feature_ms2.str() << sql_feature_ms2_transition.str();
    }

    return sql.str();
  }

  void OpenSwathOSWWriter::writeLines(const std::vector<String>& to_osw_output)
  {
    SqliteConnector conn(output_filename_);
    conn.executeStatement("BEGIN TRANSACTION");
    for (Size i = 0; i < to_osw_output.size(); i++)
    {
      conn.executeStatement(to_osw_output[i]);
    }
    conn.executeStatement("END TRANSACTION");
  }
}

