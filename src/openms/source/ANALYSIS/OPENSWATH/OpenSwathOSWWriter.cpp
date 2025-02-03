// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
  OpenSwathOSWWriter::OpenSwathOSWWriter(const String& output_filename, const UInt64 run_id, const String& input_filename, bool uis_scores) :
    output_filename_(output_filename),
    input_filename_(input_filename),
    run_id_(Internal::SqliteHelper::clearSignBit(run_id)),
    doWrite_(!output_filename.empty()),
    enable_uis_scoring_(uis_scores)
  {}

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
      "EXP_IM REAL," \
      "DELTA_IM REAL," \
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
      "EXP_IM REAL," \
      "DELTA_IM REAL," \
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
      "VAR_IM_DELTA_SCORE REAL NULL);" \

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
      "APEX_RT REAL NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "RT_FWHM REAL NOT NULL," \
      "MASSERROR_PPM REAL NULL,"
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
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL, " \
      "START_POSITION_AT_5 REAL NULL, " \
      "END_POSITION_AT_5 REAL NULL, " \
      "START_POSITION_AT_10 REAL NULL, " \
      "END_POSITION_AT_10 REAL NULL, " \
      "START_POSITION_AT_50 REAL NULL, " \
      "END_POSITION_AT_50 REAL NULL, " \
      "TOTAL_WIDTH REAL NULL, " \
      "TAILING_FACTOR REAL NULL, " \
      "ASYMMETRY_FACTOR REAL NULL, " \
      "SLOPE_OF_BASELINE REAL NULL, " \
      "BASELINE_DELTA_2_HEIGHT REAL NULL, " \
      "POINTS_ACROSS_BASELINE REAL NULL, " \
      "POINTS_ACROSS_HALF_HEIGHT REAL NULL); ";


    // Execute SQL create statement
    conn.executeStatement(create_sql);

    // Insert run_id information
    std::stringstream sql_run;
    sql_run << "INSERT INTO RUN (ID, FILENAME) VALUES ("
            << run_id_ << ", '"
            << input_filename_ << "'); ";

    // Execute SQL insert statement
    conn.executeStatement(sql_run.str());
  }

  String OpenSwathOSWWriter::getScore(const Feature& feature, const std::string& score_name) const
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

  std::vector<String> OpenSwathOSWWriter::getSeparateScore(const Feature& feature, const std::string& score_name) const
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
                                         const FeatureMap& output,
                                         const String& id) const
  {
    std::stringstream sql, sql_feature, sql_feature_ms1, sql_feature_ms1_precursor, sql_feature_ms2, sql_feature_ms2_transition, sql_feature_uis_transition;

    for (const auto& feature_it : output)
    {
      int64_t feature_id = Internal::SqliteHelper::clearSignBit(feature_it.getUniqueId()); // clear sign bit

      const auto& masserror_ppm = feature_it.metaValueExists("masserror_ppm") ? getSeparateScore(feature_it, "masserror_ppm") : std::vector<String>();

      const auto& subordinates = feature_it.getSubordinates();
      for (Size i=0; i < subordinates.size(); i++)
      {
        const auto& sub_it = subordinates[i];
        if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2")
        {
          std::string total_mi = "NULL"; // total_mi is not guaranteed to be set
          std::string masserror_ppm_query = "NULL"; // masserror_ppm is not guaranteed to be set

          if (!masserror_ppm.empty())
          {
            masserror_ppm_query = masserror_ppm[i];
          }
          if (!sub_it.getMetaValue("total_mi").isEmpty())
          {
            total_mi = sub_it.getMetaValue("total_mi").toString();
          }

          bool enable_compute_peak_shape_metrics = sub_it.metaValueExists("start_position_at_5");
          // Create sql query for storing transition level data, include peak shape metrics if they exist
          sql_feature_ms2_transition << "INSERT INTO FEATURE_TRANSITION "
                         << "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, APEX_INTENSITY, APEX_RT, RT_FWHM, MASSERROR_PPM, TOTAL_MI"
                         << (enable_compute_peak_shape_metrics ? ", START_POSITION_AT_5, END_POSITION_AT_5, "
                                         "START_POSITION_AT_10, END_POSITION_AT_10, START_POSITION_AT_50, END_POSITION_AT_50, "
                                         "TOTAL_WIDTH, TAILING_FACTOR, ASYMMETRY_FACTOR, SLOPE_OF_BASELINE, BASELINE_DELTA_2_HEIGHT, "
                                         "POINTS_ACROSS_BASELINE, POINTS_ACROSS_HALF_HEIGHT" : "")
                         << ") VALUES ("
                         << feature_id << ", "
                         << sub_it.getMetaValue("native_id") << ", "
                         << sub_it.getIntensity() << ", "
                         << sub_it.getMetaValue("total_xic") << ", "
                         << sub_it.getMetaValue("peak_apex_int") << ", "
                         << sub_it.getMetaValue("peak_apex_position") << ", "
                         << sub_it.getMetaValue("width_at_50") << ", "
                         << masserror_ppm_query << ", "
                         << total_mi;

                         if (enable_compute_peak_shape_metrics)
                         {
                            sql_feature_ms2_transition << ", "
                                          << sub_it.getMetaValue("start_position_at_5") << ", "
                                          << sub_it.getMetaValue("end_position_at_5") << ", "
                                          << sub_it.getMetaValue("start_position_at_10") << ", "
                                          << sub_it.getMetaValue("end_position_at_10") << ", "
                                          << sub_it.getMetaValue("start_position_at_50") << ", "
                                          << sub_it.getMetaValue("end_position_at_50") << ", "
                                          << sub_it.getMetaValue("total_width") << ", "
                                          << sub_it.getMetaValue("tailing_factor") << ", "
                                          << sub_it.getMetaValue("asymmetry_factor") << ", "
                                          << sub_it.getMetaValue("slope_of_baseline") << ", "
                                          << sub_it.getMetaValue("baseline_delta_2_height") << ", "
                                          << sub_it.getMetaValue("points_across_baseline") << ", "
                                          << sub_it.getMetaValue("points_across_half_height");
                         }
                         sql_feature_ms2_transition << "); ";
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
                  << run_id_ << ", "
                  << id << ", "
                  << feature_it.getRT() << ", "
                  << getScore(feature_it, "im_drift") << ", "
                  << norm_rt << ", "
                  << delta_rt << ", "
                  << feature_it.getMetaValue("leftWidth") << ", "
                  << feature_it.getMetaValue("rightWidth") << "); ";

      sql_feature_ms2 << "INSERT INTO FEATURE_MS2 " \
        "(FEATURE_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, APEX_INTENSITY, EXP_IM, DELTA_IM, TOTAL_MI, "\
        "VAR_BSERIES_SCORE, VAR_DOTPROD_SCORE, VAR_INTENSITY_SCORE, " \
        "VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE, VAR_LIBRARY_CORR,  "\
        "VAR_LIBRARY_DOTPROD, VAR_LIBRARY_MANHATTAN, VAR_LIBRARY_RMSD, VAR_LIBRARY_ROOTMEANSQUARE, "\
        "VAR_LIBRARY_SANGLE, VAR_LOG_SN_SCORE, VAR_MANHATTAN_SCORE, VAR_MASSDEV_SCORE, VAR_MASSDEV_SCORE_WEIGHTED, "\
        "VAR_MI_SCORE, VAR_MI_WEIGHTED_SCORE, VAR_MI_RATIO_SCORE, VAR_NORM_RT_SCORE, "\
        "VAR_XCORR_COELUTION,VAR_XCORR_COELUTION_WEIGHTED, VAR_XCORR_SHAPE, "\
        "VAR_XCORR_SHAPE_WEIGHTED, VAR_YSERIES_SCORE, VAR_ELUTION_MODEL_FIT_SCORE, "\
        "VAR_IM_XCORR_SHAPE, VAR_IM_XCORR_COELUTION, VAR_IM_DELTA_SCORE"
        << ") VALUES ("
                      << feature_id << ", "
                      << feature_it.getIntensity() << ", "
                      << getScore(feature_it, "total_xic") << ", "
                      << getScore(feature_it, "peak_apices_sum") << ", "
                      << getScore(feature_it, "im_drift") << ", "
                      << getScore(feature_it, "im_delta") << ", "
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
                      << getScore(feature_it, "var_im_delta_score");
      sql_feature_ms2 << "); ";

      bool enable_ms1 = feature_it.metaValueExists("var_ms1_ppm_diff");
      if (enable_ms1) // only write MS1 scores if they are present
      {
        sql_feature_ms1 << "INSERT INTO FEATURE_MS1 "\
          "(FEATURE_ID, AREA_INTENSITY, APEX_INTENSITY, EXP_IM, DELTA_IM, "\
          " VAR_MASSDEV_SCORE, VAR_IM_MS1_DELTA_SCORE, "\
          " VAR_MI_SCORE, VAR_MI_CONTRAST_SCORE, VAR_MI_COMBINED_SCORE, VAR_ISOTOPE_CORRELATION_SCORE, "\
          " VAR_ISOTOPE_OVERLAP_SCORE, VAR_XCORR_COELUTION, VAR_XCORR_COELUTION_CONTRAST, "\
          " VAR_XCORR_COELUTION_COMBINED, VAR_XCORR_SHAPE, VAR_XCORR_SHAPE_CONTRAST, VAR_XCORR_SHAPE_COMBINED "\
          ") VALUES ("
                        << feature_id << ", "
                        << getScore(feature_it, "ms1_area_intensity") << ", "
                        << getScore(feature_it, "ms1_apex_intensity") << ", "
                        << getScore(feature_it, "im_ms1_drift") << ", "
                        << getScore(feature_it, "im_ms1_delta") << ", "
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
        auto id_target_peak_apex_position = getSeparateScore(feature_it, "id_target_peak_apex_position");
        auto id_target_peak_fwhm = getSeparateScore(feature_it, "id_target_width_at_50");
        auto id_target_total_mi = getSeparateScore(feature_it, "id_target_total_mi");
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

        // check if there are compute_peak_shape_metrics scores
        auto id_target_ind_start_position_at_5 = getSeparateScore(feature_it, "id_target_ind_start_position_at_5");
        bool enable_compute_peak_shape_metrics = id_target_ind_start_position_at_5.size() > 0 && id_target_ind_start_position_at_5[0] != "0";
        
        // get scores for peak shape metrics will just be empty vector if not present
        auto start_position_at_5 = getSeparateScore(feature_it, "id_target_ind_start_position_at_5");
        auto end_position_at_5 = getSeparateScore(feature_it, "id_target_ind_end_position_at_5");
        auto start_position_at_10 = getSeparateScore(feature_it, "id_target_ind_start_position_at_10");
        auto end_position_at_10 = getSeparateScore(feature_it, "id_target_ind_end_position_at_10");
        auto start_position_at_50 = getSeparateScore(feature_it, "id_target_ind_start_position_at_50");
        auto end_position_at_50 = getSeparateScore(feature_it, "id_target_ind_end_position_at_50");
        auto total_width = getSeparateScore(feature_it, "id_target_ind_total_width");
        auto tailing_factor = getSeparateScore(feature_it, "id_target_ind_tailing_factor");
        auto asymmetry_factor = getSeparateScore(feature_it, "id_target_ind_asymmetry_factor");
        auto slope_of_baseline = getSeparateScore(feature_it, "id_target_ind_slope_of_baseline");
        auto baseline_delta_2_height = getSeparateScore(feature_it, "id_target_ind_baseline_delta_2_height");
        auto points_across_baseline = getSeparateScore(feature_it, "id_target_ind_points_across_baseline");
        auto points_across_half_height = getSeparateScore(feature_it, "id_target_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_target_num_transitions"))
        {
          int id_target_num_transitions = feature_it.getMetaValue("id_target_num_transitions");

          for (int i = 0; i < id_target_num_transitions; ++i)
          {
            sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION "\
              "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, "\
              " APEX_INTENSITY, APEX_RT, RT_FWHM, MASSERROR_PPM, TOTAL_MI, VAR_INTENSITY_SCORE, VAR_INTENSITY_RATIO_SCORE, "\
              " VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, "\
              " VAR_MASSDEV_SCORE, VAR_MI_SCORE, VAR_MI_RATIO_SCORE, "\
              " VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE "
              << (enable_compute_peak_shape_metrics ? ", START_POSITION_AT_5, END_POSITION_AT_5, "
                                         "START_POSITION_AT_10, END_POSITION_AT_10, START_POSITION_AT_50, END_POSITION_AT_50, "
                                         "TOTAL_WIDTH, TAILING_FACTOR, ASYMMETRY_FACTOR, SLOPE_OF_BASELINE, BASELINE_DELTA_2_HEIGHT, "
                                         "POINTS_ACROSS_BASELINE, POINTS_ACROSS_HALF_HEIGHT" : "")
              << ") VALUES ("
                                        << feature_id << ", "
                                        << id_target_transition_names[i] << ", "
                                        << id_target_area_intensity[i] << ", "
                                        << id_target_total_area_intensity[i] << ", "
                                        << id_target_apex_intensity[i] << ", "
                                        << id_target_peak_apex_position[i] << ", "
                                        << id_target_peak_fwhm[i] << ", "
                                        << id_target_ind_massdev_score[i] << ", "
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
                                        << id_target_ind_isotope_overlap[i];

                         if (enable_compute_peak_shape_metrics)
                         {
                            sql_feature_uis_transition << ", "
                                          << start_position_at_5[i] << ", "
                                          << end_position_at_5[i] << ", "
                                          << start_position_at_10[i] << ", "
                                          << end_position_at_10[i] << ", "
                                          << start_position_at_50[i] << ", "
                                          << end_position_at_50[i] << ", "
                                          << total_width[i] << ", "
                                          << tailing_factor[i] << ", "
                                          << asymmetry_factor[i] << ", "
                                          << slope_of_baseline[i] << ", "
                                          << baseline_delta_2_height[i] << ", "
                                          << points_across_baseline[i] << ", "
                                          << points_across_half_height[i];
                         }
                         sql_feature_uis_transition << "); ";

          }
        }

        auto id_decoy_transition_names = getSeparateScore(feature_it, "id_decoy_transition_names");
        auto id_decoy_area_intensity = getSeparateScore(feature_it, "id_decoy_area_intensity");
        auto id_decoy_total_area_intensity = getSeparateScore(feature_it, "id_decoy_total_area_intensity");
        auto id_decoy_apex_intensity = getSeparateScore(feature_it, "id_decoy_apex_intensity");
        auto id_decoy_peak_apex_position = getSeparateScore(feature_it, "id_decoy_peak_apex_position");
        auto id_decoy_peak_fwhm = getSeparateScore(feature_it, "id_decoy_width_at_50");
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

        // get scores for peak shape metrics will just be empty vector if not present
        auto decoy_start_position_at_5 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_5");
        auto decoy_end_position_at_5 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_5");
        auto decoy_start_position_at_10 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_10");
        auto decoy_end_position_at_10 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_10");
        auto decoy_start_position_at_50 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_50");
        auto decoy_end_position_at_50 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_50");
        auto decoy_total_width = getSeparateScore(feature_it, "id_decoy_ind_total_width");
        auto decoy_tailing_factor = getSeparateScore(feature_it, "id_decoy_ind_tailing_factor");
        auto decoy_asymmetry_factor = getSeparateScore(feature_it, "id_decoy_ind_asymmetry_factor");
        auto decoy_slope_of_baseline = getSeparateScore(feature_it, "id_decoy_ind_slope_of_baseline");
        auto decoy_baseline_delta_2_height = getSeparateScore(feature_it, "id_decoy_ind_baseline_delta_2_height");
        auto decoy_points_across_baseline = getSeparateScore(feature_it, "id_decoy_ind_points_across_baseline");
        auto decoy_points_across_half_height = getSeparateScore(feature_it, "id_decoy_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_decoy_num_transitions"))
        {
          int id_decoy_num_transitions = feature_it.getMetaValue("id_decoy_num_transitions");

          for (int i = 0; i < id_decoy_num_transitions; ++i)
          {
             sql_feature_uis_transition  << "INSERT INTO FEATURE_TRANSITION "\
                "(FEATURE_ID, TRANSITION_ID, AREA_INTENSITY, TOTAL_AREA_INTENSITY, "\
                " APEX_INTENSITY, APEX_RT, RT_FWHM, MASSERROR_PPM, TOTAL_MI, VAR_INTENSITY_SCORE, VAR_INTENSITY_RATIO_SCORE, "\
                " VAR_LOG_INTENSITY, VAR_XCORR_COELUTION, VAR_XCORR_SHAPE, VAR_LOG_SN_SCORE, "\
                " VAR_MASSDEV_SCORE, VAR_MI_SCORE, VAR_MI_RATIO_SCORE, "\
                " VAR_ISOTOPE_CORRELATION_SCORE, VAR_ISOTOPE_OVERLAP_SCORE "
                << (enable_compute_peak_shape_metrics ? ", START_POSITION_AT_5, END_POSITION_AT_5, "
                                         "START_POSITION_AT_10, END_POSITION_AT_10, START_POSITION_AT_50, END_POSITION_AT_50, "
                                         "TOTAL_WIDTH, TAILING_FACTOR, ASYMMETRY_FACTOR, SLOPE_OF_BASELINE, BASELINE_DELTA_2_HEIGHT, "
                                         "POINTS_ACROSS_BASELINE, POINTS_ACROSS_HALF_HEIGHT" : "")
                << ") VALUES ("
                                        << feature_id << ", "
                                        << id_decoy_transition_names[i] << ", "
                                        << id_decoy_area_intensity[i] << ", "
                                        << id_decoy_total_area_intensity[i] << ", "
                                        << id_decoy_apex_intensity[i] << ", "
                                        << id_decoy_peak_apex_position[i] << ", "
                                        << id_decoy_peak_fwhm[i] << ", "
                                        << id_decoy_ind_massdev_score[i] << ", "
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
                                        << id_decoy_ind_isotope_overlap[i];

                         if (enable_compute_peak_shape_metrics)
                         {
                            sql_feature_uis_transition << ", "
                                          << decoy_start_position_at_5[i] << ", "
                                          << decoy_end_position_at_5[i] << ", "
                                          << decoy_start_position_at_10[i] << ", "
                                          << decoy_end_position_at_10[i] << ", "
                                          << decoy_start_position_at_50[i] << ", "
                                          << decoy_end_position_at_50[i] << ", "
                                          << decoy_total_width[i] << ", "
                                          << decoy_tailing_factor[i] << ", "
                                          << decoy_asymmetry_factor[i] << ", "
                                          << decoy_slope_of_baseline[i] << ", "
                                          << decoy_baseline_delta_2_height[i] << ", "
                                          << decoy_points_across_baseline[i] << ", "
                                          << decoy_points_across_half_height[i];
                         }
                         sql_feature_uis_transition << "); ";
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

