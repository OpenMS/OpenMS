// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ML/REGRESSION/LinearRegression.h>
#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <fstream>
#include <algorithm>

#define SWATHMAPMASSCORRECTION_DEBUG

namespace OpenMS
{
  void findBestFeature(const OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType& transition_group, double& bestRT)
  {
    // Find the feature with the highest score
    bestRT = -1;
    double highest_score = -1000;
    for (const auto& mrmfeature : transition_group.getFeatures())
    {
      if (mrmfeature.getOverallQuality() > highest_score)
      {
        bestRT = mrmfeature.getRT();
        highest_score = mrmfeature.getOverallQuality();
      }
    }
  }

  std::vector<OpenSwath::SwathMap> findSwathMaps(const OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType& transition_group,
                                                 const std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    // Get the corresponding SWATH map
    std::vector<OpenSwath::SwathMap> used_maps;
    for (const auto& m : swath_maps)
    {
      if (m.lower < transition_group.getTransitions()[0].precursor_mz &&
          m.upper >= transition_group.getTransitions()[0].precursor_mz)
      {
        used_maps.push_back(m);
      }
    }
    return used_maps;
  }

  // Computes the SwathMaps for PASEF data in which windows can have the same m/z but differ by ion mobility
  // NOTE: swathMap is stored as a vector to enable compatibility with downstream function calls (as SONAR) can have multiple windows
  std::vector<OpenSwath::SwathMap> SwathMapMassCorrection::findSwathMapsPasef(const OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType& transition_group,
                                                 const std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    OPENMS_PRECONDITION(transition_group.getTransitions()[0].precursor_im != -1, "All transitions must have a valid IM value (not -1)");
    // Although theoretically there can be more than one map, for this case, just use the "best" map, best map is defined as the one in which the IM is closest to the center of the window
    std::vector<OpenSwath::SwathMap> used_maps;
    for (const auto& m : swath_maps)
    {
      // If precursor m/z and IM in Swath window
      if (m.lower < transition_group.getTransitions()[0].precursor_mz &&
          m.upper >= transition_group.getTransitions()[0].precursor_mz &&
          m.imLower < transition_group.getTransitions()[0].precursor_im &&
          m.imUpper >= transition_group.getTransitions()[0].precursor_im)
      {
        // if no other windows at this position just add it
        if (used_maps.size() == 0)
        {
          used_maps.push_back(m);
        }
        else //there is another window at this position, check if the new window found is better
        {
          double imCenterDiffOld = std::fabs(((used_maps[0].imLower + used_maps[0].imUpper) / 2) - transition_group.getTransitions()[0].precursor_im);
          double imCenterDiffNew = std::fabs(((m.imLower + m.imUpper) / 2) - transition_group.getTransitions()[0].precursor_im);
          if (imCenterDiffOld > imCenterDiffNew)
          {
            used_maps[0] = m;
          }
        }
      }
    }
    return used_maps;
  }

  SwathMapMassCorrection::SwathMapMassCorrection() :
    DefaultParamHandler("SwathMapMassCorrection")
  {
    defaults_.setValue("mz_extraction_window", -1.0, "M/z extraction window width");
    defaults_.setValue("mz_extraction_window_ppm", "false", "Whether m/z extraction is in ppm", {"advanced"});
    defaults_.setValidStrings("mz_extraction_window_ppm", {"true","false"});
    defaults_.setValue("ms1_im_calibration", "false", "Whether to use MS1 precursor data for the ion mobility calibration (default = false, uses MS2 / fragment ions for calibration)", {"advanced"});
    defaults_.setValidStrings("ms1_im_calibration", {"true","false"});
    defaults_.setValue("im_extraction_window", -1.0, "Ion mobility extraction window width");
    defaults_.setValue("mz_estimation_padding_factor", 1.3, "A padding factor to multiply the estimated m/z window by. For example, a factor of 1.3 will add a 30% padding to the estimated m/z window, so if the estimated m/z window is 18, then 5.4 will be added for a total estimated m/z window of 23.4. A factor of 1.0 will not add any padding to the estimated window.");
    defaults_.setMinFloat("mz_estimation_padding_factor", 1.0);
    defaults_.setValue("im_estimation_padding_factor", 1.3, "A padding factor to multiply the estimated ion_mobility window by. For example, a factor of 1.3 will add a 30% padding to the estimated ion_mobility window, so if the estimated ion_mobility window is 0.03, then 0.009 will be added for a total estimated ion_mobility window of 0.039. A factor of 1.0 will not add any padding to the estimated window.");
    defaults_.setMinFloat("im_estimation_padding_factor", 1.0);
    defaults_.setValue("mz_correction_function", "none", "Type of normalization function for m/z calibration.");
    defaults_.setValidStrings("mz_correction_function", {"none","regression_delta_ppm","unweighted_regression","weighted_regression","quadratic_regression","weighted_quadratic_regression","weighted_quadratic_regression_delta_ppm","quadratic_regression_delta_ppm"});
    defaults_.setValue("im_correction_function", "linear", "Type of normalization function for IM calibration.");
    defaults_.setValidStrings("im_correction_function", {"none","linear"});

    defaults_.setValue("debug_im_file", "", "Debug file for Ion Mobility calibration.");
    defaults_.setValue("debug_mz_file", "", "Debug file for m/z calibration.");

    // write defaults into Param object param_
    defaultsToParam_();
  }

  void SwathMapMassCorrection::updateMembers_()
  {
    mz_extraction_window_ = (double)param_.getValue("mz_extraction_window");
    mz_extraction_window_ppm_ = param_.getValue("mz_extraction_window_ppm") == "true";
    ms1_im_ = param_.getValue("ms1_im_calibration") == "true";
    im_extraction_window_ = (double)param_.getValue("im_extraction_window");
    mz_estimation_padding_factor_ = (double)param_.getValue("mz_estimation_padding_factor");
    im_estimation_padding_factor_ = (double)param_.getValue("im_estimation_padding_factor");
    mz_correction_function_ = param_.getValue("mz_correction_function").toString();
    im_correction_function_ = param_.getValue("im_correction_function").toString();
    debug_mz_file_ = param_.getValue("debug_mz_file").toString();
    debug_im_file_ = param_.getValue("debug_im_file").toString();
  }

  void SwathMapMassCorrection::correctIM(
    const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const bool pasef,
    TransformationDescription& im_trafo)
  {
    bool ppm = mz_extraction_window_ppm_;
    double mz_extr_window = mz_extraction_window_;
    double im_extraction_win = im_extraction_window_;
    double im_estimation_padding_factor = im_estimation_padding_factor_;

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctIM " << " window " << im_extraction_win << " mz window " << mz_extr_window << " in ppm " << ppm << std::endl;

    if (im_extraction_win < 0)
    {
      return;
    }

    if (im_correction_function_ == "none")
    {
      return;
    }
    // if it is not none, then it must be linear

    std::ofstream os_im;
    if (!debug_im_file_.empty())
    {
      std::cout.precision(16);
      os_im.open(debug_im_file_);
      os_im << "mz" << "\t" << "im" << "\t" << "theo_im" << "\t" << "RT" << "\t" << "intensity" << std::endl;
      os_im.precision(writtenDigits(double()));
    }

    std::vector<String> trgr_ids;
    std::map<std::string, double> pep_im_map;
    trgr_ids.reserve(transition_group_map.size());
    for (const auto& trgroup_it : transition_group_map)
    {
      trgr_ids.push_back(trgroup_it.first);
    }
    for (const auto& cmp : targeted_exp.getCompounds())
    {
      pep_im_map[cmp.id] = cmp.drift_time;
    }

    TransformationDescription::DataPoints data_im;
    std::vector<double> exp_im;
    std::vector<double> theo_im;

    // Collect MS1 IM pairs across all transition groups
    std::vector<double> exp_im_ms1_all, theo_im_ms1_all;

    #pragma omp parallel for
    for (SignedSize k = 0; k < (SignedSize)trgr_ids.size(); k++)
    {
      // we need at least one feature to find the best one
      auto transition_group = transition_group_map.at(trgr_ids[k]);
      if (transition_group->getFeatures().empty()) continue;

      // Find the feature with the highest score
      double bestRT;
      findBestFeature(*transition_group, bestRT);
      // Get the corresponding SWATH map(s), for SONAR there will be more than one map
      std::vector<OpenSwath::SwathMap> used_maps;
      if (!pasef)
      {
        used_maps = findSwathMaps(*transition_group, swath_maps);
      }
      // If pasef then have to check for overlap across IM
      else
      {
        used_maps = findSwathMapsPasef(*transition_group, swath_maps);
      }

      // We will collect MS1 (exp, theo) points regardless of whether ms1_im_calibration is used for fitting
      std::vector<double> exp_im_ms1_local, theo_im_ms1_local;

      std::vector<OpenSwath::SwathMap> ms1_maps;
      for (const auto& m : swath_maps) {if (m.ms1) ms1_maps.push_back(m);}

      if (used_maps.empty())
      {
        continue;
      }

      // Get the spectrum for this RT and extract raw data points for all the
      // calibrating transitions (fragment m/z values) from the spectrum
      // Note that we are not using light clones of the underlying data here,
      // so access to the data needs to be in a critical section.
      OpenSwath::SpectrumPtr sp_ms1;
      OpenSwath::SpectrumPtr sp_ms2;

      #pragma omp critical (fetch_spectrum)
      {
        RangeMobility im_range;
        if (ms1_im_)
        {
          std::vector<OpenSwath::SpectrumPtr> fetchSpectrumArr = OpenSwathScoring().fetchSpectrumSwath(ms1_maps, bestRT, 1, im_range);
          sp_ms1 = (!fetchSpectrumArr.empty()) ? fetchSpectrumArr[0] : *new(OpenSwath::SpectrumPtr);
        }
        else
        {
          std::vector<OpenSwath::SpectrumPtr> fetchSpectrumArr = OpenSwathScoring().fetchSpectrumSwath(used_maps, bestRT, 1, im_range);
          sp_ms2 = (!fetchSpectrumArr.empty()) ? fetchSpectrumArr[0] : *new(OpenSwath::SpectrumPtr);
        }
      }

      for (const auto& tr : transition_group->getTransitions())
      {
        if (ms1_im_) {continue;}
        double intensity(0), im(0), mz(0);
        RangeMZ mz_range = DIAHelpers::createMZRangePPM(tr.product_mz, mz_extr_window, ppm);

        // get drift time upper/lower offset (this assumes that all chromatograms
        // are derived from the same precursor with the same drift time)
        auto pepref = tr.getPeptideRef();
        double drift_target = pep_im_map[pepref];
        RangeMobility im_range;
        if (im_extraction_win != -1 ) // im_extraction_win is set
        {
          im_range.setMax(drift_target);
          im_range.minSpanIfSingular(im_extraction_win);
        }

        // Check that the spectrum really has a drift time array
        if (sp_ms2->getDriftTimeArray() == nullptr)
        {
          OPENMS_LOG_DEBUG << "Did not find a drift time array for peptide " << pepref << " at RT " << bestRT  << std::endl;
          for (const auto& m : used_maps)
          {
            OPENMS_LOG_DEBUG << " -- Used maps " << m.lower << " to " << m.upper << " MS1 : " << m.ms1 << true << std::endl;
          }
          continue;
        }

        DIAHelpers::integrateWindow(sp_ms2, mz, im, intensity, mz_range, im_range);

        // skip empty windows
        if (im <= 0)
        {
          continue;
        }

        #pragma omp critical (accum_points)
        {
          // store result drift time
          data_im.push_back(std::make_pair(im, drift_target));
          exp_im.push_back(im);
          theo_im.push_back(drift_target);
          if (!debug_im_file_.empty())
          {
            os_im << tr.precursor_mz << "\t" << im << "\t" << drift_target << "\t" << bestRT << "\t" << intensity << std::endl;
          }
        }
        OPENMS_LOG_DEBUG << tr.precursor_mz << "\t" << im << "\t" << drift_target << "\t" << bestRT << "\t" << intensity << std::endl;
      }

      // Always collect a few MS1 IM points for window estimation (independent of ms1_im_)
      // However, if there are no MS1 Maps, then we have to skip window estimation for MS1
      if (!transition_group->getTransitions().empty() && !ms1_maps.empty())
      {
        const auto& tr0 = transition_group->getTransitions()[0];
        const auto pepref0 = tr0.getPeptideRef();
        const double drift_target0 = pep_im_map[pepref0];

        // Define an IM window centered on the theoretical drift time of this peptide
        RangeMobility im_range_ms1(drift_target0);
        if (im_extraction_win != -1) im_range_ms1.minSpanIfSingular(im_extraction_win);

        // Fetch the MS1 spectrum at bestRT (protect raw access in critical section)
        OpenSwath::SpectrumPtr sp_ms1_collect;
        #pragma omp critical (fetch_spectrum)
        {
          std::vector<OpenSwath::SpectrumPtr> arr_ms1 =
            OpenSwathScoring().fetchSpectrumSwath(ms1_maps, bestRT, 1, im_range_ms1);
          sp_ms1_collect = (!arr_ms1.empty()) ? arr_ms1[0] : OpenSwath::SpectrumPtr();
        }

        if (sp_ms1_collect && sp_ms1_collect->getDriftTimeArray() != nullptr)
        {
          double mz = 0.0, im = 0.0, intensity = 0.0;

          // m/z range isn’t critical for IM; use a narrow window around precursor m/z
          RangeMZ dummy = DIAHelpers::createMZRangePPM(tr0.precursor_mz, mz_extr_window, ppm);

          DIAHelpers::integrateWindow(sp_ms1_collect, mz, im, intensity, dummy, im_range_ms1);

          if (im > 0.0)
          {
            exp_im_ms1_local.push_back(im);
            theo_im_ms1_local.push_back(drift_target0);
          }
        }
      }

      // Do MS1 extraction
      if (!transition_group->getTransitions().empty() && ms1_im_)
      {
        const auto& tr = transition_group->getTransitions()[0];
        double intensity(0), im(0), mz(0);
        RangeMZ mz_range = DIAHelpers::createMZRangePPM(tr.precursor_mz, mz_extr_window, ppm);

        // get drift time upper/lower offset (this assumes that all chromatograms
        // are derived from the same precursor with the same drift time)
        auto pepref = tr.getPeptideRef();
        double drift_target = pep_im_map[pepref];

        // do not need to check for IM because we are correcting IM
        RangeMobility im_range(drift_target);
        im_range.minSpanIfSingular(im_extraction_win);

        // Check that the spectrum really has a drift time array
        if (sp_ms1->getDriftTimeArray() == nullptr)
        {
          OPENMS_LOG_DEBUG << "Did not find a drift time array for peptide " << pepref << " at RT " << bestRT  << std::endl;
          for (const auto& m : used_maps)
          {
            OPENMS_LOG_DEBUG << " -- Used maps " << m.lower << " to " << m.upper << " MS1 : " << m.ms1 << true << std::endl;
          }
          continue;
        }

        DIAHelpers::integrateWindow(sp_ms1, mz, im, intensity, mz_range, im_range);

        // skip empty windows
        if (im <= 0)
        {
          continue;
        }

        #pragma omp critical (accum_points)
        {
          // store result drift time
          data_im.push_back(std::make_pair(im, drift_target));
          exp_im.push_back(im);
          theo_im.push_back(drift_target);
          if (!debug_im_file_.empty())
          {
            os_im << tr.precursor_mz << "\t" << im << "\t" << drift_target << "\t" << bestRT << "\t" << intensity << std::endl;
          }
        }
        OPENMS_LOG_DEBUG << tr.precursor_mz << "\t" << im << "\t" << drift_target << "\t" << bestRT << "\t" << intensity << std::endl;
      }

      #pragma omp critical (accum_points)
      {
        exp_im_ms1_all.insert(exp_im_ms1_all.end(),
                              exp_im_ms1_local.begin(), exp_im_ms1_local.end());
        theo_im_ms1_all.insert(theo_im_ms1_all.end(),
                               theo_im_ms1_local.begin(), theo_im_ms1_local.end());
      }
    }

    if (!debug_im_file_.empty()) {os_im.close();}

    if (exp_im.empty())
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-LinearRegression", String("Could not fit a linear model to the data (0 points)."));
    }

    // linear correction is default (none returns in the beginning of the function)
    std::vector<double> im_regression_params;
    double confidence_interval_P(0.0);
    LinearRegression lr;
    lr.computeRegression(confidence_interval_P, exp_im.begin(), exp_im.end(), theo_im.begin()); // to convert exp_im -> theoretical im
    im_regression_params.push_back(lr.getIntercept());
    im_regression_params.push_back(lr.getSlope());
    im_regression_params.push_back(0.0);

    std::cout << "# im regression parameters: Y = " << im_regression_params[0] << " + " <<
      im_regression_params[1] << " X + " << im_regression_params[2] << " X^2" << std::endl;

    // store IM transformation, using the selected model
    im_trafo.setDataPoints(data_im);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    String model_type = "linear";
    im_trafo.fitModel(model_type, model_params);

    // Estimate MS2 ion mobility window
    // Use the 0.99 quantile so the window covers ~99% of residuals, ignoring rare extremes (those that are potential outliers).
    double fragment_im_window = im_trafo.estimateWindow(0.99, false, true, im_estimation_padding_factor);
    setFragmentImWindow(fragment_im_window);

    if (!exp_im_ms1_all.empty())
    {
      TransformationDescription::DataPoints ms1_points;
      ms1_points.reserve(exp_im_ms1_all.size());
      for (Size i = 0; i < exp_im_ms1_all.size(); ++i)
      {
        // (x, y) = (experimental IM, theoretical IM)
        ms1_points.emplace_back(exp_im_ms1_all[i], theo_im_ms1_all[i]);
      }

      // Copy the fitted model; don't mutate im_trafo's datapoints
      TransformationDescription im_trafo_inv = im_trafo;
      im_trafo_inv.setDataPoints(ms1_points);
      // Use the 0.99 quantile so the window covers ~99% of residuals, ignoring rare extremes (those that are potential outliers).
      const double precursor_im_window = im_trafo_inv.estimateWindow(0.99, /*invert=*/false, /*full_window=*/true, im_estimation_padding_factor);
      setPrecursorImWindow(precursor_im_window);
    }

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctIM done." << std::endl;
  }

  void SwathMapMassCorrection::correctMZ(
    const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    const bool pasef)
  {
    bool ppm = mz_extraction_window_ppm_;
    double mz_extr_window = mz_extraction_window_;
    std::string corr_type = mz_correction_function_;
    double im_extraction = im_extraction_window_;
    double mz_estimation_padding_factor = mz_estimation_padding_factor_;

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctMZ with type " << corr_type << " and window " << mz_extr_window << " in ppm " << ppm << std::endl;

    bool is_ppm = bool(corr_type == "quadratic_regression_delta_ppm" ||
                       corr_type == "weighted_quadratic_regression_delta_ppm" ||
                       corr_type == "regression_delta_ppm");

    if (corr_type == "none")
    {
      return;
    }

    std::ofstream os;
    if (!debug_mz_file_.empty())
    {
      std::cout.precision(16);
      os.open(debug_mz_file_);
      os << "mz" << "\t" << "theo_mz" << "\t" << "drift_time" << "\t" << "diff_ppm" << "\t" << "log_intensity" << "\t" << "RT" << std::endl;
      os.precision(writtenDigits(double()));
    }

    TransformationDescription::DataPoints data_all;
    std::vector<double> weights;
    std::vector<double> exp_mz;
    std::vector<double> theo_mz;
    std::vector<double> delta_ppm;

    std::map<std::string, double> pep_im_map;
    for (const auto& cmp : targeted_exp.getCompounds())
    {
      pep_im_map[cmp.id] = cmp.drift_time;
    }

    // Collect MS1 residuals (Δppm) for precursor m/z window estimation
    std::vector<double> delta_ppm_ms1;

    // Gather MS1 maps once
    std::vector<OpenSwath::SwathMap> ms1_maps;
    for (const auto& m : swath_maps) if (m.ms1) ms1_maps.push_back(m);

    for (auto & trgroup_it : transition_group_map)
    {
      // we need at least one feature to find the best one
      auto transition_group = trgroup_it.second;

      const auto& tr = transition_group->getTransitions()[0];
      auto pepref = tr.getPeptideRef();
      double drift_target = pep_im_map[pepref];

      if (transition_group->getFeatures().empty()) continue;

      // Find the feature with the highest score
      double bestRT;
      findBestFeature(*transition_group, bestRT);
      // Get the corresponding SWATH map(s), for SONAR there will be more than one map
      std::vector<OpenSwath::SwathMap> used_maps;
      if (!pasef)
      {
        used_maps = findSwathMaps(*transition_group, swath_maps);
      }
      // If pasef then have to check for overlap across IM
      else
      {
        used_maps = findSwathMapsPasef(*transition_group, swath_maps);
      }

      if (used_maps.empty())
      {
        continue;
      }

      // if ion mobility extraction window is set than extract with ion mobility
      RangeMobility im_range;
      if (im_extraction != -1) // ion mobility extraction is set
      {
        im_range.setMax(drift_target);
        im_range.minSpanIfSingular(im_extraction);
      }

      // MS2
      // Get the spectrum for this RT and extract raw data points for all the
      // calibrating transitions (fragment m/z values) from the spectrum
      std::vector<OpenSwath::SpectrumPtr> spArr = OpenSwathScoring().fetchSpectrumSwath(used_maps, bestRT, 1, im_range);
      OpenSwath::SpectrumPtr sp = (!spArr.empty()) ? spArr[0] : *new(OpenSwath::SpectrumPtr);
      for (const auto& tr : transition_group->getTransitions())
      {
        double mz, intensity, im;
        RangeMZ mz_range = DIAHelpers::createMZRangePPM(tr.product_mz, mz_extr_window, ppm);
        bool centroided = false;

        // integrate spectrum at the position of the theoretical mass
        DIAHelpers::integrateWindow(sp, mz, im, intensity, mz_range, im_range, centroided); // Correct using the irt_im

        // skip empty windows
        if (mz == -1)
        {
          continue;
        }

        // store result masses
        data_all.push_back(std::make_pair(mz, tr.product_mz));
        // regression weight is the log2 intensity
        weights.push_back( log(intensity) / log(2.0) );
        exp_mz.push_back( mz );
        // y = target = theoretical
        theo_mz.push_back( tr.product_mz );
        double diff_ppm = (mz - tr.product_mz) * 1000000 / mz;
        // y = target = delta-ppm
        delta_ppm.push_back(diff_ppm);
        if (!debug_mz_file_.empty())
        {
          os << mz << "\t" << tr.product_mz << "\t" << drift_target << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;
        }
        OPENMS_LOG_DEBUG << mz << "\t" << tr.product_mz << "\t" << diff_ppm << "\t" << log(intensity) / log(2.0) << "\t" << bestRT << std::endl;
      }

      // MS1 precursor processing for Δppm residuals
      if (!ms1_maps.empty())
      {
        std::vector<OpenSwath::SpectrumPtr> spArr_ms1 =
          OpenSwathScoring().fetchSpectrumSwath(ms1_maps, bestRT, 1, im_range);
        OpenSwath::SpectrumPtr sp_ms1 = (!spArr_ms1.empty()) ? spArr_ms1[0] : OpenSwath::SpectrumPtr();

        if (sp_ms1)
        {
          // Use theoretical precursor m/z for the calibrant peptide
          const double theo_prec_mz = tr.precursor_mz;
          RangeMZ mz_range_ms1 = DIAHelpers::createMZRangePPM(theo_prec_mz, mz_extr_window, ppm);

          double mz{}, im{}, intensity{};
          bool centroided = false;
          DIAHelpers::integrateWindow(sp_ms1, mz, im, intensity, mz_range_ms1, im_range, centroided);

          if (mz != -1) // got a signal
          {
            const double dppm = (mz - theo_prec_mz) / theo_prec_mz * 1e6;
            delta_ppm_ms1.push_back(std::abs(dppm));
          }
        }
      }
    }

    // Estimate fragment mz window
    {
      std::ostringstream ss;
      const Size N = std::min<Size>(delta_ppm.size(), 20);
      ss << "[SwathMapMassCorrection::correctMZ] MS2 residuals (first "
         << N << " of " << delta_ppm.size() << "): ";
      for (Size i = 0; i < N; ++i) { if (i) ss << ", "; ss << delta_ppm[i]; }
      if (delta_ppm.size() > N) ss << ", ...";
      OPENMS_LOG_DEBUG << ss.str() << '\n';
    }
    // Use the 0.99 quantile so the window covers ~99% of residuals, ignoring rare extremes (those that are potential outliers).
    double fragment_mz_window = estimateWindow(delta_ppm, 0.99, true, mz_estimation_padding_factor);
    setFragmentMzWindow(fragment_mz_window);

    // Estimate precursor window from MS1 residuals (full width, ppm)
    if (!delta_ppm_ms1.empty())
    {
      std::sort(delta_ppm_ms1.begin(), delta_ppm_ms1.end());

      {
        std::ostringstream ss;
        const Size N = std::min<Size>(delta_ppm_ms1.size(), 20);
        ss << "[SwathMapMassCorrection::correctMZ] MS1 residuals (first "
           << N << " of " << delta_ppm_ms1.size() << "): ";
        for (Size i = 0; i < N; ++i) { if (i) ss << ", "; ss << delta_ppm_ms1[i]; }
        if (delta_ppm_ms1.size() > N) ss << ", ...";
        OPENMS_LOG_DEBUG << ss.str() << '\n';
      }
      // Use the 0.99 quantile so the window covers ~99% of residuals, ignoring rare extremes (those that are potential outliers).
      double precursor_mz_window = estimateWindow(delta_ppm_ms1, 0.99, true, mz_estimation_padding_factor);
      setPrecursorMzWindow(precursor_mz_window);
    }

    std::vector<double> regression_params;
    if (corr_type == "none" || data_all.size() < 3)
    {
      return;
    }
    else if (corr_type == "unweighted_regression")
    {
      double confidence_interval_P(0.0);
      LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_regression")
    {
      double confidence_interval_P(0.0);
      LinearRegression lr;
      lr.computeRegressionWeighted(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "quadratic_regression")
    {
      // Quadratic fit
      QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "weighted_quadratic_regression")
    {
      // Quadratic fit (weighted)
      QuadraticRegression qr;
      qr.computeRegressionWeighted(exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "regression_delta_ppm")
    {
      // Regression fit using ppm differences
      double confidence_interval_P(0.0);
      LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      QuadraticRegression qr;
      qr.computeRegressionWeighted(exp_mz.begin(), exp_mz.end(), delta_ppm.begin(), weights.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Unknown correction type " + corr_type);
    }

    printf("# mz regression parameters: Y = %g + %g X + %g X^2\n",
           regression_params[0],
           regression_params[1],
           regression_params[2]);

    OPENMS_LOG_DEBUG << "# mz regression parameters: Y = " << regression_params[0] << " + " <<
      regression_params[1] << " X + " << regression_params[2] << " X^2" << std::endl;

    if (!debug_mz_file_.empty()) {os.close();}

#ifdef SWATHMAPMASSCORRECTION_DEBUG
    double s_ppm_before = 0;
    double s_ppm_after = 0;
    for (TransformationDescription::DataPoints::iterator d = data_all.begin(); d != data_all.end(); ++d)
    {
      double ppm_before = (d->first - d->second) * 1000000 / d->first;
      double predict = d->first*d->first*regression_params[2] + d->first*regression_params[1]+regression_params[0];
      double ppm_after = ( predict - d->second) * 1000000 / d->first;
      if (is_ppm)
      {
        double new_mz = d->first - predict*d->first/1000000;
        ppm_after = ( new_mz - d->second) * 1000000 / d->first;
      }
      s_ppm_before += std::fabs(ppm_before);
      s_ppm_after += std::fabs(ppm_after);
    }
    std::cout <<" sum residual sq ppm before " << s_ppm_before << " / after " << s_ppm_after << std::endl;
#endif

    // Replace the swath files with a transforming wrapper.
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      swath_maps[i].sptr = std::shared_ptr<OpenSwath::ISpectrumAccess>(
        new SpectrumAccessQuadMZTransforming(swath_maps[i].sptr,
          regression_params[0], regression_params[1], regression_params[2], is_ppm));
    }

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctMZ done." << std::endl;
  }

  double SwathMapMassCorrection::estimateWindow(std::vector<double> residuals, double quantile, bool full_width, double padding_factor)
  {
    if (residuals.empty()) return 0.0;

    // Ensure residuals are absolute errors
    for (auto& d : residuals) d = std::abs(d);

    // Adaptive half-width (Tukey k=1.5, blend from robust→raw as tail density grows 1%→10%)
    // k=1.5 uses the standard Tukey upper fence (Q3 + 1.5·IQR) to cap sparse extremes proposed in Exploratory Data Analysis by John W. Tukey (1977)
    // r_sparse=0.01 means if ≤1% of |residuals| exceed the fence, treat them as outliers (favor robust quantile);
    // r_dense=0.10 means if ≥10% exceed the fence, tails are genuinely broad (favor raw quantile).
    // These values are conservative, widely used in stats.
    OpenMS::Math::AdaptiveQuantileResult adaptive_quantile_res = OpenMS::Math::adaptiveQuantile(
      residuals.begin(), residuals.end(),
      quantile);

    const double full = (full_width ? (2.0 * adaptive_quantile_res.blended) : adaptive_quantile_res.blended) * padding_factor;

    OPENMS_LOG_DEBUG
      << "[estimateWindow] n=" << residuals.size()
      << " q=" << quantile
      << " half_raw=" << adaptive_quantile_res.half_raw
      << " half_rob=" << adaptive_quantile_res.half_rob
      << " UF=" << adaptive_quantile_res.upper_fence
      << " tail_frac=" << adaptive_quantile_res.tail_fraction
      << " => half_adapt=" << adaptive_quantile_res.blended
      << " full=" << full
      << std::endl;

    return full;
  }

  double SwathMapMassCorrection::getFragmentMzWindow() const
  {
    return fragment_mz_window_;
  }

  void SwathMapMassCorrection::setFragmentMzWindow(double fragmentMzWindow)
  {
    fragment_mz_window_ = fragmentMzWindow;
  }

  double SwathMapMassCorrection::getFragmentImWindow() const
  {
    return fragment_im_window_;
  }

  void SwathMapMassCorrection::setFragmentImWindow(double fragmentImWindow)
  {
    fragment_im_window_ = fragmentImWindow;
  }

  double SwathMapMassCorrection::getPrecursorMzWindow() const
  {
    return precursor_mz_window_;
  }

  void SwathMapMassCorrection::setPrecursorMzWindow(double precursorMzWindow)
  {
    precursor_mz_window_ = precursorMzWindow;
  }

  double SwathMapMassCorrection::getPrecursorImWindow() const
  {
    return precursor_im_window_;
  }

  void SwathMapMassCorrection::setPrecursorImWindow(double precursorImWindow)
  {
    precursor_im_window_ = precursorImWindow;
  }

  }

