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
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>

// Classes
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>

// Functions
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

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
    // Get the corresponding SWATH map(s), for SONAR there will be more than one map
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

  SwathMapMassCorrection::SwathMapMassCorrection() :
    DefaultParamHandler("SwathMapMassCorrection")
  {
    defaults_.setValue("mz_extraction_window", -1.0, "M/z extraction window width");
    defaults_.setValue("mz_extraction_window_ppm", "false", "Whether m/z extraction is in ppm", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("mz_extraction_window_ppm", ListUtils::create<String>("true,false"));
    defaults_.setValue("ms1_im_calibration", "false", "Whether to use MS1 precursor data for the ion mobility calibration (default = false, uses MS2 / fragment ions for calibration)", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("ms1_im_calibration", ListUtils::create<String>("true,false"));
    defaults_.setValue("im_extraction_window", -1.0, "Ion mobility extraction window width");
    defaults_.setValue("mz_correction_function", "none", "Type of normalization function for m/z calibration.");
    defaults_.setValidStrings("mz_correction_function", ListUtils::create<String>("none,regression_delta_ppm,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));
    defaults_.setValue("im_correction_function", "linear", "Type of normalization function for IM calibration.");
    defaults_.setValidStrings("im_correction_function", ListUtils::create<String>("none,linear"));

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
    mz_correction_function_ = param_.getValue("mz_correction_function");
    im_correction_function_ = param_.getValue("im_correction_function");
    debug_mz_file_ = param_.getValue("debug_mz_file");
    debug_im_file_ = param_.getValue("debug_im_file");
  }

  void SwathMapMassCorrection::correctIM(
    const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    TransformationDescription& im_trafo,
    const OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    bool ppm = mz_extraction_window_ppm_;
    double mz_extr_window = mz_extraction_window_;
    double im_extraction_win = im_extraction_window_;

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
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (SignedSize k = 0; k < (SignedSize)trgr_ids.size(); k++)
    {
      // we need at least one feature to find the best one
      auto transition_group = transition_group_map.at(trgr_ids[k]);
      if (transition_group->getFeatures().empty()) continue;

      // Find the feature with the highest score
      double bestRT;
      findBestFeature(*transition_group, bestRT);
      // Get the corresponding SWATH map(s), for SONAR there will be more than one map
      std::vector<OpenSwath::SwathMap> used_maps = findSwathMaps(*transition_group, swath_maps);

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
      OpenSwath::SpectrumPtr sp;
      OpenSwath::SpectrumPtr sp_ms1;
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        if (ms1_im_)
        {
          sp_ms1 = OpenSwathScoring().fetchSpectrumSwath(ms1_maps, bestRT, 1, 0, 0);
        }
        else
        {
          sp = OpenSwathScoring().fetchSpectrumSwath(used_maps, bestRT, 1, 0, 0);
        }
      }

      for (const auto& tr : transition_group->getTransitions())
      {
        if (ms1_im_) {continue;}
        double intensity(0), im(0), left(tr.product_mz), right(tr.product_mz);

        // get drift time upper/lower offset (this assumes that all chromatograms
        // are derived from the same precursor with the same drift time)
        auto pepref = tr.getPeptideRef();
        double drift_target = pep_im_map[pepref];
        double drift_left(drift_target), drift_right(drift_target);
        DIAHelpers::adjustExtractionWindow(drift_right, drift_left, im_extraction_win, false);

        // Check that the spectrum really has a drift time array
        if (sp->getDriftTimeArray() == nullptr)
        {
          OPENMS_LOG_DEBUG << "Did not find a drift time array for peptide " << pepref << " at RT " << bestRT  << std::endl;
          for (const auto& m : used_maps) OPENMS_LOG_DEBUG << " -- Used maps " << m.lower << " to " << m.upper << " MS1 : " << m.ms1 << true << std::endl;
          continue;
        }

        DIAHelpers::adjustExtractionWindow(right, left, mz_extr_window, ppm);
        DIAHelpers::integrateDriftSpectrum(sp, left, right, im, intensity, drift_left, drift_right);

        // skip empty windows
        if (im <= 0)
        {
          continue;
        }

#ifdef _OPENMP
#pragma omp critical
#endif
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

      // Do MS1 extraction
      if (!transition_group->getTransitions().empty() && ms1_im_)
      {
        const auto& tr = transition_group->getTransitions()[0];
        double intensity(0), im(0), left(tr.precursor_mz), right(tr.precursor_mz);

        // get drift time upper/lower offset (this assumes that all chromatograms
        // are derived from the same precursor with the same drift time)
        auto pepref = tr.getPeptideRef();
        double drift_target = pep_im_map[pepref];
        double drift_left(drift_target), drift_right(drift_target);
        DIAHelpers::adjustExtractionWindow(drift_right, drift_left, im_extraction_win, false);

        // Check that the spectrum really has a drift time array
        if (sp_ms1->getDriftTimeArray() == nullptr)
        {
          OPENMS_LOG_DEBUG << "Did not find a drift time array for peptide " << pepref << " at RT " << bestRT  << std::endl;
          for (const auto& m : used_maps) OPENMS_LOG_DEBUG << " -- Used maps " << m.lower << " to " << m.upper << " MS1 : " << m.ms1 << true << std::endl;
          continue;
        }

        DIAHelpers::adjustExtractionWindow(right, left, mz_extr_window, ppm);
        DIAHelpers::integrateDriftSpectrum(sp, left, right, im, intensity, drift_left, drift_right);

        // skip empty windows
        if (im <= 0)
        {
          continue;
        }

#ifdef _OPENMP
#pragma omp critical
#endif
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
    }

    if (!debug_im_file_.empty()) {os_im.close();}

    // linear correction is default (none returns in the beginning of the function)
    std::vector<double> im_regression_params;
    double confidence_interval_P(0.0);
    Math::LinearRegression lr;
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

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctIM done." << std::endl;
  }

  void SwathMapMassCorrection::correctMZ(
    const std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> & transition_group_map,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    bool ppm = mz_extraction_window_ppm_;
    double mz_extr_window = mz_extraction_window_;
    std::string corr_type = mz_correction_function_;

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
      std::vector<OpenSwath::SwathMap> used_maps = findSwathMaps(*transition_group, swath_maps);

      if (used_maps.empty())
      {
        continue;
      }

      // Get the spectrum for this RT and extract raw data points for all the
      // calibrating transitions (fragment m/z values) from the spectrum
      OpenSwath::SpectrumPtr sp = OpenSwathScoring().fetchSpectrumSwath(used_maps, bestRT, 1, 0, 0);
      for (const auto& tr : transition_group->getTransitions())
      {
        double mz, intensity, left(tr.product_mz), right(tr.product_mz);
        bool centroided = false;

        // integrate spectrum at the position of the theoretical mass
        DIAHelpers::adjustExtractionWindow(right, left, mz_extr_window, ppm);
        DIAHelpers::integrateWindow(sp, left, right, mz, intensity, centroided);

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
    }

    std::vector<double> regression_params;
    if (corr_type == "none" || data_all.size() < 3)
    {
      return;
    }
    else if (corr_type == "unweighted_regression")
    {
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_regression")
    {
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegressionWeighted(confidence_interval_P, exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "quadratic_regression")
    {
      // Quadratic fit
      Math::QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), theo_mz.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "weighted_quadratic_regression")
    {
      // Quadratic fit (weighted)
      Math::QuadraticRegression qr;
      qr.computeRegressionWeighted(exp_mz.begin(), exp_mz.end(), theo_mz.begin(), weights.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      Math::QuadraticRegression qr;
      qr.computeRegression(exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(qr.getA());
      regression_params.push_back(qr.getB());
      regression_params.push_back(qr.getC());
    }
    else if (corr_type == "regression_delta_ppm")
    {
      // Regression fit using ppm differences
      double confidence_interval_P(0.0);
      Math::LinearRegression lr;
      lr.computeRegression(confidence_interval_P, exp_mz.begin(), exp_mz.end(), delta_ppm.begin());
      regression_params.push_back(lr.getIntercept());
      regression_params.push_back(lr.getSlope());
      regression_params.push_back(0.0);
    }
    else if (corr_type == "weighted_quadratic_regression_delta_ppm")
    {
      // Quadratic fit using ppm differences
      Math::QuadraticRegression qr;
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
      swath_maps[i].sptr = boost::shared_ptr<OpenSwath::ISpectrumAccess>(
        new SpectrumAccessQuadMZTransforming(swath_maps[i].sptr,
          regression_params[0], regression_params[1], regression_params[2], is_ppm));
    }

    OPENMS_LOG_DEBUG << "SwathMapMassCorrection::correctMZ done." << std::endl;
  }

}

