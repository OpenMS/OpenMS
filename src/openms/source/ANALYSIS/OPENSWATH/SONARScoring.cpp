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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/SONARScoring.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp> // for isnan

// #define DEBUG_SONAR

namespace OpenMS
{
  SONARScoring::SONARScoring() :
    DefaultParamHandler("SONARScoring")
  {
    defaults_.setValue("dia_extraction_window", 0.05, "DIA extraction window in Th or ppm.");
    defaults_.setMinFloat("dia_extraction_window", 0.0);
    defaults_.setValue("dia_extraction_unit", "Th", "DIA extraction window unit");
    defaults_.setValidStrings("dia_extraction_unit", ListUtils::create<String>("Th,ppm"));
    defaults_.setValue("dia_centroided", "false", "Use centroided DIA data.");
    defaults_.setValidStrings("dia_centroided", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
  }

  void SONARScoring::updateMembers_()
  {
    dia_extract_window_ = (double)param_.getValue("dia_extraction_window");
    dia_extraction_ppm_ = param_.getValue("dia_extraction_unit") == "ppm";
    dia_centroided_ = param_.getValue("dia_centroided").toBool();
  }

  void SONARScoring::computeXCorr_(std::vector<std::vector<double> >& sonar_profiles,
                     double& xcorr_coelution_score, double& xcorr_shape_score)
  {
    /// Cross Correlation array
    typedef OpenSwath::Scoring::XCorrArrayType XCorrArrayType;
    /// Cross Correlation matrix
    typedef std::vector<std::vector<XCorrArrayType> > XCorrMatrixType;

    XCorrMatrixType xcorr_matrix;
    xcorr_matrix.resize(sonar_profiles.size());
    for (std::size_t i = 0; i < sonar_profiles.size(); i++)
    {
      xcorr_matrix[i].resize(sonar_profiles.size());
      for (std::size_t j = i; j < sonar_profiles.size(); j++)
      {
        // compute normalized cross correlation
        xcorr_matrix[i][j] = OpenSwath::Scoring::normalizedCrossCorrelation(
                sonar_profiles[i], sonar_profiles[j], boost::numeric_cast<int>(sonar_profiles[i].size()), 1);
      }
    }

    // coelution (lag score)
    std::vector<int> deltas;
    for (std::size_t i = 0; i < xcorr_matrix.size(); i++)
    {
      for (std::size_t  j = i; j < xcorr_matrix.size(); j++)
      {
        // first is the lag value, should be an int
        deltas.push_back(std::abs(OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr_matrix[i][j])->first));
      }
    }

    {
      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(deltas.begin(), deltas.end(), msc);
      double deltas_mean = msc.mean();
      double deltas_stdv = msc.sample_stddev();
      xcorr_coelution_score = deltas_mean + deltas_stdv;
    }



    // shape score (intensity)
    std::vector<double> intensities;
    for (std::size_t i = 0; i < xcorr_matrix.size(); i++)
    {
      for (std::size_t j = i; j < xcorr_matrix.size(); j++)
      {
        // second is the Y value (intensity)
        intensities.push_back(OpenSwath::Scoring::xcorrArrayGetMaxPeak(xcorr_matrix[i][j])->second);
      }
    }
    {
      OpenSwath::mean_and_stddev msc;
      msc = std::for_each(intensities.begin(), intensities.end(), msc);
      xcorr_shape_score = msc.mean();
    }
  }

  void SONARScoring::computeSonarScores(OpenSwath::IMRMFeature* imrmfeature,
                                        const std::vector<OpenSwath::LightTransition> & transitions,
                                        std::vector<OpenSwath::SwathMap>& swath_maps,
                                        OpenSwath_Scores & scores)
  {
    if (transitions.empty()) {return;}

    double precursor_mz = transitions[0].getPrecursorMZ();

#ifdef DEBUG_SONAR
    std::ofstream debug_file;
    debug_file.open("debug_sonar_profiles.tsv",  std::fstream::in | std::fstream::out | std::fstream::app);

    String native_id = 0;
    if (transitions.size() > 0)
    {
      native_id = transitions[0].getNativeID();
    }
    debug_file << native_id << "\t" << imrmfeature->getRT() << "\tcentr";
    for (Size it = 0; it < swath_maps.size(); it++)
    {
      debug_file << "\t" << (swath_maps[it].lower + swath_maps[it].upper) / 2.0;
    }
    debug_file << "\n";


    std::cout << " doing RT " << imrmfeature->getRT() << " using maps: ";
    for (Size i  = 0; i < swath_maps.size() ; i++)
    {
      std::cout << (swath_maps[i].lower + swath_maps[i].upper) / 2 << " ";
    }
    std::cout << std::endl;

    // idea 1: check the elution profile of each SONAR scan ...
    for (Size kk = 0; kk < imrmfeature->getNativeIDs().size(); kk++)
    {
      std::vector<double> rt;
      imrmfeature->getFeature(imrmfeature->getNativeIDs()[kk])->getRT(rt);
    }
#endif


    // idea 2: check the SONAR profile (e.g. in the dimension of) of the best scan (RT apex)
    double RT = imrmfeature->getRT();

    // Aggregate sonar profiles (for each transition)
    std::vector<std::vector<double> > sonar_profiles;
    std::vector<double> sn_score;
    std::vector<double> diff_score;
    std::vector<double> trend_score;
    std::vector<double> rsq_score;
    std::vector<double> mz_median_score;
    std::vector<double> mz_stdev_score;
    for (Size k = 0; k < transitions.size(); k++)
    {
      String native_id = transitions[k].getNativeID();

      // Gather profiles across all SONAR maps
      std::vector<double> sonar_profile;
      std::vector<double> sonar_mz_profile;
      std::vector<bool> signal_exp;
      for (Size swath_idx = 0; swath_idx < swath_maps.size(); swath_idx++)
      {
        OpenSwath::SpectrumAccessPtr swath_map = swath_maps[swath_idx].sptr;

        bool expect_signal = false;
        if (swath_maps[swath_idx].ms1) {continue;} // skip MS1
        if (precursor_mz > swath_maps[swath_idx].lower && precursor_mz < swath_maps[swath_idx].upper)
        {
          expect_signal = true;
        }

        // find closest scan for current SONAR map (by retention time)
        std::vector<std::size_t> indices = swath_map->getSpectraByRT(RT, 0.0);
        if (indices.empty() )  {continue;}
        int closest_idx = boost::numeric_cast<int>(indices[0]);
        if (indices[0] != 0 &&
            std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
            std::fabs(swath_map->getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
        {
          closest_idx--;
        }
        OpenSwath::SpectrumPtr spectrum_ = swath_map->getSpectrumById(closest_idx);

        // integrate intensity within that scan
        double left = transitions[k].getProductMZ();
        double right = transitions[k].getProductMZ();
        if (dia_extraction_ppm_)
        {
          left -= left * dia_extract_window_ / 2e6;
          right += right * dia_extract_window_ / 2e6;
        }
        else
        {
          left -= dia_extract_window_ / 2.0;
          right += dia_extract_window_ / 2.0;
        }
        double mz, intensity;
        integrateWindow(spectrum_, left, right, mz, intensity, dia_centroided_);

        sonar_profile.push_back(intensity);
        sonar_mz_profile.push_back(mz);
        signal_exp.push_back(expect_signal);
      }
      sonar_profiles.push_back(sonar_profile);

#ifdef DEBUG_SONAR
      std::cout << " transition " << native_id << " at " << RT << " will analyze with " << swath_maps.size() << " maps" << std::endl;
      debug_file << native_id << "\t" << imrmfeature->getRT() << "\tint";
      for (Size it = 0; it < sonar_profile.size(); it++)
      {
        debug_file << "\t" << sonar_profile[it];
      }
      debug_file << "\n";
      debug_file << native_id << "\t" << imrmfeature->getRT() << "\tmz";
      for (Size it = 0; it < sonar_mz_profile.size(); it++)
      {
        debug_file << "\t" << sonar_mz_profile[it];
      }
      debug_file << "\n";
#endif

      // Analyze profiles
      std::vector<double> sonar_profile_pos;
      std::vector<double> sonar_mz_profile_pos;
      std::vector<double> sonar_profile_neg;
      std::vector<double> sonar_mz_profile_neg;
      for (Size it = 0; it < sonar_profile.size(); it++)
      {
        if (signal_exp[it])
        {
          sonar_profile_pos.push_back(sonar_profile[it]);
          sonar_mz_profile_pos.push_back(sonar_mz_profile[it]);
        }
        else
        {
          sonar_profile_neg.push_back(sonar_profile[it]);
          sonar_mz_profile_neg.push_back(sonar_mz_profile[it]);
        }
      }

      // try to find diff between first and last
      double sonar_trend = 1.0;
      if (sonar_profile_pos.size() > 1)
      {
        double int_end = sonar_profile_pos[sonar_profile_pos.size()-1] + sonar_profile_pos[sonar_profile_pos.size()-2];
        double int_start = sonar_profile_pos[0] + sonar_profile_pos[1];
        if (int_end > 0.0)
        {
          sonar_trend = int_start / int_end;
        }
        else
        {
          sonar_trend = 0.0;
        }
      }

      // try to find R^2 of a linear regression (optimally, there is no trend)
      std::vector<double> xvals;
      for (Size pr_idx = 0; pr_idx < sonar_profile_pos.size(); pr_idx++) {xvals.push_back(pr_idx);}
      double rsq = OpenSwath::cor_pearson( xvals.begin(), xvals.end(), sonar_profile_pos.begin() );
      if (boost::math::isnan(rsq)) rsq = 0.0; // check for nan

      // try to find largest diff
      double sonar_largediff = 0.0;
      for (Size pr_idx = 0; pr_idx < sonar_profile_pos.size()-1; pr_idx++)
      {
        double diff = std::fabs(sonar_profile_pos[pr_idx] - sonar_profile_pos[pr_idx+1]);
        if (diff > sonar_largediff)
        {
          sonar_largediff = diff;
        }
      }

      double sonar_sn = 1.0;
      double pos_med = 1.0;
      // from here on, its not sorted any more !!
      if (!sonar_profile_pos.empty() && !sonar_profile_neg.empty())
      {
        double neg_med;

        pos_med = Math::median(sonar_profile_pos.begin(), sonar_profile_pos.end());
        neg_med = Math::median(sonar_profile_neg.begin(), sonar_profile_neg.end());

        // compute the relative difference between the medians (or if the
        // medians are zero, compute the difference to the max element)
        if (neg_med > 0.0)
        {
          sonar_sn = pos_med / neg_med;
        }
        else if (*std::max_element(sonar_profile_neg.begin(), sonar_profile_neg.end()) > 0.0)
        {
          sonar_sn = pos_med / *std::max_element(sonar_profile_neg.begin(), sonar_profile_neg.end());
        }

      }

      double median_mz = 0.0;
      double mz_stdev = -1.0;
      if (!sonar_mz_profile_pos.empty())
      {
        median_mz = Math::median(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end());

        double sum = std::accumulate(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end(), 0.0);
        double mean = sum / sonar_mz_profile_pos.size();

        double sq_sum = std::inner_product(sonar_mz_profile_pos.begin(), sonar_mz_profile_pos.end(), sonar_mz_profile_pos.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / sonar_mz_profile_pos.size() - mean * mean);

        mz_stdev = stdev;
      }

#ifdef DEBUG_SONAR
      std::cout << " computed SN: " << sonar_sn  <<  "(from " << pos_med << " and neg " << neg_med <<  ")"
        << " large diff: "  << sonar_largediff << " trend " << sonar_trend << std::endl;
#endif
      sn_score.push_back(sonar_sn);
      diff_score.push_back(sonar_largediff / pos_med);
      trend_score.push_back(sonar_trend);
      rsq_score.push_back(rsq);

      mz_median_score.push_back(median_mz);
      mz_stdev_score.push_back(mz_stdev);
    }

    double xcorr_coelution_score, xcorr_shape_score;
    computeXCorr_(sonar_profiles, xcorr_coelution_score, xcorr_shape_score);

    double sn_av = std::accumulate(sn_score.begin(), sn_score.end(), 0.0) / sn_score.size();
    double diff_av = std::accumulate(diff_score.begin(), diff_score.end(), 0.0) / diff_score.size();
    double trend_av = std::accumulate(trend_score.begin(), trend_score.end(), 0.0) / trend_score.size();
    double rsq_av = std::accumulate(rsq_score.begin(), rsq_score.end(), 0.0) / rsq_score.size();

    //double mz_median = std::accumulate(mz_median_score.begin(), mz_median_score.end(), 0.0) / mz_median_score.size();
    //double mz_stdev = std::accumulate(mz_stdev_score.begin(), mz_stdev_score.end(), 0.0) / mz_stdev_score.size();

    scores.sonar_sn = sn_av;
    scores.sonar_diff = diff_av;
    scores.sonar_trend = trend_av;
    scores.sonar_rsq = rsq_av;
    scores.sonar_lag = xcorr_coelution_score;
    scores.sonar_shape = xcorr_shape_score;

#ifdef DEBUG_SONAR
    debug_file.close();
#endif
  }


}

