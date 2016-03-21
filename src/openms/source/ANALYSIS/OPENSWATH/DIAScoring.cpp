// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>

#include <numeric>
#include <algorithm>
#include <functional>

#include <boost/bind.hpp>

const double C13C12_MASSDIFF_U = 1.0033548;

namespace OpenMS
{
  DIAScoring::DIAScoring() :
    DefaultParamHandler("DIAScoring")
  {

    defaults_.setValue("dia_extraction_window", 0.05, "DIA extraction window in Th.");
    defaults_.setMinFloat("dia_extraction_window", 0.0);
    defaults_.setValue("dia_centroided", "false", "Use centroded DIA data.");
    defaults_.setValidStrings("dia_centroided", ListUtils::create<String>("true,false"));
    defaults_.setValue("dia_byseries_intensity_min", 300.0, "DIA b/y series minimum intensity to consider.");
    defaults_.setMinFloat("dia_byseries_intensity_min", 0.0);
    defaults_.setValue("dia_byseries_ppm_diff", 10.0, "DIA b/y series minimal difference in ppm to consider.");
    defaults_.setMinFloat("dia_byseries_ppm_diff", 0.0);

    defaults_.setValue("dia_nr_isotopes", 4, "DIA nr of isotopes to consider.");
    defaults_.setMinInt("dia_nr_isotopes", 0);
    defaults_.setValue("dia_nr_charges", 4, "DIA nr of charges to consider.");
    defaults_.setMinInt("dia_nr_charges", 0);

    defaults_.setValue("peak_before_mono_max_ppm_diff", 20.0, "DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.");
    defaults_.setMinFloat("peak_before_mono_max_ppm_diff", 0.0);

    // write defaults into Param object param_
    defaultsToParam_();
  }

  void DIAScoring::updateMembers_()
  {
    dia_extract_window_ = (double)param_.getValue("dia_extraction_window");
    dia_centroided_ = param_.getValue("dia_centroided").toBool();
    dia_byseries_intensity_min_ = (double)param_.getValue("dia_byseries_intensity_min");
    dia_byseries_ppm_diff_ = (double)param_.getValue("dia_byseries_ppm_diff");

    dia_nr_isotopes_ = (int)param_.getValue("dia_nr_isotopes");
    dia_nr_charges_ = (int)param_.getValue("dia_nr_charges");
    peak_before_mono_max_ppm_diff_ = (double)param_.getValue("peak_before_mono_max_ppm_diff");
  }

  void DIAScoring::set_dia_parameters(double dia_extract_window, double dia_centroided,
                                      double dia_byseries_intensity_min, double dia_byseries_ppm_diff, double dia_nr_isotopes, double dia_nr_charges)
  {
    dia_extract_window_ = dia_extract_window;
    dia_centroided_ = dia_centroided;
    dia_byseries_intensity_min_ = dia_byseries_intensity_min;
    dia_byseries_ppm_diff_ = dia_byseries_ppm_diff;

    dia_nr_isotopes_ = dia_nr_isotopes;
    dia_nr_charges_ = dia_nr_charges;
  }

  ///////////////////////////////////////////////////////////////////////////
  // DIA / SWATH scoring

  void DIAScoring::dia_isotope_scores(const std::vector<TransitionType>& transitions, SpectrumPtrType spectrum,
                                      OpenSwath::IMRMFeature* mrmfeature, double& isotope_corr, double& isotope_overlap)
  {
    isotope_corr = 0;
    isotope_overlap = 0;
    // first compute a map of relative intensities from the feature, then compute the score
    std::map<std::string, double> intensities;
    getFirstIsotopeRelativeIntensities_(transitions, mrmfeature, intensities);
    diaIsotopeScoresSub_(transitions, spectrum, intensities, isotope_corr, isotope_overlap);
  }

  void DIAScoring::dia_massdiff_score(const std::vector<TransitionType>& transitions, SpectrumPtrType spectrum,
                                      const std::vector<double>& normalized_library_intensity,
                                      double& ppm_score, double& ppm_score_weighted)
  {
    ppm_score = 0;
    ppm_score_weighted = 0;
    double mz, intensity;
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType* transition = &transitions[k];
      // Calculate the difference of the theoretical mass and the actually measured mass
      double left = transition->getProductMZ() - dia_extract_window_ / 2.0;
      double right = transition->getProductMZ() + dia_extract_window_ / 2.0;
      bool signalFound = integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!signalFound)
      {
        continue;
      }

      double diff_ppm = std::fabs(mz - transition->getProductMZ()) * 1000000 / transition->getProductMZ();
      ppm_score += diff_ppm;
      ppm_score_weighted += diff_ppm * normalized_library_intensity[k];
#ifdef MRMSCORING_TESTING
      std::cout << " weighted int of the peak is " << mz << " diff is in ppm " << diff_ppm << " thus append " << diff_ppm * diff_ppm << " or weighted " << diff_ppm * normalized_library_intensity[k] << std::endl;
#endif
    }
  }

  bool DIAScoring::dia_ms1_massdiff_score(double precursor_mz, SpectrumPtrType spectrum,
                                          double& ppm_score)
  {
    ppm_score = -1;
    double mz, intensity;
    {
      // Calculate the difference of the theoretical mass and the actually measured mass
      double left = precursor_mz - dia_extract_window_ / 2.0;
      double right = precursor_mz + dia_extract_window_ / 2.0;
      bool signalFound = integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Catch if no signal was found and replace it with the most extreme
      // value. Otherwise calculate the difference in ppm.
      if (!signalFound)
      {
        ppm_score = dia_extract_window_ / precursor_mz * 1000000;
        return false;
      }
      else
      {
        ppm_score = std::fabs(mz - precursor_mz) * 1000000 / precursor_mz;
        return true;
      }
    }
  }

  /// Precursor isotope scores
  void DIAScoring::dia_ms1_isotope_scores(double precursor_mz, SpectrumPtrType spectrum, size_t charge_state, 
                                          double& isotope_corr, double& isotope_overlap)
  {
    // collect the potential isotopes of this peak
    double max_ratio;
    int nr_occurences;
    std::vector<double> isotopes_int;
    for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
    {
      double left  = precursor_mz - dia_extract_window_ / 2.0 + iso * C13C12_MASSDIFF_U / static_cast<double>(charge_state);
      double right = precursor_mz + dia_extract_window_ / 2.0 + iso * C13C12_MASSDIFF_U / static_cast<double>(charge_state);
      double mz, intensity;
      integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      isotopes_int.push_back(intensity);
    }

    // calculate the scores:
    // isotope correlation (forward) and the isotope overlap (backward) scores
    isotope_corr = scoreIsotopePattern_(precursor_mz, isotopes_int, charge_state);
    largePeaksBeforeFirstIsotope_(spectrum, precursor_mz, isotopes_int[0], nr_occurences, max_ratio);
    isotope_overlap = max_ratio;
  }

  void DIAScoring::dia_by_ion_score(SpectrumPtrType spectrum,
                                    AASequence& sequence, int charge, double& bseries_score,
                                    double& yseries_score)
  {
    bseries_score = 0;
    yseries_score = 0;
    OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer");

    double mz, intensity, left, right;
    std::vector<double> yseries, bseries;
    OpenMS::DIAHelpers::getBYSeries(sequence, bseries, yseries, charge);
    for (Size it = 0; it < bseries.size(); it++)
    {
      left = bseries[it] - dia_extract_window_ / 2.0;
      right = bseries[it] + dia_extract_window_ / 2.0;
      bool signalFound = integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      double ppmdiff = std::fabs(bseries[it] - mz) * 1000000 / bseries[it];
      if (signalFound && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        bseries_score++;
      }
    }
    for (Size it = 0; it < yseries.size(); it++)
    {
      left = yseries[it] - dia_extract_window_ / 2.0;
      right = yseries[it] + dia_extract_window_ / 2.0;
      bool signalFound = integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      double ppmdiff = std::fabs(yseries[it] - mz) * 1000000 / yseries[it];
      if (signalFound && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        yseries_score++;
      }
    }
  }

  void DIAScoring::score_with_isotopes(SpectrumPtrType spectrum, const std::vector<TransitionType>& transitions,
                                       double& dotprod, double& manhattan)
  {
    OpenMS::DiaPrescore dp(dia_extract_window_, dia_nr_isotopes_, dia_nr_charges_);
    dp.score(spectrum, transitions, dotprod, manhattan);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private methods

  /// computes a vector of relative intensities for each feature (output to intensities)
  void DIAScoring::getFirstIsotopeRelativeIntensities_(
    const std::vector<TransitionType>& transitions,
    OpenSwath::IMRMFeature* mrmfeature, std::map<std::string, double>& intensities)
  {
    for (Size k = 0; k < transitions.size(); k++)
    {
      std::string native_id = transitions[k].getNativeID();
      double rel_intensity = mrmfeature->getFeature(native_id)->getIntensity() / mrmfeature->getIntensity();
      intensities.insert(std::pair<std::string, double>(native_id, rel_intensity));
    }
  }

  void DIAScoring::diaIsotopeScoresSub_(const std::vector<TransitionType>& transitions, SpectrumPtrType spectrum,
                                          std::map<std::string, double>& intensities, //relative intensities
                                          double& isotope_corr, double& isotope_overlap)
  {
    std::vector<double> isotopes_int;
    double max_ratio;
    int nr_occurences;
    for (Size k = 0; k < transitions.size(); k++)
    {
      isotopes_int.clear();
      String native_id = transitions[k].getNativeID();
      double rel_intensity = intensities[native_id];

      // If no charge is given, we assume it to be 1
      int putative_fragment_charge = 1;
      if (transitions[k].product_charge > 0)
      {
        putative_fragment_charge = transitions[k].product_charge;
      }

      // collect the potential isotopes of this peak
      for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
      {
        double left = transitions[k].getProductMZ() - dia_extract_window_ / 2.0 + iso * C13C12_MASSDIFF_U / static_cast<double>(putative_fragment_charge);
        double right = transitions[k].getProductMZ() + dia_extract_window_ / 2.0 + iso * C13C12_MASSDIFF_U / static_cast<double>(putative_fragment_charge);
        double mz, intensity;
        integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
        isotopes_int.push_back(intensity);
      }

      // calculate the scores:
      // isotope correlation (forward) and the isotope overlap (backward) scores
      double score = scoreIsotopePattern_(transitions[k].getProductMZ(), isotopes_int, putative_fragment_charge);
      isotope_corr += score * rel_intensity;
      largePeaksBeforeFirstIsotope_(spectrum, transitions[k].getProductMZ(), isotopes_int[0], nr_occurences, max_ratio);
      isotope_overlap += nr_occurences * rel_intensity;
    }
  }

  void DIAScoring::largePeaksBeforeFirstIsotope_(SpectrumPtrType spectrum, double mono_mz, double mono_int, int& nr_occurences, double& max_ratio)
  {
    double mz, intensity;
    nr_occurences = 0;
    max_ratio = 0.0;

    for (int ch = 1; ch <= dia_nr_charges_; ++ch)
    {
      double left = mono_mz - dia_extract_window_ / 2.0 - C13C12_MASSDIFF_U / (double) ch;
      double right = mono_mz + dia_extract_window_ / 2.0 - C13C12_MASSDIFF_U / (double) ch;
      bool signalFound = integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!signalFound)
      {
        continue;
      }

      // Compute ratio between the (presumed) monoisotopic peak intensity and the now found peak
      double ratio;
      if (mono_int != 0) { ratio = intensity / mono_int; }
      else { ratio = 0; }
      if (ratio > max_ratio) {max_ratio = ratio;}
        
      double ddiff_ppm = std::fabs(mz - (mono_mz - 1.0 / (double) ch)) * 1000000 / mono_mz;

      // FEATURE we should fit a theoretical distribution to see whether we really are a secondary peak
      if (ratio > 1 && ddiff_ppm < peak_before_mono_max_ppm_diff_)
      {
        //isotope_overlap += 1.0 * rel_intensity;

        nr_occurences += 1.0; // we count how often this happens...

#ifdef MRMSCORING_TESTING
        cout << " _ overlap diff ppm  " << ddiff_ppm << " and inten ratio " << ratio << " with " << mono_int << endl;
#endif
      }
    }
  }

  double DIAScoring::scoreIsotopePattern_(double product_mz,
                                             const std::vector<double>& isotopes_int, int putative_fragment_charge)
  {
    OPENMS_PRECONDITION(putative_fragment_charge > 0, "Charge is a positive integer");

    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;

    // create the theoretical distribution
    IsotopeDistribution d;
    TheoreticalIsotopePattern isotopes;
    d.setMaxIsotope(dia_nr_isotopes_ + 1);
    d.estimateFromPeptideWeight(product_mz * putative_fragment_charge);
    for (IsotopeDistribution::Iterator it = d.begin(); it != d.end(); ++it)
    {
      isotopes.intensity.push_back(it->second);
    }

    //FEATURE ISO pattern for peptide sequence..

    isotopes.optional_begin = 0;
    isotopes.optional_end = dia_nr_isotopes_;

    //scale the distribution to a maximum of 1
    double max = 0.0;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      if (isotopes.intensity[i] > max)
      {
        max = isotopes.intensity[i];
      }
    }
    isotopes.max = max;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      isotopes.intensity[i] /= max;
    }
    isotopes.trimmed_left = 0;

    // score the pattern against a theoretical one
    double int_score = OpenSwath::cor_pearson(isotopes_int.begin(), isotopes_int.end(), isotopes.intensity.begin());
    if (boost::math::isnan(int_score))
    {
      int_score = 0;
    }
    return int_score;

  } //end of dia_isotope_corr_sub

}
