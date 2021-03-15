// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/MATH/MISC/MathFunctions.h> // getPPM

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

    defaults_.setValue("dia_extraction_window", 0.05, "DIA extraction window in Th or ppm.");
    defaults_.setMinFloat("dia_extraction_window", 0.0);
    defaults_.setValue("dia_extraction_unit", "Th", "DIA extraction window unit");
    defaults_.setValidStrings("dia_extraction_unit", ListUtils::create<String>("Th,ppm"));
    defaults_.setValue("dia_centroided", "false", "Use centroided DIA data.");
    defaults_.setValidStrings("dia_centroided", ListUtils::create<String>("true,false"));
    defaults_.setValue("dia_byseries_intensity_min", 300.0, "DIA b/y series minimum intensity to consider.");
    defaults_.setMinFloat("dia_byseries_intensity_min", 0.0);
    defaults_.setValue("dia_byseries_ppm_diff", 10.0, "DIA b/y series minimal difference in ppm to consider.");
    defaults_.setMinFloat("dia_byseries_ppm_diff", 0.0);

    defaults_.setValue("dia_nr_isotopes", 4, "DIA number of isotopes to consider.");
    defaults_.setMinInt("dia_nr_isotopes", 0);
    defaults_.setValue("dia_nr_charges", 4, "DIA number of charges to consider.");
    defaults_.setMinInt("dia_nr_charges", 0);

    defaults_.setValue("peak_before_mono_max_ppm_diff", 20.0, "DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.");
    defaults_.setMinFloat("peak_before_mono_max_ppm_diff", 0.0);

    // write defaults into Param object param_
    defaultsToParam_();

    // for void getBYSeries
    {
      generator = new TheoreticalSpectrumGenerator();
      Param p;
      p.setValue("add_metainfo", "true",
          "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
      generator->setParameters(p);
  }

    // for simulateSpectrumFromAASequence
    //  Param p;
    //  p.setValue("add_metainfo", "false",
    //      "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    //  p.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    //  generator->setParameters(p);
  }

  DIAScoring::~DIAScoring() 
  {
    delete generator;
  }

  void DIAScoring::updateMembers_()
  {
    dia_extract_window_ = (double)param_.getValue("dia_extraction_window");
    dia_extraction_ppm_ = param_.getValue("dia_extraction_unit") == "ppm";
    dia_centroided_ = param_.getValue("dia_centroided").toBool();
    dia_byseries_intensity_min_ = (double)param_.getValue("dia_byseries_intensity_min");
    dia_byseries_ppm_diff_ = (double)param_.getValue("dia_byseries_ppm_diff");

    dia_nr_isotopes_ = (int)param_.getValue("dia_nr_isotopes");
    dia_nr_charges_ = (int)param_.getValue("dia_nr_charges");
    peak_before_mono_max_ppm_diff_ = (double)param_.getValue("peak_before_mono_max_ppm_diff");
  }

  ///////////////////////////////////////////////////////////////////////////
  // DIA / SWATH scoring

  void DIAScoring::dia_isotope_scores(const std::vector<TransitionType>& transitions, SpectrumPtrType spectrum,
                                      OpenSwath::IMRMFeature* mrmfeature, double& isotope_corr, double& isotope_overlap) const
  {
    isotope_corr = 0;
    isotope_overlap = 0;
    // first compute a map of relative intensities from the feature, then compute the score
    std::map<std::string, double> intensities;
    getFirstIsotopeRelativeIntensities_(transitions, mrmfeature, intensities);
    diaIsotopeScoresSub_(transitions, spectrum, intensities, isotope_corr, isotope_overlap);
  }

  void DIAScoring::dia_massdiff_score(const std::vector<TransitionType>& transitions,
                                      SpectrumPtrType spectrum,
                                      const std::vector<double>& normalized_library_intensity,
                                      double& ppm_score,
	                                    double& ppm_score_weighted,
                                      std::vector<double>& diff_ppm) const
  {
    ppm_score = 0;
    ppm_score_weighted = 0;
    diff_ppm.clear();
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType& transition = transitions[k];
      // Calculate the difference of the theoretical mass and the actually measured mass
      double left(transition.getProductMZ()), right(transition.getProductMZ());
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      double mz, intensity;
      bool signalFound = DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!signalFound)
      {
        continue;
      }

      double ppm = Math::getPPM(mz, transition.getProductMZ());
      diff_ppm.push_back(transition.getProductMZ());
      diff_ppm.push_back(ppm);
      ppm_score += std::fabs(ppm);
      ppm_score_weighted += std::fabs(ppm) * normalized_library_intensity[k];
#ifdef MRMSCORING_TESTING
      std::cout << " weighted int of the peak is " << mz << " diff is in ppm " << diff_ppm << " thus append " << diff_ppm * diff_ppm << " or weighted " << diff_ppm * normalized_library_intensity[k] << std::endl;
#endif
    }

    // FEATURE we should not punish so much when one transition is missing!
    ppm_score /= transitions.size();
  }

  bool DIAScoring::dia_ms1_massdiff_score(double precursor_mz, SpectrumPtrType spectrum,
                                          double& ppm_score) const
  {
    ppm_score = -1;
    double mz, intensity;
    {
      // Calculate the difference of the theoretical mass and the actually measured mass
      double left(precursor_mz), right(precursor_mz);
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      bool signalFound = DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Catch if no signal was found and replace it with the most extreme
      // value. Otherwise calculate the difference in ppm.
      if (!signalFound)
      {
        ppm_score = (right - left) / precursor_mz * 1000000;
        return false;
      }
      else
      {
        ppm_score = Math::getPPMAbs(mz, precursor_mz);
        return true;
      }
    }
  }

  /// Precursor isotope scores
  void DIAScoring::dia_ms1_isotope_scores(double precursor_mz, SpectrumPtrType spectrum,
                                          double& isotope_corr, double& isotope_overlap, const EmpiricalFormula& sum_formula) const
  {
    // although precursor_mz can be received from the empirical formula (if non-empty), the actual precursor could be
    // slightly different. And also for compounds, usually the neutral sum_formula without adducts is given.
    // Therefore calculate the isotopes based on the formula but place them at precursor_mz
    std::vector<double> isotopes_int;
    getIsotopeIntysFromExpSpec_(precursor_mz, spectrum, isotopes_int, sum_formula.getCharge());

    double max_ratio = 0;
    int nr_occurrences = 0;

    // calculate the scores:
    // isotope correlation (forward) and the isotope overlap (backward) scores
    isotope_corr = scoreIsotopePattern_(isotopes_int, sum_formula);
    largePeaksBeforeFirstIsotope_(spectrum, precursor_mz, isotopes_int[0], nr_occurrences, max_ratio);
    isotope_overlap = max_ratio;
  }

  void DIAScoring::getIsotopeIntysFromExpSpec_(double precursor_mz, SpectrumPtrType spectrum,
                            std::vector<double>& isotopes_int,
                            int charge_state) const
  {
    double abs_charge = std::fabs(static_cast<double>(charge_state));
    for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
    {
      double left  = precursor_mz + iso * C13C12_MASSDIFF_U / abs_charge;
      double right = left;
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      double mz, intensity;
      DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      isotopes_int.push_back(intensity);
    }
  }

  void DIAScoring::dia_ms1_isotope_scores_averagine(double precursor_mz, SpectrumPtrType spectrum,
                                                    double& isotope_corr, double& isotope_overlap,
                                                    int charge_state) const
  {
    std::vector<double> exp_isotopes_int;
    getIsotopeIntysFromExpSpec_(precursor_mz, spectrum, exp_isotopes_int, charge_state);
    CoarseIsotopePatternGenerator solver(dia_nr_isotopes_ + 1);
    // NOTE: this is a rough estimate of the neutral mz value since we would not know the charge carrier for negative ions
    IsotopeDistribution isotope_dist = solver.estimateFromPeptideWeight(std::fabs(precursor_mz * charge_state));

    double max_ratio;
    int nr_occurrences;
    // calculate the scores:
    // isotope correlation (forward) and the isotope overlap (backward) scores
    isotope_corr = scoreIsotopePattern_(exp_isotopes_int, isotope_dist);
    largePeaksBeforeFirstIsotope_(spectrum, precursor_mz, exp_isotopes_int[0], nr_occurrences, max_ratio);
    isotope_overlap = max_ratio;
  }

  void DIAScoring::dia_by_ion_score(SpectrumPtrType spectrum,
                                    AASequence& sequence, int charge, double& bseries_score,
                                    double& yseries_score) const
  {
    bseries_score = 0;
    yseries_score = 0;
    OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer"); // for peptides, charge should be positive

    double mz, intensity, left, right;
    std::vector<double> yseries, bseries;
    OpenMS::DIAHelpers::getBYSeries(sequence, bseries, yseries, generator, charge);
    for (const auto& b_ion_mz : bseries)
    {
      left = b_ion_mz;
      right = b_ion_mz;
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);

      bool signalFound = DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      double ppmdiff = Math::getPPMAbs(mz, b_ion_mz);
      if (signalFound && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        bseries_score++;
      }
    }
    for (const auto& y_ion_mz : yseries)
    {
      left = y_ion_mz;
      right = y_ion_mz;
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);

      bool signalFound = DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      double ppmdiff = Math::getPPMAbs(mz, y_ion_mz);
      if (signalFound && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        yseries_score++;
      }
    }
  }

  void DIAScoring::score_with_isotopes(SpectrumPtrType spectrum, const std::vector<TransitionType>& transitions,
                                       double& dotprod, double& manhattan) const
  {
    OpenMS::DiaPrescore dp(dia_extract_window_, dia_nr_isotopes_, dia_nr_charges_);
    dp.score(spectrum, transitions, dotprod, manhattan);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private methods

  /// computes a vector of relative intensities for each feature (output to intensities)
  void DIAScoring::getFirstIsotopeRelativeIntensities_(
    const std::vector<TransitionType>& transitions,
    OpenSwath::IMRMFeature* mrmfeature, std::map<std::string, double>& intensities) const
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
                                        double& isotope_corr,
                                        double& isotope_overlap) const
  {
    std::vector<double> isotopes_int;
    double max_ratio;
    int nr_occurences;
    for (Size k = 0; k < transitions.size(); k++)
    {
      isotopes_int.clear();
      const String native_id = transitions[k].getNativeID();
      double rel_intensity = intensities[native_id];

      // If no charge is given, we assume it to be 1
      int putative_fragment_charge = 1;
      if (transitions[k].fragment_charge != 0)
      {
        putative_fragment_charge = transitions[k].fragment_charge;
      }

      // collect the potential isotopes of this peak
      double abs_charge = std::fabs(static_cast<double>(putative_fragment_charge));
      for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
      {
        double left = transitions[k].getProductMZ() + iso * C13C12_MASSDIFF_U / abs_charge;
        double right = left;
        DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
        double mz, intensity;
        DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);
        isotopes_int.push_back(intensity);
      }

      // calculate the scores:
      // isotope correlation (forward) and the isotope overlap (backward) scores
      double score = scoreIsotopePattern_(isotopes_int, transitions[k].getProductMZ(), putative_fragment_charge);
      isotope_corr += score * rel_intensity;
      largePeaksBeforeFirstIsotope_(spectrum, transitions[k].getProductMZ(), isotopes_int[0], nr_occurences, max_ratio);
      isotope_overlap += nr_occurences * rel_intensity;
    }
  }

  void DIAScoring::largePeaksBeforeFirstIsotope_(SpectrumPtrType spectrum, double mono_mz, double mono_int, int& nr_occurences, double& max_ratio) const
  {
    double mz, intensity;
    nr_occurences = 0;
    max_ratio = 0.0;

    for (int ch = 1; ch <= dia_nr_charges_; ++ch)
    {
      double center = mono_mz - C13C12_MASSDIFF_U / (double) ch;
      double left = center;
      double right = center;
      DIAHelpers::adjustExtractionWindow(right, left, dia_extract_window_, dia_extraction_ppm_);
      bool signalFound = DIAHelpers::integrateWindow(spectrum, left, right, mz, intensity, dia_centroided_);

      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!signalFound)
      {
        continue;
      }

      // Compute ratio between the (presumed) monoisotopic peak intensity and the now found peak
      double ratio;
      if (mono_int != 0)
      {
        ratio = intensity / mono_int;
      }
      else
      {
        ratio = 0;
      }
      if (ratio > max_ratio) {max_ratio = ratio;}

      double ddiff_ppm = std::fabs(mz - center) * 1e6 / center;

      // FEATURE we should fit a theoretical distribution to see whether we really are a secondary peak
      if (ratio > 1 && ddiff_ppm < peak_before_mono_max_ppm_diff_)
      {
        //isotope_overlap += 1.0 * rel_intensity;

        nr_occurences += 1; // we count how often this happens...

#ifdef MRMSCORING_TESTING
        cout << " _ overlap diff ppm  " << ddiff_ppm << " and inten ratio " << ratio << " with " << mono_int << endl;
#endif
      }
    }
  }

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          double product_mz,
                                          int putative_fragment_charge) const
  {
    OPENMS_PRECONDITION(putative_fragment_charge != 0, "Charge needs to be set to != 0"); // charge can be positive and negative

    IsotopeDistribution isotope_dist;

    // create the theoretical distribution from the peptide weight
    CoarseIsotopePatternGenerator solver(dia_nr_isotopes_ + 1);
    // NOTE: this is a rough estimate of the neutral mz value since we would not know the charge carrier for negative ions
    isotope_dist = solver.estimateFromPeptideWeight(std::fabs(product_mz * putative_fragment_charge));

    return scoreIsotopePattern_(isotopes_int, isotope_dist);
  } //end of dia_isotope_corr_sub

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          const EmpiricalFormula& empf) const
  {
    return scoreIsotopePattern_(isotopes_int,
                                empf.getIsotopeDistribution(CoarseIsotopePatternGenerator(dia_nr_isotopes_ + 1)));
  }

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          const IsotopeDistribution& isotope_dist) const
  {
    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;

    TheoreticalIsotopePattern isotopes;

    for (IsotopeDistribution::ConstIterator it = isotope_dist.begin(); it != isotope_dist.end(); ++it)
    {
      isotopes.intensity.push_back(it->getIntensity());
    }
    isotopes.optional_begin = 0;
    isotopes.optional_end = dia_nr_isotopes_;

    // scale the distribution to a maximum of 1
    double max = 0.0;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      if (isotopes.intensity[i] > max)
      {
        max = isotopes.intensity[i];
      }
    }
    isotopes.max = max;
    if (max == 0.) max = 1.;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      isotopes.intensity[i] /= max;
    }
    isotopes.trimmed_left = 0;

    // score the pattern against a theoretical one
    OPENMS_POSTCONDITION(isotopes_int.size() == isotopes.intensity.size(), "Vectors for pearson correlation do not have the same size.");
    double int_score = OpenSwath::cor_pearson(isotopes_int.begin(), isotopes_int.end(), isotopes.intensity.begin());
    if (boost::math::isnan(int_score))
    {
      int_score = 0;
    }
    return int_score;
  } //end of dia_isotope_corr_sub
}
