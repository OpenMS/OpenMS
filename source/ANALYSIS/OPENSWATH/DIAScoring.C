// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/corr.h"
#include <OpenMS/ANALYSIS/OPENSWATH/OpenMSHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DiaPrescoring.h>

#include <numeric>
#include <algorithm>
#include <functional>

#include <boost/bind.hpp>

// const double C13C12_MASSDIFF_U = 1.0033548; // ou
namespace OpenSwath
{
  ///////////////////////////////////////////////////////////////////////////
  // DIA / SWATH scoring

  void DIAScoring::dia_isotope_scores(const std::vector<TransitionType> & transitions, SpectrumType spectrum,
    OpenSwath::IMRMFeature * mrmfeature, double & isotope_corr, double & isotope_overlap)
  {
    // first compute a map of relative intensities from the feature, then compute the score
    std::map<std::string, double> intensities;
    getFirstIsotopeRelativeIntensities(transitions, mrmfeature, intensities);
    dia_isotope_scores_sub(transitions, spectrum, intensities, isotope_corr, isotope_overlap);
  } 

  void DIAScoring::dia_massdiff_score(const std::vector<TransitionType> & transitions, SpectrumType spectrum,
    const std::vector<double> & normalized_library_intensity,
    double & ppm_score, double & ppm_score_weighted)
  {
    ppm_score = 0;
    ppm_score_weighted = 0;
    double mz, intensity, left, right;
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType * transition = &transitions[k];
      // Calculate the difference of the theoretical mass and the actually measured mass
      left = transition->getProductMZ() - dia_extract_window_ / 2.0;
      right = transition->getProductMZ() + dia_extract_window_ / 2.0;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      if (mz == -1)
      {
        mz = (left + right) / 2.0;
      }

      //double diff = std::fabs( mz - transition->getProductMZ() );
      double diff_ppm = std::fabs(mz - transition->getProductMZ()) * 1000000 / transition->getProductMZ();
      //ppm_score += diff_ppm * diff_ppm;
      ppm_score += diff_ppm;
      ppm_score_weighted += diff_ppm * normalized_library_intensity[k];
#ifdef MRMSCORING_TESTING
      std::cout << " weighted int of the peak is " << mz << " diff is in ppm " << diff_ppm << " thus append " << diff_ppm * diff_ppm << " or weighted " << diff_ppm * normalized_library_intensity[k] << std::endl;
#endif
    }
  }

  void DIAScoring::dia_by_ion_score(SpectrumType spectrum,
    AASequence & sequence, int charge, double & bseries_score,
    double & yseries_score)
  {
    OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer");

    double mz, intensity, left, right;
    std::vector<double> yseries, bseries;
    std::vector<double> extr_yseries, extr_bseries;
    OpenMS::getBYSeries(sequence, bseries, yseries, charge);
    double ppmdiff;
    for (Size it = 0; it < bseries.size(); it++)
    {
      left = bseries[it] - dia_extract_window_ / 2.0;
      right = bseries[it] + dia_extract_window_ / 2.0;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      ppmdiff = std::fabs(bseries[it] - mz) * 1000000 / bseries[it];
      if (mz != -1 && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        bseries_score++;
      }
    }
    for (Size it = 0; it < yseries.size(); it++)
    {
      left = yseries[it] - dia_extract_window_ / 2.0;
      right = yseries[it] + dia_extract_window_ / 2.0;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      ppmdiff = std::fabs(yseries[it] - mz) * 1000000 / yseries[it];
      if (mz != -1 && ppmdiff < dia_byseries_ppm_diff_ && intensity > dia_byseries_intensity_min_)
      {
        yseries_score++;
      }
    }
  }

  void DIAScoring::score_with_isotopes(SpectrumType spectrum, const std::vector<TransitionType> & transitions,
    double & dotprod, double & manhattan)
  {
    OpenMS::DiaPrescore2 dp;
    dp.score(spectrum, transitions, dotprod, manhattan);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private methods

  /// computes a vector of relative intensities for each feature (output to intensities) 
  void DIAScoring::getFirstIsotopeRelativeIntensities(const std::vector<TransitionType> & transitions,
    OpenSwath::IMRMFeature * mrmfeature, std::map<std::string, double> & intensities)
  {
    double rel_intensity;
    for (Size k = 0; k < transitions.size(); k++)
    {
      std::string native_id = transitions[k].getNativeID(); 
      rel_intensity = mrmfeature->getFeature(native_id)->getIntensity() / mrmfeature->getIntensity();
      intensities.insert(std::pair<std::string, double>(native_id, rel_intensity));
    }
  }

#if 0

  void DIAScoring::getSpectrumIntensities(const std::vector<TransitionType> & transitions, SpectrumType spectrum, double extractWindow,
    std::vector<double> & mzv, std::vector<double> & intensityv)
  {
    for (std::size_t k = 0; k < transitions.size(); ++k)
    {
      // Calculate the difference of the theoretical mass and the actually measured mass
      double left = transitions[k].getProductMZ() - extractWindow / 2.0;
      double right = transitions[k].getProductMZ() + extractWindow / 2.0;
      double mz, intensity;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      if (mz == -1)
      {
        mz = (left + right) / 2.0;
      }
      mzv.push_back(mz);
      intensityv.push_back(intensity);
    }
  }
  //gets the relative intensities from spectrum
  void DIAScoring::getFirstIsotopeRelativeIntensities(
    const std::vector<TransitionType> & transitions,
    SpectrumType spectrum,
    std::map<std::string, double> & intensities   //experimental intensities of transitions
    )
  {
    typedef std::pair<std::string, double> TPair;
    typedef std::vector<TransitionType> Transitiontypevec;
    Transitiontypevec::const_iterator beg = transitions.begin();
    Transitiontypevec::const_iterator end = transitions.end();

    std::vector<double> tmpint;
    std::vector<std::string> refs;
    double mz, intensity;
    for (; beg != end; ++beg)
    {
      std::string native_id = beg->getNativeID();
      double left = beg->getProductMZ() - dia_extract_window_ / 2.0;
      double right = beg->getProductMZ() + dia_extract_window_ / 2.0;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      tmpint.push_back(intensity);
      refs.push_back(native_id);
    }

    //compute total intensities
    double totalInt = std::accumulate(tmpint.begin(), tmpint.end(), 0.0);
    //normalize intensities
    if (totalInt > 0)
    {
      std::transform(tmpint.begin(), tmpint.end(), tmpint.begin(),
                     std::bind2nd(std::divides<double>(), totalInt));
    }
    //add to result
    for (uint32_t i = 0; i < tmpint.size(); ++i)
    {
      intensities.insert(TPair(refs[i], tmpint[i]));
    }
  }
#endif

  void DIAScoring::dia_isotope_scores_sub(const std::vector<TransitionType> & transitions, SpectrumType spectrum,
    std::map<std::string, double> & intensities,  //relative intensities
    double & isotope_corr, double & isotope_overlap)
  {
    std::vector<double> isotopes_int;
    double max_ppm_diff = 20.0; // TODO (hroest) make this a proper parameter

    for (Size k = 0; k < transitions.size(); k++)
    {
      isotopes_int.clear();
      String native_id = transitions[k].getNativeID(); 
      double rel_intensity = intensities[native_id];

      // If no charge is given, we assume it to be 1
      int putative_fragment_charge = 1; 
      if (transitions[k].charge > 0)
      {
        putative_fragment_charge = transitions[k].charge; 
      }

      // collect the potential isotopes of this peak
      for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
      {
        // TODO (hroest) multiply with C13C12_MASSDIFF_U -->  +iso * C13C12_MASSDIFF_U / charge
        // TODO (hroest) get fragment charge from transition --> add to Transition interface
        double left = transitions[k].getProductMZ() - dia_extract_window_ / 2.0 + iso / static_cast<DoubleReal>(putative_fragment_charge);
        double right = transitions[k].getProductMZ() + dia_extract_window_ / 2.0 + iso / static_cast<DoubleReal>(putative_fragment_charge);
        double mz, intensity;
        getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
        isotopes_int.push_back(intensity);
      }

      // calculate the scores:
      // isotope correlation (forward) and the isotope overlap (backward) scores
      double score = scoreIsotopePattern(transitions[k].getProductMZ(), isotopes_int, putative_fragment_charge);
      isotope_corr += score * rel_intensity;
      score = largePeaksBeforeFirstIsotope(transitions[k].getProductMZ(), spectrum, max_ppm_diff, isotopes_int[0]);
      isotope_overlap += score * rel_intensity;
    }
  }

  /// Search for a large peak _before_ (lower m/z) the current peak 
  DoubleReal DIAScoring::largePeaksBeforeFirstIsotope(double product_mz,
    SpectrumType & spectrum, double max_ppm_diff, double main_peak)
  {
    double result = 0;
    double mz, intensity, left, right, ratio;

    for (int ch = 1; ch <= dia_nr_charges_; ++ch)
    {
      // TODO (hroest) multiply with C13C12_MASSDIFF_U -->  + 1.0 * C13C12_MASSDIFF_U / charge
      left = product_mz - dia_extract_window_ / 2.0 - 1.0 / (DoubleReal) ch;
      right = product_mz + dia_extract_window_ / 2.0 - 1.0 / (DoubleReal) ch;
      getIntensePeakInWindow(spectrum, left, right, mz, intensity, dia_centroided_);
      if (mz == -1)
      {
        mz = (left + right) / 2.0;
      }

      ratio = intensity / main_peak;
      if (main_peak == 0)
      {
        ratio = 0;
      }
      double ddiff_ppm = std::fabs(mz - (product_mz - 1.0 / (DoubleReal) ch)) * 1000000 / product_mz;

      // FEATURE we should fit a theoretical distribution to see whether we really are a secondary peak
      if (ratio > 1 && ddiff_ppm < max_ppm_diff)
      {
        //isotope_overlap += 1.0 * rel_intensity;

        result += 1.0; // we count how often this happens...

#ifdef MRMSCORING_TESTING
        cout << " _ overlap diff ppm  " << ddiff_ppm << " and inten ratio " << ratio << " with " << main_peak << endl;
#endif
      }
    }
    return result;
  }

  /// Compare an experimental isotope pattern to a theoretical one 
  DoubleReal DIAScoring::scoreIsotopePattern(double product_mz,
    const std::vector<double> & isotopes_int, int putative_fragment_charge)
  {
    OPENMS_PRECONDITION(putative_fragment_charge > 0, "Charge is a positive integer");

    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;
    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern IsotopePattern;

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
    DoubleReal max = 0.0;
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
    DoubleReal int_score = OpenSwath::cor_pearson(isotopes_int.begin(), isotopes_int.end(), isotopes.intensity.begin());
    if (boost::math::isnan(int_score))
    {
      int_score = 0;
    }
    return int_score;

  } //end of dia_isotope_corr_sub

  void DIAScoring::getIntensePeakInWindow(const SpectrumType spectrum, double mz_start, 
      double mz_end, double & mz, double & intensity, bool centroided)
  {
    OPENMS_PRECONDITION(
        std::adjacent_find(spectrum->getMZArray()->data.begin(), spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
      "MZ vector needs to be sorted!");

    intensity = 0;
    if (!centroided)
    {
      // get the weighted average for noncentroided data.
      // TODO this is not optimal if there are two peaks in this window (e.g. if
      // the window is too large)
      mz = 0;
      intensity = 0;

      std::vector<double>::const_iterator mz_arr_end =
        spectrum->getMZArray()->data.end();
      std::vector<double>::const_iterator int_it =
        spectrum->getIntensityArray()->data.begin();

      // this assumes that the spectra are sorted!
      std::vector<double>::const_iterator mz_it = std::lower_bound(
        spectrum->getMZArray()->data.begin(),
        spectrum->getMZArray()->data.end(), mz_start);
      std::vector<double>::const_iterator mz_it_end = std::lower_bound(mz_it,
                                                                       mz_arr_end, mz_end);

      // also advance intensity iterator now
      int iterator_pos =
        std::distance(
          (std::vector<double>::const_iterator)spectrum->getMZArray()->data.begin(),
          mz_it);
      std::advance(int_it, iterator_pos);

      for (; mz_it != mz_it_end; ++mz_it, ++int_it)
      {
        intensity += (*int_it);
        mz += (*int_it) * (*mz_it);
      }

      // This is equivalent to (but faster)
      // * std::vector<double>::const_iterator mz_arr_end = spectrum->getMZArray()->data.end();
      // * std::vector<double>::const_iterator int_it = spectrum->getIntensityArray()->data.begin();
      // * std::vector<double>::const_iterator mz_it = spectrum->getMZArray()->data.begin();
      // * for (; mz_it != mz_arr_end; mz_it++, int_it++)
      // * {
      // *   if ((*mz_it) < mz_start) continue;
      // *   if ((*mz_it) > mz_end) break;
      // *   std::cout << " add " <<  (*int_it) << " @ " << (*mz_it) << std::endl;
      // *   intensity += (*int_it);
      // *   mz += (*int_it) * (*mz_it);
      // * }

      if (intensity != 0)
      {
        mz /= intensity;
      }
      else
      {
        mz = -1;
        intensity = 0;
      }

    }
    else
    {
      // not implemented
      throw "Not implemented";
    }
  }

//mrmfeature
}
