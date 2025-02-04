// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <boost/math/distributions/binomial.hpp>
#include <numeric>

using namespace std;

namespace OpenMS
{

  float XQuestScores::preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta)
  {

    if ( (matched_alpha <= 0 && matched_beta <= 0) || ions_alpha <= 0 || ions_beta <= 0)
    {
      return 0.0;
    }

    // avoid 0 values in multiplication, adds a "dynamic range" among candidates with no matching linear peaks to one of the peptides
    float matched_alpha_float = matched_alpha;
    if (matched_alpha <= 0)
    {
      matched_alpha_float = 0.1f;
    }
    float matched_beta_float = matched_beta;
    if (matched_beta <= 0)
    {
      matched_beta_float = 0.1f;
    }

      float result = sqrt((static_cast<float>(matched_alpha_float) / static_cast<float>(ions_alpha)) * (static_cast<float>(matched_beta_float) / static_cast<float>(ions_beta)));
      return result;
  }

  float XQuestScores::preScore(Size matched_alpha, Size ions_alpha)
  {
    if (ions_alpha <= 0)
    {
      return 0.0;
    }

    float result = static_cast<float>(matched_alpha) / static_cast<float>(ions_alpha);
    return result;
  }

  double XQuestScores::matchOddsScore(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, Size n_charges)
  {
    using boost::math::binomial;
    Size theo_size = theoretical_spec.size();

    if (matched_size < 1 || theo_size < 1)
    {
      return 0;
    }

    double range = theoretical_spec[theo_size-1].getMZ() -  theoretical_spec[0].getMZ();

    // Compute fragment tolerance in Da for the mean of MZ values, if tolerance in ppm (rough approximation)
    double mean = 0.0;
    for (Size i = 0; i < theo_size; ++i)
    {
      mean += theoretical_spec[i].getMZ();
    }
    mean = mean / theo_size;
    double tolerance_Th = fragment_mass_tolerance_unit_ppm ? mean * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;

    // A priori probability of a random match given info about the theoretical spectrum
    double a_priori_p = 0;

    if (is_xlink_spectrum)
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  (static_cast<double>(theo_size) / static_cast<double>(n_charges)))));
    }
    else
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  static_cast<int>(theo_size))));
    }

    double match_odds = 0;

    binomial flip(theo_size, a_priori_p);
    // min double number to avoid 0 values, causing scores with the value "inf"
    match_odds = -log(cdf(complement(flip, matched_size)) + std::numeric_limits<double>::min());

    // score lower than 0 does not make sense, but can happen if cfd = 0, -log( 1 + min() ) < 0
    if (match_odds >= 0.0)
    {
      return match_odds;
    }
    else
    {
      return 0;
    }
  }

  double XQuestScores::matchOddsScoreSimpleSpec(const std::vector< SimpleTSGXLMS::SimplePeak >& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, Size n_charges)
  {
    using boost::math::binomial;
    Size theo_size = theoretical_spec.size();

    if (matched_size < 1 || theo_size < 1)
    {
      return 0;
    }

    double range = theoretical_spec[theo_size-1].mz - theoretical_spec[0].mz;

    // Compute fragment tolerance in Da for the mean of MZ values, if tolerance in ppm (rough approximation)
    double mean = 0.0;
    for (Size i = 0; i < theo_size; ++i)
    {
      mean += theoretical_spec[i].mz;
    }
    mean = mean / theo_size;
    double tolerance_Th = fragment_mass_tolerance_unit_ppm ? mean * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;

    // A priori probability of a random match given info about the theoretical spectrum
    double a_priori_p = 0;

    if (is_xlink_spectrum)
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  (static_cast<double>(theo_size) / static_cast<double>(n_charges)))));
    }
    else
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  static_cast<int>(theo_size))));
    }

    double match_odds = 0;

    binomial flip(theo_size, a_priori_p);
    // min double number to avoid 0 values, causing scores with the value "inf"
    match_odds = -log(cdf(complement(flip, matched_size)) + std::numeric_limits<double>::min());

    // score lower than 0 does not make sense, but can happen if cfd = 0, -log( 1 + min() ) < 0
    if (match_odds >= 0.0)
    {
      return match_odds;
    }
    else
    {
      return 0;
    }
  }

  double XQuestScores::logOccupancyProb(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
  {
    using boost::math::binomial;
    Size theo_size = theoretical_spec.size();

    if (matched_size < 1 || theo_size < 1)
    {
      return 0;
    }

    double range;
    double used_tolerance;

    if (fragment_mass_tolerance_unit_ppm)
    {
      range = std::log(theoretical_spec.back().getMZ()) - std::log(theoretical_spec[0].getMZ());
      used_tolerance = fragment_mass_tolerance / 1e6;
    }
    else
    {
      range = theoretical_spec.back().getMZ() - theoretical_spec[0].getMZ();
      used_tolerance = fragment_mass_tolerance;
    }

    // A priori probability of a random match given info about the theoretical spectrum
    double a_priori_p = 0;
    a_priori_p = 1 - pow(1 - 2 * used_tolerance / range,  static_cast<double>(theo_size));

    double log_occu_prob = 0;
    binomial flip(theo_size, a_priori_p);
    // min double number to avoid 0 values, causing scores with the value "inf"
    log_occu_prob = -log(cdf(complement(flip, matched_size)) + std::numeric_limits<double>::min());

    // score lower than 0 does not make sense, but can happen, if cfd = 0, then -log( 1 + <double>::min() ) < 0
    if (log_occu_prob >= 0.0)
    {
      return log_occu_prob;
    }
    else // underflow warning?
    {
      return 0;
    }
  }

  double XQuestScores::weightedTICScoreXQuest(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
  {
    // maxdigestlength and mindigestlength from standard settings of xQuest
    double maxdigestlength = 50;
    double mindigestlength = 5;
    if (!type_is_cross_link)
    {
      beta_size = ( maxdigestlength + mindigestlength ) - alpha_size;
    }

    double aatotal = alpha_size + beta_size;

    double invMax = 1 / (mindigestlength / (mindigestlength + maxdigestlength));
    double invFrac_alpha = 1 / (alpha_size / aatotal);
    double invFrac_beta = 1 / (beta_size / aatotal);
    double TIC_weight_alpha = invFrac_alpha / invMax;
    double TIC_weight_beta = invFrac_beta / invMax;

    double wTIC = TIC_weight_alpha * (intsum_alpha / total_current ) + TIC_weight_beta * (intsum_beta / total_current);
    return wTIC;
  }

  double XQuestScores::weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
  {
    if (!type_is_cross_link)
    {
      beta_size = alpha_size;
    }

    double aatotal = alpha_size + beta_size;

    // deviation from xQuest algorithm: invMax is not a constant anymore
    // and scales by the actual length difference between alpha and beta, rather than the maximal possible difference between any two peptides
    // this results in a local scaling, rather than a global one
    double invMax = 1 / (min(alpha_size, beta_size) / aatotal);

    double invFrac_alpha = 1 / (alpha_size / aatotal);
    double invFrac_beta = 1 / (beta_size / aatotal);
    double TIC_weight_alpha = invFrac_alpha / invMax;
    double TIC_weight_beta = invFrac_beta / invMax;

    double wTIC = TIC_weight_alpha * (intsum_alpha / total_current ) + TIC_weight_beta * (intsum_beta / total_current);
    return wTIC;
  }

  double XQuestScores::matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_linear, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks)
  {
    double intsum = 0;
    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_linear.size()); ++j)
    {
      intsum += spectrum_linear_peaks[matched_spec_linear[j].second].getIntensity();
    }
    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks.size()); ++j)
    {
      intsum += spectrum_xlink_peaks[matched_spec_xlinks[j].second].getIntensity();
    }
    return intsum;
  }

  double XQuestScores::totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_linear_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks)
  {
    // make vectors of matched peak indices
    double intsum(0);
    std::vector< Size > indices_linear;
    std::vector< Size > indices_xlinks;
    for (Size j = 0; j < matched_spec_linear_alpha.size(); ++j)
    {
      indices_linear.push_back(matched_spec_linear_alpha[j].second);
    }
    for (Size j = 0; j < matched_spec_linear_beta.size(); ++j)
    {
      indices_linear.push_back(matched_spec_linear_beta[j].second);
    }
    for (Size j = 0; j < matched_spec_xlinks_alpha.size(); ++j)
    {
      indices_xlinks.push_back(matched_spec_xlinks_alpha[j].second);
    }
    for (Size j = 0; j < matched_spec_xlinks_beta.size(); ++j)
    {
      indices_xlinks.push_back(matched_spec_xlinks_beta[j].second);
    }

    // make the indices in the vectors unique, to not sum up peak intensities multiple times
    sort(indices_linear.begin(), indices_linear.end());
    sort(indices_xlinks.begin(), indices_xlinks.end());
    std::vector< Size >::iterator last_unique_linear = unique(indices_linear.begin(), indices_linear.end());
    std::vector< Size >::iterator last_unique_xlinks = unique(indices_xlinks.begin(), indices_xlinks.end());
    indices_linear.erase(last_unique_linear, indices_linear.end());
    indices_xlinks.erase(last_unique_xlinks, indices_xlinks.end());

    // sum over intensities under the unique indices
    for (Size j = 0; j < indices_linear.size(); ++j)
    {
      intsum += spectrum_linear_peaks[indices_linear[j]].getIntensity();
    }
    for (Size j = 0; j < indices_xlinks.size(); ++j)
    {
      intsum += spectrum_xlink_peaks[indices_xlinks[j]].getIntensity();
    }
    return intsum;
  }

  std::vector< double > XQuestScores::xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance)
  {
    // generate vector of results, filled with zeroes
    std::vector< double > results(maxshift * 2 + 1, 0);

    // return 0 = no correlation, when one of the spectra is empty
    if (spec1.empty() || spec2.empty()) {
      return results;
    }

    double maxionsize = std::max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
    Int table_size = ceil(maxionsize / tolerance)+1;
    std::vector< double > ion_table1(table_size, 0);
    std::vector< double > ion_table2(table_size, 0);

    // Build tables of the same size, each bin has the size of the tolerance
    for (Size i = 0; i < spec1.size(); ++i)
    {
      Size pos = static_cast<Size>(ceil(spec1[i].getMZ() / tolerance));
      ion_table1[pos] = 10.0;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
      ion_table2[pos] = 10.0;
    }

    // Compute means
    double mean1 = (std::accumulate(ion_table1.begin(), ion_table1.end(), 0.0)) / table_size;
    double mean2 = (std::accumulate(ion_table2.begin(), ion_table2.end(), 0.0)) / table_size;

    // Compute denominator
    double s1 = 0;
    double s2 = 0;
    for (Int i = 0; i < table_size; ++i)
    {
      s1 += pow((ion_table1[i] - mean1), 2);
      s2 += pow((ion_table2[i] - mean2), 2);
    }
    double denom = sqrt(s1 * s2);

    // Calculate correlation for each shift
    for (Int shift = -maxshift; shift <= maxshift; ++shift)
    {
      double s = 0;
      for (Int i = 0; i < table_size; ++i)
      {
        Int j = i + shift;
        if ( (j >= 0) && (j < table_size))
        {
          s += (ion_table1[i] - mean1) * (ion_table2[j] - mean2);
        }
      }
      if (denom > 0)
      {
        results[shift + maxshift] = s / denom;
      }
    }
    return results;
  }

  double XQuestScores::xCorrelationPrescore(const PeakSpectrum & spec1, const PeakSpectrum & spec2, double tolerance)
  {
    // return 0 = no correlation, when one of the spectra is empty
    if (spec1.empty() || spec2.empty()) {
      return 0.0;
    }

    double maxionsize = std::max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
    Int table_size = ceil(maxionsize / tolerance)+1;
    std::vector< double > ion_table1(table_size, 0);
    std::vector< double > ion_table2(table_size, 0);

    // Build tables of the same size, each bin has the size of the tolerance
    for (Size i = 0; i < spec1.size(); ++i)
    {
      Size pos = static_cast<Size>(ceil(spec1[i].getMZ() / tolerance));
      ion_table1[pos] = 1;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
      ion_table2[pos] = 1;

    }

    double dot_product = 0.0;
    for (Size i = 0; i < ion_table1.size(); ++i)
    {
      dot_product += ion_table1[i] * ion_table2[i];
    }

    // determine the smaller spectrum and normalize by the number of peaks in it
    double peaks = std::min(spec1.size(), spec2.size());
    return dot_product / peaks;
  }

}
