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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
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

    // avoid 0 values in multiplication, adds a "dynamic range" among candidates with no matching common peaks to one of the peptides
    float matched_alpha_float = matched_alpha;
    if (matched_alpha <= 0)
    {
//      matched_alpha_float = std::numeric_limits<float>::min();
      matched_alpha_float = 0.1;
    }
    float matched_beta_float = matched_beta;
    if (matched_beta <= 0)
    {
//      matched_beta_float = std::numeric_limits<float>::min();
      matched_beta_float = 0.1;
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

    // Size matched_size = matched_spec.size();
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
    match_odds = -log(1 - cdf(flip, matched_size) + std::numeric_limits<double>::min());

    //     cout << "TEST a_priori_prob: " << a_priori_p << " | tolerance: " << tolerance_Th << " | theo_size: " << theo_size << " | matched_size: " << matched_size << " | cumul_binom: " << cumulativeBinomial_(theo_size, matched_size, a_priori_p)
    //              << " | match_odds: " << match_odds << endl;

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
      // log transform theoretical spectrum for more accurate a_priori_p estimation
      //
      // vector<double> log_theo_spec;
      // for (auto peak : theoretical_spec)
      // {
      //   log_theo_spec.push_back(std::log(peak.getMZ()));
      // }
      // range = log_theo_spec.back() - log_theo_spec[0];
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
    log_occu_prob = -log(1 - cdf(flip, matched_size) + std::numeric_limits<double>::min());

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
      // TODO what to do for mono-links, does this work?
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

  double XQuestScores::matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
  {
    double intsum = 0;
    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_common.size()); ++j)
    {
      intsum += spectrum_common_peaks[matched_spec_common[j].second].getIntensity();
    }
    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks.size()); ++j)
    {
      intsum += spectrum_xlink_peaks[matched_spec_xlinks[j].second].getIntensity();
    }
    return intsum;
  }

  double XQuestScores::totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
  {
    // make vectors of matched peak indices
    double intsum = 0;
    std::vector< Size > indices_common;
    std::vector< Size > indices_xlinks;
    for (Size j = 0; j < matched_spec_common_alpha.size(); ++j)
    {
      indices_common.push_back(matched_spec_common_alpha[j].second);
    }
    for (Size j = 0; j < matched_spec_common_beta.size(); ++j)
    {
      indices_common.push_back(matched_spec_common_beta[j].second);
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
    sort(indices_common.begin(), indices_common.end());
    sort(indices_xlinks.begin(), indices_xlinks.end());
    std::vector< Size >::iterator last_unique_common = unique(indices_common.begin(), indices_common.end());
    std::vector< Size >::iterator last_unique_xlinks = unique(indices_xlinks.begin(), indices_xlinks.end());
    indices_common.erase(last_unique_common, indices_common.end());
    indices_xlinks.erase(last_unique_xlinks, indices_xlinks.end());

    // sum over intensities under the unique indices
    for (Size j = 0; j < indices_common.size(); ++j)
    {
      intsum += spectrum_common_peaks[indices_common[j]].getIntensity();
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
    if (spec1.size() == 0 || spec2.size() == 0) {
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
      // with this line, use real intensities
//      ion_table1[pos] = spec1[i].getIntensity();
      // with this line, use intensities normalized to 10
      ion_table1[pos] = 10.0;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
      // with this line, use real intensities
//      ion_table2[pos] = spec2[i].getIntensity();
      // with this line, use intensities normalized to 10
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

}
