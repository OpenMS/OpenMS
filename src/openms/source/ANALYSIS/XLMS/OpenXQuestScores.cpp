// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/ANALYSIS/XLMS/OpenXQuestScores.h>
#include <boost/math/special_functions/binomial.hpp>
//#include <numeric>

using namespace std;

namespace OpenMS
{
  // fast pre-Score for cross-links
  // required: numbers of peaks for each chain, and how many of them were matched
  float OpenXQuestScores::preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta)
  {
    if ( (ionsAlpha > 0) && (ionsBeta > 0) )
    {
      float result = sqrt((static_cast<float>(matchedAlpha) / static_cast<float>(ionsAlpha)) * (static_cast<float>(matchedBeta) / static_cast<float>(ionsBeta)));
      return result;
    } else
    {
      return 0.0;
    }
  }

  // fast pre-Score for mono-links and loop-links
  float OpenXQuestScores::preScore(Size matchedAlpha, Size ionsAlpha)
  {
    if (ionsAlpha > 0)
    {
      float result = static_cast<float>(matchedAlpha) / static_cast<float>(ionsAlpha);
      return result;
    } else
    {
      return 0.0;
    }
  }

  // Statistics/Combinatorics functions for match-odds score
  // Standard cumulative binomial distribution
  double OpenXQuestScores::cumulativeBinomial(Size n, Size k, double p)
  {
    double p_cumul = 0.0;
    if (p < 1e-99) return static_cast<double>(k == 0); //  (not true/false, but 1/0 as probability)
    if (1 - p < 1e-99) return static_cast<double>(k != n); //
    if (k > n)  return 1.0;

    //cout << "TEST cumulBinom, passed if's, p = " << p << endl;

    for (Size j = 0; j < k; j++)
    {
      double coeff = boost::math::binomial_coefficient<double>(static_cast<unsigned int>(n), static_cast<unsigned int>(j));
      p_cumul += coeff * pow(p,  j) * pow((1-p), (n-j));
      //cout << "TEST coeff: " << coeff << " | first pow: " << pow(p,  j) << " | second pow: " <<  pow((1-p), (n-j)) << " | just added: " << coeff * pow(p,  j) * pow((1-p), (n-j)) << " | new p_cumul: " << p_cumul << endl;
    }

    // match-odds score becomes INFINITY for p_cumul >= 1, p_cumul might reach 1 because of insufficient precision, solved by using largest value smaller than 1
    if (p_cumul >= 1.0)
    {
      p_cumul = nexttoward(1.0, 0.0);
    }

    return p_cumul;
  }

  // match odds score, spectra must be sorted by position
  double OpenXQuestScores::match_odds_score(const RichPeakSpectrum& theoretical_spec,  const std::vector< std::pair< Size, Size > >& matched_spec, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, Size n_charges)
  {
    // TODO Theoretical spectra for cross-links contain 1. and 2. isotopic peak, is mostly one of them matched making theo_size = 2 * matched_size in the best case?
    // how does that skew the statistics?
    // should we use theo_size / 2 for cross-links?
    // or n_charges * 2?
    Size matched_size = matched_spec.size();
    Size theo_size = theoretical_spec.size();
    double range = theoretical_spec[theo_size-1].getMZ() -  theoretical_spec[0].getMZ();

    // Compute fragment tolerance for the middle of the range / mean of MZ values, if ppm
    // TODO mean should be used, so sum over MZs and devide by number, if ppm
    double mean = 0.0;
    for (Size i = 0; i < theo_size; ++i)
    {
      mean += theoretical_spec[i].getMZ();
    }
    mean = mean / theo_size;
    double tolerance_Th = fragment_mass_tolerance_unit_ppm ? mean * 1e-6 * fragment_mass_tolerance : fragment_mass_tolerance;

    // A priori probability of a random match given info about the theoretical spectrum
    //    double a_priori_p = a_priori_probability(tolerance_Th, theo_size, range, 3);
    double a_priori_p = 0;

    if (is_xlink_spectrum)
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  (theo_size / n_charges))));
    }
    else
    {
      a_priori_p = (1 - ( pow( (1 - 2 * tolerance_Th / (0.5 * range)),  (theo_size))));
    }

    double match_odds = 0;
    match_odds = -log(1 - cumulativeBinomial(theo_size, matched_size, a_priori_p) + 1e-5);

//     cout << "TEST a_priori_prob: " << a_priori_p << " | tolerance: " << tolerance_Th << " | theo_size: " << theo_size << " | matched_size: " << matched_size << " | cumul_binom: " << cumulativeBinomial(theo_size, matched_size, a_priori_p)
//              << " | match_odds: " << match_odds << endl;

    // score lower than 0 does not make sense, but can happen if cumBinom = 0, -log( 1 + 1e5 ) < 0
    if (match_odds >= 0.0)
    {
      return match_odds;
    }
    else
    {
      return 0;
    }
  }

//   // Cross-correlation, with shifting the second spectrum from -maxshift to +maxshift of tolerance bins (Tolerance in Da, a constant binsize)
//  template <typename SpectrumType1, typename SpectrumType2>
//  std::vector< double > OpenXQuestScores::xCorrelation(const SpectrumType1 & spec1, const SpectrumType2 & spec2, Int maxshift, double tolerance)
//  {
//    // generate vector of results, filled with zeroes
//    std::vector< double > results(maxshift * 2 + 1, 0);

//    // return 0 = no correlation, either positive nor negative, when one of the spectra is empty (e.g. when no common ions or xlink ions could be matched between light and heavy spectra)
//    if (spec1.size() == 0 || spec2.size() == 0) {
//      return results;
//    }

//    double maxionsize = max(spec1[spec1.size()-1].getMZ(), spec2[spec2.size()-1].getMZ());
//    Int table_size = ceil(maxionsize / tolerance)+1;
//    std::vector< double > ion_table1(table_size, 0);
//    std::vector< double > ion_table2(table_size, 0);

//    // Build tables of the same size, each bin has the size of the tolerance
//    for (Size i = 0; i < spec1.size(); ++i)
//    {
//      Size pos = static_cast<Size>(ceil(spec1[i].getMZ() / tolerance));
//      // TODO this line leads to using real intensities
////      ion_table1[pos] = spec1[i].getIntensity();
//      // TODO this line leads to using intensities normalized to 10
//      ion_table1[pos] = 10.0;
//    }
//    for (Size i = 0; i < spec2.size(); ++i)
//    {
//      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
//      // TODO this line leads to using real intensities
////      ion_table2[pos] = spec2[i].getIntensity();
//      // TODO this line leads to using intensities normalized to 10
//      ion_table2[pos] = 10.0;
//    }

//    // Compute means for real intensities
//    double mean1 = (std::accumulate(ion_table1.begin(), ion_table1.end(), 0.0)) / table_size;
//    double mean2 = (std::accumulate(ion_table2.begin(), ion_table2.end(), 0.0)) / table_size;

//    // Compute denominator
//    double s1 = 0;
//    double s2 = 0;
//    for (Int i = 0; i < table_size; ++i)
//    {
//      s1 += pow((ion_table1[i] - mean1), 2);
//      s2 += pow((ion_table2[i] - mean2), 2);
//    }
//    double denom = sqrt(s1 * s2);

//    // Calculate correlation for each shift
//    for (Int shift = -maxshift; shift <= maxshift; ++shift)
//    {
//      double s = 0;
//      for (Int i = 0; i < table_size; ++i)
//      {
//        Int j = i + shift;
//        if ( (j >= 0) && (j < table_size))
//        {
//          s += (ion_table1[i] - mean1) * (ion_table2[j] - mean2);
//        }
//      }
//      if (denom > 0)
//      {
//        results[shift + maxshift] = s / denom;
//      }
//    }
//    return results;
//  }

  // weigthed TIC score, using standard max- and mindigestlength, TODO remove digestlength from equation after benchmarking scores against xQuest?
  double OpenXQuestScores::weighted_TIC_score(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double intsum, double total_current, bool type_is_cross_link)
  {
    // TODO from xquest.def, but not used in this program aside from this calculation
    double maxdigestlength = 50;
    double mindigestlength = 5;
    if (!type_is_cross_link)
    {
      beta_size = ( maxdigestlength + mindigestlength ) - alpha_size;
      // this should already be the case
      intsum_beta = 0;
      intsum_alpha = intsum;
    }

    double aatotal = alpha_size + beta_size;

    // TODO maybe use this alternative version, based only on the lengths of the sequences?
    //double invMax = 1 / (min(alpha_size, beta_size) / aatotal);
    double invMax = 1 / (mindigestlength / (mindigestlength + maxdigestlength));
    double invFrac_alpha = 1 / (alpha_size / aatotal);
    double invFrac_beta = 1 / (beta_size / aatotal);
    double TIC_weight_alpha = invFrac_alpha / invMax;
    double TIC_weight_beta = invFrac_beta / invMax;

    double wTIC = TIC_weight_alpha * (intsum_alpha / total_current ) + TIC_weight_beta * (intsum_beta / total_current);
    return wTIC;
  }

//  // Sum of matched ion intesity, for Intsum score and %TIC score
//  double OpenXQuestScores::matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
//  {
//    double intsum = 0;
//    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_common.size()); ++j)
//    {
//      intsum += spectrum_common_peaks[matched_spec_common[j].second].getIntensity();
//    }
//    for (SignedSize j = 0; j < static_cast<SignedSize>(matched_spec_xlinks.size()); ++j)
//    {
//      intsum += spectrum_xlink_peaks[matched_spec_xlinks[j].second].getIntensity();
//    }
//    return intsum;
//  }

//  double OpenXQuestScores::total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
//  {
//    // make vectors of matched peak indices
//    double intsum = 0;
//    std::vector< Size > indices_common;
//    std::vector< Size > indices_xlinks;
//    for (Size j = 0; j < matched_spec_common_alpha.size(); ++j)
//    {
//      indices_common.push_back(matched_spec_common_alpha[j].second);
//    }
//    for (Size j = 0; j < matched_spec_common_beta.size(); ++j)
//    {
//      indices_common.push_back(matched_spec_common_beta[j].second);
//    }
//    for (Size j = 0; j < matched_spec_xlinks_alpha.size(); ++j)
//    {
//      indices_xlinks.push_back(matched_spec_xlinks_alpha[j].second);
//    }
//    for (Size j = 0; j < matched_spec_xlinks_beta.size(); ++j)
//    {
//      indices_xlinks.push_back(matched_spec_xlinks_beta[j].second);
//    }

//    // make the indices in the vectors unique
//    sort(indices_common.begin(), indices_common.end());
//    sort(indices_xlinks.begin(), indices_xlinks.end());
//    std::vector< Size >::iterator last_unique_common = unique(indices_common.begin(), indices_common.end());
//    std::vector< Size >::iterator last_unique_xlinks = unique(indices_xlinks.begin(), indices_xlinks.end());
//    indices_common.erase(last_unique_common, indices_common.end());
//    indices_xlinks.erase(last_unique_xlinks, indices_xlinks.end());

//    // sum over intensities under the unique indices
//    for (Size j = 0; j < indices_common.size(); ++j)
//    {
//      intsum += spectrum_common_peaks[indices_common[j]].getIntensity();
//    }
//    for (Size j = 0; j < indices_xlinks.size(); ++j)
//    {
//      intsum += spectrum_xlink_peaks[indices_xlinks[j]].getIntensity();
//    }

//    return intsum;
//  }


  // check whether the candidate pair is within the given tolerance to at least one precursor mass in the spectra data
  void filter_and_add_candidate (vector<OpenXQuestScores::XLPrecursor>& mass_to_candidates, vector< double >& spectrum_precursors, int charge_min, int charge_max, bool precursor_mass_tolerance_unit_ppm, double precursor_mass_tolerance, OpenXQuestScores::XLPrecursor precursor)
  {
    bool found_matching_precursors = false;
    // loop over all considered ion charges;
    // TODO: maybe precompute uncharged masses from precursor m/z values instead? don't forget to filter by charge then
    for (int charge = charge_min; charge <= charge_max; ++charge)
    {
      // if a precursor with the previous charge was found, there is no need to continue searching, stop loop
      if (found_matching_precursors)
      {
        break;
      }

      // use candidate mass and current charge to compute m/z
      double cross_link_mz = (precursor.precursor_mass + (static_cast<double>(charge) * Constants::PROTON_MASS_U)) / static_cast<double>(charge);

      vector< double >::const_iterator low_it;
      vector< double >::const_iterator up_it;

      // compute absolute tolerance from relative, if necessary
      double allowed_error = 0;
      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = cross_link_mz * precursor_mass_tolerance * 1e-6;
      }
      else // Dalton
      {
        allowed_error = precursor_mass_tolerance;
      }

      // find precursor with m/z >= low end of range
      low_it = lower_bound(spectrum_precursors.begin(), spectrum_precursors.end(), cross_link_mz - allowed_error);
      // find precursor with m/z > (not equal to) high end of range
      up_it =  upper_bound(spectrum_precursors.begin(), spectrum_precursors.end(), cross_link_mz + allowed_error);
      // if these two are equal, there is no precursor within the range

      if (low_it != up_it) // if they are not equal, there are matching precursors in the data
      {
        found_matching_precursors = true;
      }
    }

    // if precursors were found in the above for-loop, add candidate to results vector
    if (found_matching_precursors)
    {
// don't access this vector from two processing threads at the same time
#ifdef _OPENMP
#pragma omp critical
#endif
      mass_to_candidates.push_back(precursor);
    }
  }


  // Enumerate all pairs of peptides from the searched database and calculate their masses (inlcuding mono-links and loop-links)
  vector<OpenXQuestScores::XLPrecursor> OpenXQuestScores::enumerateCrossLinksAndMasses_(const vector<AASequence>&  peptide_AASeqs, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, int charge_min, int charge_max)
  {
    // initialize empty vector for the results
    vector<OpenXQuestScores::XLPrecursor> mass_to_candidates;
    // initialize progress counter
    Size countA = 0;

// Multithreading options: schedule: static, dynamic, guided
// use OpenMP to run this for-loop on multiple CPU cores
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (Size p1 = 0; p1 < peptide_AASeqs.size(); ++p1)
    {
      // get the amino acid sequence of this peptide as a character string
      String seq_first = peptide_AASeqs[p1].toUnmodifiedString();

      // every 500 peptides print current progress to console
      countA += 1;
      if (countA % 500 == 0)
      {
        cout << "Enumerating pairs with sequence " << countA << " of " << peptide_AASeqs.size() << ";\t Current pair count: " << mass_to_candidates.size() << endl;
      }

      // generate mono-links: one cross-linker with one peptide attached to one side
      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
      {
        // Monoisotopic weight of the peptide + cross-linker
        double cross_linked_pair_mass = peptide_AASeqs[p1].getMonoWeight() + cross_link_mass_mono_link[i];

        // Make sure it is clear only one peptide is considered here. Use NULL value for the second peptide.
        // to check: if(precursor.beta_index) returns "false" for NULL, "true" for any other value
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = NULL;

        // call function to compare with spectrum precursor masses
        // will only add this candidate, if the mass is within the given tolerance to any precursor in the spectra data
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, charge_min, charge_max, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }

      // test if this peptide could have loop-links: one cross-link with both sides attached to the same peptide
      // TODO check for distance between the two linked residues
      bool first_res = false; // is there a residue the first side of the linker can attach to?
      bool second_res = false; // is there a residue the second side of the linker can attach to?
      for (Size k = 0; k < seq_first.size()-1; ++k)
      {
        for (Size i = 0; i < cross_link_residue1.size(); ++i)
        {
          if (seq_first.substr(k, 1) == cross_link_residue1[i])
          {
            first_res = true;
          }
        }
        for (Size i = 0; i < cross_link_residue2.size(); ++i)
        {
          if (seq_first.substr(k, 1) == cross_link_residue2[i])
          {
            second_res = true;
          }
        }
      }

      // If both sides of a cross-linker can link to this peptide, generate the loop-link
      if (first_res && second_res)
      {
        // Monoisotopic weight of the peptide + cross-linker
        double cross_linked_pair_mass = peptide_AASeqs[p1].getMonoWeight() + cross_link_mass;

        // also only one peptide
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = NULL;

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, charge_min, charge_max, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }

      // Generate cross-links: one cross-linker linking two separate peptides, the most important case
      // Loop over all p2 peptide candidates, that come after p1 in the list
      for (Size p2 = p1; p2 < peptide_AASeqs.size(); ++p2)
      {
        // Monoisotopic weight of the first peptide + the second peptide + cross-linker
        double cross_linked_pair_mass = peptide_AASeqs[p1].getMonoWeight() + peptide_AASeqs[p2].getMonoWeight() + cross_link_mass;

        // this time both peptides have valid indices
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = p2;

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, charge_min, charge_max, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }
    }
    return mass_to_candidates;
  }

//  // Enumerates all possible combinations containing a cross-link, without specific cross-link positions. (There are cases where multiple positions are possible, but they have the same precursor mass)
//  // At this point the only difference between mono-links and loop-links is the added cross-link mass
//  multimap<double, pair<const AASequence*, const AASequence*> > OpenXQuestScores::enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2)
//  {
//    multimap<double, pair<const AASequence*, const AASequence*> > mass_to_candidates;
//    Size countA = 0;

//    vector<const StringView*> peptide_SVs;
//    vector<const AASequence*> peptide_AASeqs;
//    // preparing vectors compatible with openmp multi-threading, TODO this should only be a temporary fix (too much overhead?)
//    for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a)
//    {
//      peptide_SVs.push_back(&(a->first));
//      peptide_AASeqs.push_back(&(a->second));
//    }

//    //for (map<StringView, AASequence>::const_iterator a = peptides.begin(); a != peptides.end(); ++a) // old loop version

//// Multithreading options: schedule: static, dynamic, guided with chunk size
//#ifdef _OPENMP
//#pragma omp parallel for schedule(guided)
//#endif
//    for (Size p1 = 0; p1 < peptide_AASeqs.size(); ++p1)
//    {
//      String seq_first = peptide_AASeqs[p1]->toUnmodifiedString();

//      countA += 1;
//      if (countA % 500 == 0)
//      {
//        //LOG_DEBUG << "Enumerating pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << endl;
//        cout << "Enumerating pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << endl;

//      }

//      // generate mono-links
//      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
//      {
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + cross_link_mass_mono_link[i];
//        // Make sure it is clear this is a monolink, (is a NULL pointer a good idea?)
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], NULL)));
//      }

//      // generate loop-links
//      bool first_res = false;
//      bool second_res = false;
//      for (Size k = 0; k < seq_first.size()-1; ++k)
//      {
//        for (Size i = 0; i < cross_link_residue1.size(); ++i)
//        {
//          if (seq_first.substr(k, 1) == cross_link_residue1[i])
//          {
//            first_res = true;
//          }
//        }
//        for (Size i = 0; i < cross_link_residue2.size(); ++i)
//        {
//          if (seq_first.substr(k, 1) == cross_link_residue2[i])
//          {
//            second_res = true;
//          }
//        }
//      }
//      // If both sides of a homo- or heterobifunctional cross-linker can link to this peptide, generate the loop-link
//      if (first_res && second_res)
//      {
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + cross_link_mass;
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], NULL)));
//      }

//      // Generate cross-link between two peptides
//      //for (map<StringView, AASequence>::const_iterator b = a; b != peptides.end(); ++b)
//      for (Size p2 = p1; p2 < peptide_AASeqs.size(); ++p2)
//      {
//        // mass peptide1 + mass peptide2 + cross linker mass - cross link loss
//        double cross_linked_pair_mass = peptide_AASeqs[p1]->getMonoWeight() + peptide_AASeqs[p2]->getMonoWeight() + cross_link_mass;
//#ifdef _OPENMP
//#pragma omp critical
//#endif
//        mass_to_candidates.insert(make_pair(cross_linked_pair_mass, make_pair<const AASequence*, const AASequence*>(peptide_AASeqs[p1], peptide_AASeqs[p2])));
//      }
//    }

//    return mass_to_candidates;
//  }

}
