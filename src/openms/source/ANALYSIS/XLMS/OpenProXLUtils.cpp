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

#include <OpenMS/ANALYSIS/XLMS/OpenProXLUtils.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/Base64.h>
#include <boost/math/special_functions/binomial.hpp>
#include <limits>
#include <fstream>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
//#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
//#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

using namespace std;

namespace OpenMS
{
  // fast pre-Score for cross-links
  // required: numbers of peaks for each chain, and how many of them were matched
  float OpenProXLUtils::preScore(Size matchedAlpha, Size ionsAlpha, Size matchedBeta, Size ionsBeta)
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
  float OpenProXLUtils::preScore(Size matchedAlpha, Size ionsAlpha)
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
  double OpenProXLUtils::cumulativeBinomial(Size n, Size k, double p)
  {
    double p_cumul = 0.0;
    if (p < 1e-99) return static_cast<double>(k == 0); //  (not true/false, but 1/0 as probability)
    if (1 - p < 1e-99) return static_cast<double>(k != n); //
    if (k > n)  return 1.0;

    for (Size j = 0; j < k; ++j)
    {
      double coeff = 0;

      try
      {
        coeff = boost::math::binomial_coefficient<double>(static_cast<unsigned int>(n), static_cast<unsigned int>(j));
      }
      catch (std::overflow_error const& e)
      {
        cout << "Warning: Binomial coefficient for match-odds score has overflowed! Setting value to the maximal double value." << endl;
        cout << "binomial_coefficient was called with N = " << n << " and k = " << j << std::endl;
        coeff = std::numeric_limits<double>::max();
      }

      p_cumul += coeff * pow(p,  j) * pow((1-p), (n-j));
    }

    // match-odds score becomes INFINITY for p_cumul >= 1, p_cumul might reach 1 because of insufficient precision, solved by using largest value smaller than 1
    if (p_cumul >= 1.0)
    {
      p_cumul = nexttoward(1.0, 0.0);
    }

    return p_cumul;
  }

  // match odds score, spectra must be sorted by position
  double OpenProXLUtils::match_odds_score(const PeakSpectrum& theoretical_spec,  const std::vector< std::pair< Size, Size > >& matched_spec, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, Size n_charges)
  {
    // TODO Theoretical spectra for cross-links contain 1. and 2. isotopic peak, is mostly one of them matched making theo_size = 2 * matched_size in the best case?
    // how does that skew the statistics?
    // should we use theo_size / 2 for cross-links?
    // or n_charges * 2?
    Size matched_size = matched_spec.size();
    Size theo_size = theoretical_spec.size();

    if (matched_size < 1 || theo_size < 1)
    {
      return 0;
    }

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
    //recompile this shit fo shure

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
//  std::vector< double > OpenProXLUtils::xCorrelation(const SpectrumType1 & spec1, const SpectrumType2 & spec2, Int maxshift, double tolerance)
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
  double OpenProXLUtils::weighted_TIC_score(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double intsum, double total_current, bool type_is_cross_link)
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
//  double OpenProXLUtils::matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
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

//  double OpenProXLUtils::total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
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
  void filter_and_add_candidate (vector<OpenProXLUtils::XLPrecursor>& mass_to_candidates, vector< double >& spectrum_precursors, bool precursor_mass_tolerance_unit_ppm, double precursor_mass_tolerance, OpenProXLUtils::XLPrecursor precursor)
  {
    bool found_matching_precursors = false;
    // loop over all considered ion charges;
    // TODO: maybe precompute uncharged masses from precursor m/z values instead? don't forget to filter by charge then

    // use candidate mass and current charge to compute m/z
    //double cross_link_mz = (precursor.precursor_mass + (static_cast<double>(charge) * Constants::PROTON_MASS_U)) / static_cast<double>(charge);

    vector< double >::const_iterator low_it;
    vector< double >::const_iterator up_it;

    // compute absolute tolerance from relative, if necessary
    double allowed_error = 0;
    if (precursor_mass_tolerance_unit_ppm) // ppm
    {
      allowed_error = precursor.precursor_mass * precursor_mass_tolerance * 1e-6;
    }
    else // Dalton
    {
      allowed_error = precursor_mass_tolerance;
    }

    // find precursor with m/z >= low end of range
    low_it = lower_bound(spectrum_precursors.begin(), spectrum_precursors.end(), precursor.precursor_mass - allowed_error);
    // find precursor with m/z > (not equal to) high end of range
    up_it =  upper_bound(spectrum_precursors.begin(), spectrum_precursors.end(), precursor.precursor_mass + allowed_error);
    // if these two are equal, there is no precursor within the range

    if (low_it != up_it) // if they are not equal, there are matching precursors in the data
    {
      found_matching_precursors = true;
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
  vector<OpenProXLUtils::XLPrecursor> OpenProXLUtils::enumerateCrossLinksAndMasses_(const vector<OpenProXLUtils::PeptideMass>&  peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, vector< double >& spectrum_precursors, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
  {
    // initialize empty vector for the results
    vector<OpenProXLUtils::XLPrecursor> mass_to_candidates;
    // initialize progress counter
    Size countA = 0;

    double min_precursor = spectrum_precursors[0];
    double max_precursor = spectrum_precursors[spectrum_precursors.size()-1];

// Multithreading options: schedule: static, dynamic, guided
// use OpenMP to run this for-loop on multiple CPU cores
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize p1 = 0; p1 < static_cast<SignedSize>(peptides.size()); ++p1)
    {
      // get the amino acid sequence of this peptide as a character string
      String seq_first = peptides[p1].peptide_seq.toUnmodifiedString();

      // every 500 peptides print current progress to console
      countA += 1;
      if (countA % 500 == 0)
      {
        cout << "Enumerating pairs with sequence " << countA << " of " << peptides.size() << ";\t Current pair count: " << mass_to_candidates.size() << endl;
      }

      // generate mono-links: one cross-linker with one peptide attached to one side
      for (Size i = 0; i < cross_link_mass_mono_link.size(); i++)
      {
        // Monoisotopic weight of the peptide + cross-linker
        double cross_linked_pair_mass = peptides[p1].peptide_mass + cross_link_mass_mono_link[i];

        // Make sure it is clear only one peptide is considered here. Use NULL value for the second peptide.
        // to check: if(precursor.beta_index) returns "false" for NULL, "true" for any other value
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = NULL;

        // call function to compare with spectrum precursor masses
        // will only add this candidate, if the mass is within the given tolerance to any precursor in the spectra data
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
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
        double cross_linked_pair_mass = peptides[p1].peptide_mass + cross_link_mass;

        // also only one peptide
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = NULL;

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }

      // check for minimal mass of second peptide, jump farther than current peptide if possible
      double allowed_error = 0;
      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = min_precursor * precursor_mass_tolerance * 1e-6;
      }
      else // Dalton
      {
        allowed_error = precursor_mass_tolerance;
      }
      double min_second_peptide_mass = min_precursor - cross_link_mass - peptides[p1].peptide_mass - allowed_error;

      if (precursor_mass_tolerance_unit_ppm) // ppm
      {
        allowed_error = max_precursor * precursor_mass_tolerance * 1e-6;
      }
      double max_second_peptide_mass = max_precursor - cross_link_mass - peptides[p1].peptide_mass + allowed_error;

      // Generate cross-links: one cross-linker linking two separate peptides, the most important case
      // Loop over all p2 peptide candidates, that come after p1 in the list
      for (Size p2 = p1; p2 < peptides.size(); ++p2)
      {
        // skip peptides, that are too small in any case
        if (peptides[p2].peptide_mass < min_second_peptide_mass)
        {
          continue;
        } else if (peptides[p2].peptide_mass > max_second_peptide_mass)
        {
          break;
        }

        // Monoisotopic weight of the first peptide + the second peptide + cross-linker
        double cross_linked_pair_mass = peptides[p1].peptide_mass + peptides[p2].peptide_mass + cross_link_mass;

        // this time both peptides have valid indices
        XLPrecursor precursor;
        precursor.precursor_mass = cross_linked_pair_mass;
        precursor.alpha_index = p1;
        precursor.beta_index = p2;

        // call function to compare with spectrum precursor masses
        filter_and_add_candidate(mass_to_candidates, spectrum_precursors, precursor_mass_tolerance_unit_ppm, precursor_mass_tolerance, precursor);
      }
    }
    return mass_to_candidates;
  }

//  // Enumerates all possible combinations containing a cross-link, without specific cross-link positions. (There are cases where multiple positions are possible, but they have the same precursor mass)
//  // At this point the only difference between mono-links and loop-links is the added cross-link mass
//  multimap<double, pair<const AASequence*, const AASequence*> > OpenProXLUtils::enumerateCrossLinksAndMasses_(const multimap<StringView, AASequence>&  peptides, double cross_link_mass, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2)
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

 // Write xQuest.xml output
  void  OpenProXLUtils::writeXQuestXML(String out_file, String base_name, const std::vector< PeptideIdentification >& peptide_ids, const std::vector< std::vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra,
                                                String precursor_mass_tolerance_unit, String fragment_mass_tolerance_unit, double precursor_mass_tolerance, double fragment_mass_tolerance, double fragment_mass_tolerance_xlinks, String cross_link_name,
                                                double cross_link_mass_light, DoubleList cross_link_mass_mono_link, String in_fasta, String in_decoy_fasta, StringList cross_link_residue1, StringList cross_link_residue2, double cross_link_mass_iso_shift, String enzyme_name, Size missed_cleavages)
  {
    String spec_xml_name = base_name + "_matched";

    cout << "Writing xquest.xml to " << out_file << endl;
    ofstream xml_file;
    xml_file.open(out_file.c_str(), ios::trunc);
    // XML Header
    xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
    xml_file << "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>" << endl;

    // TODO!!! write actual experiment data
    // original date/time format: Fri Dec 18 12:28:23 2015
    DateTime time= DateTime::now();
    String timestring = time.getDate() + " " + time.getTime();

//    String precursor_mass_tolerance_unit = getStringOption_("precursor:mass_tolerance_unit");
//    String fragment_mass_tolerance_unit = getStringOption_("fragment:mass_tolerance_unit");
//    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
//    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
//    double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");

//    String cross_link_name = getStringOption_("cross_linker:name");
//    double cross_link_mass_light = getDoubleOption_("cross_linker:mass_light");
//    DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");
    String mono_masses;
    for (Size k = 0; k < cross_link_mass_mono_link.size()-1; ++k)
    {
      mono_masses += String(cross_link_mass_mono_link[k]) + ", ";
    }
    mono_masses += cross_link_mass_mono_link[cross_link_mass_mono_link.size()-1];

//    const string in_fasta(getStringOption_("database"));
//    const string in_decoy_fasta(getStringOption_("decoy_database"));
//    StringList cross_link_residue1 = getStringList_("cross_linker:residue1");
//    StringList cross_link_residue2 = getStringList_("cross_linker:residue2");
    String aarequired1, aarequired2;
    for (Size k= 0; k < cross_link_residue1.size()-1; ++k)
    {
      aarequired1 += cross_link_residue1[k] + ",";
    }
    aarequired1 += cross_link_residue1[cross_link_residue1.size()-1];
    for (Size k= 0; k < cross_link_residue2.size()-1; ++k)
    {
      aarequired2 += cross_link_residue2[k] + ",";
    }
    aarequired2 += cross_link_residue2[cross_link_residue2.size()-1];

//    double cross_link_mass_iso_shift = 0;
//    // This parameter is only available for the algorithm for labeled linkers
//    try
//    {
//      cross_link_mass_iso_shift = getDoubleOption_("cross_linker:mass_iso_shift");
//    }
//    catch (...)
//    {
//    }


//    String enzyme_name = getStringOption_("peptide:enzyme");
//    Size missed_cleavages = getIntOption_("peptide:missed_cleavages");

    xml_file << "<xquest_results xquest_version=\"OpenProXL 1.0\" date=\"" << timestring <<
             "\" author=\"Eugen Netz, Timo Sachsenberg\" tolerancemeasure_ms1=\"" << precursor_mass_tolerance_unit  <<
             "\" tolerancemeasure_ms2=\"" << fragment_mass_tolerance_unit << "\" ms1tolerance=\"" << precursor_mass_tolerance <<
             "\" ms2tolerance=\"" << fragment_mass_tolerance << "\" xlink_ms2tolerance=\"" << fragment_mass_tolerance_xlinks <<
             "\" crosslinkername=\"" << cross_link_name << "\" xlinkermw=\"" << cross_link_mass_light <<
             "\" monolinkmw=\"" << mono_masses << "\" database=\"" << in_fasta << "\" database_dc=\"" << in_decoy_fasta <<
             "\" xlinktypes=\"1111\" AArequired1=\"" << aarequired1 << "\" AArequired2=\"" << aarequired2 <<  "\" cp_isotopediff=\"" << cross_link_mass_iso_shift <<
             "\" enzyme_name=\"" << enzyme_name << "\" outputpath=\"" << spec_xml_name <<
             "\" Iontag_charges_for_index=\"1\" missed_cleavages=\"" << missed_cleavages <<
             "\" ntermxlinkable=\"0\" CID_match2ndisotope=\"1" <<
             "\" variable_mod=\"TODO\" nocutatxlink=\"1\" xcorrdelay=\"5\" >" << endl;



    for (vector< vector< CrossLinkSpectrumMatch > >::const_iterator top_csms_spectrum = all_top_csms.begin(); top_csms_spectrum != all_top_csms.end(); ++top_csms_spectrum)
    {
      vector< CrossLinkSpectrumMatch > top_vector = (*top_csms_spectrum);

      if (top_vector.empty())
      {
        continue;
      }
      // Spectrum Data, for each spectrum
      Size scan_index_light = top_vector[0].scan_index_light;
      Size scan_index_heavy = scan_index_light;
      if (cross_link_mass_iso_shift > 0)
      {
        scan_index_heavy = top_vector[0].scan_index_heavy;
      }
      const PeakSpectrum& spectrum_light = spectra[scan_index_light];
      double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();

      double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();
      double precursor_rt = spectrum_light.getRT();
      double precursor_mass = precursor_mz * static_cast<double>(precursor_charge) - static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U;

      double precursor_mz_heavy = spectra[scan_index_heavy].getPrecursors()[0].getMZ();
      double precursor_rt_heavy = spectra[scan_index_heavy].getRT();

      // print information about new peak to file (starts with <spectrum_search..., ends with </spectrum_search>
      String spectrum_light_name = base_name + ".light." + scan_index_light;
      String spectrum_heavy_name = base_name + ".heavy." + scan_index_heavy;

      String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;
      String rt_scans = String(precursor_rt) + ":" + String(precursor_rt_heavy);
      String mz_scans = String(precursor_mz) + ":" + String(precursor_mz_heavy);

      // Mean ion intensity (light spectrum, TODO add heavy spectrum?)
      double mean_intensity= 0;
      if (cross_link_mass_iso_shift > 0)
      {
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) mean_intensity += spectrum_light[j].getIntensity();
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectra[scan_index_heavy].size()); ++j) mean_intensity += spectra[scan_index_heavy][j].getIntensity();
        mean_intensity = mean_intensity / (spectrum_light.size() + spectra[scan_index_heavy].size());
      }
      else
      {
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum_light.size()); ++j) mean_intensity += spectrum_light[j].getIntensity();
        mean_intensity = mean_intensity / spectrum_light.size();
      }

      xml_file << "<spectrum_search spectrum=\"" << spectrum_name << "\" mean_ionintensity=\"" << mean_intensity << "\" ionintensity_stdev=\"" << "TODO" << "\" addedMass=\"" << "TODO" << "\" iontag_ncandidates=\"" << "TODO"
          << "\"  apriori_pmatch_common=\"" << "TODO" << "\" apriori_pmatch_xlink=\"" << "TODO" << "\" ncommonions=\"" << "TODO" << "\" nxlinkions=\"" << "TODO" << "\" mz_precursor=\"" << precursor_mz
          << "\" scantype=\"" << "light_heavy" << "\" charge_precursor=\"" << precursor_charge << "\" Mr_precursor=\"" << precursor_mass <<  "\" rtsecscans=\"" << rt_scans << "\" mzscans=\"" << mz_scans << "\" >" << endl;


      for (vector< CrossLinkSpectrumMatch>::const_iterator top_csm = top_csms_spectrum->begin(); top_csm != top_csms_spectrum->end(); ++top_csm)
      {
        String xltype = "monolink";
        String structure = top_csm->cross_link.alpha.toUnmodifiedString();
        String letter_first = structure.substr(top_csm->cross_link.cross_link_position.first, 1);


         // TODO track or otherwise find out, which kind of mono-link it was (if there are several possibilities for the weigths)
        double weight = top_csm->cross_link.alpha.getMonoWeight() + top_csm->cross_link.cross_linker_mass;
//          bool is_monolink = (top_csm->cross_link.cross_link_position.second == -1);
        int alpha_pos = top_csm->cross_link.cross_link_position.first + 1;
        int beta_pos = top_csm->cross_link.cross_link_position.second + 1;

        String topology = String("a") + alpha_pos;
        String id = structure + String("-") + letter_first + alpha_pos + String("-") + static_cast<int>(top_csm->cross_link.cross_linker_mass);

        if (top_csm->cross_link.getType() == ProteinProteinCrossLink::CROSS)
        {
          xltype = "xlink";
          structure += "-" + top_csm->cross_link.beta.toUnmodifiedString();
          topology += String("-b") + beta_pos;
          weight += top_csm->cross_link.beta.getMonoWeight();
          id = structure + "-" + topology;
        }
        else if (top_csm->cross_link.getType() == ProteinProteinCrossLink::LOOP)
        {
          xltype = "intralink";
          topology += String("-b") + beta_pos;
          String letter_second = structure.substr(top_csm->cross_link.cross_link_position.second, 1);
          id = structure + String("-") + letter_first + alpha_pos + String("-") + letter_second + beta_pos;
        }

         // Error calculation
        double cl_mz = (weight + (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / static_cast<double>(precursor_charge);
        double error = precursor_mz - cl_mz;
        double rel_error = (error / cl_mz) / 1e-6;

        PeptideIdentification pep_id = peptide_ids[top_csm->peptide_id_index];
        vector< PeptideHit > pep_hits = pep_id.getHits();

        String prot_alpha = pep_hits[0].getPeptideEvidences()[0].getProteinAccession();
        if (pep_hits[0].getPeptideEvidences().size() > 1)
        {
          for (Size i = 1; i < pep_hits[0].getPeptideEvidences().size(); ++i)
          {
            prot_alpha = prot_alpha + "," + pep_hits[0].getPeptideEvidences()[i].getProteinAccession();
          }
        }

        String prot_beta = "";

        if (pep_hits.size() > 1)
        {
          prot_beta= pep_hits[1].getPeptideEvidences()[0].getProteinAccession();
          if (pep_hits[1].getPeptideEvidences().size() > 1)
          {
            for (Size i = 1; i < pep_hits[1].getPeptideEvidences().size(); ++i)
            {
              prot_alpha = prot_alpha + "," + pep_hits[1].getPeptideEvidences()[i].getProteinAccession();
            }
          }
        }
        // Hit Data, for each cross-link to Spectrum Hit (e.g. top 5 per spectrum)
        xml_file << "<search_hit search_hit_rank=\"" <<top_csm->rank << "\" id=\"" << id << "\" type=\"" << xltype << "\" structure=\"" << structure << "\" seq1=\"" << top_csm->cross_link.alpha.toUnmodifiedString() << "\" seq2=\"" << top_csm->cross_link.beta.toUnmodifiedString()
              << "\" prot1=\"" << prot_alpha << "\" prot2=\"" << prot_beta << "\" topology=\"" << topology << "\" xlinkposition=\"" << (top_csm->cross_link.cross_link_position.first+1) << "," << (top_csm->cross_link.cross_link_position.second+1)
              << "\" Mr=\"" << weight << "\" mz=\"" << cl_mz << "\" charge=\"" << precursor_charge << "\" xlinkermass=\"" << top_csm->cross_link.cross_linker_mass << "\" measured_mass=\"" << precursor_mass << "\" error=\"" << error
              << "\" error_rel=\"" << rel_error << "\" xlinkions_matched=\"" << (top_csm->matched_xlink_alpha + top_csm->matched_xlink_beta) << "\" backboneions_matched=\"" << (top_csm->matched_common_alpha + top_csm->matched_common_beta)
              << "\" weighted_matchodds_mean=\"" << "TODO" << "\" weighted_matchodds_sum=\"" << "TODO" << "\" match_error_mean=\"" << "TODO" << "\" match_error_stdev=\"" << "TODO" << "\" xcorrx=\"" << top_csm->xcorrx_max << "\" xcorrb=\"" << top_csm->xcorrc_max << "\" match_odds=\"" <<top_csm->match_odds << "\" prescore=\"" << top_csm->pre_score
              << "\" prescore_alpha=\"" << "TODO" << "\" prescore_beta=\"" << "TODO" << "\" match_odds_alphacommon=\"" << "TODO" << "\" match_odds_betacommon=\"" << "TODO" << "\" match_odds_alphaxlink=\"" << "TODO"
              << "\" match_odds_betaxlink=\"" << "TODO" << "\" num_of_matched_ions_alpha=\"" << (top_csm->matched_common_alpha + top_csm->matched_xlink_alpha) << "\" num_of_matched_ions_beta=\"" << (top_csm->matched_common_beta + top_csm->matched_xlink_beta) << "\" num_of_matched_common_ions_alpha=\"" << top_csm->matched_common_alpha
              << "\" num_of_matched_common_ions_beta=\"" << top_csm->matched_common_beta << "\" num_of_matched_xlink_ions_alpha=\"" << top_csm->matched_xlink_alpha << "\" num_of_matched_xlink_ions_beta=\"" << top_csm->matched_xlink_beta << "\" xcorrall=\"" << "TODO" << "\" TIC=\"" << top_csm->percTIC
              << "\" TIC_alpha=\"" << "TODO" << "\" TIC_beta=\"" << "TODO" << "\" wTIC=\"" << top_csm->wTIC << "\" intsum=\"" << top_csm->int_sum * 100 << "\" apriori_match_probs=\"" << "TODO" << "\" apriori_match_probs_log=\"" << "TODO"
              << "\" HyperCommon=\"" << top_csm->HyperCommon << "\" HyperXLink=\"" << top_csm->HyperXlink << "\" HyperBoth=\"" << top_csm->HyperBoth << "\" PScoreCommon=\"" << top_csm->PScoreCommon << "\" PScoreXLink=\"" << top_csm->PScoreXlink << "\" PScoreBoth=\"" << top_csm->PScoreBoth
              << "\" series_score_mean=\"" << "TODO" << "\" annotated_spec=\"" << "" << "\" score=\"" << top_csm->score << "\" >" << endl;
        xml_file << "</search_hit>" << endl;
      }
      // Closing tag for Spectrum
      xml_file << "</spectrum_search>" << endl;
    }

    // Closing tag for results (end of file)
    xml_file << "</xquest_results>" << endl;
    xml_file.close();

    return;
  }

  PeakSpectrum OpenProXLUtils::mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum)
  {
    // merge peaks: create new spectrum, insert peaks from first and then from second spectrum
    PeakSpectrum resulting_spectrum;
    resulting_spectrum.insert(resulting_spectrum.end(), first_spectrum.begin(), first_spectrum.end());
    resulting_spectrum.insert(resulting_spectrum.end(), second_spectrum.begin(), second_spectrum.end());

    // merge DataArrays in a similar way
    for (Size i = 0; i < first_spectrum.getFloatDataArrays().size(); i++)
    {
      // TODO instead of this "if", get second array by name if available.  would not be dependent on order.
      if (second_spectrum.getFloatDataArrays().size() > i)
      {
        PeakSpectrum::FloatDataArray float_array;
        float_array.insert(float_array.end(), first_spectrum.getFloatDataArrays()[i].begin(), first_spectrum.getFloatDataArrays()[i].end());
        float_array.insert(float_array.end(), second_spectrum.getFloatDataArrays()[i].begin(), second_spectrum.getFloatDataArrays()[i].end());
        resulting_spectrum.getFloatDataArrays().push_back(float_array);
      }
    }

    for (Size i = 0; i < first_spectrum.getStringDataArrays().size(); i++)
    {
      if (second_spectrum.getStringDataArrays().size() > i)
      {
        PeakSpectrum::StringDataArray string_array;
        string_array.insert(string_array.end(), first_spectrum.getStringDataArrays()[i].begin(), first_spectrum.getStringDataArrays()[i].end());
        string_array.insert(string_array.end(), second_spectrum.getStringDataArrays()[i].begin(), second_spectrum.getStringDataArrays()[i].end());
        resulting_spectrum.getStringDataArrays().push_back(string_array);
      }
    }

    for (Size i = 0; i < first_spectrum.getIntegerDataArrays().size(); i++)
    {
      if (second_spectrum.getIntegerDataArrays().size() > i)
      {
        PeakSpectrum::IntegerDataArray integer_array;
        integer_array.insert(integer_array.end(), first_spectrum.getIntegerDataArrays()[i].begin(), first_spectrum.getIntegerDataArrays()[i].end());
        integer_array.insert(integer_array.end(), second_spectrum.getIntegerDataArrays()[i].begin(), second_spectrum.getIntegerDataArrays()[i].end());
        resulting_spectrum.getIntegerDataArrays().push_back(integer_array);
      }
    }

    // Spectra were simply concatenated, so they are not sorted by position anymore
    resulting_spectrum.sortByPosition();
    return resulting_spectrum;
  }

  void OpenProXLUtils::nLargestSpectrumFilter(PeakSpectrum spectrum, int peak_count)
  {
    if (spectrum.size() <= peak_count) return;

    // sort by reverse intensity
    spectrum.sortByIntensity(true);

    // keep the n largest peaks if more than n are present
    spectrum.resize(peak_count);

    // also resize DataArrays
    for (Size i = 0; i < spectrum.getFloatDataArrays().size(); i++)
    {
      spectrum.getFloatDataArrays()[i].resize(peak_count);
    }
    for (Size i = 0; i < spectrum.getStringDataArrays().size(); i++)
    {
      spectrum.getStringDataArrays()[i].resize(peak_count);
    }
    for (Size i = 0; i < spectrum.getIntegerDataArrays().size(); i++)
    {
      spectrum.getIntegerDataArrays()[i].resize(peak_count);
    }
  }

  void OpenProXLUtils::wrap_(const String& input, Size width, String & output)
  {
    Size start = 0;

    while (start + width < input.size())
    {
      output += input.substr(start, width) + "\n";
      start += width;
    }

    if (start < input.size())
    {
      output += input.substr(start, input.size() - start) + "\n";
    }
  }

  String OpenProXLUtils::getxQuestBase64EncodedSpectrum(const PeakSpectrum& spec, String header)
  {
    vector<String> in_strings;
    StringList sl;

    double precursor_mz = spec.getPrecursors()[0].getMZ();
    double precursor_z = spec.getPrecursors()[0].getCharge();

    // header lines
    if (!header.empty()) // common or xlinker spectrum will be reported
    {
      sl.push_back(header + "\n"); // e.g. GUA1372-S14-A-LRRK2_DSS_1A3.03873.03873.3.dta,GUA1372-S14-A-LRRK2_DSS_1A3.03863.03863.3.dta
      sl.push_back(String(precursor_mz) + "\n");
      sl.push_back(String(precursor_z) + "\n");
    }
    else // light or heavy spectrum will be reported
    {
      sl.push_back(String(precursor_mz) + "\t" + String(precursor_z) + "\n");
    }

    // write peaks
    for (Size i = 0; i != spec.size(); ++i)
    {
      String s;
      s += String(spec[i].getMZ()) + "\t";
      s += String(spec[i].getIntensity()) + "\t";

      // add fragment charge if meta value exists (must be present for 'common' and 'xlinker'.
//      if (spec[i].metaValueExists("z"))
//      {
//        s += String(spec[i].getMetaValue("z"));
//      }
      s += "0";

      s += "\n";

      sl.push_back(s);
    }

    String out;
    out.concatenate(sl.begin(), sl.end(), "");
    in_strings.push_back(out);

    String out_encoded;
    Base64().encodeStrings(in_strings, out_encoded, false, false);
    String out_wrapped;
    wrap_(out_encoded, 76, out_wrapped);
    return out_wrapped;
  }

  vector<ResidueModification> OpenProXLUtils::getModificationsFromStringList(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
  }

  void OpenProXLUtils::preprocessSpectraLabeled(PeakMap& exp, double fragment_mass_tolerance_xlinks, bool fragment_mass_tolerance_unit_ppm)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);
    // TODO perl code filters by dynamic range (1000), meaning everything below max_intensity / 1000 is filtered out additionally to 0 int, before scaling / normalizing

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);
    // TODO perl code scales to 0-100: int / max_int * 100

    // sort by rt
    exp.sortSpectra(false);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      // sort by mz and deisotope
      exp[exp_index].sortByPosition();
      exp[exp_index] = OpenProXLUtils::deisotopeAndSingleChargeMSSpectrum(exp[exp_index] , 1, 7, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm);
     }
  }

  void OpenProXLUtils::getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const PeakSpectrum & s1, const PeakSpectrum & s2, double tolerance, bool relative_tolerance, double intensity_cutoff)
  {
    if (!s1.isSorted() || !s2.isSorted())
    {
		throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
    }

    // clear result
    alignment.clear();

    // the weight of the intensity ratio term of the alignment score, and the gap penalty (since the ratio will always be between 0 and 1)
    // TODO this weight should somehow be proportional to the tolerance, since the tolerance penalty varies, maybe 1/2 of absolute tol?
    double intensity_weight = 0.1;

//    if (!(relative_tolerance))
//    {
      std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
      std::map<Size, std::map<Size, double> > matrix;

      // initialize to tolerance, will not change if tolerance is absolute (Da), updated for each position if tol is relative (ppm)
      double max_dist = tolerance;

      // init the matrix with "gap costs" tolerance + 1 for worst intensity ratio
      matrix[0][0] = 0;
      for (Size i = 1; i <= s1.size(); ++i)
      {
        // update relative max_dist at new position
        if (relative_tolerance)
        {
          max_dist = s1[i-1].getMZ() * tolerance * 1e-6;
        }
        matrix[i][0] = i * max_dist + i;
        traceback[i][0]  = std::make_pair(i - 1, 0);
      }
      for (Size j = 1; j <= s2.size(); ++j)
      {
        if (relative_tolerance)
        {
          max_dist = s2[j-1].getMZ() * tolerance * 1e-6;
        }
        matrix[0][j] = j * max_dist + j;
        traceback[0][j] = std::make_pair(0, j - 1);
      }

      // fill in the matrix
      Size left_ptr(1);
      Size last_i(0), last_j(0);

      //Size off_band_counter(0);
      for (Size i = 1; i <= s1.size(); ++i)
      {
        double pos1(s1[i - 1].getMZ());

        // update relative max_dist at new position
        if (relative_tolerance)
        {
          max_dist = pos1 * tolerance * 1e-6;
        }

        for (Size j = left_ptr; j <= s2.size(); ++j)
        {
          bool off_band(false);
          // find min of the three possible directions
          double pos2(s2[j - 1].getMZ());
          double diff_align = fabs(pos1 - pos2);

          // running off the right border of the band?
          if (pos2 > pos1 && diff_align >= max_dist)
          {
            if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
            {
              off_band = true;
            }
          }

          // can we tighten the left border of the band?
          if (pos1 > pos2 && diff_align >= max_dist && j > left_ptr + 1)
          {
            ++left_ptr;
          }

          double score_align = diff_align;

          if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j - 1) != matrix[i - 1].end())
          {
            score_align += matrix[i - 1][j - 1];
          }
          else
          {
            score_align += (i - 1 + j - 1) * max_dist + (i - 1 + j - 1) * intensity_weight;
          }

          double score_up = max_dist + intensity_weight;
          if (matrix.find(i) != matrix.end() && matrix[i].find(j - 1) != matrix[i].end())
          {
            score_up += matrix[i][j - 1];
          }
          else
          {
            score_up += (i + j - 1) * max_dist + (i + j - 1) / 10;
          }

          double score_left = max_dist + intensity_weight;
          if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j) != matrix[i - 1].end())
          {
            score_left += matrix[i - 1][j];
          }
          else
          {
            score_left += (i - 1 + j) * max_dist + (i - 1 + j) * intensity_weight;
          }

          // check for similar intensity values
          double intensity1(s1[i - 1].getIntensity());
          double intensity2(s2[j - 1].getIntensity());
          double int_ratio = min(intensity1, intensity2) / max(intensity1, intensity2);
          bool diff_int_clear = int_ratio > intensity_cutoff;

          // check for same charge (loose restriction, only excludes matches if both charges are known but unequal)
          bool charge_fits = true;
          if (s1.getIntegerDataArrays().size() > 0 && s2.getIntegerDataArrays().size() > 0)
          {
            int s1_charge = s1.getIntegerDataArrays()[0][i - 1];
            int s2_charge = s2.getIntegerDataArrays()[0][j - 1];
            charge_fits = s1_charge == s2_charge || s1_charge == 0 || s2_charge == 0;
//          LOG_DEBUG << "s1 charge: " << s1_charges[i - 1] << " | s2 charge: " << s2_charges[j - 1] << endl;
          }

          // int_ratio is between 0 and 1, multiply with intensity_weight for penalty
          score_align += (1 - int_ratio) * intensity_weight;

          if (score_align <= score_up && score_align <= score_left && diff_align < max_dist && diff_int_clear && charge_fits)
          {
//            cout << "Aligning peaks | score_align = " << score_align << "\t| int_ratio = " << int_ratio << "\t| score_up = " << score_up << "\t| score_left = " << score_left << endl;
            matrix[i][j] = score_align;
            traceback[i][j] = std::make_pair(i - 1, j - 1);
            last_i = i;
            last_j = j;
          }
          else
          {
            if (score_up <= score_left)
            {
              matrix[i][j] = score_up;
              traceback[i][j] = std::make_pair(i, j - 1);
            }
            else
            {
              matrix[i][j] = score_left;
              traceback[i][j] = std::make_pair(i - 1, j);
            }
          }

          if (off_band)
          {
            break;
          }
        }
      }

      // do traceback
      Size i = last_i;
      Size j = last_j;

      while (i >= 1 && j >= 1)
      {
        if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
        {
          alignment.push_back(std::make_pair(i - 1, j - 1));
        }
        Size new_i = traceback[i][j].first;
        Size new_j = traceback[i][j].second;

        i = new_i;
        j = new_j;
      }

      std::reverse(alignment.begin(), alignment.end());
  }

    PeakSpectrum OpenProXLUtils::deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped, Size min_isopeaks, Size max_isopeaks, bool make_single_charged)
    {
      PeakSpectrum out;
      PeakSpectrum::IntegerDataArray charge_array;

      vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
      vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
      if (old_spectrum.empty())
      {
        return out;
      }

      // determine charge seeds and extend them
      vector<Int> features(old_spectrum.size(), -1);
      Int feature_number = 0;

      for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
      {
        double current_mz = old_spectrum[current_peak].getMZ();
        mono_iso_peak_intensity[current_peak] = old_spectrum[current_peak].getIntensity();

        for (Int q = max_charge; q >= min_charge; --q)   // important: test charge hypothesis from high to low
        {
          // try to extend isotopes from mono-isotopic peak
          // if extension larger then min_isopeaks possible:
          //   - save charge q in mono_isotopic_peak[]
          //   - annotate all isotopic peaks with feature number
          if (features[current_peak] == -1)   // only process peaks which have no assigned feature number
          {
            bool has_min_isopeaks = true;
            vector<Size> extensions;
            for (Size i = 0; i < max_isopeaks; ++i)
            {
              double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
              Size p = old_spectrum.findNearest(expected_mz);
              double tolerance_dalton = fragment_tolerance_unit_ppm ? fragment_tolerance * old_spectrum[p].getMZ() * 1e-6 : fragment_tolerance;
              if (fabs(old_spectrum[p].getMZ() - expected_mz) > tolerance_dalton)   // test for missing peak
              {
                if (i < min_isopeaks)
                {
                  has_min_isopeaks = false;
                }
                break;
              }
              else
              {
                // TODO: include proper averagine model filtering. assuming the intensity gets lower for heavier peaks does not work for the high masses of cross-linked peptides
//                Size n_extensions = extensions.size();
//                if (n_extensions != 0)
//                {
//                  if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
//                  {
//                    if (i < min_isopeaks)
//                    {
//                      has_min_isopeaks = false;
//                    }
//                    break;
//                  }
//                }

                // averagine check passed
                extensions.push_back(p);
                mono_iso_peak_intensity[current_peak] += old_spectrum[p].getIntensity();
              }
            }

            if (has_min_isopeaks)
            {
              //LOG_DEBUG << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
              mono_isotopic_peak[current_peak] = q;
              for (Size i = 0; i != extensions.size(); ++i)
              {
                features[extensions[i]] = feature_number;
              }
              feature_number++;
            }
          }
        }
      }


      // creating PeakSpectrum containing charges
      //out.clear(false);

      for (Size i = 0; i != old_spectrum.size(); ++i)
      {
        Int z = mono_isotopic_peak[i];
        if (keep_only_deisotoped)
        {
          if (z == 0)
          {
            continue;
          }

          // if already single charged or no decharging selected keep peak as it is
          if (!make_single_charged)
          {
            RichPeak1D p;
            p.setMZ(old_spectrum[i].getMZ());
//            p.setIntensity(old_spectrum[i].getIntensity());
            p.setIntensity(mono_iso_peak_intensity[i]);
//            p.setMetaValue("z", z);
            charge_array.push_back(z);
            out.push_back(p);
          }
          else
          {
            RichPeak1D p;
//            p.setIntensity(old_spectrum[i].getIntensity());
            p.setIntensity(mono_iso_peak_intensity[i]);
            p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
//            p.setMetaValue("z", 1);
            charge_array.push_back(1);
            out.push_back(p);
          }
        }
        else
        {
          // keep all unassigned peaks
          if (features[i] < 0)
          {
            RichPeak1D p;
            p.setMZ(old_spectrum[i].getMZ());
            p.setIntensity(old_spectrum[i].getIntensity());
//            p.setMetaValue("z", 0);
            charge_array.push_back(0);
            out.push_back(p);
            continue;
          }

          // convert mono-isotopic peak with charge assigned by deisotoping
          if (z != 0)
          {
            if (!make_single_charged)
            {
              RichPeak1D p;
              p.setMZ(old_spectrum[i].getMZ());
              p.setIntensity(mono_iso_peak_intensity[i]);
//              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
              charge_array.push_back(z);
              out.push_back(p);
            }
            else
            {
              RichPeak1D p;
              p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
              p.setIntensity(mono_iso_peak_intensity[i]);
//              p.setIntensity(old_spectrum[i].getIntensity());
//              p.setMetaValue("z", z);
              charge_array.push_back(z);
              out.push_back(p);
            }
          }
        }
      }
      out.setPrecursors(old_spectrum.getPrecursors());
      out.setRT(old_spectrum.getRT());

      out.setNativeID(old_spectrum.getNativeID());
      out.setInstrumentSettings(old_spectrum.getInstrumentSettings());
      out.setAcquisitionInfo(old_spectrum.getAcquisitionInfo());
      out.setSourceFile(old_spectrum.getSourceFile());
      out.setDataProcessing(old_spectrum.getDataProcessing());
      out.setType(old_spectrum.getType());
      out.setMSLevel(old_spectrum.getMSLevel());
      out.setName(old_spectrum.getName());

      out.getIntegerDataArrays().push_back(charge_array);

//      out.sortByPosition();
      return out;
    }

  std::vector<OpenProXLUtils::PeptideMass> OpenProXLUtils::digestDatabase(vector<FASTAFile::FASTAEntry> fasta_db, EnzymaticDigestion digestor, Size min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, std::vector<ResidueModification> fixed_modifications, std::vector<ResidueModification> variable_modifications, Size max_variable_mods_per_peptide, Size count_proteins, Size count_peptides, bool n_term_linker, bool c_term_linker)
  {
    multimap<StringView, AASequence> processed_peptides;
    vector<OpenProXLUtils::PeptideMass> peptide_masses;
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
    // digest and filter database
    for (SignedSize fasta_index = 0; fasta_index < static_cast<SignedSize>(fasta_db.size()); ++fasta_index)
    {
//#ifdef _OPENMP
//#pragma omp atomic
//#endif
      ++count_proteins;

      // store vector of substrings pointing in fasta database (bounded by pairs of begin, end iterators)
      vector<StringView> current_digest;
      digestor.digestUnmodifiedString(fasta_db[fasta_index].sequence, current_digest, min_peptide_length);

      for (vector<StringView>::iterator cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        // skip peptides with invalid AAs // TODO is this necessary?
        if (cit->getString().has('B') || cit->getString().has('O') || cit->getString().has('U') || cit->getString().has('X') || cit->getString().has('Z')) continue;

        OpenProXLUtils::PeptidePosition position = OpenProXLUtils::INTERNAL;
        if (fasta_db[fasta_index].sequence.hasPrefix(cit->getString()))
        {
          position = OpenProXLUtils::N_TERM;
        } else if (fasta_db[fasta_index].sequence.hasSuffix(cit->getString()))
        {
          position = OpenProXLUtils::C_TERM;
        }

        // skip if no cross-linked residue
        bool skip = true;
        for (Size k = 0; k < cross_link_residue1.size(); k++)
        {
          if (cit->getString().find(cross_link_residue1[k]) < cit->getString().size()-1)
          {
            skip = false;
          }
          if (n_term_linker && position == OpenProXLUtils::N_TERM)
          {
            skip = false;
          }
          if (c_term_linker && position == OpenProXLUtils::C_TERM)
          {
            skip = false;
          }
        }
        for (Size k = 0; k < cross_link_residue2.size(); k++)
        {
          if (cit->getString().find(cross_link_residue2[k]) < cit->getString().size()-1)
          {
            skip = false;
          }
          if (n_term_linker && position == OpenProXLUtils::N_TERM)
          {
            skip = false;
          }
          if (c_term_linker && position == OpenProXLUtils::C_TERM)
          {
            skip = false;
          }
        }
        if (skip) continue;

        bool already_processed = false;
//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
        {
          if (processed_peptides.find(*cit) != processed_peptides.end())
          {
            // peptide (and all modified variants) already processed so skip it
            already_processed = true;
          }
        }

        if (already_processed)
        {
          continue;
        }
//        if (cit->getString().find('K') >= cit->getString().size()-1)
//        {
//          continue;
//        }



//#ifdef _OPENMP
//#pragma omp atomic
//#endif
        ++count_peptides;

        vector<AASequence> all_modified_peptides;

        // generate all modified variants of a peptide
        // Note: no critial section is needed despite ResidueDB not beeing thread sage.
        //       It is only written to on introduction of novel modified residues. These residues have been already added above (single thread context).
        {
          AASequence aas = AASequence::fromString(cit->getString());
          ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications.begin(), fixed_modifications.end(), aas);
          ModifiedPeptideGenerator::applyVariableModifications(variable_modifications.begin(), variable_modifications.end(), aas, max_variable_mods_per_peptide, all_modified_peptides);
        }

        for (SignedSize mod_pep_idx = 0; mod_pep_idx < static_cast<SignedSize>(all_modified_peptides.size()); ++mod_pep_idx)
        {
          const AASequence& candidate = all_modified_peptides[mod_pep_idx];
          OpenProXLUtils::PeptideMass pep_mass;
          pep_mass.peptide_mass = candidate.getMonoWeight();
          pep_mass.peptide_seq = candidate;
          pep_mass.position = position;

//#ifdef _OPENMP
//#pragma omp critical (processed_peptides_access)
//#endif
          {
            processed_peptides.insert(pair<StringView, AASequence>(*cit, candidate));
            peptide_masses.push_back(pep_mass);
          }
        }
      }
    }
    sort(peptide_masses.begin(), peptide_masses.end());
    return peptide_masses;
  }

  vector <ProteinProteinCrossLink> OpenProXLUtils::buildCandidates(const std::vector< OpenProXLUtils::XLPrecursor > & candidates, const std::vector<OpenProXLUtils::PeptideMass> & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, const DoubleList & cross_link_mass_mono_link, double precursor_mass, double allowed_error, String cross_link_name, bool n_term_linker, bool c_term_linker)
  {
    vector <ProteinProteinCrossLink> cross_link_candidates;
    for (Size i = 0; i < candidates.size(); ++i)
    {
      OpenProXLUtils::XLPrecursor candidate = candidates[i];
      vector <SignedSize> link_pos_first;
      vector <SignedSize> link_pos_second;
      AASequence peptide_first = peptide_masses[candidate.alpha_index].peptide_seq;
      OpenProXLUtils::PeptidePosition peptide_pos_first = peptide_masses[candidate.alpha_index].position;
      AASequence peptide_second;
      OpenProXLUtils::PeptidePosition peptide_pos_second = OpenProXLUtils::INTERNAL;
      if (candidate.beta_index)
      {
        peptide_second = peptide_masses[candidate.beta_index].peptide_seq;
        peptide_pos_second = peptide_masses[candidate.beta_index].position;
      }
      String seq_first = peptide_first.toUnmodifiedString();
      String seq_second =  peptide_second.toUnmodifiedString();

      // TODO mono-links and loop-links with different masses can be generated for the same precursor mass, but only one of them can be valid each time.
      // Find out which is the case. But it should not happen often enough to slow down the tool significantly.
      bool is_loop = abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass)) <= allowed_error;

      for (Size k = 0; k < seq_first.size()-1; ++k)
      {
        for (Size x = 0; x < cross_link_residue1.size(); ++x)
        {
          if (seq_first.substr(k, 1) == cross_link_residue1[x]) link_pos_first.push_back(k);
        }
      }
      if (candidate.beta_index)
      {
        for (Size k = 0; k < seq_second.size()-1; ++k)
        {
          for (Size x = 0; x < cross_link_residue2.size(); ++x)
          {
            if (seq_second.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
          }
        }
      } else
      {
        // Second position defining a mono-link and the second positions on the same peptide for loop links (only one of these two is valid for any specific precursor)
        if (!is_loop)
        {
          link_pos_second.push_back(-1);
        }
        else
        {
          for (Size k = 0; k < seq_first.size()-1; ++k)
          {
            for (Size x = 0; x < cross_link_residue2.size(); ++x)
            {
              if (seq_first.substr(k, 1) == cross_link_residue2[x]) link_pos_second.push_back(k);
            }
          }
        }
      }

        // Determine larger peptide (alpha) by sequence length, use mass as tie breaker
      bool alpha_first = true;

      if (seq_second.size() > seq_first.size())
      {
        alpha_first = false;
      } else if (seq_second.size() == seq_first.size() && peptide_second.getMonoWeight() > peptide_first.getMonoWeight())
      {
        alpha_first = false;
      }

      // TODO remodel this, there should be a simpler way, e.g. the peptides were sorted so "second" is always heavier?
      // generate cross_links for all valid combinations
      for (Size x = 0; x < link_pos_first.size(); ++x)
      {
        for (Size y = 0; y < link_pos_second.size(); ++y)
        {
          ProteinProteinCrossLink cross_link_candidate;
          cross_link_candidate.cross_linker_name = cross_link_name;
          // if loop link, and the positions are the same, then it is linking the same residue with itself,  skip this combination, also pos1 > pos2 would be the same link as pos1 < pos2
          if (((seq_second.size() == 0) && (link_pos_first[x] >= link_pos_second[y])) && (link_pos_second[y] != -1))
          {
            continue;
          }
          if (alpha_first)
          {
            cross_link_candidate.alpha = peptide_first;
            cross_link_candidate.beta = peptide_second;
            cross_link_candidate.cross_link_position.first = link_pos_first[x];
            cross_link_candidate.cross_link_position.second = link_pos_second[y];
            cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
            cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
          }
          else
          {
            cross_link_candidate.alpha = peptide_second;
            cross_link_candidate.beta = peptide_first;
            cross_link_candidate.cross_link_position.first = link_pos_second[y];
            cross_link_candidate.cross_link_position.second = link_pos_first[x];
            cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
            cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
          }
          // Cross-linker mass is only one of the mono-link masses, if there is no second position (second == -1), otherwise the normal linker mass
          if (link_pos_second[y] != -1)
          {
            cross_link_candidate.cross_linker_mass = cross_link_mass;
            cross_link_candidates.push_back(cross_link_candidate);
          }
          else
          {
            for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
            {
              // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
              if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
              {
                cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];;
                cross_link_candidates.push_back(cross_link_candidate);
              }
            }
          }
        }
      }

      if (peptide_pos_second != OpenProXLUtils::INTERNAL)
      {
        ResidueModification::TermSpecificity second_spec;
        Size mod_pos;
        bool compatible = false;
        if (n_term_linker && (peptide_pos_second == OpenProXLUtils::N_TERM))
        {
          second_spec = ResidueModification::N_TERM;
          mod_pos = 0;
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_second == OpenProXLUtils::C_TERM))
        {
          second_spec = ResidueModification::C_TERM;
          mod_pos = peptide_second.size()-1;
          compatible = true;
        }
        if (compatible)
        {
          for (Size x = 0; x < link_pos_first.size(); ++x)
          {
            ProteinProteinCrossLink cross_link_candidate;
            if (alpha_first)
            {
              cross_link_candidate.alpha = peptide_first;
              cross_link_candidate.beta = peptide_second;
              cross_link_candidate.cross_link_position.first = link_pos_first[x];
              cross_link_candidate.cross_link_position.second = mod_pos;
              cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;
              cross_link_candidate.term_spec_beta = second_spec;
            }
            else
            {
              cross_link_candidate.alpha = peptide_second;
              cross_link_candidate.beta = peptide_first;
              cross_link_candidate.cross_link_position.first = mod_pos;
              cross_link_candidate.cross_link_position.second = link_pos_first[x];
              cross_link_candidate.term_spec_alpha = second_spec;
              cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;
            }
            // If second peptide has a term specificity, there must be a second peptide, so we don't have to consider mono or loop-links
            cross_link_candidate.cross_linker_mass = cross_link_mass;
            cross_link_candidate.cross_linker_name = cross_link_name;
            cross_link_candidates.push_back(cross_link_candidate);

          }
        }
      }

      if (peptide_pos_first != OpenProXLUtils::INTERNAL)
      {
        ResidueModification::TermSpecificity first_spec;
        Size mod_pos;
        bool compatible = false;
        if (n_term_linker && (peptide_pos_first == OpenProXLUtils::N_TERM))
        {
          first_spec = ResidueModification::N_TERM;
          mod_pos = 0;
          compatible = true;
        }
        if (c_term_linker && (peptide_pos_first == OpenProXLUtils::C_TERM))
        {
          first_spec = ResidueModification::C_TERM;
          mod_pos = peptide_first.size()-1;
          compatible = true;
        }
        if (compatible)
        {
          for (Size x = 0; x < link_pos_second.size(); ++x)
          {
            ProteinProteinCrossLink cross_link_candidate;
            cross_link_candidate.cross_linker_name = cross_link_name;
            if (alpha_first)
            {
              cross_link_candidate.alpha = peptide_first;
              cross_link_candidate.beta = peptide_second;
              cross_link_candidate.cross_link_position.first = mod_pos;
              cross_link_candidate.cross_link_position.second = link_pos_second[x];
              cross_link_candidate.term_spec_alpha = first_spec;
              cross_link_candidate.term_spec_beta = ResidueModification::ANYWHERE;;
            }
            else
            {
              cross_link_candidate.alpha = peptide_second;
              cross_link_candidate.beta = peptide_first;
              cross_link_candidate.cross_link_position.first = link_pos_second[x];
              cross_link_candidate.cross_link_position.second = mod_pos;
              cross_link_candidate.term_spec_alpha = ResidueModification::ANYWHERE;;
              cross_link_candidate.term_spec_beta = first_spec;
            }
            // Cross-linker mass is only one of the mono-link masses, if there is no second position (second == -1), otherwise the normal linker mass
            if (link_pos_second[x] != -1)
            {
              cross_link_candidate.cross_linker_mass = cross_link_mass;
              cross_link_candidates.push_back(cross_link_candidate);
            }
            else
            {
              for (Size k = 0; k < cross_link_mass_mono_link.size(); ++k)
              {
                // only use the correct mono-links (at this point we know it is a mono-link, but not which one)
                if (abs(precursor_mass - (peptide_first.getMonoWeight() + cross_link_mass_mono_link[k])) <= allowed_error)
                {
                  cross_link_candidate.cross_linker_mass = cross_link_mass_mono_link[k];
                  cross_link_candidates.push_back(cross_link_candidate);
                }
              }
            }
          }
        }
      }
    }
    return cross_link_candidates;
  }

  void OpenProXLUtils::buildFragmentAnnotations(std::vector<PeptideHit::FragmentAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum)
  {
    if (theoretical_spectrum.empty() || experiment_spectrum.empty())
    {
      return;
    }
    PeakSpectrum::IntegerDataArray charges = theoretical_spectrum.getIntegerDataArrays()[0];
    PeakSpectrum::StringDataArray names = theoretical_spectrum.getStringDataArrays()[0];
    for (Size k = 0; k < matching.size(); ++k)
    {
      PeptideHit::FragmentAnnotation frag_anno;
      frag_anno.mz = experiment_spectrum[matching[k].second].getMZ();
      frag_anno.intensity = experiment_spectrum[matching[k].second].getIntensity();

      frag_anno.charge = charges[matching[k].first];
      frag_anno.annotation = names[matching[k].first];
      frag_annotations.push_back(frag_anno);
    }
  }

  void OpenProXLUtils::buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy)
  {
    for (Size i = 0; i < top_csms_spectrum.size(); ++i)
    {
      PeptideIdentification peptide_id;

      const PeakSpectrum& spectrum_light = spectra[scan_index];
      double precursor_charge = spectrum_light.getPrecursors()[0].getCharge();
      double precursor_mz = spectrum_light.getPrecursors()[0].getMZ();

      String xltype = "cross-link";
      SignedSize alpha_pos = top_csms_spectrum[i].cross_link.cross_link_position.first;
      SignedSize beta_pos = top_csms_spectrum[i].cross_link.cross_link_position.second;

      if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::MONO)
      {
        xltype = "mono-link";
      }
      else if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::LOOP)
      {
        xltype = "loop-link";
      }

      PeptideHit ph_alpha, ph_beta;
      // Set monolink as a modification or add MetaValue for cross-link identity and mass
      AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
      ResidueModification::TermSpecificity alpha_term_spec = top_csms_spectrum[i].cross_link.term_spec_alpha;
      if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::MONO)
      {
        //AASequence seq_alpha = top_csms_spectrum[i].cross_link.alpha;
        vector< String > mods;
        const String residue = seq_alpha[alpha_pos].getOneLetterCode();
        LOG_DEBUG << "Searching mono-link for " << residue << " | " << alpha_pos << endl;
        ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001, residue, ResidueModification::ANYWHERE);
        LOG_DEBUG << "number of modifications fitting the diff mass: " << mods.size() << endl;
        bool mod_set = false;
        if (mods.size() > 0) // If several mods have the same diff mass, try to resolve ambiguity by cross-linker name (e.g. DSS and BS3 are different reagents, but have the same result after the reaction)
        {
          for (Size s = 0; s < mods.size(); ++s)
          {
            if (mods[s].hasSubstring(top_csms_spectrum[i].cross_link.cross_linker_name))
            {
              LOG_DEBUG << "applied modification: " << mods[s] << endl;
              seq_alpha.setModification(alpha_pos, mods[s]);
              mod_set = true;
              break;
            }
          }
        }
        else if (mods.size() == 0 && (alpha_pos == 0 || alpha_pos == seq_alpha.size()-1))
        {
          LOG_DEBUG << "No residue specific mono-link found, searching for terminal mods..." << endl;
          ModificationsDB::getInstance()->searchModificationsByDiffMonoMass(mods, top_csms_spectrum[i].cross_link.cross_linker_mass, 0.001, "", alpha_term_spec);
          if (mods.size() > 0)
          {
            Size mod_index = 0;
            for (Size s = 0; s < mods.size(); ++s)
            {
              if (mods[s].hasSubstring(top_csms_spectrum[i].cross_link.cross_linker_name))
              {
                mod_index = s;
              }
            }
//            cout << "Terminal Mod Test; mod: " << mod_index <<  " | " << mods[mod_index] << " | term_spec: " << alpha_term_spec << endl;
            if (alpha_term_spec == ResidueModification::N_TERM)
            {
              LOG_DEBUG << "Setting N-term mono-link: " << mods[mod_index] << endl;
              seq_alpha.setNTerminalModification(mods[mod_index]);
            }
            else
            {
              LOG_DEBUG << "Setting C-term mono-link: " << mods[mod_index] << endl;
              seq_alpha.setCTerminalModification(mods[mod_index]);
            }
            mod_set = true;
          }
        }

        if ( (mods.size() > 0) && (!mod_set) ) // If resolving by name did not work, use any with matching diff mass
        {
          seq_alpha.setModification(alpha_pos, mods[0]);
          mod_set = true;
        }
        if (!mod_set) // If no equivalent mono-link exists in the UNIMOD or XLMOD databases, use the given name to construct a placeholder
        {
          String mod_name = String("unknown mono-link " + top_csms_spectrum[i].cross_link.cross_linker_name + " mass " + String(top_csms_spectrum[i].cross_link.cross_linker_mass));
          //seq_alpha.setModification(alpha_pos, mod_name);
          LOG_DEBUG << "unknown mono-link" << endl;
          ph_alpha.setMetaValue("xl_mod", mod_name);
          ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
        }
      }
      else
      {
        ph_alpha.setMetaValue("xl_mod", top_csms_spectrum[i].cross_link.cross_linker_name);
        ph_alpha.setMetaValue("xl_mass", DataValue(top_csms_spectrum[i].cross_link.cross_linker_mass));
      }


      if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::LOOP)
      {
        ph_alpha.setMetaValue("xl_pos2", DataValue(beta_pos));
      }



      String alpha_term = "ANYWHERE";
      if (alpha_term_spec == ResidueModification::N_TERM)
      {
        alpha_term = "N_TERM";
      }
      else if (alpha_term_spec == ResidueModification::C_TERM)
      {
        alpha_term = "C_TERM";
      }

      ResidueModification::TermSpecificity beta_term_spec = top_csms_spectrum[i].cross_link.term_spec_beta;
      String beta_term = "ANYWHERE";
      if (beta_term_spec == ResidueModification::N_TERM)
      {
        beta_term = "N_TERM";
      }
      else if (beta_term_spec == ResidueModification::C_TERM)
      {
        beta_term = "C_TERM";
      }

      vector<PeptideHit> phs;

      ph_alpha.setSequence(seq_alpha);
      ph_alpha.setCharge(precursor_charge);
      ph_alpha.setScore(top_csms_spectrum[i].score);
      ph_alpha.setRank(DataValue(i+1));
      ph_alpha.setMetaValue("xl_chain", "MS:1002509");  // donor (longer, heavier, alphabetically earlier)
      ph_alpha.setMetaValue("xl_pos", DataValue(alpha_pos));
      ph_alpha.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
      ph_alpha.setMetaValue("xl_type", xltype);
      ph_alpha.setMetaValue("xl_rank", DataValue(i + 1));
      ph_alpha.setMetaValue("xl_term_spec", alpha_term);

      if (scan_index_heavy != scan_index)
      {
        ph_alpha.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
        ph_alpha.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
        ph_alpha.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
      }

      ph_alpha.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
      ph_alpha.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
      ph_alpha.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
      ph_alpha.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
      ph_alpha.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

      ph_alpha.setMetaValue("OpenProXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
      ph_alpha.setMetaValue("OpenProXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
      ph_alpha.setMetaValue("OpenProXL:HyperAlpha", top_csms_spectrum[i].HyperAlpha);
      ph_alpha.setMetaValue("OpenProXL:HyperBeta", top_csms_spectrum[i].HyperBeta);
      ph_alpha.setMetaValue("OpenProXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

      ph_alpha.setMetaValue("OpenProXL:PScoreCommon",top_csms_spectrum[i].PScoreCommon);
      ph_alpha.setMetaValue("OpenProXL:PScoreXlink",top_csms_spectrum[i].PScoreXlink);
      ph_alpha.setMetaValue("OpenProXL:PScoreAlpha",top_csms_spectrum[i].PScoreAlpha);
      ph_alpha.setMetaValue("OpenProXL:PScoreBeta",top_csms_spectrum[i].PScoreBeta);
      ph_alpha.setMetaValue("OpenProXL:PScoreBoth",top_csms_spectrum[i].PScoreBoth);

      ph_alpha.setMetaValue("selected", "false");

      ph_alpha.setFragmentAnnotations(top_csms_spectrum[i].frag_annotations);
      LOG_DEBUG << "Annotations of size " << ph_alpha.getFragmentAnnotations().size() << endl;
      phs.push_back(ph_alpha);

      if (top_csms_spectrum[i].cross_link.getType() == ProteinProteinCrossLink::CROSS)
      {
        ph_beta.setSequence(top_csms_spectrum[i].cross_link.beta);
        ph_beta.setCharge(precursor_charge);
        ph_beta.setScore(top_csms_spectrum[i].score);
        ph_beta.setRank(DataValue(i+1));
        ph_beta.setMetaValue("xl_chain", "MS:1002510"); // receiver
        ph_beta.setMetaValue("xl_pos", DataValue(beta_pos));
        ph_beta.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
        ph_beta.setMetaValue("xl_term_spec", beta_term);

        if (scan_index_heavy != scan_index)
        {
          ph_beta.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
          ph_beta.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
          ph_beta.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
        }

        ph_beta.setMetaValue("OpenXQuest:xcorr xlink", top_csms_spectrum[i].xcorrx_max);
        ph_beta.setMetaValue("OpenXQuest:xcorr common", top_csms_spectrum[i].xcorrc_max);
        ph_beta.setMetaValue("OpenXQuest:match-odds", top_csms_spectrum[i].match_odds);
        ph_beta.setMetaValue("OpenXQuest:intsum", top_csms_spectrum[i].int_sum);
        ph_beta.setMetaValue("OpenXQuest:wTIC", top_csms_spectrum[i].wTIC);

        ph_beta.setMetaValue("OpenProXL:HyperCommon",top_csms_spectrum[i].HyperCommon);
        ph_beta.setMetaValue("OpenProXL:HyperXlink",top_csms_spectrum[i].HyperXlink);
        ph_beta.setMetaValue("OpenProXL:HyperAlpha",top_csms_spectrum[i].HyperAlpha);
        ph_beta.setMetaValue("OpenProXL:HyperBeta",top_csms_spectrum[i].HyperBeta);
        ph_beta.setMetaValue("OpenProXL:HyperBoth",top_csms_spectrum[i].HyperBoth);

        ph_beta.setMetaValue("OpenProXL:PScoreCommon",top_csms_spectrum[i].PScoreCommon);
        ph_beta.setMetaValue("OpenProXL:PScoreXlink",top_csms_spectrum[i].PScoreXlink);
        ph_beta.setMetaValue("OpenProXL:PScoreAlpha",top_csms_spectrum[i].PScoreAlpha);
        ph_beta.setMetaValue("OpenProXL:PScoreBeta",top_csms_spectrum[i].PScoreBeta);
        ph_beta.setMetaValue("OpenProXL:PScoreBoth",top_csms_spectrum[i].PScoreBoth);

        ph_beta.setMetaValue("selected", "false");

        phs.push_back(ph_beta);
      }

      peptide_id.setRT(spectrum_light.getRT());
      peptide_id.setMZ(precursor_mz);
      String specIDs;
      if (scan_index_heavy != scan_index)
      {
        specIDs = spectra[scan_index].getNativeID() + "," + spectra[scan_index_heavy].getNativeID();
      }
      else
      {
        specIDs = spectra[scan_index].getNativeID();
      }

      peptide_id.setMetaValue("spectrum_reference", specIDs);
//      peptide_id.setMetaValue("spec_heavy_RT", spectra[scan_index_heavy].getRT());
//      peptide_id.setMetaValue("spec_heavy_MZ", spectra[scan_index_heavy].getPrecursors()[0].getMZ());
//      peptide_id.setMetaValue("spectrum_reference", spectra[scan_index].getNativeID());
//      peptide_id.setMetaValue("spectrum_reference_heavy", spectra[scan_index_heavy].getNativeID());
//      peptide_id.setMetaValue("xl_type", xltype); // TODO: needs CV term
//      peptide_id.setMetaValue("xl_rank", DataValue(i + 1));

      peptide_id.setHits(phs);
      peptide_id.setScoreType("OpenXQuest:combined score");

#ifdef _OPENMP
#pragma omp critical (peptides_ids_access)
#endif
      {
        peptide_ids.push_back(peptide_id);
        all_top_csms[all_top_csms_current_index][i].peptide_id_index = peptide_ids.size()-1;
      }
    }
  }

  void OpenProXLUtils::writeXQuestXMLSpec(String out_file, String base_name, const PreprocessedPairSpectra& preprocessed_pair_spectra, const vector< pair<Size, Size> >& spectrum_pairs, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
  {
    //String spec_xml_filename = base_name + "_matched.spec.xml";
    // XML Header
    ofstream spec_xml_file;
    cout << "Writing spec.xml to " << out_file << endl;
    spec_xml_file.open(out_file.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
    // TODO write actual data
    spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra compare_peaks_version=\"3.4\" date=\"Tue Nov 24 12:41:18 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"aleitner_M1012_004_matched\" deffile=\"xquest.def\" >" << endl;

    for (Size i = 0; i < spectrum_pairs.size(); ++i)
    {
      if (!all_top_csms[i].empty())
      {
        Size scan_index_light = spectrum_pairs[i].first;
        Size scan_index_heavy = spectrum_pairs[i].second;
        // TODO more correct alternative
        String spectrum_light_name = base_name + ".light." + scan_index_light;
        String spectrum_heavy_name = base_name + ".heavy." + scan_index_heavy;
//        String spectrum_light_name = String("spectrumlight") + scan_index_light;
//        String spectrum_heavy_name = String("spectrumheavy") + scan_index_heavy;
        String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

        // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
        spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[scan_index_light], String(""));
        spec_xml_file << "</spectrum>" << endl;

        spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[scan_index_heavy], String(""));
        spec_xml_file << "</spectrum>" << endl;

        String spectrum_common_name = spectrum_name + String("_common.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(preprocessed_pair_spectra.spectra_common_peaks[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << endl;

        String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(preprocessed_pair_spectra.spectra_xlink_peaks[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << endl;
      }
    }

    spec_xml_file << "</xquest_spectra>" << endl;
    spec_xml_file.close();

    return;
  }

  void OpenProXLUtils::writeXQuestXMLSpec(String out_file, String base_name, const vector< vector< CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra)
  {
    // String spec_xml_filename = base_name + "_matched.spec.xml";
    // XML Header
    ofstream spec_xml_file;
    cout << "Writing spec.xml to " << out_file << endl;
    spec_xml_file.open(out_file.c_str(), ios::trunc); // ios::app = append to file, ios::trunc = overwrites file
    // TODO write actual data
    spec_xml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><xquest_spectra compare_peaks_version=\"3.4\" date=\"Tue Nov 24 12:41:18 2015\" author=\"Thomas Walzthoeni,Oliver Rinner\" homepage=\"http://proteomics.ethz.ch\" resultdir=\"aleitner_M1012_004_matched\" deffile=\"xquest.def\" >" << endl;

    for (Size i = 0; i < spectra.size(); ++i)
    {
      if (!all_top_csms[i].empty())
      {
        String spectrum_light_name = base_name + ".light." + i;
        String spectrum_heavy_name = base_name + ".heavy." + i;

        String spectrum_name = spectrum_light_name + String("_") + spectrum_heavy_name;

        // 4 Spectra resulting from a light/heavy spectra pair.  Write for each spectrum, that is written to xquest.xml (should be all considered pairs, or better only those with at least one sensible Hit, meaning a score was computed)
        spec_xml_file << "<spectrum filename=\"" << spectrum_light_name << ".dta" << "\" type=\"light\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[i], String(""));
        spec_xml_file << "</spectrum>" << endl;

        spec_xml_file << "<spectrum filename=\"" << spectrum_heavy_name << ".dta" << "\" type=\"heavy\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[i], String(""));
        spec_xml_file << "</spectrum>" << endl;

        String spectrum_common_name = spectrum_name + String("_common.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_common_name << "\" type=\"common\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << endl;

        String spectrum_xlink_name = spectrum_name + String("_xlinker.txt");
        spec_xml_file << "<spectrum filename=\"" << spectrum_xlink_name << "\" type=\"xlinker\">" << endl;
        spec_xml_file << getxQuestBase64EncodedSpectrum(spectra[i], spectrum_light_name + ".dta," + spectrum_heavy_name + ".dta");
        spec_xml_file << "</spectrum>" << endl;
      }
    }

    spec_xml_file << "</xquest_spectra>" << endl;
    spec_xml_file.close();

    return;
  }

  double OpenProXLUtils::matched_current_chain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
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

  double OpenProXLUtils::total_matched_current(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks)
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

    // make the indices in the vectors unique
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

  std::vector< double > OpenProXLUtils::xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance)
  {
    // generate vector of results, filled with zeroes
    std::vector< double > results(maxshift * 2 + 1, 0);

    // return 0 = no correlation, either positive nor negative, when one of the spectra is empty (e.g. when no common ions or xlink ions could be matched between light and heavy spectra)
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
      // TODO this line leads to using real intensities
//      ion_table1[pos] = spec1[i].getIntensity();
      // TODO this line leads to using intensities normalized to 10
      ion_table1[pos] = 10.0;
    }
    for (Size i = 0; i < spec2.size(); ++i)
    {
      Size pos =static_cast<Size>(ceil(spec2[i].getMZ() / tolerance));
      // TODO this line leads to using real intensities
//      ion_table2[pos] = spec2[i].getIntensity();
      // TODO this line leads to using intensities normalized to 10
      ion_table2[pos] = 10.0;
    }

    // Compute means for real intensities
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
