// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Petra Gutenbrunner $
// $Authors: David Wojnar, Timo Sachsenberg, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/ID/PScore.h>

#include <limits>
#include <vector>

namespace OpenMS
{
  class PeptideHit;
  class AASequence;
  
  class ProbablePhosphoSites
  {
  public:

    Size first;
    Size second;
    Size seq_1; ///< index of best permutation with site in phosphorylated state
    Size seq_2; ///< index of permutation with site in unphosphorylated state
    Size peak_depth; ///< filtering level that gave rise to maximum discriminatory score
    Size AScore;
  };
  
  /**
      @brief Implementation of the Ascore
      For a given peptide sequence and its MS/MS spectrum it identifies the most probable phosphorylation-site(s).
      For each phosphorylation site a probability score is calculated.
      The algorithm is implemented according to Beausoleil et al. (Nat. Biotechnol. 2006).

      @htmlinclude OpenMS_AScore.parameters
  */
  class OPENMS_DLLAPI AScore: public DefaultParamHandler
  {
    friend struct PScore;
    
  public:
    ///Default constructor
    AScore();

    ///Destructor
    ~AScore() override;

    /**
        @brief Computes the AScore and returns all computed phospho-sites. The saved sequences contain only phospho information. All other modifications are dropped due to simplicity.

        @param	hit a PeptideHit
        @param real_spectrum spectrum mapped to hit

        @note the original sequence is saved in the PeptideHits as MetaValue Search_engine_sequence.
    */
    PeptideHit compute(const PeptideHit& hit, PeakSpectrum& real_spectrum);

  protected:
    int compareMZ_(double mz1, double mz2) const;
    
    /// getSpectrumDifference_ works similar as the method std::set_difference (http://en.cppreference.com/w/cpp/algorithm/set_difference). 
    /// set_difference was reimplemented, because it was necessary to overwrite the compare operator to be able to compare the m/z values.
    /// not implemented as "operator<", because using tolerances for comparison does not imply total ordering    
    template <class InputIterator1, class InputIterator2, class OutputIterator>
    OutputIterator getSpectrumDifference_(InputIterator1 first1, InputIterator1 last1,
      InputIterator2 first2, InputIterator2 last2, OutputIterator result) const
    {
      while (first1 != last1 && first2 != last2)
      {
        double mz1 = first1->getMZ();
        double mz2 = first2->getMZ();
        int val = compareMZ_(mz1, mz2);
        
        if (val == -1)
        { 
          *result = *first1; 
          ++result; 
          ++first1; 
        }
        else if (val == 1)
        {
          ++first2;
        }
        else // check if more ions are within the same tolerance. If so, these can not be site determining ions
        {
          //check mz2 until no match
          ++first2;
          if (first2 != last2)
          {
            int ret = compareMZ_(mz1, first2->getMZ());
            while (ret == 0 && first2 != last2)
            {
              ++first2;
              ret = compareMZ_(mz1, first2->getMZ());
            }
          }
          
          //check mz1 until no match
          ++first1;
          if (first1 != last1)
          {
            int ret = compareMZ_(first1->getMZ(), mz2);
            while (ret == 0 && first1 != last1)
            {
              ++first1;
              ret = compareMZ_(first1->getMZ(), mz2);
            }
          }
        }
      }
      return std::copy(first1, last1, result);
    }
    
    ///Computes the site determining_ions for the given AS and sequences in candidates
    void computeSiteDeterminingIons_(const std::vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, std::vector<PeakSpectrum>& site_determining_ions) const;

    /// return all phospho sites
    std::vector<Size> getSites_(const String& unmodified_sequence) const;

    /// calculate all n_phosphorylation_events sized sets of phospho sites (all versions of the peptides with exactly n_phosphorylation_events)
    std::vector<std::vector<Size>> computePermutations_(const std::vector<Size>& sites, Int n_phosphorylation_events) const;

    /// Computes number of matched ions between windows and the given spectrum. All spectra have to be sorted by position!
    Size numberOfMatchedIons_(const PeakSpectrum& th, const PeakSpectrum& windows, Size depth) const;

    /// Computes the peptide score according to Beausoleil et al. page 1291
    double peptideScore_(const std::vector<double>& scores) const;

    /**
        @brief Finds the peptides with the highest PeptideScores and outputs all information for computing the AScore
        @note This function assumes that there are more permutations than the assumed number of phosphorylations!
    */
    void determineHighestScoringPermutations_(const std::vector<std::vector<double>>& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const std::vector<std::vector<Size>>& permutations, std::multimap<double, Size>& ranking) const;

    /// Computes probability for a peak depth of one given spectra and mass_tolerance variables
    double computeBaseProbability_(double ppm_reference_mz) const;

    /// Computes the cumulative binomial probabilities.
    double computeCumulativeScore_(Size N, Size n, double p) const;
    
    /// Computes number of phospho events in a sequence
    Size numberOfPhosphoEvents_(const String& sequence) const;
    
    /// Create variant of the peptide with all phosphorylations removed
    AASequence removePhosphositesFromSequence_(const String& sequence) const;
    
    /// Create theoretical spectra with all combinations with the number of phosphorylation events
    std::vector<PeakSpectrum> createTheoreticalSpectra_(const std::vector<std::vector<Size>>& permutations, const AASequence& seq_without_phospho) const;
    
    /// Pick top 10 intensity peaks for each 100 Da windows
    std::vector<PeakSpectrum> peakPickingPerWindowsInSpectrum_(PeakSpectrum& real_spectrum) const;
    
    /// Create 10 scores for each theoretical spectrum (permutation), according to Beausoleil et al. Figure 3 b
    std::vector<std::vector<double>> calculatePermutationPeptideScores_(std::vector<PeakSpectrum>& th_spectra, const std::vector<PeakSpectrum>& windows_top10) const;
    
    /// Rank weighted permutation scores ascending
    std::multimap<double, Size> rankWeightedPermutationPeptideScores_(const std::vector<std::vector<double>>& peptide_site_scores) const;

    /// Reimplemented from @ref DefaultParamHandler
    void updateMembers_() override;

    // variables:
    double fragment_mass_tolerance_; ///< Fragment mass tolerance for spectrum comparisons
    bool fragment_tolerance_ppm_; ///< Is fragment mass tolerance given in ppm (or Da)?
    Size max_peptide_length_; ///< Limit for peptide lengths that can be analyzed
    Size max_permutations_; ///< Limit for number of sequence permutations that can be handled
    double unambiguous_score_; ///< Score for unambiguous assignments (all sites phosphorylated)
    double base_match_probability_; ///< Probability of a match at a peak depth of 1
  };

} // namespace OpenMS
