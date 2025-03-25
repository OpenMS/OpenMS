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
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

#include <limits>
#include <vector>

namespace OpenMS
{
  class PeptideHit;
  class AASequence;
  
  struct ProbablePhosphoSites
  {
    /// Index of the site evaluated in the best scoring permutation
    Size first{};
    
    /// Index of the site evaluated in the second best scoring permutation
    Size second{};
    
    /// Index of best permutation with site in phosphorylated state
    Size seq_1{};
    
    /// Index of permutation with site in unphosphorylated state
    Size seq_2{};
    
    /// Filtering level (peak depth) that gave rise to maximum discriminatory score
    Size peak_depth{};
    double ascore{}; ///< The score difference between phosphorylated and unphosphorylated state for that site
  };
  
  /**
      @brief Implementation of the Ascore
      For a given peptide sequence and its MS/MS spectrum it identifies the most probable phosphorylation-site(s).
      For each phosphorylation site a probability score is calculated.
      
      The algorithm is implemented according to Beausoleil et al. (Nat. Biotechnol. 2006):
      "A probability-based approach for high-throughput protein phosphorylation analysis and site localization"
      
      The AScore is a probabilistic measure of phosphorylation site localization based on the presence
      and intensity of site-determining ions in MS/MS spectra. The algorithm:
      
      1. Generates all possible phosphorylation site permutations for a given peptide
      2. Creates theoretical spectra for each permutation
      3. Compares these theoretical spectra with the experimental spectrum
      4. Calculates a score for each permutation based on matched peaks
      5. Determines the most likely phosphorylation site(s) based on score differences

      In addition to standard phosphorylation on S, T, and Y residues, this implementation also supports
      PhosphoDecoy modifications on A, G, and L residues which can be used for false localization rate (FLR) calculations.
      
      The AScore value represents the probability of correct phosphorylation site localization. Higher scores
      indicate greater confidence in site localization. Typically, scores above 19 correspond to >99% certainty
      in site localization.

      @htmlinclude OpenMS_AScore.parameters
  */
  class OPENMS_DLLAPI AScore: public DefaultParamHandler
  {
    friend struct PScore;
    
  public:
    ///Default constructor
    AScore();

    ///Destructor
    ~AScore() override = default;

    /**
        @brief Computes the AScore and returns all computed phospho-sites
        
        This is the main method for phosphorylation site localization. It takes a peptide hit and its
        corresponding MS/MS spectrum and determines the most probable phosphorylation site(s).
        
        The method:
        1. Extracts phosphorylation sites from the peptide sequence
        2. Generates all possible phosphorylation site permutations
        3. Creates theoretical spectra for each permutation
        4. Compares these with the experimental spectrum using a windowed peak picking approach
        5. Calculates scores for each permutation and determines the most likely sites
        6. Returns a modified PeptideHit with AScore values for each phosphorylation site
        
        The saved sequences contain only phospho information. All other modifications are dropped for simplicity.
        The original sequence is saved in the PeptideHits as MetaValue Search_engine_sequence.

        @param hit A PeptideHit containing the peptide sequence with potential phosphorylation sites
        @param real_spectrum The experimental MS/MS spectrum mapped to the hit
        @return A modified PeptideHit with AScore values and phosphorylation site localization information
                The score field contains the best AScore value (lower is better)
    */
    PeptideHit compute(const PeptideHit& hit, PeakSpectrum& real_spectrum);

    /// return is char is a phospho decoy site
    static bool isPhosphoDecoySite(char residue);

    /// return is char is a phospho site
    static bool isPhosphoSite(char residue);

  protected:
    /**
     * @brief Compares two m/z values using the fragment mass tolerance
     * 
     * @param mz1 First m/z value
     * @param mz2 Second m/z value
     * @return int -1 if mz1 < mz2, 1 if mz1 > mz2, 0 if they are within tolerance
     */
    int compareMZ_(double mz1, double mz2) const;
    
    /**
     * @brief Finds the difference between two spectra based on m/z values
     * 
     * This method works similar to std::set_difference but is specialized for comparing
     * mass spectra with a tolerance window. It was reimplemented because it was necessary
     * to overwrite the compare operator to be able to compare the m/z values with tolerance.
     * It is not implemented as "operator<" because using tolerances for comparison does not
     * imply total ordering.
     */
    template <class InputIterator1, class InputIterator2, class OutputIterator>
    OutputIterator getSpectrumDifference_(InputIterator1 first1, InputIterator1 last1,
      InputIterator2 first2, InputIterator2 last2, OutputIterator result) const
    {
      // Compare spectra and find peaks that are unique to the first spectrum
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
    
    /**
     * @brief Computes the site-determining ions for the given phosphorylation site candidates
     * 
     * Site-determining ions are fragment ions that can distinguish between different phosphorylation
     * site localizations. This method identifies ions that are unique to each of the two best-scoring
     * permutations for a given phosphorylation site.
     * 
     * @param th_spectra Vector of theoretical spectra for all permutations
     * @param candidates The phosphorylation site candidates to evaluate
     * @param site_determining_ions Output vector to store the site-determining ions
     */
    void computeSiteDeterminingIons_(const std::vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, std::vector<PeakSpectrum>& site_determining_ions) const;

    /**
     * @brief Identifies all potential phosphorylation sites in a peptide sequence
     * 
     * Returns the positions of all residues that can be phosphorylated (S, T, Y)
     * and optionally PhosphoDecoy sites (A, G, L) if enabled.
     * 
     * @param unmodified_sequence The unmodified peptide sequence
     * @return Vector of positions (0-based) of potential phosphorylation sites
     */
    std::vector<Size> getSites_(const String& unmodified_sequence) const;

    /**
     * @brief Generates all possible combinations of phosphorylation sites
     * 
     * Calculates all possible ways to place n_phosphorylation_events phosphorylations
     * on the available phosphorylation sites. Uses an efficient combinatorial algorithm
     * with early termination for large numbers of permutations.
     * 
     * @param sites Vector of positions of potential phosphorylation sites
     * @param n_phosphorylation_events Number of phosphorylation events to place
     * @return Vector of all possible phosphorylation site combinations
     */
    std::vector<std::vector<Size>> computePermutations_(const std::vector<Size>& sites, Int n_phosphorylation_events) const;

    /**
     * @brief Counts the number of matched ions between a theoretical spectrum and experimental spectrum window
     * 
     * Compares a theoretical spectrum with an experimental spectrum window and counts
     * the number of matching peaks within the fragment mass tolerance.
     * 
     * @param th Theoretical spectrum (must be sorted by m/z)
     * @param windows Experimental spectrum window (must be sorted by m/z)
     * @param depth Maximum number of peaks to consider from the window (peak depth)
     * @return Number of matched ions
     */
    Size numberOfMatchedIons_(const PeakSpectrum& th, const PeakSpectrum& windows, Size depth) const;

    /**
     * @brief Computes the weighted peptide score according to Beausoleil et al.
     * 
     * Calculates a weighted average of scores at different peak depths according to
     * the formula described in Beausoleil et al. (page 1291). The weights are:
     * - 0.5 for depth 1
     * - 0.75 for depth 2
     * - 1.0 for depths 3-6
     * - 0.75 for depth 7
     * - 0.5 for depth 8
     * - 0.25 for depths 9-10
     * The sum is divided by 7.0 to normalize.
     * 
     * @param scores Vector of 10 scores at different peak depths
     * @return Weighted peptide score
     */
    double peptideScore_(const std::vector<double>& scores) const;

    /**
        @brief Finds the peptides with the highest PeptideScores and determines site-specific AScores
        
        For each phosphorylation site in the highest-scoring permutation, this function:
        1. Finds the next best permutation where this site is not phosphorylated
        2. Determines the peak depth that maximizes the score difference between these permutations
        3. Stores this information for AScore calculation
        
        @param peptide_site_scores Scores for each permutation at each peak depth
        @param sites Output vector to store phosphorylation site information
        @param permutations Vector of all phosphorylation site permutations
        @param ranking Ranked permutations by their weighted peptide scores
        
        @note This function assumes that there are more permutations than the number of phosphorylations!
    */
    void determineHighestScoringPermutations_(const std::vector<std::vector<double>>& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const std::vector<std::vector<Size>>& permutations, std::multimap<double, Size>& ranking) const;

    /**
     * @brief Computes the base match probability for a given reference m/z
     * 
     * Calculates the probability of a random peak match based on the fragment mass tolerance.
     * For Da tolerance, this is 2 * tolerance / 100.
     * For ppm tolerance, this is 2 * tolerance * reference_mz * 1e-6 / 100.
     * 
     * @param ppm_reference_mz Reference m/z value for ppm calculations
     * @return Base probability of a random peak match
     */
    double computeBaseProbability_(double ppm_reference_mz) const;

    /**
     * @brief Computes the cumulative binomial probability P(X ≥ n)
     * 
     * Calculates the probability of observing at least n successes in N trials
     * with individual success probability p. This is used to determine the
     * statistical significance of peak matches.
     * 
     * @param N Total number of trials (theoretical peaks)
     * @param n Number of successes (matched peaks)
     * @param p Probability of success for a single trial
     * @return Cumulative binomial probability P(X ≥ n)
     */
    double computeCumulativeScore_(Size N, Size n, double p) const;
    
    /**
     * @brief Counts the number of phosphorylation events in a peptide sequence
     * 
     * Counts both regular phosphorylation events (Phospho) and decoy phosphorylation
     * events (PhosphoDecoy) if enabled.
     * 
     * @param sequence The peptide sequence as a string
     * @return Number of phosphorylation events
     */
    Size numberOfPhosphoEvents_(const String& sequence) const;
    
    /**
     * @brief Creates a variant of the peptide with all phosphorylations removed
     * 
     * Removes both regular Phospho and PhosphoDecoy modifications from the sequence.
     * 
     * @param sequence The peptide sequence as a string
     * @return AASequence object without phosphorylations
     */
    AASequence removePhosphositesFromSequence_(const String& sequence) const;
    
    /**
     * @brief Creates theoretical spectra for all phosphorylation site permutations
     * 
     * For each permutation of phosphorylation sites, this method:
     * 1. Creates a peptide sequence with phosphorylations at the specified positions
     * 2. Generates a theoretical spectrum with b- and y-ions
     * 3. Sets the spectrum name to the peptide sequence
     * 
     * @param permutations Vector of all phosphorylation site permutations
     * @param seq_without_phospho Base peptide sequence without phosphorylations
     * @return Vector of theoretical spectra for all permutations
     */
    std::vector<PeakSpectrum> createTheoreticalSpectra_(const std::vector<std::vector<Size>>& permutations, const AASequence& seq_without_phospho) const;
    
    /**
     * @brief Picks the top 10 intensity peaks for each 100 Da window in the spectrum
     * 
     * Divides the spectrum into 100 Da windows and selects the 10 most intense peaks
     * from each window. This approach normalizes for intensity variations across
     * the m/z range and focuses on the most informative peaks.
     * 
     * @param real_spectrum The experimental MS/MS spectrum
     * @return Vector of spectra containing the top 10 peaks for each window
     */
    std::vector<PeakSpectrum> peakPickingPerWindowsInSpectrum_(PeakSpectrum& real_spectrum) const;
    
    /**
     * @brief Calculates scores for each permutation at 10 different peak depths
     * 
     * For each theoretical spectrum (permutation), this method:
     * 1. Counts matched ions at each peak depth (1-10)
     * 2. Calculates the cumulative binomial probability
     * 3. Converts to a score using -10*log10(probability)
     * 
     * This corresponds to Figure 3b in Beausoleil et al.
     * 
     * @param th_spectra Vector of theoretical spectra for all permutations
     * @param windows_top10 Vector of experimental spectra windows with top 10 peaks
     * @return Vector of scores for each permutation at each peak depth
     */
    std::vector<std::vector<double>> calculatePermutationPeptideScores_(std::vector<PeakSpectrum>& th_spectra, const std::vector<PeakSpectrum>& windows_top10) const;
    
    /**
     * @brief Ranks permutations by their weighted peptide scores
     * 
     * Calculates the weighted peptide score for each permutation and
     * creates a multimap ranking them in ascending order.
     * 
     * @param peptide_site_scores Scores for each permutation at each peak depth
     * @return Multimap of weighted scores to permutation indices
     */
    std::multimap<double, Size> rankWeightedPermutationPeptideScores_(const std::vector<std::vector<double>>& peptide_site_scores) const;

    /**
     * @brief Generates a ProForma-like string with phosphorylation site localization scores
     * 
     * Creates a string representation of the peptide with phosphorylation sites
     * annotated with their localization scores. The format is:
     * X[Phospho|score=0.9999] where X is the amino acid and 0.9999 is the localization probability.
     * 
     * @param peptide The peptide sequence
     * @param ascores Map of positions to AScore values
     * @return ProForma-like string with phosphorylation site scores
     */
    String generateProFormaString_(const AASequence& peptide, const std::map<Size, double>& ascores) const;

    /// Reimplemented from @ref DefaultParamHandler
    void updateMembers_() override;

    // variables:
    /// Fragment mass tolerance for spectrum comparisons
    double fragment_mass_tolerance_;
    
    /// Is fragment mass tolerance given in ppm (true) or Da (false)?
    bool fragment_tolerance_ppm_;
    
    /// Maximum peptide length that can be analyzed (0 for no limit)
    Size max_peptide_length_;
    
    /// Maximum number of sequence permutations that can be handled (0 for no limit)
    Size max_permutations_;
    
    /// Score for unambiguous assignments (all sites phosphorylated)
    const double unambiguous_score_ = 1e6;
    
    /// Probability of a random peak match at a peak depth of 1
    double base_match_probability_;
    
    /// Include PhosphoDecoy sites (A, G, L) in phosphorylation site analysis
    bool add_decoys_;
    
    /// Generator for theoretical spectra
    TheoreticalSpectrumGenerator spectrum_generator_;
  };

} // namespace OpenMS
