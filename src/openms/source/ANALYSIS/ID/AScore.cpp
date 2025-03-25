// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Petra Gutenbrunner $
// $Authors: David Wojnar, Timo Sachsenberg, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AScore.h>

#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <sstream>
#include <iomanip>

using namespace std;

namespace OpenMS
{

  /**
   * @brief Constructor for AScore
   * 
   * Initializes the AScore object with default parameters:
   * - fragment_mass_tolerance: 0.05 Da
   * - Theoretical spectrum generator configured for phosphorylation analysis
   */
  AScore::AScore():
    DefaultParamHandler("AScore") // Initialize base class with algorithm name
  {
    defaults_.setValue("fragment_mass_tolerance", 0.05, "Fragment mass tolerance for spectrum comparisons");
    defaults_.setMinFloat("fragment_mass_tolerance", 0.0);

    defaults_.setValue("fragment_mass_unit", "Da", "Unit of fragment mass tolerance");
    defaults_.setValidStrings("fragment_mass_unit", {"Da","ppm"});

    vector<std::string> advanced(1, "advanced"); // tag for advanced parameters

    defaults_.setValue("max_peptide_length", 40, "Restrict scoring to peptides with a length no greater than this value ('0' for 'no restriction')", advanced);
    defaults_.setMinInt("max_peptide_length", 0);

    defaults_.setValue("max_num_perm", 16384, "Maximum number of permutations a sequence can have to be processed ('0' for 'no restriction')", advanced);
    defaults_.setMinInt("max_num_perm", 0);

    defaults_.setValue("unambiguous_score", 1000, "Score to use for unambiguous assignments, where all sites on a peptide are phosphorylated. (Note: If a peptide is not phosphorylated at all, its score is set to '-1'.)", advanced);

    defaults_.setValue("add_decoys", "false", "Include PhosphoDecoy sites (A, G, L) in phosphorylation site analysis for FLR calculation", advanced);
    defaults_.setValidStrings("add_decoys", {"true", "false"});

    // Apply default parameters to the current parameter set
    defaultsToParam_();

    Param p = spectrum_generator_.getDefaults();
    p.setValue("isotope_model", "none");
    p.setValue("add_first_prefix_ion", "true");
    spectrum_generator_.setParameters(p);
  }

  /**
   * @brief Main method to compute phosphorylation site localization scores
   * 
   * This method implements the AScore algorithm for phosphorylation site localization:
   * 1. Extracts phosphorylation sites from the peptide sequence
   * 2. Generates all possible phosphorylation site permutations
   * 3. Creates theoretical spectra for each permutation
   * 4. Compares these with the experimental spectrum using a windowed peak picking approach
   * 5. Calculates scores for each permutation and determines the most likely sites
   * 
   * @param hit The peptide hit containing the sequence with potential phosphorylation sites
   * @param real_spectrum The experimental MS/MS spectrum
   * @return A modified PeptideHit with AScore values for each phosphorylation site
   *         The score field contains the best AScore value
   */
  PeptideHit AScore::compute(const PeptideHit& hit, PeakSpectrum& real_spectrum)
  {
    PeptideHit phospho = hit;
    
    // reset phospho
    phospho.setScore(-1);
    phospho.setMetaValue("search_engine_sequence", hit.getSequence().toString());

    // Early termination for empty spectra
    if (real_spectrum.empty()) 
    {
      OPENMS_LOG_WARN << "Empty spectrum provided to AScore::compute. Returning original hit." << std::endl;
      return phospho;
    }
        
    String sequence_str = phospho.getSequence().toString();
    String unmodified_sequence_str = phospho.getSequence().toUnmodifiedString();

    // Count phosphorylation events in the sequence (both normal and decoy phosphorylations)
    // and remove phosphorylations to get the base peptide sequence
    Size number_of_phosphorylation_events = numberOfPhosphoEvents_(sequence_str); 
    AASequence seq_without_phospho = removePhosphositesFromSequence_(sequence_str); // remove Phospho and PhosphoDecoy from sequence

    // Initialize counters for regular and decoy phosphorylation events
    // These will be used for metadata in the output
    // We need to count them separately from the total number of phosphorylation events
    // to provide detailed information about the types of modifications
    Size regular_phospho_count = 0;
    Size decoy_phospho_count = 0;
    
    // Count regular phosphorylation events
    Size pos = 0;
    while ((pos = sequence_str.find("(Phospho)", pos)) != std::string::npos)
    {
      ++regular_phospho_count;
      pos += 9; // Move past "(Phospho)"
    }
    
    // Count PhosphoDecoy events
    pos = 0;
    while ((pos = sequence_str.find("(PhosphoDecoy)", pos)) != std::string::npos)
    {
      ++decoy_phospho_count;
      pos += 14; // Move past "(PhosphoDecoy)"
    }
    
    if (add_decoys_)
    {
      OPENMS_LOG_DEBUG << "Found " << regular_phospho_count << " regular phosphorylation events and " 
                      << decoy_phospho_count << " PhosphoDecoy events in sequence: " << sequence_str << std::endl;
    }
    else if (decoy_phospho_count > 0)
    {
      OPENMS_LOG_WARN << "PhosphoDecoy sites found in sequence, but add_decoys is set to false. "
                      << "Please enable add_decoys to include PhosphoDecoy sites in the analysis." 
                      << "Returning original hit." << std::endl;
      return phospho;
    }

    // Check if peptide is too long for analysis (if max_peptide_length_ is set)
    if ((max_peptide_length_ > 0) && (unmodified_sequence_str.size() > max_peptide_length_))
    {
      OPENMS_LOG_DEBUG << "\tcalculation aborted: peptide too long: " << seq_without_phospho.toString() << std::endl;
      return phospho;
    }

    // Get all potential phosphorylation sites (S,T,Y and optionally A,G,L for decoys)
    // determine all phospho sites (including decoy sites if add_decoys is true)
    vector<Size> sites = getSites_(unmodified_sequence_str);
    Size number_of_sites = sites.size();

    if (number_of_phosphorylation_events == 0 || number_of_sites == 0)
    {
      return phospho;
    }


    // Add metadata about phosphorylation types
    if (add_decoys_)
    {
      phospho.setMetaValue("regular_phospho_count", regular_phospho_count);
      phospho.setMetaValue("phospho_decoy_count", decoy_phospho_count);
    }
    phospho.setMetaValue("phospho_sites", number_of_phosphorylation_events);

    // Special case: If all possible sites are phosphorylated, the localization is unambiguous
    // Set a high score to indicate certainty
    if (number_of_sites == number_of_phosphorylation_events) 
    {
      phospho.setScore(unambiguous_score_);
      
      // Create a map with scores of 1.0 for all phosphorylation sites
      std::map<Size, double> site_scores;
      for (Size site : sites)  // 'sites' contains all phosphorylation site positions
      {
        site_scores[site] = unambiguous_score_;  // This will be converted to probability ~1.0
      }
      
      // Create ProForma string with scores for all sites
      String proforma = generateProFormaString_(phospho.getSequence(), site_scores);
      phospho.setMetaValue("ProForma", proforma);
      phospho.setMetaValue("AScore_pep_score", unambiguous_score_); // initialize score with unambiguous score  
      return phospho;
    } 

    // Generate all possible permutations of phosphorylation sites
    vector<vector<Size>> permutations = computePermutations_(sites, (Int)number_of_phosphorylation_events);
    OPENMS_LOG_DEBUG << "\tnumber of permutations: " << permutations.size() << std::endl;

    // Check if permutations is empty (early termination) or exceeds the maximum allowed
    // Early termination can happen if the number of permutations is too large to compute
    // or if there are more phosphorylation events than available sites
    // TODO: using a heuristic to calculate the best phospho sites if the number of permutations exceeds the maximum.
    // A heuristic could be to calculate the best site for the first phosphorylation and based on this the best site for the second
    // phosphorylation and so on until every site is determined.
    if (permutations.empty() || ((max_permutations_ > 0) && (permutations.size() > max_permutations_)))
    {
      OPENMS_LOG_DEBUG << "\tcalculation aborted: number of permutations exceeded or early termination" << std::endl;
      return phospho;
    }

    // Create theoretical spectra for all permutations
    vector<PeakSpectrum> th_spectra = createTheoreticalSpectra_(permutations, seq_without_phospho);

    // prepare real spectrum windows
    if (!real_spectrum.isSorted())
    {
      real_spectrum.sortByPosition();
    }
    vector<PeakSpectrum> windows_top10 = peakPickingPerWindowsInSpectrum_(real_spectrum);
    
    // compute match probability for a peak depth of 1
    // This is the probability of a random match between a theoretical and experimental peak
    base_match_probability_ = computeBaseProbability_(real_spectrum.back().getMZ());

    // calculate peptide score for each possible phospho site permutation
    // This generates scores at 10 different peak depths for each permutation
    vector<vector<double>> peptide_site_scores = calculatePermutationPeptideScores_(th_spectra, windows_top10);

    // Sort peptide permutations by (weighted) peptide score. 
    // Lower scores are better (represent lower p-values)
    multimap<double, Size> ranking = rankWeightedPermutationPeptideScores_(peptide_site_scores);

    
    // Get the best scoring permutation and set it as the peptide sequence
    multimap<double, Size>::reverse_iterator rev = ranking.rbegin();

    // Note: there might be multiple permutations with the highest score. If the original
    // sequence is one of them, we take this one. Otherwise, we take the first one.
    double best_score = ranking.rbegin()->first;
    auto best_range = ranking.equal_range(best_score);
    for (auto it = best_range.first; it != best_range.second; ++it)
    {
      auto index = it->second; // Get the index of the permutation with the best score

      // get sequence of this best-scoring permutation
      String best_sequence = th_spectra[index].getName(); // Retrieve the sequence of the best scoring permutation
      // if the sequence is the same as the original sequence choose this one
      if (best_sequence == sequence_str)
      {
        std::swap(rev->second, it->second); // Swap the indices to ensure the original sequence is first
        break;
      }
    }

    // The name of the spectrum contains the peptide sequence
    String seq1 = th_spectra[rev->second].getName();
    phospho.setSequence(AASequence::fromString(seq1));

    double peptide1_score = rev->first;
    phospho.setMetaValue("AScore_pep_score", peptide1_score); // initialize score with highest peptide score (aka highest weighted score)

    ++rev;
    String seq2 = th_spectra[rev->second].getName();
    double peptide2_score = rev->first;

    // Determine the highest scoring permutations for each phosphorylation site
    vector<ProbablePhosphoSites> phospho_sites;
    determineHighestScoringPermutations_(peptide_site_scores, phospho_sites, permutations, ranking);

    // Calculate AScore for each phosphorylation site
    // AScore is the difference in scores between the best and second-best permutations
    Int rank = 1;
    double best_Ascore = std::numeric_limits<double>::max(); // the lower the better
    map<Size, double> site2score; // map to store the scores for each phospho site
    for (auto& phospho_site : phospho_sites)
    {
      double Ascore = 0.0;
      if (peptide1_score == peptide2_score) // set Ascore = 0 for each phosphorylation site
      {
        OPENMS_LOG_DEBUG << "\tscore of best (" << seq1 << ") and second best peptide (" << seq2 << ") are equal (" << peptide1_score << ")" << std::endl;
      }
      else
      {
        // Calculate site-determining ions - these are ions that can distinguish between phosphorylation sites
        vector<PeakSpectrum> site_determining_ions;

        computeSiteDeterminingIons_(th_spectra, phospho_site, site_determining_ions);
        Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
        double p = static_cast<double>(phospho_site.peak_depth) * base_match_probability_;

        Size n_first = 0; // number of matching peaks for first peptide
        // Count matched ions across all windows for the first permutation
        for (Size window_idx = 0; window_idx != windows_top10.size(); ++window_idx) // for each 100 m/z window
        {
          n_first += numberOfMatchedIons_(site_determining_ions[0], windows_top10[window_idx], phospho_site.peak_depth);        
        }
        double P_first = computeCumulativeScore_(N, n_first, p);

        Size n_second = 0; // number of matching peaks for second peptide
        // Count matched ions across all windows for the second permutation
        for (Size window_idx = 0; window_idx < windows_top10.size(); ++window_idx) // each 100 m/z window
        {
          n_second += numberOfMatchedIons_(site_determining_ions[1], windows_top10[window_idx], phospho_site.peak_depth);        
        }
        Size N2 = site_determining_ions[1].size(); // all possibilities have the same number so take the first one
        double P_second = computeCumulativeScore_(N2, n_second, p);

        // Convert probabilities to scores: -10 * log10(P)
        //abs is used to avoid -0 score values
        double score_first = abs(-10. * log10(P_first));
        double score_second = abs(-10. * log10(P_second));

        OPENMS_LOG_DEBUG << "\tfirst - N: " << N << ",p: " << p << ",n: " << n_first << ", score: " << score_first << std::endl;
        OPENMS_LOG_DEBUG << "\tsecond - N: " << N2 << ",p: " << p << ",n: " << n_second << ", score: " << score_second << std::endl;

        // AScore is the difference between the two scores
        Ascore = score_first - score_second;
        OPENMS_LOG_DEBUG << "\tAscore_" << rank << ": " << Ascore << std::endl;
      }

      if (Ascore < best_Ascore)
      {
        best_Ascore = Ascore;
      }
      phospho.setMetaValue("AScore_" + String(rank), Ascore);
      //std::cout << "Rank:" << phospho_site.first << " " << phospho_site.second << " " << phospho_site.peak_depth << " " << Ascore << std::endl;
      site2score[phospho_site.first] = Ascore;
      ++rank;      
    }

    // Generate a ProForma-like string with phosphorylation site scores
    // Add ProForma string as meta value
    String proforma = generateProFormaString_(phospho.getSequence(), site2score);
    phospho.setMetaValue("ProForma", proforma);
    phospho.setScore(best_Ascore);
    return phospho;
  }

  /**
   * @brief Computes the base probability of a random peak match
   * 
   * The base probability is calculated as:
   * - For Da tolerance: 2 * tolerance / 100
   * - For ppm tolerance: 2 * tolerance * reference_mz * 1e-6 / 100
   * 
   * @param ppm_reference_mz Reference m/z value for ppm calculations
   * @return Base probability of a random peak match
   */
  double AScore::computeBaseProbability_(double ppm_reference_mz) const
  {
    double base_match_probability = 2. * fragment_mass_tolerance_ / 100.;
    if (fragment_tolerance_ppm_)
    {
      base_match_probability *= ppm_reference_mz * 1e-6; // 1e-6 converts fragment_mass_tolerance_ to ppm
    }
    return base_match_probability;
  }

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
   */
  double AScore::computeCumulativeScore_(Size N, Size n, double p) const
  {
    OPENMS_PRECONDITION(n <= N, "The number of matched ions (n) can be at most as large as the number of trials (N).");
    OPENMS_PRECONDITION(p >= 0 && p <= 1.0, "p must be a probability [0,1].");

    // Use the numerically stable implementation from MathFunctions
    // This calculates P(X ≥ n) for a binomial distribution with parameters N and p
    // which is exactly what we need for the AScore calculation
    return Math::binomial_cdf_complement(N, n, p);
  }

  /**
   * @brief Determines the highest scoring permutations for each phosphorylation site
   * 
   * For each phosphorylation site in the highest-scoring permutation, this function:
   * 1. Finds the next best permutation where this site is not phosphorylated
   * 2. Determines the peak depth that maximizes the score difference between these permutations
   * 3. Stores this information for AScore calculation
   * 
   * The algorithm is complex because it needs to find permutations that differ by exactly
   * one phosphorylation site while keeping all other sites the same.
   * 
   * @param peptide_site_scores Scores for each permutation at each peak depth
   * @param sites Output vector to store phosphorylation site information
   * @param permutations Vector of all phosphorylation site permutations
   * @param ranking Ranked permutations by their weighted peptide scores
   */
  void AScore::determineHighestScoringPermutations_(const std::vector<std::vector<double>>& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const vector<vector<Size>>& permutations, std::multimap<double, Size>& ranking) const
  {
    // For every phospho site of the highest (weighted) scoring phospho site assignment:
    // 1. determine the next best (weighted) score assignment with this site in unphosporylated state.
    // 2. determine the filtering level (peak depths) that maximizes the (unweighted) score difference between these two assignments

    sites.clear();
    // take first set of phospho site assignments
    sites.resize(permutations[0].size());    
    const vector<Size> & best_peptide_sites = permutations[ranking.rbegin()->second]; // sites of the assignment that achieved the highest weighted score

    for (Size i = 0; i < best_peptide_sites.size(); ++i)  // for each phosphorylated site
    {
      // Start with the highest scoring permutation
      multimap<double, Size>::reverse_iterator rev = ranking.rbegin();      
      sites[i].first = best_peptide_sites[i]; // store the site
      sites[i].seq_1 = rev->second; // and permutation
      bool peptide_not_found = true;
      
      // iterate from best scoring peptide to the first peptide that doesn't contain the current phospho site
      do
      {
        // Move to the next permutation in the ranking
        ++rev;
        // Check if this permutation has the same phosphorylation sites as the best one,
        // except for the current site (i) which should be different
        for (Size j = 0; j < best_peptide_sites.size(); ++j)
        {
          if (j == i)
          {
            if (find(permutations[rev->second].begin(), permutations[rev->second].end(), best_peptide_sites[j]) != permutations[rev->second].end())
            {
              peptide_not_found = true;
              break;
            }
            else
            {
              peptide_not_found = false;
            }
          }
          else
          {
            if (find(permutations[rev->second].begin(), permutations[rev->second].end(), best_peptide_sites[j]) == permutations[rev->second].end())
            {
              peptide_not_found = true;
              break;
            }
            else
            {
              peptide_not_found = false;
            }
          }
        }
      } while (peptide_not_found);

      // store permutation of peptide without the phospho site i (seq_2)
      sites[i].seq_2 = rev->second;

      // store phospho site location that is not contained in the best scoring (seq_1) but in seq_2.
      for (Size j = 0; j < permutations[sites[i].seq_2].size(); ++j)
      {
        if (find(permutations[sites[i].seq_1].begin(), permutations[sites[i].seq_1].end(), permutations[sites[i].seq_2][j]) == permutations[sites[i].seq_1].end())
        {
          sites[i].second = permutations[sites[i].seq_2][j];
          break;
        }
      }
    }

    // Find the peak depth that maximizes the score difference for each phosphorylation site
    // store peak depth that achieves maximum score difference between best and runner up for every phospho site.
    for (Size i = 0; i < sites.size(); ++i)
    {
      double maximum_score_difference = 0.0;
      sites[i].peak_depth = 1;
      vector<double>::const_iterator first_it = peptide_site_scores[sites[i].seq_1].begin();
      vector<double>::const_iterator second_it = peptide_site_scores[sites[i].seq_2].begin();
      
      for (Size depth = 1; second_it != peptide_site_scores[sites[i].seq_2].end(); ++second_it, ++first_it, ++depth)
      {
        double phospho_at_site_score = *first_it;
        double no_phospho_at_site_score = *second_it;
        double score_difference = phospho_at_site_score - no_phospho_at_site_score;
        
        if (score_difference > maximum_score_difference)
        {
          maximum_score_difference = score_difference;
          sites[i].peak_depth = depth;
          sites[i].ascore = maximum_score_difference;
        }
      }
    }
  }
  
  /**
   * @brief Computes the site-determining ions for phosphorylation site candidates
   * 
   * Site-determining ions are fragment ions that can distinguish between different phosphorylation
   * site localizations. This method identifies ions that are unique to each of the two best-scoring
   * permutations for a given phosphorylation site.
   * 
   * @param th_spectra Vector of theoretical spectra for all permutations
   * @param candidates The phosphorylation site candidates to evaluate
   * @param site_determining_ions Output vector to store the site-determining ions
   */
  void AScore::computeSiteDeterminingIons_(const vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, vector<PeakSpectrum>& site_determining_ions) const
  {
    site_determining_ions.clear();
    site_determining_ions.resize(2);
    
    PeakSpectrum spectrum_first = th_spectra[candidates.seq_1];
    PeakSpectrum spectrum_second = th_spectra[candidates.seq_2];
    
    // Find peaks that are unique to the first spectrum (not in the second)
    PeakSpectrum spectrum_first_diff;
    AScore::getSpectrumDifference_(
      spectrum_first.begin(), spectrum_first.end(),
      spectrum_second.begin(), spectrum_second.end(),
      std::inserter(spectrum_first_diff, spectrum_first_diff.begin()));
      
    PeakSpectrum spectrum_second_diff;
    // Find peaks that are unique to the second spectrum (not in the first)
    AScore::getSpectrumDifference_(
      spectrum_second.begin(), spectrum_second.end(),
      spectrum_first.begin(), spectrum_first.end(),
      std::inserter(spectrum_second_diff, spectrum_second_diff.begin()));
      
    OPENMS_LOG_DEBUG << spectrum_first_diff << std::endl;
    OPENMS_LOG_DEBUG << spectrum_second_diff << std::endl;
      
    site_determining_ions[0] = spectrum_first_diff;
    site_determining_ions[1] = spectrum_second_diff;
    
    site_determining_ions[0].sortByPosition();
    site_determining_ions[1].sortByPosition(); 
  }

  /**
   * @brief Counts the number of matched ions between a theoretical spectrum and experimental spectrum window
   * 
   * Compares a theoretical spectrum with an experimental spectrum window and counts
   * the number of matching peaks within the fragment mass tolerance.
   * 
   * @param th Theoretical spectrum (must be sorted by m/z)
   * @param window Experimental spectrum window (must be sorted by m/z)
   * @param depth Maximum number of peaks to consider from the window (peak depth)
   */
  Size AScore::numberOfMatchedIons_(const PeakSpectrum& th, const PeakSpectrum& window, Size depth) const
  {
    PeakSpectrum window_reduced = window;
    if (window_reduced.size() > depth)
    {
      window_reduced.resize(depth);
    }
    
    window_reduced.sortByPosition();
    Size matched_peaks(0);
    // Use different matching iterators based on the tolerance type (ppm or Da)
    if (fragment_tolerance_ppm_)
    {
      MatchedIterator<PeakSpectrum, PpmTrait> it(th, window_reduced, fragment_mass_tolerance_);
      for (; it != it.end(); ++it) ++matched_peaks;
    }
    else
    {
      MatchedIterator<PeakSpectrum, DaTrait> it(th, window_reduced, fragment_mass_tolerance_);
      for (; it != it.end(); ++it) ++matched_peaks;
    }
    
    return matched_peaks;
  }

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
   */
  double AScore::peptideScore_(const std::vector<double>& scores) const
  {
    OPENMS_PRECONDITION(scores.size() == 10, "Scores vector must contain a score for every peak level."); 
    return (scores[0] * 0.5
            + scores[1] * 0.75
            + scores[2]
            + scores[3]
            + scores[4]
            + scores[5]
            + scores[6] * 0.75
            + scores[7] * 0.5
            + scores[8] * 0.25
            + scores[9] * 0.25)
           / 7.0;
  }

  /**
   * @brief Checks if a residue is a PhosphoDecoy site (A, G, L)
   */
  bool AScore::isPhosphoDecoySite(char residue)
  {
    return (residue == 'A' || residue == 'G' || residue == 'L');
  }

  /**
   * @brief Checks if a residue is a standard phosphorylation site (S, T, Y)
   * 
   * These are the standard amino acids that can be phosphorylated in biological systems:
   * - Serine (S)
   * - Threonine (T)
   * - Tyrosine (Y) */
  bool AScore::isPhosphoSite(char residue)
  {
    return (residue == 'S' || residue == 'T' || residue == 'Y');
  }

  vector<Size> AScore::getSites_(const String& unmodified) const
  {
    vector<Size> ret;
    for (Size i = 0; i < unmodified.size(); ++i)
    {
      // Always include standard phosphorylation sites (S, T, Y)
      if (isPhosphoSite(unmodified[i]))
      {
        ret.push_back(i);
      }
      // Include PhosphoDecoy sites (A, G, L) if enabled
      else if (add_decoys_ && isPhosphoDecoySite(unmodified[i]))
      {
        ret.push_back(i);
      }
    }
    return ret;
  }

  /**
   * @brief Counts the number of phosphorylation events in a peptide sequence
   * 
   * Counts both regular phosphorylation events (Phospho) and decoy phosphorylation
   * events (PhosphoDecoy) if enabled.
   * 
   * @param sequence The peptide sequence as a string
   */
  Size AScore::numberOfPhosphoEvents_(const String& sequence) const 
  {
    Size cnt_phospho_events = 0;
    Size pos = 0;
    
    // Count regular phosphorylation events
    while ((pos = sequence.find("(Phospho)", pos)) != std::string::npos)
    {
      ++cnt_phospho_events;
      pos += 9; // Move past "(Phospho)"
    }
    
    // Count PhosphoDecoy events
    pos = 0; // Reset position for PhosphoDecoy search
    while ((pos = sequence.find("(PhosphoDecoy)", pos)) != std::string::npos)
    {
      ++cnt_phospho_events;
      pos += 14; // Move past "(PhosphoDecoy)"
    }
    
    OPENMS_LOG_DEBUG << "Found " << cnt_phospho_events << " phosphorylation events in sequence: " << sequence << std::endl;
    
    return cnt_phospho_events;
  }

  /**
   * @brief Computes all possible combinations of phosphorylation sites.
   *
   * This method generates all possible combinations of n_phosphorylation_events
   * chosen from the given sites vector. It uses an iterative bitmask approach
   * to generate combinations without recursion, which is more efficient and
   * avoids stack overflow for large inputs.
   *
   * Algorithm:
   * 1. Create a bitmask with k 1's followed by (n-k) 0's
   * 2. For each permutation of the bitmask:
   *    a. Select sites corresponding to 1's in the bitmask
   *    b. Add the selected sites as a new permutation
   * 3. Continue until all permutations are generated
   *
   * Early termination conditions:
   * - If the estimated number of permutations exceeds max_permutations_
   * - If more phosphorylation events are requested than available sites
   * - If the actual number of permutations exceeds max_permutations_ during iteration
   *
   * Special cases:
   * - n_phosphorylation_events = 0: Returns empty permutations
   * - n_phosphorylation_events = 1: Returns each site individually
   * - n_phosphorylation_events = sites.size(): Returns one permutation with all sites
   *
   * Time Complexity: O(C(n,k)) where C(n,k) is the binomial coefficient (n choose k)
   * Space Complexity: O(k * C(n,k)) where k is n_phosphorylation_events
   */
  vector<vector<Size>> AScore::computePermutations_(const vector<Size>& sites, Int n_phosphorylation_events) const
  {
    vector<vector<Size>> permutations;
    // Early termination if we can estimate that the number of permutations will exceed the maximum
    if (max_permutations_ > 0 && n_phosphorylation_events >= 1)
    {
      // Estimate the number of permutations using the binomial coefficient (n choose k)
      double estimated_permutations = 0.0;
      
      // Check if k > n (more phosphorylation events than sites)
      if ((Size)n_phosphorylation_events > sites.size()) // TODO: check how this can happen
      {
        OPENMS_LOG_DEBUG << "\tEarly termination: more phosphorylation events (" << n_phosphorylation_events 
                         << ") than available sites (" << sites.size() << ")" << std::endl;
        return permutations; // Return empty permutations to signal that the calculation should be aborted
      }
      else
      {
        try
        {
          estimated_permutations = boost::math::binomial_coefficient<double>((unsigned int)sites.size(), (unsigned int)n_phosphorylation_events);
        }
        catch (std::exception const& /*e*/) // Catch any exception, not just overflow
        {
          estimated_permutations = std::numeric_limits<double>::max();
        }
      }
      if (estimated_permutations > max_permutations_)
      {
        OPENMS_LOG_DEBUG << "\tEarly termination: estimated permutations (" << estimated_permutations << ") exceeds maximum (" << max_permutations_ << ")" << std::endl;
        return permutations; // Return empty permutations to signal that the calculation should be aborted
      }
    }
    
    if (n_phosphorylation_events == 0)
    {
      return permutations;
    }
    else if (n_phosphorylation_events == 1)
    {
      for (Size i = 0; i < sites.size(); ++i)
      {
        vector<Size> temp;
        temp.push_back(sites[i]);
        permutations.push_back(temp);
      }
      return permutations;
    }
    // All sites are phosphorylated? Return one permutation containing all sites at once.
    else if (sites.size() == (Size)n_phosphorylation_events)
    {
      permutations.push_back(sites);
      return permutations;
    }
    else
    {
      // Generate all n_phosphorylation_events sized combinations from sites using an iterative approach
      Size n = sites.size();
      Size k = n_phosphorylation_events;
      
      // Create a bitmask: k 1's followed by (n-k) 0's
      // This represents selecting k elements from n elements
      vector<bool> bitmask(n, false);
      fill(bitmask.begin(), bitmask.begin() + k, true);
      
      // Generate all combinations using the bitmask
      // std::prev_permutation generates all permutations in lexicographically decreasing order
      do 
      {
        // For each permutation of the bitmask, create a combination of sites
        vector<Size> combination;
        for (Size i = 0; i < n; ++i)
        {
          if (bitmask[i])  // If this position is selected in the bitmask
          {
            combination.push_back(sites[i]);  // Add the corresponding site
          }
        }
        permutations.push_back(combination);
        
        // Early termination check to avoid excessive computation
        if (max_permutations_ > 0 && permutations.size() > max_permutations_)
        {
          OPENMS_LOG_DEBUG << "\tEarly termination during iteration: current permutations (" << permutations.size() << ") exceeds maximum (" << max_permutations_ << ")" << std::endl;
          permutations.clear();  // Clear permutations to signal abortion
          return permutations;
        }
      }
      while (std::prev_permutation(bitmask.begin(), bitmask.end()));
      
      return permutations;
    }
  }
  
    
  /**
   * @brief Creates a variant of the peptide with all phosphorylations removed
   * 
   * Removes both regular Phospho and PhosphoDecoy modifications from the sequence.
   * 
   * @param sequence The peptide sequence as a string
   */
  AASequence AScore::removePhosphositesFromSequence_(const String& sequence) const 
  {
    String seq(sequence);
    // Remove both regular and decoy phosphorylation markers
    seq.substitute("(Phospho)", "");
    seq.substitute("(PhosphoDecoy)", "");
    AASequence without_phospho = AASequence::fromString(seq);
    
    return without_phospho;
  }
  
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
   */
  vector<PeakSpectrum> AScore::createTheoreticalSpectra_(const vector<vector<Size>>& permutations, const AASequence& seq_without_phospho) const
  {
    vector<PeakSpectrum> th_spectra;
    
    
    // Get the unmodified sequence as a string to check residue types
    String seq_string = seq_without_phospho.toUnmodifiedString();

    th_spectra.resize(permutations.size());
    for (Size i = 0; i < permutations.size(); ++i)
    {
      AASequence seq(seq_without_phospho);
      Size permu = 0;
      
      for (Size as = 0; as < seq.size(); ++as)
      {
        if (as == permutations[i][permu])
        {
          // Apply the appropriate modification based on residue type
          char residue = seq_string[as];
          
          // Determine which modification to apply based on residue type
          if (add_decoys_ && isPhosphoDecoySite(residue))
          {
            // This is a PhosphoDecoy site (A, G, L)
            seq.setModification(as, "PhosphoDecoy");
            OPENMS_LOG_DEBUG << "Set PhosphoDecoy modification at position " << as 
                            << " (residue " << residue << ")" << std::endl;
          }
          else if (residue == 'S' || residue == 'T' || residue == 'Y')
          {
            // This is a standard phosphorylation site (S, T, Y)
            seq.setModification(as, "Phospho");
            OPENMS_LOG_DEBUG << "Set Phospho modification at position " << as 
                            << " (residue " << residue << ")" << std::endl;
          }
          else
          {
            // This is neither a standard phosphorylation site nor a PhosphoDecoy site
            // This should not happen, but log it just in case
            OPENMS_LOG_WARN << "Attempted to set phosphorylation on non-phosphorylatable residue " 
                           << residue << " at position " << as << std::endl;
          }
          
          ++permu;
        }
        
        if (permu == permutations[i].size()) 
        {
          break;
        }
      }

      // we mono-charge spectra, generating b- and y-ions is the default behavior of the TSG
      spectrum_generator_.getSpectrum(th_spectra[i], seq, 1, 1);
      th_spectra[i].setName(seq.toString());
    }
    return th_spectra;
  }

  /**
   * @brief Picks the top 10 intensity peaks for each 100 Da window in the spectrum
   * 
   * Divides the spectrum into 100 Da windows and selects the 10 most intense peaks
   * from each window. This approach normalizes for intensity variations across
   * the m/z range and focuses on the most informative peaks.
   * 
   * @param real_spectrum The experimental MS/MS spectrum
   */
  std::vector<PeakSpectrum> AScore::peakPickingPerWindowsInSpectrum_(PeakSpectrum& real_spectrum) const
  {
    vector<PeakSpectrum> windows_top10;
    
    double spect_lower_bound = floor(real_spectrum.front().getMZ() / 100) * 100;
    double spect_upper_bound = ceil(real_spectrum.back().getMZ() / 100) * 100;
    
    Size number_of_windows = static_cast<Size>(ceil((spect_upper_bound - spect_lower_bound) / 100));
    windows_top10.resize(number_of_windows);
    
    PeakSpectrum::Iterator it_current_peak = real_spectrum.begin();
    Size window_upper_bound(spect_lower_bound + 100);
    
    for (Size current_window = 0; current_window < number_of_windows; ++current_window)
    {
      PeakSpectrum real_window;
      while ((it_current_peak < real_spectrum.end()) && ((*it_current_peak).getMZ() <= window_upper_bound))
      {
        real_window.push_back(*it_current_peak);
        ++it_current_peak;
      }
      
      real_window.sortByIntensity(true);
      for (Size i = 0; (i < 10) & (i < real_window.size()); ++i)
      {
        windows_top10[current_window].push_back(real_window[i]);
      }
      
      window_upper_bound += 100;
    }
    return windows_top10;
  }
  
  /**
   * @brief Calculates scores for each permutation at 10 different peak depths
   * 
   * For each theoretical spectrum (permutation), this method:
   * 1. Counts matched ions at each peak depth (1-10)
   * 2. Calculates the cumulative binomial probability
   * 3. Converts to a score using -10*log10(probability)
   * 
   * @param th_spectra Vector of theoretical spectra for all permutations
   * @param windows_top10 Vector of experimental spectra windows with top 10 peaks
   */
  std::vector<std::vector<double>> AScore::calculatePermutationPeptideScores_(vector<PeakSpectrum>& th_spectra, const vector<PeakSpectrum>& windows_top10) const
  {
    //prepare peak depth for all windows in the actual spectrum
    vector<vector<double>> permutation_peptide_scores(th_spectra.size());
    vector<vector<double>>::iterator site_score = permutation_peptide_scores.begin();
    
    // for each phospho site assignment
    for (vector<PeakSpectrum>::iterator it = th_spectra.begin(); it != th_spectra.end(); ++it, ++site_score)
    {
      // the number of theoretical peaks (all b- and y-ions) correspond to the number of trials N
      Size N = it->size();
      site_score->resize(10);
      for (Size i = 1; i <= 10; ++i)
      {
        Size n = 0;
        for (Size current_win = 0; current_win < windows_top10.size(); ++current_win) // count matched ions over all 100 Da windows
        {
          n += numberOfMatchedIons_(*it, windows_top10[current_win], i);
        }
        double p = static_cast<double>(i) * base_match_probability_;
        double cumulative_score = computeCumulativeScore_(N, n, p);

        //abs is used to avoid -0 score values
        (*site_score)[i - 1] = abs((-10.0 * log10(cumulative_score)));
      }
    }
    return permutation_peptide_scores;
  }
  
  /**
   * @brief Ranks permutations by their weighted peptide scores
   * 
   * Calculates the weighted peptide score for each permutation and
   * creates a multimap ranking them in ascending order.
   * @param peptide_site_scores Scores for each permutation at each peak depth
   */
  std::multimap<double, Size> AScore::rankWeightedPermutationPeptideScores_(const vector<vector<double>>& peptide_site_scores) const
  {
    multimap<double, Size> ranking;
    
    for (Size i = 0; i != peptide_site_scores.size(); ++i)
    {
      double weighted_score = peptideScore_(peptide_site_scores[i]);
      ranking.insert(pair<double, Size>(weighted_score, i));
    }
    
    return ranking;
  }
  
  /**
   * @brief Compares two m/z values using the fragment mass tolerance
   * 
   * Determines if two m/z values are within the specified tolerance:
   * - Returns -1 if mz1 < mz2 (outside tolerance)
   * - Returns 1 if mz1 > mz2 (outside tolerance)
   * - Returns 0 if they are within tolerance
   * 
   * @param mz1 First m/z value
   * @param mz2 Second m/z value
   */
  int AScore::compareMZ_(double mz1, double mz2) const
  {
    double tolerance = fragment_mass_tolerance_;        
    double error = mz1 - mz2;
    
    if (fragment_tolerance_ppm_)
    {
      double avg_mass = (mz1 + mz2) / 2.0;
      tolerance = tolerance * avg_mass / 1.e6;
    }
    
    if (error < -tolerance)
    { 
      return -1;
    }
    else if (error > tolerance)
    {
      return 1;
    }
    else 
    { 
      return 0;
    }
  }

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
  String AScore::generateProFormaString_(const AASequence& peptide, const std::map<Size, double>& ascores) const
  {
    // Get the unmodified sequence as a string
    String unmodified_str = peptide.toUnmodifiedString();
    
    // Create a map to store scores for each position
    std::map<Size, double> position_scores;
    
    // Convert AScores to probabilities (0-1 range) for each phosphorylation site
    for (const auto& site_pair : ascores)
    {
      Size position = site_pair.first;
      double ascore = site_pair.second;


      // AScore is a -log10 scale, convert to probability
      // Higher AScore means higher confidence in site localization
      // double probability = 1.0 - pow(10.0, -ascore); this would be the actual probability 
      // but will lead to all numbers very close to 1.0.
      double probability = 1.0 - (1.0 / (1.0 + ascore));

      // Cap probability between 0 and 1
      probability = std::max(0.0, std::min(1.0, probability));
      
      // Store the probability for this position
      position_scores[position] = probability;
    }
    
    // Build the ProForma string
    String result;
    for (Size i = 0; i < unmodified_str.size(); ++i)
    {
      result += unmodified_str[i];
      
      // Check if this position has a phosphorylation
      auto it = position_scores.find(i);
      if (it != position_scores.end())
      {
        // Determine the modification type based on the residue
        char residue = unmodified_str[i];
        String mod_name;
        
        if (isPhosphoSite(residue))
        {
          // Standard phosphorylation site (S, T, Y)
          mod_name = "Phospho";
        }
        else if (add_decoys_ && isPhosphoDecoySite(residue))
        {
          // PhosphoDecoy site (A, G, L)
          mod_name = "PhosphoDecoy";
        }
        else
        {
          continue; // Skip if not a valid phosphorylation site
        }
        
        // Format the score with 4 decimal places
        std::stringstream ss;
        ss << std::fixed << std::setprecision(4) << it->second;
        // Add the ProForma-like annotation
        result += "[" + mod_name + "|score=" + ss.str() + "]";
      }
    }
    
    return result;
  }
  
  /**
   * @brief Updates member variables from the parameter object
   * 
   * Called when parameters are changed to update the internal state of the object
   */
  void AScore::updateMembers_()
  {
    fragment_mass_tolerance_ = param_.getValue("fragment_mass_tolerance");
    fragment_tolerance_ppm_ = (param_.getValue("fragment_mass_unit") == "ppm");
    max_peptide_length_ = param_.getValue("max_peptide_length");
    max_permutations_ = param_.getValue("max_num_perm");
    add_decoys_ = param_.getValue("add_decoys").toBool();
  }
  
} // namespace OpenMS
