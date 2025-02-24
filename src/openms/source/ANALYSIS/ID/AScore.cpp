// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Petra Gutenbrunner $
// $Authors: David Wojnar, Timo Sachsenberg, Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AScore.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
#include <OpenMS/KERNEL/RangeUtils.h>

#include <boost/math/special_functions/binomial.hpp>

using namespace std;

namespace OpenMS
{

  AScore::AScore():
    DefaultParamHandler("AScore")
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

    defaultsToParam_();
  }

  AScore::~AScore() = default;

  PeptideHit AScore::compute(const PeptideHit& hit, PeakSpectrum& real_spectrum)
  {
    PeptideHit phospho = hit;
    
    //reset phospho
    phospho.setScore(-1);
    if (real_spectrum.empty())
    {
      return phospho;
    }
    
    String sequence_str = phospho.getSequence().toString();
    String unmodified_sequence_str = phospho.getSequence().toUnmodifiedString();
    
    Size number_of_phosphorylation_events = numberOfPhosphoEvents_(sequence_str);
    AASequence seq_without_phospho = removePhosphositesFromSequence_(sequence_str);

    if ((max_peptide_length_ > 0) && (unmodified_sequence_str.size() > max_peptide_length_))
    {
      OPENMS_LOG_DEBUG << "\tcalculation aborted: peptide too long: " << seq_without_phospho.toString() << std::endl;
      return phospho;
    }

    // determine all phospho sites
    vector<Size> sites = getSites_(unmodified_sequence_str);
    Size number_of_STY = sites.size();

    if (number_of_phosphorylation_events == 0 || number_of_STY == 0)
    {
      return phospho;
    }
    if (number_of_STY == number_of_phosphorylation_events)
    {
      phospho.setScore(unambiguous_score_);
      return phospho;
    } 

    vector<vector<Size>> permutations = computePermutations_(sites, (Int)number_of_phosphorylation_events);
    OPENMS_LOG_DEBUG << "\tnumber of permutations: " << permutations.size() << std::endl;

    // TODO: using a heuristic to calculate the best phospho sites if the number of permutations are exceeding the maximum.
    // A heuristic could be to calculate the best site for the first phosphorylation and based on this the best site for the second 
    // phosphorylation and so on until every site is determined
    if ((max_permutations_ > 0) && (permutations.size() > max_permutations_))
    {
      OPENMS_LOG_DEBUG << "\tcalculation aborted: number of permutations exceeded" << std::endl;
      return phospho;
    }

    vector<PeakSpectrum> th_spectra = createTheoreticalSpectra_(permutations, seq_without_phospho);

    // prepare real spectrum windows
    if (!real_spectrum.isSorted())
    {
      real_spectrum.sortByPosition();
    }
    vector<PeakSpectrum> windows_top10 = peakPickingPerWindowsInSpectrum_(real_spectrum);

    // compute match probability for a peak depth of 1
    base_match_probability_ = computeBaseProbability_(real_spectrum.back().getMZ());

    // calculate peptide score for each possible phospho site permutation
    vector<vector<double>> peptide_site_scores = calculatePermutationPeptideScores_(th_spectra, windows_top10);

    // rank peptide permutations ascending
    multimap<double, Size> ranking = rankWeightedPermutationPeptideScores_(peptide_site_scores);

    multimap<double, Size>::reverse_iterator rev = ranking.rbegin();
    String seq1 = th_spectra[rev->second].getName();
    phospho.setSequence(AASequence::fromString(seq1));
    phospho.setMetaValue("search_engine_sequence", hit.getSequence().toString());

    double peptide1_score = rev->first;
    phospho.setMetaValue("AScore_pep_score", peptide1_score); // initialize score with highest peptide score (aka highest weighted score)

    ++rev;
    String seq2 = th_spectra[rev->second].getName();
    double peptide2_score = rev->first;

    vector<ProbablePhosphoSites> phospho_sites;
    determineHighestScoringPermutations_(peptide_site_scores, phospho_sites, permutations, ranking);

    Int rank = 1;
    double best_Ascore = std::numeric_limits<double>::max(); // the lower the better
    for (vector<ProbablePhosphoSites>::iterator s_it = phospho_sites.begin(); s_it != phospho_sites.end(); ++s_it)
    {
      double Ascore = 0;
      if (peptide1_score == peptide2_score) // set Ascore = 0 for each phosphorylation site
      {
        OPENMS_LOG_DEBUG << "\tscore of best (" << seq1 << ") and second best peptide (" << seq2 << ") are equal (" << peptide1_score << ")" << std::endl;
      }
      else
      {
        vector<PeakSpectrum> site_determining_ions;

        computeSiteDeterminingIons_(th_spectra, *s_it, site_determining_ions);
        Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
        double p = static_cast<double>(s_it->peak_depth) * base_match_probability_;

        Size n_first = 0; // number of matching peaks for first peptide
        for (Size window_idx = 0; window_idx != windows_top10.size(); ++window_idx) // for each 100 m/z window
        {
          n_first += numberOfMatchedIons_(site_determining_ions[0], windows_top10[window_idx], s_it->peak_depth);        
        }
        double P_first = computeCumulativeScore_(N, n_first, p);

        Size n_second = 0; // number of matching peaks for second peptide
        for (Size window_idx = 0; window_idx <  windows_top10.size(); ++window_idx) //each 100 m/z window
        {
          n_second += numberOfMatchedIons_(site_determining_ions[1], windows_top10[window_idx], s_it->peak_depth);        
        }
        Size N2 = site_determining_ions[1].size(); // all possibilities have the same number so take the first one
        double P_second = computeCumulativeScore_(N2, n_second, p);

        //abs is used to avoid -0 score values
        double score_first = abs(-10 * log10(P_first));
        double score_second = abs(-10 * log10(P_second));

        OPENMS_LOG_DEBUG << "\tfirst - N: " << N << ",p: " << p << ",n: " << n_first << ", score: " << score_first << std::endl;
        OPENMS_LOG_DEBUG << "\tsecond - N: " << N2 << ",p: " << p << ",n: " << n_second << ", score: " << score_second << std::endl;

        Ascore = score_first - score_second;
        OPENMS_LOG_DEBUG << "\tAscore_" << rank << ": " << Ascore << std::endl;
      }
      if (Ascore < best_Ascore)
      {
        best_Ascore = Ascore;
      }
      phospho.setMetaValue("AScore_" + String(rank), Ascore);
      ++rank;      
    }
    phospho.setScore(best_Ascore);
    return phospho;
  }

  double AScore::computeBaseProbability_(double ppm_reference_mz) const
  {
    double base_match_probability = 2. * fragment_mass_tolerance_ / 100.;
    if (fragment_tolerance_ppm_)
    {
      base_match_probability *= ppm_reference_mz * 1e-6; // 1e-6 converts fragment_mass_tolerance_ to ppm
    }
    return base_match_probability;
  }

  double AScore::computeCumulativeScore_(Size N, Size n, double p) const
  {
    OPENMS_PRECONDITION(n <= N, "The number of matched ions (n) can be at most as large as the number of trials (N).");
    OPENMS_PRECONDITION(p >= 0 && p <= 1.0, "p must be a probability [0,1].");

    // return bad p value if none has been matched (see Beausoleil et al.)
    if (n == 0) return 1.0;

    double score = 0.0;
    // score = sum_{k=n..N}(\choose{N}{k}p^k(1-p)^{N-k})
    for (Size k = n; k <= N; ++k)
    {
      double coeff = 0;

      try
      {
        coeff = boost::math::binomial_coefficient<double>((unsigned int)N, (unsigned int)k);
      }
      catch (std::overflow_error const& /*e*/)
      {
        // not sure if a warning is appropriate here, since if it happens, it will happen very often for the same spectrum and flood the stdout
//        std::cout << "Warning: Binomial coefficient for AScore has overflowed! Setting value to the maximal double value." << std::endl;
//        std::cout << "binomial_coefficient was called with N = " << N << " and k = " << k << std::endl;
        coeff = std::numeric_limits<double>::max();
      }

      double pow1 = pow((double)p, (int)k);
      double pow2 = pow(double(1 - p), double(N - k));

      score += coeff * pow1 * pow2;
    }

    return score;
  }

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
      multimap<double, Size>::reverse_iterator rev = ranking.rbegin();      
      sites[i].first = best_peptide_sites[i]; // store the site
      sites[i].seq_1 = rev->second; // and permutation
      bool peptide_not_found = true;
      
      // iterate from best scoring peptide to the first peptide that doesn't contain the current phospho site
      do
      {      
        ++rev;
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
        }
      }
    }
  }
  
  // calculation of the number of different speaks between the theoretical spectra of the two best scoring peptide permutations, respectively
  void AScore::computeSiteDeterminingIons_(const vector<PeakSpectrum>& th_spectra, const ProbablePhosphoSites& candidates, vector<PeakSpectrum>& site_determining_ions) const
  {
    site_determining_ions.clear();
    site_determining_ions.resize(2);
    
    PeakSpectrum spectrum_first = th_spectra[candidates.seq_1];
    PeakSpectrum spectrum_second = th_spectra[candidates.seq_2];
    
    PeakSpectrum spectrum_first_diff;
    AScore::getSpectrumDifference_(
      spectrum_first.begin(), spectrum_first.end(),
      spectrum_second.begin(), spectrum_second.end(),
      std::inserter(spectrum_first_diff, spectrum_first_diff.begin()));
      
    PeakSpectrum spectrum_second_diff;
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

  Size AScore::numberOfMatchedIons_(const PeakSpectrum& th, const PeakSpectrum& window, Size depth) const
  {
    PeakSpectrum window_reduced = window;
    if (window_reduced.size() > depth)
    {
      window_reduced.resize(depth);
    }
    
    window_reduced.sortByPosition();
    Size matched_peaks(0);
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

  vector<Size> AScore::getSites_(const String& unmodified) const
  {
    vector<Size> tupel;
    for (Size i = 0; i < unmodified.size(); ++i)
    {
      if (unmodified[i] == 'Y' || unmodified[i] == 'T' || unmodified[i] == 'S')
      {
        tupel.push_back(i);
      }
    }
    return tupel;
  }

  vector<vector<Size>> AScore::computePermutations_(const vector<Size>& sites, Int n_phosphorylation_events) const
  {
    vector<vector<Size>  > permutations;
    
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
    // Generate all n_phosphorylation_events sized sets from sites
    {
      vector<Size> head;
      vector<vector<Size>> tail;
      
      // all permutations with first site selected
      head.push_back(sites[0]);
      vector<Size> tupel_left(++sites.begin(), sites.end());
      Int tail_phospho_sites = n_phosphorylation_events - 1;
      
      tail = computePermutations_(tupel_left, tail_phospho_sites);
      
      for (vector<vector<Size>>::iterator it = tail.begin(); it != tail.end(); ++it)
      {
        vector<Size> temp(head);
        temp.insert(temp.end(), it->begin(), it->end());
        permutations.push_back(temp);
      }

      // all permutations with first site not selected
      vector<vector<Size>> other_possibilities(computePermutations_(tupel_left, n_phosphorylation_events));
      permutations.insert(permutations.end(), other_possibilities.begin(), other_possibilities.end());
      return permutations;
    }
  }
  
  /// Computes number of phospho events in a sequence
  Size AScore::numberOfPhosphoEvents_(const String& sequence) const 
  {
    Size cnt_phospho_events = 0;
    
    for (Size i = sequence.find("Phospho"); i != std::string::npos; i = sequence.find("Phospho", i + 7))
    {
      ++cnt_phospho_events;
    }
    
    return cnt_phospho_events;
  }
    
  /// Create variant of the peptide with all phosphorylations removed
  AASequence AScore::removePhosphositesFromSequence_(const String& sequence) const 
  {
    String seq(sequence);
    seq.substitute("(Phospho)", "");
    AASequence without_phospho = AASequence::fromString(seq);
    
    return without_phospho;
  }
  
  /// Create theoretical spectra
  vector<PeakSpectrum> AScore::createTheoreticalSpectra_(const vector<vector<Size>>& permutations, const AASequence& seq_without_phospho) const
  {
    vector<PeakSpectrum> th_spectra;
    TheoreticalSpectrumGenerator spectrum_generator;
    
    th_spectra.resize(permutations.size());
    for (Size i = 0; i < permutations.size(); ++i)
    {
      AASequence seq(seq_without_phospho);
      Size permu = 0;
      
      for (Size as = 0; as < seq.size(); ++as)
      {
        if (as == permutations[i][permu])
        {
          seq.setModification(as, "Phospho");
          ++permu;
        }
        
        if (permu == permutations[i].size()) 
        {
          break;
        }
      }

      // we mono-charge spectra, generating b- and y-ions is the default behavior of the TSG
      spectrum_generator.getSpectrum(th_spectra[i], seq, 1, 1);
      th_spectra[i].setName(seq.toString());
    }
    return th_spectra;
  }

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
  
  int AScore::compareMZ_(double mz1, double mz2) const
  {
    double tolerance = fragment_mass_tolerance_;        
    double error = mz1 - mz2;
    
    if (fragment_tolerance_ppm_)
    {
      double avg_mass = (mz1 + mz2) / 2;
      tolerance = tolerance * avg_mass / 1e6;
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

  void AScore::updateMembers_()
  {
    fragment_mass_tolerance_ = param_.getValue("fragment_mass_tolerance");
    fragment_tolerance_ppm_ = (param_.getValue("fragment_mass_unit") == "ppm");
    max_peptide_length_ = param_.getValue("max_peptide_length");
    max_permutations_ = param_.getValue("max_num_perm");
    unambiguous_score_ = param_.getValue("unambiguous_score");
  }
  
} // namespace OpenMS
