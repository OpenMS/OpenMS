// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar, Timo Sachsenberg $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <cmath>
#include <algorithm> //find
#include <boost/math/special_functions/binomial.hpp>

#include <iostream>

using namespace std;

namespace OpenMS
{

  AScore::AScore()
  {
  }

  AScore::~AScore()
  {
  }

  PeptideHit AScore::compute(PeptideHit & hit, RichPeakSpectrum & real_spectrum, double fragment_mass_tolerance, bool fragment_mass_unit_ppm) const
  {
    String with_phospho_string = hit.getSequence().toString();

    // count number of phosphorylation events
    Size number_of_phosphorylation_events = 0;
    for (Size i = with_phospho_string.find("Phospho"); i != std::string::npos; i = with_phospho_string.find("Phospho", i + 7))
    {
      ++number_of_phosphorylation_events;
    }

    // create variant of the peptide with all phosphorylations removed
    String without_phospho_str(hit.getSequence().toString());
    without_phospho_str.substitute("(Phospho)", "");
    AASequence without_phospho = AASequence::fromString(without_phospho_str);

    // create unmodified variant
    String umps = without_phospho.toUnmodifiedString();

    // count number of potential phospho sites
    Size number_of_STY = std::count(umps.begin(), umps.end(), 'S') + std::count(umps.begin(), umps.end(), 'T') + std::count(umps.begin(), umps.end(), 'Y');

    if (real_spectrum.empty() || number_of_phosphorylation_events == 1 || number_of_STY == 0)
    {
      return hit;
    }

    TheoreticalSpectrumGenerator spectrum_generator;
    vector<RichPeakSpectrum> th_spectra;

    // produce theoretical spectra with all combinations with the number of phosphorylation events
    if (number_of_STY < number_of_phosphorylation_events)
    {
      number_of_phosphorylation_events = number_of_STY;
    }

    // determine all phospho sites
    vector<Size> sites(getSites(without_phospho));

    vector<vector<Size> > permutations(computePermutations(sites, number_of_phosphorylation_events));

    th_spectra.resize(permutations.size());
    for (Size i = 0; i < permutations.size(); ++i)
    {
      AASequence temp(without_phospho);
      Size permu = 0;
      for (Size as = 0; as < temp.size(); ++as)
      {
        if (as == permutations[i][permu])
        {
          temp.setModification(as, "Phospho");
          ++permu;
        }
        if (permu == permutations[i].size())
          break;
      }

      // we mono-charge spectra
      spectrum_generator.addPeaks(th_spectra[i], temp, Residue::BIon, 1);
      spectrum_generator.addPeaks(th_spectra[i], temp, Residue::YIon, 1);
      th_spectra[i].setName(temp.toString());
    }

    if (!real_spectrum.isSorted())
    {
      real_spectrum.sortByPosition();
    }

    vector<vector<double> > peptide_site_scores(th_spectra.size());
    vector<RichPeakSpectrum> windows_top10;

    // prepare peak depth for all windows in the actual spectrum
    {
      double biggest_window = real_spectrum.back().getMZ();
      Size number_of_windows = static_cast<Size>(ceil(biggest_window / 100));

      windows_top10.resize(number_of_windows);

      // for each of the 100 Da windows, retain the top 10 intensity peaks
      RichPeakSpectrum::Iterator begin_window = real_spectrum.MZBegin(0);
      RichPeakSpectrum::Iterator end_window = real_spectrum.MZBegin(100);
      for (Size current_window = 0; current_window < number_of_windows; ++current_window)
      {
        RichPeakSpectrum real_window;
        RichPeakSpectrum::Iterator runner = begin_window;
        while (runner <= end_window)
        {
          real_window.push_back(*runner);
          ++runner;
        }

        real_window.sortByIntensity(true);
        for (Size i = 0; i < 10; ++i)
        {
          if (i < real_window.size())
            windows_top10[current_window].push_back(real_window[i]);
        }
        begin_window = end_window + 1.0;
        end_window = real_spectrum.MZBegin((current_window + 1.0) * 100.0);
      }
    }

    //prepare peak depth for all windows in the actual spectrum
    vector<vector<double> >::iterator site_score = peptide_site_scores.begin();

    // for each phospho site assignment
    for (vector<RichPeakSpectrum>::iterator it = th_spectra.begin(); it < th_spectra.end(); ++it, ++site_score)
    {     
      // the number of theoretical peaks (all b- and y-ions) correspond to the number of trials N
      Size N = it->size();
      site_score->resize(10);
      for (Size i = 1; i <= 10; ++i)
      {
        UInt n = 0;
        for (Size depth = 0; depth < windows_top10.size(); ++depth) // count matched ions over all 100 Da windows
        {
          n += numberOfMatchedIons(*it, windows_top10[depth], i, fragment_mass_tolerance, fragment_mass_unit_ppm);
        }
        double p = static_cast<double>(i) / 100.0;
        double cumulative_score = computeCumulativeScore(N, n, p);
        (*site_score)[i - 1] = (-10.0 * log10(cumulative_score));
      }
    }

    vector<ProbablePhosphoSites> phospho_sites;

    PeptideHit phospho;

    // (permutations[0].size() == permutations.size() is the special case, that there are as many phosphorylation events as phospho sites.
    // No selection needs to be performed in this case.
    if (permutations[0].size() < permutations.size())
    {
      determineHighestScoringPermutations(peptide_site_scores, phospho_sites, permutations);

      // initialize score
      phospho.setScore(peptideScore(peptide_site_scores[phospho_sites[0].seq_1]));
    }

    {
      // initialize score with highest peptide score (aka highest weighted score)
      multimap<double, Size> ranking;
      for (Size i = 0; i != peptide_site_scores.size(); ++i)
      {
        double weighted_score = peptideScore(peptide_site_scores[i]);
        ranking.insert(pair<double, Size>(weighted_score, i));
      }
      phospho.setScore(ranking.rbegin()->first);
      phospho.setSequence(AASequence::fromString(th_spectra[ranking.rbegin()->second].getName()));
    }

    phospho.setCharge(hit.getCharge());
    phospho.setMetaValue("Search_engine_sequence", hit.getSequence().toString());

    Int rank = 1;
    for (vector<ProbablePhosphoSites>::iterator s_it = phospho_sites.begin(); s_it < phospho_sites.end(); ++s_it)
    {
      vector<RichPeakSpectrum> site_determining_ions;
      // Previously, the precursor charge was used here. This is clearly wrong and it is better to use charge 1 here.
      computeSiteDeterminingIons(th_spectra, *s_it, 1, site_determining_ions);
      Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
      double p = static_cast<double>(s_it->peak_depth) / 100.0;

      Size n_first = 0;
      for (Size depth = 0; depth != windows_top10.size(); ++depth) // for each 100 m/z window
      {
        n_first += numberOfMatchedIons(site_determining_ions[0], windows_top10[depth], depth, fragment_mass_tolerance, fragment_mass_unit_ppm);
      }

      double P_first = computeCumulativeScore(N, n_first, p);

      Size n_second = 0;
      for (Size depth = 0; depth <  windows_top10.size(); ++depth) //each 100 m/z window
      {
        n_second += numberOfMatchedIons(site_determining_ions[1], windows_top10[depth], depth, fragment_mass_tolerance, fragment_mass_unit_ppm);
      }

      double P_second = computeCumulativeScore(N, n_second, p);
      double score_first = -10 * log10(P_first);
      double score_second = -10 * log10(P_second);
      double AScore_first = score_first - score_second;
      phospho.setMetaValue("AScore_" + String(rank), AScore_first);
      ++rank;
    }
    return phospho;
  }

  double AScore::computeCumulativeScore(Size N, Size n, double p) const
  {
    OPENMS_PRECONDITION(n <= N, "The number of matched ions (n) can be at most as large as the number of trials (N).");
    OPENMS_PRECONDITION(p >= 0 && p <= 1.0, "p must be a probability [0,1].");

    // return bad p value if none has been matched (see Beausoleil et al.)
    if (n == 0) return 1.0;

    double score = 0.0;
    // score = sum_{k=n..N}(\choose{N}{k}p^k(1-p)^{N-k})
    for (Size k = n; k <= N; ++k)
    {
      double coeff = boost::math::binomial_coefficient<double>((double)N, k);
      double pow1 = pow((double)p, (int)k);
      double pow2 = pow(double(1 - p), double(N - k));
      score += coeff * pow1 * pow2;
    }

    return score;
  }

  void AScore::determineHighestScoringPermutations(const std::vector<std::vector<double> >& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const vector<vector<Size> >& permutations) const
  {
    // For every phospho site of the highest (weighted) scoring phospho site assignment:
    // 1. determine the next best (weighted) score assignment with this site in unphosporylated state.
    // 2. determine the filtering level (peak depths) that maximizes the (unweighted) score difference between these two assignments

    sites.clear();

    // take first set of phospho site assignments
    sites.resize(permutations[0].size());

    multimap<double, Size> ranking;
    for (Size permutation_idx = 0; permutation_idx != peptide_site_scores.size(); ++permutation_idx)  // for each scored phospho site assignment (=permutation)
    {
      // calculate weighted score by linearly weighting all depths (top 1 - top 10 intensity levels).
      double current_score = peptideScore(peptide_site_scores[permutation_idx]);

      // store weighted score and permutation
      ranking.insert(pair<double, Size>(current_score, permutation_idx));
    }

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
      }
      while (peptide_not_found);

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

      for (Size depth = 1; second_it < peptide_site_scores[sites[i].seq_2].end(); ++second_it, ++first_it, ++depth)
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

  void AScore::computeSiteDeterminingIons(vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, Int charge, vector<RichPeakSpectrum> & site_determining_ions) const
  {
    site_determining_ions.clear();
    site_determining_ions.resize(2);
    TheoreticalSpectrumGenerator spectrum_generator;
    AASequence pref, suf, pref_with_phospho_first, pref_with_phospho_second, suf_with_phospho_first, suf_with_phospho_second;
    RichPeakSpectrum prefix, suffix, prefix_with_phospho_first, prefix_with_phospho_second, suffix_with_phospho_second, suffix_with_phospho_first;
    Size permutation_with_site_phosphorylated = candidates.seq_1;
    Size permutation_with_site_not_phosphorylated = candidates.seq_2;
    AASequence first(AASequence::fromString(th_spectra[permutation_with_site_phosphorylated].getName()));
    AASequence second(AASequence::fromString(th_spectra[permutation_with_site_not_phosphorylated].getName()));

    if (candidates.first < candidates.second)
    {
      pref = AASequence::fromString(first.getPrefix(candidates.first + 1).toString());
      suf = AASequence::fromString(second.getSuffix(second.size() - candidates.second - 1).toString());
      pref_with_phospho_first = AASequence::fromString(first.getPrefix(candidates.second + 1).toString());
      pref_with_phospho_second = AASequence::fromString(second.getPrefix(candidates.second + 1).toString());
      suf_with_phospho_first = AASequence::fromString(first.getSuffix(first.size() - candidates.first).toString());
      suf_with_phospho_second = AASequence::fromString(second.getSuffix(second.size() - candidates.first).toString());
    }
    else
    {
      pref = AASequence::fromString(second.getPrefix(candidates.second + 1).toString());
      suf = AASequence::fromString(first.getSuffix(first.size() - candidates.first - 1).toString());
      pref_with_phospho_first = AASequence::fromString(first.getPrefix(candidates.first + 1).toString());
      pref_with_phospho_second = AASequence::fromString(second.getPrefix(candidates.first + 1).toString());
      suf_with_phospho_first = AASequence::fromString(first.getSuffix(first.size() - candidates.second).toString());
      suf_with_phospho_second = AASequence::fromString(second.getSuffix(first.size() - candidates.second).toString());
    }

    spectrum_generator.addPeaks(prefix, pref, Residue::BIon, charge);
    spectrum_generator.addPeaks(suffix, suf, Residue::YIon, charge);
    spectrum_generator.addPeaks(prefix_with_phospho_first, pref_with_phospho_first, Residue::BIon, charge);
    spectrum_generator.addPeaks(prefix_with_phospho_second, pref_with_phospho_second, Residue::BIon, charge);
    spectrum_generator.addPeaks(suffix_with_phospho_first, suf_with_phospho_first, Residue::YIon, charge);
    spectrum_generator.addPeaks(suffix_with_phospho_second, suf_with_phospho_second, Residue::YIon, charge);

    if (!prefix.empty())
    {
      for (RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
      {
        if (it->getMZ() > prefix.back().getMZ())
        {
          site_determining_ions[0].push_back(*it);
        }
      }
      for (RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
      {
        if (it->getMZ() > prefix.back().getMZ())
        {
          site_determining_ions[1].push_back(*it);
        }
      }
    }
    else
    {
      for (RichPeakSpectrum::iterator it = prefix_with_phospho_first.begin(); it < prefix_with_phospho_first.end(); ++it)
      {
        site_determining_ions[0].push_back(*it);
      }
      for (RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
      {
        site_determining_ions[1].push_back(*it);
      }
    }
    if (!suffix.empty())
    {
      for (RichPeakSpectrum::iterator it = suffix_with_phospho_first.begin(); it < suffix_with_phospho_first.end(); ++it)
      {
        if (it->getMZ() > suffix.back().getMZ())
        {
          site_determining_ions[0].push_back(*it);
        }
      }
      for (RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
      {
        if (it->getMZ() > suffix.back().getMZ())
        {
          site_determining_ions[1].push_back(*it);
        }
      }
    }
    else
    {
      RichPeakSpectrum::iterator it1 = suffix_with_phospho_first.begin();
      RichPeakSpectrum::iterator it2 = suffix_with_phospho_second.begin();
      if (!suf.empty())
      {
        ++it1;
        ++it2;
      }
      for (; it1 < suffix_with_phospho_first.end(); ++it1)
      {
        site_determining_ions[0].push_back(*it1);
      }
      for (; it2 < suffix_with_phospho_second.end(); ++it2)
      {
        site_determining_ions[1].push_back(*it2);
      }
    }
    site_determining_ions[0].sortByPosition();
    site_determining_ions[1].sortByPosition();
  }

  Size AScore::numberOfMatchedIons(const RichPeakSpectrum & th, const RichPeakSpectrum & windows, Size depth, double fragment_mass_tolerance, bool fragment_mass_tolerance_ppm) const
  {
    Size n = 0;
    for (Size i = 0; i < windows.size() && i <= depth; ++i)
    {
      Size nearest_peak = th.findNearest(windows[i].getMZ());
      if (nearest_peak < th.size())
      {
        double theo_mz = th[nearest_peak].getMZ();
        double error = abs(theo_mz - windows[i].getMZ());

        if (fragment_mass_tolerance_ppm)
        {
          error = error / theo_mz * 1e6;
        }

        if (error < fragment_mass_tolerance)
        {
          ++n;
        }
      }
    }
    return n;
  }

  double AScore::peptideScore(const std::vector<double> & scores) const
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
           / 10.0;
  }

  vector<Size> AScore::getSites(AASequence& without_phospho) const
  {
    vector<Size> tupel;
    String unmodified = without_phospho.toUnmodifiedString();
    for (Size i = 0; i < unmodified.size(); ++i)
    {
      if (unmodified[i] == 'Y' || unmodified[i] == 'T' || unmodified[i] == 'S')
      {
        tupel.push_back(i);
      }
    }
    return tupel;
  }

  vector<vector<Size> > AScore::computePermutations(vector<Size> sites, Int n_phosphorylation_events) const
  {
    // Only one phosphorylation event? Then return every site is a permutation.
    if (n_phosphorylation_events == 1)
    {
      vector<vector<Size>  > permutations;
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
      vector<vector<Size> > permutations;
      permutations.push_back(sites);
      return permutations;
    }
    else
    // Generate all n_phosphorylation_events sized sets from sites
    {
      vector<vector<Size> > permutations;
      vector<Size> head;
      vector<vector<Size> > tail;
      // all permutations with first site selected
      head.push_back(sites[0]);
      vector<Size> tupel_left(++sites.begin(), sites.end());
      Int tail_phospho_sites = n_phosphorylation_events - 1;
      tail = computePermutations(tupel_left, tail_phospho_sites);
      for (vector<vector<Size> >::iterator it = tail.begin(); it < tail.end(); ++it)
      {
        vector<Size> temp(head);
        temp.insert(temp.end(), it->begin(), it->end());
        permutations.push_back(temp);
      }

      // all permutations with first site not selected
      vector<vector<Size> > other_possibilities(computePermutations(tupel_left, n_phosphorylation_events));
      permutations.insert(permutations.end(), other_possibilities.begin(), other_possibilities.end());
      return permutations;
    }
  }
} // namespace OpenMS
