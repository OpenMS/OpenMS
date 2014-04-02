// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
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

  PeptideHit AScore::compute(PeptideHit & hit, RichPeakSpectrum & real_spectrum, double fmt, Int number_of_phospho_sites)
  {
    String without_phospho_str(hit.getSequence().toString());
    size_t found = without_phospho_str.find("(Phospho)");
    while (found != string::npos)
    {
      without_phospho_str.erase(found, String("(Phospho)").size());
      found = without_phospho_str.find("(Phospho)");
    }
    AASequence without_phospho(without_phospho_str);
    Int number_of_STY = Int(without_phospho.getNumberOf("S") + without_phospho.getNumberOf("T") + without_phospho.getNumberOf("Y"));
    if (real_spectrum.empty() || number_of_phospho_sites < 1 || number_of_STY == 0)
    {
      return PeptideHit(-1, 0, hit.getCharge(), without_phospho);
    }
    TheoreticalSpectrumGenerator spectrum_generator;
    vector<RichPeakSpectrum> th_spectra;     //typedef MSSpectrum<RichPeak1D> RichPeakSpectrum;
    ///produce theoretical spectra
    if (number_of_STY < number_of_phospho_sites)
    {
      number_of_phospho_sites = number_of_STY;
    }
    //Int number_of_permutations = (Int)boost::math::binomial_coefficient<double>((double)number_of_STY, number_of_phospho_sites);
    vector<Size> tupel(computeTupel_(without_phospho));
    vector<vector<Size> > permutations(computePermutations_(tupel, number_of_phospho_sites));
    //cout<<"number of permutations "<<permutations.size();//<< " - " << number_of_permutations;
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
      spectrum_generator.addPeaks(th_spectra[i], temp, Residue::BIon, hit.getCharge());
      spectrum_generator.addPeaks(th_spectra[i], temp, Residue::YIon, hit.getCharge());
      th_spectra[i].setName(temp.toString());
    }
    ///produce theoretical spectra - END

    if (!real_spectrum.isSorted())
    {
      real_spectrum.sortByPosition();
    }
    vector<vector<double> > peptide_site_scores(th_spectra.size());
    vector<RichPeakSpectrum> windows_with_all_peak_depths;
    ///prepare peak depth for all windows in the actual spectrum
    {
      double biggest_window = real_spectrum[real_spectrum.size() - 1].getMZ();
      Size number_of_windows = 1;
      if (biggest_window > 100)
        number_of_windows = ceil(biggest_window / 100);
      windows_with_all_peak_depths.resize(number_of_windows);
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
            windows_with_all_peak_depths[current_window].push_back(real_window[i]);
        }
        begin_window = ++end_window;
        end_window = real_spectrum.MZBegin((current_window + 1) * 100);
      }
    }
    ///prepare peak depth for all windows in the actual spectrum - END
    UInt N;
    vector<vector<double> >::iterator side_scores = peptide_site_scores.begin();
    for (vector<RichPeakSpectrum>::iterator it = th_spectra.begin(); it < th_spectra.end(); ++it, ++side_scores)    //each theoretical spectrum
    {
      N = UInt(it->size());       //real or theo!!
      UInt i = 1;
      side_scores->resize(10);
      while (i <= 10)
      {
        //Auslagern
        UInt n = 0;
        for (Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth)         //each 100 m/z window
        {
          n += numberOfMatchedIons_(*it, windows_with_all_peak_depths[depth], i, fmt);
        }
        //auslagern_end
        double p = (double)i / 100;
        double cum_socre = computeCumulativeScore(N, n, p);
        (*side_scores)[i - 1] = (-10 * log10(cum_socre));    //computeCumulativeScore(N,n,p);
        ++i;
      }
    }
    vector<ProbablePhosphoSites> highest_peptides;
    PeptideHit phospho;
    if (permutations[0].size() < permutations.size())
    {
      computeHighestPeptides(peptide_site_scores, highest_peptides, permutations);
      phospho.setScore(peptideScore_(peptide_site_scores[highest_peptides[0].seq_1]));
      //phospho.setSequence(AASequence(th_spectra[highest_peptides[0].seq_1].getName()));
    }
    {
      multimap<double, Size> ranking;
      for (Size it = 0; it < peptide_site_scores.size(); ++it)
      {
        double current_score = peptideScore_(peptide_site_scores[it]);
        ranking.insert(pair<double, Size>(current_score, it));
      }
      phospho.setScore(ranking.rbegin()->first);
      phospho.setSequence(AASequence(th_spectra[ranking.rbegin()->second].getName()));
    }
    phospho.setCharge(hit.getCharge());
    phospho.setMetaValue("Search_engine_sequence", hit.getSequence().toString());
    ///Calculate AScore
    //Int add_ascores = (Int)highest_peptides.size()/permutations[0].size();
    //Int count_ascores = 0;
    //vector<ProbablePhosphoSites>::iterator winner = highest_peptides.begin();
    //vector<ProbablePhosphoSites>::iterator actual = winner;
    //double ascore_addition = 0.0;
    //double highest_addition = 0.0;
    Int rank = 1;
    for (vector<ProbablePhosphoSites>::iterator hp = highest_peptides.begin(); hp < highest_peptides.end(); ++hp)
    {
      vector<RichPeakSpectrum> site_determining_ions;
      compute_site_determining_ions(th_spectra, *hp, hit.getCharge(), site_determining_ions);
      N = UInt(site_determining_ions[0].size());       // all possiblities have the same number
      double p = (double)hp->peak_depth / 100;
      Int n = 0;
      //Auslagern
      for (Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth)       //each 100 m/z window
      {
        n += numberOfMatchedIons_(site_determining_ions[0], windows_with_all_peak_depths[depth], (Size)p, fmt);
      }
      //auslagern_end
      double P_first = computeCumulativeScore(N, n, p);
      Int n2 = 0;
      //Auslagern
      for (Size depth = 0; depth <  windows_with_all_peak_depths.size(); ++depth)       //each 100 m/z window
      {
        n2 += numberOfMatchedIons_(site_determining_ions[1], windows_with_all_peak_depths[depth], (Size)p, fmt);
      }
      //auslagern_end
      double P_second = computeCumulativeScore(N, n2, p);
      double score_first = -10 * log10(P_first);
      double score_second = -10 * log10(P_second);
      double AScore_first = score_first - score_second;
      phospho.setMetaValue("AScore_" + String(rank), AScore_first);
      ++rank;
      /*    hp->AScore = AScore_first;
          ++count_ascores;
          ascore_addition += AScore_first;
          if(count_ascores == add_ascores)
          {

              count_ascores = 0;
              if(ascore_addition > highest_addition)
              {
                  highest_addition = ascore_addition;
                  winner = actual;
              }
              actual = hp;
              ascore_addition = 0;
          }*/
    }
    return phospho;
  }

  double AScore::computeCumulativeScore(UInt N, UInt n, double p)
  {
    if (n > N)
    {
      return -1.0;
    }
    double score = 0.0;
    for (UInt k = n; k <= N; ++k)
    {
      double coeff = boost::math::binomial_coefficient<double>((double)N, k);
      double pow1 = pow((double)p, (int)k);
      double pow2 = pow(double(1 - p), double(N - k));
      score += coeff * pow1 * pow2;
    }
    if (score == 0.0)
    {
      return 1.0;
    }
    return score;
  }

  void AScore::computeHighestPeptides(std::vector<std::vector<double> > & peptide_site_scores, std::vector<ProbablePhosphoSites> & sites, vector<vector<Size> > & permutations)
  {
    sites.clear();
    sites.resize(permutations[0].size());
    multimap<double, Size> ranking;
    for (Size it = 0; it < peptide_site_scores.size(); ++it)
    {
      double current_score = peptideScore_(peptide_site_scores[it]);
      ranking.insert(pair<double, Size>(current_score, it));
    }
    vector<Size> & hps = permutations[ranking.rbegin()->second /*it->second*/];       //highest peptide score
    for (Size i = 0; i < hps.size(); ++i)
    {
      multimap<double, Size>::reverse_iterator rev = ranking.rbegin();
      sites[i].first = hps[i];
      sites[i].seq_1 = rev->second;          //it->second;
      bool peptide_not_found = true;
      do
      {
        ++rev;
        for (Size j = 0; j < hps.size(); ++j)
        {
          if (j == i)
          {
            if (find(permutations[rev->second].begin(), permutations[rev->second].end(), hps[j]) != permutations[rev->second].end())
            {
              peptide_not_found = true;
              break;
            }
            else
            {
              peptide_not_found = false;
            }
          }
          else if (find(permutations[rev->second].begin(), permutations[rev->second].end(), hps[j]) == permutations[rev->second].end())
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
      while (peptide_not_found);
      sites[i].seq_2 = rev->second;
      for (Size j = 0; j < permutations[sites[i].seq_2].size(); ++j)
      {
        if (find(permutations[sites[i].seq_1].begin(), permutations[sites[i].seq_1].end(), permutations[sites[i].seq_2][j]) == permutations[sites[i].seq_1].end())
        {
          sites[i].second = permutations[sites[i].seq_2][j];
          break;
        }
      }
    }
    for (Size i = 0; i < sites.size(); ++i)
    {
      double current_peak_depth = 0.0;
      sites[i].peak_depth = 1;
      vector<double>::iterator first_it = peptide_site_scores[sites[i].seq_1].begin();
      Size depth = 1;
      for (vector<double>::iterator second_it = peptide_site_scores[sites[i].seq_2].begin(); second_it < peptide_site_scores[sites[i].seq_2].end(); ++second_it, ++first_it)
      {
        if ((*first_it - *second_it) > current_peak_depth)
        {
          current_peak_depth = *first_it - *second_it;
          sites[i].peak_depth = depth;
        }
        ++depth;
      }
    }
  }

  void AScore::compute_site_determining_ions(vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, Int charge, vector<RichPeakSpectrum> & site_determining_ions)
  {
    site_determining_ions.clear();
    site_determining_ions.resize(2);
    TheoreticalSpectrumGenerator spectrum_generator;
    AASequence pref, suf, pref_with_phospho_first, pref_with_phospho_second, suf_with_phospho_first, suf_with_phospho_second;
    RichPeakSpectrum prefix, suffix, prefix_with_phospho_first, prefix_with_phospho_second, suffix_with_phospho_second, suffix_with_phospho_first;
    AASequence first(th_spectra[candidates.seq_1].getName());
    AASequence second(th_spectra[candidates.seq_2].getName());
    if (candidates.first < candidates.second)
    {
      pref.setStringSequence(first.getPrefix(candidates.first + 1).toString());
      suf.setStringSequence(second.getSuffix(second.size() - candidates.second - 1).toString());
      pref_with_phospho_first.setStringSequence(first.getPrefix(candidates.second + 1).toString());
      pref_with_phospho_second.setStringSequence(second.getPrefix(candidates.second + 1).toString());
      suf_with_phospho_first.setStringSequence(first.getSuffix(first.size() - candidates.first).toString());
      suf_with_phospho_second.setStringSequence(second.getSuffix(second.size() - candidates.first).toString());
    }
    else
    {
      pref.setStringSequence(second.getPrefix(candidates.second + 1).toString());
      suf.setStringSequence(first.getSuffix(first.size() - candidates.first - 1).toString());
      pref_with_phospho_first.setStringSequence(first.getPrefix(candidates.first + 1).toString());
      pref_with_phospho_second.setStringSequence(second.getPrefix(candidates.first + 1).toString());
      suf_with_phospho_first.setStringSequence(first.getSuffix(first.size() - candidates.second).toString());
      suf_with_phospho_second.setStringSequence(second.getSuffix(first.size() - candidates.second).toString());
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
        if (it->getMZ() > prefix[prefix.size() - 1].getMZ())
          site_determining_ions[0].push_back(*it);
      }
      for (RichPeakSpectrum::iterator it = prefix_with_phospho_second.begin(); it < prefix_with_phospho_second.end(); ++it)
      {
        if (it->getMZ() > prefix[prefix.size() - 1].getMZ())
          site_determining_ions[1].push_back(*it);
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
        if (it->getMZ() > suffix[suffix.size() - 1].getMZ())
          site_determining_ions[0].push_back(*it);
      }
      for (RichPeakSpectrum::iterator it = suffix_with_phospho_second.begin(); it < suffix_with_phospho_second.end(); ++it)
      {
        if (it->getMZ() > suffix[suffix.size() - 1].getMZ())
          site_determining_ions[1].push_back(*it);
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

  Int AScore::numberOfMatchedIons_(const RichPeakSpectrum & th, const RichPeakSpectrum & windows, Size depth, double fmt)
  {
    Int n = 0;
    for (Size i = 0; i < windows.size() && i <= depth; ++i)
    {
      Size nearest_peak = th.findNearest(windows[i].getMZ());
      if (nearest_peak < th.size() && fabs(th[nearest_peak].getMZ() - windows[i].getMZ()) < fmt)
        ++n;
    }
    return n;
  }

  double AScore::peptideScore_(std::vector<double> & scores)
  {
    return (scores[0] * 0.5
            + scores[1] * 0.75
            + scores[2]           //*1
            + scores[3]           //*1
            + scores[4]           //*1
            + scores[5]           //*1
            + scores[6] * 0.75
            + scores[7] * 0.5
            + scores[8] * 0.25
            + scores[9] * 0.25)
           / 10;
  }

  vector<Size> AScore::computeTupel_(AASequence & without_phospho)
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

  vector<vector<Size> > AScore::computePermutations_(vector<Size> tupel, Int number_of_phospho_sites)
  {
    if (number_of_phospho_sites == 1)
    {
      vector<vector<Size>  > permutations;
      for (Size i = 0; i < tupel.size(); ++i)
      {
        vector<Size> temp;
        temp.push_back(tupel[i]);
        permutations.push_back(temp);
      }
      return permutations;
    }
    else if (tupel.size() == (Size)number_of_phospho_sites)
    {
      vector<vector<Size> > permutations;
      permutations.push_back(tupel);
      return permutations;
    }
    else
    {
      vector<vector<Size> > permutations;
      vector<Size> head;
      vector<vector<Size> > tail;
      head.push_back(tupel[0]);
      vector<Size> tupel_left(++tupel.begin(), tupel.end());
      Int tail_phospho_sites = number_of_phospho_sites - 1;
      tail = computePermutations_(tupel_left, tail_phospho_sites);
      for (vector<vector<Size> >::iterator it = tail.begin(); it < tail.end(); ++it)
      {
        vector<Size> temp(head);
        temp.insert(temp.end(), it->begin(), it->end());
        permutations.push_back(temp);
      }
      vector<vector<Size> > other_possibilities(computePermutations_(tupel_left, number_of_phospho_sites));
      permutations.insert(permutations.end(), other_possibilities.begin(), other_possibilities.end());
      return permutations;
    }
  }

} // namespace OpenMS
