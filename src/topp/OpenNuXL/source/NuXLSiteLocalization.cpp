// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "NuXLSiteLocalization.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>

using namespace std;

namespace OpenMS
{
  // static
  pair<int, int> NuXLSiteLocalization::supportedRange(const vector<double>& n_shifts,
                                                      const vector<double>& n_noshifts,
                                                      const vector<double>& c_shifts,
                                                      const vector<double>& c_noshifts)
  {
    const int length = static_cast<int>(n_shifts.size());
    const double shifted_intensity = accumulate(n_shifts.begin(), n_shifts.end(), 0.0) + accumulate(c_shifts.begin(), c_shifts.end(), 0.0);
    if (length == 0 || shifted_intensity <= 0.0) { return {-1, -1}; } // no fragment evidence for any site

    // intensity of the ions consistent with a cross-link at residue k minus the inconsistent ones
    vector<double> support(length, 0.0);
    for (int k = 0; k != length; ++k)
    {
      for (int i = 0; i != length; ++i)
      {
        support[k] += (i >= k ? n_shifts[i] : -n_shifts[i]);     // prefix ending at i contains k: must be shifted
        support[k] += (i < k ? n_noshifts[i] : -n_noshifts[i]);  // prefix ending at i doesn't contain k: must be unshifted
        support[k] += (i <= k ? c_shifts[i] : -c_shifts[i]);     // suffix starting at i contains k: must be shifted
        support[k] += (i > k ? c_noshifts[i] : -c_noshifts[i]);  // suffix starting at i doesn't contain k: must be unshifted
      }
    }

    // residues not separated by any ion get exactly the same terms: allow for rounding only.
    // E.g. without b1 and y(n-1) ions, the first two residues always tie.
    const double best = *max_element(support.begin(), support.end());
    const double total_intensity = shifted_intensity
      + accumulate(n_noshifts.begin(), n_noshifts.end(), 0.0) + accumulate(c_noshifts.begin(), c_noshifts.end(), 0.0);
    const double tolerance = 1e-9 * total_intensity;
    int begin = -1, end = -1;
    for (int k = 0; k != length; ++k)
    {
      if (support[k] < best - tolerance) { continue; }
      if (begin == -1) { begin = k; }
      end = k;
    }
    return {begin, end};
  }

  // static
  NuXLSiteLocalization::Site NuXLSiteLocalization::localizeNHS(const string& sequence,
                                                               bool reactive_protein_n_term,
                                                               const vector<double>& n_shifts,
                                                               const vector<double>& n_noshifts,
                                                               const vector<double>& c_shifts,
                                                               const vector<double>& c_noshifts)
  {
    Site site;
    const auto [begin, end] = supportedRange(n_shifts, n_noshifts, c_shifts, c_noshifts);
    if (begin < 0 || end >= static_cast<int>(sequence.size())) { return site; }
    site.range_begin = begin;
    site.range_end = end;

    // 1. free protein N-terminus
    if (reactive_protein_n_term && begin == 0)
    {
      site.position = 0;
      site.method = "N-term";
      return site;
    }

    // residues of the range ordered by distance to its middle (ties: closer to the N-terminus)
    const int middle = begin + (end - begin + 2) / 2 - 1;
    vector<int> by_distance(end - begin + 1);
    iota(by_distance.begin(), by_distance.end(), begin);
    stable_sort(by_distance.begin(), by_distance.end(), [middle](int a, int b) { return abs(a - middle) < abs(b - middle); });

    // 2. lysine
    for (int p : by_distance)
    {
      if (sequence[p] == 'K')
      {
        site.position = p;
        site.method = "K";
        return site;
      }
    }

    // 3. serine, threonine or tyrosine with a histidine within four residues
    for (int p : by_distance)
    {
      if (sequence[p] != 'S' && sequence[p] != 'T' && sequence[p] != 'Y') { continue; }
      for (int h = max(0, p - 4); h <= min(static_cast<int>(sequence.size()) - 1, p + 4); ++h)
      {
        if (h != p && sequence[h] == 'H')
        {
          site.position = p;
          site.method = "STY+H";
          return site;
        }
      }
    }

    // 4. no residue an NHS ester reacts with
    return site;
  }

  // static
  vector<string> NuXLSiteLocalization::proteinNTermini(const vector<string_view>& protein_sequences)
  {
    vector<string> n_termini;
    n_termini.reserve(2 * protein_sequences.size());
    for (const string_view& s : protein_sequences)
    {
      n_termini.emplace_back(s);
      if (!s.empty() && s[0] == 'M') { n_termini.emplace_back(s.substr(1)); } // initiator methionine removed
    }
    sort(n_termini.begin(), n_termini.end());
    return n_termini;
  }

  // static
  bool NuXLSiteLocalization::isProteinNTerminal(const vector<string>& protein_n_termini, string_view peptide)
  {
    // the smallest entry not less than the peptide is the one to start with it, if any does
    const auto it = lower_bound(protein_n_termini.begin(), protein_n_termini.end(), peptide,
      [](const string& entry, string_view value) { return string_view(entry) < value; });
    return it != protein_n_termini.end() && string_view(*it).substr(0, peptide.size()) == peptide;
  }
} // namespace OpenMS
