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
  vector<int> NuXLSiteLocalization::supportedResidues(const vector<double>& n_shifts,
                                                      const vector<double>& n_noshifts,
                                                      const vector<double>& c_shifts,
                                                      const vector<double>& c_noshifts,
                                                      double tolerance)
  {
    const int length = static_cast<int>(n_shifts.size());
    // unshifted a1/b1 ions are not used: b1 of lysine (m/z 129.10) is also a common lysine ion of any peptide
    // (e.g. y1-H2O of a C-terminal lysine), whether the lysine is cross-linked or not
    const double total_intensity = accumulate(n_shifts.begin(), n_shifts.end(), 0.0)
      + (length > 1 ? accumulate(n_noshifts.begin() + 1, n_noshifts.end(), 0.0) : 0.0)
      + accumulate(c_shifts.begin(), c_shifts.end(), 0.0) + accumulate(c_noshifts.begin(), c_noshifts.end(), 0.0);
    if (length == 0 || total_intensity <= 0.0) { return {}; } // no fragment evidence for any site

    // intensity of the ions consistent with a cross-link at residue k minus the inconsistent ones.
    // Unshifted ions alone also locate it: they can't contain the cross-linked residue.
    vector<double> support(length, 0.0);
    for (int k = 0; k != length; ++k)
    {
      for (int i = 0; i != length; ++i)
      {
        support[k] += (i >= k ? n_shifts[i] : -n_shifts[i]);     // prefix ending at i contains k: must be shifted
        if (i > 0) { support[k] += (i < k ? n_noshifts[i] : -n_noshifts[i]); } // prefix ending at i doesn't contain k: must be unshifted
        support[k] += (i <= k ? c_shifts[i] : -c_shifts[i]);     // suffix starting at i contains k: must be shifted
        support[k] += (i > k ? c_noshifts[i] : -c_noshifts[i]);  // suffix starting at i doesn't contain k: must be unshifted
      }
    }

    // Residues not separated by any ion get exactly the same support (e.g. the first two residues without a1/b1
    // and y(n-1) ions). Support differences below the tolerance (a fraction of all fragment intensity) are
    // not conclusive either: a single weak or misassigned ion would otherwise decide.
    const double best = *max_element(support.begin(), support.end());
    const double margin = max(tolerance, 1e-9) * total_intensity;
    vector<int> supported;
    for (int k = 0; k != length; ++k)
    {
      if (support[k] >= best - margin) { supported.push_back(k); }
    }
    return supported;
  }

  // static
  NuXLSiteLocalization::Site NuXLSiteLocalization::localizeNHS(const string& sequence,
                                                               bool reactive_protein_n_term,
                                                               const vector<double>& n_shifts,
                                                               const vector<double>& n_noshifts,
                                                               const vector<double>& c_shifts,
                                                               const vector<double>& c_noshifts,
                                                               double tolerance,
                                                               const vector<bool>& blocked)
  {
    Site site;
    const auto reactive = [&blocked](int p) { return p >= static_cast<int>(blocked.size()) || !blocked[p]; };
    // contradicting ions can support residues on both sides of unsupported ones: only supported residues are candidates
    vector<int> by_distance = supportedResidues(n_shifts, n_noshifts, c_shifts, c_noshifts, tolerance);
    if (by_distance.empty() || by_distance.back() >= static_cast<int>(sequence.size())) { return site; }
    const int begin = by_distance.front(), end = by_distance.back();
    site.range_begin = begin;
    site.range_end = end;

    // 1. free protein N-terminus
    if (reactive_protein_n_term && begin == 0)
    {
      site.position = 0;
      site.method = "N-term";
      return site;
    }

    // supported residues ordered by distance to the middle of their range (ties: closer to the N-terminus)
    const int middle = begin + (end - begin + 2) / 2 - 1;
    stable_sort(by_distance.begin(), by_distance.end(), [middle](int a, int b) { return abs(a - middle) < abs(b - middle); });

    // 2. lysine
    for (int p : by_distance)
    {
      if (sequence[p] == 'K' && reactive(p))
      {
        site.position = p;
        site.method = "K";
        return site;
      }
    }

    // 3. serine, threonine or tyrosine with a histidine within four residues
    for (int p : by_distance)
    {
      if ((sequence[p] != 'S' && sequence[p] != 'T' && sequence[p] != 'Y') || !reactive(p)) { continue; }
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
  NuXLSiteLocalization::LossAwareSite NuXLSiteLocalization::localizeLossAware(const vector<double>& n_shifts,
                                                                             const vector<double>& n_noshifts,
                                                                             const vector<double>& c_shifts,
                                                                             const vector<double>& c_noshifts,
                                                                             const vector<double>& residue_evidence,
                                                                             const vector<bool>& blocked,
                                                                             double unshifted_weight)
  {
    LossAwareSite site;
    const int length = static_cast<int>(n_shifts.size());
    if (length == 0) { return site; }

    // votes[k]: prefix rungs via a difference array (a prefix ending at i contains residues 0..i,
    // a suffix starting at i contains residues i..length-1)
    vector<double> scores(length, 0.0);
    double shifted_total = 0.0;
    for (int i = 0; i != length; ++i)
    {
      if (n_shifts[i] > 0.0) // shifted prefix: for 0..i, against i+1..
      {
        const double w = sqrt(n_shifts[i]);
        shifted_total += w;
        for (int k = 0; k != length; ++k) { scores[k] += (k <= i ? w : -w); }
      }
      else if (n_noshifts[i] > 0.0) // unshifted prefix only: against 0..i
      {
        const double w = unshifted_weight * sqrt(n_noshifts[i]);
        for (int k = 0; k <= i; ++k) { scores[k] -= w; }
      }

      if (c_shifts[i] > 0.0) // shifted suffix: for i.., against 0..i-1
      {
        const double w = sqrt(c_shifts[i]);
        shifted_total += w;
        for (int k = 0; k != length; ++k) { scores[k] += (k >= i ? w : -w); }
      }
      else if (c_noshifts[i] > 0.0) // unshifted suffix only: against i..
      {
        const double w = unshifted_weight * sqrt(c_noshifts[i]);
        for (int k = i; k != length; ++k) { scores[k] -= w; }
      }
    }

    double evidence_max = 0.0;
    for (int k = 0; k != length && k < static_cast<int>(residue_evidence.size()); ++k)
    {
      scores[k] += residue_evidence[k];
      evidence_max = max(evidence_max, residue_evidence[k]);
    }
    if (shifted_total <= 0.0 && evidence_max <= 0.0) { return site; } // nothing carries the cross-link

    const double norm = shifted_total > 0.0 ? shifted_total : evidence_max;
    for (double& s : scores) { s /= norm; }

    const auto allowed = [&blocked](int k) { return k >= static_cast<int>(blocked.size()) || !blocked[k]; };
    bool found = false;
    double best = 0.0;
    for (int k = 0; k != length; ++k)
    {
      if (allowed(k) && (!found || scores[k] > best)) { best = scores[k]; found = true; }
    }
    if (!found) { return site; }

    const double eps = 1e-9;
    double second = best;
    bool has_second = false;
    for (int k = 0; k != length; ++k)
    {
      if (!allowed(k)) { continue; }
      if (scores[k] >= best - eps) { site.best.push_back(k); }
      else if (!has_second || scores[k] > second) { second = scores[k]; has_second = true; }
    }
    site.position = site.best.front();
    site.range_begin = site.best.front();
    site.range_end = site.best.back();
    site.score = best;
    site.delta = has_second ? best - second : 0.0;
    site.scores = std::move(scores);
    return site;
  }

  // static
  void NuXLSiteLocalization::digestProtein(const EnzymaticDigestion& digestion,
                                           string_view protein,
                                           size_t min_length,
                                           size_t max_length,
                                           bool clip_nterm_methionine,
                                           vector<string_view>& peptides)
  {
    digestion.digestUnmodified(protein, peptides, min_length, max_length);
    if (!clip_nterm_methionine || protein.size() < 2 || protein[0] != 'M') { return; }

    // initial methionine removed: the N-terminal peptides of the rest of the protein
    vector<string_view> clipped;
    digestion.digestUnmodified(protein.substr(1), clipped, min_length, max_length);
    for (const string_view& p : clipped)
    {
      // views into the protein: N-terminal products of the clipped protein start right after the methionine
      if (p.data() == protein.data() + 1) { peptides.push_back(p); }
    }
  }

  // static
  NuXLSiteLocalization::ProteinTermini NuXLSiteLocalization::proteinTermini(const vector<string_view>& proteins,
                                                                            const EnzymaticDigestion& digestion,
                                                                            size_t min_length,
                                                                            size_t max_length,
                                                                            bool clip_nterm_methionine)
  {
    ProteinTermini termini;
    vector<string_view> peptides;
    for (const string_view& protein : proteins)
    {
      digestProtein(digestion, protein, min_length, max_length, clip_nterm_methionine, peptides);
      for (const string_view& p : peptides)
      {
        const bool clipped_start = clip_nterm_methionine && protein[0] == 'M' && p.data() == protein.data() + 1;
        if (p.data() == protein.data() || clipped_start) { termini.n_terminal.emplace_back(p); }
        if (p.data() + p.size() == protein.data() + protein.size()) { termini.c_terminal.emplace_back(p); }
      }
    }
    for (vector<string>* v : {&termini.n_terminal, &termini.c_terminal})
    {
      sort(v->begin(), v->end());
      v->erase(unique(v->begin(), v->end()), v->end());
    }
    return termini;
  }

  // static
  bool NuXLSiteLocalization::isProteinNTerminal(const ProteinTermini& termini, string_view peptide)
  {
    return binary_search(termini.n_terminal.begin(), termini.n_terminal.end(), peptide,
      [](const auto& a, const auto& b) { return string_view(a) < string_view(b); });
  }

  // static
  bool NuXLSiteLocalization::isProteinCTerminal(const ProteinTermini& termini, string_view peptide)
  {
    return binary_search(termini.c_terminal.begin(), termini.c_terminal.end(), peptide,
      [](const auto& a, const auto& b) { return string_view(a) < string_view(b); });
  }
} // namespace OpenMS
