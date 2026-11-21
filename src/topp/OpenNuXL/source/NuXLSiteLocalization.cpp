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
                                                                             double unshifted_weight,
                                                                             double secondary_min_gain)
  {
    LossAwareSite site;
    const int length = static_cast<int>(n_shifts.size());
    if (length == 0) { return site; }

    // the rungs that vote: a prefix ending at i contains residues 0..i, a suffix starting at i contains i..length-1;
    // an unshifted ion on a rung with a shifted one doesn't vote (loss of the cross-link)
    struct Rung
    {
      int lo;       // first residue the ion contains
      int hi;       // last residue the ion contains
      double vote;  // > 0
      bool shifted; // true: for the residues inside, against those outside; false: against those inside
    };
    vector<Rung> rungs;
    for (int i = 0; i != length; ++i)
    {
      if (n_shifts[i] > 0.0) { rungs.push_back({0, i, sqrt(n_shifts[i]), true}); }
      else if (n_noshifts[i] > 0.0 && unshifted_weight > 0.0) { rungs.push_back({0, i, unshifted_weight * sqrt(n_noshifts[i]), false}); }
      if (c_shifts[i] > 0.0) { rungs.push_back({i, length - 1, sqrt(c_shifts[i]), true}); }
      else if (c_noshifts[i] > 0.0 && unshifted_weight > 0.0) { rungs.push_back({i, length - 1, unshifted_weight * sqrt(c_noshifts[i]), false}); }
    }
    const auto contains = [](const Rung& r, int k) { return r.lo <= k && k <= r.hi; };

    vector<double> scores(length, 0.0);
    double shifted_total = 0.0;
    for (const Rung& r : rungs)
    {
      if (r.shifted) { shifted_total += r.vote; }
      for (int k = 0; k != length; ++k)
      {
        if (contains(r, k)) { scores[k] += r.shifted ? r.vote : -r.vote; }
        else if (r.shifted) { scores[k] -= r.vote; }
      }
    }

    const auto evidence = [&residue_evidence](int k) { return k < static_cast<int>(residue_evidence.size()) ? residue_evidence[k] : 0.0; };
    const auto same = [](double a, double b) { return std::fabs(a - b) <= 1e-9 * max(1.0, max(std::fabs(a), std::fabs(b))); };
    double evidence_max = 0.0;
    for (int k = 0; k != length; ++k)
    {
      scores[k] += evidence(k);
      evidence_max = max(evidence_max, evidence(k));
    }
    if (shifted_total <= 0.0 && evidence_max <= 0.0) { return site; } // nothing carries the cross-link

    const double norm = shifted_total > 0.0 ? shifted_total : evidence_max;
    for (double& s : scores) { s /= norm; }

    const auto allowed = [&blocked](int k) { return k >= static_cast<int>(blocked.size()) || !blocked[k]; };
    const double eps = 1e-9;
    int top = -1;
    for (int k = 0; k != length; ++k)
    {
      if (allowed(k) && (top < 0 || scores[k] > scores[top] + eps)) { top = k; } // N-terminal-most of the best
    }
    if (top < 0) { return site; }

    // tied with the localized residue: no rung contains one but not the other, and the same residue evidence
    const auto tied = [&](int k)
    {
      for (const Rung& r : rungs) { if (contains(r, k) != contains(r, top)) { return false; } }
      return same(evidence(k), evidence(top));
    };
    double second = 0.0;
    bool has_second = false;
    for (int k = 0; k != length; ++k)
    {
      if (!allowed(k)) { continue; }
      if (tied(k)) { site.best.push_back(k); }
      else if (!has_second || scores[k] > second) { second = scores[k]; has_second = true; }
    }
    site.position = top;
    site.range_begin = site.best.front();
    site.range_end = site.best.back();
    site.score = scores[top];
    site.delta = has_second ? scores[top] - second : 0.0;
    site.scores = std::move(scores);

    if (secondary_min_gain <= 0.0) { return site; }

    // a second cross-linked form at j, given the one at top: shifted rungs that contain j but not top count for
    // instead of against, unshifted-only rungs that contain top but not j no longer count against, and the
    // residue evidence counts at the higher of the two. j needs evidence of its own to be told apart from top.
    int secondary = -1;
    double secondary_gain = 0.0;
    for (int j = 0; j != length; ++j)
    {
      if (!allowed(j) || tied(j)) { continue; }
      double explained = 0.0;
      double relieved = 0.0;
      for (const Rung& r : rungs)
      {
        if (r.shifted && contains(r, j) && !contains(r, top)) { explained += r.vote; }
        else if (!r.shifted && contains(r, top) && !contains(r, j)) { relieved += r.vote; }
      }
      const double extra_evidence = same(evidence(j), evidence(top)) ? 0.0 : max(0.0, evidence(j) - evidence(top));
      if (explained <= 0.0 && extra_evidence <= 0.0) { continue; }
      const double gain = (2.0 * explained + relieved + extra_evidence) / norm;
      if (secondary < 0 || gain > secondary_gain + eps) { secondary = j; secondary_gain = gain; } // N-terminal-most of the best
    }
    if (secondary < 0 || secondary_gain < secondary_min_gain - eps) { return site; }

    // tied with the secondary site, given top: the rungs that count for the pair (shifted ones without top,
    // unshifted-only ones with top) contain both or neither, and the residue evidence of the pair is the same
    const auto secondary_tied = [&](int k)
    {
      for (const Rung& r : rungs)
      {
        const bool counts_for_pair = r.shifted ? !contains(r, top) : contains(r, top);
        if (counts_for_pair && contains(r, k) != contains(r, secondary)) { return false; }
      }
      return same(max(evidence(k), evidence(top)), max(evidence(secondary), evidence(top)));
    };
    for (int k = 0; k != length; ++k)
    {
      if (allowed(k) && !tied(k) && secondary_tied(k)) { site.secondary_best.push_back(k); }
    }
    site.secondary = secondary;
    site.secondary_range_begin = site.secondary_best.front();
    site.secondary_range_end = site.secondary_best.back();
    site.secondary_gain = secondary_gain;
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
