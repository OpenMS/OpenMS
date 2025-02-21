// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/ConvolutionBasedProteinFilter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <utility>

namespace OpenMS
{
inline const Size max_node_cntr = 300;

std::vector<Residue> FLASHTaggerAlgorithm::getAA_(double l, double r, double tol, int consider_ion_diff, int mode) const
{
  std::vector<Residue> ret;
  if (l == r) return ret;

  const auto& shifts = mode == 0 ? common_shifts_ : (mode == 1 ? n_term_shifts_ : c_term_shifts_);

  bool connected = false;
  for (const auto shift : shifts)
  {
    if (consider_ion_diff == 0 && shift != 0) continue;
    if (consider_ion_diff != 0 && shift == 0) continue;
    double diff1 = std::abs(std::abs(r - l) - shift);
    double diff2 = std::abs(std::abs(r - l) + shift);
    double abs_tol = std::max(l, r) * tol / 1e6 * 2;
    auto iter = aa_mass_map_.lower_bound(diff1 - abs_tol);

    while (iter != aa_mass_map_.end())
    {
      if (std::abs(diff1 - iter->first) < abs_tol || std::abs(diff2 - iter->first) < abs_tol)
      {
        for (const auto& aa : iter->second)
        {
          ret.push_back(aa);
          connected = true;
        }
      }
      else if (iter->first - diff2 > abs_tol) { break; }
      iter++;
    }
    if (connected) break;
  }
  return ret;
}

std::vector<std::vector<Residue>> FLASHTaggerAlgorithm::getGap_(double l, double r, double tol, int iso_offset) const
{
  std::vector<std::vector<Residue>> ret;
  if (l == r) return ret;
  double iso_mass = std::abs(iso_offset * Constants::C13C12_MASSDIFF_U);

  double diff1 = std::abs(std::abs(r - l) - iso_mass);
  double diff2 = std::abs(std::abs(r - l) + iso_mass);
  double abs_tol = 2 * std::abs(std::max(l, r) * tol / 1e6);
  auto iter = gap_mass_map_.lower_bound(diff1 - abs_tol);

  while (iter != gap_mass_map_.end())
  {
    if (std::abs(diff1 - iter->first) < abs_tol || std::abs(diff2 - iter->first) < abs_tol)
    {
      for (const auto& aa : iter->second)
        ret.push_back(aa);
    }
    else if (iter->first - diff2 > abs_tol) { break; }
    iter++;
  }
  return ret;
}

void FLASHTaggerAlgorithm::updateEdgeMasses_()
{
  aa_mass_map_.clear();
  gap_mass_map_.clear();

  for (const auto& aa : aas_)
  {
    double aa_mass = aa->getMonoWeight(Residue::Internal);
    if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end()) aa_mass_map_[aa_mass] = std::vector<Residue>();
    aa_mass_map_[aa_mass].push_back(*aa);
  }

  if (max_gap_count_ > 0)
  {
    std::map<double, std::vector<std::vector<Residue>>> prev_gap_mass_map_;
    prev_gap_mass_map_[.0] = std::vector<std::vector<Residue>>(1, std::vector<Residue>());
    for (int i = 0; i < max_aa_in_gap_; i++)
    {
      for (const auto& prev : prev_gap_mass_map_)
      {
        for (const auto& current : aa_mass_map_)
        {
          if (gap_mass_map_.find(prev.first + current.first) == gap_mass_map_.end())
            gap_mass_map_[prev.first + current.first] = std::vector<std::vector<Residue>>();

          for (const auto& aa_vec : prev.second)
          {
            for (const auto& aa : current.second)
            {
              auto new_aa_vec(aa_vec);
              new_aa_vec.push_back(aa);
              gap_mass_map_[prev.first + current.first].push_back(new_aa_vec);
            }
          }
        }
      }
      prev_gap_mass_map_ = gap_mass_map_;
    }

    std::map<double, std::vector<std::vector<Residue>>> tmp_gap_mass_map;
    for (const auto& e : gap_mass_map_)
    {
      std::vector<std::vector<Residue>> tmp_e;
      for (const auto& f : e.second)
      {
        if (f.size() <= 1) continue;
        tmp_e.push_back(f);
      }
      if (tmp_e.empty()) continue;
      tmp_gap_mass_map[e.first] = tmp_e;
    }
    gap_mass_map_ = tmp_gap_mass_map;
  }
}

inline Size FLASHTaggerAlgorithm::getVertex_(int index, int path_score, int level, int diff_ion_jump, int gap_level) const
{
  return (((index * (max_tag_length_ + 1) + level) * ((consider_diff_ion_jumps_ ? 1 : 0) + 1) + diff_ion_jump) * (max_gap_count_ + 1) + gap_level)
           * (max_path_score_ - min_path_score_ + 1)
         + (path_score - min_path_score_);
}

inline int FLASHTaggerAlgorithm::getIndex_(Size vertex) const
{
  return (((int)vertex / (max_path_score_ - min_path_score_ + 1)) / (max_gap_count_ + 1)) / ((consider_diff_ion_jumps_ ? 1 : 0) + 1)
         / (max_tag_length_ + 1);
}

void FLASHTaggerAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         int length,
                                         double tol,
                                         int mode)
{
  // from source to sink, connect but the edge direction is from sink to source.
  edge_aa_map_.clear();

  int start_index = 1; // zeroth = source.
  int end_index = 1;
  boost::dynamic_bitset<> visited(dag.size());
  visited[getVertex_(0, 0, 0, 0, 0)] = true;

  while (end_index < (int)mzs.size())
  {
    auto r = mzs[end_index];

    // first, make edge from r to source.
    Size vertex1 = getVertex_(end_index, scores[end_index], 0, 0, 0);
    Size vertex2 = getVertex_(0, 0, 0, 0, 0);

    bool connected = dag.addEdge(vertex1, vertex2, visited); // move to DAG?
    if (! connected)
    {
      end_index++;
      continue;
    }
    // from an edge i, j to class edge.  for each i, j make a unique key. key to an edge.

    while (start_index < end_index && r - mzs[start_index] > max_edge_mass_)
      start_index++;

    for (int g = 0; g < (max_gap_count_ > 0 ? 2 : 1); g++) // 0 for only a.a 1 for with gap.
    {
      for (int n = 0; n < (consider_diff_ion_jumps_ ? 2 : 1); n++) // 0 for all a.a 1 for diff ion type jumps. Allow only one isotope errors.
      {
        for (int current_index = start_index; current_index < end_index; current_index++)
        {
          auto l = mzs[current_index];
          int edge_score = scores[end_index];

          // make edge from r to l if they make an a.a. mass.
          std::vector<Residue> aas;
          std::vector<std::vector<Residue>> gaps;
          if (g == 0)
          {
            aas = getAA_(l, r, tol, n, mode);
            if (aas.empty()) continue;
          }
          else
          {
            gaps = getGap_(l, r, tol, 0);
            if (gaps.empty()) continue;
          }
          // end_index, current_index to amino acid strings.
          if (edge_aa_map_.find(end_index) == edge_aa_map_.end()) { edge_aa_map_[end_index] = std::map<int, std::vector<String>>(); }
          auto& e = edge_aa_map_[end_index];

          if (e.find(current_index) == e.end()) { e[current_index] = std::vector<String>(); }

          if (g == 0)
          {
            for (const auto& aa : aas)
            {
              auto aaStr = n == 0 ? aa.toString() : aa.toString().toLower();
              e[current_index].push_back(aaStr);
            }
          }
          else
          {
            for (const auto& gap : gaps)
            {
              std::stringstream ss;
              for (const auto& aa : gap)
                ss << aa.toString().toLower();
              e[current_index].emplace_back(ss.str());
            }
          }

          int gap_diff = g == 0 ? 0 : 1;

          for (int d = 0; d + gap_diff <= max_gap_count_; d++)
          {
            for (int iso = 0; iso + n <= (consider_diff_ion_jumps_ ? 1 : 0); iso++)
            {
              int t_edge_score = edge_score -  (iso == 0 ? 0 : 1); // if different ion types, subtract 1
              for (int lvl = 0; lvl < length; lvl++)
              {
                for (int score = min_path_score_; score <= max_path_score_; score++)
                {
                  if (score - t_edge_score < min_path_score_) continue;
                  if (score - t_edge_score > max_path_score_) break;

                  vertex1 = getVertex_(end_index, score, lvl + 1, iso + n, d + gap_diff);
                  vertex2 = getVertex_(current_index, score - t_edge_score, lvl, iso, d);
                  dag.addEdge(vertex1, vertex2, visited);
                }
              }
            }
          }
        }
      }
    }
    //  make edge from sink to r
    if (end_index < (int)mzs.size() - 1)
    {
      for (int d = 0; d <= max_gap_count_; d++)
      {
        for (int iso = 0; iso <= (consider_diff_ion_jumps_ ? 1 : 0); iso++)
        {
          for (int score = min_path_score_; score <= max_path_score_; score++)
          {
            vertex1 = getVertex_((int)mzs.size() - 1, score, length, iso, d);
            vertex2 = getVertex_(end_index, score, length, iso, d);
            dag.addEdge(vertex1, vertex2, visited);
          }
        }
      }
    }
    end_index++;
  }
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(): DefaultParamHandler("FLASHTaggerAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(const FLASHTaggerAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHTaggerAlgorithm& FLASHTaggerAlgorithm::operator=(const FLASHTaggerAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHTaggerAlgorithm::setDefaultParams_()
{
//  defaults_.setValue("max_count", 400,
//                     "Maximum number of tags per length, defined by the -min_length and -max_length options. Tags with different amino acid combinations but identical masses are counted once (e.g., TII, TIL, TLI, and TLL are distinct but counted as one).");
//  defaults_.setMinInt("max_count", 0);

  defaults_.setValue(
    "min_length", 3,
    "Minimum tag length, where each mass gap contributes one unit of length (even if represented by multiple amino acids). For instance, the length of TA[255]G is 4.");
  defaults_.setMaxInt("min_length", 9);
  defaults_.setMinInt("min_length", 3);

  defaults_.setValue(
    "max_length", 8,
    "Maximum tag length, where each mass gap contributes one unit of length (even if represented by multiple amino acids). For instance, the length of TA[255]G is 4.");
  defaults_.setMaxInt("max_length", 30);
  defaults_.setMinInt("max_length", 3);

//  defaults_.setValue("flanking_mass_tol", -1,
//                     "Flanking mass tolerance (the flanking mass minus protein mass up to the matching amino acid) in Da. This limits the possible "
//                     "terminal modification mass. Set to a positive value to activate.");
//  defaults_.addTag("flanking_mass_tol", "advanced");

  //defaults_.setValue("max_iso_error_count", 0, "Maximum isotope error count allowed per tag.");

  // defaults_.setValue("allow_iso_error", "false", "Allow up to one isotope error in each tag.");
  // defaults_.setValidStrings("allow_iso_error", {"true", "false"});
  // defaults_.addTag("allow_iso_error", "advanced");

//  defaults_.setValue("min_matched_aa", 3, "Minimum number of amino acids in matched proteins covered by tags.");
//  defaults_.addTag("min_matched_aa", "advanced");

  defaults_.setValue("allow_gap", "false", "Allows mass gaps (representing multiple consecutive amino acids) in tags.");
  defaults_.setValidStrings("allow_gap", {"true", "false"});
  defaults_.addTag("allow_gap", "advanced");

  defaults_.setValue("max_aa_in_gap", 2, "Maximum number of amino acids in a mass gap.");
  defaults_.setMaxInt("max_aa_in_gap", 3);
  defaults_.setMinInt("max_aa_in_gap", 2);
  defaults_.addTag("max_aa_in_gap", "advanced");

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaultsToParam_();
}

bool areSetsEqualWithinTolerance(const std::set<double>& set1, const std::set<double>& set2, double tolerance) {
  if (set1.size() != set2.size()) {
    return false;  // Early exit if sizes are different
  }

  auto it1 = set1.begin();
  auto it2 = set2.begin();

  while (it1 != set1.end() && it2 != set2.end()) {
    if (std::fabs(*it1 - *it2) > tolerance) {
      return false;  // Elements are not within tolerance
    }
    ++it1;
    ++it2;
  }

  return true;  // All elements matched within tolerance
}

void FLASHTaggerAlgorithm::updateMembers_()
{
  min_tag_length_ = param_.getValue("min_length");
  max_tag_length_ = param_.getValue("max_length");

  std::vector<double> prefix_shifts, suffix_shifts;
  common_shifts_.clear();
  n_term_shifts_.clear();
  c_term_shifts_.clear();
  for (const auto& ion_str : param_.getValue("ion_type").toStringVector())
  {
    if (ion_str == "a") { prefix_shifts.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b") { prefix_shifts.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c") { prefix_shifts.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x") { suffix_shifts.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y") { suffix_shifts.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z") { suffix_shifts.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else if (ion_str == "zp1") { suffix_shifts.push_back(Residue::getInternalToZp1Ion().getMonoWeight()); }
    else if (ion_str == "zp2") { suffix_shifts.push_back(Residue::getInternalToZp2Ion().getMonoWeight()); }
    else { continue; }
  }
  for (const auto& m1 : prefix_shifts)
  {
    for (const auto& m2 : prefix_shifts)
    {
      if (m1 < m2) continue;
      n_term_shifts_.insert(m1 - m2);
    }
  }

  for (const auto& m1 : suffix_shifts)
  {
    for (const auto& m2 : suffix_shifts)
    {
      if (m1 < m2) continue;
      c_term_shifts_.insert(m1 - m2);
    }
  }

  if (areSetsEqualWithinTolerance(n_term_shifts_, c_term_shifts_, .0001)) common_shifts_ = c_term_shifts_;

  consider_diff_ion_jumps_ = common_shifts_.size() > 1 || n_term_shifts_.size() > 1 || c_term_shifts_.size() > 1;
  //min_cov_aa_ = (int)param_.getValue("min_matched_aa");
  max_aa_in_gap_ = param_.getValue("max_aa_in_gap");
  max_gap_count_ = param_.getValue("allow_gap").toString() == "true" ? 1 : 0;
  // prsm_fdr_ = param_.getValue("fdr");
  //flanking_mass_tol_ = param_.getValue("flanking_mass_tol");
  updateEdgeMasses_();
  max_edge_mass_ = consider_diff_ion_jumps_ ? (EmpiricalFormula("H4O2").getMonoWeight()) : 0;
  max_edge_mass_ += gap_mass_map_.empty() ? aa_mass_map_.rbegin()->first : std::max(aa_mass_map_.rbegin()->first, gap_mass_map_.rbegin()->first);
}

const MSSpectrum& FLASHTaggerAlgorithm::getSpectrum() const
{
  return spec_;
}

void FLASHTaggerAlgorithm::run(const DeconvolvedSpectrum& deconvolved_spectrum, double ppm)
{
  // setLogType(CMD);
  if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy()) return;

  auto tags = std::vector<FLASHHelperClasses::Tag>();
  tags.reserve(max_tag_length_ * max_tag_counts_.back());

  spec_.setRT(deconvolved_spectrum.getOriginalSpectrum().getRT());
  if (deconvolved_spectrum.getPrecursorPeakGroup().getMonoMass() > 0)
  {
    spec_.setMetaValue("PrecursorMass", deconvolved_spectrum.getPrecursorPeakGroup().getMonoMass());
  }

  for (int mode = 0; mode < (common_shifts_.empty() ? 3 : 1); mode++)
    getTags_(deconvolved_spectrum, ppm, mode);

  std::sort(tags_.rbegin(), tags_.rend());

  int index = 0;
  for (auto& tag : tags_)
  {
    tag.setIndex(index++);
    tag.setRetentionTime(deconvolved_spectrum.getOriginalSpectrum().getRT());
  }
}

void FLASHTaggerAlgorithm::indexFasta_(const std::vector<FASTAFile::FASTAEntry>& fasta_entry)
{
  if (!indexed_fasta_.empty()) return;
  std::vector<String> keys = {""};
  keys.reserve( std::pow(aas_.size(), min_tag_length_) + 1);
  for (Size d = 0; d < min_tag_length_; d++)
  {
    auto t_keys = keys;
    for (const auto& aa : aas_)
    {
      const auto& aaStr = aa->toString();
      for (const auto& key : t_keys)
      {
        if (key.length() != d) continue;
        keys.push_back(key + aaStr);
      }
    }
  }

#pragma omp parallel for default(none) shared(keys, fasta_entry, std::cout)
  for (int i = 0 ; i < (int) keys.size(); i++)
  {
    const auto& key = keys[i];
    if (key.length() < min_tag_length_) continue;
    std::set<Size> value;
    for (Size j = 0 ; j < fasta_entry.size(); j++)
    {
      if (!fasta_entry[j].sequence.hasSubstring(key)) continue;
      value.insert(j);
    }
    //std::cout<<key<< " " << value.size()<<std::endl;
#pragma omp critical
    indexed_fasta_[key] = value;
  }
}

int FLASHTaggerAlgorithm::getPeakGroupScore(const PeakGroup& peak_group)
{
  // (int)round(5 * log10(std::max(1e-1, pg.getQscore() / std::max(1e-1, (1.0 - pg.getQscore())))));
  return std::max(1, (int)round(max_peak_group_score * peak_group.getQscore()));
}

void FLASHTaggerAlgorithm::getTags_(const DeconvolvedSpectrum& dspec, double ppm, int mode) // mode 0 : common 1 : n 2 : c term
{
  if (mode == 0 && common_shifts_.empty()) return;
  if (mode == 1 && n_term_shifts_.empty()) return;
  if (mode == 2 && c_term_shifts_.empty()) return;

  std::vector<double> mzs;
  std::vector<int> scores;
  mzs.reserve(dspec.size());
  scores.reserve(dspec.size());

  for (const auto& pg : dspec)
  {
    int score = getPeakGroupScore(pg);
    // if (score <= -5) continue;
    scores.push_back(score);
    mzs.push_back(pg.getMonoMass());
  }
  getTags_(mzs, scores, dspec.getScanNumber(), ppm, mode);
}

void FLASHTaggerAlgorithm::updateTagSet_(std::set<FLASHHelperClasses::Tag>& tag_set,
                                         std::map<String, std::vector<FLASHHelperClasses::Tag>>& seq_tag,
                                         const std::vector<Size>& path,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         int scan,
                                         double ppm,
                                         int mode)
{
  double flanking_mass = -1;
  std::vector<String> seqs {""};
  std::vector<double> tag_mzs;
  std::vector<int> tag_scores;
  tag_mzs.reserve(path.size() - 1);
  tag_scores.reserve(path.size() - 1);

  for (int j = 1; j < (int)path.size(); j++)
  {
    int i1 = getIndex_(path[j - 1]); // c term side
    int i2 = getIndex_(path[j]);     // n term side

    if (edge_aa_map_.find(i1) != edge_aa_map_.end() && edge_aa_map_[i1].find(i2) != edge_aa_map_[i1].end())
    {
      auto& edge_aa = edge_aa_map_[i1];
      std::vector<String> tmp_seqs;
      tmp_seqs.reserve(seqs.size() * edge_aa[i2].size());
      for (const auto& tmp_seq : seqs)
      {
        for (const auto& seq : edge_aa[i2])
        {
          tmp_seqs.emplace_back(seq + tmp_seq);
        }
      }
      seqs = tmp_seqs;
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
    }
    else if (i2 == 0) // nterm
    {
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
      flanking_mass = mzs[i1];
    }
  }

  std::vector<double> rev_tag_mzs;
  rev_tag_mzs.reserve(tag_mzs.size());
  std::vector<int> rev_tag_scores;
  rev_tag_scores.reserve(tag_scores.size());

  for (int i = (int)tag_mzs.size() - 1; i >= 0; i--)
  {
    rev_tag_mzs.push_back(tag_mzs[i]);
    rev_tag_scores.push_back(tag_scores[i]);
  }

  for (const auto& seq : seqs)
  {
    String rev_seq = String(seq).reverse();
    if (rev_seq == seq) continue; // exclude palindromic tags.
    auto iter = seq_tag.find(seq);
    bool pass = (mode != 2);
    if (pass && iter != seq_tag.end()) // remove overlapping tags.
    {
      for (const auto& pt : iter->second)
      {
        if (pt.getNtermMass() < 0) continue;
        if (abs(pt.getNtermMass() - flanking_mass) / std::max(pt.getNtermMass(), flanking_mass) * 1e6 > ppm) continue;
        pass = false;
        break;
      }
    }
    if (pass)
    {
      const auto direct_tag = FLASHHelperClasses::Tag(seq, flanking_mass, -1, tag_mzs, tag_scores, scan);
      tag_set.insert(direct_tag);
      seq_tag[seq].push_back(direct_tag);
    }

    pass = (mode != 1);

    iter = seq_tag.find(rev_seq);
    if (pass && iter != seq_tag.end()) // remove overlapping tags.
    {
      for (const auto& pt : iter->second)
      {
        if (pt.getCtermMass() < 0) continue;
        if (abs(pt.getCtermMass() - flanking_mass) / std::max(pt.getCtermMass(), flanking_mass) * 1e6 > ppm) continue;
        pass = false;
        break;
      }
    }
    if (pass)
    {
      const auto reverse_tag = FLASHHelperClasses::Tag(rev_seq, -1, flanking_mass, rev_tag_mzs, rev_tag_scores, scan);
      tag_set.insert(reverse_tag);
      seq_tag[rev_seq].push_back(reverse_tag);
    }
  }
}

void FLASHTaggerAlgorithm::getScoreAndMatchCount_(const std::vector<int>& spec_vec,
                                                  const boost::dynamic_bitset<>& pro_vec,
                                                  // const boost::dynamic_bitset<>& mask_pro_vec,
                                                  const std::set<int>& spec_pro_diffs,
                                                  const std::vector<int>& spec_scores,
                                                  int& max_score,
                                                  int& match_cntr) const
{
  max_score = 0;
  match_cntr = 0;
  for (int d : spec_pro_diffs)
  {
    int score = 0;
    int cntr = 0;
    for (int i = 0; i < spec_vec.size(); i++)
    {
      int index = d + spec_vec[i];
      if (index < 0) continue;
      if (index >= pro_vec.size()) break;
      if (! pro_vec[index]) continue;
      cntr++;
      score += spec_scores[i];
    }
    if (cntr >= match_cntr)
    {
      max_score = std::max(max_score, score);
      match_cntr = cntr;
    }
  }
}


void FLASHTaggerAlgorithm::getTags_(const std::vector<double>& mzs, const std::vector<int>& scores, int scan, double ppm, int mode)
{
  std::vector<int> _scores;
  std::vector<double> _mzs;
  int threshold;

  if (mzs.size() >= max_node_cntr)
  {
    _scores = scores;
    std::sort(_scores.rbegin(), _scores.rend());
    threshold = _scores[max_node_cntr - 1];
    _scores.clear();

    _mzs.reserve(max_node_cntr + 1);
    _scores.reserve(max_node_cntr + 1);
  }
  else
  {
    _mzs.reserve(mzs.size() + 1);
    _scores.reserve(mzs.size() + 1);
    threshold = *std::min_element(scores.begin(), scores.end());
  }

  _mzs.push_back(.0);
  _scores.push_back(0);
  spec_.reserve(_mzs.size());
  for (int i = 0; i < (int)mzs.size(); i++)
  {
    if (scores[i] < threshold) continue;
    _mzs.push_back(mzs[i]);
    _scores.push_back(scores[i]);
    spec_.emplace_back(mzs[i], scores[i]);
  }

  // filtration of top 500 masses is done

  int max_vertex_score = *std::max_element(_scores.begin(), _scores.end());
  int min_vertex_score = *std::min_element(_scores.begin(), _scores.end());

  max_path_score_ = std::max(max_vertex_score, max_vertex_score) * (max_tag_length_ + 2);
  min_path_score_ = std::max(min_vertex_score, min_vertex_score) * (max_tag_length_ + 2);

  max_path_score_ = std::max(max_path_score_, std::max(max_vertex_score, max_vertex_score) * (min_tag_length_ - 2));
  min_path_score_ = std::min(min_path_score_, std::max(min_vertex_score, min_vertex_score) * (min_tag_length_ - 2));
  min_path_score_ = std::max(0, min_path_score_);

  std::set<FLASHHelperClasses::Tag> tagSet;
  std::map<String, std::vector<FLASHHelperClasses::Tag>> seq_tag;

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    auto tag_count = max_tag_counts_[std::min((int)max_tag_counts_.size() - 1, length)];
    FLASHHelperClasses::DAG dag(_mzs.size() * (1 + max_tag_length_) * (1 + max_gap_count_) * (1 + (consider_diff_ion_jumps_ ? 1 : 0))
                                * (1 + max_path_score_ - min_path_score_));
    constructDAG_(dag, _mzs, _scores, length, ppm, mode);

    std::set<FLASHHelperClasses::Tag> _tagSet;
    for (int score = max_path_score_; score >= min_path_score_ && (int)_tagSet.size() < tag_count; score--)
    {
      std::vector<std::vector<Size>> all_paths;
      all_paths.reserve(tag_count);
      for (int d = 0; d <= max_gap_count_; d++)
      {
        for (int diff_ion_index = 0; diff_ion_index <= (consider_diff_ion_jumps_ ? 1 : 0); diff_ion_index++)
        {
          dag.findAllPaths(getVertex_((int)_mzs.size() - 1, score, length, diff_ion_index, d), getVertex_(0, 0, 0, 0, 0), all_paths, tag_count);
        }
      }

      for (const auto& path : all_paths)
      {
        updateTagSet_(_tagSet, seq_tag, path, _mzs, _scores, scan, ppm, mode);
      }
    }
    tagSet.insert(_tagSet.begin(), _tagSet.end());
  }

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    // int count = 0;
    for (const auto& tag : tagSet)
    {
      if ((int)tag.getLength() != length) continue;
      tags_.push_back(tag);
    }
  }

  std::sort(tags_.begin(), tags_.end(), [](const FLASHHelperClasses::Tag& a, const FLASHHelperClasses::Tag& b) {
    return a.getLength() == b.getLength()
             ? (a.getScore() == b.getScore()
                  ? (a.getUppercaseSequence() == b.getUppercaseSequence()
                       ? (a.getNtermMass() == b.getNtermMass() ? a.getCtermMass() > b.getCtermMass() : a.getNtermMass() > b.getNtermMass())
                       : a.getUppercaseSequence() > b.getUppercaseSequence())
                  : a.getScore() > b.getScore())
             : a.getLength() > b.getLength();
  });
}

Size FLASHTaggerAlgorithm::find_with_X_(const std::string_view& A, const String& B, Size pos) // allow a single X. pos is in A
{
  if (A.length() <= B.length()) return String::npos;
  for (size_t i = pos; i <= A.length() - B.length(); ++i)
  {
    bool match = true;
    int x_cntr = 0;
    for (size_t j = 0; j < B.length(); ++j)
    {
      if (A[i + j] == 'X') x_cntr++;
      if ((A[i + j] != B[j] && A[i + j] != 'X') || x_cntr > 1)
      {
        match = false;
        break;
      }
    }
    if (match) { return i; }
  }
  return String::npos;
}

// Make output struct containing all information about matched entries and tags, coverage, score etc.
void FLASHTaggerAlgorithm::runMatching(const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                                       const DeconvolvedSpectrum& deconvolved_spectrum,
                                       const std::vector<boost::dynamic_bitset<>>& vectorized_fasta_entry,
                                       const std::vector<boost::dynamic_bitset<>>& reversed_vectorized_fasta_entry,
                                       const std::vector<std::map<int, double>>& mass_map,
                                       const std::vector<std::map<int, double>>& rev_mass_map,
                                       const double max_mod_mass)
{
  //indexFasta_(fasta_entry);
  std::vector<std::pair<ProteinHit, std::vector<int>>> pairs;
  //std::set<Size> matched_fasta_entries;
  int scan = deconvolved_spectrum.getScanNumber();
  protein_hits_.clear();

  // int min_hit_tag_score = max_path_score_;
  startProgress(0, (SignedSize)fasta_entry.size(), "Running FLASHTagger: searching database");

  if (decoy_factor_ < 0)
  {
    double taget_count = 0;
    double decoy_count = 0;
    for (const auto& fe : fasta_entry)
    {
      if (fe.identifier.hasPrefix("DECOY")) { decoy_count++; }
      else { taget_count++; }
      decoy_factor_ = decoy_count / taget_count;
    }
  }

  std::vector<int> spec_vec;
  std::vector<int> spec_scores;

  spec_vec.reserve(deconvolved_spectrum.size() + 1);
  spec_scores.reserve(deconvolved_spectrum.size() + 1);

  spec_vec.push_back(0);
  spec_scores.push_back(1);
  for (const auto& pg : deconvolved_spectrum)
  {
    Size mn = round(ConvolutionBasedProteinFilter::multi_factor_for_vectorization * pg.getMonoMass());
    spec_vec.push_back(mn);
    spec_scores.push_back(FLASHTaggerAlgorithm::getPeakGroupScore(pg));
  }

#pragma omp parallel for default(none)                                                                                                     \
  shared(mass_map, rev_mass_map, spec_vec, spec_scores, deconvolved_spectrum, pairs, fasta_entry, \
           scan, max_mod_mass, vectorized_fasta_entry, reversed_vectorized_fasta_entry, std::cout)
  for (int i = 0; i < fasta_entry.size(); i++)
  {
    const auto& fe = fasta_entry[i];
    nextProgress();
    std::vector<int> matched_tag_indices;
    String seq = fe.sequence;
    seq.erase(remove(seq.begin(), seq.end(), 'X'), seq.end()); // remove all X
    // auto x_pos = fe.sequence.find('X');
    const auto& pro_vec = vectorized_fasta_entry[i];
    const auto& rev_pro_vec = reversed_vectorized_fasta_entry[i];
    const auto& pro_masses = mass_map[i];
    const auto& rev_pro_masses = rev_mass_map[i];

    // find range, match allowing X.
    std::set<double> matched_masses;
    std::set<int> n_spec_pro_diffs, c_spec_pro_diffs;

    for (int j = 0; j < (int)tags_.size(); j++)
    {
      const auto& tag = tags_[j];
      const auto& tag_mzs = tag.getMzs();
      bool isSubset = std::includes(matched_masses.begin(), matched_masses.end(), tag_mzs.begin(), tag_mzs.end());
      if (isSubset) { continue; }

      bool isNterm = tag.getNtermMass() > 0;

      const auto& uppercase_tag_seq = tag.getUppercaseSequence();
      std::vector<int> positions;
      Size tpos = 0;
      while (true)
      {
        tpos = seq.find(uppercase_tag_seq, tpos);
        if (tpos == std::string_view::npos) break;
        positions.push_back((int)tpos);
        tpos++;
      }

      bool matched = false;
      double n_term_mass, c_term_mass;

      for (const auto& pos : positions)
      {
        if (isNterm && pos >= 0)
        {
          n_term_mass = pro_masses.at(pos);
          double flanking_mass = n_term_mass - tag.getNtermMass();
          if (max_mod_mass > 0 && flanking_mass < -max_mod_mass) continue;
          n_spec_pro_diffs.insert(round(ConvolutionBasedProteinFilter::multi_factor_for_vectorization * flanking_mass));
        }
        else if (! isNterm && pos + tag.getSequence().length() < seq.length())
        {
          c_term_mass = rev_pro_masses.at(pos + tag.getSequence().length());
          double flanking_mass = c_term_mass - tag.getCtermMass();
          if (max_mod_mass > 0 && flanking_mass < -max_mod_mass) continue;
          c_spec_pro_diffs.insert(round(ConvolutionBasedProteinFilter::multi_factor_for_vectorization * flanking_mass));
        }
        else
          continue;

        matched = true;
      }
      if (matched)
      {
        matched_masses.insert(tag.getMzs().begin(), tag.getMzs().end());
        matched_tag_indices.push_back(j); // tag indices
      }
    }
    if (matched_tag_indices.empty()) continue;
    int match_cntr;
    int match_score1 = 0, match_score2 = 0;
    int match_cntr1 = 0, match_cntr2 = 0;

    if (! n_spec_pro_diffs.empty()) { getScoreAndMatchCount_(spec_vec, pro_vec, n_spec_pro_diffs, spec_scores, match_score1, match_cntr1); }
    if (! c_spec_pro_diffs.empty()) { getScoreAndMatchCount_(spec_vec, rev_pro_vec, c_spec_pro_diffs, spec_scores, match_score2, match_cntr2); }

    match_cntr = match_cntr1 + match_cntr2;
    if (match_cntr < min_cov_aa_) continue;
    ProteinHit hit(std::max(match_score1, match_score2), 0, fe.identifier, fe.sequence); //
    hit.setDescription(fe.description);
    hit.setMetaValue("Scan", scan);
    hit.setMetaValue("MatchedAA", match_cntr);
    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length());
#pragma omp critical
    {
      pairs.emplace_back(hit, matched_tag_indices);
    }
  }

  endProgress();
  protein_hits_.reserve(pairs.size());

  for (auto& [hit, indices] : pairs)
  {
    hit.setMetaValue("TagIndices", indices);
    protein_hits_.push_back(hit);
  }
  std::sort(protein_hits_.begin(), protein_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore()
             ? (left.getCoverage() == right.getCoverage() ? (left.getDescription() > right.getDescription()) : (left.getCoverage() > right.getCoverage()))
             : (left.getScore() > right.getScore());
  });
}

void FLASHTaggerAlgorithm::getProteinHits(std::vector<ProteinHit>& hits, int max_target_out) const
{
  hits.reserve(protein_hits_.size());
  int count = 0;

  count = 0;
  double prev_score = 0;
  for (const auto& hit : protein_hits_)
  {

    if (hit.getScore() <= 0) break;
    hits.push_back(hit);
    //if (hit.getAccession().hasPrefix("DECOY")) continue;
    if (max_target_out > 0 && count++ >= max_target_out)
    {
      if (hit.getScore() < prev_score) break; // keep adding if the score is the same
    }
    prev_score = hit.getScore();
  }
}

void FLASHTaggerAlgorithm::getTags(std::vector<FLASHHelperClasses::Tag>& tags, int tag_length) const
{
  for (const auto& tag : tags_)
  {
    if (tag_length > 0 && tag.getLength() != tag_length) continue;
    tags.push_back(tag);
  }
}

void FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(std::vector<int>& positions,
                                                                   std::vector<double>& masses,
                                                                   double flanking_mass_tol,
                                                                   const String& seq, //
                                                                   const FLASHHelperClasses::Tag& tag)
{
  std::vector<int> indices;
  const auto& tagseq = tag.getUppercaseSequence();
  Size pos = find_with_X_(seq, tagseq);
  while (pos != String::npos)
  {
    double delta_mass = - flanking_mass_tol - 1;

    auto x_pos = seq.find('X');
    if (tag.getNtermMass() > 0)
    {
      auto nterm = seq.substr(0, pos);
      if (! nterm.empty())
      {
        if (x_pos != String::npos) { nterm.erase(remove(nterm.begin(), nterm.end(), 'X'), nterm.end()); }
        delta_mass = tag.getNtermMass() - (nterm.empty() ? 0 : AASequence::fromString(nterm).getMonoWeight(Residue::Internal));
      }
    }
    else
    {
      auto cterm = seq.substr(pos + tag.getSequence().length());
      if (! cterm.empty())
      {
        if (x_pos != String::npos) cterm.erase(remove(cterm.begin(), cterm.end(), 'X'), cterm.end());
        delta_mass = tag.getCtermMass() - (cterm.empty() ? 0 : AASequence::fromString(cterm).getMonoWeight(Residue::Internal));
      }
    }

    if (delta_mass != - flanking_mass_tol - 1 && (flanking_mass_tol < 0 || std::abs(delta_mass) <= flanking_mass_tol))
    {
      masses.push_back(delta_mass);
      positions.push_back((int)pos);
    }
    pos = find_with_X_(seq, tagseq, pos + 1);
  }
}
} // namespace OpenMS
