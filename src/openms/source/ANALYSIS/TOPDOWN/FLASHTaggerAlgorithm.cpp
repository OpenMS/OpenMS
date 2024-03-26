// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <utility>

namespace OpenMS
{
inline const Size max_node_cntr = 500;
class FLASHTaggerAlgorithm::DAC_
{
private:
  int vertex_count_;
  // 0, 1, 2, ... ,vertex_count - 1
  std::vector<std::vector<int>> adj_list_; //

public:
  DAC_(int vertice_count): vertex_count_(vertice_count), adj_list_(vertice_count)
  {
  }

  int size()
  {
    return vertex_count_;
  }

  void addEdge(int src, int dest)
  {
    adj_list_[src].push_back(dest); //
  }

  void findAllPaths(int source, int sink, std::vector<std::vector<int>>& all_paths, int max_count)
  {
    boost::dynamic_bitset<> visited(vertex_count_);
    std::vector<int> path;

    findAllPaths_(source, sink, visited, path, all_paths, max_count); // reverse traveling
  }

private:
  void findAllPaths_(int current,
                     int destination,
                     boost::dynamic_bitset<>& visited,
                     std::vector<int>& path,
                     std::vector<std::vector<int>>& all_paths,
                     int max_count)
  {
    if ((int)all_paths.size() >= max_count) return;
    visited[current] = true;
    path.push_back(current);

    if (current == destination)
    {
      // add the current path
      all_paths.push_back(path);
    }
    else
    {
      // Recursively explore neighbors
      for (const auto& neighbor : adj_list_[current])
      {
        if (! visited[neighbor]) { findAllPaths_(neighbor, destination, visited, path, all_paths, max_count); }
      }
    }

    // Backtrack
    visited[current] = false;
    path.pop_back();
  }
};

std::vector<Residue> FLASHTaggerAlgorithm::getAA_(double l, double r, double tol, int iso_offset) const
{
  std::vector<Residue> ret;
  if (l == r) return ret;
  double iso_mass = std::abs(iso_offset * Constants::C13C12_MASSDIFF_U);
  double diff1 = std::abs(std::abs(r - l) - iso_mass);
  double diff2 = std::abs(std::abs(r - l) + iso_mass);
  double abs_tol = std::max(l, r) * tol / 1e6 * 2;
  auto iter = aa_mass_map_.lower_bound(diff1 - abs_tol);

  while (iter != aa_mass_map_.end())
  {
    if (std::abs(diff1 - iter->first) < abs_tol || std::abs(diff2 - iter->first) < abs_tol)
    {
      for (auto& aa : iter->second)
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
  // gap_mass_map_.clear();

  for (const auto& aa : aas_)
  {
    double aa_mass = aa->getMonoWeight(Residue::Internal);
    if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end()) aa_mass_map_[aa_mass] = std::vector<Residue>();
    aa_mass_map_[aa_mass].push_back(*aa);
  }
}

int FLASHTaggerAlgorithm::getVertex_(int index, int path_score, int level, int iso_level) const
{
  return ((index * (max_tag_length_ + 1) + level) * (max_iso_in_tag_ + 1) + iso_level) * (max_path_score_ - min_path_score_ + 1)
         + (path_score - min_path_score_);
}

int FLASHTaggerAlgorithm::getIndex_(int vertex) const
{
  return ((vertex / (max_path_score_ - min_path_score_ + 1)) / (max_iso_in_tag_ + 1)) / (max_tag_length_ + 1);
}

bool FLASHTaggerAlgorithm::connectEdge_(FLASHTaggerAlgorithm::DAC_& dac, int vertex1, int vertex2, boost::dynamic_bitset<>& visited)
{
  if (vertex1 < 0 || vertex2 < 0 || vertex1 >= (int)visited.size() || vertex2 >= (int)visited.size()) return false;
  if (! visited[vertex2]) return false;

  dac.addEdge(vertex1, vertex2);
  return visited[vertex1] = true;
}


void FLASHTaggerAlgorithm::constructDAC_(FLASHTaggerAlgorithm::DAC_& dac,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         int length,
                                         double tol)
{
  // from source to sink, connect but the edge direction is from sink to source.
  edge_aa_map_.clear();
  int start_index = 1; // zeroth = source.
  int end_index = 1;
  boost::dynamic_bitset<> visited(dac.size());
  visited[getVertex_(0, 0, 0, 0)] = true;

  while (end_index < (int)mzs.size())
  {
    auto r = mzs[end_index];

    // first, make edge from r to source and sink to r.
    int vertex1 = getVertex_(end_index, scores[end_index], 0, 0);
    int vertex2 = getVertex_(0, 0, 0, 0);

    connectEdge_(dac, vertex1, vertex2, visited);

    // from an edge i, j to class edge.  for each i, j make a unique key. key to an edge.

    while (start_index < end_index && r - mzs[start_index] > max_edge_mass_)
      start_index++;

    for (int n = 0; n < 2; n++) // 0 for all a.a 1 for isotope errors. Allow only one isotope errors.
    {
      for (int current_index = start_index; current_index < end_index; current_index++)
      {
        auto l = mzs[current_index];
        int edge_score = scores[end_index];

        // make edge from r to l if they make an a.a. mass.
        std::vector<Residue> aas;

        aas = getAA_(l, r, tol, n);
        if (aas.empty()) continue;


        // end_index, current_index to amino acid strings.
        if (edge_aa_map_.find(end_index) == edge_aa_map_.end()) { edge_aa_map_[end_index] = std::map<int, std::vector<String>>(); }
        auto& e = edge_aa_map_[end_index];

        if (e.find(current_index) == e.end()) { e[current_index] = std::vector<String>(); }

        for (auto& aa : aas)
        {
          auto aaStr = n == 0 ? aa.toString() : aa.toString().toLower();
          e[current_index].push_back(aaStr);
        }

        for (int g = 0; g + n <= max_iso_in_tag_; g++)
        {
          for (int lvl = 0; lvl < length; lvl++)
          {
            for (int score = min_path_score_; score <= max_path_score_; score++)
            {
              if (score - edge_score < min_path_score_) continue;
              if (score - edge_score > max_path_score_) break;

              int vertex1 = getVertex_(end_index, score, lvl + 1, g + n);
              int vertex2 = getVertex_(current_index, score - edge_score, lvl, g);
              connectEdge_(dac, vertex1, vertex2, visited);
            }
          }
        }
      }
      if (max_iso_in_tag_ == 0) break;
    }

    if (end_index < (int)mzs.size() - 1)
    {
      for (int g = 0; g <= max_iso_in_tag_; g++)
      {
        for (int score = min_path_score_; score <= max_path_score_; score++)
        {
          int edge_score = scores[mzs.size() - 1];
          if (score - edge_score < min_path_score_) continue;
          if (score - edge_score > max_path_score_) break;

          int vertex1 = (int)getVertex_((int)mzs.size() - 1, score, length, g);
          int vertex2 = (int)getVertex_(end_index, score - edge_score, length, g);
          connectEdge_(dac, vertex1, vertex2, visited);
        }
      }
    }
    end_index++;
  }
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(): DefaultParamHandler("FLASHTaggerAlgorithm")
{
  setDefaultParams_();
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(const FLASHTaggerAlgorithm& other): DefaultParamHandler(other)
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
  defaults_.setValue("max_tag_count", 500,
                     "Maximum number of the tags per length (lengths set by -min_length and -max_length options). The tags with different amino acid "
                     "combinations are all treated separately. E.g., "
                     "TII, TIL, TLI, TLL are distinct tags even though they have the same mass differences. "
                     "but are counted as four different tags. ");
  defaults_.setMinInt("max_tag_count", 0);

  defaults_.setValue(
    "min_length", 4,
    "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("min_length", 30);
  defaults_.setMinInt("min_length", 3);

  defaults_.setValue(
    "max_length", 10,
    "Maximum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("max_length", 30);
  defaults_.setMinInt("max_length", 3);

  defaults_.setValue("flanking_mass_tol", 500.0, "Flanking mass tolerance in Da.");
  defaults_.setValue("max_iso_error_count", 0, "Maximum isotope error count per tag.");
  defaults_.setMaxInt("max_iso_error_count", 2);
  defaults_.setMinInt("max_iso_error_count", 0);
  defaults_.addTag("max_iso_error_count", "advanced");
  defaults_.setValue("min_matched_aa", 5, "Minimum number of amino acids in matched proteins, covered by tags.");

  defaults_.setValue("fdr", 1.0, "Protein FDR threshold.");
  defaults_.setMaxFloat("fdr", 1.0);
  defaults_.setMinFloat("fdr", 0.01);
  defaults_.setValue("keep_decoy", "false", "Keep decoy proteins.");
  defaults_.addTag("keep_decoy", "advanced");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaultsToParam_();
}

void FLASHTaggerAlgorithm::updateMembers_()
{
  max_tag_count_ = param_.getValue("max_tag_count");
  min_tag_length_ = param_.getValue("min_length");
  max_tag_length_ = param_.getValue("max_length");
  max_iso_in_tag_ = param_.getValue("max_iso_error_count");
  min_cov_aa_ = (int)param_.getValue("min_matched_aa");
  fdr_ = param_.getValue("fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  updateEdgeMasses_();
  max_edge_mass_ = aa_mass_map_.rbegin()->first + max_iso_in_tag_ * Constants::C13C12_MASSDIFF_U;
}

void FLASHTaggerAlgorithm::run(const std::vector<DeconvolvedSpectrum>& deconvolved_spectra, double ppm)
{
  DeconvolvedSpectrum dspec_for_tagging;
  for (const auto& dspec : deconvolved_spectra)
  {
    if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() == 1) continue;
    for (const auto& pg : dspec)
    {
      dspec_for_tagging.push_back(pg);
    }
  }
  if (deconvolved_spectra.size() > 1)
  {
    dspec_for_tagging.sort();
    SpectralDeconvolution::removeOverlappingPeakGroups(dspec_for_tagging, ppm * 1e-6); // merged peak groups have scan number information!
  }
  run(dspec_for_tagging, ppm);
}

void FLASHTaggerAlgorithm::run(const DeconvolvedSpectrum& dspec, double ppm)
{
  std::vector<double> mzs;
  std::vector<int> scores;
  std::vector<int> scans;
  mzs.reserve(dspec.size());
  scores.reserve(dspec.size());
  std::vector<double> qscores;
  qscores.reserve(dspec.size());

  for (auto& pg : dspec)
  {
    qscores.push_back(pg.getQscore2D());
  }
  std::sort(qscores.rbegin(), qscores.rend());
  auto end = std::min(qscores.end(), qscores.begin() + max_node_cntr);
  double random_hit_prob = std::accumulate(qscores.begin(), end, .0);
  random_hit_prob /= (double)std::distance(qscores.begin(), end);

  for (auto& pg : dspec)
  {
    mzs.push_back(pg.getMonoMass());
    int score = (int)round(10 * log10(std::max(1e-6, pg.getQscore2D() / std::max(1e-6, (1.0 - random_hit_prob)))));

    scores.push_back(score); //
    scans.push_back(pg.getScanNumber());
  }
  run(mzs, scores, scans, ppm);
}

void FLASHTaggerAlgorithm::updateTagSet_(std::set<FLASHDeconvHelperStructs::Tag>& tag_set,
                                         std::map<String, std::vector<FLASHDeconvHelperStructs::Tag>>& seq_tag,
                                         const std::vector<int>& path,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         const std::vector<int>& scans,
                                         double ppm)
{
  double flanking_mass = -1;

  std::vector<String> seqs {""};
  std::vector<double> tag_mzs;
  std::vector<int> tag_scores;
  std::vector<int> tag_scans;
  tag_mzs.reserve(path.size() - 1);
  tag_scores.reserve(path.size() - 1);
  tag_scans.reserve(path.size() - 1);

  for (int j = 1; j < (int)path.size(); j++)
  {
    int i1 = getIndex_(path[j - 1]); // c term size
    int i2 = getIndex_(path[j]);     // n term side

    if (edge_aa_map_.find(i1) != edge_aa_map_.end() && edge_aa_map_[i1].find(i2) != edge_aa_map_[i1].end())
    {
      std::vector<String> tmp_seqs;
      tmp_seqs.reserve(seqs.size());
      for (const auto& tmp_seq : seqs)
      {
        for (const auto& seq : edge_aa_map_[i1][i2])
        {
          tmp_seqs.push_back(seq + tmp_seq);
        }
      }
      seqs = tmp_seqs;
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
      tag_scans.push_back(scans[i1]);
    }
    else if (i2 == 0) // nterm
    {
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
      tag_scans.push_back(scans[i1]);
      flanking_mass = mzs[i1];
    }
  }

  std::vector<double> rev_tag_mzs(tag_mzs);
  std::reverse(rev_tag_mzs.begin(), rev_tag_mzs.end());

  std::vector<int> rev_tag_scores(tag_scores);
  std::reverse(rev_tag_scores.begin(), rev_tag_scores.end());

  std::vector<int> rev_tag_scans(tag_scans);
  std::reverse(rev_tag_scans.begin(), rev_tag_scans.end());

  for (const auto& seq : seqs)
  {
    auto iter = seq_tag.find(seq);
    bool pass = true;
    if (iter != seq_tag.end()) // remove overlapping tags.
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
      auto direct_tag = FLASHDeconvHelperStructs::Tag(seq, flanking_mass, -1, tag_mzs, tag_scores, tag_scans);
      tag_set.insert(direct_tag);
      seq_tag[seq].push_back(direct_tag);
    }

    pass = true;
    String rev_seq = String(seq).reverse();
    iter = seq_tag.find(rev_seq);
    if (iter != seq_tag.end()) // remove overlapping tags.
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
      auto reverse_tag = FLASHDeconvHelperStructs::Tag(rev_seq, -1, flanking_mass, rev_tag_mzs, rev_tag_scores, rev_tag_scans);
      tag_set.insert(reverse_tag);
      seq_tag[rev_seq].push_back(reverse_tag);
    }
  }
}

void FLASHTaggerAlgorithm::run(const std::vector<double>& mzs, const std::vector<int>& scores, const std::vector<int>& scans, double ppm)
{
  if (max_tag_count_ == 0) return;

  std::vector<double> _mzs;
  std::vector<int> _scores, _scans;
  int threshold;

  if (mzs.size() >= max_node_cntr)
  {
    _scores = scores;
    std::sort(_scores.rbegin(), _scores.rend());
    threshold = _scores[max_node_cntr - 1];
    _scores.clear();

    _mzs.reserve(max_node_cntr + 1);
    _scores.reserve(max_node_cntr + 1);
    _scans.reserve(max_node_cntr + 1);
  }
  else
  {
    _mzs.reserve(mzs.size() + 1);
    _scores.reserve(mzs.size() + 1);
    _scans.reserve(mzs.size() + 1);
    threshold = *std::min_element(scores.begin(), scores.end());
  }

  _mzs.push_back(.0);
  _scores.push_back(0);
  _scans.push_back(0);
  for (int i = 0; i < (int)mzs.size(); i++)
  {
    if (scores[i] < threshold) continue;
    _mzs.push_back(mzs[i]);
    _scores.push_back(scores[i]);
    _scans.push_back(scans[i]);
  }
  // filtration of top 500 masses is done

  int max_vertex_score = *std::max_element(_scores.begin(), _scores.end());
  int min_vertex_score = *std::min_element(_scores.begin(), _scores.end());

  max_path_score_ = std::max(max_vertex_score, max_vertex_score) * (max_tag_length_ + 2);
  min_path_score_ = std::max(min_vertex_score, min_vertex_score) * (max_tag_length_ + 2);

  max_path_score_ = std::max(max_path_score_, std::max(max_vertex_score, max_vertex_score) * (min_tag_length_ - 2));
  min_path_score_ = std::min(min_path_score_, std::max(min_vertex_score, min_vertex_score) * (min_tag_length_ - 2));

  std::set<FLASHDeconvHelperStructs::Tag> tagSet;
  std::map<String, std::vector<FLASHDeconvHelperStructs::Tag>> seq_tag;

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    FLASHTaggerAlgorithm::DAC_ dac(_mzs.size() * (1 + max_tag_length_) * (1 + max_iso_in_tag_) * (1 + max_path_score_ - min_path_score_));
    constructDAC_(dac, _mzs, _scores, length, ppm);

    std::set<FLASHDeconvHelperStructs::Tag> _tagSet;
    for (int score = max_path_score_; score >= min_path_score_ && (int)_tagSet.size() < max_tag_count_; score--)
    {
      std::vector<std::vector<int>> all_paths;
      all_paths.reserve(max_tag_count_);
      for (int g = 0; g <= max_iso_in_tag_; g++)
      {
        dac.findAllPaths(getVertex_(_mzs.size() - 1, score, length, g), getVertex_(0, 0, 0, 0), all_paths, max_tag_count_);
      }
      for (const auto& path : all_paths)
      {
        updateTagSet_(_tagSet, seq_tag, path, _mzs, _scores, _scans, ppm);
      }
    }
    tagSet.insert(_tagSet.begin(), _tagSet.end());
  }

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    int count = 0;
    for (const auto& tag : tagSet)
    {
      if ((int)tag.getLength() != length) continue;
      tags_.push_back(tag);
      if (++count == max_tag_count_) break;
    }
    OPENMS_LOG_INFO << "Tag count with length " << length << ": " << count << std::endl;
  }

  std::sort(tags_.begin(), tags_.end(),
            [](const FLASHDeconvHelperStructs::Tag& a, const FLASHDeconvHelperStructs::Tag& b) { return a.getScore() > b.getScore(); });
}

Size FLASHTaggerAlgorithm::find_with_X_(const std::string_view& A, const String& B, Size pos) // allow a single X. pos is in A
{
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
void FLASHTaggerAlgorithm::runMatching(const String& fasta_file)
{
  if (tags_.empty()) return;
  std::vector<FASTAFile::FASTAEntry> fasta_entry;
  FASTAFile ffile;
  ffile.load(fasta_file, fasta_entry);

  std::vector<std::pair<ProteinHit, std::vector<int>>> pairs;
  std::vector<int> start_loc(tags_.size(), 0);
  std::vector<int> end_loc(tags_.size(), 0);

  // for each tag, find the possible start and end locations in the protein sequence. If C term, they are negative values to specify values are from
  // the end of the protein
#pragma omp parallel for default(none) shared(end_loc, start_loc)
  for (int i = 0; i < (int)tags_.size(); i++)
  {
    const auto& tag = tags_[i];
    auto flanking_mass = std::max(tag.getNtermMass(), tag.getCtermMass());

    start_loc[i] = std::max(0, int(floor(flanking_mass - flanking_mass_tol_) / aa_mass_map_.rbegin()->first));
    end_loc[i] = int(ceil(flanking_mass + flanking_mass_tol_) / aa_mass_map_.begin()->first) + (int)tag.getLength() + 1;
  }

  int min_hit_tag_score = max_path_score_;
  double decoy_mul = .0;

  for (int n = 0; n < 2; n++)
  {
#pragma omp parallel for default(none) shared(pairs, fasta_entry, start_loc, end_loc, decoy_mul, min_hit_tag_score, n)
    for (int i = 0; i < (int)fasta_entry.size(); i++)
    {
      const auto& fe = fasta_entry[i];
      bool is_decoy = false;
      if (fe.identifier.hasPrefix("DECOY")) { is_decoy = true; }

      if (is_decoy && n == 0)
      {
#pragma omp critical
        decoy_mul++;
        continue;
      }
      if (! is_decoy && n != 0) continue;

      std::vector<int> matched_tag_indices;
      auto x_pos = fe.sequence.find('X');
      std::map<Size, int> matched_pos_score;
      // find range, match allowing X.
      for (int j = 0; j < (int)tags_.size(); j++)
      {
        auto& tag = tags_[j];
        if (is_decoy && tag.getScore() < min_hit_tag_score) break;
        bool isNterm = tag.getNtermMass() > 0;

        int s, n;
        if (isNterm) { s = start_loc[j]; }
        else { s = std::max(0, int(fe.sequence.length()) - 1 - end_loc[j]); }
        n = std::min(end_loc[j] - start_loc[j], int(fe.sequence.length()) - s);
        if (n < (int)tag.getLength()) continue;
        const auto sub_seq = std::string_view(fe.sequence.data() + s, n);

        auto uppercase_tag_seq = tag.getSequence().toUpper();
        std::vector<int> positions;
        Size tpos = 0;
        while (true)
        {
          tpos = sub_seq.find(uppercase_tag_seq, tpos);
          if (tpos == std::string_view::npos) break;
          positions.push_back((int)tpos + s);
          tpos++;
        }

        if (positions.empty() && (int)x_pos >= s && (int)x_pos <= s + n) // only if perfect hits are not found and X exists
        {
          tpos = 0;
          while (true)
          {
            tpos = find_with_X_(sub_seq, uppercase_tag_seq, tpos);
            if (tpos == std::string_view::npos) break;
            positions.push_back((int)tpos + s);
            tpos++;
          }
        }

        bool matched = false;
        for (const auto& pos : positions)
        {
          if (tag.getNtermMass() > 0 && pos >= 0)
          {
            auto nterm = fe.sequence.substr(0, pos);
            if (x_pos != String::npos) { nterm.erase(remove(nterm.begin(), nterm.end(), 'X'), nterm.end()); }
            double aamass = nterm.empty() ? 0 : AASequence::fromString(nterm).getMonoWeight();
            if (std::abs(tag.getNtermMass() - aamass) > flanking_mass_tol_) continue;
          }

          if (tag.getCtermMass() > 0 && pos + tag.getSequence().length() < fe.sequence.length())
          {
            auto cterm = fe.sequence.substr(pos + tag.getSequence().length());
            if (x_pos != String::npos) cterm.erase(remove(cterm.begin(), cterm.end(), 'X'), cterm.end());

            double aamass = cterm.empty() ? 0 : AASequence::fromString(cterm).getMonoWeight();
            if (std::abs(tag.getCtermMass() - aamass) > flanking_mass_tol_) continue;
          }

          for (int off = 0; off < (int)tag.getLength(); off++)
          {
            int score = tag.getScore(off);
            auto iter = matched_pos_score.find(pos + off);
            if (iter != matched_pos_score.end()) score = std::max(score, iter->second);
            matched_pos_score[pos + off] = score;
            matched = true;
          }
        }
        if (matched)
        {
          matched_tag_indices.push_back(j); // tag indices
        }
        else
          continue;
#pragma omp critical
        if (! is_decoy) min_hit_tag_score = std::min(min_hit_tag_score, tag.getScore());
      }
      if (matched_tag_indices.empty()) continue;

      int match_cntr = 0;
      int match_score = 0;
      for (const auto& ps : matched_pos_score)
      {
        if (fe.sequence[ps.first] == 'X') continue;
        match_cntr++;
        match_score += ps.second;
      }

      if (match_cntr < min_cov_aa_) continue;
      //(double score, UInt rank, String accession, String sequence)
      ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
      hit.setDescription(fe.description);
      hit.setMetaValue("MatchedAA", match_cntr);
      hit.setMetaValue("IsDecoy", is_decoy ? 1 : 0);
      hit.setCoverage(double(match_cntr) / fe.sequence.length());
      hit.setScore(match_score);
#pragma omp critical
      {
        pairs.emplace_back(hit, matched_tag_indices);
      }
    }
  }
  if (pairs.empty()) return;

  protein_hits_.reserve(pairs.size());

  std::sort(pairs.begin(), pairs.end(),
            [](const std::pair<ProteinHit, std::vector<int>>& left, const std::pair<ProteinHit, std::vector<int>>& right) {
              return left.first.getScore() > right.first.getScore();
            });

  // FDR calculation
  double cum_target_count = 0;
  double cum_decoy_count = 0;

  decoy_mul /= fasta_entry.size() - decoy_mul;

  for (auto& [hit, indices] : pairs)
  {
    bool is_decoy = (int)hit.getMetaValue("IsDecoy") > 0;
    if (is_decoy) { cum_decoy_count += 1.0 / decoy_mul; }
    else { cum_target_count++; }

    double qvalue = decoy_mul != 0 ? (cum_decoy_count / (cum_target_count + cum_decoy_count)) : -1.0;

    hit.setMetaValue("qvalue", qvalue);
  }

  double min_qvalue = 1;
  for (auto iter = pairs.rbegin(); iter != pairs.rend(); iter++)
  {
    min_qvalue = std::min(min_qvalue, (double)iter->first.getMetaValue("qvalue"));
    iter->first.setMetaValue("qvalue", min_qvalue);
  }

  matching_tags_indices_.reserve(pairs.size());
  matching_hits_indices_ = std::vector<std::vector<int>>(tags_.size());

  for (const auto& [hit, indices] : pairs)
  {
    if ((double)hit.getMetaValue("qvalue") > fdr_) continue;
    if ((int)hit.getMetaValue("IsDecoy") > 0 && ! keep_decoy_) continue;

    protein_hits_.push_back(hit);
    matching_tags_indices_.push_back(indices);
    for (const auto& index : indices)
    {
      matching_hits_indices_[index].push_back(protein_hits_.size() - 1);
    }
  }
}

int FLASHTaggerAlgorithm::getProteinIndex(const ProteinHit& hit) const
{
  auto iter = std::find(protein_hits_.begin(), protein_hits_.end(), hit);
  if (iter == protein_hits_.end()) return -1;
  return std::distance(protein_hits_.begin(), iter);
}

int FLASHTaggerAlgorithm::getTagIndex(const FLASHDeconvHelperStructs::Tag& tag) const
{
  auto iter = std::find(tags_.begin(), tags_.end(), tag);
  if (iter == tags_.end()) return -1;
  return std::distance(tags_.begin(), iter);
}


const std::vector<ProteinHit>& FLASHTaggerAlgorithm::getProteinHits() const
{
  return protein_hits_;
}

const std::vector<ProteinHit> FLASHTaggerAlgorithm::getProteinHits(const FLASHDeconvHelperStructs::Tag& tag) const
{
  std::vector<ProteinHit> hits;
  int index = getTagIndex(tag);
  if (index < 0) return hits;
  for (auto i : matching_hits_indices_[index])
  {
    hits.push_back(protein_hits_[i]);
  }
  return hits;
}

const std::vector<FLASHDeconvHelperStructs::Tag>& FLASHTaggerAlgorithm::getTags() const
{
  return tags_;
}

std::vector<int> FLASHTaggerAlgorithm::getMatchedPositions(const ProteinHit& hit, const FLASHDeconvHelperStructs::Tag& tag) const
{
  Size pos = 0;
  std::vector<int> indices;
  auto seq = hit.getSequence();
  auto tagseq = tag.getSequence().toUpper();
  while (true)
  {
    pos = find_with_X_(seq, tagseq, pos + 1);
    if (pos == String::npos) break;
    indices.push_back((int)pos);
  }
  return indices;
}

std::vector<FLASHDeconvHelperStructs::Tag> FLASHTaggerAlgorithm::getTags(const ProteinHit& hit) const
{
  std::vector<FLASHDeconvHelperStructs::Tag> tags;
  int index = getProteinIndex(hit);
  if (index < 0) return tags;
  for (auto i : matching_tags_indices_[index])
  {
    tags.push_back(tags_[i]);
  }
  return tags;
}

} // namespace OpenMS