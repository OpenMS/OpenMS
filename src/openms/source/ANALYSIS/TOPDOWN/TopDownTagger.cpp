// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>
#include <utility>

namespace OpenMS
{
  class TopDownTagger::DAC_
  {
  private:
    int vertex_count_;
    std::vector<std::vector<int>> adj_list_;

  public:
    DAC_(int vertice_count) : vertex_count_(vertice_count), adj_list_(vertice_count)
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

    void findAllPaths(int source, int sink, std::vector<std::vector<int>>& all_paths)
    {
      boost::dynamic_bitset<> visited(vertex_count_);
      std::vector<int> path;

      findAllPaths_(source, sink, visited, path, all_paths); // reverse traveling
    }

  private:
    void findAllPaths_(int current, int destination, boost::dynamic_bitset<>& visited, std::vector<int>& path, std::vector<std::vector<int>>& all_paths)
    {
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
          if (!visited[neighbor])
          {
            findAllPaths_(neighbor, destination, visited, path, all_paths);
          }
        }
      }

      // Backtrack
      visited[current] = false;
      path.pop_back();
    }
  };

  std::vector<Residue> TopDownTagger::getAA_(double l, double r, double tol, int z, int iso_offset) const
  {
    std::vector<Residue> ret;
    if (l == r)
      return ret;
    double abs_iso_mass = std::abs(iso_offset * Constants::C13C12_MASSDIFF_U);
    double diff = std::abs(std::abs((r - l)) - abs_iso_mass) / z;
    double abs_tol = 2 * std::max(l, r) * tol / 1e6;
    auto iter = aa_mass_map_.lower_bound(diff - abs_tol);

    while (iter != aa_mass_map_.end())
    {
      if (std::abs(diff - iter->first) < abs_tol)
      {
        for (auto& aa : iter->second)
          ret.push_back(aa);
      }
      else if (iter->first - diff > abs_tol)
      {
        break;
      }
      iter++;
    }

    return ret;
  }

  void TopDownTagger::updateEdgeMasses_()
  {
    aa_mass_map_.clear();
    //gap_mass_map_.clear();

    for (const auto& aa : aas_)
    {
      double aa_mass = aa->getMonoWeight(Residue::Internal);
      if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end())
        aa_mass_map_[aa_mass] = std::vector<Residue>();
      aa_mass_map_[aa_mass].push_back(*aa);
    }
  }

  int TopDownTagger::getVertex_(int index, int path_score, int level, int iso_level) const
  {
    return ((index * (max_tag_length_ + 1) + level) * (max_iso_in_tag_ + 1) + iso_level) * (max_path_score_ - min_path_score_ + 1) + (path_score - min_path_score_);
  }

  int TopDownTagger::getIndex_(int vertex) const
  {
    return ((vertex / (max_path_score_ - min_path_score_ + 1)) / (max_iso_in_tag_ + 1)) / (max_tag_length_ + 1);
  }

  int TopDownTagger::edgeScore_(int vertex_score1, int vertex_score2)
  {
    return vertex_score1 + vertex_score2;
  }

  bool TopDownTagger::connectEdge_(TopDownTagger::DAC_& dac, int vertex1, int vertex2, boost::dynamic_bitset<>& visited)
  {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= visited.size() || vertex2 >= visited.size())
      return false;
    if (!visited[vertex2])
      return false;

    dac.addEdge(vertex1, vertex2);
    return visited[vertex1] = true;
  }


  void TopDownTagger::constructDAC_(TopDownTagger::DAC_& dac, const std::vector<double>& mzs, const std::vector<int>& scores, int z, int length, double tol)
  {
    // from source to sink, connect but the edge direction is from sink to source.
    edge_aa_map_.clear();
    int start_index = 1; // zeroth = source.
    int end_index = 1;
    boost::dynamic_bitset<> visited(dac.size());
    visited[getVertex_(0, 0, 0, 0)] = true;

    while (end_index < mzs.size())
    {
      auto r = mzs[end_index];

      // first, make edge from r to source and sink to r.
      int vertex1 = getVertex_(end_index, edge_score_(scores[end_index], scores[0]), 0, 0);
      int vertex2 = getVertex_(0, 0, 0, 0);

      connectEdge_(dac, vertex1, vertex2, visited);

      // from an edge i, j to class edge.  for each i, j make a unique key. key to an edge.

      while (start_index < end_index && r - mzs[start_index] > max_edge_mass_)
        start_index++;

      for (int n = 0; n < 2; n++) // 0 for all a.a 1 for isotope errors
      {
        for (int current_index = start_index; current_index < end_index; current_index++)
        {
          auto l = mzs[current_index];
          int edge_score = edge_score_(scores[end_index], scores[current_index]);

          // make edge from r to l if they make an a.a. mass.
          std::vector<Residue> aas;

          aas = getAA_(l, r, tol, z, n);
          if (aas.empty())
            continue;


          // end_index, current_index to amino acid strings.
          if (edge_aa_map_.find(end_index) == edge_aa_map_.end())
          {
            edge_aa_map_[end_index] = std::map<int, std::vector<String>>();
          }
          auto& e = edge_aa_map_[end_index];

          if (e.find(current_index) == e.end())
          {
            e[current_index] = std::vector<String>();
          }

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
                if (score - edge_score < min_path_score_)
                  continue;
                if (score - edge_score > max_path_score_)
                  break;

                int vertex1 = getVertex_(end_index, score, lvl + 1, g + n);
                int vertex2 = getVertex_(current_index, score - edge_score, lvl, g);
                connectEdge_(dac, vertex1, vertex2, visited);
              }
            }
          }
        }
        if (max_iso_in_tag_ == 0)
          break;
      }

      if (end_index < mzs.size() - 1)
      {
        for (int g = 0; g <= max_iso_in_tag_; g++)
        {
          for (int score = min_path_score_; score <= max_path_score_; score++)
          {
            int edge_score = edge_score_(scores[mzs.size() - 1], scores[end_index]);
            if (score - edge_score < min_path_score_)
              continue;
            if (score - edge_score > max_path_score_)
              break;

            int vertex1 = getVertex_(mzs.size() - 1, score, length, g);
            int vertex2 = getVertex_(end_index, score - edge_score, length, g);
            connectEdge_(dac, vertex1, vertex2, visited);
          }
        }
      }

      end_index++;
    }
  }

  TopDownTagger::TopDownTagger() : DefaultParamHandler("TopDownTagger")
  {
    setDefaultParams_();
  }

  TopDownTagger::TopDownTagger(const TopDownTagger& other) : DefaultParamHandler(other)
  {
  }

  TopDownTagger& TopDownTagger::operator=(const TopDownTagger& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void TopDownTagger::setDefaultParams_()
  {
    defaults_.setValue("max_tag_count", 0,
                       "Maximum number of the tags per length (lengths set by -min_length and -max_length options). The tags with different amino acid combinations are all treated separately. E.g., "
                       "TII, TIL, TLI, TLL are distinct tags even though they have the same mass differences. "
                       "but are counted as four different tags. ");
    defaults_.setMinInt("max_tag_count", 0);

    defaults_.setValue("min_length", 4, "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
    defaults_.setMaxInt("min_length", 30);
    defaults_.setMinInt("min_length", 3);

    defaults_.setValue("max_length", 10, "Maximum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
    defaults_.setMaxInt("max_length", 30);
    defaults_.setMinInt("max_length", 3);

    defaults_.setValue("flanking_mass_tol", 200.0, "Flanking mass tolerance in Da.");

    defaults_.setValue("min_charge", 1, "Minimum charge state of the tags (can be negative for negative mode)");
    defaults_.setValue("max_charge", 1, "Maximum charge state of the tags (can be negative for negative mode)");

    defaults_.setValue("max_iso_error_count", 0, "Maximum isotope error count per tag.");
    defaults_.setMaxInt("max_iso_error_count", 2);
    defaults_.setMinInt("max_iso_error_count", 0);
    defaults_.addTag("max_iso_error_count", "advanced");

    defaultsToParam_();
  }

  void TopDownTagger::updateMembers_()
  {
    max_tag_count_ = param_.getValue("max_tag_count");
    min_tag_length_ = param_.getValue("min_length");
    max_tag_length_ = param_.getValue("max_length");
    min_charge_ = param_.getValue("min_charge");
    max_charge_ = param_.getValue("max_charge");
    max_iso_in_tag_ = param_.getValue("max_iso_error_count");

    updateEdgeMasses_();
    max_edge_mass_ = aa_mass_map_.rbegin()->first + max_iso_in_tag_ * Constants::C13C12_MASSDIFF_U;
  }

  void TopDownTagger::run(const DeconvolvedSpectrum& dspec, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags)
  {
    std::vector<double> mzs;
    std::vector<int> scores;
    mzs.reserve(dspec.size());
    scores.reserve(dspec.size());

    for (auto& pg : dspec)
    {
      mzs.push_back(pg.getMonoMass());
      int score = (int)round(100 * log10(1e-6 + pg.getQscore2D()));
      scores.push_back(score); //
    }
    run(mzs, scores, ppm, tags);
  }

  void TopDownTagger::run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags)
  {
    run(mzs, scores, ppm, tags, edgeScore_);
  }

  void TopDownTagger::updateTagSet_(std::set<FLASHDeconvHelperStructs::Tag>& tag_set, const std::vector<int>& path, const std::vector<double>& mzs, int z, int score)
  {
    double nmass = -1, cmass = -1;
    bool is_positive = min_charge_ > 0;

    std::vector<String> seqs {""};
    std::vector<double> tag_mzs;
    tag_mzs.reserve(path.size() - 1);

    for (int j = 1; j < path.size(); j++)
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
      }
      else if (i2 == 0) // nterm
      {
        tag_mzs.push_back(mzs[i1]);
        nmass = mzs[i1] / z;
      }
    }

    std::vector<double> rev_tag_mzs(tag_mzs);
    std::reverse(tag_mzs.begin(), tag_mzs.end());

    for (const auto& seq : seqs)
    {
      auto direct_tag = FLASHDeconvHelperStructs::Tag(seq, nmass, cmass, is_positive ? z : -z, std::pow(10.0, score / 100.0), tag_mzs);
      auto reverse_tag = FLASHDeconvHelperStructs::Tag(String(seq).reverse(), cmass, nmass, is_positive ? z : -z, std::pow(10.0, score / 100.0), rev_tag_mzs);

      auto iter = tag_set.find(direct_tag);
      if (iter == tag_set.end())
      {
        tag_set.insert(direct_tag);
      }
      else if (iter->getScore() < direct_tag.getScore())
      {
        tag_set.erase(iter);
        tag_set.insert(direct_tag);
      }

      iter = tag_set.find(reverse_tag);
      if (iter == tag_set.end())
      {
        tag_set.insert(reverse_tag);
      }
      else if (iter->getScore() < reverse_tag.getScore())
      {
        tag_set.erase(iter);
        tag_set.insert(reverse_tag);
      }
    }
  }

  void TopDownTagger::run(const std::vector<double>& mzs, const std::vector<int>& scores, double ppm, std::vector<FLASHDeconvHelperStructs::Tag>& tags, const std::function<int(int, int)>& edge_score)
  {
    const Size max_node_cntr = 500;
    if (max_tag_count_ == 0)
      return;

    edge_score_ = edge_score;

    std::vector<double> _mzs;
    std::vector<int> _scores;
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
      _scores.reserve(_scores.size() + 1);
      threshold = *std::min_element(scores.begin(), scores.end());
    }

    _mzs.push_back(.0);
    _scores.push_back(0);

    for (int i = 0; i < mzs.size(); i++)
    {
      if (scores[i] < threshold)
        continue;
      _mzs.push_back(mzs[i]);
      _scores.push_back(scores[i]);
    }
    int min_abs_charge = abs(min_charge_);
    int max_abs_charge = abs(max_charge_);

    if (min_abs_charge > max_abs_charge)
    {
      int tmp = min_abs_charge;
      min_abs_charge = max_abs_charge;
      max_abs_charge = tmp;
    }

    int max_vertex_score = *std::max_element(_scores.begin(), _scores.end());
    int min_vertex_score = *std::min_element(_scores.begin(), _scores.end());

    max_path_score_ = edge_score_(max_vertex_score, max_vertex_score) * (max_tag_length_ + 2);
    min_path_score_ = edge_score_(min_vertex_score, min_vertex_score) * (max_tag_length_ + 2);

    max_path_score_ = std::max(max_path_score_, edge_score_(max_vertex_score, max_vertex_score) * (min_tag_length_ - 2));
    min_path_score_ = std::min(min_path_score_, edge_score_(min_vertex_score, min_vertex_score) * (min_tag_length_ - 2));

    std::set<FLASHDeconvHelperStructs::Tag> tagSet;

    for (int z = min_abs_charge; z <= max_abs_charge; z++)
    {
      for (int length = min_tag_length_; length <= max_tag_length_; length++)
      {
        TopDownTagger::DAC_ dac(_mzs.size() * (1 + max_tag_length_) * (1 + max_iso_in_tag_) * (1 + max_path_score_ - min_path_score_));
        constructDAC_(dac, _mzs, _scores, z, length, ppm);
        std::vector<std::vector<int>> all_paths;
        std::set<FLASHDeconvHelperStructs::Tag> _tagSet;
        for (int score = max_path_score_; score >= min_path_score_ && _tagSet.size() < max_tag_count_; score--)
        {
          for (int g = 0; g <= max_iso_in_tag_; g++)
          {
            dac.findAllPaths(getVertex_(_mzs.size() - 1, score, length, g), getVertex_(0, 0, 0, 0), all_paths);
          }
          for (const auto& path : all_paths)
          {
            updateTagSet_(_tagSet, path, _mzs, z, score);
          }
        }
        tagSet.insert(_tagSet.begin(), _tagSet.end());
      }
    }

    for (int length = min_tag_length_; length <= max_tag_length_; length++)
    {
      int count = 0;
      for (const auto& tag : tagSet)
      {
        if (tag.getLength() != length)
          continue;
        tags.push_back(tag);
        if (++count == max_tag_count_)
          break;
      }
      std::sort(tags.begin(), tags.end(), [](const FLASHDeconvHelperStructs::Tag& a, const FLASHDeconvHelperStructs::Tag& b) {
        return a.getLength() == b.getLength() ? a.getScore() > b.getScore() : a.getLength() < b.getLength();
      });
      OPENMS_LOG_INFO << "Tag count with length " << length << ": " << count << std::endl;
    }
  }

  std::vector<std::pair<FASTAFile::FASTAEntry, std::vector<FLASHDeconvHelperStructs::Tag>>> TopDownTagger::runMatching(const std::vector<FLASHDeconvHelperStructs::Tag>& tags,
                                                                                                                       const std::vector<FASTAFile::FASTAEntry>& fasta_entry) const
  {
    std::vector<std::pair<FASTAFile::FASTAEntry, std::vector<FLASHDeconvHelperStructs::Tag>>> pairs;
#pragma omp parallel for default(none) shared(pairs, fasta_entry, tags)
    for (int i = 0; i < fasta_entry.size(); i++)
    {
      auto& fe = fasta_entry[i];
      std::vector<FLASHDeconvHelperStructs::Tag> matched_tags;
      auto seq = fe.sequence;
      for (auto& tag : tags)
      {
        bool matched = !seq.empty() && seq.hasSubstring(tag.getSequence().toUpper());

        if (matched)
        {
          auto pos = seq.find(tag.getSequence().toUpper());
          if (tag.getNtermMass() > 0)
          {
            auto nterm = seq.substr(0, pos);
            double aamass = nterm.empty() ? 0 : AASequence::fromString(nterm).getMonoWeight();
            if (std::abs(tag.getNtermMass() - aamass) > flanking_mass_tol_)
              matched = false;
          }
          if (matched && tag.getCtermMass() > 0)
          {
            auto cterm = seq.substr(pos + tag.getSequence().length());
            double aamass = cterm.empty() ? 0 : AASequence::fromString(cterm).getMonoWeight();
            if (std::abs(tag.getCtermMass() - aamass) > flanking_mass_tol_)
              matched = false;
          }
        }

        if (matched)
        {
          matched_tags.push_back(tag);
        }
      }
      if (matched_tags.empty())
        continue;
#pragma omp critical
      pairs.emplace_back(fe, matched_tags);
    }

    std::sort(pairs.begin(), pairs.end(),
              [](const std::pair<FASTAFile::FASTAEntry, std::vector<FLASHDeconvHelperStructs::Tag>>& left, const std::pair<FASTAFile::FASTAEntry, std::vector<FLASHDeconvHelperStructs::Tag>>& right) {
                return left.second.size() > right.second.size();
              });
    return pairs;
  }
} // namespace OpenMS