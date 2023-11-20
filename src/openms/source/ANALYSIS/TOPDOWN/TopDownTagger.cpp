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

  std::vector<Residue> TopDownTagger::getAA_(double l, double r, double tol, int z) const
  {
    std::vector<Residue> ret;
    if (l == r)
      return ret;

    double diff = std::abs((r - l) * z);
    double abs_tol = 2 * std::abs(std::max(l, r) * tol / 1e6);
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

  std::vector<std::vector<Residue>> TopDownTagger::getGap_(double l, double r, double tol, int z) const
  {
    std::vector<std::vector<Residue>> ret;
    if (l == r)
      return ret;

    double diff = std::abs((r - l) * z);
    double abs_tol = 2 * std::abs(std::max(l, r) * tol / 1e6);
    auto iter = gap_mass_map_.lower_bound(diff - abs_tol);

    while (iter != gap_mass_map_.end())
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
    gap_mass_map_.clear();

    for (const auto& aa : aas_)
    {
      double aa_mass = aa->getMonoWeight(Residue::Internal);
      if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end())
        aa_mass_map_[aa_mass] = std::vector<Residue>();
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
    }

    if (allowed_isotope_error_ > 0)
    {
      for (int i = -allowed_isotope_error_; i <= allowed_isotope_error_; i++)
      {
        if (i == 0)
          continue;
        for (const auto& aa : aas_)
        {
          double aa_mass = aa->getMonoWeight(Residue::Internal) + i * Constants::C13C12_MASSDIFF_U;
          if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end())
            aa_mass_map_[aa_mass] = std::vector<Residue>();
          aa_mass_map_[aa_mass].push_back(*aa);
        }
      }
    }

    if (max_gap_count_ > 0)
    {
      std::map<double, std::vector<std::vector<Residue>>> tmp_gap_mass_map;
      for (auto& e : gap_mass_map_)
      {
        std::vector<std::vector<Residue>> tmp_e;
        for (auto&f : e.second)
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

  int TopDownTagger::getVertex_(int index, int path_score, int level, int gap_level) const
  {
    return ((index * (tag_length_ + 1) + level) * (max_gap_count_ + 1) + gap_level) * (max_path_score_ - min_path_score_ + 1) + (path_score - min_path_score_);
  }

  int TopDownTagger::getIndex_(int vertex) const
  {
    return ((vertex / (max_path_score_ - min_path_score_ + 1)) / (max_gap_count_ + 1)) / (tag_length_ + 1);
  }

  int TopDownTagger::edgeScore_(int vertex_score1, int vertex_score2)
  {
    return vertex_score1; // std::max(vertex_score1, vertex_score2);
  }

  bool TopDownTagger::connectEdge_(TopDownTagger::DAC_& dac, int vertex1, int vertex2, boost::dynamic_bitset<>& visited)
  {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= visited.size() || vertex2 >= visited.size())
      return false;
    if (!visited[vertex1] && !visited[vertex2])
      return false;

    dac.addEdge(vertex1, vertex2);
    return visited[vertex1] = visited[vertex2] = true;
  }


  void TopDownTagger::constructDAC_(TopDownTagger::DAC_& dac, const std::vector<double>& mzs, const std::vector<int>& scores, int z, double tol)
  {
    // from source to sink, connect but the edge direction is from sink to source.
    edge_aa_map_.clear();
    int start_index = 0;
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

      bool connected = false;
      for (int n = 0; n < 2; n++)
      {
        for (int current_index = start_index; current_index < end_index; current_index++)
        {
          auto l = mzs[current_index];
          int edge_score = edge_score_(scores[end_index], scores[current_index]);

          // make edge from r to l if they make an a.a. mass.
          std::vector<Residue> aas;
          std::vector<std::vector<Residue>> gaps;
          if (n == 0)
          {
            aas = getAA_(l, r, tol, z);
            if (aas.empty())
              continue;
          }
          else
          {
            gaps = getGap_(l, r, tol, z);
            if (gaps.empty())
              continue;
          }

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

          if (n == 0)
          {
            for (auto& aa : aas)
            {
              e[current_index].push_back(aa.toString());
            }
          }
          else
          {
            for (auto& gap : gaps)
            {
              std::stringstream ss;
              for (auto& aa : gap)
                ss << aa.toString().toLower();
              e[current_index].emplace_back(ss.str());
            }
          }
          int gap_diff = n == 0 ? 0 : 1;
          for (int g = 0; g + gap_diff <= max_gap_count_; g++)
          {
            for (int lvl = 0; lvl <= tag_length_; lvl++)
            {
              for (int score = min_path_score_; score <= max_path_score_; score++)
              {
                if (score - edge_score < min_path_score_)
                  continue;
                if (score - edge_score > max_path_score_)
                  break;

                int vertex1 = getVertex_(end_index, score, std::min(tag_length_, lvl + 1), g + gap_diff);
                int vertex2 = getVertex_(current_index, score - edge_score, lvl, g);
                bool con = connectEdge_(dac, vertex1, vertex2, visited);
                connected |= con;
              }
            }
          }
        }
        if (max_gap_count_ == 0)
          break;
      }

      if (end_index < mzs.size() - 1)
      {
        for (int g = 0; g <= max_gap_count_; g++)
        {
          for (int score = min_path_score_; score <= max_path_score_; score++)
          {
            int edge_score = edge_score_(scores[mzs.size() - 1], scores[end_index]);
            if (score - edge_score < min_path_score_)
              continue;
            if (score - edge_score > max_path_score_)
              break;

            int vertex1 = getVertex_(mzs.size() - 1, score, tag_length_, g);
            int vertex2 = getVertex_(end_index, score - edge_score, tag_length_, g);

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
                       "Maximum number of the tags. The tags with different amino acid and mass gap combinations are all treated separately. E.g., TII, TIL, TLI, TLL have the same mass differences "
                       "but are counted as four different tags. ");
    defaults_.setMinInt("max_tag_count", 0);

    defaults_.setValue("allowed_isotope_error", 0, "Allowed_isotope_error for tag generation. It only applies to amino acids, not mass gaps.");
    defaults_.setMaxInt("allowed_isotope_error", 1);
    defaults_.setMinInt("allowed_isotope_error", 0);

    defaults_.setValue("min_length", 4, "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids).");
    defaults_.setMaxInt("min_length", 20);
    defaults_.setMinInt("min_length", 3);


    defaults_.setValue("min_charge", 1, "Minimum charge state of the tags (can be negative for negative mode)");
    defaults_.setValue("max_charge", 1, "Maximum charge state of the tags (can be negative for negative mode)");

    defaults_.setValue("max_gap_count", 0, "Maximum mass gap count per tag.");
    defaults_.setMaxInt("max_gap_count", 2);
    defaults_.setMinInt("max_gap_count", 0);

    defaults_.setValue("max_aa_in_gap", 2, "Maximum amino acid count in a mass gap.");
    defaults_.setMaxInt("max_aa_in_gap", 3);
    defaults_.setMinInt("max_aa_in_gap", 2);


    defaultsToParam_();
  }

  void TopDownTagger::updateMembers_()
  {
    max_tag_count_ = param_.getValue("max_tag_count");
    tag_length_ = param_.getValue("min_length");
    min_charge_ = param_.getValue("min_charge");
    max_charge_ = param_.getValue("max_charge");
    max_aa_in_gap_ = param_.getValue("max_aa_in_gap");
    max_gap_count_ = param_.getValue("max_gap_count");
    allowed_isotope_error_ = param_.getValue("allowed_isotope_error");

    updateEdgeMasses_();
    max_edge_mass_ = gap_mass_map_.empty() ? aa_mass_map_.rbegin()->first : std::max(aa_mass_map_.rbegin()->first, gap_mass_map_.rbegin()->first);
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
      int score = (int)round(100 * log10(1e-6 + pg.getQscore()));
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
    double nmass = 0, cmass = 0;
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
    if (max_tag_count_ == 0)
      return;

    edge_score_ = edge_score;

    std::vector<double> _mzs;
    std::vector<int> _scores;
    _mzs.reserve(mzs.size() + 1);
    _scores.reserve(_scores.size() + 1);
    _mzs.push_back(.0);
    _scores.push_back(0);

    for (int i = 0; i < mzs.size(); i++)
    {
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

    max_path_score_ = edge_score_(max_vertex_score, max_vertex_score) * (tag_length_ + 3);
    min_path_score_ = edge_score_(min_vertex_score, min_vertex_score) * (tag_length_ + 3);

    std::vector<std::vector<int>> all_paths;
    for (int z = min_abs_charge; z <= max_abs_charge; z++)
    {
      TopDownTagger::DAC_ dac(_mzs.size() * (1 + tag_length_) * (1 + max_gap_count_) * (1 + max_path_score_ - min_path_score_));
      constructDAC_(dac, _mzs, _scores, z, ppm);
      std::set<FLASHDeconvHelperStructs::Tag> _tagSet;

      for (int score = max_path_score_; score >= min_path_score_ && _tagSet.size() < max_tag_count_; score--)
      {
        for (int g = 0; g <= max_gap_count_; g++)
        {
          dac.findAllPaths(getVertex_(_mzs.size() - 1, score, tag_length_, g), getVertex_(0, 0, 0, 0), all_paths);
        }
        for (const auto& path : all_paths)
        {
          updateTagSet_(_tagSet, path, _mzs, z, score);
        }
      }

      for (const auto& tag : _tagSet)
      {
        tags.push_back(tag);
      }
    }

    std::sort(tags.begin(), tags.end(), [](FLASHDeconvHelperStructs::Tag& a, FLASHDeconvHelperStructs::Tag& b) { return a.getScore() > b.getScore(); });

    while (tags.size() > max_tag_count_)
    {
      tags.pop_back();
    }
  }
} // namespace OpenMS