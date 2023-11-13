// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>
#include <boost/dynamic_bitset.hpp>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
  /*
  class DAC
  {
  private:
    int vertice_count_;
    std::vector<std::vector<int>> adj_list_;
    std::vector<float> vertex_scores_;
  public:
    DAC(int vertice_count) : vertice_count_(vertice_count), adj_list_(vertice_count), vertex_scores_(vertice_count, .0f)
    {
    }

    void addEdge(int src, int dest, float src_score)
    {
      adj_list_[src].push_back(dest);
      vertex_scores_[src] = src_score;
    }

    void findAllPaths(int source, int destination, float min_score)
    {
      boost::dynamic_bitset<> visited(vertice_count_);
      std::vector<int> path;
      float score = 0.0; // where to start?
      std::vector<std::vector<int>> all_paths;
      findAllPaths_(source, destination, visited, path, score, min_score, all_paths);
    }

  private:
    float getScore_(int current, int destination)
    {
      return std::max(vertex_scores_[current], vertex_scores_[destination]);
    }

    void findAllPaths_(int current, int destination, boost::dynamic_bitset<>& visited, std::vector<int>& path, float score, float min_score, std::vector<std::vector<int>>& all_paths)
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
            float new_score = score + getScore_(current, neighbor);
            if (new_score >= min_score)
              findAllPaths_(neighbor, destination, visited, path, new_score, min_score, all_paths);
          }
        }
      }

      // Backtrack
      visited[current] = false;
      path.pop_back();
    }
  };


  struct edge {
    double mass;
    String aa;    // maybe empty for mass gaps and edges from sink or to source.
    double score; //
  };

  inline int getVertex(int index, int level, int max_level)
  {
    return index * (max_level + 1) + level;
  }

  const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

  bool isAA(PeakGroup& l, PeakGroup& r, double tol)
  {
    return true;
  }

  bool isGap(PeakGroup& l, PeakGroup& r, double tol)
  {
    return true;
  }

  void constructDAC(DAC& dac, DeconvolvedSpectrum& dspec, int min_tag_length, float min_score, double tol) // dac(dspec.size() * (min_tag_length + 1));
  {
    double max_edge_mass_ = 200; // TODO
    // from source to sink, connect but the edge direction is from sink to source.
    int start_index = 0;
    int end_index = 1;

    while (end_index < dspec.size())
    {
      auto r = dspec[end_index];
      float score = r.getQscore();
      // first, make edge from r to source and sink to r.
      dac.addEdge(getVertex(end_index, 0, min_tag_length), getVertex(0, 0, min_tag_length), score);
      // from an edge i, j to class edge.  for each i, j make a unique key. key to an edge.
      if (end_index < dspec.size() - 1) dac.addEdge(getVertex(dspec.size() - 1, 0, min_tag_length), getVertex(end_index, 0, min_tag_length), score);

      while (r.getMonoMass() - dspec[start_index].getMonoMass() > max_edge_mass_)
        start_index++;

      bool connected = false;
      for (int i = start_index; i < end_index; i++)
      {
        auto l = dspec[i];
        // make edge from r to l if they make an a.a. mass.
        if (isAA(l, r ,tol))
        {
          for (int lvl = 0; lvl < min_tag_length; lvl++)
          {
            dac.addEdge(getVertex(end_index, lvl + 1, min_tag_length), getVertex(i, lvl, min_tag_length), score);
          }
          dac.addEdge(getVertex(end_index, min_tag_length, min_tag_length), getVertex(i, min_tag_length, min_tag_length), score);
          connected = true;
        }
      }
      if (!connected)
      {
        for (int i = start_index; i < end_index; i++)
        {
          auto l = dspec[i];
          // make edge from r to l if they make a gap mass.
          if (isGap(l, r ,tol))
          {
            for (int lvl = 0; lvl < min_tag_length; lvl++)
            {
              dac.addEdge(getVertex(end_index, lvl + 1, min_tag_length), getVertex(i, lvl, min_tag_length), score);
            }
            dac.addEdge(getVertex(end_index, min_tag_length, min_tag_length), getVertex(i, min_tag_length, min_tag_length), score);
          }
        }
      }
      end_index++;
    }
  }
*/

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
    defaults_.setValue("min_length", 3, "Minimum length of the tags.");
    defaults_.setMinInt("min_length", 3);

    defaults_.setValue("max_length", 1000, "Maximum length of the tags.");
    defaults_.setMaxInt("max_length", 1000);

    defaults_.setValue("tol", DoubleList {10.0, 10.0}, "ppm tolerance for tag generation.");
    defaultsToParam_();
  }

  void TopDownTagger::updateMembers_()
  {
    min_tag_length_ = param_.getValue("min_length");
    max_tag_length_ = param_.getValue("max_length");
    max_tag_length_ = max_tag_length_ < min_tag_length_ ? min_tag_length_ : max_tag_length_;
    ppm_ = param_.getValue("tol");
  }

  void TopDownTagger::run(DeconvolvedSpectrum& dspec, std::vector<std::string>& tags)
  {
    double tol = ppm_[dspec.getOriginalSpectrum().getMSLevel() - 1];
    auto spec = dspec.toSpectrum(1, 10, tol);
    auto tagger = Tagger(min_tag_length_, tol, max_tag_length_, 1, 1);
    // tagger.setUseAbsoluteMzForTol();
    tagger.getTag(spec, tags);
  }
} // namespace OpenMS