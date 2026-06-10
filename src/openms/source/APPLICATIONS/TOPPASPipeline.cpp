// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPASPipeline.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>

namespace OpenMS
{
  int TOPPASPipeline::addNode(const Node& n)
  {
    nodes_.push_back(n);
    // do not trust caller-provided edge lists for a freshly added node
    nodes_.back().in_edges.clear();
    nodes_.back().out_edges.clear();
    return (int)nodes_.size() - 1;
  }

  int TOPPASPipeline::addEdge(int from, int to, Int source_out_param, Int target_in_param)
  {
    Edge e;
    e.from = from;
    e.to = to;
    e.source_out_param = source_out_param;
    e.target_in_param = target_in_param;
    int idx = (int)edges_.size();
    edges_.push_back(e);
    nodes_[from].out_edges.push_back(idx);
    nodes_[to].in_edges.push_back(idx);
    return idx;
  }

  bool TOPPASPipeline::topoSort(std::string& error_msg)
  {
    const Size n = nodes_.size();
    std::vector<int> indeg(n, 0);
    for (const Edge& e : edges_)
    {
      if (e.to >= 0) ++indeg[e.to];
    }
    // process zero-indegree nodes in index order (deterministic)
    std::vector<int> queue;
    for (Size i = 0; i < n; ++i)
    {
      if (indeg[i] == 0) queue.push_back((int)i);
    }
    Size head = 0;
    Int topo = 1;
    Size processed = 0;
    while (head < queue.size())
    {
      int u = queue[head++];
      nodes_[u].topo_nr = topo++;
      ++processed;
      for (int ei : nodes_[u].out_edges)
      {
        int v = edges_[ei].to;
        if (--indeg[v] == 0) queue.push_back(v);
      }
    }
    if (processed != n)
    {
      error_msg = "TOPPASPipeline: workflow graph contains a cycle.";
      return false;
    }
    return true;
  }

  const TOPPASPipeline::Filenames& TOPPASPipeline::getFileNames(int node_index, Int param_index, int round) const
  {
    const Node& node = nodes_[node_index];
    if ((Size)round >= node.output_files.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, round, node.output_files.size());
    }
    const RoundPackage& rp = node.output_files[round];
    RoundPackage::const_iterator it = rp.find(param_index);
    if (it == rp.end())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, param_index, rp.size());
    }
    return it->second.filenames;
  }

  bool TOPPASPipeline::buildRoundPackages(int node_index, RoundPackages& pkg, std::string& error_msg) const
  {
    const Node& node = nodes_[node_index];
    if (node.in_edges.empty())
    {
      error_msg = "buildRoundPackages() called on vertex with no input edges!\n";
      return false;
    }

    // -- determine number of rounds from incoming, NON-recycling edges (they must all agree)
    int round_common = -1;
    int no_recycle_count = 0;
    for (int ei : node.in_edges)
    {
      const Node& up = nodes_[edges_[ei].from];
      if (up.recycle_output) continue;
      ++no_recycle_count;
      if (round_common == -1) round_common = up.round_total; // first non-recycler sets the pace
      if (round_common != up.round_total)
      {
        error_msg = "Number of rounds for incoming edges of node #" + std::to_string(node.topo_nr) +
                    " are not equal. No idea on how to combine them! Did you want to recycle its input?\n";
        return false;
      }
    }

    // -- at least one non-recycling input is required to anchor the round count
    if (no_recycle_count == 0)
    {
      error_msg = "Number of rounds of node #" + std::to_string(node.topo_nr) +
                  " cannot be determined since all input nodes have recycling enabled. Disable for at least one input!\n";
      return false;
    }

    // -- recyclers must tile evenly into round_common
    for (int ei : node.in_edges)
    {
      const Node& up = nodes_[edges_[ei].from];
      if (!up.recycle_output) continue;
      // (the round_total == 0 guard avoids a modulo-by-zero on a degenerate recycler that the
      //  historical TOPPASVertex::buildRoundPackages lacked; it turns UB into a clean error)
      if (up.round_total == 0 || round_common % up.round_total != 0)
      {
        error_msg = std::to_string(up.round_total) + " rounds for incoming edges of node #" + std::to_string(node.topo_nr) +
                    " are recycled to meet a total of " + std::to_string(round_common) +
                    " rounds. But modulo is not 0. No idea on how to combine them! Adapt the number of input files?\n";
        return false;
      }
    }

    if (round_common <= 0)
    {
      error_msg = "Number of input rounds is 0 or negative. This cannot be! Aborting!\n";
      return false;
    }

    pkg.clear();
    pkg.resize(round_common);

    for (int ei : node.in_edges)
    {
      const Edge& e = edges_[ei];
      const Node& up = nodes_[e.from];
      Int param_index_src_out = e.source_out_param;
      Int param_index_tgt_in = e.target_in_param;
      for (int round = 0; round < round_common; ++round)
      {
        VertexRoundPackage rpg;
        rpg.edge_index = ei;
        int upstream_round = round;
        if (up.recycle_output && upstream_round >= up.round_total)
        {
          upstream_round %= up.round_total;
        }
        rpg.filenames = getFileNames(e.from, param_index_src_out, upstream_round);

        // mergers have several incoming edges all keyed -1: spread over -1, -2, -3, ...
        while (pkg[round].count(param_index_tgt_in))
        {
          --param_index_tgt_in;
        }
        pkg[round][param_index_tgt_in] = rpg;
      }
    }
    return true;
  }

  void TOPPASPipeline::runInput_(Node& node) const
  {
    // each input file becomes one round (keyed by -1)
    node.output_files.clear();
    node.output_files.resize(node.input_files.size());
    for (Size f = 0; f < node.input_files.size(); ++f)
    {
      node.output_files[f][-1].filenames.push_back(node.input_files[f]);
    }
    node.round_total = (int)node.input_files.size();
    node.finished = true;
  }

  void TOPPASPipeline::runMerger_(Node& node, const RoundPackages& pkg) const
  {
    const Size input_rounds = pkg.size();
    node.round_total = node.round_based_merge ? (int)input_rounds : 1;
    node.output_files.clear();
    node.output_files.resize(node.round_total);
    for (Size round = 0; round < input_rounds; ++round)
    {
      Filenames files;
      for (RoundPackage::const_iterator it = pkg[round].begin(); it != pkg[round].end(); ++it)
      {
        files.insert(files.end(), it->second.filenames.begin(), it->second.filenames.end());
      }
      const Size round_index = node.round_based_merge ? round : 0;
      Filenames& dest = node.output_files[round_index][-1].filenames;
      dest.insert(dest.end(), files.begin(), files.end());
    }
    node.finished = true;
  }

  void TOPPASPipeline::runSplitter_(Node& node, const RoundPackages& pkg) const
  {
    // 1 round of N files becomes N rounds of 1 file
    node.output_files.clear();
    for (const RoundPackage& rp : pkg)
    {
      // there can only be one upstream (input) node
      const Filenames& files = rp.begin()->second.filenames;
      for (const std::string& f : files)
      {
        RoundPackage new_pkg;
        new_pkg[-1].filenames.push_back(f);
        node.output_files.push_back(new_pkg);
      }
    }
    node.round_total = (int)node.output_files.size();
    node.finished = true;
  }

  bool TOPPASPipeline::computeDataFlow(std::string& error_msg)
  {
    // visit nodes in topological order so every upstream is processed first
    std::vector<int> order(nodes_.size());
    for (Size i = 0; i < nodes_.size(); ++i) order[i] = (int)i;
    std::sort(order.begin(), order.end(), [this](int a, int b) { return nodes_[a].topo_nr < nodes_[b].topo_nr; });

    for (int ni : order)
    {
      Node& node = nodes_[ni];
      switch (node.kind)
      {
        case NodeKind::INPUT:
          runInput_(node);
          break;
        case NodeKind::MERGER:
        {
          RoundPackages pkg;
          if (!buildRoundPackages(ni, pkg, error_msg)) return false;
          runMerger_(node, pkg);
          break;
        }
        case NodeKind::SPLITTER:
        {
          RoundPackages pkg;
          if (!buildRoundPackages(ni, pkg, error_msg)) return false;
          runSplitter_(node, pkg);
          break;
        }
        case NodeKind::TOOL:
        {
          RoundPackages pkg;
          if (!buildRoundPackages(ni, pkg, error_msg)) return false;
          node.round_total = (int)pkg.size();
          // placeholder: the actual output files are produced by a PipelineExecutor
          node.output_files.assign(node.round_total, RoundPackage());
          node.finished = true;
          break;
        }
        case NodeKind::OUTPUT:
        case NodeKind::OUTPUT_FOLDER:
        {
          RoundPackages pkg;
          if (!buildRoundPackages(ni, pkg, error_msg)) return false;
          node.round_total = (int)pkg.size();
          node.finished = true;
          break;
        }
      }
    }
    return true;
  }

} // namespace OpenMS
