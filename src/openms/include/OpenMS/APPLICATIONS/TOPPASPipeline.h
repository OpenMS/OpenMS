// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Qt-free data-flow ("round package") engine for TOPPAS workflows.

    This class implements the round-based data-flow routing that TOPPAS uses to decide,
    for every node and every execution round, which files arrive on which input parameter.
    It is the Qt-free extraction of the routing logic previously living in the
    @c TOPPAS*Vertex classes (@c buildRoundPackages, the merger/splitter/input semantics and
    output recycling), so that a pipeline can be scheduled and executed without the GUI.

    The engine operates on an explicit node/edge graph (integer parameter indices, plain
    @c std types) and is deliberately decoupled from how nodes are parsed (from a @em .toppas
    document) and how tools are executed (a PipelineExecutor built on top).

    @par Data model
    A node's data flow is described by RoundPackages: one RoundPackage per execution round,
    each mapping a (signed) parameter index to the files arriving for it. The signed key
    mirrors the historic TOPPAS convention: @c -1 means "no specific parameter" (used by
    input/output nodes and mergers); mergers walk to @c -2, @c -3, ... when several incoming
    edges share the @c -1 slot.

    @ingroup Applications
  */
  class OPENMS_DLLAPI TOPPASPipeline
  {
  public:
    /// File list arriving along one edge for one round.
    using Filenames = std::vector<std::string>;

    /// Files for one (node, round, parameter) plus the edge they came from.
    struct VertexRoundPackage
    {
      Filenames filenames;
      int edge_index = -1; ///< index into edges_ of the producing edge (-1 if none)
    };

    /// One round's inputs/outputs of a node: parameter index -> files. Key may be negative.
    using RoundPackage = std::map<Int, VertexRoundPackage>;
    /// One entry per execution round.
    using RoundPackages = std::vector<RoundPackage>;

    /// Kind of a workflow node.
    enum class NodeKind
    {
      INPUT,         ///< seeds input files (one round per file)
      TOOL,          ///< a TOPP tool invocation
      MERGER,        ///< merges incoming file lists (round-based or collect)
      SPLITTER,      ///< splits one round of N files into N rounds of 1 file
      OUTPUT,        ///< collects output files
      OUTPUT_FOLDER  ///< collects output into a folder
    };

    /// A workflow node with its topology and data-flow state.
    struct Node
    {
      NodeKind kind = NodeKind::TOOL;
      std::string id;                  ///< node id (as in the .toppas file)
      Int topo_nr = 0;                 ///< topological number (1-based, like TOPPAS)
      bool recycle_output = false;     ///< recycle this node's output to match downstream round counts
      bool round_based_merge = true;   ///< merger only: round-based ("Merge") vs collect ("Collect")

      Filenames input_files;           ///< INPUT node: the seed files

      std::vector<int> in_edges;       ///< indices into edges_ (incoming)
      std::vector<int> out_edges;      ///< indices into edges_ (outgoing)

      // --- data-flow state (filled by the engine) ---
      RoundPackages output_files;      ///< files this node produces, per round and output parameter
      int round_total = -1;            ///< number of rounds this node runs
      bool finished = false;
    };

    /// A directed edge connecting an output parameter of one node to an input parameter of another.
    struct Edge
    {
      int from = -1;             ///< index into nodes_ (source)
      int to = -1;               ///< index into nodes_ (target)
      Int source_out_param = -1; ///< output-parameter index on the source node (-1 if none)
      Int target_in_param = -1;  ///< input-parameter index on the target node (-1 if none)
    };

    /// Add a node and return its index.
    int addNode(const Node& n);
    /// Add an edge (by node indices) and wire it into both endpoints' in/out lists; returns its index.
    int addEdge(int from, int to, Int source_out_param = -1, Int target_in_param = -1);

    /// Access nodes/edges.
    std::vector<Node>& nodes() { return nodes_; }
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Edge>& edges() const { return edges_; }

    /**
      @brief Topologically order the nodes and assign @c topo_nr (1-based).
      @return false if the graph contains a cycle (then @p error_msg is set).
    */
    bool topoSort(std::string& error_msg);

    /**
      @brief Build the per-round input packages for @p node_index from its incoming edges.

      Faithful port of @c TOPPASVertex::buildRoundPackages: the round count is set by the
      non-recycling upstream nodes (all must agree), recycling upstream nodes must tile evenly
      into it (and are wrapped via modulo), and several incoming edges sharing the @c -1 slot
      (mergers) are spread over @c -1, @c -2, ... .
      @return false on an inconsistent graph (then @p error_msg is set).
    */
    bool buildRoundPackages(int node_index, RoundPackages& pkg, std::string& error_msg) const;

    /**
      @brief Compute the data flow for all non-tool nodes in topological order.

      Seeds INPUT nodes, performs the (virtual) merge/split, and propagates recycling, filling
      each node's @c output_files / @c round_total. TOOL nodes get their @c round_total and a
      placeholder @c output_files (one empty RoundPackage per round); a PipelineExecutor is
      responsible for actually producing a tool's output files. @c topoSort() must have run.
      @return false on an inconsistent graph (then @p error_msg is set).
    */
    bool computeDataFlow(std::string& error_msg);

    /// Files produced by @p node_index for output parameter @p param_index in @p round.
    /// Throws Exception::IndexOverflow if the round/parameter is absent.
    const Filenames& getFileNames(int node_index, Int param_index, int round) const;

  protected:
    /// Seed an INPUT node: one round per input file, keyed by -1.
    void runInput_(Node& node) const;
    /// Virtual merge of @p node from its built round packages.
    void runMerger_(Node& node, const RoundPackages& pkg) const;
    /// Virtual split of @p node (one round of N files becomes N rounds of one file).
    void runSplitter_(Node& node, const RoundPackages& pkg) const;

    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
  };

} // namespace OpenMS
