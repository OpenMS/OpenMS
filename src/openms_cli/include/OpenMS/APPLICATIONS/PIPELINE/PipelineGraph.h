// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/OpenMS_CLIConfig.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <string>
#include <vector>

namespace OpenMS
{
/**
  @brief Qt-independent workflow definition shared by ExecutePipeline and TOPPAS.

  Node identifiers are stable references, not vector offsets. Node and edge order
  are significant for legacy TOPPAS numbering and merger input ordering. Runtime
  status, processes and generated filenames belong to the executor.
*/
class OPENMS_CLI_DLLAPI PipelineGraph
{
public:
  /// A workflow vertex, including optional GUI layout metadata.
  struct Node
  {
    enum class Kind
    {
      INPUT,
      TOOL,
      MERGER,
      SPLITTER,
      OUTPUT,
      OUTPUT_DIRECTORY
    };
    Size id = 0;
    Kind kind = Kind::INPUT;
    std::string tool_name;
    std::string tool_type;
    Param parameters;
    std::vector<std::string> files;
    bool recycle_output = false;
    bool round_based = true;
    std::string resource_key;
    std::string output_folder;
    double x = 0;
    double y = 0;
    Size topo_number = 0;
  };

  using Kind = Node::Kind;

  /// A directed connection; empty ports represent the unnamed ports of structural nodes.
  struct Edge
  {
    Size source = 0;
    Size target = 0;
    std::string source_port;
    std::string target_port;
  };

  std::vector<Node> nodes;
  std::vector<Edge> edges;
  std::string version;
  std::string description;
  std::string filename;
  /// Older files without info:version encode tool ports as numeric indices.
  bool legacy_port_indices = false;

  /// Find a node by identifier; throws Exception::InvalidParameter if absent.
  const Node& node(Size id) const;
  Node& node(Size id);

  /// Return node identifiers in the historical TOPPAS stable scan order; reject cycles.
  std::vector<Size> topologicalOrder() const;

  /// Validate identifiers, topology and structural port bindings without launching tools.
  void validate() const;

  /// Assign one-based topological numbers and update default (but not custom) resource keys.
  void assignTopologicalNumbers();
};
} // namespace OpenMS
