// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineGraph.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <algorithm>
#include <map>
#include <set>
#include <tuple>

namespace OpenMS
{
namespace
{
  [[noreturn]] void invalidGraph(const std::string& message)
  { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message); }

  bool isOutput(PipelineGraph::Kind kind)
  { return kind == PipelineGraph::Kind::OUTPUT || kind == PipelineGraph::Kind::OUTPUT_DIRECTORY; }
} // namespace

const PipelineGraph::Node& PipelineGraph::node(Size id) const
{
  const auto it = std::find_if(nodes.begin(), nodes.end(), [id](const Node& value) { return value.id == id; });
  if (it == nodes.end()) { invalidGraph("Unknown pipeline node " + StringUtils::toStr(id)); }
  return *it;
}

PipelineGraph::Node& PipelineGraph::node(Size id)
{
  const auto it = std::find_if(nodes.begin(), nodes.end(), [id](const Node& value) { return value.id == id; });
  if (it == nodes.end()) { invalidGraph("Unknown pipeline node " + StringUtils::toStr(id)); }
  return *it;
}

std::vector<Size> PipelineGraph::topologicalOrder() const
{
  std::map<Size, Size> index;
  for (Size i = 0; i < nodes.size(); ++i)
  {
    if (! index.emplace(nodes[i].id, i).second) { invalidGraph("Duplicate pipeline node identifier " + StringUtils::toStr(nodes[i].id)); }
  }
  std::vector<Size> incoming(nodes.size(), 0);
  std::vector<std::vector<Size>> successors(nodes.size());
  for (const auto& edge : edges)
  {
    const auto source = index.find(edge.source);
    const auto target = index.find(edge.target);
    if (source == index.end() || target == index.end())
    {
      invalidGraph("Pipeline edge references an unknown node: " + StringUtils::toStr(edge.source) + "/" + StringUtils::toStr(edge.target));
    }
    ++incoming[target->second];
    successors[source->second].push_back(target->second);
  }

  std::vector<Size> result;
  std::vector<bool> visited(nodes.size(), false);
  result.reserve(nodes.size());
  while (result.size() < nodes.size())
  {
    const Size previous_size = result.size();
    // Deliberately match TOPPAS' scan order, including nodes made ready within a scan.
    // A priority queue or a breadth-first traversal would change legacy resource keys.
    for (Size i = 0; i < nodes.size(); ++i)
    {
      if (visited[i] || incoming[i] != 0) { continue; }
      visited[i] = true;
      result.push_back(nodes[i].id);
      for (const Size successor : successors[i])
      {
        --incoming[successor];
      }
    }
    if (result.size() == previous_size) { invalidGraph("Pipeline contains a directed cycle."); }
  }
  return result;
}

void PipelineGraph::validate() const
{
  topologicalOrder(); // Also checks unique IDs and all edge endpoints.
  std::set<std::string> resource_keys;
  for (const auto& vertex : nodes)
  {
    switch (vertex.kind)
    {
      case Kind::INPUT:
        if (! vertex.resource_key.empty())
        {
          if (vertex.resource_key.find(':') != std::string::npos) { invalidGraph("Input resource key must not contain ':': " + vertex.resource_key); }
          if (! resource_keys.insert(vertex.resource_key).second) { invalidGraph("Duplicate input resource key: " + vertex.resource_key); }
        }
        break;
      case Kind::TOOL:
        if (vertex.tool_name.empty()) { invalidGraph("Pipeline tool node has no tool name."); }
        break;
      case Kind::MERGER:
      case Kind::SPLITTER:
      case Kind::OUTPUT:
      case Kind::OUTPUT_DIRECTORY:
        break;
      default:
        invalidGraph("Unknown pipeline node kind.");
    }
  }

  std::set<std::tuple<Size, Size, std::string, std::string>> bindings;
  std::set<std::pair<Size, std::string>> tool_inputs;
  std::map<Size, Size> incoming;
  for (const auto& edge : edges)
  {
    const auto& source = node(edge.source);
    const auto& target = node(edge.target);
    if (target.kind == Kind::INPUT || isOutput(source.kind)) { invalidGraph("Pipeline edge has an invalid node direction."); }
    if (isOutput(target.kind) && source.kind != Kind::TOOL) { invalidGraph("An output node must receive a tool output directly."); }
    if (++incoming[edge.target] > 1 && (isOutput(target.kind) || target.kind == Kind::SPLITTER))
    {
      invalidGraph("An output or splitter node must have only one incoming edge.");
    }
    if (! bindings.emplace(edge.source, edge.target, edge.source_port, edge.target_port).second) { invalidGraph("Duplicate pipeline edge."); }
    if (target.kind == Kind::TOOL && ! tool_inputs.emplace(edge.target, edge.target_port).second)
    {
      invalidGraph("Multiple edges bind the same tool input: " + target.tool_name + ":" + edge.target_port);
    }
    if (source.kind != Kind::TOOL && ! edge.source_port.empty()) { invalidGraph("A structural source node cannot have a named output port."); }
    if (target.kind != Kind::TOOL && ! edge.target_port.empty()) { invalidGraph("A structural target node cannot have a named input port."); }
  }
}

void PipelineGraph::assignTopologicalNumbers()
{
  const auto order = topologicalOrder();
  Size number = 0;
  for (const Size id : order)
  {
    auto& vertex = node(id);
    ++number;
    if (vertex.kind == Kind::INPUT && (vertex.resource_key.empty() || vertex.resource_key == StringUtils::toStr(vertex.topo_number)))
    {
      vertex.resource_key = StringUtils::toStr(number);
    }
    vertex.topo_number = number;
  }
}
} // namespace OpenMS
