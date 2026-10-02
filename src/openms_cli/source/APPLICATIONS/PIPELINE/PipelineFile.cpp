// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <charconv>
#include <filesystem>
#include <map>
#include <set>

namespace OpenMS
{
namespace
{
  [[noreturn]] void invalidFile(const std::string& message)
  { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message); }

  std::string pathString(const std::filesystem::path& path)
  {
    const auto value = path.generic_u8string();
    return {value.begin(), value.end()};
  }

  std::filesystem::path utf8Path(const std::string& value)
  { return std::filesystem::path(std::u8string(value.begin(), value.end())); }

  const ParamValue& requiredValue(const Param& parameters, const std::string& key)
  {
    if (! parameters.exists(key)) { invalidFile("Missing TOPPAS field '" + key + "'."); }
    return parameters.getValue(key);
  }

  Size indexValue(const std::string& value, const std::string& context)
  {
    Size result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (value.empty() || parsed.ec != std::errc() || parsed.ptr != value.data() + value.size())
    {
      invalidFile("Invalid nonnegative integer for " + context + ": '" + value + "'.");
    }
    return result;
  }

  bool boolValue(const Param& parameters, const std::string& key, bool default_value)
  {
    if (! parameters.exists(key)) { return default_value; }
    const std::string value = parameters.getValue(key).toString();
    if (value == "true") { return true; }
    if (value == "false") { return false; }
    invalidFile("Invalid boolean for TOPPAS field '" + key + "'.");
  }

  std::vector<std::string> sectionNames(const Param& parameters)
  {
    std::vector<std::string> names;
    std::set<std::string> seen;
    for (auto it = parameters.begin(); it != parameters.end(); ++it)
    {
      const std::string name = it.getName();
      const auto colon = name.find(':');
      if (colon == std::string::npos) { invalidFile("Expected TOPPAS section, found '" + name + "'."); }
      const auto section = name.substr(0, colon);
      if (seen.insert(section).second) { names.push_back(section); }
    }
    return names;
  }

  std::string portName(const std::string& value, bool legacy)
  {
    if (value.empty() || value == "__no_name__" || (legacy && value == "-1")) { return ""; }
    if (legacy) { indexValue(value, "legacy port index"); }
    return value;
  }

  PipelineGraph::Kind nodeKind(const std::string& name)
  {
    if (name == "input file list") { return PipelineGraph::Kind::INPUT; }
    if (name == "tool") { return PipelineGraph::Kind::TOOL; }
    if (name == "merger") { return PipelineGraph::Kind::MERGER; }
    if (name == "splitter") { return PipelineGraph::Kind::SPLITTER; }
    if (name == "output file list") { return PipelineGraph::Kind::OUTPUT; }
    if (name == "output folder") { return PipelineGraph::Kind::OUTPUT_DIRECTORY; }
    invalidFile("Unknown TOPPAS node type '" + name + "'.");
  }

  std::string nodeKindName(PipelineGraph::Kind kind)
  {
    switch (kind)
    {
      case PipelineGraph::Kind::INPUT:
        return "input file list";
      case PipelineGraph::Kind::TOOL:
        return "tool";
      case PipelineGraph::Kind::MERGER:
        return "merger";
      case PipelineGraph::Kind::SPLITTER:
        return "splitter";
      case PipelineGraph::Kind::OUTPUT:
        return "output file list";
      case PipelineGraph::Kind::OUTPUT_DIRECTORY:
        return "output folder";
    }
    invalidFile("Unknown pipeline node kind.");
  }

  std::string decodeFileURL(const std::string& url)
  {
    const auto colon = url.find(':');
    if (colon == std::string::npos || StringUtils::toLowered(url.substr(0, colon)) != "file")
    {
      invalidFile("Unsupported pipeline resource URL; expected file: URL: " + url);
    }
    std::string path = url.substr(colon + 1);
    if (path.find_first_of("?#") != std::string::npos) { invalidFile("File resource URLs must escape literal '?' and '#': " + url); }
    // file:///path has an empty authority; retain //host/path for UNC resources.
    if (path.starts_with("///")) { path.erase(0, 2); }
    std::string decoded;
    for (Size i = 0; i < path.size(); ++i)
    {
      if (path[i] != '%')
      {
        decoded.push_back(path[i]);
        continue;
      }
      if (i + 2 >= path.size()) { invalidFile("Malformed percent escape in resource URL: " + url); }
      unsigned int byte = 0;
      const auto parsed = std::from_chars(path.data() + i + 1, path.data() + i + 3, byte, 16);
      if (parsed.ec != std::errc() || parsed.ptr != path.data() + i + 3 || byte == 0)
      {
        invalidFile("Malformed percent escape in resource URL: " + url);
      }
      decoded.push_back(static_cast<char>(byte));
      i += 2;
    }
#ifdef OPENMS_WINDOWSPLATFORM
    // Standard file:///C:/... URLs name local drive paths on Windows.
    if (decoded.size() >= 3 && decoded[0] == '/' && decoded[2] == ':') { decoded.erase(0, 1); }
#endif
    if (decoded.empty()) { invalidFile("Empty file resource URL."); }
    // Do not normalize away a UNC authority (//server/share) on POSIX hosts.
    // As with QUrl::toLocalFile, relative file: URLs remain relative to the caller.
    return decoded;
  }
} // namespace

void PipelineFile::load(const std::string& filename, PipelineGraph& graph) const
{
  Param parameters;
  ParamXMLFile().load(filename, parameters);
  loadParam(parameters, graph, filename);
}

void PipelineFile::store(const std::string& filename, const PipelineGraph& graph) const
{ ParamXMLFile().store(filename, storeParam(graph, filename)); }

void PipelineFile::loadParam(const Param& parameters, PipelineGraph& graph, const std::string& filename) const
{
  PipelineGraph loaded;
  loaded.filename = filename;
  loaded.legacy_port_indices = ! parameters.exists("info:version");
  if (! loaded.legacy_port_indices) { loaded.version = parameters.getValue("info:version").toString(); }
  if (parameters.exists("info:description"))
  {
    loaded.description = parameters.getValue("info:description").toString();
    if (loaded.description.starts_with("<![CDATA[") && loaded.description.ends_with("]]>"))
    {
      loaded.description = loaded.description.substr(9, loaded.description.size() - 12);
    }
  }
  const Size vertex_count = indexValue(requiredValue(parameters, "info:num_vertices").toString(), "info:num_vertices");
  const Size edge_count = indexValue(requiredValue(parameters, "info:num_edges").toString(), "info:num_edges");
  const Param vertices = parameters.copy("vertices:", true);
  const auto vertex_ids = sectionNames(vertices);
  if (vertex_ids.size() != vertex_count) { invalidFile("TOPPAS vertex count does not match its vertex records."); }
  const auto base = filename.empty() ? std::filesystem::path() : std::filesystem::absolute(utf8Path(filename)).parent_path();
  for (const auto& key : vertex_ids)
  {
    PipelineGraph::Node vertex;
    vertex.id = indexValue(key, "vertex ID");
    if (vertex.id >= vertex_count) { invalidFile("TOPPAS vertex ID exceeds info:num_vertices: " + key); }
    const std::string prefix = key + ":";
    vertex.kind = nodeKind(requiredValue(vertices, prefix + "toppas_type").toString());
    if (vertices.exists(prefix + "x_pos")) { vertex.x = vertices.getValue(prefix + "x_pos"); }
    if (vertices.exists(prefix + "y_pos")) { vertex.y = vertices.getValue(prefix + "y_pos"); }
    vertex.recycle_output = boolValue(vertices, prefix + "recycle_output", false);
    if (vertex.kind == PipelineGraph::Kind::INPUT)
    {
      vertex.files = static_cast<std::vector<std::string>>(requiredValue(vertices, prefix + "file_names"));
      for (auto& file : vertex.files)
      {
        auto path = utf8Path(file);
        if (path.is_relative() && ! base.empty()) { path = base / path; }
        // A lexical cleanup of symlink/../file can silently select a different input.
        file = pathString(path);
      }
      if (vertices.exists(prefix + "resource_key")) { vertex.resource_key = vertices.getValue(prefix + "resource_key").toString(); }
    }
    else if (vertex.kind == PipelineGraph::Kind::TOOL)
    {
      vertex.tool_name = requiredValue(vertices, prefix + "tool_name").toString();
      if (vertices.exists(prefix + "tool_type")) { vertex.tool_type = vertices.getValue(prefix + "tool_type").toString(); }
      vertex.parameters = vertices.copy(prefix + "parameters:", true);
    }
    else if (vertex.kind == PipelineGraph::Kind::MERGER) { vertex.round_based = boolValue(vertices, prefix + "round_based", true); }
    else if (vertex.kind == PipelineGraph::Kind::OUTPUT || vertex.kind == PipelineGraph::Kind::OUTPUT_DIRECTORY)
    {
      if (vertices.exists(prefix + "output_folder_name")) { vertex.output_folder = vertices.getValue(prefix + "output_folder_name").toString(); }
    }
    loaded.nodes.push_back(std::move(vertex));
  }

  const Param edges = parameters.copy("edges:", true);
  const auto edge_ids = sectionNames(edges);
  if (edge_ids.size() != edge_count) { invalidFile("TOPPAS edge count does not match its edge records."); }
  for (const auto& key : edge_ids)
  {
    indexValue(key, "edge ID");
    const std::string prefix = key + ":";
    const std::string endpoints = requiredValue(edges, prefix + "source/target:").toString();
    const auto slash = endpoints.find('/');
    if (slash == std::string::npos) { invalidFile("Invalid TOPPAS edge endpoints: " + endpoints); }
    PipelineGraph::Edge edge;
    edge.source = indexValue(endpoints.substr(0, slash), "edge source");
    edge.target = indexValue(endpoints.substr(slash + 1), "edge target");
    edge.source_port = portName(requiredValue(edges, prefix + "source_out_param:").toString(), loaded.legacy_port_indices);
    edge.target_port = portName(requiredValue(edges, prefix + "target_in_param:").toString(), loaded.legacy_port_indices);
    loaded.edges.push_back(std::move(edge));
  }
  loaded.assignTopologicalNumbers();
  loaded.validate();
  graph = std::move(loaded);
}

Param PipelineFile::storeParam(const PipelineGraph& graph, const std::string& filename) const
{
  PipelineGraph numbered = graph;
  numbered.assignTopologicalNumbers();
  numbered.validate();
  const auto order = numbered.topologicalOrder();
  std::map<Size, Size> persisted_ids;
  Param parameters;
  if (! numbered.legacy_port_indices) { parameters.setValue("info:version", std::string(VersionInfo::getVersion())); }
  parameters.setValue("info:num_vertices", static_cast<int>(numbered.nodes.size()));
  parameters.setValue("info:num_edges", static_cast<int>(numbered.edges.size()));
  parameters.setValue("info:description", "<![CDATA[" + numbered.description + "]]>");
  const auto base = filename.empty() ? std::filesystem::path() : std::filesystem::absolute(utf8Path(filename)).parent_path();
  for (Size i = 0; i < order.size(); ++i)
  {
    const auto& vertex = numbered.node(order[i]);
    persisted_ids.emplace(vertex.id, i);
    const std::string prefix = "vertices:" + StringUtils::toStr(i) + ":";
    parameters.setValue(prefix + "toppas_type", nodeKindName(vertex.kind));
    parameters.setValue(prefix + "x_pos", vertex.x);
    parameters.setValue(prefix + "y_pos", vertex.y);
    parameters.setValue(prefix + "recycle_output", vertex.recycle_output ? "true" : "false");
    switch (vertex.kind)
    {
      case PipelineGraph::Kind::INPUT: {
        auto files = vertex.files;
        if (! base.empty())
        {
          for (auto& file : files)
          {
            const auto path = std::filesystem::absolute(utf8Path(file));
            const auto relative = path.lexically_relative(base);
            file = pathString(path);
            if (! relative.empty())
            {
              // Going up from a symlinked workflow directory may resolve somewhere
              // different from its lexical parent. Retain the absolute path then.
              std::error_code source_error, relative_error;
              const auto source = std::filesystem::weakly_canonical(path, source_error);
              const auto reloaded = std::filesystem::weakly_canonical(base / relative, relative_error);
              if (! source_error && ! relative_error && source == reloaded) { file = pathString(relative); }
            }
          }
        }
        parameters.setValue(prefix + "file_names", files);
        parameters.setValue(prefix + "resource_key", vertex.resource_key);
        break;
      }
      case PipelineGraph::Kind::TOOL:
        parameters.setValue(prefix + "tool_name", vertex.tool_name);
        parameters.setValue(prefix + "tool_type", vertex.tool_type);
        parameters.insert(prefix + "parameters:", vertex.parameters);
        break;
      case PipelineGraph::Kind::MERGER:
        parameters.setValue(prefix + "round_based", vertex.round_based ? "true" : "false");
        break;
      case PipelineGraph::Kind::OUTPUT:
      case PipelineGraph::Kind::OUTPUT_DIRECTORY:
        parameters.setValue(prefix + "output_folder_name", vertex.output_folder);
        break;
      case PipelineGraph::Kind::SPLITTER:
        break;
    }
  }
  for (Size i = 0; i < numbered.edges.size(); ++i)
  {
    const auto& edge = numbered.edges[i];
    if (numbered.legacy_port_indices)
    {
      portName(edge.source_port, true);
      portName(edge.target_port, true);
    }
    const std::string prefix = "edges:" + StringUtils::toStr(i) + ":";
    parameters.setValue(prefix + "source/target:",
                        StringUtils::toStr(persisted_ids.at(edge.source)) + "/" + StringUtils::toStr(persisted_ids.at(edge.target)));
    const std::string unnamed = numbered.legacy_port_indices ? "-1" : "__no_name__";
    parameters.setValue(prefix + "source_out_param:", edge.source_port.empty() ? unnamed : edge.source_port);
    parameters.setValue(prefix + "target_in_param:", edge.target_port.empty() ? unnamed : edge.target_port);
  }
  return parameters;
}

void PipelineFile::loadResources(const std::string& filename, PipelineGraph& graph) const
{
  Param parameters;
  ParamXMLFile().load(filename, parameters);
  loadResourceParam(parameters, graph);
}

void PipelineFile::loadResourceParam(const Param& parameters, PipelineGraph& graph) const
{
  std::map<std::string, std::vector<std::string>> resources;
  for (auto it = parameters.begin(); it != parameters.end(); ++it)
  {
    const std::string name = it.getName();
    const auto colon = name.find(':');
    if (colon == std::string::npos || colon == 0 || name.substr(colon) != ":url_list" || it->value.valueType() != ParamValue::STRING_LIST)
    {
      invalidFile("Invalid resource entry '" + name + "'; expected <input key>:url_list.");
    }
    std::vector<std::string> files;
    for (const auto& url : static_cast<std::vector<std::string>>(it->value))
    {
      files.push_back(decodeFileURL(url));
    }
    resources.emplace(name.substr(0, colon), std::move(files));
  }
  // Parse everything before modifying the graph: a bad late entry cannot leave partial bindings.
  PipelineGraph loaded = graph;
  loaded.assignTopologicalNumbers();
  for (auto& vertex : loaded.nodes)
  {
    if (vertex.kind != PipelineGraph::Kind::INPUT) { continue; }
    const auto resource = resources.find(vertex.resource_key);
    vertex.files = resource == resources.end() ? std::vector<std::string>() : resource->second;
  }
  graph = std::move(loaded);
}
} // namespace OpenMS
