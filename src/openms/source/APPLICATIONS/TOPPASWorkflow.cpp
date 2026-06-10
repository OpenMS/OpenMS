// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPASWorkflow.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

namespace OpenMS
{
  const std::string TOPPASWorkflow::NO_NAME = "__no_name__";

  namespace
  {
    /// Split @p s at every occurrence of @p delim (keeps empty fields; never returns an empty vector).
    std::vector<std::string> splitStr(const std::string& s, char delim)
    {
      std::vector<std::string> out;
      std::string::size_type start = 0, pos;
      while ((pos = s.find(delim, start)) != std::string::npos)
      {
        out.push_back(s.substr(start, pos - start));
        start = pos + 1;
      }
      out.push_back(s.substr(start));
      return out;
    }
  } // namespace

  TOPPASWorkflow::NodeType TOPPASWorkflow::typeFromString(const std::string& s)
  {
    if (s == "input file list") return NodeType::INPUT_FILE_LIST;
    if (s == "output file list") return NodeType::OUTPUT_FILE_LIST;
    if (s == "output folder") return NodeType::OUTPUT_FOLDER;
    if (s == "tool") return NodeType::TOOL;
    if (s == "merger") return NodeType::MERGER;
    if (s == "splitter") return NodeType::SPLITTER;
    return NodeType::UNKNOWN;
  }

  std::string TOPPASWorkflow::typeToString(NodeType t)
  {
    switch (t)
    {
      case NodeType::INPUT_FILE_LIST: return "input file list";
      case NodeType::OUTPUT_FILE_LIST: return "output file list";
      case NodeType::OUTPUT_FOLDER: return "output folder";
      case NodeType::TOOL: return "tool";
      case NodeType::MERGER: return "merger";
      case NodeType::SPLITTER: return "splitter";
      default: return "";
    }
  }

  void TOPPASWorkflow::fromParam(const Param& param)
  {
    version.clear();
    description.clear();
    nodes.clear();
    edges.clear();

    if (param.exists("info:version"))
    {
      version = param.getValue("info:version").toString();
    }
    if (param.exists("info:description"))
    {
      std::string text = param.getValue("info:description").toString();
      StringUtils::substitute(text, "<![CDATA[", "");
      StringUtils::substitute(text, "]]>", "");
      description = StringUtils::trimmed(text);
    }

    // --- nodes ---
    Param vertices_param = param.copy("vertices:", true);
    for (Param::ParamIterator it = vertices_param.begin(); it != vertices_param.end(); ++it)
    {
      std::vector<std::string> parts = splitStr(std::string(it.getName()), ':');
      if (parts.back() != "toppas_type") // every node is introduced by its 'toppas_type' entry
      {
        continue;
      }

      Node n;
      n.id = parts[0];
      n.type = typeFromString(it->value.toString());
      if (n.type == NodeType::UNKNOWN) // unknown/forward-compatible type: skip (as the old TOPPASScene did)
      {
        OPENMS_LOG_WARN << "Unknown vertex type '" << it->value.toString() << "' in TOPPAS workflow -- skipping node." << std::endl;
        continue;
      }
      const std::string base = n.id + ":";
      if (vertices_param.exists(base + "x_pos")) n.x_pos = vertices_param.getValue(base + "x_pos");
      if (vertices_param.exists(base + "y_pos")) n.y_pos = vertices_param.getValue(base + "y_pos");
      if (vertices_param.exists(base + "recycle_output"))
      {
        n.recycle_output = (vertices_param.getValue(base + "recycle_output").toString() == "true");
      }

      switch (n.type)
      {
        case NodeType::INPUT_FILE_LIST:
          if (vertices_param.exists(base + "file_names"))
          {
            n.file_names = ListUtils::toStringList<std::string>(vertices_param.getValue(base + "file_names"));
          }
          break;
        case NodeType::OUTPUT_FILE_LIST:
        case NodeType::OUTPUT_FOLDER:
          if (vertices_param.exists(base + "output_folder_name"))
          {
            n.output_folder_name = vertices_param.getValue(base + "output_folder_name").toString();
          }
          break;
        case NodeType::TOOL:
          if (vertices_param.exists(base + "tool_name")) n.tool_name = vertices_param.getValue(base + "tool_name").toString();
          if (vertices_param.exists(base + "tool_type")) n.tool_type = vertices_param.getValue(base + "tool_type").toString();
          n.parameters = vertices_param.copy(base + "parameters:", true);
          break;
        case NodeType::MERGER:
          // default to round-based if the flag is absent
          n.round_based = vertices_param.exists(base + "round_based")
                            ? (vertices_param.getValue(base + "round_based").toString() == "true")
                            : true;
          break;
        default:
          break;
      }
      nodes.push_back(n);
    }

    // --- edges ---
    // each edge contributes three consecutive entries: source/target, source_out_param, target_in_param
    Param edges_param = param.copy("edges:", true);
    for (Param::ParamIterator it = edges_param.begin(); it != edges_param.end(); ++it)
    {
      std::vector<std::string> st = splitStr(it->value.toString(), '/');
      if (st.size() != 2) // malformed source/target: stop (as the old TOPPASScene did)
      {
        break;
      }
      Edge e;
      e.source_id = st[0];
      e.target_id = st[1];
      if (++it == edges_param.end()) break;
      e.source_out_param = it->value.toString();
      if (++it == edges_param.end()) break;
      e.target_in_param = it->value.toString();
      edges.push_back(e);
    }
  }

  Param TOPPASWorkflow::toParam() const
  {
    Param param;
    param.setValue("info:version", VersionInfo::getVersion());
    param.setValue("info:num_vertices", (int)nodes.size());
    param.setValue("info:num_edges", (int)edges.size());
    param.setValue("info:description", std::string("<![CDATA[") + description + std::string("]]>"));

    for (const Node& n : nodes)
    {
      const std::string base = "vertices:" + n.id + ":";
      param.setValue(base + "toppas_type", typeToString(n.type));
      param.setValue(base + "x_pos", n.x_pos);
      param.setValue(base + "y_pos", n.y_pos);
      param.setValue(base + "recycle_output", n.recycle_output ? "true" : "false");

      switch (n.type)
      {
        case NodeType::INPUT_FILE_LIST:
          param.setValue(base + "file_names", n.file_names);
          break;
        case NodeType::OUTPUT_FILE_LIST:
        case NodeType::OUTPUT_FOLDER:
          param.setValue(base + "output_folder_name", n.output_folder_name);
          break;
        case NodeType::TOOL:
          param.setValue(base + "tool_name", n.tool_name);
          param.setValue(base + "tool_type", n.tool_type);
          param.insert(base + "parameters:", n.parameters);
          break;
        case NodeType::MERGER:
          param.setValue(base + "round_based", n.round_based ? "true" : "false");
          break;
        default:
          break;
      }
    }

    int counter = 0;
    for (const Edge& e : edges)
    {
      const std::string base = "edges:" + std::to_string(counter) + ":";
      param.setValue(base + "source/target:", e.source_id + "/" + e.target_id);
      param.setValue(base + "source_out_param:", e.source_out_param);
      param.setValue(base + "target_in_param:", e.target_in_param);
      ++counter;
    }

    return param;
  }

  void TOPPASWorkflow::load(const std::string& filename)
  {
    Param p;
    ParamXMLFile().load(filename, p);
    fromParam(p);
  }

  void TOPPASWorkflow::store(const std::string& filename) const
  {
    ParamXMLFile().store(filename, toParam());
  }

} // namespace OpenMS
