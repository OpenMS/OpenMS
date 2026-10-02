// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineFile.h>
#include <OpenMS/APPLICATIONS/PIPELINE/PipelineGraph.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <filesystem>
#include <fstream>

using namespace OpenMS;

namespace
{
PipelineGraph exampleGraph()
{
  PipelineGraph graph;
  PipelineGraph::Node input;
  input.id = 10;
  input.files = {"sample.mzML"};
  input.resource_key = "samples";
  input.x = -30.5;
  input.y = 2.75;
  PipelineGraph::Node tool;
  tool.id = 30;
  tool.kind = PipelineGraph::Kind::TOOL;
  tool.tool_name = "FileInfo";
  tool.parameters.setValue("in", "", "Input file", {"input file"});
  tool.parameters.setValue("out", "", "Output file", {"output file"});
  tool.parameters.setValue("nested:threshold", 2.5);
  PipelineGraph::Node output;
  output.id = 50;
  output.kind = PipelineGraph::Kind::OUTPUT;
  output.output_folder = "results";
  graph.nodes = {input, tool, output};
  graph.edges = {{10, 30, "", "in"}, {30, 50, "out", ""}};
  graph.description = "Workflow <description> & details";
  return graph;
}
} // namespace

START_TEST(PipelineGraph, "$Id$")

START_SECTION((const Node& node(Size id) const))
{
  auto graph = exampleGraph();
  TEST_EQUAL(graph.node(30).tool_name, "FileInfo")
  TEST_EXCEPTION(Exception::InvalidParameter, graph.node(1))
  graph.node(10).files = {"changed.mzML"};
  const auto& const_graph = graph;
  TEST_EQUAL(const_graph.node(10).files.front(), "changed.mzML")
}
END_SECTION

START_SECTION((std::vector<Size> topologicalOrder() const))
{
  // A ready node later in a scan must run before one skipped earlier in that scan.
  PipelineGraph graph;
  for (const Size id : {1, 0, 2, 3})
  {
    PipelineGraph::Node vertex;
    vertex.id = id;
    vertex.kind = id == 0 ? PipelineGraph::Kind::INPUT : PipelineGraph::Kind::MERGER;
    graph.nodes.push_back(vertex);
  }
  graph.edges = {{0, 1, "", ""}, {0, 2, "", ""}, {1, 3, "", ""}};
  const auto order = graph.topologicalOrder();
  TEST_EQUAL(order[0], 0)
  TEST_EQUAL(order[1], 2)
  TEST_EQUAL(order[2], 1)
  TEST_EQUAL(order[3], 3)
  graph.edges.push_back({3, 1, "", ""});
  TEST_EXCEPTION(Exception::InvalidParameter, graph.topologicalOrder())
}
END_SECTION

START_SECTION((void assignTopologicalNumbers()))
{
  PipelineGraph graph;
  for (const Size id : {0, 1, 10, 11, 2, 3, 4, 5, 6, 7, 8, 9})
  {
    PipelineGraph::Node vertex;
    vertex.id = id;
    graph.nodes.push_back(vertex);
  }
  graph.node(11).resource_key = "custom input";
  graph.assignTopologicalNumbers();
  TEST_EQUAL(graph.node(10).topo_number, 3)
  TEST_EQUAL(graph.node(10).resource_key, "3")
  TEST_EQUAL(graph.node(11).resource_key, "custom input")
  TEST_EQUAL(graph.node(9).resource_key, "12")
  const Param stored = PipelineFile().storeParam(graph);
  PipelineGraph loaded;
  PipelineFile().loadParam(stored, loaded);
  TEST_EQUAL(loaded.nodes.size(), 12)
  TEST_EQUAL(loaded.node(2).resource_key, "3")
  TEST_EQUAL(loaded.node(3).resource_key, "custom input")
  TEST_EQUAL(loaded.node(11).resource_key, "12")
}
END_SECTION

START_SECTION((void validate() const))
{
  auto graph = exampleGraph();
  graph.validate();
  graph.edges.push_back({10, 30, "", "in"});
  TEST_EXCEPTION(Exception::InvalidParameter, graph.validate())
  graph = exampleGraph();
  graph.nodes.push_back(graph.nodes.front());
  TEST_EXCEPTION(Exception::InvalidParameter, graph.validate())
  graph = exampleGraph();
  graph.edges.front().source = 999;
  TEST_EXCEPTION(Exception::InvalidParameter, graph.validate())
  graph = exampleGraph();
  graph.edges.front().target = 10;
  TEST_EXCEPTION(Exception::InvalidParameter, graph.validate())
  graph = exampleGraph();
  graph.node(10).resource_key = "invalid:key";
  TEST_EXCEPTION(Exception::InvalidParameter, graph.validate())
  // Disconnected editing leftovers, including input-free tools, are valid graph definitions.
  graph = exampleGraph();
  graph.edges.clear();
  graph.validate();
}
END_SECTION

START_SECTION((PipelineFile workflow round trip and sparse identifier remapping))
{
  const auto graph = exampleGraph();
  const auto stored = PipelineFile().storeParam(graph);
  PipelineGraph loaded;
  PipelineFile().loadParam(stored, loaded);
  TEST_EQUAL(loaded.nodes.size(), 3)
  TEST_EQUAL(loaded.edges[0].source, 0)
  TEST_EQUAL(loaded.edges[0].target, 1)
  TEST_EQUAL(loaded.edges[1].target, 2)
  TEST_EQUAL(loaded.node(0).resource_key, "samples")
  TEST_REAL_SIMILAR(loaded.node(0).x, -30.5)
  TEST_REAL_SIMILAR(loaded.node(0).y, 2.75)
  TEST_EQUAL(loaded.node(1).parameters, graph.node(30).parameters)
  TEST_EQUAL(loaded.node(2).output_folder, "results")
  TEST_EQUAL(loaded.description, graph.description)
  TEST_FALSE(loaded.legacy_port_indices)
}
END_SECTION

START_SECTION((PipelineFile preserves all node kinds and flags))
{
  PipelineGraph graph;
  for (const auto kind : {PipelineGraph::Kind::INPUT, PipelineGraph::Kind::TOOL, PipelineGraph::Kind::MERGER, PipelineGraph::Kind::SPLITTER,
                          PipelineGraph::Kind::OUTPUT, PipelineGraph::Kind::OUTPUT_DIRECTORY})
  {
    PipelineGraph::Node vertex;
    vertex.id = graph.nodes.size();
    vertex.kind = kind;
    vertex.tool_name = "FileInfo";
    vertex.round_based = false;
    vertex.recycle_output = true;
    vertex.output_folder = "folder";
    graph.nodes.push_back(vertex);
  }
  PipelineGraph loaded;
  PipelineFile().loadParam(PipelineFile().storeParam(graph), loaded);
  for (Size i = 0; i < graph.nodes.size(); ++i)
  {
    TEST_EQUAL(static_cast<int>(loaded.node(i).kind), static_cast<int>(graph.node(i).kind))
    TEST_TRUE(loaded.node(i).recycle_output)
  }
  TEST_FALSE(loaded.node(2).round_based)
  TEST_EQUAL(loaded.node(5).output_folder, "folder")
}
END_SECTION

START_SECTION((PipelineFile reads edge records by field name and loads transactionally))
{
  const auto graph = exampleGraph();
  auto parameters = PipelineFile().storeParam(graph);
  parameters.removeAll("edges:");
  parameters.setValue("edges:1:target_in_param:", "__no_name__");
  parameters.setValue("edges:1:source_out_param:", "out");
  parameters.setValue("edges:1:source/target:", "1/2");
  parameters.setValue("edges:0:target_in_param:", "in");
  parameters.setValue("edges:0:source/target:", "0/1");
  parameters.setValue("edges:0:source_out_param:", "__no_name__");
  PipelineGraph loaded;
  PipelineFile().loadParam(parameters, loaded);
  TEST_EQUAL(loaded.edges[0].source_port, "out")
  TEST_EQUAL(loaded.edges[1].target_port, "in")
  parameters.setValue("edges:0:source/target:", "-1/1");
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadParam(parameters, loaded))
  TEST_EQUAL(loaded.nodes.size(), 3)
  TEST_EQUAL(loaded.edges[1].source, 0)
  parameters.setValue("edges:0:source/target:", "999/1");
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadParam(parameters, loaded))
  parameters.setValue("edges:0:source/target:", "0/1");
  parameters.setValue("vertices:1:toppas_type", "future unknown kind");
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadParam(parameters, loaded))
}
END_SECTION

START_SECTION((PipelineFile preserves legacy index ports without launching tools))
{
  auto parameters = PipelineFile().storeParam(exampleGraph());
  parameters.remove("info:version");
  parameters.setValue("edges:0:source_out_param:", "-1");
  parameters.setValue("edges:0:target_in_param:", "0");
  parameters.setValue("edges:1:source_out_param:", "1");
  parameters.setValue("edges:1:target_in_param:", "-1");
  PipelineGraph loaded;
  PipelineFile().loadParam(parameters, loaded);
  TEST_TRUE(loaded.legacy_port_indices)
  TEST_EQUAL(loaded.edges[0].target_port, "0")
  const auto stored = PipelineFile().storeParam(loaded);
  TEST_FALSE(stored.exists("info:version"))
  TEST_EQUAL(stored.getValue("edges:0:source_out_param:").toString(), "-1")
  parameters.setValue("edges:0:target_in_param:", "-2");
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadParam(parameters, loaded))
}
END_SECTION

START_SECTION((PipelineFile resolves embedded paths relative to the workflow))
{
  std::string filename;
  NEW_TMP_FILE(filename)
  auto graph = exampleGraph();
  const auto input = std::filesystem::absolute(std::filesystem::path(filename)).parent_path() / "input data.mzML";
  graph.node(10).files = {input.generic_string()};
  PipelineFile().store(filename, graph);
  Param stored;
  ParamXMLFile().load(filename, stored);
  const auto embedded = static_cast<std::vector<std::string>>(stored.getValue("vertices:0:file_names"));
  TEST_EQUAL(embedded.front(), "input data.mzML")
  PipelineGraph loaded;
  PipelineFile().load(filename, loaded);
  TEST_EQUAL(loaded.node(0).files.front(), input.generic_string())
}
END_SECTION

START_SECTION([EXTRA] Workflow input paths retain symlink and parent directory semantics)
{
  namespace fs = std::filesystem;
  std::string temporary;
  NEW_TMP_FILE(temporary)
  const auto directory = fs::absolute(to_path(temporary + "_paths"));
  struct Cleanup
  {
    fs::path directory;
    ~Cleanup()
    {
      std::error_code ignored;
      fs::remove_all(directory, ignored);
    }
  } cleanup {directory};
  fs::create_directories(directory / "actual" / "subdir");
  std::ofstream(directory / "sample.mzML") << "outer";
  std::ofstream(directory / "actual" / "sample.mzML") << "intended";
  std::error_code error;
  fs::create_directory_symlink(directory / "actual" / "subdir", directory / "shortcut", error);
  // Windows can require extra privileges to create directory symlinks.
  if (! error)
  {
    auto path_string = [](const fs::path& path) {
      const auto text = path.generic_u8string();
      return std::string(text.begin(), text.end());
    };
    auto graph = exampleGraph();
    graph.node(10).files = {"shortcut/../sample.mzML"};
    const auto workflow = path_string(directory / "workflow.toppas");
    PipelineGraph loaded;
    PipelineFile().loadParam(PipelineFile().storeParam(graph), loaded, workflow);
    TEST_TRUE(fs::equivalent(to_path(loaded.node(0).files.front()), directory / "actual" / "sample.mzML"))

    graph.node(10).files = {path_string(directory / "shortcut" / ".." / "sample.mzML")};
    PipelineFile().store(workflow, graph);
    PipelineFile().load(workflow, loaded);
    TEST_TRUE(fs::equivalent(to_path(loaded.node(0).files.front()), directory / "actual" / "sample.mzML"))

    // A workflow opened through a directory symlink cannot use a lexical ../
    // path to refer to the symlink's sibling: that would select actual/sample.
    const auto linked_workflow = path_string(directory / "shortcut" / "workflow.toppas");
    graph.node(10).files = {path_string(directory / "sample.mzML")};
    PipelineFile().store(linked_workflow, graph);
    PipelineFile().load(linked_workflow, loaded);
    TEST_TRUE(fs::equivalent(to_path(loaded.node(0).files.front()), directory / "sample.mzML"))
  }
}
END_SECTION

START_SECTION((PipelineFile replaces all resource bindings and rejects invalid URLs transactionally))
{
  auto graph = exampleGraph();
  PipelineGraph::Node other;
  other.id = 70;
  other.files = {"old.mzML"};
  other.resource_key = "missing";
  graph.nodes.push_back(other);
  Param resources;
  resources.setValue("samples:url_list", std::vector<std::string> {"file:sample%20one.mzML", "file:sample%23two.mzML"});
  PipelineFile().loadResourceParam(resources, graph);
  TEST_EQUAL(graph.node(10).files[0], "sample one.mzML")
  TEST_EQUAL(graph.node(10).files[1], "sample#two.mzML")
  TEST_TRUE(graph.node(70).files.empty())
  resources.setValue("samples:url_list", std::vector<std::string> {"file:replacement.mzML"});
  resources.setValue("extra:url_list", std::vector<std::string> {"https://example.org/not-supported.mzML"});
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadResourceParam(resources, graph))
  TEST_EQUAL(graph.node(10).files[0], "sample one.mzML")
  resources.removeAll("extra:");
  resources.setValue("samples:url_list", std::vector<std::string> {"file:bad%GG.mzML"});
  TEST_EXCEPTION(Exception::InvalidParameter, PipelineFile().loadResourceParam(resources, graph))
  resources.setValue("samples:url_list", std::vector<std::string> {"file://server/share/sample.mzML"});
  PipelineFile().loadResourceParam(resources, graph);
  TEST_EQUAL(graph.node(10).files.front(), "//server/share/sample.mzML")
}
END_SECTION

END_TEST
