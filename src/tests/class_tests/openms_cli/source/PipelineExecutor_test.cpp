// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineExecutor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

using namespace OpenMS;
namespace fs = std::filesystem;

namespace
{
using Kind = PipelineGraph::Kind;
using State = PipelineExecutor::State;

std::string readFile(const fs::path& filename)
{
  std::ifstream stream(filename, std::ios::binary);
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

std::string readFile(const std::string& filename)
{ return readFile(to_path(filename)); }

std::string pathString(const fs::path& path)
{
  const auto utf8 = path.generic_u8string();
  return {utf8.begin(), utf8.end()};
}

// Each section has an independently owned workspace, including paths with spaces.
struct Fixture
{
  fs::path directory = fs::temp_directory_path() / (File::getUniqueName() + " pipeline paths");
  PipelineExecutor::Options options;

  Fixture()
  {
    fs::create_directories(directory / "temporary");
    fs::create_directories(directory / "results");
    options.output_directory = pathString(directory / "results");
    options.temp_directory = pathString(directory / "temporary");
    options.keep_temporary_files = true;
  }

  ~Fixture()
  {
    std::error_code error;
    fs::remove_all(directory, error);
  }

  std::string input(const std::string& name, const std::string& sequence = "PEPTIDE") const
  {
    const auto filename = directory / to_path(name + ".fasta");
    std::ofstream stream(filename, std::ios::binary);
    stream << ">" << name << "\n" << sequence << "\n";
    return pathString(filename);
  }

  std::vector<fs::path> writtenFiles() const
  {
    std::vector<fs::path> result;
    const auto published = to_path(options.output_directory) / "TOPPAS_out";
    if (! fs::is_directory(published)) return result;
    for (const auto& item : fs::recursive_directory_iterator(published))
    {
      if (item.is_regular_file()) result.push_back(item.path());
    }
    std::sort(result.begin(), result.end());
    return result;
  }
};

PipelineGraph::Node node(Size id, Kind kind)
{
  PipelineGraph::Node result;
  result.id = id;
  result.kind = kind;
  if (kind == Kind::TOOL) result.tool_name = "PipelineTestTool";
  if (kind == Kind::OUTPUT || kind == Kind::OUTPUT_DIRECTORY) result.output_folder = "collected";
  return result;
}

PipelineGraph singleToolGraph(const std::vector<std::string>& inputs)
{
  PipelineGraph graph;
  graph.nodes = {node(10, Kind::INPUT), node(20, Kind::TOOL), node(30, Kind::OUTPUT)};
  graph.node(10).files = inputs;
  graph.edges = {{10, 20, "", "in"}, {20, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  return graph;
}

std::vector<std::string> toolContents(const PipelineExecutor::Result& result, Size id = 20, const std::string& port = "out")
{
  std::vector<std::string> contents;
  for (const auto& round : result.nodes.at(id).outputs)
  {
    for (const auto& filename : round.at(port))
      contents.push_back(readFile(to_path(filename)));
  }
  return contents;
}

PipelineGraph mergerGraph(const std::vector<std::string>& first, const std::vector<std::string>& second)
{
  auto graph = singleToolGraph({});
  graph.nodes.insert(graph.nodes.begin() + 1, node(11, Kind::INPUT));
  graph.nodes.insert(graph.nodes.begin() + 2, node(15, Kind::MERGER));
  graph.node(10).files = first;
  graph.node(11).files = second;
  graph.edges = {{10, 15, "", ""}, {11, 15, "", ""}, {15, 20, "", "in"}, {20, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  return graph;
}
} // namespace

START_TEST(PipelineExecutor, "$Id$")

START_SECTION([EXTRA] Executes real TOPP subprocesses with named ports and preserves round order)
{
  Fixture fixture;
  const auto first = fixture.input("first sample", "PEPTIDE");
  const auto second = fixture.input("second sample \xE4\xB8\xAD\xE6\x96\x87 \xC3\xA4", "SEQUENCE");
  auto graph = singleToolGraph({first, second});
  const auto caller = std::this_thread::get_id();
  bool callbacks_on_caller = true;
  Size completed = 0;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options, [&](const PipelineExecutor::Event& event) {
    callbacks_on_caller = callbacks_on_caller && std::this_thread::get_id() == caller;
    if (event.type == PipelineExecutor::Event::Type::ROUND_COMPLETED && event.node_id == 20) ++completed;
  });
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  TEST_TRUE(callbacks_on_caller)
  TEST_EQUAL(completed, 2)
  TEST_TRUE(result.nodes.at(20).state == State::SUCCEEDED)
  TEST_TRUE(result.nodes.at(30).state == State::SUCCEEDED)
  TEST_EQUAL(result.nodes.at(20).completed_rounds, 2)
  const auto contents = toolContents(result);
  TEST_EQUAL(contents.size(), 2)
  TEST_EQUAL(contents.at(0), readFile(first))
  TEST_EQUAL(contents.at(1), readFile(second))
  const auto written = fixture.writtenFiles();
  TEST_EQUAL(written.size(), 2)
  for (const auto& file : written)
    TEST_EQUAL(file.extension().string(), ".fasta")
}
END_SECTION

START_SECTION([EXTRA] Round merger preserves edge order and repeats recycled inputs)
{
  Fixture fixture;
  const auto a = fixture.input("a");
  const auto b = fixture.input("b");
  const auto c = fixture.input("c");
  const auto d = fixture.input("d");
  const auto x = fixture.input("x");
  const auto y = fixture.input("y");
  auto graph = mergerGraph({a, b, c, d}, {x, y});
  graph.node(11).recycle_output = true;
  fixture.options.num_jobs = 2;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  const auto contents = toolContents(result);
  TEST_EQUAL(contents.size(), 4)
  TEST_EQUAL(contents.at(0), readFile(x) + readFile(a))
  TEST_EQUAL(contents.at(1), readFile(y) + readFile(b))
  TEST_EQUAL(contents.at(2), readFile(x) + readFile(c))
  TEST_EQUAL(contents.at(3), readFile(y) + readFile(d))
}
END_SECTION

START_SECTION([EXTRA] Content based output naming is propagated to downstream tools)
{
  Fixture fixture;
  const auto input = fixture.input("input");
  auto graph = singleToolGraph({input});
  graph.nodes.insert(graph.nodes.end() - 1, node(40, Kind::TOOL));
  graph.edges = {{10, 20, "", "in"}, {20, 40, "untyped_out", "single_in"}, {40, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  TEST_EQUAL(toolContents(result, 40).at(0), readFile(input))
  const auto& renamed = result.nodes.at(20).outputs.at(0).at("untyped_out").at(0);
  TEST_EQUAL(to_path(renamed).extension().string(), ".fasta")
  TEST_EQUAL(readFile(renamed), readFile(input))
}
END_SECTION

START_SECTION([EXTRA] Collect and splitter convert between lists and processing rounds)
{
  Fixture fixture;
  const auto a = fixture.input("a");
  const auto b = fixture.input("b");
  auto graph = singleToolGraph({a, b});
  graph.nodes.insert(graph.nodes.begin() + 1, node(15, Kind::MERGER));
  graph.node(15).round_based = false;
  graph.edges = {{10, 15, "", ""}, {15, 20, "", "in"}, {20, 30, "out_list", ""}};
  graph.assignTopologicalNumbers();
  PipelineExecutor collector;
  const auto collected = collector.run(graph, fixture.options);
  TEST_EQUAL(collected.exit_code, 0)
  TEST_EQUAL(collected.error_message, "")
  TEST_EQUAL(collected.nodes.at(20).completed_rounds, 1)
  const auto contents = toolContents(collected, 20, "out_list");
  TEST_EQUAL(contents.size(), 2)
  for (const auto& value : contents)
    TEST_EQUAL(value, readFile(a) + readFile(b))

  graph.nodes.insert(graph.nodes.begin() + 2, node(16, Kind::SPLITTER));
  graph.edges = {{10, 15, "", ""}, {15, 16, "", ""}, {16, 20, "", "single_in"}, {20, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  fixture.options.output_directory = pathString(fixture.directory / "split results");
  PipelineExecutor splitter;
  const auto split = splitter.run(graph, fixture.options);
  TEST_EQUAL(split.exit_code, 0)
  TEST_EQUAL(split.error_message, "")
  TEST_EQUAL(split.nodes.at(20).completed_rounds, 2)
  const auto split_contents = toolContents(split);
  TEST_EQUAL(split_contents.at(0), readFile(a))
  TEST_EQUAL(split_contents.at(1), readFile(b))
}
END_SECTION

START_SECTION([EXTRA] Smart filenames distinguish equal basenames from future upstream output folders)
{
  Fixture fixture;
  const auto input = fixture.input("sample");
  PipelineGraph graph;
  graph.nodes
    = {node(10, Kind::INPUT), node(20, Kind::TOOL), node(15, Kind::MERGER), node(16, Kind::SPLITTER), node(40, Kind::TOOL), node(30, Kind::OUTPUT)};
  graph.node(10).files = {input};
  graph.node(15).round_based = false;
  graph.edges
    = {{10, 20, "", "in"}, {20, 15, "out", ""}, {20, 15, "untyped_out", ""}, {15, 16, "", ""}, {16, 40, "", "single_in"}, {40, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  const auto contents = toolContents(result, 40);
  TEST_EQUAL(contents.size(), 2)
  TEST_EQUAL(contents.at(0), readFile(input))
  TEST_EQUAL(contents.at(1), readFile(input))
  const auto files = fixture.writtenFiles();
  TEST_EQUAL(files.size(), 2)
  std::vector<std::string> names;
  for (const auto& file : files)
    names.push_back(file.filename().string());
  std::sort(names.begin(), names.end());
  TEST_EQUAL(names.at(0), "out.fasta")
  TEST_EQUAL(names.at(1), "untyped_out.fasta")
}
END_SECTION

START_SECTION([EXTRA] Output alias checks follow the destination filesystem)
{
  Fixture fixture;
  fs::create_directory(fixture.directory / "CaseProbe");
  const bool case_sensitive = fs::create_directory(fixture.directory / "caseprobe");
  fs::create_directories(fixture.directory / "upper");
  fs::create_directories(fixture.directory / "lower");
  const auto graph = singleToolGraph({fixture.input("upper/Sample"), fixture.input("lower/sample")});
  Size started = 0;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options, [&](const PipelineExecutor::Event& event) {
    if (event.type == PipelineExecutor::Event::Type::NODE_STARTED) ++started;
  });
  if (case_sensitive)
  {
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(fixture.writtenFiles().size(), 2)
  }
  else
  {
    TEST_NOT_EQUAL(result.exit_code, 0)
    TEST_FALSE(result.error_message.empty())
    TEST_EQUAL(started, 0)
    TEST_TRUE(fixture.writtenFiles().empty())
  }
}
END_SECTION

START_SECTION([EXTRA] Independent branches cannot overwrite each others published outputs)
{
  Fixture fixture;
  fs::create_directories(fixture.directory / "first");
  fs::create_directories(fixture.directory / "second");
  auto graph = singleToolGraph({fixture.input("first/sample", "FIRST")});
  graph.nodes.push_back(node(11, Kind::INPUT));
  graph.nodes.push_back(node(21, Kind::TOOL));
  graph.nodes.push_back(node(31, Kind::OUTPUT));
  graph.node(11).files = {fixture.input("second/sample", "SECOND")};
  graph.edges.push_back({11, 21, "", "in"});
  graph.edges.push_back({21, 31, "out", ""});
  Size started = 0;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options, [&](const PipelineExecutor::Event& event) {
    if (event.type == PipelineExecutor::Event::Type::NODE_STARTED) ++started;
  });
  TEST_NOT_EQUAL(result.exit_code, 0)
  TEST_TRUE(result.error_message.find("collide") != std::string::npos)
  TEST_EQUAL(started, 0)
  TEST_TRUE(fixture.writtenFiles().empty())

  // Untyped output names become known only after execution. The second branch
  // must still fail without replacing the first branch's completed publication.
  graph.edges[1].source_port = "untyped_out";
  graph.edges.back().source_port = "untyped_out";
  PipelineExecutor runtime;
  const auto runtime_result = runtime.run(graph, fixture.options);
  TEST_NOT_EQUAL(runtime_result.exit_code, 0)
  TEST_TRUE(runtime_result.error_message.find("collide") != std::string::npos)
  const auto published = fixture.writtenFiles();
  TEST_EQUAL(published.size(), 1)
  TEST_EQUAL(readFile(published.at(0)), readFile(graph.node(10).files.at(0)))
}
END_SECTION

START_SECTION([EXTRA] Output symlinks cannot redirect publication outside TOPPAS_out)
{
  Fixture fixture;
  const auto outside = fixture.directory / "outside";
  fs::create_directories(outside);
  fs::create_directories(to_path(fixture.options.output_directory) / "TOPPAS_out");
  std::error_code error;
  fs::create_directory_symlink(outside, to_path(fixture.options.output_directory) / "TOPPAS_out" / "collected", error);
  // Creating symlinks can require privileges on Windows.
  if (! error)
  {
    PipelineExecutor executor;
    const auto result = executor.run(singleToolGraph({fixture.input("sample")}), fixture.options);
    TEST_NOT_EQUAL(result.exit_code, 0)
    TEST_TRUE(result.error_message.find("outside TOPPAS_out") != std::string::npos)
    TEST_TRUE(fs::is_empty(outside))
  }
}
END_SECTION

START_SECTION([EXTRA] A symlinked temporary parent preserves valid output paths)
{
  Fixture fixture;
  const auto alias = fixture.directory / "temp alias";
  std::error_code error;
  fs::create_directory_symlink(to_path(fixture.options.temp_directory), alias, error);
  if (! error)
  {
    fixture.options.temp_directory = pathString(alias);
    const auto input = fixture.input("sample");
    PipelineExecutor executor;
    const auto result = executor.run(singleToolGraph({input}), fixture.options);
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(result.error_message, "")
    TEST_EQUAL(fixture.writtenFiles().size(), 1)
    TEST_EQUAL(readFile(fixture.writtenFiles().at(0)), readFile(input))
  }
}
END_SECTION

#ifndef _WIN32
START_SECTION([EXTRA] Long valid basenames survive content renaming and publication)
{
  for (const std::string port : {"out", "untyped_out"})
  {
    Fixture fixture;
    const auto input = fixture.input(std::string(230, 'x'));
    auto graph = singleToolGraph({input});
    graph.edges.back().source_port = port;
    PipelineExecutor executor;
    const auto result = executor.run(graph, fixture.options);
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(result.error_message, "")
    const auto published = fixture.writtenFiles();
    TEST_EQUAL(published.size(), 1)
    TEST_EQUAL(pathString(published.at(0).filename()), pathString(to_path(input).filename()))
    TEST_EQUAL(readFile(published.at(0)), readFile(input))
  }
}
END_SECTION
#endif

START_SECTION([EXTRA] Relative input names beginning with a dash are not child options)
{
  Fixture fixture;
  const auto input = fixture.input("-sample");
  struct WorkingDirectoryGuard
  {
    fs::path saved = fs::current_path();
    ~WorkingDirectoryGuard()
    { fs::current_path(saved); }
  } guard;
  fs::current_path(fixture.directory);
  PipelineExecutor executor;
  const auto result = executor.run(singleToolGraph({"-sample.fasta"}), fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  TEST_EQUAL(toolContents(result).at(0), readFile(input))
}
END_SECTION

START_SECTION([EXTRA] Rejects incompatible round counts and invalid recycling)
{
  Fixture fixture;
  const auto a = fixture.input("a");
  const auto b = fixture.input("b");
  const auto c = fixture.input("c");
  auto graph = mergerGraph({a, b, c}, {a, b});
  PipelineExecutor unequal;
  TEST_NOT_EQUAL(unequal.run(graph, fixture.options).exit_code, 0)
  graph.node(11).recycle_output = true;
  PipelineExecutor nondivisible;
  TEST_NOT_EQUAL(nondivisible.run(graph, fixture.options).exit_code, 0)
  graph.node(10).recycle_output = true;
  PipelineExecutor all_recycled;
  TEST_NOT_EQUAL(all_recycled.run(graph, fixture.options).exit_code, 0)
  TEST_TRUE(fixture.writtenFiles().empty())
}
END_SECTION

START_SECTION([EXTRA] Large file lists and nested external tool ports are passed through INI files)
{
  // Twelve files exercise the long-list INI branch independently of the ETool prefix.
  for (const std::string port : {"in", "ETool:in"})
  {
    Fixture fixture;
    std::vector<std::string> files;
    std::string expected;
    for (Size i = 0; i < 12; ++i)
    {
      files.push_back(fixture.input("sample \xC3\xBC " + std::to_string(i)));
      expected += readFile(files.back());
    }
    auto graph = singleToolGraph(files);
    graph.nodes.insert(graph.nodes.begin() + 1, node(15, Kind::MERGER));
    graph.node(15).round_based = false;
    graph.edges = {{10, 15, "", ""}, {15, 20, "", port}, {20, 30, "out_list", ""}};
    graph.assignTopologicalNumbers();
    PipelineExecutor executor;
    const auto result = executor.run(graph, fixture.options);
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(result.error_message, "")
    const auto contents = toolContents(result, 20, "out_list");
    TEST_EQUAL(contents.size(), 12)
    for (const auto& value : contents)
      TEST_EQUAL(value, expected)
    TEST_EQUAL(fixture.writtenFiles().size(), 12)
  }
}
END_SECTION

START_SECTION([EXTRA] Numeric child failures are preserved including exit code 14)
{
  for (const int exit_code : {9, 14})
  {
    Fixture fixture;
    auto graph = singleToolGraph({fixture.input("input")});
    graph.node(20).parameters.setValue("mode", "fail");
    graph.node(20).parameters.setValue("exit_code", exit_code);
    PipelineExecutor executor;
    const auto result = executor.run(graph, fixture.options);
    TEST_EQUAL(result.exit_code, exit_code)
    TEST_TRUE(result.nodes.at(20).state == State::FAILED)
    TEST_FALSE(result.nodes.at(30).state == State::SUCCEEDED)
    TEST_FALSE(result.error_message.empty())
    TEST_TRUE(readFile(to_path(fixture.options.output_directory) / "TOPPAS.log").find(result.error_message) != std::string::npos)
    TEST_TRUE(fixture.writtenFiles().empty())
  }
}
END_SECTION

START_SECTION([EXTRA] A successful process must still produce its connected outputs)
{
  Fixture fixture;
  auto graph = singleToolGraph({fixture.input("input")});
  graph.node(20).parameters.setValue("mode", "skip_output");
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, TOPPBase::UNEXPECTED_RESULT)
  TEST_TRUE(result.nodes.at(20).state == State::FAILED)
  TEST_FALSE(result.error_message.empty())
  TEST_TRUE(fixture.writtenFiles().empty())
}
END_SECTION

START_SECTION([EXTRA] Missing executables and malformed graphs fail before executing jobs)
{
  Fixture fixture;
  auto graph = singleToolGraph({fixture.input("input")});
  graph.node(20).tool_name = "PipelineTestTool_missing_" + File::getUniqueName();
  PipelineExecutor missing;
  const auto missing_result = missing.run(graph, fixture.options);
  TEST_EQUAL(missing_result.exit_code, TOPPBase::EXTERNAL_PROGRAM_NOTFOUND)
  TEST_FALSE(missing_result.error_message.empty())

  graph.node(20).tool_name = "PipelineTestTool";
  graph.nodes.push_back(node(40, Kind::TOOL));
  graph.edges.push_back({20, 40, "out", "in"});
  graph.edges.push_back({40, 20, "out", "single_in"});
  PipelineExecutor cyclic;
  const auto cyclic_result = cyclic.run(graph, fixture.options);
  TEST_NOT_EQUAL(cyclic_result.exit_code, 0)
  TEST_FALSE(cyclic_result.error_message.empty())
  TEST_TRUE(fixture.writtenFiles().empty())
}
END_SECTION

START_SECTION([EXTRA] Parallel jobs respect their configured concurrency bound)
{
  for (const Size jobs : {Size(1), Size(2)})
  {
    Fixture fixture;
    auto graph = singleToolGraph({fixture.input("a"), fixture.input("b"), fixture.input("c"), fixture.input("d")});
    const auto trace_directory = fixture.directory / "traces";
    graph.node(20).parameters.setValue("delay_ms", 200);
    graph.node(20).parameters.setValue("trace_directory", pathString(trace_directory));
    fixture.options.num_jobs = jobs;
    PipelineExecutor executor;
    const auto result = executor.run(graph, fixture.options);
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(result.error_message, "")
    std::vector<std::pair<long long, int>> events;
    for (const auto& item : fs::directory_iterator(trace_directory))
    {
      if (item.path().extension() != ".trace") continue;
      long long begin = 0;
      long long end = 0;
      std::ifstream trace(item.path());
      trace >> begin >> end;
      TEST_TRUE(begin > 0 && end >= begin)
      events.emplace_back(begin, 1);
      events.emplace_back(end, -1);
    }
    TEST_EQUAL(events.size(), 8)
    std::sort(events.begin(), events.end());
    int active = 0;
    int peak = 0;
    for (const auto& event : events)
    {
      active += event.second;
      peak = std::max(peak, active);
    }
    TEST_EQUAL(active, 0)
    TEST_TRUE(peak >= 1)
    TEST_TRUE(peak <= static_cast<int>(jobs))
    TEST_EQUAL(result.nodes.at(20).completed_rounds, 4)
  }
}
END_SECTION

START_SECTION([EXTRA] Unavailable tools in wholly disconnected editing branches are ignored)
{
  Fixture fixture;
  const auto input = fixture.input("input");
  auto graph = singleToolGraph({input});
  graph.nodes.push_back(node(40, Kind::TOOL));
  graph.nodes.push_back(node(50, Kind::OUTPUT));
  graph.node(40).tool_name = "PipelineTestTool_missing_" + File::getUniqueName();
  graph.node(50).output_folder = "disconnected";
  graph.edges.push_back({40, 50, "out", ""});
  graph.assignTopologicalNumbers();
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  TEST_TRUE(result.nodes.at(30).state == State::SUCCEEDED)
  TEST_TRUE(result.nodes.at(40).state == State::BLOCKED)
  TEST_TRUE(result.nodes.at(50).state == State::BLOCKED)
  const auto files = fixture.writtenFiles();
  TEST_EQUAL(files.size(), 1)
  TEST_EQUAL(readFile(files.at(0)), readFile(input))
}
END_SECTION

START_SECTION([EXTRA] A join with a live input and a disconnected required branch fails preflight)
{
  Fixture fixture;
  auto graph = singleToolGraph({fixture.input("input")});
  graph.nodes.push_back(node(40, Kind::TOOL));
  graph.node(40).tool_name = "PipelineTestTool_missing_" + File::getUniqueName();
  graph.edges.push_back({40, 20, "out", "single_in"});
  graph.assignTopologicalNumbers();
  Size started = 0;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options, [&](const PipelineExecutor::Event& event) {
    if (event.type == PipelineExecutor::Event::Type::NODE_STARTED) ++started;
  });
  TEST_EQUAL(result.exit_code, TOPPBase::INCOMPATIBLE_INPUT_DATA)
  TEST_TRUE(result.error_message.find("disconnected") != std::string::npos)
  TEST_EQUAL(started, 0)
  TEST_TRUE(fixture.writtenFiles().empty())
}
END_SECTION

START_SECTION([EXTRA] Output sinks copy every upstream round even when the tool permits recycling)
{
  for (const auto kind : {Kind::OUTPUT, Kind::OUTPUT_DIRECTORY})
  {
    Fixture fixture;
    const auto a = fixture.input("first");
    const auto b = fixture.input("second");
    auto graph = singleToolGraph({a, b});
    graph.node(20).recycle_output = true;
    graph.node(30).kind = kind;
    if (kind == Kind::OUTPUT_DIRECTORY) graph.edges.back().source_port = "out_dir";
    PipelineExecutor executor;
    const auto result = executor.run(graph, fixture.options);
    TEST_EQUAL(result.exit_code, 0)
    TEST_EQUAL(result.error_message, "")
    TEST_TRUE(result.nodes.at(30).state == State::SUCCEEDED)
    TEST_EQUAL(result.nodes.at(30).completed_rounds, 2)
    const auto files = fixture.writtenFiles();
    TEST_EQUAL(files.size(), 2)
    std::vector<std::string> contents;
    for (const auto& file : files)
      contents.push_back(readFile(file));
    std::sort(contents.begin(), contents.end());
    std::vector<std::string> expected {readFile(a), readFile(b)};
    std::sort(expected.begin(), expected.end());
    TEST_EQUAL(contents.at(0), expected.at(0))
    TEST_EQUAL(contents.at(1), expected.at(1))
  }
}
END_SECTION

START_SECTION([EXTRA] Cancellation from another thread terminates an active child)
{
  Fixture fixture;
  auto graph = singleToolGraph({fixture.input("input")});
  const auto trace_directory = fixture.directory / "traces";
  graph.node(20).parameters.setValue("delay_ms", 1000);
  graph.node(20).parameters.setValue("trace_directory", pathString(trace_directory));
  std::atomic_bool finished {false};
  std::atomic_bool observed_child {false};
  PipelineExecutor executor;
  std::thread cancellation([&] {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (! finished && std::chrono::steady_clock::now() < deadline)
    {
      std::error_code error;
      if (fs::is_directory(trace_directory, error) && ! fs::is_empty(trace_directory, error))
      {
        observed_child = true;
        executor.cancel();
        return;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
  });
  struct JoinGuard
  {
    std::thread& thread;
    std::atomic_bool& finished;
    ~JoinGuard()
    {
      finished = true;
      if (thread.joinable()) { thread.join(); }
    }
  } join_guard {cancellation, finished};
  const auto begin = std::chrono::steady_clock::now();
  const auto result = executor.run(graph, fixture.options);
  finished = true;
  cancellation.join();
  TEST_TRUE(observed_child.load())
  TEST_NOT_EQUAL(result.exit_code, 0)
  TEST_TRUE(result.nodes.at(20).state == State::CANCELLED)
  TEST_FALSE(result.nodes.at(30).state == State::SUCCEEDED)
  TEST_TRUE(std::chrono::steady_clock::now() - begin < std::chrono::seconds(5))
  TEST_TRUE(fixture.writtenFiles().empty())
  Size completed_children = 0;
  for (const auto& item : fs::directory_iterator(trace_directory))
  {
    if (item.path().extension() == ".trace") ++completed_children;
  }
  TEST_EQUAL(completed_children, 0)

  // Cancellation belongs to the completed run; it must not poison a later run.
  graph.node(20).parameters.setValue("delay_ms", 0);
  const auto retry = executor.run(graph, fixture.options);
  TEST_EQUAL(retry.exit_code, 0)
  TEST_EQUAL(retry.error_message, "")
  TEST_TRUE(retry.nodes.at(30).state == State::SUCCEEDED)
  TEST_EQUAL(fixture.writtenFiles().size(), 1)
}
END_SECTION

START_SECTION([EXTRA] A pre - start cancellation is consumed by one run and the executor can be reused)
{
  Fixture fixture;
  const auto input = fixture.input("input");
  const auto graph = singleToolGraph({input});
  PipelineExecutor executor;
  executor.cancel();
  Size started = 0;
  const auto cancelled = executor.run(graph, fixture.options, [&](const PipelineExecutor::Event& event) {
    if (event.type == PipelineExecutor::Event::Type::NODE_STARTED) ++started;
  });
  TEST_NOT_EQUAL(cancelled.exit_code, 0)
  TEST_EQUAL(started, 0)
  TEST_TRUE(fixture.writtenFiles().empty())

  const auto retry = executor.run(graph, fixture.options);
  TEST_EQUAL(retry.exit_code, 0)
  TEST_EQUAL(retry.error_message, "")
  TEST_EQUAL(toolContents(retry).at(0), readFile(input))
}
END_SECTION

START_SECTION([EXTRA] Rerun reuses retained upstream outputs and rejects missing cached files)
{
  Fixture fixture;
  const auto input = fixture.input("input");
  const auto original = readFile(input);
  auto graph = singleToolGraph({input});
  graph.nodes.insert(graph.nodes.end() - 1, node(40, Kind::TOOL));
  graph.edges = {{10, 20, "", "in"}, {20, 40, "out", "in"}, {40, 30, "out", ""}};
  graph.assignTopologicalNumbers();
  PipelineExecutor initial;
  const auto previous = initial.run(graph, fixture.options);
  TEST_EQUAL(previous.exit_code, 0)
  TEST_EQUAL(previous.error_message, "")
  fixture.input("input", "CHANGED");
  // Any accidental upstream execution now fails, even if output contents looked plausible.
  graph.node(20).parameters.setValue("mode", "fail");
  Size upstream_starts = 0;
  PipelineExecutor rerun;
  const auto result = rerun.run(
    graph, fixture.options,
    [&](const PipelineExecutor::Event& event) {
      if (event.type == PipelineExecutor::Event::Type::NODE_STARTED && event.node_id == 20) ++upstream_starts;
    },
    &previous, 40);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  TEST_EQUAL(upstream_starts, 0)
  TEST_EQUAL(toolContents(result, 40).at(0), original)
  fs::remove(to_path(previous.nodes.at(20).outputs.at(0).at("out").at(0)));
  Size missing_cache_starts = 0;
  PipelineExecutor missing;
  const auto missing_result = missing.run(
    graph, fixture.options,
    [&](const PipelineExecutor::Event& event) {
      if (event.type == PipelineExecutor::Event::Type::NODE_STARTED) ++missing_cache_starts;
    },
    &previous, 40);
  TEST_NOT_EQUAL(missing_result.exit_code, 0)
  TEST_FALSE(missing_result.error_message.empty())
  TEST_EQUAL(missing_cache_starts, 0)
}
END_SECTION

START_SECTION([EXTRA] Directory output copies immediate files and temporary files are cleaned)
{
  Fixture fixture;
  const auto input = fixture.input("input");
  auto graph = singleToolGraph({input});
  graph.node(30).kind = Kind::OUTPUT_DIRECTORY;
  graph.edges.back().source_port = "out_dir";
  fixture.options.keep_temporary_files = false;
  PipelineExecutor executor;
  const auto result = executor.run(graph, fixture.options);
  TEST_EQUAL(result.exit_code, 0)
  TEST_EQUAL(result.error_message, "")
  const auto files = fixture.writtenFiles();
  TEST_EQUAL(files.size(), 1)
  TEST_EQUAL(files.at(0).filename().string(), "direct.fasta")
  TEST_EQUAL(readFile(files.at(0)), readFile(input))
  TEST_TRUE(fs::is_empty(fixture.directory / "temporary"))
}
END_SECTION

END_TEST
