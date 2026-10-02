// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// $Maintainer: Timo Sachsenberg $

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineFile.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <QApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QTemporaryDir>
#include <fstream>
#include <thread>

using namespace OpenMS;

namespace
{
bool finish(TOPPASScene& scene)
{
  QElapsedTimer timer;
  timer.start();
  while (scene.isPipelineRunning() && timer.elapsed() < 10000)
  {
    QApplication::processEvents();
    std::this_thread::yield();
  }
  return ! scene.isPipelineRunning();
}
} // namespace

START_TEST(TOPPASScene, "$Id$")

// No display server is needed for the scene adapter tests.
qputenv("QT_QPA_PLATFORM", "offscreen");
QApplication application(argc, argv);
QTemporaryDir directory;
const std::string root = directory.path().toStdString();
const std::string input_file = root + "/input.txt";
std::ofstream(input_file) << "test input\n";
PipelineGraph graph;
PipelineGraph::Node input;
input.id = 4;
input.kind = PipelineGraph::Kind::INPUT;
input.resource_key = "samples";
input.files = {input_file};
input.x = 20.0;
input.y = 40.0;
PipelineGraph::Node merger;
merger.id = 17;
merger.kind = PipelineGraph::Kind::MERGER;
merger.round_based = true;
merger.recycle_output = true;
graph.nodes = {input, merger};
graph.edges = {{4, 17, "", ""}};
graph.description = "Workflow <b>description</b>";
graph.assignTopologicalNumbers();
const std::string workflow = root + "/input.toppas";
PipelineFile().store(workflow, graph);

START_SECTION((shared serialization preserves input keys, coordinates, and graph structure))
{
  TOPPASScene scene(nullptr, directory.path(), false);
  scene.load(workflow);
  scene.setDescription("Updated description");
  TEST_TRUE(scene.store(root + "/saved.toppas"))
  PipelineGraph saved;
  PipelineFile().load(root + "/saved.toppas", saved);
  TEST_EQUAL(saved.nodes.size(), 2)
  TEST_EQUAL(saved.edges.size(), 1)
  TEST_EQUAL(saved.nodes[0].resource_key, "samples")
  TEST_EQUAL(saved.nodes[0].files.front(), input_file)
  TEST_REAL_SIMILAR(saved.nodes[0].x, 20.0)
  TEST_REAL_SIMILAR(saved.nodes[0].y, 40.0)
  TEST_EQUAL(saved.description, "Updated description")
  TEST_TRUE(saved.nodes[1].recycle_output)
}
END_SECTION

START_SECTION((shared execution, layout edits, and definitions - only include))
{
  TOPPASScene scene(nullptr, directory.path(), false);
  scene.load(workflow);
  scene.setOutDir(directory.path());
  int completed = 0;
  int failed = 0;
  QObject::connect(&scene, &TOPPASScene::entirePipelineFinished, [&]() { ++completed; });
  QObject::connect(&scene, &TOPPASScene::pipelineExecutionFailed, [&](int) { ++failed; });
  scene.runPipeline();
  (*scene.verticesBegin())->setSelected(true);
  scene.moveSelectedItems(25, 10); // layout edits must not cancel the immutable run
  TEST_TRUE(finish(scene))
  TEST_EQUAL(completed, 1)
  TEST_EQUAL(failed, 0)
  for (auto it = scene.verticesBegin(); it != scene.verticesEnd(); ++it)
  {
    TEST_TRUE((*it)->isFinished())
  }

  TOPPASScene included(nullptr, directory.path(), false);
  included.include(&scene);
  auto* copied_input = qobject_cast<TOPPASInputFileListVertex*>(*included.verticesBegin());
  TEST_NOT_EQUAL(copied_input, nullptr)
  TEST_EQUAL(copied_input->getFileNames().size(), 1)
  TEST_EQUAL(fromQString(copied_input->getFileNames().front()), input_file)
  for (auto it = included.verticesBegin(); it != included.verticesEnd(); ++it)
  {
    TEST_FALSE((*it)->isFinished())
  }
  auto* copied_merger = qobject_cast<TOPPASMergerVertex*>(*(included.verticesBegin() + 1));
  TEST_NOT_EQUAL(copied_merger, nullptr)
  TEST_TRUE(copied_merger->getFileNames().empty())
}
END_SECTION

START_SECTION((abort and restart reject stale queued results; scenes are independent))
{
  TOPPASScene first(nullptr, directory.path(), false);
  TOPPASScene second(nullptr, directory.path(), false);
  first.load(workflow);
  second.load(workflow);
  first.setOutDir(toQString(root + "/first"));
  second.setOutDir(toQString(root + "/second"));
  TEST_NOT_EQUAL(first.getTempDir().toStdString(), second.getTempDir().toStdString())
  int completed = 0;
  QObject::connect(&first, &TOPPASScene::entirePipelineFinished, [&]() { ++completed; });
  first.runPipeline();
  first.abortPipeline();
  TEST_FALSE(first.isPipelineRunning())
  TEST_EQUAL(QDir(first.getTempDir()).entryList({"TOPPAS_run_*"}, QDir::Dirs | QDir::NoDotAndDotDot).size(), 0)
  first.runPipeline();
  second.runPipeline();
  TEST_TRUE(finish(first))
  TEST_TRUE(finish(second))
  TEST_EQUAL(completed, 1)
}
END_SECTION

START_SECTION((tool rerun retains upstream outputs and parameter edits cancel the old snapshot))
{
  const std::string fasta = root + "/input.fasta";
  std::ofstream(fasta) << ">test\nPEPTIDE\n";
  PipelineGraph tools;
  PipelineGraph::Node source = input;
  source.files = {fasta};
  source.id = 0;
  PipelineGraph::Node first_tool;
  first_tool.id = 1;
  first_tool.kind = PipelineGraph::Kind::TOOL;
  first_tool.tool_name = "PipelineTestTool";
  PipelineGraph::Node second_tool = first_tool;
  second_tool.id = 2;
  PipelineGraph::Node output;
  output.id = 3;
  output.kind = PipelineGraph::Kind::OUTPUT;
  tools.nodes = {source, first_tool, second_tool, output};
  tools.edges = {{0, 1, "", "in"}, {1, 2, "out", "in"}, {2, 3, "out", ""}};
  tools.assignTopologicalNumbers();
  const std::string tools_file = root + "/tools.toppas";
  PipelineFile().store(tools_file, tools);

  TOPPASScene scene(nullptr, directory.path(), false);
  scene.load(tools_file);
  scene.setOutDir(toQString(root + "/tool_results"));
  auto* first = qobject_cast<TOPPASToolVertex*>(*(scene.verticesBegin() + 1));
  auto* second = qobject_cast<TOPPASToolVertex*>(*(scene.verticesBegin() + 2));
  TEST_NOT_EQUAL(first, nullptr)
  TEST_NOT_EQUAL(second, nullptr)
  TEST_TRUE(first->getParam().exists("in"))
  TEST_TRUE(first->getParam().exists("out"))
  Param invalid_parameters;
  invalid_parameters.setValue("unknown_parameter", 1);
  TEST_EXCEPTION(Exception::InvalidParameter, first->setParam(invalid_parameters))
  TEST_TRUE(first->getParam().exists("out"))
  scene.refreshParameters(); // parameter refresh must work before the first run creates job directories
  TEST_TRUE(first->isToolReady())
  TEST_TRUE(second->isToolReady())
  int completed = 0;
  int failed = 0;
  std::string diagnostics;
  QObject::connect(&scene, &TOPPASScene::messageReady, [&](const QString& message) { diagnostics += fromQString(message) + "\n"; });
  QObject::connect(&scene, &TOPPASScene::entirePipelineFinished, [&]() { ++completed; });
  QObject::connect(&scene, &TOPPASScene::pipelineExecutionFailed, [&](int) {
    ++failed;
    std::cerr << diagnostics;
  });
  scene.runPipeline();
  TEST_TRUE(finish(scene))
  TEST_EQUAL(completed, 1)
  TEST_EQUAL(failed, 0)
  TEST_EQUAL(first->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
  TEST_EQUAL(second->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
  const auto upstream_files = first->getFileNames();
  const auto original_downstream_files = second->getFileNames();
  TEST_FALSE(upstream_files.empty())

  // Cancellation before run() initializes must retain the upstream cache too.
  scene.resumePipeline(second);
  scene.abortPipeline();
  QApplication::processEvents();
  TEST_EQUAL(first->getFileNames() == upstream_files, true)

  Param parameters = second->getParam();
  parameters.setValue("delay_ms", 5000);
  second->setParam(parameters);
  second->parameterChanged(true);
  scene.resumePipeline(second);
  QElapsedTimer timer;
  timer.start();
  while (scene.isPipelineRunning() && second->getStatus() != TOPPASToolVertex::TOOL_RUNNING && timer.elapsed() < 10000)
  {
    QApplication::processEvents();
    std::this_thread::yield();
  }
  TEST_EQUAL(second->getStatus(), TOPPASToolVertex::TOOL_RUNNING)
  parameters.setValue("delay_ms", 0);
  second->setParam(parameters);
  second->parameterChanged(true); // cancels and joins before invalidating the downstream cache
  TEST_FALSE(scene.isPipelineRunning())
  QApplication::processEvents(); // discard delayed completion from the cancelled generation
  TEST_FALSE(second->isFinished())
  TEST_EQUAL(first->getFileNames() == upstream_files, true)

  scene.resumePipeline(second);
  TEST_TRUE(finish(scene))
  TEST_EQUAL(completed, 2)
  TEST_EQUAL(failed, 0)
  TEST_EQUAL(first->getFileNames() == upstream_files, true)
  TEST_TRUE(second->getFileNames() != original_downstream_files)
  TEST_EQUAL(second->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
  for (const auto& file : upstream_files)
  {
    TEST_TRUE(File::exists(fromQString(file)))
  }

  // Each downstream rerun replaces its own temporary files while preserving
  // the original run directory that still supplies the upstream cache.
  const auto obsolete_downstream_files = second->getFileNames();
  scene.resumePipeline(second);
  TEST_TRUE(finish(scene))
  TEST_EQUAL(completed, 3)
  TEST_EQUAL(failed, 0)
  TEST_EQUAL(QDir(scene.getTempDir()).entryList({"TOPPAS_run_*"}, QDir::Dirs | QDir::NoDotAndDotDot).size(), 2)
  for (const auto& file : upstream_files)
  {
    TEST_TRUE(File::exists(fromQString(file)))
  }
  for (const auto& file : obsolete_downstream_files)
  {
    TEST_FALSE(File::exists(fromQString(file)))
  }
  const auto last_downstream_files = second->getFileNames();
  // Matching a naming pattern does not make an unrelated directory ours.
  const QString unrelated = scene.getTempDir() + "/TOPPAS_run_unowned";
  TEST_TRUE(QDir().mkpath(unrelated))
  scene.runPipeline();
  TEST_TRUE(finish(scene))
  TEST_EQUAL(completed, 4)
  TEST_EQUAL(failed, 0)
  TEST_EQUAL(QDir(scene.getTempDir()).entryList({"TOPPAS_run_*"}, QDir::Dirs | QDir::NoDotAndDotDot).size(), 2)
  TEST_TRUE(QDir(unrelated).exists())
  for (const auto& file : upstream_files)
  {
    TEST_FALSE(File::exists(fromQString(file)))
  }
  for (const auto& file : last_downstream_files)
  {
    TEST_FALSE(File::exists(fromQString(file)))
  }
  TEST_EQUAL(first->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
  TEST_EQUAL(second->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
}
END_SECTION

START_SECTION((editing input files during a run preserves the new selection))
{
  const std::string original = root + "/original.fasta";
  const std::string replacement = root + "/replacement.fasta";
  std::ofstream(original) << ">original\nPEPTIDE\n";
  std::ofstream(replacement) << ">replacement\nREPLACEMENT\n";
  PipelineGraph tools;
  PipelineGraph::Node source = input;
  source.id = 0;
  source.files = {original};
  PipelineGraph::Node tool_node;
  tool_node.id = 1;
  tool_node.kind = PipelineGraph::Kind::TOOL;
  tool_node.tool_name = "PipelineTestTool";
  tool_node.parameters.setValue("delay_ms", 5000);
  PipelineGraph::Node output;
  output.id = 2;
  output.kind = PipelineGraph::Kind::OUTPUT;
  tools.nodes = {source, tool_node, output};
  tools.edges = {{0, 1, "", "in"}, {1, 2, "out", ""}};
  const std::string tools_file = root + "/input_edit.toppas";
  PipelineFile().store(tools_file, tools);

  TOPPASScene scene(nullptr, directory.path(), false);
  scene.load(tools_file);
  scene.setOutDir(toQString(root + "/input_edit_results"));
  auto* source_vertex = qobject_cast<TOPPASInputFileListVertex*>(*scene.verticesBegin());
  auto* tool = qobject_cast<TOPPASToolVertex*>(*(scene.verticesBegin() + 1));
  TEST_NOT_EQUAL(source_vertex, nullptr)
  TEST_NOT_EQUAL(tool, nullptr)
  scene.runPipeline();
  QElapsedTimer timer;
  timer.start();
  while (scene.isPipelineRunning() && tool->getStatus() != TOPPASToolVertex::TOOL_RUNNING && timer.elapsed() < 10000)
  {
    QApplication::processEvents();
    std::this_thread::yield();
  }
  TEST_EQUAL(tool->getStatus(), TOPPASToolVertex::TOOL_RUNNING)
  // Match the input dialog: update the definition before emitting its change.
  source_vertex->setFilenames({toQString(replacement)});
  source_vertex->parameterChanged(true);
  TEST_FALSE(scene.isPipelineRunning())
  TEST_EQUAL(source_vertex->getFileNames().size(), 1)
  TEST_EQUAL(fromQString(source_vertex->getFileNames().front()), replacement)
  QApplication::processEvents();
  TEST_EQUAL(fromQString(source_vertex->getFileNames().front()), replacement)

  Param parameters = tool->getParam();
  parameters.setValue("delay_ms", 0);
  tool->setParam(parameters);
  tool->parameterChanged(true);
  scene.runPipeline();
  TEST_TRUE(finish(scene))
  TEST_EQUAL(tool->getStatus(), TOPPASToolVertex::TOOL_SUCCESS)
  TEST_EQUAL(fromQString(source_vertex->getFileNames().front()), replacement)
  const auto published = (*(scene.verticesBegin() + 2))->getFileNames();
  TEST_EQUAL(published.size(), 1)
  if (! published.empty())
  {
    std::ifstream stream(fromQString(published.front()));
    std::string first_line;
    std::getline(stream, first_line);
    TEST_EQUAL(first_line, ">replacement")
  }
}
END_SECTION

START_SECTION((destruction joins the worker and drops queued GUI callbacks))
{
  {
    TOPPASScene scene(nullptr, directory.path(), false);
    scene.load(workflow);
    scene.setOutDir(directory.path());
    scene.runPipeline();
  }
  QApplication::processEvents();
  TEST_TRUE(true)
}
END_SECTION

END_TEST
