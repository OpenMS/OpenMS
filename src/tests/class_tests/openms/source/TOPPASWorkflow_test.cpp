// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/TOPPASWorkflow.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

using NodeType = TOPPASWorkflow::NodeType;

// compare two workflows field-by-field (ignoring 'version', which store()/toParam() always rewrite to the current OpenMS version)
static void compareWorkflows(const TOPPASWorkflow& a, const TOPPASWorkflow& b)
{
  TEST_EQUAL(a.description, b.description)
  TEST_EQUAL(a.nodes.size(), b.nodes.size())
  TEST_EQUAL(a.edges.size(), b.edges.size())
  for (Size i = 0; i < a.nodes.size() && i < b.nodes.size(); ++i)
  {
    const TOPPASWorkflow::Node& x = a.nodes[i];
    const TOPPASWorkflow::Node& y = b.nodes[i];
    TEST_EQUAL(x.id, y.id)
    TEST_EQUAL(x.type == y.type, true)
    TEST_REAL_SIMILAR(x.x_pos, y.x_pos)
    TEST_REAL_SIMILAR(x.y_pos, y.y_pos)
    TEST_EQUAL(x.recycle_output, y.recycle_output)
    TEST_EQUAL(x.tool_name, y.tool_name)
    TEST_EQUAL(x.tool_type, y.tool_type)
    TEST_EQUAL(x.parameters == y.parameters, true)
    TEST_EQUAL(x.file_names == y.file_names, true)
    TEST_EQUAL(x.output_folder_name, y.output_folder_name)
    TEST_EQUAL(x.round_based, y.round_based)
  }
  for (Size i = 0; i < a.edges.size() && i < b.edges.size(); ++i)
  {
    const TOPPASWorkflow::Edge& x = a.edges[i];
    const TOPPASWorkflow::Edge& y = b.edges[i];
    TEST_EQUAL(x.source_id, y.source_id)
    TEST_EQUAL(x.target_id, y.target_id)
    TEST_EQUAL(x.source_out_param, y.source_out_param)
    TEST_EQUAL(x.target_in_param, y.target_in_param)
  }
}

START_TEST(TOPPASWorkflow, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOPPASWorkflow* ptr = nullptr;
TOPPASWorkflow* null_ptr = nullptr;
START_SECTION(TOPPASWorkflow())
{
  ptr = new TOPPASWorkflow();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TOPPASWorkflow())
{
  delete ptr;
}
END_SECTION

START_SECTION((static NodeType typeFromString(const std::string& s)))
{
  TEST_EQUAL(TOPPASWorkflow::typeFromString("input file list") == NodeType::INPUT_FILE_LIST, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("output file list") == NodeType::OUTPUT_FILE_LIST, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("output folder") == NodeType::OUTPUT_FOLDER, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("tool") == NodeType::TOOL, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("merger") == NodeType::MERGER, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("splitter") == NodeType::SPLITTER, true)
  TEST_EQUAL(TOPPASWorkflow::typeFromString("nonsense") == NodeType::UNKNOWN, true)
}
END_SECTION

START_SECTION((static std::string typeToString(NodeType t)))
{
  TEST_EQUAL(TOPPASWorkflow::typeToString(NodeType::INPUT_FILE_LIST), "input file list")
  TEST_EQUAL(TOPPASWorkflow::typeToString(NodeType::TOOL), "tool")
  TEST_EQUAL(TOPPASWorkflow::typeToString(NodeType::MERGER), "merger")
  // string <-> enum should round-trip for all real types
  const NodeType all[] = {NodeType::INPUT_FILE_LIST, NodeType::OUTPUT_FILE_LIST, NodeType::OUTPUT_FOLDER,
                          NodeType::TOOL, NodeType::MERGER, NodeType::SPLITTER};
  for (NodeType t : all)
  {
    TEST_EQUAL(TOPPASWorkflow::typeFromString(TOPPASWorkflow::typeToString(t)) == t, true)
  }
}
END_SECTION

START_SECTION((Param toParam() const) + (void fromParam(const Param& param)))
{
  // build a representative workflow covering all node types and named/unnamed edges
  TOPPASWorkflow wf;
  wf.description = "my <test> pipeline";

  TOPPASWorkflow::Node in;
  in.id = "0";
  in.type = NodeType::INPUT_FILE_LIST;
  in.x_pos = -140.0;
  in.y_pos = -160.0;
  in.recycle_output = false;
  in.file_names = {"../a.mzML", "../b.mzML"};
  wf.nodes.push_back(in);

  TOPPASWorkflow::Node tool;
  tool.id = "1";
  tool.type = NodeType::TOOL;
  tool.x_pos = 10.0;
  tool.y_pos = 20.0;
  tool.recycle_output = true;
  tool.tool_name = "FileFilter";
  tool.tool_type = "";
  Param tp;
  tp.setValue("in", "");
  tp.setValue("out", "");
  tp.setValue("threads", 1);
  tool.parameters = tp;
  wf.nodes.push_back(tool);

  TOPPASWorkflow::Node mg;
  mg.id = "2";
  mg.type = NodeType::MERGER;
  mg.x_pos = 5.0;
  mg.y_pos = 30.0;
  mg.round_based = false;
  wf.nodes.push_back(mg);

  TOPPASWorkflow::Node out;
  out.id = "3";
  out.type = NodeType::OUTPUT_FILE_LIST;
  out.x_pos = 0.0;
  out.y_pos = 40.0;
  out.output_folder_name = "results";
  wf.nodes.push_back(out);

  TOPPASWorkflow::Edge e1;
  e1.source_id = "0";
  e1.target_id = "1";
  e1.source_out_param = TOPPASWorkflow::NO_NAME;
  e1.target_in_param = "in";
  wf.edges.push_back(e1);

  TOPPASWorkflow::Edge e2;
  e2.source_id = "1";
  e2.target_id = "2";
  e2.source_out_param = "out";
  e2.target_in_param = TOPPASWorkflow::NO_NAME;
  wf.edges.push_back(e2);

  // in-memory round-trip via Param
  TOPPASWorkflow wf_param;
  wf_param.fromParam(wf.toParam());
  compareWorkflows(wf, wf_param);

  // the serialized Param uses the documented .toppas layout
  Param p = wf.toParam();
  TEST_EQUAL(p.exists("info:num_vertices"), true)
  TEST_EQUAL((int)p.getValue("info:num_vertices"), 4)
  TEST_EQUAL((int)p.getValue("info:num_edges"), 2)
  TEST_EQUAL(p.getValue("vertices:1:tool_name").toString(), "FileFilter")
  TEST_EQUAL(p.getValue("vertices:2:round_based").toString(), "false")
  TEST_EQUAL(p.getValue("edges:0:source/target:").toString(), "0/1")
  TEST_EQUAL(p.getValue("edges:0:target_in_param:").toString(), "in")
}
END_SECTION

START_SECTION((void store(const std::string& filename) const) + (void load(const std::string& filename)))
{
  TOPPASWorkflow wf;
  wf.description = "round trip";
  TOPPASWorkflow::Node tool;
  tool.id = "0";
  tool.type = NodeType::TOOL;
  tool.tool_name = "FeatureFinderCentroided";
  Param tp;
  tp.setValue("in", "");
  tp.setValue("out", "");
  tool.parameters = tp;
  wf.nodes.push_back(tool);

  std::string tmp;
  NEW_TMP_FILE(tmp)
  wf.store(tmp);
  TOPPASWorkflow wf_loaded;
  wf_loaded.load(tmp);
  compareWorkflows(wf, wf_loaded);
}
END_SECTION

START_SECTION([EXTRA] load/store round-trip on a real .toppas file)
{
  TOPPASWorkflow wf1;
  wf1.load(OPENMS_GET_TEST_DATA_PATH("FileHandler_toppas.toppas"));
  TEST_EQUAL(wf1.nodes.empty(), false)

  std::string tmp;
  NEW_TMP_FILE(tmp)
  wf1.store(tmp);
  TOPPASWorkflow wf2;
  wf2.load(tmp);
  compareWorkflows(wf1, wf2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
