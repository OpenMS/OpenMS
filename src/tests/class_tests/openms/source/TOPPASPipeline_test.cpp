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
#include <OpenMS/APPLICATIONS/TOPPASPipeline.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

using NodeKind = TOPPASPipeline::NodeKind;
using RoundPackages = TOPPASPipeline::RoundPackages;

namespace
{
  // helper: make a node of a given kind
  TOPPASPipeline::Node mkNode(NodeKind kind, const std::string& id)
  {
    TOPPASPipeline::Node n;
    n.kind = kind;
    n.id = id;
    return n;
  }

  // helper: an input node seeded with the given files
  TOPPASPipeline::Node mkInput(const std::string& id, const std::vector<std::string>& files, bool recycle = false)
  {
    TOPPASPipeline::Node n = mkNode(NodeKind::INPUT, id);
    n.input_files = files;
    n.recycle_output = recycle;
    return n;
  }

  // collect the files of a node's output for a given round and parameter index
  std::vector<std::string> outFiles(const TOPPASPipeline& p, int node, Int param, int round)
  {
    return p.getFileNames(node, param, round);
  }
}

START_TEST(TOPPASPipeline, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TOPPASPipeline* ptr = nullptr;
TOPPASPipeline* null_ptr = nullptr;
START_SECTION(TOPPASPipeline())
{
  ptr = new TOPPASPipeline();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TOPPASPipeline())
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] input node seeding: one round per file)
{
  TOPPASPipeline p;
  int in = p.addNode(mkInput("0", {"a.mzML", "b.mzML", "c.mzML"}));
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)

  TEST_EQUAL(p.nodes()[in].round_total, 3)
  TEST_EQUAL(p.nodes()[in].output_files.size(), 3)
  TEST_EQUAL(outFiles(p, in, -1, 0) == std::vector<std::string>{"a.mzML"}, true)
  TEST_EQUAL(outFiles(p, in, -1, 1) == std::vector<std::string>{"b.mzML"}, true)
  TEST_EQUAL(outFiles(p, in, -1, 2) == std::vector<std::string>{"c.mzML"}, true)
}
END_SECTION

START_SECTION([EXTRA] linear input -> tool routing)
{
  TOPPASPipeline p;
  int in = p.addNode(mkInput("0", {"a.mzML", "b.mzML"}));
  int tool = p.addNode(mkNode(NodeKind::TOOL, "1"));
  p.addEdge(in, tool, /*src_out*/ -1, /*tgt_in*/ 0); // input -> tool's input param 0
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)

  // tool runs 2 rounds
  TEST_EQUAL(p.nodes()[tool].round_total, 2)

  // its input packages: each round, param 0 carries the matching input file
  RoundPackages pkg;
  TEST_EQUAL(p.buildRoundPackages(tool, pkg, err), true)
  TEST_EQUAL(pkg.size(), 2)
  TEST_EQUAL(pkg[0].count(0), 1)
  TEST_EQUAL(pkg[0][0].filenames == std::vector<std::string>{"a.mzML"}, true)
  TEST_EQUAL(pkg[1][0].filenames == std::vector<std::string>{"b.mzML"}, true)
}
END_SECTION

START_SECTION([EXTRA] merger round-based combines per round and keeps the round count)
{
  TOPPASPipeline p;
  int a = p.addNode(mkInput("0", {"a1", "a2"}));
  int b = p.addNode(mkInput("1", {"b1", "b2"}));
  int m = p.addNode(mkNode(NodeKind::MERGER, "2")); // round_based_merge defaults true
  p.addEdge(a, m, -1, -1);
  p.addEdge(b, m, -1, -1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)

  TEST_EQUAL(p.nodes()[m].round_total, 2)
  TEST_EQUAL(p.nodes()[m].output_files.size(), 2)
  // Each round merges the matching files of both inputs. The merger keys the first incoming edge
  // at -1 and the second at -2 (collision walk); iterating the round map in ascending key order
  // (std::map) therefore visits -2 before -1, i.e. the SECOND edge (b) before the first (a).
  // This faithfully reproduces the historical TOPPASMergerVertex behavior.
  TEST_EQUAL((outFiles(p, m, -1, 0) == std::vector<std::string>{"b1", "a1"}), true)
  TEST_EQUAL((outFiles(p, m, -1, 1) == std::vector<std::string>{"b2", "a2"}), true)
}
END_SECTION

START_SECTION([EXTRA] merger (collect): flattens everything into a single round)
{
  TOPPASPipeline p;
  int a = p.addNode(mkInput("0", {"a1", "a2"}));
  int b = p.addNode(mkInput("1", {"b1"}));
  TOPPASPipeline::Node mn = mkNode(NodeKind::MERGER, "2");
  mn.round_based_merge = false; // collect
  int m = p.addNode(mn);
  p.addEdge(a, m, -1, -1);
  p.addEdge(b, m, -1, -1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  // collect requires the non-recyclers to agree on round count -> a has 2, b has 1 -> mismatch
  TEST_EQUAL(p.computeDataFlow(err), false)
}
END_SECTION

START_SECTION([EXTRA] merger (collect) with equal round counts)
{
  TOPPASPipeline p;
  int a = p.addNode(mkInput("0", {"a1", "a2"}));
  int b = p.addNode(mkInput("1", {"b1", "b2"}));
  TOPPASPipeline::Node mn = mkNode(NodeKind::MERGER, "2");
  mn.round_based_merge = false; // collect
  int m = p.addNode(mn);
  p.addEdge(a, m, -1, -1);
  p.addEdge(b, m, -1, -1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)
  TEST_EQUAL(p.nodes()[m].round_total, 1)
  // everything collapsed into round 0; per input round the second edge (b) comes before the first (a)
  // (see ordering note above), so: round0 -> b1,a1 then round1 -> b2,a2, appended
  TEST_EQUAL((outFiles(p, m, -1, 0) == std::vector<std::string>{"b1", "a1", "b2", "a2"}), true)
}
END_SECTION

START_SECTION([EXTRA] splitter: 1 round of N files becomes N rounds of 1 file)
{
  TOPPASPipeline p;
  int in = p.addNode(mkInput("0", {"x1", "x2", "x3"}));
  TOPPASPipeline::Node cn = mkNode(NodeKind::MERGER, "1");
  cn.round_based_merge = false; // collect 3 rounds into 1 round of 3 files
  int collect = p.addNode(cn);
  int split = p.addNode(mkNode(NodeKind::SPLITTER, "2"));
  p.addEdge(in, collect, -1, -1);
  p.addEdge(collect, split, -1, -1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)

  TEST_EQUAL(p.nodes()[collect].round_total, 1)
  TEST_EQUAL((outFiles(p, collect, -1, 0) == std::vector<std::string>{"x1", "x2", "x3"}), true)

  TEST_EQUAL(p.nodes()[split].round_total, 3)
  TEST_EQUAL((outFiles(p, split, -1, 0) == std::vector<std::string>{"x1"}), true)
  TEST_EQUAL((outFiles(p, split, -1, 1) == std::vector<std::string>{"x2"}), true)
  TEST_EQUAL((outFiles(p, split, -1, 2) == std::vector<std::string>{"x3"}), true)
}
END_SECTION

START_SECTION([EXTRA] output recycling: divisible round counts wrap via modulo)
{
  TOPPASPipeline p;
  int a = p.addNode(mkInput("0", {"r1", "r2"}, /*recycle*/ true)); // 2 files, recycled
  int b = p.addNode(mkInput("1", {"f1", "f2", "f3", "f4"}));       // 4 files, sets the pace
  int tool = p.addNode(mkNode(NodeKind::TOOL, "2"));
  p.addEdge(a, tool, -1, 0); // recycled input on param 0
  p.addEdge(b, tool, -1, 1); // non-recycled input on param 1
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)
  TEST_EQUAL(p.nodes()[tool].round_total, 4)

  RoundPackages pkg;
  TEST_EQUAL(p.buildRoundPackages(tool, pkg, err), true)
  TEST_EQUAL(pkg.size(), 4)
  // param 1 = b straight through; param 0 = a wrapped modulo 2
  TEST_EQUAL((pkg[0][1].filenames == std::vector<std::string>{"f1"}), true)
  TEST_EQUAL((pkg[1][1].filenames == std::vector<std::string>{"f2"}), true)
  TEST_EQUAL((pkg[2][1].filenames == std::vector<std::string>{"f3"}), true)
  TEST_EQUAL((pkg[3][1].filenames == std::vector<std::string>{"f4"}), true)
  TEST_EQUAL((pkg[0][0].filenames == std::vector<std::string>{"r1"}), true)
  TEST_EQUAL((pkg[1][0].filenames == std::vector<std::string>{"r2"}), true)
  TEST_EQUAL((pkg[2][0].filenames == std::vector<std::string>{"r1"}), true) // wrapped
  TEST_EQUAL((pkg[3][0].filenames == std::vector<std::string>{"r2"}), true) // wrapped
}
END_SECTION

START_SECTION([EXTRA] output recycling: non-divisible round counts are an error)
{
  TOPPASPipeline p;
  int a = p.addNode(mkInput("0", {"r1", "r2", "r3"}, /*recycle*/ true)); // 3 files, recycled
  int b = p.addNode(mkInput("1", {"f1", "f2", "f3", "f4"}));             // 4 files
  int tool = p.addNode(mkNode(NodeKind::TOOL, "2"));
  p.addEdge(a, tool, -1, 0);
  p.addEdge(b, tool, -1, 1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  // 4 % 3 != 0 -> cannot recycle
  TEST_EQUAL(p.computeDataFlow(err), false)
  TEST_NOT_EQUAL(err.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] diamond: fan-out from one source and rejoin at a merger)
{
  TOPPASPipeline p;
  int in = p.addNode(mkInput("0", {"d1", "d2"}));
  int m1 = p.addNode(mkNode(NodeKind::MERGER, "1")); // pass-through (single input)
  int m2 = p.addNode(mkNode(NodeKind::MERGER, "2")); // pass-through (single input)
  int join = p.addNode(mkNode(NodeKind::MERGER, "3"));
  p.addEdge(in, m1, -1, -1);
  p.addEdge(in, m2, -1, -1);
  p.addEdge(m1, join, -1, -1);
  p.addEdge(m2, join, -1, -1);
  std::string err;
  TEST_EQUAL(p.topoSort(err), true)
  TEST_EQUAL(p.computeDataFlow(err), true)

  // both paths preserve 2 rounds; the join concatenates them per round (m1 before m2)
  TEST_EQUAL(p.nodes()[join].round_total, 2)
  TEST_EQUAL((outFiles(p, join, -1, 0) == std::vector<std::string>{"d1", "d1"}), true)
  TEST_EQUAL((outFiles(p, join, -1, 1) == std::vector<std::string>{"d2", "d2"}), true)
}
END_SECTION

START_SECTION([EXTRA] cycle detection)
{
  TOPPASPipeline p;
  int t1 = p.addNode(mkNode(NodeKind::TOOL, "0"));
  int t2 = p.addNode(mkNode(NodeKind::TOOL, "1"));
  p.addEdge(t1, t2, 0, 0);
  p.addEdge(t2, t1, 0, 0); // cycle
  std::string err;
  TEST_EQUAL(p.topoSort(err), false)
  TEST_NOT_EQUAL(err.size(), 0)
}
END_SECTION

START_SECTION((bool buildRoundPackages(int node_index, RoundPackages& pkg, std::string& error_msg) const))
{
  // node without input edges -> error
  TOPPASPipeline p;
  int tool = p.addNode(mkNode(NodeKind::TOOL, "0"));
  RoundPackages pkg;
  std::string err;
  TEST_EQUAL(p.buildRoundPackages(tool, pkg, err), false)
  TEST_NOT_EQUAL(err.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
