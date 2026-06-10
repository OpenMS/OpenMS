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
#include <OpenMS/APPLICATIONS/PipelineExecutor.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(PipelineExecutor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PipelineExecutor* ptr = nullptr;
PipelineExecutor* null_ptr = nullptr;
START_SECTION(PipelineExecutor())
{
  ptr = new PipelineExecutor();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PipelineExecutor())
{
  delete ptr;
}
END_SECTION

START_SECTION((static std::vector<IOParam> collectIO(const Param& param, bool input)))
{
  // build a tool-like Param with input/output file params (the indices collectIO() assigns must match
  // what TOPPASToolVertex::getParameters_ produced, since .toppas edges reference them)
  Param p;
  p.setValue("in", "", "single input file");
  p.addTag("in", "input file");
  p.setValidStrings("in", ListUtils::create<std::string>("*.mzML"));

  p.setValue("database", "", "another single input file");
  p.addTag("database", "input file");
  p.setValidStrings("database", ListUtils::create<std::string>("*.fasta"));

  std::vector<std::string> empty_list;
  p.setValue("id_list", empty_list, "an input file list");
  p.addTag("id_list", "input file");
  p.setValidStrings("id_list", ListUtils::create<std::string>("*.idXML"));

  p.setValue("out", "", "single output file");
  p.addTag("out", "output file");
  p.setValidStrings("out", ListUtils::create<std::string>("*.featureXML"));

  p.setValue("out_dir", "", "output directory");
  p.addTag("out_dir", "output dir");

  p.setValue("threads", 1, "not an io param");

  // --- inputs ---
  std::vector<PipelineExecutor::IOParam> ins = PipelineExecutor::collectIO(p, true);
  TEST_EQUAL(ins.size(), 3)
  // single files sort before lists; within a type, by name -> database, in, then id_list
  TEST_EQUAL(ins[0].name, "database")
  TEST_EQUAL(ins[0].is_list, false)
  TEST_EQUAL(ins[0].ext, "fasta")
  TEST_EQUAL(ins[1].name, "in")
  TEST_EQUAL(ins[1].is_list, false)
  TEST_EQUAL(ins[1].ext, "mzML")
  TEST_EQUAL(ins[2].name, "id_list")
  TEST_EQUAL(ins[2].is_list, true)
  TEST_EQUAL(ins[2].ext, "idXML")

  // --- outputs ---
  std::vector<PipelineExecutor::IOParam> outs = PipelineExecutor::collectIO(p, false);
  TEST_EQUAL(outs.size(), 2)
  // output file sorts before output directory
  TEST_EQUAL(outs[0].name, "out")
  TEST_EQUAL(outs[0].is_list, false)
  TEST_EQUAL(outs[0].ext, "featureXML")
  TEST_EQUAL(outs[1].name, "out_dir")
  TEST_EQUAL(outs[1].is_list, false)
  TEST_EQUAL(outs[1].ext, "") // no restriction

  // a Param without any io tags yields nothing
  Param empty;
  empty.setValue("foo", 1);
  TEST_EQUAL(PipelineExecutor::collectIO(empty, true).empty(), true)
  TEST_EQUAL(PipelineExecutor::collectIO(empty, false).empty(), true)
}
END_SECTION

START_SECTION((static std::map<std::string, std::vector<std::string>> loadResourceFile(const std::string& filename)))
{
  // write a .trf-style resource file: each input node keyed by its topo number carries a 'url_list'
  Param trf;
  std::vector<std::string> n1 = {"file:///abs/a.mzML", "file:/abs/b.mzML"};
  std::vector<std::string> n3 = {"/abs/c.mzML"}; // no URL scheme
  trf.setValue("1:url_list", n1);
  trf.setValue("3:url_list", n3);

  std::string tmp;
  NEW_TMP_FILE(tmp)
  ParamXMLFile().store(tmp, trf);

  std::map<std::string, std::vector<std::string>> m = PipelineExecutor::loadResourceFile(tmp);
  TEST_EQUAL(m.size(), 2)
  TEST_EQUAL(m.count("1"), 1)
  TEST_EQUAL(m.count("3"), 1)
  // 'file://' and 'file:' schemes are stripped; a plain path is kept verbatim
  TEST_EQUAL(m["1"].size(), 2)
  TEST_EQUAL(m["1"][0], "/abs/a.mzML")
  TEST_EQUAL(m["1"][1], "/abs/b.mzML")
  TEST_EQUAL(m["3"].size(), 1)
  TEST_EQUAL(m["3"][0], "/abs/c.mzML")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
