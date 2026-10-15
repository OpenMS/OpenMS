// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParamCWLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <sstream>

///////////////////////////

using namespace OpenMS;

START_TEST(ParamCWLFile, "$Id")

Param p;
p.setValue("ParamCWLFileTest:1:in", "", "input file", {"input file"});
p.setValue("ParamCWLFileTest:1:value", 1, "a number");

ToolInfo info;
info.name_ = "ParamCWLFileTest";
info.version_ = "1.0";

START_SECTION((static bool isSupported()))
{
#if defined(ENABLE_TDL)
  TEST_EQUAL(ParamCWLFile::isSupported(), true)
#else
  TEST_EQUAL(ParamCWLFile::isSupported(), false)
#endif
}
END_SECTION

START_SECTION((void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const))
{
  ParamCWLFile paramFile;
  std::string filename;
  NEW_TMP_FILE(filename)

  if (ParamCWLFile::isSupported())
  {
    // OpenMS exceptions only: TOPPBase (the caller) does not catch std ones
    TEST_EXCEPTION(Exception::UnableToCreateFile, paramFile.store("/does/not/exist/FileDoesNotExist.cwl", p, info))
    paramFile.store(filename, p, info);
    TEST_EQUAL(File::empty(filename), false)
  }
  else
  {
    // a build without TDL used to throw std::runtime_error (crashing every TOPP tool given
    // '-write_cwl') after opening, and thus truncating, the target file
    TEST_EXCEPTION(Exception::NotImplemented, paramFile.store(filename, p, info))
    TEST_EQUAL(File::exists(filename), false)
  }
}
END_SECTION

START_SECTION((void writeCWLToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const))
{
  ParamCWLFile paramFile;
  paramFile.flatHierarchy = true;
  std::stringstream ss;
  if (ParamCWLFile::isSupported())
  {
    paramFile.writeCWLToStream(&ss, p, info);
    TEST_EQUAL(ss.str().find("cwlVersion") != std::string::npos, true)
  }
  else
  {
    TEST_EXCEPTION(Exception::NotImplemented, paramFile.writeCWLToStream(&ss, p, info))
  }
}
END_SECTION

END_TEST
