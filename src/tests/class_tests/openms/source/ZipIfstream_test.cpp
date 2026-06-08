// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ZipIfstream.h>
///////////////////////////

#include <OpenMS/FORMAT/GzipIfstream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

START_TEST(ZipIfstream, "$Id$")

/////////////////////////////////////////////////////////////

ZipIfstream* ptr = nullptr;
ZipIfstream* null_ptr = nullptr;

START_SECTION(ZipIfstream())
  ptr = new ZipIfstream();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(~ZipIfstream())
  delete ptr;
END_SECTION

START_SECTION(void open(const char* filename))
  // Valid single-entry ZIP
  ZipIfstream zif;
  zif.open(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"));
  TEST_EQUAL(zif.isOpen(), true)
  TEST_EQUAL(zif.streamEnd(), false)
  zif.close();
END_SECTION

START_SECTION(size_t read(char* s, size_t n))
  // Read some bytes and verify they match the uncompressed file
  ZipIfstream zif(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"));
  char buf[256];
  size_t n = zif.read(buf, 255);
  TEST_EQUAL(n > 0, true)
  buf[n] = '\0';
  // First line of an mzML file should contain <?xml
  string content(buf, n);
  TEST_EQUAL(content.find("<?xml") != string::npos, true)
END_SECTION

START_SECTION(bool streamEnd() const)
  ZipIfstream zif(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"));
  TEST_EQUAL(zif.streamEnd(), false)
  // Read until EOF
  char buf[4096];
  while (!zif.streamEnd())
  {
    zif.read(buf, sizeof(buf));
  }
  TEST_EQUAL(zif.streamEnd(), true)
END_SECTION

START_SECTION(bool isOpen() const)
  ZipIfstream zif;
  TEST_EQUAL(zif.isOpen(), false)
  zif.open(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"));
  TEST_EQUAL(zif.isOpen(), true)
  zif.close();
  TEST_EQUAL(zif.isOpen(), false)
END_SECTION

START_SECTION(void close())
  ZipIfstream zif(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"));
  TEST_EQUAL(zif.isOpen(), true)
  zif.close();
  TEST_EQUAL(zif.isOpen(), false)
  TEST_EQUAL(zif.streamEnd(), true)
END_SECTION

START_SECTION([EXTRA] open with directory entry plus one file entry)
  // ZIP with dir + one file should work
  ZipIfstream zif;
  zif.open(OPENMS_GET_TEST_DATA_PATH("DirEntry.mzML.zip"));
  TEST_EQUAL(zif.isOpen(), true)
  char buf[256];
  size_t n = zif.read(buf, 255);
  TEST_EQUAL(n > 0, true)
  zif.close();
END_SECTION

START_SECTION([EXTRA] open empty ZIP throws FileEmpty)
  ZipIfstream zif;
  TEST_EXCEPTION(Exception::FileEmpty, zif.open(OPENMS_GET_TEST_DATA_PATH("Empty.zip")))
END_SECTION

START_SECTION([EXTRA] open multi-entry ZIP throws ParseError)
  ZipIfstream zif;
  TEST_EXCEPTION(Exception::ParseError, zif.open(OPENMS_GET_TEST_DATA_PATH("MultiEntry.mzML.zip")))
END_SECTION

START_SECTION([EXTRA] open nonexistent file throws FileNotFound)
  ZipIfstream zif;
  TEST_EXCEPTION(Exception::FileNotFound, zif.open("/nonexistent/path/to/file.zip"))
END_SECTION

START_SECTION([EXTRA] FileHandler::getTypeByFileName strips .zip extension)
  TEST_EQUAL(FileHandler::getTypeByFileName("sample.mzML.zip"), FileTypes::MZML)
  TEST_EQUAL(FileHandler::getTypeByFileName("data.featureXML.zip"), FileTypes::FEATUREXML)
  TEST_EQUAL(FileHandler::getTypeByFileName("ids.idXML.zip"), FileTypes::IDXML)
END_SECTION

START_SECTION([EXTRA] End-to-end: load mzML from ZIP via FileHandler)
  MSExperiment exp_zip, exp_plain;

  FileHandler().loadExperiment(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip"), exp_zip);
  FileHandler().loadExperiment(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML"), exp_plain);

  TEST_EQUAL(exp_zip.size(), exp_plain.size())
  for (Size i = 0; i < exp_zip.size(); ++i)
  {
    TEST_EQUAL(exp_zip[i].size(), exp_plain[i].size())
    TEST_EQUAL(exp_zip[i].getMSLevel(), exp_plain[i].getMSLevel())
  }
END_SECTION

END_TEST
