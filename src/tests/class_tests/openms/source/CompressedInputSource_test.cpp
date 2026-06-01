// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
using namespace OpenMS;


///////////////////////////

START_TEST(CompressedInputSource, "$Id$")

xercesc::XMLPlatformUtils::Initialize();
CompressedInputSource* ptr = nullptr;
CompressedInputSource* nullPointer = nullptr;

START_SECTION(CompressedInputSource(const String& file_path, const char * header, xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager))
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  ptr = new CompressedInputSource(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"), bz);
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~CompressedInputSource()))
  delete ptr;
END_SECTION

START_SECTION(CompressedInputSource(const XMLCh *const file_path, const char *header, xercesc::MemoryManager *const manager=xercesc::XMLPlatformUtils::fgMemoryManager))
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  String filename(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  ptr = new CompressedInputSource(Internal::StringManager().convert(filename.c_str()).c_str(), bz);
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
END_SECTION


START_SECTION(virtual xercesc::BinInputStream* makeStream() const)
  char header[3];
  header[0] = 'B';
  header[1] = 'Z';
  header[2] = '\0';
  String bz = String(header);
  CompressedInputSource source(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist"), bz);
  TEST_EXCEPTION(Exception::FileNotFound,source.makeStream())
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
