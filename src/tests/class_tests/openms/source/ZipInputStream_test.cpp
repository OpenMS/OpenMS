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

#include <OpenMS/FORMAT/ZipInputStream.h>

#include <string>

///////////////////////////

START_TEST(ZipInputStream, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace xercesc;
using namespace std;

// A committed fixture: a ZIP archive containing a single mzML (XML) file.
const std::string zip_file = OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML.zip");

ZipInputStream* ptr = nullptr;
ZipInputStream* null_ptr = nullptr;

START_SECTION((explicit ZipInputStream(const std::string& file_name)))
{
  ptr = new ZipInputStream(zip_file);
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getIsOpen(), true)
}
END_SECTION

START_SECTION((~ZipInputStream()))
{
  delete ptr;
}
END_SECTION

START_SECTION((explicit ZipInputStream(const char* const file_name)))
{
  ZipInputStream zis(zip_file.c_str());
  TEST_EQUAL(zis.getIsOpen(), true)
}
END_SECTION

START_SECTION((bool getIsOpen() const))
{
  ZipInputStream zis(zip_file);
  TEST_EQUAL(zis.getIsOpen(), true)
}
END_SECTION

START_SECTION((XMLFilePos curPos() const))
{
  ZipInputStream zis(zip_file);
  TEST_EQUAL(zis.curPos(), 0)

  XMLByte buf[16];
  zis.readBytes(buf, 16);
  TEST_EQUAL(zis.curPos(), 16) // advances by the number of bytes read
}
END_SECTION

START_SECTION((XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read)))
{
  ZipInputStream zis(zip_file);

  // first read decompresses the start of the contained mzML (XML declaration)
  XMLByte buf[32];
  XMLSize_t n = zis.readBytes(buf, 32);
  TEST_EQUAL(n, 32)
  std::string head((const char*)buf, (size_t)n);
  TEST_EQUAL(head.substr(0, 5), "<?xml")
  TEST_EQUAL(zis.curPos(), 32)

  // read the rest of the stream; curPos must track the total read
  XMLSize_t total = n;
  XMLByte big[4096];
  XMLSize_t r;
  while ((r = zis.readBytes(big, 4096)) > 0)
  {
    total += r;
  }
  // the fixture decompresses to a fixed number of bytes
  TEST_EQUAL(total, 28008)
  TEST_EQUAL(zis.curPos(), 28008)

  // further reads at end-of-stream return 0 and do not advance
  TEST_EQUAL(zis.readBytes(big, 4096), 0)
  TEST_EQUAL(zis.curPos(), 28008)
}
END_SECTION

START_SECTION((const XMLCh* getContentType() const))
{
  ZipInputStream zis(zip_file);
  TEST_EQUAL(zis.getContentType() == nullptr, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
