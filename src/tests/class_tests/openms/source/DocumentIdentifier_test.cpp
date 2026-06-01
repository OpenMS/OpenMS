// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/FORMAT/FileTypes.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DocumentIdentifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DocumentIdentifier* ptr = nullptr;
DocumentIdentifier* nullPointer = nullptr;
START_SECTION(DocumentIdentifier())
{
	ptr = new DocumentIdentifier();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DocumentIdentifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((DocumentIdentifier(const DocumentIdentifier &source)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	di1.setLoadedFilePath( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
	di1.setLoadedFileType( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));

	DocumentIdentifier di2(di1);
	TEST_EQUAL(di2.getIdentifier(), "this is a test");
	TEST_EQUAL(di2.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"))
  TEST_EQUAL(FileTypes::typeToName(di2.getLoadedFileType()) == "unknown", true)
}
END_SECTION

START_SECTION((DocumentIdentifier& operator=(const DocumentIdentifier &source)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	di1.setLoadedFilePath( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
	di1.setLoadedFileType( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));

	DocumentIdentifier di2 = di1;
	TEST_EQUAL(di2.getIdentifier(), "this is a test");
	TEST_EQUAL(di2.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"))
  TEST_EQUAL(FileTypes::typeToName(di2.getLoadedFileType()) == "unknown", true)
}
END_SECTION

START_SECTION((bool operator==(const DocumentIdentifier &rhs) const))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	DocumentIdentifier di2(di1);
	TEST_TRUE(di1 == di2)
}
END_SECTION

START_SECTION((void setIdentifier(const String &id)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	TEST_EQUAL(di1.getIdentifier(), "this is a test")
}
END_SECTION

START_SECTION((const String& getIdentifier() const))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setLoadedFileType(const String &file_name)))
{
  DocumentIdentifier di1;
	di1.setLoadedFileType( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
  TEST_EQUAL(FileTypes::typeToName(di1.getLoadedFileType()), "unknown")
}
END_SECTION

START_SECTION((const FileTypes::Type& getLoadedFileType() const))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setLoadedFilePath(const String &file_name)))
{
  DocumentIdentifier di1;
	di1.setLoadedFilePath( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
	TEST_EQUAL(di1.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"))
}
END_SECTION

START_SECTION((const String& getLoadedFilePath() const))
{
	// tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void swap(DocumentIdentifier& from)))
{
  DocumentIdentifier di1;
	di1.setIdentifier("this is a test");
	di1.setLoadedFilePath( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
	di1.setLoadedFileType( OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"));
	DocumentIdentifier di2;
	di1.swap(di2);
	TEST_EQUAL(di1.getIdentifier().empty(), true)
	TEST_EQUAL(di1.getIdentifier().empty(), true)
	TEST_EQUAL(di1.getIdentifier().empty(), true)
	TEST_EQUAL(di2.getIdentifier() == "this is a test", true)
  TEST_EQUAL(di2.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt"))
  TEST_EQUAL(FileTypes::typeToName(di2.getLoadedFileType()) == "unknown", true)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



