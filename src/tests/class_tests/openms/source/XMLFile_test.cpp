// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

///////////////////////////

START_TEST(XMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

XMLFile* ptr = nullptr;
XMLFile* nullPointer = nullptr;

START_SECTION(XMLFile())
	ptr = new XMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer)
END_SECTION

START_SECTION(~XMLFile())
	delete ptr;
END_SECTION

START_SECTION(XMLFile(const String &schema_location, const String &version))
	NOT_TESTABLE
END_SECTION

START_SECTION(bool isValid(const String &filename,  std::ostream& os = std::cerr) )
	XMLFile f("","");
	TEST_EXCEPTION(Exception::NotImplemented, f.isValid("", std::cerr))
END_SECTION

START_SECTION(const String& getVersion() const)
	XMLFile f("","1.567");
	TEST_EQUAL( f.getVersion(),"1.567")
END_SECTION


START_SECTION(([EXTRA] String writeXMLEscape(const String& to_escape)))
  String s1("nothing_to_escape. Just a regular string...");
  String s2("This string contains an ampersand, &, which must be escaped.");
  String s3("This string also contains characters which is not allowed, and must be escaped; the characters are '>' and \"<\"");

  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s1), "nothing_to_escape. Just a regular string...")
  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s2), "This string contains an ampersand, &amp;, which must be escaped.")
  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s3), "This string also contains characters which is not allowed, and must be escaped; the characters are &apos;&gt;&apos; and &quot;&lt;&quot;");
END_SECTION

START_SECTION(static DataValue fromXSDString(const String& type, const String& value))
{
  TEST_EQUAL((Int64)XMLHandler::fromXSDString("xsd:int", "2147483647"), 2147483647)
  TEST_EQUAL((Int64)XMLHandler::fromXSDString("xsd:long", "9223372036854775807"), 9223372036854775807)
  TEST_EQUAL((double)XMLHandler::fromXSDString("xsd:decimal", "123.45"), 123.45)
  TEST_EQUAL((Int64)XMLHandler::fromXSDString("xsd:unsignedLong", "9223372036854775807"), 9223372036854775807)

  // input exceeds valid range
  TEST_EXCEPTION(Exception::ConversionError, XMLHandler::fromXSDString("xsd:int", "2147483648"))           // +1 larger than 2^31-1
  TEST_EXCEPTION(Exception::ConversionError, XMLHandler::fromXSDString("xsd:long", "9223372036854775808")) // +1 larger than 2^63-1

  // things we SHOULD support, but don't, due to using a signed 64bit type in DataValue
  TEST_EXCEPTION(Exception::ConversionError, XMLHandler::fromXSDString("xsd:unsignedLong", "9223372036854775808")) // +1 larger than 2^63-1; 'xsd:unsignedLong' up to 2^64-1

  // things which are really hard to support (arbitrarily large numbers)
  TEST_EXCEPTION(Exception::ConversionError, XMLHandler::fromXSDString("xsd:integer", "9223372036854775808")) // +1 larger than 2^63-1; 'xsd:integer' can be any number... hard to support :)
  TEST_EXCEPTION(Exception::ConversionError,
                 XMLHandler::fromXSDString("xsd:negativeInteger", "-9223372036854775809")) // -1 smaller than 2^63; 'xsd:negativeInteger' can be any negative number... hard to support :)
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
