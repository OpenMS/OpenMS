// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <sstream>

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

START_SECTION(XMLFile(const std::string &schema_location, const std::string &version))
	NOT_TESTABLE
END_SECTION

START_SECTION(bool isValid(const std::string &filename,  std::ostream& os = std::cerr) )
	XMLFile f("","");
	TEST_EXCEPTION(Exception::NotImplemented, f.isValid("", std::cerr))
END_SECTION

START_SECTION(const std::string& getVersion() const)
	XMLFile f("","1.567");
	TEST_EQUAL( f.getVersion(),"1.567")
END_SECTION


START_SECTION(([EXTRA] std::string writeXMLEscape(const std::string& to_escape)))
  std::string s1("nothing_to_escape. Just a regular string...");
  std::string s2("This string contains an ampersand, &, which must be escaped.");
  std::string s3("This string also contains characters which is not allowed, and must be escaped; the characters are '>' and \"<\"");

  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s1), "nothing_to_escape. Just a regular string...")
  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s2), "This string contains an ampersand, &amp;, which must be escaped.")
  TEST_STRING_EQUAL(XMLHandler::writeXMLEscape(s3), "This string also contains characters which is not allowed, and must be escaped; the characters are &apos;&gt;&apos; and &quot;&lt;&quot;");
END_SECTION

START_SECTION(static DataValue fromXSDString(const std::string& type, const std::string& value))
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

START_SECTION(([EXTRA] remote XML error diagnostics do not reopen URIs as local files))
{
  for (const std::string uri : {"http://unused.invalid/test.mzML", "https://unused.invalid/test.mzML",
                                "ftp://unused.invalid/test.mzML", "s3://unused/test.mzML"})
  {
    XMLHandler handler(uri, "1.0");
    TEST_EXCEPTION(Exception::ParseError, handler.fatalError(XMLHandler::LOAD, "offline parser error"))
  }
  class DiagnosticHandler : public XMLHandler
  {
  public:
    using XMLHandler::XMLHandler;
    const std::string& diagnosticFile() const { return file_; }
    void throwDirect() const
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_, "offline direct handler error");
    }
  };
  DiagnosticHandler handler(std::string("s3://") + "user:password@" + "bucket/key?signature=private-value", "1.0");
  TEST_EQUAL(handler.diagnosticFile(), "s3://bucket/key")
  try
  {
    handler.fatalError(XMLHandler::LOAD, "offline parser error");
    TEST_TRUE(false)
  }
  catch (const Exception::ParseError& error)
  {
    const std::string message = error.getMessage();
    TEST_TRUE(message.find("offline parser error") != std::string::npos)
    TEST_TRUE(message.find("password") == std::string::npos)
    TEST_TRUE(message.find("signature") == std::string::npos)
    TEST_TRUE(message.find("private-value") == std::string::npos)
  }
  try
  {
    handler.throwDirect();
    TEST_TRUE(false)
  }
  catch (const Exception::ParseError& error)
  {
    const std::string message = error.getMessage();
    TEST_TRUE(message.find("offline direct handler error") != std::string::npos)
    TEST_TRUE(message.find("password") == std::string::npos)
    TEST_TRUE(message.find("signature") == std::string::npos)
    TEST_TRUE(message.find("private-value") == std::string::npos)
  }
  std::ostringstream logged;
  getThreadLocalLogError().insert(logged);
  getThreadLocalLogWarn().insert(logged);
  getThreadLocalLogDebug().insert(logged);
  handler.error(XMLHandler::LOAD, "offline non-fatal error");
  handler.warning(XMLHandler::LOAD, "offline warning");
  getThreadLocalLogError().remove(logged);
  getThreadLocalLogWarn().remove(logged);
  getThreadLocalLogDebug().remove(logged);
  TEST_TRUE(logged.str().find("offline non-fatal error") != std::string::npos)
  TEST_TRUE(logged.str().find("offline warning") != std::string::npos)
  TEST_TRUE(logged.str().find("password") == std::string::npos)
  TEST_TRUE(logged.str().find("signature") == std::string::npos)
  TEST_TRUE(logged.str().find("private-value") == std::string::npos)
}
END_SECTION

#ifndef WITH_S3
START_SECTION(([EXTRA] XML readers preserve remote metadata before reporting disabled S3))
{
  std::string uri = "s3://unused-bucket/object.featureXML?region=eu-west-1";
  FeatureMap features;
  TEST_EXCEPTION(Exception::ParseError, FeatureXMLFile().load(uri, features))
  TEST_EQUAL(features.getLoadedFilePath(), "s3://unused-bucket/object.featureXML")
  TEST_EQUAL(features.getLoadedFileType(), FileTypes::FEATUREXML)
  uri = "s3://unused-bucket/object.consensusXML?region=eu-west-1";
  ConsensusMap consensus;
  TEST_EXCEPTION(Exception::ParseError, ConsensusXMLFile().load(uri, consensus))
  TEST_EQUAL(consensus.getLoadedFilePath(), "s3://unused-bucket/object.consensusXML")
  TEST_EQUAL(consensus.getLoadedFileType(), FileTypes::CONSENSUSXML)
  uri = "s3://unused-bucket/object.mzXML?region=eu-west-1";
  MSExperiment experiment;
  TEST_EXCEPTION(Exception::ParseError, MzXMLFile().load(uri, experiment))
  TEST_EQUAL(experiment.getLoadedFilePath(), "s3://unused-bucket/object.mzXML")
  TEST_EQUAL(experiment.getLoadedFileType(), FileTypes::MZXML)
  uri = "s3://unused-bucket/object.mzData?region=eu-west-1";
  TEST_EXCEPTION(Exception::ParseError, MzDataFile().load(uri, experiment))
  TEST_EQUAL(experiment.getLoadedFilePath(), "s3://unused-bucket/object.mzData")
  TEST_EQUAL(experiment.getLoadedFileType(), FileTypes::MZDATA)
  uri = "s3://unused-bucket/object?filename=object.featureXML";
  TEST_EXCEPTION(Exception::ParseError, FeatureXMLFile().load(uri, features))
  TEST_EQUAL(features.getLoadedFilePath(), "s3://unused-bucket/object")
  TEST_EQUAL(features.getLoadedFileType(), FileTypes::UNKNOWN)
}
END_SECTION
#endif
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
