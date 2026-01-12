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
#include <OpenMS/METADATA/NativeIDParser.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NativeIDParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static bool isNativeID(const String& id)))
{
  // Test recognized native ID prefixes
  TEST_EQUAL(NativeIDParser::isNativeID("scan=123"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("scanID=456"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("controllerType=0 controllerNumber=1 scan=100"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("function=2 process=1 scan=100"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("sample=1 period=1 cycle=42 experiment=1"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("index=789"), true);
  TEST_EQUAL(NativeIDParser::isNativeID("spectrum=101112"), true);

  // Test non-native IDs
  TEST_EQUAL(NativeIDParser::isNativeID("123"), false);
  TEST_EQUAL(NativeIDParser::isNativeID("file=42"), false);  // file= is not a recognized prefix for isNativeID
  TEST_EQUAL(NativeIDParser::isNativeID("scanId=42"), false);  // case-sensitive: scanId != scanID
  TEST_EQUAL(NativeIDParser::isNativeID(""), false);
  TEST_EQUAL(NativeIDParser::isNativeID("some_random_string"), false);
  TEST_EQUAL(NativeIDParser::isNativeID("SCAN=123"), false);  // case-sensitive
}
END_SECTION


START_SECTION((static std::string getRegExFromNativeID(const String& native_id)))
{
  // Test Thermo format: "controllerType=0 controllerNumber=1 scan=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("controllerType=0 controllerNumber=1 scan=100"), R"(scan=(?<GROUP>\d+))");

  // Test Waters format: "function= process= scan=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("function=2 process=1 scan=100"), R"(scan=(?<GROUP>\d+))");

  // Test simple scan format: "scan=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("scan=123"), R"(scan=(?<GROUP>\d+))");

  // Test index format: "index=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("index=456"), R"(index=(?<GROUP>\d+))");

  // Test Agilent MassHunter format: "scanId=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("scanId=789"), R"(scanId=(?<GROUP>\d+))");

  // Test spectrum format: "spectrum=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("spectrum=101112"), R"(spectrum=(?<GROUP>\d+))");

  // Test Bruker FID format: "file=NUMBER"
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("file=42"), R"(file=(?<GROUP>\d+))");

  // Test plain number (fallback)
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("123"), R"((?<GROUP>\d+))");

  // Test unknown format falls back to plain number extraction
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID("unknown=123"), R"((?<GROUP>\d+))");

  // Test empty string falls back to plain number extraction
  TEST_EQUAL(NativeIDParser::getRegExFromNativeID(""), R"((?<GROUP>\d+))");
}
END_SECTION


START_SECTION((static Int extractScanNumber(const String& native_id, const boost::regex& scan_regexp, bool no_error = false)))
{
  // Test successful extraction with spectrum= format
  boost::regex re_spectrum("spectrum=(?<SCAN>\\d+)");
  TEST_EQUAL(NativeIDParser::extractScanNumber("spectrum=42", re_spectrum), 42);
  TEST_EQUAL(NativeIDParser::extractScanNumber("spectrum=0", re_spectrum), 0);
  TEST_EQUAL(NativeIDParser::extractScanNumber("spectrum=99999", re_spectrum), 99999);

  // Test successful extraction with scan= format
  boost::regex re_scan("scan=(?<SCAN>\\d+)");
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=123", re_scan), 123);
  TEST_EQUAL(NativeIDParser::extractScanNumber("controllerType=0 controllerNumber=1 scan=456", re_scan), 456);

  // Test no_error=true returns -1 instead of throwing
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", re_spectrum, true), -1);
  TEST_EQUAL(NativeIDParser::extractScanNumber("no_match_here", re_spectrum, true), -1);
  TEST_EQUAL(NativeIDParser::extractScanNumber("", re_spectrum, true), -1);

  // Test exception when no_error=false (default) and no match
  TEST_EXCEPTION(Exception::ParseError, NativeIDParser::extractScanNumber("scan=42", re_spectrum));
  TEST_EXCEPTION(Exception::ParseError, NativeIDParser::extractScanNumber("no_match", re_spectrum, false));

  // Test multiple matches - should return last one
  boost::regex re_multi("(?<SCAN>\\d+)");
  TEST_EQUAL(NativeIDParser::extractScanNumber("1 2 3 42", re_multi), 42);

  // Test edge cases
  boost::regex re_index("index=(?<SCAN>\\d+)");
  TEST_EQUAL(NativeIDParser::extractScanNumber("index=0", re_index), 0);
}
END_SECTION


START_SECTION((static Int extractScanNumber(const String& native_id, const String& native_id_type_accession)))
{
  // Test Thermo nativeID format (MS:1000768) - scan=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:1000768"), 42);
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=0", "MS:1000768"), 0);
  TEST_EQUAL(NativeIDParser::extractScanNumber("controllerType=0 controllerNumber=1 scan=123", "MS:1000768"), 123);

  // Test Waters nativeID format (MS:1000769) - function=X process=Y scan=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:1000769"), 42);
  TEST_EQUAL(NativeIDParser::extractScanNumber("function=1 process=1 scan=99", "MS:1000769"), 99);

  // Test WIFF nativeID format (MS:1000770) - cycle * 1000 + experiment
  TEST_EQUAL(NativeIDParser::extractScanNumber("sample=1 period=1 cycle=42 experiment=1", "MS:1000770"), 42001);
  TEST_EQUAL(NativeIDParser::extractScanNumber("sample=1 period=1 cycle=1 experiment=999", "MS:1000770"), 1999);
  TEST_EQUAL(NativeIDParser::extractScanNumber("sample=1 period=1 cycle=100 experiment=50", "MS:1000770"), 100050);

  // Test Bruker BAF nativeID format (MS:1000771) - scan=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:1000771"), 42);

  // Test Bruker U2 nativeID format (MS:1000772) - scan=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:1000772"), 42);

  // Test Bruker FID nativeID format (MS:1000773) - file=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("file=42", "MS:1000773"), 42);
  TEST_EQUAL(NativeIDParser::extractScanNumber("file=0", "MS:1000773"), 0);

  // Test multiple spectra per file (MS:1000774) - index=NUMBER -> returns index+1
  TEST_EQUAL(NativeIDParser::extractScanNumber("index=42", "MS:1000774"), 43);
  TEST_EQUAL(NativeIDParser::extractScanNumber("index=0", "MS:1000774"), 1);

  // Test single peak list nativeID format (MS:1000775) - file=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("file=42", "MS:1000775"), 42);

  // Test Thermo/Bruker TDF nativeID format (MS:1000776) - scan=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:1000776"), 42);

  // Test spectrum identifier nativeID format (MS:1000777) - spectrum=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("spectrum=42", "MS:1000777"), 42);

  // Test Agilent MassHunter nativeID format (MS:1001508) - scanId=NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("scanId=42", "MS:1001508"), 42);

  // Test mzML unique identifier (MS:1001530) - plain NUMBER
  TEST_EQUAL(NativeIDParser::extractScanNumber("42", "MS:1001530"), 42);
  TEST_EQUAL(NativeIDParser::extractScanNumber("0", "MS:1001530"), 0);

  // Test invalid accession returns -1
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "MS:9999999"), -1);
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", "invalid"), -1);
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=42", ""), -1);

  // Test invalid native_id for given accession returns -1
  TEST_EQUAL(NativeIDParser::extractScanNumber("spectrum=42", "MS:1000768"), -1);  // expects scan=
  TEST_EQUAL(NativeIDParser::extractScanNumber("", "MS:1000768"), -1);
  TEST_EQUAL(NativeIDParser::extractScanNumber("malformed", "MS:1000768"), -1);

  // Test merged spectra - should return last scan number
  TEST_EQUAL(NativeIDParser::extractScanNumber("controllerType=0 controllerNumber=1 scan=100 merged controllerType=0 controllerNumber=1 scan=200", "MS:1000768"), 200);
}
END_SECTION


// Additional edge case tests
START_SECTION([EXTRA] Edge cases and error handling)
{
  // Test very large numbers
  TEST_EQUAL(NativeIDParser::extractScanNumber("scan=999999999", "MS:1000768"), 999999999);

  // Test whitespace handling in WIFF format
  TEST_EQUAL(NativeIDParser::extractScanNumber("sample=1 period=1 cycle=42   experiment=1", "MS:1000770"), 42001);

  // Test isNativeID with various whitespace
  TEST_EQUAL(NativeIDParser::isNativeID(" scan=123"), false);  // leading space
  TEST_EQUAL(NativeIDParser::isNativeID("scan =123"), false);  // space before =

  // Test getRegExFromNativeID preserves regex structure
  std::string regex = NativeIDParser::getRegExFromNativeID("scan=123");
  boost::regex re(regex);
  boost::smatch match;
  std::string test = "scan=456";
  TEST_EQUAL(boost::regex_search(test, match, re), true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
