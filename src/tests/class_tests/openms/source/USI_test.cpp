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
#include <OpenMS/METADATA/USI.h>
///////////////////////////

#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(USI, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

USI* ptr = nullptr;
USI* null_ptr = nullptr;

START_SECTION((USI()))
{
  ptr = new USI();
  TEST_NOT_EQUAL(ptr, null_ptr);
  TEST_EQUAL(ptr->isValid(), false);
}
END_SECTION

START_SECTION((~USI()))
{
  delete ptr;
}
END_SECTION

START_SECTION((USI(const String& collection, const String& ms_run, IndexType index_type, const String& index, const String& interpretation)))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345", "PEPTIDEK/2");
  TEST_EQUAL(usi.isValid(), true);
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000561");
  TEST_STRING_EQUAL(usi.getMSRun(), "sample.mzML");
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::SCAN);
  TEST_STRING_EQUAL(usi.getIndex(), "12345");
  TEST_STRING_EQUAL(usi.getInterpretation(), "PEPTIDEK/2");
}
END_SECTION

START_SECTION((explicit USI(const String& usi_string)))
{
  // Test valid USI parsing
  USI usi("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2");
  TEST_EQUAL(usi.isValid(), true);
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000561");
  TEST_STRING_EQUAL(usi.getMSRun(), "sample.mzML");
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::SCAN);
  TEST_STRING_EQUAL(usi.getIndex(), "12345");
  TEST_STRING_EQUAL(usi.getInterpretation(), "PEPTIDEK/2");

  // Test USI with ProForma interpretation containing colon (e.g., UNIMOD accessions)
  USI usi3("mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2");
  TEST_EQUAL(usi3.isValid(), true);
  TEST_STRING_EQUAL(usi3.getInterpretation(), "EM[UNIMOD:35]K/2");
  TEST_STRING_EQUAL(usi3.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2");

  // Test USI without interpretation
  USI usi4("mzspec:PXD000561:sample.mzML:scan:12345");
  TEST_EQUAL(usi4.isValid(), true);
  TEST_STRING_EQUAL(usi4.getInterpretation(), "");

  // Test invalid USI
  TEST_EXCEPTION(Exception::ParseError, USI("invalid_usi_string"));
  TEST_EXCEPTION(Exception::ParseError, USI("mzspec:PXD000561"));  // too few parts
}
END_SECTION

START_SECTION((USI(const USI&)))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2(usi1);
  TEST_EQUAL(usi1 == usi2, true);
}
END_SECTION

START_SECTION((USI(USI&&)))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2(std::move(usi1));
  TEST_STRING_EQUAL(usi2.getCollection(), "PXD000561");
}
END_SECTION

START_SECTION((USI& operator=(const USI&)))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2;
  usi2 = usi1;
  TEST_EQUAL(usi1 == usi2, true);
}
END_SECTION

START_SECTION((USI& operator=(USI&&)))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2;
  usi2 = std::move(usi1);
  TEST_STRING_EQUAL(usi2.getCollection(), "PXD000561");
}
END_SECTION

START_SECTION((bool operator==(const USI& rhs) const))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi3("PXD000562", "sample.mzML", USI::IndexType::SCAN, "12345");
  
  TEST_EQUAL(usi1 == usi2, true);
  TEST_EQUAL(usi1 == usi3, false);
}
END_SECTION

START_SECTION((bool operator!=(const USI& rhs) const))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2("PXD000562", "sample.mzML", USI::IndexType::SCAN, "12345");
  
  TEST_EQUAL(usi1 != usi2, true);
}
END_SECTION

START_SECTION((bool operator<(const USI& rhs) const))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  USI usi2("PXD000562", "sample.mzML", USI::IndexType::SCAN, "12345");
  
  TEST_EQUAL(usi1 < usi2, true);
}
END_SECTION

START_SECTION((bool isValid() const))
{
  USI usi1;
  TEST_EQUAL(usi1.isValid(), false);
  
  USI usi2("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_EQUAL(usi2.isValid(), true);
  
  USI usi3("", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_EQUAL(usi3.isValid(), false);
}
END_SECTION

START_SECTION((static bool isValidUSI(const String& usi_string)))
{
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345"), true);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"), true);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"), true);
  TEST_EQUAL(USI::isValidUSI("invalid_string"), false);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561"), false);
  TEST_EQUAL(USI::isValidUSI(""), false);
}
END_SECTION

START_SECTION((const String& getCollection() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000561");
}
END_SECTION

START_SECTION((void setCollection(const String& collection)))
{
  USI usi;
  usi.setCollection("PXD000562");
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000562");
}
END_SECTION

START_SECTION((const String& getMSRun() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getMSRun(), "sample.mzML");
}
END_SECTION

START_SECTION((void setMSRun(const String& ms_run)))
{
  USI usi;
  usi.setMSRun("test.mzML");
  TEST_STRING_EQUAL(usi.getMSRun(), "test.mzML");
}
END_SECTION

START_SECTION((IndexType getIndexType() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::INDEX, "12345");
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::INDEX);
}
END_SECTION

START_SECTION((void setIndexType(IndexType index_type)))
{
  USI usi;
  usi.setIndexType(USI::IndexType::NATIVEID);
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::NATIVEID);
}
END_SECTION

START_SECTION((const String& getIndex() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getIndex(), "12345");
}
END_SECTION

START_SECTION((void setIndex(const String& index)))
{
  USI usi;
  usi.setIndex("99999");
  TEST_STRING_EQUAL(usi.getIndex(), "99999");
}
END_SECTION

START_SECTION((const String& getInterpretation() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345", "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi.getInterpretation(), "PEPTIDEK/2");
}
END_SECTION

START_SECTION((void setInterpretation(const String& interpretation)))
{
  USI usi;
  usi.setInterpretation("SEQUENCE/3");
  TEST_STRING_EQUAL(usi.getInterpretation(), "SEQUENCE/3");
}
END_SECTION

START_SECTION((bool hasInterpretation() const))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345", "PEPTIDEK/2");
  TEST_EQUAL(usi1.hasInterpretation(), true);
  
  USI usi2("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_EQUAL(usi2.hasInterpretation(), false);
}
END_SECTION

START_SECTION((String toString() const))
{
  USI usi1("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345", "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi1.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2");
  
  USI usi2("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:sample.mzML:scan:12345");
  
  USI usi3("PXD000561", "sample.mzML", USI::IndexType::INDEX, "0");
  TEST_STRING_EQUAL(usi3.toString(), "mzspec:PXD000561:sample.mzML:index:0");
  
  USI usi_empty;
  TEST_STRING_EQUAL(usi_empty.toString(), "");
}
END_SECTION

START_SECTION((bool fromString(const String& usi_string)))
{
  USI usi;
  
  TEST_EQUAL(usi.fromString("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"), true);
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000561");
  TEST_STRING_EQUAL(usi.getMSRun(), "sample.mzML");
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::SCAN);
  TEST_STRING_EQUAL(usi.getIndex(), "12345");
  TEST_STRING_EQUAL(usi.getInterpretation(), "PEPTIDEK/2");
  
  // Test with index type
  TEST_EQUAL(usi.fromString("mzspec:PXD000561:sample.mzML:index:0"), true);
  TEST_EQUAL(usi.getIndexType(), USI::IndexType::INDEX);
  
  // Test invalid formats
  TEST_EQUAL(usi.fromString("invalid"), false);
  TEST_EQUAL(usi.fromString("mzspec:PXD000561"), false);
  TEST_EQUAL(usi.fromString(""), false);
}
END_SECTION

START_SECTION((static USI createFromScanNumber(const String& dataset_id, const String& filename, int scan_number, const String& interpretation)))
{
  USI usi1 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345);
  TEST_STRING_EQUAL(usi1.toString(), "mzspec:PXD000561:sample.mzML:scan:12345");

  USI usi2 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345, "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2");

  USI usi3 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345, "EM[UNIMOD:35]K/2");
  TEST_STRING_EQUAL(usi3.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2");
}
END_SECTION

START_SECTION((static USI createFromNativeID(const String& dataset_id, const String& filename, const String& native_id)))
{
  // Test with extractable scan number
  USI usi1 = USI::createFromNativeID("PXD000561", "sample.mzML", "scan=12345");
  TEST_STRING_EQUAL(usi1.toString(), "mzspec:PXD000561:sample.mzML:scan:12345");
  
  USI usi2 = USI::createFromNativeID("PXD000561", "sample.mzML", "controllerType=0 controllerNumber=1 scan=54321");
  TEST_STRING_EQUAL(usi2.getIndex(), "54321");
  
  // Test with non-extractable native ID (falls back to nativeId type)
  USI usi3 = USI::createFromNativeID("PXD000561", "sample.mzML", "custom_id_format");
  TEST_EQUAL(usi3.getIndexType(), USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(usi3.getIndex(), "custom_id_format");
}
END_SECTION

START_SECTION((static std::optional<int> extractScanNumberFromNativeID(const String& native_id)))
{
  // Test various native ID formats
  auto result1 = USI::extractScanNumberFromNativeID("scan=12345");
  TEST_EQUAL(result1.has_value(), true);
  TEST_EQUAL(result1.value(), 12345);
  
  auto result2 = USI::extractScanNumberFromNativeID("controllerType=0 controllerNumber=1 scan=54321");
  TEST_EQUAL(result2.has_value(), true);
  TEST_EQUAL(result2.value(), 54321);
  
  auto result3 = USI::extractScanNumberFromNativeID("scanId=99999");
  TEST_EQUAL(result3.has_value(), true);
  TEST_EQUAL(result3.value(), 99999);
  
  auto result4 = USI::extractScanNumberFromNativeID("index=100");
  TEST_EQUAL(result4.has_value(), true);
  TEST_EQUAL(result4.value(), 100);
  
  auto result5 = USI::extractScanNumberFromNativeID("12345");
  TEST_EQUAL(result5.has_value(), true);
  TEST_EQUAL(result5.value(), 12345);
  
  auto result6 = USI::extractScanNumberFromNativeID("custom_format_no_number");
  TEST_EQUAL(result6.has_value(), false);
}
END_SECTION

START_SECTION((static String indexTypeToString(IndexType index_type)))
{
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::SCAN), "scan");
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::INDEX), "index");
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::NATIVEID), "nativeId");
}
END_SECTION

START_SECTION((static IndexType indexTypeFromString(const String& type_string)))
{
  TEST_EQUAL(USI::indexTypeFromString("scan"), USI::IndexType::SCAN);
  TEST_EQUAL(USI::indexTypeFromString("SCAN"), USI::IndexType::SCAN);
  TEST_EQUAL(USI::indexTypeFromString("index"), USI::IndexType::INDEX);
  TEST_EQUAL(USI::indexTypeFromString("INDEX"), USI::IndexType::INDEX);
  TEST_EQUAL(USI::indexTypeFromString("nativeId"), USI::IndexType::NATIVEID);
  TEST_EQUAL(USI::indexTypeFromString("NATIVEID"), USI::IndexType::NATIVEID);
  
  TEST_EXCEPTION(Exception::InvalidValue, USI::indexTypeFromString("invalid_type"));
}
END_SECTION

START_SECTION((static String extractBasename(const String& filepath)))
{
  // Test Unix-style paths
  TEST_STRING_EQUAL(USI::extractBasename("/path/to/sample.mzML"), "sample.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("/data/experiments/run1.mzML"), "run1.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("/sample.mzML"), "sample.mzML");
  
  // Test Windows-style paths
  TEST_STRING_EQUAL(USI::extractBasename("C:\\data\\sample.mzML"), "sample.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("C:\\Users\\test\\data\\run1.mzML"), "run1.mzML");
  
  // Test file:// URIs (common in mzML and consensusXML)
  TEST_STRING_EQUAL(USI::extractBasename("file:///path/to/sample.mzML"), "sample.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("file:///C:/Users/bielow/data/sample.mzML"), "sample.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("file:///C:/data/experiments/ES-0014b_2.mzML"), "ES-0014b_2.mzML");
  
  // Test basename only (no path)
  TEST_STRING_EQUAL(USI::extractBasename("sample.mzML"), "sample.mzML");
  TEST_STRING_EQUAL(USI::extractBasename("run1.raw"), "run1.raw");
  
  // Test empty string
  TEST_STRING_EQUAL(USI::extractBasename(""), "");
  
  // Test paths with mixed separators
  TEST_STRING_EQUAL(USI::extractBasename("/data/sub\\folder/sample.mzML"), "sample.mzML");
}
END_SECTION

START_SECTION((static const String& getCVAccession()))
{
  TEST_STRING_EQUAL(USI::getCVAccession(), "MS:1003063");
}
END_SECTION

START_SECTION((static const String& getCVName()))
{
  TEST_STRING_EQUAL(USI::getCVName(), "universal spectrum identifier");
}
END_SECTION

START_SECTION((std::ostream& operator<<(std::ostream& os, const USI& usi)))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  stringstream ss;
  ss << usi;
  TEST_STRING_EQUAL(ss.str(), "mzspec:PXD000561:sample.mzML:scan:12345");
}
END_SECTION

// Additional tests for file path handling in USI construction
START_SECTION(([EXTRA] USI construction with file paths))
{
  // Test createFromScanNumber with full path - should use path as-is
  USI usi1 = USI::createFromScanNumber("PXD000561", "/path/to/sample.mzML", 12345);
  TEST_STRING_EQUAL(usi1.getMSRun(), "/path/to/sample.mzML");
  
  // Test creating USI with extracted basename
  String full_path = "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML";
  String basename = USI::extractBasename(full_path);
  USI usi2 = USI::createFromScanNumber("PXD000561", basename, 219);
  TEST_STRING_EQUAL(usi2.getMSRun(), "ES-0014b_2.mzML");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:ES-0014b_2.mzML:scan:219");
  
  // Test with Windows path
  String win_path = "C:\\data\\experiments\\sample.mzML";
  USI usi3 = USI::createFromScanNumber("PXD000561", USI::extractBasename(win_path), 100);
  TEST_STRING_EQUAL(usi3.getMSRun(), "sample.mzML");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
