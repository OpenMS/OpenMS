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

START_SECTION((USI(const std::string& collection, const std::string& ms_run, IndexType index_type, const std::string& index, const std::string& interpretation)))
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

START_SECTION((explicit USI(const std::string& usi_string)))
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

START_SECTION((static bool isValidUSI(const std::string& usi_string)))
{
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345"), true);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2"), true);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2"), true);
  TEST_EQUAL(USI::isValidUSI("invalid_string"), false);
  TEST_EQUAL(USI::isValidUSI("mzspec:PXD000561"), false);
  TEST_EQUAL(USI::isValidUSI(""), false);
}
END_SECTION

START_SECTION((const std::string& getCollection() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000561");
}
END_SECTION

START_SECTION((void setCollection(const std::string& collection)))
{
  USI usi;
  usi.setCollection("PXD000562");
  TEST_STRING_EQUAL(usi.getCollection(), "PXD000562");
}
END_SECTION

START_SECTION((const std::string& getMSRun() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getMSRun(), "sample.mzML");
}
END_SECTION

START_SECTION((void setMSRun(const std::string& ms_run)))
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

START_SECTION((const std::string& getIndex() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345");
  TEST_STRING_EQUAL(usi.getIndex(), "12345");
}
END_SECTION

START_SECTION((void setIndex(const std::string& index)))
{
  USI usi;
  usi.setIndex("99999");
  TEST_STRING_EQUAL(usi.getIndex(), "99999");
}
END_SECTION

START_SECTION((const std::string& getInterpretation() const))
{
  USI usi("PXD000561", "sample.mzML", USI::IndexType::SCAN, "12345", "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi.getInterpretation(), "PEPTIDEK/2");
}
END_SECTION

START_SECTION((void setInterpretation(const std::string& interpretation)))
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

START_SECTION((std::string toString() const))
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

START_SECTION((bool fromString(const std::string& usi_string)))
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

START_SECTION((static USI createFromScanNumber(const std::string& dataset_id, const std::string& filename, int scan_number, const std::string& interpretation)))
{
  USI usi1 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345);
  TEST_STRING_EQUAL(usi1.toString(), "mzspec:PXD000561:sample.mzML:scan:12345");

  USI usi2 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345, "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2");

  USI usi3 = USI::createFromScanNumber("PXD000561", "sample.mzML", 12345, "EM[UNIMOD:35]K/2");
  TEST_STRING_EQUAL(usi3.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2");
}
END_SECTION

START_SECTION((static USI createFromNativeID(const std::string& dataset_id, const std::string& filename, const std::string& native_id)))
{
  // Thermo standard → scan:
  USI usi1 = USI::createFromNativeID("PXD000561", "sample", "scan=12345");
  TEST_STRING_EQUAL(usi1.toString(), "mzspec:PXD000561:sample:scan:12345");

  USI usi2 = USI::createFromNativeID("PXD000561", "sample", "controllerType=0 controllerNumber=1 scan=54321");
  TEST_STRING_EQUAL(usi2.getIndex(), "54321");
  TEST_EQUAL(usi2.getIndexType(), USI::IndexType::SCAN);

  // Bruker TDF → nativeId: comma-separated
  USI usi3 = USI::createFromNativeID("PXD000561", "sample", "frame=120 scan=475");
  TEST_EQUAL(usi3.getIndexType(), USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(usi3.getIndex(), "120,475");

  // index= → index:
  USI usi4 = USI::createFromNativeID("PXD000561", "sample", "index=2074");
  TEST_EQUAL(usi4.getIndexType(), USI::IndexType::INDEX);
  TEST_STRING_EQUAL(usi4.getIndex(), "2074");

  // Unknown format → nativeId: passthrough
  USI usi5 = USI::createFromNativeID("PXD000561", "sample", "custom_id_format");
  TEST_EQUAL(usi5.getIndexType(), USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(usi5.getIndex(), "custom_id_format");

  // Empty → invalid
  USI usi6 = USI::createFromNativeID("PXD000561", "sample", "");
  TEST_EQUAL(usi6.isValid(), false);
}
END_SECTION

START_SECTION((static std::optional<int> extractScanNumberFromNativeID(const std::string& native_id)))
{
  // Only standard Thermo (controllerType=0) and bare scan= should return a scan number
  auto result1 = USI::extractScanNumberFromNativeID("scan=12345");
  TEST_EQUAL(result1.has_value(), true);
  TEST_EQUAL(result1.value(), 12345);

  auto result2 = USI::extractScanNumberFromNativeID("controllerType=0 controllerNumber=1 scan=54321");
  TEST_EQUAL(result2.has_value(), true);
  TEST_EQUAL(result2.value(), 54321);

  // Plain number
  auto result5 = USI::extractScanNumberFromNativeID("12345");
  TEST_EQUAL(result5.has_value(), true);
  TEST_EQUAL(result5.value(), 12345);

  // controllerType=0 but controllerNumber!=1 → not a standard scan (§3.6.4)
  TEST_EQUAL(USI::extractScanNumberFromNativeID("controllerType=0 controllerNumber=2 scan=7").has_value(), false);

  // These should NOT return scan numbers per USI spec
  TEST_EQUAL(USI::extractScanNumberFromNativeID("controllerType=5 controllerNumber=1 scan=7").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("frame=120 scan=475").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("function=10 process=1 scan=345").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("sample=1 period=1 cycle=2740 experiment=10").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("index=100").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("scanId=99999").has_value(), false);
  TEST_EQUAL(USI::extractScanNumberFromNativeID("custom_format_no_number").has_value(), false);
}
END_SECTION

START_SECTION((static std::optional<std::pair<USI::IndexType, std::string>> nativeIdToUSIIndex(const std::string& native_id)))
{
  // Standard Thermo → scan:
  auto r1 = USI::nativeIdToUSIIndex("controllerType=0 controllerNumber=1 scan=43920");
  TEST_EQUAL(r1.has_value(), true);
  TEST_EQUAL(r1->first, USI::IndexType::SCAN);
  TEST_STRING_EQUAL(r1->second, "43920");

  // Bare scan= → scan:
  auto r2 = USI::nativeIdToUSIIndex("scan=12345");
  TEST_EQUAL(r2.has_value(), true);
  TEST_EQUAL(r2->first, USI::IndexType::SCAN);
  TEST_STRING_EQUAL(r2->second, "12345");

  // controllerType=0 but controllerNumber!=1 → nativeId: (§3.6.4)
  auto r2b = USI::nativeIdToUSIIndex("controllerType=0 controllerNumber=2 scan=7");
  TEST_EQUAL(r2b.has_value(), true);
  TEST_EQUAL(r2b->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r2b->second, "0,2,7");

  // Non-standard Thermo controllerType → nativeId: comma-separated
  auto r3 = USI::nativeIdToUSIIndex("controllerType=5 controllerNumber=1 scan=7");
  TEST_EQUAL(r3.has_value(), true);
  TEST_EQUAL(r3->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r3->second, "5,1,7");

  // Bruker TDF → nativeId:
  auto r4 = USI::nativeIdToUSIIndex("frame=120 scan=475");
  TEST_EQUAL(r4.has_value(), true);
  TEST_EQUAL(r4->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r4->second, "120,475");

  // Bruker DIA adds a window-group number; it is not part of the USI tuple.
  auto r4_dia = USI::nativeIdToUSIIndex("frame=120 windowGroup=3 scan=475");
  TEST_EQUAL(r4_dia->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r4_dia->second, "120,475");

  // DDA IDs may carry a precursor number, which is also excluded.
  auto r4_dda = USI::nativeIdToUSIIndex("frame=120 scan=475 precursor=2");
  TEST_EQUAL(r4_dda->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r4_dda->second, "120,475");

  // Waters → nativeId:
  auto r5 = USI::nativeIdToUSIIndex("function=10 process=1 scan=345");
  TEST_EQUAL(r5.has_value(), true);
  TEST_EQUAL(r5->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r5->second, "10,1,345");

  // SCIEX WIFF → nativeId:
  auto r6 = USI::nativeIdToUSIIndex("sample=1 period=1 cycle=2740 experiment=10");
  TEST_EQUAL(r6.has_value(), true);
  TEST_EQUAL(r6->first, USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(r6->second, "1,1,2740,10");

  // index= → index:
  auto r7 = USI::nativeIdToUSIIndex("index=2074");
  TEST_EQUAL(r7.has_value(), true);
  TEST_EQUAL(r7->first, USI::IndexType::INDEX);
  TEST_STRING_EQUAL(r7->second, "2074");

  // Plain number → scan:
  auto r8 = USI::nativeIdToUSIIndex("42");
  TEST_EQUAL(r8.has_value(), true);
  TEST_EQUAL(r8->first, USI::IndexType::SCAN);
  TEST_STRING_EQUAL(r8->second, "42");

  // Empty → nullopt
  TEST_EQUAL(USI::nativeIdToUSIIndex("").has_value(), false);
}
END_SECTION

START_SECTION((static std::string indexTypeToString(IndexType index_type)))
{
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::SCAN), "scan");
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::INDEX), "index");
  TEST_STRING_EQUAL(USI::indexTypeToString(USI::IndexType::NATIVEID), "nativeId");
}
END_SECTION

START_SECTION((static IndexType indexTypeFromString(const std::string& type_string)))
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

START_SECTION((static std::string extractBasename(const std::string& filepath)))
{
  // Test Unix-style paths (extension stripped per USI spec 3.3.4)
  TEST_STRING_EQUAL(USI::extractBasename("/path/to/sample.mzML"), "sample");
  TEST_STRING_EQUAL(USI::extractBasename("/data/experiments/run1.mzML"), "run1");
  TEST_STRING_EQUAL(USI::extractBasename("/sample.mzML"), "sample");

  // Test Windows-style paths
  TEST_STRING_EQUAL(USI::extractBasename("C:\\data\\sample.mzML"), "sample");
  TEST_STRING_EQUAL(USI::extractBasename("C:\\Users\\test\\data\\run1.mzML"), "run1");

  // Test file:// URIs (common in mzML and consensusXML)
  TEST_STRING_EQUAL(USI::extractBasename("file:///path/to/sample.mzML"), "sample");
  TEST_STRING_EQUAL(USI::extractBasename("file:///C:/Users/bielow/data/sample.mzML"), "sample");
  TEST_STRING_EQUAL(USI::extractBasename("file:///C:/data/experiments/ES-0014b_2.mzML"), "ES-0014b_2");

  // Test basename only (no path)
  TEST_STRING_EQUAL(USI::extractBasename("sample.mzML"), "sample");
  TEST_STRING_EQUAL(USI::extractBasename("run1.raw"), "run1");

  // Compound extension (.mzML.gz)
  TEST_STRING_EQUAL(USI::extractBasename("/path/to/sample.mzML.gz"), "sample");

  // Unknown dotted suffix stays intact (e.g. pepXML base_name without extension)
  TEST_STRING_EQUAL(USI::extractBasename("/data/HeLa_1.5ug"), "HeLa_1.5ug");

  // Test without extension
  TEST_STRING_EQUAL(USI::extractBasename("sample"), "sample");

  // Test empty string
  TEST_STRING_EQUAL(USI::extractBasename(""), "");

  // Test paths with mixed separators
  TEST_STRING_EQUAL(USI::extractBasename("/data/sub\\folder/sample.mzML"), "sample");
}
END_SECTION

START_SECTION((static const std::string& getCVAccession()))
{
  TEST_STRING_EQUAL(USI::getCVAccession(), "MS:1003063");
}
END_SECTION

START_SECTION((static const std::string& getCVName()))
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

  // Test creating USI with extracted basename (extension stripped)
  std::string full_path = "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML";
  std::string basename = USI::extractBasename(full_path);
  USI usi2 = USI::createFromScanNumber("PXD000561", basename, 219);
  TEST_STRING_EQUAL(usi2.getMSRun(), "ES-0014b_2");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:ES-0014b_2:scan:219");

  // Test with Windows path
  std::string win_path = "C:\\data\\experiments\\sample.mzML";
  USI usi3 = USI::createFromScanNumber("PXD000561", USI::extractBasename(win_path), 100);
  TEST_STRING_EQUAL(usi3.getMSRun(), "sample");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
