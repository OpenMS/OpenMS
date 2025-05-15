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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

START_TEST(FileHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static FileTypes::Type getTypeByFileName(const String &filename)))
FileHandler tmp;
TEST_EQUAL(tmp.getTypeByFileName("test.bla"), FileTypes::UNKNOWN)
TEST_EQUAL(tmp.getTypeByFileName("test.dta"), FileTypes::DTA)
TEST_EQUAL(tmp.getTypeByFileName("test.DTA2D"), FileTypes::DTA2D)
TEST_EQUAL(tmp.getTypeByFileName("test.MzData"), FileTypes::MZDATA)
TEST_EQUAL(tmp.getTypeByFileName("test.MZXML"), FileTypes::MZXML)
TEST_EQUAL(tmp.getTypeByFileName("test.featureXML"), FileTypes::FEATUREXML)
TEST_EQUAL(tmp.getTypeByFileName("test.idXML"), FileTypes::IDXML)
TEST_EQUAL(tmp.getTypeByFileName("test.consensusXML"), FileTypes::CONSENSUSXML)
TEST_EQUAL(tmp.getTypeByFileName("test.mGf"), FileTypes::MGF)
TEST_EQUAL(tmp.getTypeByFileName("test.ini"), FileTypes::INI)
TEST_EQUAL(tmp.getTypeByFileName("test.toPPas"), FileTypes::TOPPAS)
TEST_EQUAL(tmp.getTypeByFileName("test.TraFoXML"), FileTypes::TRANSFORMATIONXML)
TEST_EQUAL(tmp.getTypeByFileName("test.MzML"), FileTypes::MZML)
TEST_EQUAL(tmp.getTypeByFileName(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")), FileTypes::MZML)
TEST_EQUAL(tmp.getTypeByFileName(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")), FileTypes::MZML)
TEST_EQUAL(tmp.getTypeByFileName("test.mS2"), FileTypes::MS2)
TEST_EQUAL(tmp.getTypeByFileName("test.pepXML"), FileTypes::PEPXML)
TEST_EQUAL(tmp.getTypeByFileName("test.pep.xml"), FileTypes::PEPXML)
TEST_EQUAL(tmp.getTypeByFileName("test.protXML"), FileTypes::PROTXML)
TEST_EQUAL(tmp.getTypeByFileName("test.prot.xml"), FileTypes::PROTXML)
TEST_EQUAL(tmp.getTypeByFileName("test.mzid"), FileTypes::MZIDENTML)
TEST_EQUAL(tmp.getTypeByFileName("test.GELML"), FileTypes::GELML)
TEST_EQUAL(tmp.getTypeByFileName("test.TRAML"), FileTypes::TRAML)
TEST_EQUAL(tmp.getTypeByFileName("test.MSP"), FileTypes::MSP)
TEST_EQUAL(tmp.getTypeByFileName("test.OMSSAXML"), FileTypes::OMSSAXML)
TEST_EQUAL(tmp.getTypeByFileName("test.png"), FileTypes::PNG)
TEST_EQUAL(tmp.getTypeByFileName("./foo.bar/XMass/fid"), FileTypes::XMASS);
TEST_EQUAL(tmp.getTypeByFileName("test.TSV"), FileTypes::TSV)
TEST_EQUAL(tmp.getTypeByFileName("test.PEPLIST"), FileTypes::PEPLIST)
TEST_EQUAL(tmp.getTypeByFileName("test.HARDKLOER"), FileTypes::HARDKLOER)
TEST_EQUAL(tmp.getTypeByFileName("test.fasta"), FileTypes::FASTA)
TEST_EQUAL(tmp.getTypeByFileName("test.EDTA"), FileTypes::EDTA)
TEST_EQUAL(tmp.getTypeByFileName("test.csv"), FileTypes::CSV)
TEST_EQUAL(tmp.getTypeByFileName("test.txt"), FileTypes::TXT)
END_SECTION

START_SECTION((static bool hasValidExtension(const String& filename, const FileTypes::Type type)))
TEST_EQUAL(FileHandler::hasValidExtension("test.bla", FileTypes::UNKNOWN), true)
TEST_EQUAL(FileHandler::hasValidExtension("test.idXML", FileTypes::IDXML), true)
TEST_EQUAL(FileHandler::hasValidExtension("test.consensusXML", FileTypes::CONSENSUSXML), true)

// tmp (UNKNOWN)
TEST_EQUAL(FileHandler::hasValidExtension("test.tmp", FileTypes::UNKNOWN), true)
TEST_EQUAL(FileHandler::hasValidExtension("test.tmp", FileTypes::IDXML), true)
TEST_EQUAL(FileHandler::hasValidExtension("test.tmp", FileTypes::CONSENSUSXML), true)

// known other file type
TEST_EQUAL(FileHandler::hasValidExtension("test.consensusXML", FileTypes::IDXML), false)
TEST_EQUAL(FileHandler::hasValidExtension("test.idXML", FileTypes::CONSENSUSXML), false)
END_SECTION

START_SECTION((static FileTypes::Type getTypeByContent(const String &filename)))
  FileHandler tmp;
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData")), FileTypes::MZDATA)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML")), FileTypes::MZXML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML")), FileTypes::MZML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")), FileTypes::MZML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")), FileTypes::MZML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid")), FileTypes::MZIDENTML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML")), FileTypes::FEATUREXML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML")), FileTypes::CONSENSUSXML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML")), FileTypes::IDXML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTAFile_test.dta")), FileTypes::DTA)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d")), FileTypes::DTA2D)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_2.dta2d")), FileTypes::DTA2D)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("class_test_infile.txt")), FileTypes::UNKNOWN)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML")), FileTypes::TRANSFORMATIONXML)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")), FileTypes::FASTA)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FileHandler_toppas.toppas")), FileTypes::TOPPAS)
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FileHandler_MGFbyContent1.mgf")), FileTypes::MGF) // detect via 'FORMAT=Mascot generic\n'
  TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FileHandler_MGFbyContent2.mgf")), FileTypes::MGF) // detect via 'BEGIN IONS\n'

  TEST_EXCEPTION(Exception::FileNotFound, tmp.getTypeByContent("/bli/bla/bluff"))
END_SECTION

START_SECTION((static FileTypes::Type getType(const String &filename)))
  FileHandler tmp;
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("header_file.h")), FileTypes::UNKNOWN)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("class_test_infile.txt")), FileTypes::TXT)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML")), FileTypes::IDXML)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile.consensusXML")), FileTypes::CONSENSUSXML)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML")), FileTypes::TRANSFORMATIONXML)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("FileHandler_toppas.toppas")), FileTypes::TOPPAS)
  TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("pepnovo.txt")), FileTypes::TXT)

  TEST_EXCEPTION(Exception::FileNotFound, tmp.getType("/bli/bla/bluff"))
END_SECTION



START_SECTION((static String stripExtension(const String& file)))
  TEST_STRING_EQUAL(FileHandler::stripExtension(""), "")
  TEST_STRING_EQUAL(FileHandler::stripExtension(".unknown"), "")
  TEST_STRING_EQUAL(FileHandler::stripExtension(".idXML"), "")
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home/doe/file"), "/home/doe/file")
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home/doe/file.txt"), "/home/doe/file")
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home/doe/file.mzML.gz"), "/home/doe/file") // special extension, known to OpenMS
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home/doe/file.txt.tgz"), "/home/doe/file.txt") // not special to us... just strip the last one
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home/doe/file.unknown"), "/home/doe/file")
  TEST_STRING_EQUAL(FileHandler::stripExtension("/home.with.dot/file"), "/home.with.dot/file")
  TEST_STRING_EQUAL(FileHandler::stripExtension("c:\\home.with.dot\\file"), "c:\\home.with.dot\\file")
  TEST_STRING_EQUAL(FileHandler::stripExtension("./filename"), "./filename")
END_SECTION

START_SECTION((static String swapExtension(const String& filename, const FileTypes::Type new_type)))
  TEST_STRING_EQUAL(FileHandler::swapExtension("", FileTypes::UNKNOWN), ".unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension(".unknown", FileTypes::UNKNOWN), ".unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension(".idXML", FileTypes::UNKNOWN), ".unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home/doe/file", FileTypes::UNKNOWN), "/home/doe/file.unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home/doe/file.txt", FileTypes::FEATUREXML), "/home/doe/file.featureXML")
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home/doe/file.mzML.gz", FileTypes::IDXML), "/home/doe/file.idXML") // special extension, known to OpenMS
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home/doe/file.txt.tgz", FileTypes::UNKNOWN), "/home/doe/file.txt.unknown") // not special to us... just strip the last one
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home/doe/file.unknown", FileTypes::UNKNOWN), "/home/doe/file.unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension("/home.with.dot/file", FileTypes::UNKNOWN), "/home.with.dot/file.unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension("c:\\home.with.dot\\file", FileTypes::UNKNOWN), "c:\\home.with.dot\\file.unknown")
  TEST_STRING_EQUAL(FileHandler::swapExtension("./filename", FileTypes::UNKNOWN), "./filename.unknown")
END_SECTION

START_SECTION((FileTypes::Type FileHandler::getConsistentOutputfileType(const String& output_filename, const String& requested_type)))
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("", ""), FileTypes::UNKNOWN)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("a.unknown", "weird"), FileTypes::UNKNOWN)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("a.idXML", ""), FileTypes::IDXML)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.txt", ""), FileTypes::TXT)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.txt", "featureXML"), FileTypes::UNKNOWN)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.mzML.gz", ""), FileTypes::MZML)          // special extension, known to OpenMS
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.mzML.gz", "mzML"), FileTypes::MZML)      // special extension, known to OpenMS
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.txt.tgz", ""), FileTypes::UNKNOWN)       // not special to us... 
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home/doe/file.unknown", "idxml"), FileTypes::IDXML)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("/home.with.dot/file", "mzML"), FileTypes::MZML)
  TEST_EQUAL(FileHandler::getConsistentOutputfileType("c:\\home.with.dot\\file", "mzML"), FileTypes::MZML)
END_SECTION



START_SECTION((template < class PeakType > bool loadExperiment(const String &filename, MSExperiment< PeakType > &exp, FileTypes::Type force_type=FileTypes::UNKNOWN, ProgressLogger::LogType log=ProgressLogger::NONE, const bool compute_hash=true)))
FileHandler tmp;
PeakMap exp;
TEST_EXCEPTION(Exception::FileNotFound, tmp.loadExperiment("test.bla", exp))

tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), exp);
TEST_REAL_SIMILAR(exp[1][0].getPosition()[0], 110)
TEST_REAL_SIMILAR(exp[1][1].getPosition()[0], 120)
TEST_REAL_SIMILAR(exp[1][2].getPosition()[0], 130)

// starts with 110, so this one should skip the first
tmp.getOptions().setMZRange(DRange<1>(115, 1000));
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), exp);
TEST_REAL_SIMILAR(exp[1][0].getPosition()[0], 120)
TEST_REAL_SIMILAR(exp[1][1].getPosition()[0], 130)

tmp.getOptions() = PeakFileOptions();
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"), exp);
TEST_REAL_SIMILAR(exp[2][0].getPosition()[0], 100)
TEST_REAL_SIMILAR(exp[2][1].getPosition()[0], 110)
TEST_REAL_SIMILAR(exp[2][2].getPosition()[0], 120)

tmp.getOptions().setMZRange(DRange<1>(115, 1000));
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"), exp);
TEST_REAL_SIMILAR(exp[2][0].getPosition()[0], 120)
TEST_REAL_SIMILAR(exp[2][1].getPosition()[0], 130)
TEST_REAL_SIMILAR(exp[2][2].getPosition()[0], 140)

tmp.getOptions() = PeakFileOptions();
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp, {FileTypes::MZML}, OpenMS::ProgressLogger::NONE, true, true);
TEST_EQUAL(exp.size(), 4)
TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "36007593dbca0ba59a1f4fc32fb970f0e8991fa6")
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::SHA1)

tmp.getOptions() = PeakFileOptions();
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp);
TEST_REAL_SIMILAR(exp[0][0].getPosition()[0], 230.02)
TEST_REAL_SIMILAR(exp[0][1].getPosition()[0], 430.02)
TEST_REAL_SIMILAR(exp[0][2].getPosition()[0], 630.02)

tmp.getOptions().setMZRange(DRange<1>(300, 1000));
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp, {FileTypes::DTA2D}, OpenMS::ProgressLogger::NONE, true, true);
TEST_REAL_SIMILAR(exp[0][0].getPosition()[0], 430.02)
TEST_REAL_SIMILAR(exp[0][1].getPosition()[0], 630.02)
TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "d50d5144cc3805749b9e8d16f3bc8994979d8142")
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::SHA1)

// disable hash computation
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp, {}, ProgressLogger::NONE, true, false);
TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "")
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::UNKNOWN_CHECKSUM)
// Test that we fail if given a known bogus type restriction
TEST_EXCEPTION(Exception::ParseError, tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp, {FileTypes::SIZE_OF_TYPE}, ProgressLogger::NONE, true, false))

END_SECTION

START_SECTION((static String computeFileHash(const String& filename)))
PeakMap exp;
FileHandler tmp;
// Test that we load with the correct file type restriction
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp, {FileTypes::DTA2D}, ProgressLogger::NONE, true, true);
// compute hash
TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "d50d5144cc3805749b9e8d16f3bc8994979d8142")
END_SECTION


START_SECTION((static bool isSupported(FileTypes::Type type)))
FileHandler tmp;
TEST_EQUAL(false, tmp.isSupported(FileTypes::UNKNOWN));
TEST_EQUAL(true, tmp.isSupported(FileTypes::DTA));
TEST_EQUAL(true, tmp.isSupported(FileTypes::DTA2D));
TEST_EQUAL(true, tmp.isSupported(FileTypes::MZDATA));
TEST_EQUAL(true, tmp.isSupported(FileTypes::MZML));
TEST_EQUAL(true, tmp.isSupported(FileTypes::MZXML));
TEST_EQUAL(true, tmp.isSupported(FileTypes::XMASS));
TEST_EQUAL(true, tmp.isSupported(FileTypes::FEATUREXML));
END_SECTION

START_SECTION((const PeakFileOptions &getOptions() const))
FileHandler a;
TEST_EQUAL(a.getOptions().hasMSLevels(), false)
END_SECTION

START_SECTION((PeakFileOptions & getOptions()))
FileHandler a;
a.getOptions().addMSLevel(1);
TEST_EQUAL(a.getOptions().hasMSLevels(), true);
END_SECTION

START_SECTION((template <class FeatureType> bool loadFeatures(const String &filename, FeatureMap<FeatureType>&map, FileTypes::Type force_type = FileTypes::UNKNOWN)))
FileHandler tmp;
FeatureMap map;
TEST_EXCEPTION(Exception::FileNotFound, tmp.loadFeatures("test.bla", map))
tmp.loadFeatures(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), map);
TEST_EQUAL(map.size(), 7);
tmp.loadFeatures(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), map);
TEST_EQUAL(map.size(), 7);
END_SECTION

START_SECTION((void storeExperiment(const String &filename, const MSExperiment<>&exp, ProgressLogger::LogType log = ProgressLogger::NONE)))
FileHandler fh;
PeakMap exp;
fh.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);

//test mzML
String filename, filename2;
NEW_TMP_FILE_EXT(filename, ".mzML");
fh.storeExperiment(filename, exp, {FileTypes::MZML}, ProgressLogger::NONE);
TEST_EQUAL(fh.getTypeByContent(filename), FileTypes::MZML)

//Test that we throw an exception when given a bogus type restriction
TEST_EXCEPTION(Exception::InvalidFileType, fh.storeExperiment(filename, exp, {FileTypes::SIZE_OF_TYPE}))

//other types cannot be tested, because the NEW_TMP_FILE template does not support file extensions...
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
