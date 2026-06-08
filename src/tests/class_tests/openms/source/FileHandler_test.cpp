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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <set>

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
TEST_EQUAL(tmp.getTypeByFileName("test.peff"), FileTypes::PEFF)
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
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::ChecksumType::SHA1)

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
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::ChecksumType::SHA1)

// disable hash computation
tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"), exp, {}, ProgressLogger::NONE, true, false);
TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(), "")
TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(), SourceFile::ChecksumType::UNKNOWN_CHECKSUM)
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

// Regression test: computeFileHash must work with non-ASCII paths
{
  namespace fs = std::filesystem;
  // Use u8"" for portable path construction. Only use Latin-1 characters (ä, ü) that are
  // representable in Windows-1252 so that path::string() round-trips correctly on Windows.
  // CJK characters would fail because Windows path::string() uses the Active Code Page.
  fs::path nonascii_dir = fs::path(std::string(File::getTempDirectory())) / u8"openms_t\u00e4st_\u00fc";
  std::error_code ec;
  fs::create_directories(nonascii_dir, ec);
  if (!ec)
  {
    fs::path nonascii_file = nonascii_dir / "test.txt";
    {
      std::ofstream ofs{nonascii_file, std::ios::binary};
      ofs << "hello";
    }
    String hash = FileHandler::computeFileHash(nonascii_file.string());
    TEST_EQUAL(hash.empty(), false) // must succeed, not return ""
    // Clean up
    fs::remove(nonascii_file, ec);
    fs::remove(nonascii_dir, ec);
  }
  else
  {
    STATUS("Skipping non-ASCII path test: could not create directory")
  }
}
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

START_SECTION(([EXTRA] storeIdentifications_loadIdentifications_idparquet_round_trip))
{
  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;

  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  prot_ids.push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(0.5);
  pid.getHits().push_back(hit);
  pep_ids.push_back(pid);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".idparquet";

  FileHandler().storeIdentifications(dir, prot_ids, pep_ids, {FileTypes::IDPARQUET});

  std::vector<ProteinIdentification> prot_ids_in;
  PeptideIdentificationList pep_ids_in;
  FileHandler().loadIdentifications(dir, prot_ids_in, pep_ids_in, {FileTypes::IDPARQUET});

  TEST_EQUAL(prot_ids_in.size(), 1);
  TEST_EQUAL(pep_ids_in.size(), 1);

  // Verify key fields actually round-trip rather than just counting containers.
  // Identifier is synthesized on load per IdXMLFile.cpp:530 parity; stored "run_1"
  // becomes `<search_engine>_<date>_<UniqueIdGenerator>`. Pep_ids re-stamp in lock-step.
  TEST_NOT_EQUAL(prot_ids_in[0].getIdentifier(), "");
  TEST_NOT_EQUAL(prot_ids_in[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(prot_ids_in[0].getScoreType(), "score");
  TEST_EQUAL(prot_ids_in[0].isHigherScoreBetter(), true);

  TEST_STRING_EQUAL(pep_ids_in[0].getIdentifier(), prot_ids_in[0].getIdentifier());
  TEST_STRING_EQUAL(pep_ids_in[0].getScoreType(), "score");
  TEST_EQUAL(pep_ids_in[0].isHigherScoreBetter(), true);
  TEST_EQUAL(pep_ids_in[0].getHits().size(), 1);
  const PeptideHit& h = pep_ids_in[0].getHits()[0];
  TEST_STRING_EQUAL(h.getSequence().toString(), "PEPTIDE");
  TEST_EQUAL(h.getCharge(), 2);
  TEST_REAL_SIMILAR(h.getScore(), 0.5);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] storeFeatures_loadFeatures_featureparquet_round_trip))
{
  FeatureMap fm;
  Feature f;
  f.setUniqueId(42);
  f.setRT(123.4);
  f.setMZ(567.89);
  f.setIntensity(1000.0f);
  f.setCharge(2);
  f.setOverallQuality(0.9f);

  // Minimal ID sidecar: ProteinIdentification on the map, PeptideIdentification
  // on the feature, linked via shared identifier.
  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  fm.getProteinIdentifications().push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(2);
  hit.setScore(0.8);
  pid.getHits().push_back(hit);
  f.getPeptideIdentifications().push_back(pid);

  fm.push_back(f);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".featureparquet";

  FileHandler().storeFeatures(dir, fm, {FileTypes::FEATUREPARQUET});

  FeatureMap fm_in;
  FileHandler().loadFeatures(dir, fm_in, {FileTypes::FEATUREPARQUET});

  TEST_EQUAL(fm_in.size(), 1);
  TEST_EQUAL(fm_in[0].getUniqueId(), 42);
  TEST_REAL_SIMILAR(fm_in[0].getRT(), 123.4);
  TEST_REAL_SIMILAR(fm_in[0].getMZ(), 567.89);
  TEST_REAL_SIMILAR(fm_in[0].getIntensity(), 1000.0f);
  TEST_EQUAL(fm_in[0].getCharge(), 2);
  TEST_REAL_SIMILAR(fm_in[0].getOverallQuality(), 0.9f);

  // ID sidecar round-trip — identifier is synthesized on load (IdXMLFile.cpp:530
  // parity); pep_ids re-stamp in lock-step.
  TEST_EQUAL(fm_in.getProteinIdentifications().size(), 1);
  TEST_NOT_EQUAL(fm_in.getProteinIdentifications()[0].getIdentifier(), "");
  TEST_NOT_EQUAL(fm_in.getProteinIdentifications()[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(fm_in.getProteinIdentifications()[0].getScoreType(), "score");
  TEST_EQUAL(fm_in.getProteinIdentifications()[0].isHigherScoreBetter(), true);
  TEST_EQUAL(fm_in[0].getPeptideIdentifications().size(), 1);
  TEST_STRING_EQUAL(fm_in[0].getPeptideIdentifications()[0].getIdentifier(),
                    fm_in.getProteinIdentifications()[0].getIdentifier());
  TEST_EQUAL(fm_in[0].getPeptideIdentifications()[0].getHits().size(), 1);
  const PeptideHit& h_in = fm_in[0].getPeptideIdentifications()[0].getHits()[0];
  TEST_STRING_EQUAL(h_in.getSequence().toString(), "PEPTIDE");
  TEST_EQUAL(h_in.getCharge(), 2);
  TEST_REAL_SIMILAR(h_in.getScore(), 0.8);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] consensusparquet_round_trip_ProteomicsLFQ_real_output))
{
  // Round-trip a real ProteomicsLFQ consensusXML through .consensusparquet
  // (load -> store as parquet bundle -> load back) and verify the bits that
  // matter for downstream consumers: map references (column headers <->
  // FeatureHandle.map_index), run references (PeptideIdentification.identifier
  // -> ProteinIdentification.identifier), id_merge_index linking PSMs to their
  // original input file via spectra_data, inference info (indistinguishable
  // proteins), and per-PSM hit content.
  ConsensusMap cmap_ref;
  FileHandler().loadConsensusFeatures(
    OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_ProteomicsLFQ_1_subset_out.consensusXML"),
    cmap_ref, {FileTypes::CONSENSUSXML});

  // Sanity-anchor the reference: stable counts for this test fixture.
  TEST_EQUAL(cmap_ref.size(), 39);
  TEST_EQUAL(cmap_ref.getColumnHeaders().size(), 5);
  TEST_EQUAL(cmap_ref.getProteinIdentifications().size(), 1);
  TEST_EQUAL(cmap_ref.getUnassignedPeptideIdentifications().size(), 65);
  TEST_EQUAL(cmap_ref.getProteinIdentifications()[0].getHits().size(), 25);
  TEST_EQUAL(cmap_ref.getProteinIdentifications()[0].getIndistinguishableProteins().size(), 18);
  // ConsensusXMLHandler rebuilds the run identifier as "<engine>_<date>_<hash>"
  // on load — the file-level XML id="PI_0" is only an internal cross-reference.
  TEST_EQUAL(cmap_ref.getProteinIdentifications()[0].getIdentifier().hasPrefix("OMSSA_"), true);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".consensusparquet";

  FileHandler().storeConsensusFeatures(dir, cmap_ref, {FileTypes::CONSENSUSPARQUET});

  ConsensusMap cmap_in;
  FileHandler().loadConsensusFeatures(dir, cmap_in, {FileTypes::CONSENSUSPARQUET});

  // ---- Top-level container counts ----
  TEST_EQUAL(cmap_in.size(), cmap_ref.size());
  TEST_EQUAL(cmap_in.getColumnHeaders().size(), cmap_ref.getColumnHeaders().size());
  TEST_EQUAL(cmap_in.getProteinIdentifications().size(), cmap_ref.getProteinIdentifications().size());
  TEST_EQUAL(cmap_in.getUnassignedPeptideIdentifications().size(), cmap_ref.getUnassignedPeptideIdentifications().size());

  // ---- ColumnHeaders (map references): filename, label, unique_id per map ----
  for (const auto& [map_idx, header_ref] : cmap_ref.getColumnHeaders())
  {
    auto it = cmap_in.getColumnHeaders().find(map_idx);
    TEST_EQUAL(it != cmap_in.getColumnHeaders().end(), true);
    if (it == cmap_in.getColumnHeaders().end()) continue;
    TEST_STRING_EQUAL(it->second.filename, header_ref.filename);
    TEST_STRING_EQUAL(it->second.label, header_ref.label);
    TEST_EQUAL(it->second.unique_id, header_ref.unique_id);
  }

  // ---- FeatureHandle.map_index distribution ----
  // Build per-map_index handle multisets keyed on (uid,map_idx) so any
  // permutation of features still compares equal.
  auto handle_keys = [](const ConsensusMap& m) {
    std::multiset<std::pair<UInt64, UInt64>> keys;
    for (const auto& cf : m)
      for (const auto& fh : cf.getFeatures())
        keys.emplace(fh.getUniqueId(), fh.getMapIndex());
    return keys;
  };
  TEST_EQUAL(handle_keys(cmap_in) == handle_keys(cmap_ref), true);

  // Per-map_index handle counts must match (covers the 5-map distribution).
  auto handles_per_map = [](const ConsensusMap& m) {
    std::map<UInt64, Size> counts;
    for (const auto& cf : m)
      for (const auto& fh : cf.getFeatures())
        ++counts[fh.getMapIndex()];
    return counts;
  };
  TEST_EQUAL(handles_per_map(cmap_in) == handles_per_map(cmap_ref), true);

  // ---- Run references: every PeptideIdentification points to a known run ----
  // Both lanes synthesize fresh identifiers on load (IdXMLFile.cpp:530 parity);
  // identifier SUFFIXES differ between lanes by design, but each lane is internally
  // consistent (no dangling pep_id->prot_id references) and the set sizes match.
  std::set<String> run_ids_ref, run_ids_in;
  for (const auto& p : cmap_ref.getProteinIdentifications()) run_ids_ref.insert(p.getIdentifier());
  for (const auto& p : cmap_in.getProteinIdentifications()) run_ids_in.insert(p.getIdentifier());
  TEST_EQUAL(run_ids_in.size(), run_ids_ref.size());

  Size dangling_ref = 0, dangling_in = 0;
  for (const auto& cf : cmap_ref)
    for (const auto& pid : cf.getPeptideIdentifications())
      if (!run_ids_ref.count(pid.getIdentifier())) ++dangling_ref;
  for (const auto& pid : cmap_ref.getUnassignedPeptideIdentifications())
    if (!run_ids_ref.count(pid.getIdentifier())) ++dangling_ref;
  for (const auto& cf : cmap_in)
    for (const auto& pid : cf.getPeptideIdentifications())
      if (!run_ids_in.count(pid.getIdentifier())) ++dangling_in;
  for (const auto& pid : cmap_in.getUnassignedPeptideIdentifications())
    if (!run_ids_in.count(pid.getIdentifier())) ++dangling_in;
  TEST_EQUAL(dangling_ref, 0);
  TEST_EQUAL(dangling_in, 0);

  // ---- ProteinIdentification (run-level) round-trip ----
  const auto& prot_ref = cmap_ref.getProteinIdentifications()[0];
  const auto& prot_in  = cmap_in.getProteinIdentifications()[0];
  // Identifier prefix (search_engine_date) must match; the UniqueIdGenerator
  // suffix differs between lanes by design.
  TEST_EQUAL(prot_in.getIdentifier().hasPrefix("OMSSA_"), true);
  TEST_EQUAL(prot_ref.getIdentifier().hasPrefix("OMSSA_"), true);
  TEST_STRING_EQUAL(prot_in.getSearchEngine(), prot_ref.getSearchEngine());
  TEST_STRING_EQUAL(prot_in.getSearchEngineVersion(), prot_ref.getSearchEngineVersion());
  TEST_STRING_EQUAL(prot_in.getScoreType(), prot_ref.getScoreType());
  TEST_EQUAL(prot_in.isHigherScoreBetter(), prot_ref.isHigherScoreBetter());
  TEST_EQUAL(prot_in.getHits().size(), prot_ref.getHits().size());

  // spectra_data (the per-merge-index file list)
  StringList spectra_ref, spectra_in;
  prot_ref.getPrimaryMSRunPath(spectra_ref);
  prot_in.getPrimaryMSRunPath(spectra_in);
  TEST_EQUAL(spectra_in == spectra_ref, true);
  TEST_EQUAL(spectra_ref.size(), 5);

  // ProteinHit accession set must round-trip exactly.
  auto accessions = [](const ProteinIdentification& p) {
    std::set<String> acc;
    for (const auto& h : p.getHits()) acc.insert(h.getAccession());
    return acc;
  };
  TEST_EQUAL(accessions(prot_in) == accessions(prot_ref), true);

  // ---- Inference info: IndistinguishableProteins groups ----
  // Compare as a set of (probability, sorted-accession-list) tuples so group
  // ordering doesn't matter.
  auto group_signature = [](const std::vector<ProteinIdentification::ProteinGroup>& gs) {
    std::set<std::pair<double, std::vector<String>>> sig;
    for (const auto& g : gs) {
      std::vector<String> accs(g.accessions.begin(), g.accessions.end());
      std::sort(accs.begin(), accs.end());
      sig.emplace(g.probability, std::move(accs));
    }
    return sig;
  };
  TEST_EQUAL(group_signature(prot_in.getIndistinguishableProteins())
             == group_signature(prot_ref.getIndistinguishableProteins()), true);
  TEST_EQUAL(group_signature(prot_in.getProteinGroups())
             == group_signature(prot_ref.getProteinGroups()), true);

  // ---- id_merge_index distribution (PSM -> original-run linkage) ----
  auto merge_index_hist = [](const ConsensusMap& m) {
    std::map<int, Size> hist;
    auto count = [&](const PeptideIdentification& pid) {
      if (pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
        ++hist[(int)pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)];
      else
        ++hist[-1];
    };
    for (const auto& cf : m)
      for (const auto& pid : cf.getPeptideIdentifications()) count(pid);
    for (const auto& pid : m.getUnassignedPeptideIdentifications()) count(pid);
    return hist;
  };
  auto hist_ref = merge_index_hist(cmap_ref);
  auto hist_in  = merge_index_hist(cmap_in);
  TEST_EQUAL(hist_in == hist_ref, true);
  TEST_EQUAL(hist_ref.count(-1), 0); // every PSM in this fixture must have id_merge_index
  TEST_EQUAL(hist_ref.size(), 5);    // values 0..4 (one per spectra_data entry)

  // ---- PSM hit content fidelity (sequence/charge/score) ----
  // Build (id_merge_index, RT, MZ, sequence, charge, score) tuples for unassigned
  // PSMs (best hit only). RT/MZ make the tuple stable across order. The pid
  // identifier is omitted because both lanes synthesize independently on load
  // (IdXMLFile.cpp:530 parity); each lane is internally consistent (dangling
  // checks above prove that) but the UniqueIdGenerator suffixes differ.
  using PSMSig = std::tuple<int, double, double, String, int, double>;
  auto psm_sigs = [](const auto& pids) {
    std::vector<PSMSig> sigs;
    for (const auto& pid : pids) {
      if (pid.getHits().empty()) continue;
      const auto& h = pid.getHits()[0];
      int mi = pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX)
               ? (int)pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX) : -1;
      sigs.emplace_back(mi, pid.getRT(), pid.getMZ(),
                        h.getSequence().toString(), h.getCharge(), h.getScore());
    }
    std::sort(sigs.begin(), sigs.end());
    return sigs;
  };
  TEST_EQUAL(psm_sigs(cmap_in.getUnassignedPeptideIdentifications())
             == psm_sigs(cmap_ref.getUnassignedPeptideIdentifications()), true);

  // Same comparison for feature-attached PSMs, flattened across all features.
  auto attached_sigs = [&](const ConsensusMap& m) {
    std::vector<PeptideIdentification> all;
    for (const auto& cf : m)
      for (const auto& pid : cf.getPeptideIdentifications()) all.push_back(pid);
    return psm_sigs(all);
  };
  TEST_EQUAL(attached_sigs(cmap_in) == attached_sigs(cmap_ref), true);

  File::removeDirRecursively(dir);
}
END_SECTION

START_SECTION(([EXTRA] storeConsensusFeatures_loadConsensusFeatures_consensusparquet_round_trip))
{
  ConsensusMap cmap;
  ConsensusFeature cf;
  cf.setUniqueId(7);
  cf.setRT(50.0);
  cf.setMZ(300.5);
  cf.setIntensity(2000.0f);
  cf.setCharge(3);
  cf.setQuality(0.8f);

  // Minimal ID sidecar: ProteinIdentification on the map, PeptideIdentification
  // on the consensus feature, linked via shared identifier.
  ProteinIdentification prot;
  prot.setIdentifier("run_1");
  prot.setScoreType("score");
  prot.setHigherScoreBetter(true);
  cmap.getProteinIdentifications().push_back(prot);

  PeptideIdentification pid;
  pid.setIdentifier("run_1");
  pid.setScoreType("score");
  pid.setHigherScoreBetter(true);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setCharge(3);
  hit.setScore(0.7);
  pid.getHits().push_back(hit);
  cf.getPeptideIdentifications().push_back(pid);

  cmap.push_back(cf);

  String dir;
  NEW_TMP_FILE(dir)
  dir += ".consensusparquet";

  FileHandler().storeConsensusFeatures(dir, cmap, {FileTypes::CONSENSUSPARQUET});

  ConsensusMap cmap_in;
  FileHandler().loadConsensusFeatures(dir, cmap_in, {FileTypes::CONSENSUSPARQUET});

  TEST_EQUAL(cmap_in.size(), 1);
  TEST_EQUAL(cmap_in[0].getUniqueId(), 7);
  TEST_REAL_SIMILAR(cmap_in[0].getRT(), 50.0);
  TEST_REAL_SIMILAR(cmap_in[0].getMZ(), 300.5);
  TEST_REAL_SIMILAR(cmap_in[0].getIntensity(), 2000.0f);
  TEST_EQUAL(cmap_in[0].getCharge(), 3);
  TEST_REAL_SIMILAR(cmap_in[0].getQuality(), 0.8f);

  // ID sidecar round-trip — identifier is synthesized on load (IdXMLFile.cpp:530
  // parity); pep_ids re-stamp in lock-step.
  TEST_EQUAL(cmap_in.getProteinIdentifications().size(), 1);
  TEST_NOT_EQUAL(cmap_in.getProteinIdentifications()[0].getIdentifier(), "");
  TEST_NOT_EQUAL(cmap_in.getProteinIdentifications()[0].getIdentifier(), "run_1");
  TEST_STRING_EQUAL(cmap_in.getProteinIdentifications()[0].getScoreType(), "score");
  TEST_EQUAL(cmap_in.getProteinIdentifications()[0].isHigherScoreBetter(), true);
  TEST_EQUAL(cmap_in[0].getPeptideIdentifications().size(), 1);
  TEST_STRING_EQUAL(cmap_in[0].getPeptideIdentifications()[0].getIdentifier(),
                    cmap_in.getProteinIdentifications()[0].getIdentifier());
  TEST_EQUAL(cmap_in[0].getPeptideIdentifications()[0].getHits().size(), 1);
  const PeptideHit& h_in = cmap_in[0].getPeptideIdentifications()[0].getHits()[0];
  TEST_STRING_EQUAL(h_in.getSequence().toString(), "PEPTIDE");
  TEST_EQUAL(h_in.getCharge(), 3);
  TEST_REAL_SIMILAR(h_in.getScore(), 0.7);

  File::removeDirRecursively(dir);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
