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

#include <OpenMS/FORMAT/FileInfo.h>

#include <string>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FileInfo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FileInfo* ptr = nullptr;
FileInfo* null_ptr = nullptr;

START_SECTION((FileInfo()))
  ptr = new FileInfo();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION((~FileInfo()))
  delete ptr;
END_SECTION

START_SECTION((Result run(const std::string& filename, const Options& options)) - featureXML)
{
  FileInfo fi;
  FileInfo::Result r = fi.runAll(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"));
  TEST_EQUAL(r.meta.file_type == FileTypes::FEATUREXML, true)
  TEST_STRING_EQUAL(r.meta.file_type_name, "featureXML")
  TEST_EQUAL(r.feature.has_value(), true)
  TEST_EQUAL(r.feature->is_consensus, false)
  TEST_EQUAL(r.feature->num_features > 0, true)
  TEST_EQUAL(r.ranges.is_experiment, false)
  // the human-readable / TSV renderings are produced and cached
  TEST_EQUAL(r.text.find("-- General information --") != std::string::npos, true)
  TEST_EQUAL(r.tsv.find("general: number of features") != std::string::npos, true)
}
END_SECTION

START_SECTION((Result run ...) - consensusXML)
{
  FileInfo fi;
  FileInfo::Result r = fi.runAll(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.consensusXML"));
  TEST_EQUAL(r.meta.file_type == FileTypes::CONSENSUSXML, true)
  TEST_EQUAL(r.feature.has_value(), true)
  TEST_EQUAL(r.feature->is_consensus, true)
}
END_SECTION

START_SECTION((Result run ...) - mzML peaks)
{
  FileInfo fi;
  FileInfo::Result r = fi.runAll(OPENMS_GET_TEST_DATA_PATH("MzMLFile_2_minimal.mzML"));
  TEST_EQUAL(r.meta.file_type == FileTypes::MZML, true)
  TEST_EQUAL(r.peak.has_value(), true)
  TEST_EQUAL(r.ranges.is_experiment, true)
  // spectra-per-MS-level is consistent with the spectrum count
  UInt64 sum = 0;
  for (const auto& kv : r.peak->spectra_per_ms_level) { sum += kv.second; }
  TEST_EQUAL(sum, r.peak->num_spectra)
}
END_SECTION

START_SECTION((Result run ...) - FASTA)
{
  FileInfo fi;
  FileInfo::Result r = fi.runAll(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(r.meta.file_type == FileTypes::FASTA, true)
  TEST_EQUAL(r.fasta.has_value(), true)
  TEST_EQUAL(r.fasta->num_sequences > 0, true)
  TEST_EQUAL(r.fasta->total_residues > 0, true)
  // length_stats must be populated for any non-empty FASTA (also for 1-2 sequences)
  TEST_EQUAL(r.fasta->length_stats.count, r.fasta->num_sequences)
}
END_SECTION

START_SECTION((static std::string toText(const Result& r)) and toTSV)
{
  FileInfo fi;
  FileInfo::Result r = fi.runAll(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"));
  // static renderers return the cached CLI text / tsv
  TEST_EQUAL(FileInfo::toText(r) == r.text, true)
  TEST_EQUAL(FileInfo::toTSV(r) == r.tsv, true)
  TEST_EQUAL(r.text.empty(), false)
}
END_SECTION

START_SECTION((Options gating))
{
  FileInfo fi;
  const std::string in = OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML");

  FileInfo::Result r_default = fi.run(in); // no -s
  TEST_EQUAL(r_default.text.find("-- Statistics --") != std::string::npos, false)

  FileInfo::Options opt;
  opt.statistics = true; // == CLI -s
  FileInfo::Result r_stats = fi.run(in, opt);
  TEST_EQUAL(r_stats.text.find("-- Statistics --") != std::string::npos, true)
}
END_SECTION

START_SECTION((forced type bypasses auto-detection))
{
  FileInfo fi;
  FileInfo::Options opt;
  opt.forced_type = FileTypes::FEATUREXML;
  FileInfo::Result r = fi.run(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), opt);
  TEST_EQUAL(r.meta.file_type == FileTypes::FEATUREXML, true)
  TEST_EQUAL(r.feature.has_value(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
