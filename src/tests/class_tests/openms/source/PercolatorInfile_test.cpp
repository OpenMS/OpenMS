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
#include <OpenMS/FORMAT/PercolatorInfile.h>
///////////////////////////

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <fstream>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(PercolatorInfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PercolatorInfile* ptr = nullptr;
PercolatorInfile* null_pointer = nullptr;

START_SECTION(PercolatorInfile())
{
  ptr = new PercolatorInfile();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~PercolatorInfile())
{
  delete ptr;
}
END_SECTION

START_SECTION(PeptideIdentificationList PercolatorInfile::load(const std::string& pin_file, bool higher_score_better, const std::string& score_name, std::string decoy_prefix))
{
  StringList filenames;
  // test loading of pin file with automatic update of target/decoy annotation based on decoy prefix in protein accessions

  // test some extra scores
  StringList extra_scores = {"ln(delta_next)", "ln(delta_best)", "matched_peaks"};

  auto pids = PercolatorInfile::load(OPENMS_GET_TEST_DATA_PATH("sage.pin"), 
    true, 
    "ln(hyperscore)", 
    extra_scores, 
    filenames, 
    "DECOY_");
  TEST_EQUAL(pids.size(), 9)
  TEST_EQUAL(filenames.size(), 2)
  TEST_EQUAL(pids[0].getSpectrumReference(), "30381")
  TEST_EQUAL(pids[6].getSpectrumReference(), "spectrum=2041")
  TEST_EQUAL(pids[7].getHits()[0].getMetaValue("target_decoy"),"decoy") // 8th entry is annotated as target in pin file but only maps to decoy proteins with prefix "DECOY_" -> set to decoy
}
END_SECTION

START_SECTION((static StringList getStandardFeatureSet(int min_charge, int max_charge)))
{
  // The standard Percolator feature set is the three mandatory header columns
  // (SpecId, Label, ScanNr), then the mass/length features, then one "charge<c>"
  // column per charge state in [min,max], then the enzyme / mass-delta features.
  // Callers append any extra_features and finally Peptide, Proteins before writing
  // the .pin header that Percolator expects (#9195). This validates that contract.

  // exact ordered feature set for charges 2..4: three mandatory columns, then
  // mass/length features, then one "charge<c>" column per charge state, then the
  // enzyme / mass-delta features.
  StringList fs = PercolatorInfile::getStandardFeatureSet(2, 4);
  StringList expected = {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen",
                         "charge2", "charge3", "charge4",
                         "enzN", "enzC", "enzInt", "dm", "absdm"};
  TEST_EQUAL(fs.size(), expected.size())
  for (Size i = 0; i < expected.size(); ++i)
  {
    TEST_STRING_EQUAL(fs[i], expected[i])
  }

  // the three mandatory leading columns
  TEST_STRING_EQUAL(fs[0], "SpecId")
  TEST_STRING_EQUAL(fs[1], "Label")
  TEST_STRING_EQUAL(fs[2], "ScanNr")

  // a single charge state yields exactly one "charge3" column and no other charge column
  StringList single = PercolatorInfile::getStandardFeatureSet(3, 3);
  StringList expected_single = {"SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "mass", "peplen",
                                "charge3", "enzN", "enzC", "enzInt", "dm", "absdm"};
  TEST_EQUAL(single.size(), expected_single.size())
  for (Size i = 0; i < expected_single.size(); ++i)
  {
    TEST_STRING_EQUAL(single[i], expected_single[i])
  }

  // the assembled .pin header (features + Peptide + Proteins) begins with SpecId
  // and ends with Peptide, Proteins -- the exact shape Percolator parses (#9195).
  StringList header = fs;
  header.push_back("Peptide");
  header.push_back("Proteins");
  TEST_STRING_EQUAL(header.front(), "SpecId")
  TEST_STRING_EQUAL(header[header.size() - 2], "Peptide")
  TEST_STRING_EQUAL(header.back(), "Proteins")
}
END_SECTION

START_SECTION((static void store(const std::string& pin_file, const PeptideIdentificationList& peptide_ids, const StringList& feature_set, const std::string& enz, int min_charge, int max_charge)))
{
  // Validates the -out_pin output that PercolatorAdapter/ProSE/OpenNuXL emit:
  // the header line must be the exact column contract Percolator parses
  // (SpecId Label ScanNr ... Peptide Proteins) and store() must compute and
  // write the standard per-PSM features (stampPinFeaturesOnHits). Build one
  // fully-annotated target PSM so it is not skipped and a data row is written.
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("SAMPLER"));
  hit.setCharge(2);
  hit.setScore(1.0);
  hit.setTargetDecoyType(PeptideHit::TargetDecoyType::TARGET);
  PeptideEvidence ev;
  ev.setProteinAccession("PROT1");
  ev.setAABefore('K'); // tryptic flanks so enzN/enzC are true
  ev.setAAAfter('S');
  hit.setPeptideEvidences(std::vector<PeptideEvidence>{ev});

  PeptideIdentification pid;
  pid.setMZ(500.25);
  pid.setRT(123.4);
  pid.setSpectrumReference("scan=529");
  pid.setHits(std::vector<PeptideHit>{hit});

  PeptideIdentificationList pids;
  pids.push_back(pid);

  // the real .pin column contract: standard features + Peptide + Proteins
  StringList feature_set = PercolatorInfile::getStandardFeatureSet(2, 3);
  feature_set.push_back("Peptide");
  feature_set.push_back("Proteins");

  std::string pin_file;
  NEW_TMP_FILE(pin_file);
  PercolatorInfile::store(pin_file, pids, feature_set, "trypsin", 2, 3);

  // read the written .pin back
  std::ifstream is(pin_file.c_str());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(is, line)) lines.push_back(line);

  // header + exactly one data row (the PSM must not be skipped)
  TEST_EQUAL(lines.size(), 2)
  ABORT_IF(lines.empty())

  // --- header: column names and order Percolator expects ---
  StringList header;
  StringUtils::split(lines[0], '\t', header);
  TEST_EQUAL(header.size(), feature_set.size())
  TEST_STRING_EQUAL(header[0], "SpecId")
  TEST_STRING_EQUAL(header[1], "Label")
  TEST_STRING_EQUAL(header[2], "ScanNr")
  TEST_STRING_EQUAL(header[header.size() - 2], "Peptide")
  TEST_STRING_EQUAL(header.back(), "Proteins")

  auto has = [&header](const std::string& n) {
    return std::find(header.begin(), header.end(), n) != header.end();
  };
  // mass + enzyme feature columns are present
  TEST_EQUAL(has("ExpMass"), true)
  TEST_EQUAL(has("CalcMass"), true)
  TEST_EQUAL(has("mass"), true)
  TEST_EQUAL(has("peplen"), true)
  TEST_EQUAL(has("charge2"), true)
  TEST_EQUAL(has("charge3"), true)
  TEST_EQUAL(has("enzN"), true)
  TEST_EQUAL(has("enzC"), true)
  TEST_EQUAL(has("enzInt"), true)
  TEST_EQUAL(has("dm"), true)
  TEST_EQUAL(has("absdm"), true)

  // --- data row: same column count + the deterministic computed feature values ---
  ABORT_IF(lines.size() < 2)
  StringList row;
  StringUtils::split(lines[1], '\t', row);
  TEST_EQUAL(row.size(), header.size())

  auto col = [&header](const std::string& n) -> Size {
    return (Size)(std::find(header.begin(), header.end(), n) - header.begin());
  };
  TEST_STRING_EQUAL(row[col("SpecId")], "scan=529")
  TEST_STRING_EQUAL(row[col("ScanNr")], "529")
  TEST_STRING_EQUAL(row[col("Label")], "1")          // target -> 1
  TEST_STRING_EQUAL(row[col("peplen")], "7")         // SAMPLER
  TEST_STRING_EQUAL(row[col("charge2")], "1")        // charge == 2 (one-hot)
  TEST_STRING_EQUAL(row[col("charge3")], "0")
  TEST_STRING_EQUAL(row[col("enzN")], "1")           // tryptic N-terminus
  TEST_STRING_EQUAL(row[col("enzC")], "1")           // tryptic C-terminus
  TEST_STRING_EQUAL(row[col("Peptide")], "K.SAMPLER.S")
  TEST_STRING_EQUAL(row[col("Proteins")], "PROT1")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
