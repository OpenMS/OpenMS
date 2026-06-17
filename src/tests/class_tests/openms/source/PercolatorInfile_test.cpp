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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
