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

START_SECTION(vector<PeptideIdentification> PercolatorInfile::load(const String& pin_file, bool higher_score_better, const String& score_name, String decoy_prefix))
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
