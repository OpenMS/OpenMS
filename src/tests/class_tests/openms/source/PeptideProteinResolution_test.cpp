// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeptideProteinResolution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideProteinResolution* ptr = nullptr;
PeptideProteinResolution* null_ptr = nullptr;
START_SECTION(PeptideProteinResolution())
{
	ptr = new PeptideProteinResolution();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(static void PeptideProteinResolution::run(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides))
{
  vector<ProteinIdentification> prots;
  vector<PeptideIdentification> peps;
  IdXMLFile idf;
  idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in.idXML"), prots, peps);  
  PeptideProteinResolution::run(prots, peps);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  IdXMLFile().store(tmp_filename, prots, peps);
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out.idXML"), tmp_filename);

  prots.clear();
  peps.clear();
  tmp_filename.clear();
  NEW_TMP_FILE(tmp_filename);
  idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in2.idXML"), prots, peps);  
  PeptideProteinResolution::run(prots, peps);
  IdXMLFile().store(tmp_filename, prots, peps);
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out2.idXML"), tmp_filename);
}
END_SECTION

START_SECTION(~PeptideProteinResolution())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



