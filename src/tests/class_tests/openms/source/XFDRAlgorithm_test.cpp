// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>



///////////////////////////
#include <OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(XFDRAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

XFDRAlgorithm* ptr = 0;
XFDRAlgorithm* null_ptr = 0;
START_SECTION(XFDRAlgorithm())
{
  ptr = new XFDRAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~XFDRAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION(ExitCodes run(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id))

std::vector<PeptideIdentification> peptide_ids;
std::vector<ProteinIdentification> protein_ids;
ProteinIdentification protein_id;

XQuestResultXMLFile xquest_file;
xquest_file.load(OPENMS_GET_TEST_DATA_PATH("XFDRAlgorithm_input.xquest.xml"), peptide_ids, protein_ids);
protein_id = protein_ids[0];


XFDRAlgorithm fdr_algorithm;
Param algo_param = fdr_algorithm.getParameters();
algo_param.setValue("binsize", 0.1);
fdr_algorithm.setParameters(algo_param);

// run algorithm
XFDRAlgorithm::ExitCodes exit_code = fdr_algorithm.run(peptide_ids, protein_id);

TEST_EQUAL(exit_code, XFDRAlgorithm::EXECUTION_OK)
TEST_EQUAL(protein_ids.size(), 1)
TEST_EQUAL(peptide_ids.size(), 310)


for (Size i = 0; i < peptide_ids.size(); i += 30)
{
  auto pep_hits = peptide_ids[i].getHits();
  // the first hit is always the alpha chain
  TEST_EQUAL(pep_hits[0].metaValueExists("xl_target_decoy_alpha"), true)
  TEST_EQUAL(pep_hits[0].metaValueExists("XFDR:FDR"), true)
  TEST_EQUAL(pep_hits[0].metaValueExists("XFDR:used_for_FDR"), true)
  TEST_EQUAL(pep_hits[0].metaValueExists("XFDR:fdr_type"), true)
  if (pep_hits[0].getMetaValue("xl_type") == "cross-link")
  {
    TEST_EQUAL(pep_hits[0].metaValueExists("BetaPepEv:pre"), true)
  }
}

TEST_EQUAL(peptide_ids[50].getHits()[0].getMetaValue("XFDR:FDR"), -0.025)
TEST_EQUAL(peptide_ids[100].getHits()[0].getMetaValue("XFDR:FDR"), 0.934782608695652)
TEST_EQUAL(peptide_ids[250].getHits()[0].getMetaValue("XFDR:FDR"), 0.934782608695652)
TEST_EQUAL(peptide_ids[300].getHits()[0].getMetaValue("XFDR:FDR"), 0.934782608695652)
TEST_EQUAL(peptide_ids[309].getHits()[0].getMetaValue("XFDR:FDR"), -0.025)
TEST_EQUAL(peptide_ids[25].getHits()[0].getMetaValue("XFDR:FDR"), 0.020618556701031)
TEST_EQUAL(peptide_ids[75].getHits()[0].getMetaValue("XFDR:FDR"), 0.934782608695652)
TEST_EQUAL(peptide_ids[275].getHits()[0].getMetaValue("XFDR:FDR"), 0.01063829787234)
TEST_EQUAL(peptide_ids[276].getHits()[0].getMetaValue("XFDR:FDR"), -0.025)

END_SECTION

END_TEST
