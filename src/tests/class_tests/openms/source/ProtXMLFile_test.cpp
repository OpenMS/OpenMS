// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ProtXMLFile.h>
///////////////////////////
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

using namespace OpenMS;
using namespace std;

START_TEST(ProtXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProtXMLFile* ptr = nullptr;
ProtXMLFile* nullPointer = nullptr;
ProtXMLFile file;
START_SECTION(ProtXMLFile())
  ptr = new ProtXMLFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ProtXMLFile())
  delete ptr;
END_SECTION

START_SECTION(void load(const String &filename, ProteinIdentification &protein_ids, PeptideIdentification &peptide_ids))
{
  ProtXMLFile f;
  ProteinIdentification proteins;
  PeptideIdentification peptides;
  String prot_file;

  StringList ids = ListUtils::create<String>("16627578304933075941,13229490167902618598");
  // we do this twice, just to check that members are correctly reset etc..
  for (Int i=0;i<2;++i)
  {
    prot_file = OPENMS_GET_TEST_DATA_PATH("ProtXMLFile_input_1.protXML");
    f.load(prot_file, proteins, peptides);

    // groups
    TEST_EQUAL(proteins.getProteinGroups().size(), 7);
    TEST_EQUAL(proteins.getProteinGroups()[0].probability, 0.9990);
    TEST_EQUAL(proteins.getProteinGroups()[0].accessions.size(), 1);
    TEST_EQUAL(proteins.getProteinGroups()[3].accessions.size(), 2);
    TEST_EQUAL(proteins.getProteinGroups()[3].accessions[0],
               "P01876|IGHA1_HUMAN");
    TEST_EQUAL(proteins.getProteinGroups()[3].accessions[1],
               "P01877|IGHA2_HUMAN");
    TEST_EQUAL(proteins.getProteinGroups()[6].probability, 0.2026);
    TEST_EQUAL(proteins.getProteinGroups()[6].accessions.size(), 1);

    TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 7);
    TEST_EQUAL(proteins.getIndistinguishableProteins()[0].accessions.size(), 1);
    TEST_EQUAL(proteins.getIndistinguishableProteins()[3].accessions.size(), 2);
    TEST_EQUAL(proteins.getIndistinguishableProteins()[3].accessions[0],
               "P01876|IGHA1_HUMAN");
    TEST_EQUAL(proteins.getIndistinguishableProteins()[3].accessions[1],
               "P01877|IGHA2_HUMAN");
    TEST_EQUAL(proteins.getIndistinguishableProteins()[6].accessions.size(), 1);

    // proteins
    TEST_EQUAL(proteins.getHits().size(), 9);
    TEST_EQUAL(proteins.getHits()[0].getAccession(), "P02787|TRFE_HUMAN");
    TEST_EQUAL(proteins.getHits()[0].getCoverage(), 8.6);
    TEST_EQUAL(proteins.getHits()[0].getScore(), 0.9990);
    // this one is indistinguishable... therefore no coverage (but the score
    // got transferred from the "leader" protein):
    TEST_EQUAL(proteins.getHits()[6].getAccession(), "P00739|HPTR_HUMAN");
    TEST_EQUAL(proteins.getHits()[6].getCoverage(), -1);
    TEST_EQUAL(proteins.getHits()[6].getScore(), 0.2663);

    TEST_EQUAL(proteins.getHits()[8].getAccession(), "P04217|A1BG_HUMAN");
    TEST_EQUAL(proteins.getHits()[8].getCoverage(), 2.0);
    TEST_EQUAL(proteins.getHits()[8].getScore(), 0.2026);

    // peptides
    TEST_EQUAL(peptides.getHits().size(), 16);
    AASequence aa_seq = AASequence::fromString("MYLGYEYVTAIR");
    TEST_EQUAL(peptides.getHits()[0].getSequence(), aa_seq);
    TEST_EQUAL(peptides.getHits()[0].getCharge(), 2);
    TEST_EQUAL(peptides.getHits()[0].getScore(), 0.8633);
    set<String> protein_accessions = peptides.getHits()[0].extractProteinAccessionsSet();
    TEST_EQUAL(protein_accessions.size(), 1);
    TEST_EQUAL(*protein_accessions.begin(), "P02787|TRFE_HUMAN");
    TEST_EQUAL(peptides.getHits()[0].getMetaValue("is_unique"), true);
    TEST_EQUAL(peptides.getHits()[0].getMetaValue("is_contributing"), true);

    // load 2 nd file and
    prot_file = OPENMS_GET_TEST_DATA_PATH("ProtXMLFile_input_2.protXML");

  }
}
END_SECTION

START_SECTION(void store(const String &filename, const ProteinIdentification &protein_ids, const PeptideIdentification &peptide_ids, const String &document_id=""))
{
  ProtXMLFile f;
  ProteinIdentification proteins;
  PeptideIdentification peptides;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("notimplemented.protXML", proteins, peptides ))
}
END_SECTION

// Test that probability=0 (subsumable) proteins are filtered out during parsing
// while preserving normal protein groups and indistinguishable proteins.
// Regression test for GitHub issue #6038.
START_SECTION([EXTRA] probability zero protein filtering)
{
  ProtXMLFile f;
  ProteinIdentification proteins;
  PeptideIdentification peptides;
  String prot_file = OPENMS_GET_TEST_DATA_PATH("ProtXMLFile_input_3.protXML");
  f.load(prot_file, proteins, peptides);

  // 3 protein_groups in input, all should still exist
  TEST_EQUAL(proteins.getProteinGroups().size(), 3);

  // Group 1: leader only (subsumable protein was filtered)
  TEST_EQUAL(proteins.getProteinGroups()[0].accessions.size(), 1);
  TEST_EQUAL(proteins.getProteinGroups()[0].accessions[0], "LEADER_PROT");
  TEST_REAL_SIMILAR(proteins.getProteinGroups()[0].probability, 0.9998);

  // Group 2: protein + indistinguishable sibling (unaffected)
  TEST_EQUAL(proteins.getProteinGroups()[1].accessions.size(), 2);
  TEST_EQUAL(proteins.getProteinGroups()[1].accessions[0], "PROT_A");
  TEST_EQUAL(proteins.getProteinGroups()[1].accessions[1], "PROT_B");

  // Group 3: single protein (unaffected)
  TEST_EQUAL(proteins.getProteinGroups()[2].accessions.size(), 1);
  TEST_EQUAL(proteins.getProteinGroups()[2].accessions[0], "SINGLE_PROT");

  // 3 indistinguishable groups (NOT 4 -- the subsumable protein should not create one)
  TEST_EQUAL(proteins.getIndistinguishableProteins().size(), 3);

  // Indist group 1: leader only
  TEST_EQUAL(proteins.getIndistinguishableProteins()[0].accessions.size(), 1);
  TEST_EQUAL(proteins.getIndistinguishableProteins()[0].accessions[0], "LEADER_PROT");

  // Indist group 2: PROT_A + PROT_B
  TEST_EQUAL(proteins.getIndistinguishableProteins()[1].accessions.size(), 2);
  TEST_EQUAL(proteins.getIndistinguishableProteins()[1].accessions[0], "PROT_A");
  TEST_EQUAL(proteins.getIndistinguishableProteins()[1].accessions[1], "PROT_B");

  // Indist group 3: single
  TEST_EQUAL(proteins.getIndistinguishableProteins()[2].accessions.size(), 1);
  TEST_EQUAL(proteins.getIndistinguishableProteins()[2].accessions[0], "SINGLE_PROT");

  // 4 protein hits (LEADER_PROT, PROT_A, PROT_B, SINGLE_PROT -- no SUBSUMABLE_PROT)
  TEST_EQUAL(proteins.getHits().size(), 4);
  TEST_EQUAL(proteins.getHits()[0].getAccession(), "LEADER_PROT");
  TEST_REAL_SIMILAR(proteins.getHits()[0].getScore(), 0.9998);
  TEST_EQUAL(proteins.getHits()[1].getAccession(), "PROT_A");
  TEST_REAL_SIMILAR(proteins.getHits()[1].getScore(), 0.9000);
  TEST_EQUAL(proteins.getHits()[2].getAccession(), "PROT_B");
  TEST_REAL_SIMILAR(proteins.getHits()[2].getScore(), 0.9000); // inherited from leader
  TEST_EQUAL(proteins.getHits()[3].getAccession(), "SINGLE_PROT");
  TEST_REAL_SIMILAR(proteins.getHits()[3].getScore(), 0.8000);

  // 4 peptide hits (2 from leader, 1 from group 2, 1 from group 3 -- NOT the subsumable's peptide)
  TEST_EQUAL(peptides.getHits().size(), 4);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

