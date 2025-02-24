// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <QtCore/QStringList>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


std::vector<FASTAFile::FASTAEntry> toFASTAVec(const QStringList& sl_prot, const QStringList& identifier = QStringList())
{
  std::vector<FASTAFile::FASTAEntry> proteins;
  for (int i = 0; i < sl_prot.size(); ++i)
  {
    String id = i < identifier.size() ? identifier[i] : String(i); // use identifier if given; or create automatically
    proteins.push_back(FASTAFile::FASTAEntry(id, "", sl_prot[int(i)]));
  }
  return proteins;
}


std::vector<PeptideIdentification> toPepVec(const QStringList& sl_pep)
{
  std::vector<PeptideIdentification> pep_vec;
  for (int i = 0; i < sl_pep.size(); ++i)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(sl_pep[int(i)]));
    std::vector<PeptideHit> hits;
    hits.push_back(hit);
    PeptideIdentification pi;
    pi.setHits(hits);
    pep_vec.push_back(pi);
  }
  return pep_vec;
}

START_TEST(PeptideIndexing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideIndexing* ptr = 0;
PeptideIndexing* null_ptr = 0;
START_SECTION(PeptideIndexing())
{
  ptr = new PeptideIndexing();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~PeptideIndexing())
{
  delete ptr;
}
END_SECTION


START_SECTION((ExitCodes run(std::vector<FASTAFile::FASTAEntry>& proteins, std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids)))
{
  // regression test: https://github.com/OpenMS/OpenMS/issues/3447
  {
    PeptideIndexing indexer;
    Param p = indexer.getParameters();
    p.setValue("decoy_string", "DECOY_");
    indexer.setParameters(p);
    std::vector<FASTAFile::FASTAEntry> proteins = toFASTAVec(QStringList() << "AAAKEEEKTTTK");
    std::vector<ProteinIdentification> prot_ids;
    std::vector<PeptideIdentification> pep_ids = toPepVec(QStringList() << "EEEK(Label:13C(6))");
    indexer.run(proteins, prot_ids, pep_ids);
    TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessionsSet().size(), 1); // one exact hit
    indexer.run(proteins, prot_ids, pep_ids);
    TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessionsSet().size(), 1); // one exact hit
  }

  PeptideIndexing pi;
  Param p = pi.getParameters();
  PeptideIndexing::ExitCodes r;

  // easy case:
  std::vector<FASTAFile::FASTAEntry> proteins = toFASTAVec(QStringList() << "*MLT*EAXK"); // 1 X!!  ; extra * chars (should be ignored)
  std::vector<ProteinIdentification> prot_ids;
  std::vector<PeptideIdentification> pep_ids = toPepVec(QStringList() << "MLTEAEK"); // requires 1 ambAA
  p.setValue("aaa_max", 0);
  p.setValue("decoy_string", "DECOY_");
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessionsSet().size(), 0); // no hit or one hit!
  p.setValue("aaa_max", 1);
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessionsSet().size(), 1); // one hit! -- no ambAA's to spare
  p.setValue("aaa_max", 10);
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessionsSet().size(), 1); // one hit! -- plenty of ambAA's to spare

  // 2 AmbAA's...
  proteins = toFASTAVec(QStringList() << "B*EBE*"); // DB with 2 ambiguous AA's; and extra * chars (should be ignored)
  pep_ids = toPepVec(QStringList() << "NENE" << "NEDE" << "DENE" << "DEDE"); // each is a hit, if >= 2 ambAA's are allowed;
 
  for (int i_aa = 0; i_aa < 5; ++i_aa)
  {
    p.setValue("aaa_max", i_aa);
    pi.setParameters(p);
    std::vector<FASTAFile::FASTAEntry> proteins_local = proteins;
    std::vector<PeptideIdentification> pep_ids_local = pep_ids;
    r = pi.run(proteins_local, prot_ids, pep_ids_local);
    for (Size i = 0; i < pep_ids.size(); ++i)
    {
      set<String> protein_accessions = pep_ids_local[i].getHits()[0].extractProteinAccessionsSet();
      TEST_EQUAL(protein_accessions.size(), i_aa >= 2 ? 1 : 0); // no hit or one hit!
    }
  }
 

  std::cerr << "\n\n testing larger protein with CASIQK...\n\n";
  proteins = toFASTAVec(QStringList() << "SSLDIVLHDTYYVVAHFHYVLSMGAVFAIMGGFIHWFPLFSGYTLDQTYAKIHFTIMFIGVNLTFFPXXXXXXXXXXRRXSDYPDAYTTWNILSSVGSFISLTAVMLMIFXIXEXXASXXKXLMXXXXSXXXXXXXXXXXXXHTFEEPVYMKS");
  //                                                                                                                                                         ^
  //                                   exists      does not exist......                                                                                      CASIQK
  pep_ids = toPepVec(QStringList() << "CASIQK" << "ASIQKFGER" << "KDAVAASIQK" << "KPASIQKR");
  p.setValue("enzyme:specificity", "none");
  p.setValue("missing_decoy_action", "warn");
  pi.setParameters(p);
  for (int i_aa = 0; i_aa < 5; ++i_aa)
  {
    p.setValue("aaa_max", i_aa);
    pi.setParameters(p);
    std::vector<FASTAFile::FASTAEntry> proteins_local = proteins;
    std::vector<PeptideIdentification> pep_ids_local = pep_ids;
    pi.run(proteins_local, prot_ids, pep_ids_local);
    for (Size i = 0; i < pep_ids.size(); ++i)
    {
      set<String> protein_accessions = pep_ids_local[i].getHits()[0].extractProteinAccessionsSet();
      bool is_CASIQK = (i == 0);
      bool allow_at_least_3_ambAA = (i_aa >= 3);
      std::cerr << "TEST: ambAA=" << i_aa << ", hit#:" << i << " ==> prots: " << protein_accessions.size() << "==" << (is_CASIQK & allow_at_least_3_ambAA ? 1 : 0) << "?\n";
      TEST_EQUAL(protein_accessions.size(), is_CASIQK & allow_at_least_3_ambAA ? 1 : 0);
    }
  }

  // empty FASTA (proteins) --> FAIL
  proteins = toFASTAVec(QStringList());
  pep_ids = toPepVec(QStringList() << "SOME" << "PEPTIDES");
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(r, PeptideIndexing::DATABASE_EMPTY);

  // empty idXML (peptides) --> FAIL
  proteins = toFASTAVec(QStringList("PROTEINSEQ"));
  pep_ids = toPepVec(QStringList());
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(r, PeptideIndexing::PEPTIDE_IDS_EMPTY);

  // duplicate accession -- will not be detected and the peptide will have two protein hits.
  // However, extractProteinAccessionsSet() returns a set<>, i.e. only one hit.
  p.setValue("aaa_max", 2);
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBE" << "PROTEIN" << "BEBE", QStringList() << "P_BEBE" << "P_PROTEIN" << "P_BEBE"); //
  pep_ids = toPepVec(QStringList() << "NENE" << "NEDE" << "DENE" << "DEDE"); // 4 hits;
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(proteins.size(), 3) // all three present
   
  // I/L conversion
  p.setValue("aaa_max", 2); // testing I / L conversion, with additional ambAA's to saturate the max_aaa = 2 constraint to ensure that internally 'J' is not used for 'I' or 'L'
  p.setValue("IL_equivalent", "false"); // NOT default
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBEI" << "BEBEL"); //
  pep_ids = toPepVec(QStringList() << "NENEL" << "NEDEL" << "DENEI" << "DEDEI"); // each PSM hits either one or two proteins, depending on I/L setting;
  r = pi.run(proteins, prot_ids, pep_ids);
  for (Size i = 0; i < pep_ids.size(); ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessionsSet().size(), 1); // one hit!
  // ... separate
  p.setValue("IL_equivalent", "true"); // default
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  for (Size i = 0; i < 4; ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessionsSet().size(), 2); // two hits!
  TEST_EQUAL(pep_ids[0].getHits()[0].getSequence().toUnmodifiedString(), "NENEL"); // make sure the PEPTIDE(!) sequence itself is unchanged
  TEST_EQUAL(pep_ids[2].getHits()[0].getSequence().toUnmodifiedString(), "DENEI"); // make sure the PEPTIDE(!) sequence itself is unchanged

  // insertion / deletion
  p.setValue("aaa_max", 2);
  p.setValue("IL_equivalent", "true"); // default
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBE"); //
  pep_ids = toPepVec(QStringList() << "NEKNE" << "NEE"); // 1 insertion, 1 deletion;
  r = pi.run(proteins, prot_ids, pep_ids);
  for (Size i = 0; i < pep_ids.size(); ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessionsSet().size(), 0); // no hits

  // auto mode for decoy strings and position
  std::vector<ProteinIdentification> prot_ids_2;
  std::vector<PeptideIdentification> pep_ids_2;

  {
  // simple prefix
  PeptideIndexing pi_2;
  Param p_2 = pi_2.getParameters();
  std::vector<FASTAFile::FASTAEntry> proteins_2 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "Protein1" << "DECOY_Protein2");
  pi_2.run(proteins_2, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_2.getDecoyString(), "DECOY_");
  TEST_EQUAL(pi_2.isPrefix(), true);
  }

  {
  // simple prefix without special characters
  PeptideIndexing pi_3;
  Param p_3 = pi_3.getParameters();
  std::vector<FASTAFile::FASTAEntry> proteins_3 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "Protein1" << "DECOYProtein2");
  pi_3.run(proteins_3, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_3.getDecoyString(), "DECOY");
  TEST_EQUAL(pi_3.isPrefix(), true);
  }

  {
  // wrong suffix
  PeptideIndexing pi_4;
  Param p_4 = pi_4.getParameters();
  std::vector<FASTAFile::FASTAEntry> proteins_4 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "Protein1" << "Protein2DECOY_");
  pi_4.run(proteins_4, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_4.getDecoyString(), "DECOY_"); //here DECOY_ is the default when finding an affix fails
  TEST_EQUAL(pi_4.isPrefix(), true); // prefix is default too
  }

  {
  // simple suffix
  PeptideIndexing pi_42;
  Param p_42 = pi_42.getParameters();
  std::vector<FASTAFile::FASTAEntry> proteins_42 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "Protein1" << "Protein2_DECOY");
  pi_42.run(proteins_42, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_42.getDecoyString(), "_DECOY");
  TEST_EQUAL(pi_42.isPrefix(), false);
  }

  {
  // complex prefix with one false friend
  PeptideIndexing pi_5;
  Param p_5 = pi_5.getParameters();
  std::vector<FASTAFile::FASTAEntry> proteins_5 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX" << "PEPTLDEXXX" << "PEPTLDEXXX" << "PEPTLDEXXX" << "PEPTLDEXXX",
                                                             QStringList() << "Protein1" << "__id_decoy__Protein2" << "Protein3" <<"Protein4rev" << "__id_decoy__Protein5" << "__id_decoy__Protein6");
  pi_5.run(proteins_5, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_5.getDecoyString(), "__id_decoy__");
  TEST_EQUAL(pi_5.isPrefix(), true);
  }

  {
  // test for self containing decoys: rev vs reverse should output the longer decoy -> reverse?
  PeptideIndexing pi_6;
  std::vector<FASTAFile::FASTAEntry> proteins_6 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "Protein1" << "reverse_Protein");
  pi_6.run(proteins_6, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_6.getDecoyString(), "reverse_");
  TEST_EQUAL(pi_6.isPrefix(), true);
  }

  {
  // impossible to determine automatically -> exit code: DECOYSTRING_EMPTY?
  PeptideIndexing pi_7;
  std::vector<FASTAFile::FASTAEntry> proteins_7 = toFASTAVec(QStringList() << "PEPTIDEXXX" << "PEPTLDEXXX", QStringList() << "rev_Protein1" << "reverse_Protein");
  pi_7.run(proteins_7, prot_ids_2, pep_ids_2);
  TEST_STRING_EQUAL(pi_7.getDecoyString(), "DECOY_");
  TEST_EQUAL(pi_7.isPrefix(), true);
  }

  {
  // test if ambiguous AA's can occur in peptides and are matched without using AAAs or MMs
  PeptideIndexing pi_8;
  Param p_8 = pi_8.getParameters();
  p_8.setValue("aaa_max", 0);
  p_8.setValue("mm_max", 0);
  pi_8.setParameters(p_8);
  std::vector<FASTAFile::FASTAEntry> proteins_8 = toFASTAVec(QStringList() << "PEPTIDERXXXBEBEAR"
                                                                           << "PEPTLDEXXXXBEEEAR",
                                                             QStringList() << "Protein1" << "otherProtein");
  pep_ids = toPepVec(QStringList() << "PEPTIDER"    // matches Protein1;
                                   << "XXXBEBEAR"); // matches Protein1;
  pi_8.run(proteins_8, prot_ids_2, pep_ids);
  for (auto& pep : pep_ids)
  {
    TEST_EQUAL(pep.getHits().size(), 1)
    const auto r = pep.getHits()[0].extractProteinAccessionsSet();
    TEST_EQUAL(r.size(), 1); // one hit!
    TEST_EQUAL(*r.begin(), "Protein1"); // one hit!
  }
  }

  {
  // test no-cleavage (e.g. matching a peptide FASTA DB exactly)
  PeptideIndexing pi_8;
  Param p_8 = pi_8.getParameters();
  p_8.setValue("aaa_max", 0);
  p_8.setValue("mm_max", 0);
  p_8.setValue("enzyme:name", "no cleavage");
  p_8.setValue("allow_nterm_protein_cleavage", "false");
  pi_8.setParameters(p_8);
  std::vector<FASTAFile::FASTAEntry> proteins_8 = toFASTAVec(QStringList() << "MKDPLMMLK"
                                                                           << "KDPLMMLK",
                                                             QStringList() << "Protein1" << "otherProtein");
  pep_ids = toPepVec(QStringList() << "KDPLMMLK" << "MKD"); 
  pi_8.run(proteins_8, prot_ids_2, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits().size(), 1)
  const auto r = pep_ids[0].getHits()[0].extractProteinAccessionsSet();
  TEST_EQUAL(r.size(), 1); // one hit!
  TEST_EQUAL(*r.begin(), "otherProtein"); // one hit!

  TEST_EQUAL(pep_ids[1].getHits()[0].extractProteinAccessionsSet().size(), 0) // no hit for "MKD"
  }

  {
    // test no-cleavage (e.g. matching a peptide FASTA DB exactly) but with ASP/PRO cleavage enabled
    PeptideIndexing pi_8;
    Param p_8 = pi_8.getParameters();
    p_8.setValue("aaa_max", 0);
    p_8.setValue("mm_max", 0);
    p_8.setValue("enzyme:name", "no cleavage");
    p_8.setValue("allow_nterm_protein_cleavage", "false");
    pi_8.setParameters(p_8);
    std::vector<FASTAFile::FASTAEntry> proteins_8 = toFASTAVec(QStringList() << "MKDPLMMLK" // should not hit, due to !allow_nterm_protein_cleavage
                                                                             << "KDPLMMLK", // target
                                                               QStringList() << "Protein1"
                                                                             << "otherProtein");
    prot_ids_2.resize(1);
    prot_ids_2[0].setSearchEngine("XTANDEM"); // enable random ASP/PRO cleavage in PeptideIndexing
    pep_ids = toPepVec(QStringList() << "KDPLMMLK"  // one hit
                                     << "KD");      // one hit due to D/P cleavage
    pi_8.run(proteins_8, prot_ids_2, pep_ids);
    for (auto& pep : pep_ids)
    {
      TEST_EQUAL(pep.getHits().size(), 1)
      const auto r = pep.getHits()[0].extractProteinAccessionsSet();
      TEST_EQUAL(r.size(), 1);                // one hit!
      TEST_EQUAL(*r.begin(), "otherProtein"); // one hit!
    }
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



