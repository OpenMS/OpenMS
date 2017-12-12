// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <Qt/qstringlist.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


std::vector<FASTAFile::FASTAEntry> toFASTAVec(const QStringList& sl_prot, const QStringList& identifier = QStringList())
{
  std::vector<FASTAFile::FASTAEntry> proteins;
  for (Size i = 0; i < sl_prot.size(); ++i)
  {
    String id = i < identifier.size() ? identifier[i] : String(i); // use identifier if given; or create automatically
    proteins.push_back(FASTAFile::FASTAEntry(id, "", sl_prot[int(i)]));
  }
  return proteins;
}


std::vector<PeptideIdentification> toPepVec(const QStringList& sl_pep)
{
  std::vector<PeptideIdentification> pep_vec;
  for (Size i = 0; i < sl_pep.size(); ++i)
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
  PeptideIndexing pi;
  Param p = pi.getParameters();
  PeptideIndexing::ExitCodes r;

  // easy case:
  std::vector<FASTAFile::FASTAEntry> proteins = toFASTAVec(QStringList() << "*MLT*EAXK"); // 1 X!!  ; extra * chars (should be ignored)
  std::vector<ProteinIdentification> prot_ids;
  std::vector<PeptideIdentification> pep_ids = toPepVec(QStringList() << "MLTEAEK"); // requires 1 ambAA
  p.setValue("aaa_max", 0);
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessions().size(), 0); // no hit or one hit!
  p.setValue("aaa_max", 1);
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessions().size(), 1); // one hit! -- no ambAA's to spare
  p.setValue("aaa_max", 10);
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(pep_ids[0].getHits()[0].extractProteinAccessions().size(), 1); // one hit! -- plenty of ambAA's to spare

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
      set<String> protein_accessions = pep_ids_local[i].getHits()[0].extractProteinAccessions();
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
    PeptideIndexing::ExitCodes r = pi.run(proteins_local, prot_ids, pep_ids_local);
    for (Size i = 0; i < pep_ids.size(); ++i)
    {
      set<String> protein_accessions = pep_ids_local[i].getHits()[0].extractProteinAccessions();
      bool is_CASIQK = (i == 0);
      bool allow_at_least_3_ambAA = (i_aa >= 3);
      std::cerr << "TEST: ambAA=" << i_aa << ", hit#:" << i << " ==> prots: " << protein_accessions.size() << "==" << (is_CASIQK & allow_at_least_3_ambAA ? 1 : 0) << "?\n";
      TEST_EQUAL(protein_accessions.size(), is_CASIQK & allow_at_least_3_ambAA ? 1 : 0);
    }
  }

  // empty FASTA (proteins --> FAIL
  proteins = toFASTAVec(QStringList());
  pep_ids = toPepVec(QStringList() << "SOME" << "PEPTIDES");
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(r, PeptideIndexing::DATABASE_EMPTY);

  // empty idXML (peptides) --> FAIL
  proteins = toFASTAVec(QStringList("PROTEINSEQ"));
  pep_ids = toPepVec(QStringList());
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(r, PeptideIndexing::PEPTIDE_IDS_EMPTY);

  // duplicate accession -- just a warning, nothing else; second protein will be removed
  p.setValue("aaa_max", 2);
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBE" << "PROTEIN" << "BEBE", QStringList() << "P_BEBE" << "P_PROTEIN" << "P_BEBE"); //
  pep_ids = toPepVec(QStringList() << "NENE" << "NEDE" << "DENE" << "DEDE"); // 4 hits;
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(proteins.size(), 2) // one removed!
  for (int i = 0; i < pep_ids.size(); ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessions().size(), 1); // one hit!
  // ... however, if sequences are not equal: bail out
  proteins = toFASTAVec(QStringList() << "BEBE" << "PROTEIN" << "NOT*BEBE", QStringList() << "P_BEBE" << "P_PROTEIN" << "P_BEBE"); //
  pep_ids = toPepVec(QStringList() << "NENE" << "NEDE" << "DENE" << "DEDE"); // 4 hits;
  r = pi.run(proteins, prot_ids, pep_ids);
  TEST_EQUAL(r, PeptideIndexing::DATABASE_CONTAINS_MULTIPLES)
   
  // I/L conversion
  p.setValue("aaa_max", 2); // testing I / L conversion, with additional ambAA's to saturate the max_aaa = 2 constraint to ensure that internally 'J' is not used for 'I' or 'L', since 'J' is unknown to SeqAn and will get converted to 'X' 
  p.setValue("IL_equivalent", "false"); // NOT default
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBEI" << "BEBEL"); //
  pep_ids = toPepVec(QStringList() << "NENEL" << "NEDEL" << "DENEI" << "DEDEI"); // each PSM hits either one or two proteins, depending on I/L setting;
  r = pi.run(proteins, prot_ids, pep_ids);
  for (int i = 0; i < pep_ids.size(); ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessions().size(), 1); // one hit!
  // ... separate
  p.setValue("IL_equivalent", "true"); // default
  pi.setParameters(p);
  r = pi.run(proteins, prot_ids, pep_ids);
  for (int i = 0; i < 4; ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessions().size(), 2); // two hits!
  TEST_EQUAL(pep_ids[0].getHits()[0].getSequence().toUnmodifiedString(), "NENEL"); // make sure the PEPTIDE(!) sequence itself is unchanged
  TEST_EQUAL(pep_ids[2].getHits()[0].getSequence().toUnmodifiedString(), "DENEI"); // make sure the PEPTIDE(!) sequence itself is unchanged

  // insertion / deletion
  p.setValue("aaa_max", 2);
  p.setValue("IL_equivalent", "true"); // default
  pi.setParameters(p);
  proteins = toFASTAVec(QStringList() << "BEBE"); //
  pep_ids = toPepVec(QStringList() << "NEKNE" << "NEE"); // 1 insertion, 1 deletion;
  r = pi.run(proteins, prot_ids, pep_ids);
  for (int i = 0; i < pep_ids.size(); ++i) TEST_EQUAL(pep_ids[i].getHits()[0].extractProteinAccessions().size(), 0); // no hits

 

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



