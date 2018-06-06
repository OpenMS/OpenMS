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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

