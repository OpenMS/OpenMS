// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>


using namespace OpenMS;
using namespace std;

START_TEST(MzIdentMLFile, "$Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MzIdentMLFile* ptr = 0;
MzIdentMLFile* nullPointer = 0;
START_SECTION((MzIdentMLFile()))
    ptr = new MzIdentMLFile;
    TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MzIdentMLFile()))
    delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids) )
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),1)
  TEST_EQUAL(peptide_ids.size(),3)

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getScoreType(),"Mascot:score")
  TEST_EQUAL(protein_ids[0].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[0].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),"2.1.0")
  TEST_EQUAL(protein_ids[0].getDateTime().getDate(),"2006-01-12")
  TEST_EQUAL(protein_ids[0].getDateTime().getTime(),"00:00:00")
  TEST_EQUAL(protein_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
  TEST_EQUAL(protein_ids[0].getSearchParameters().enzyme,ProteinIdentification::TRYPSIN)
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"+1, +2")
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::AVERAGE)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().peak_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_tolerance,1.0)
  TEST_EQUAL((String)(protein_ids[0].getMetaValue("name")),"ProteinIdentification")

  //ProteinGroups not nupported yet
  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 0);
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 0);

  TEST_EQUAL(protein_ids[0].getHits().size(),2)
  //protein hit 1
  TEST_REAL_SIMILAR(protein_ids[0].getHits()[0].getScore(),34.4)
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),"PROT1")
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),"ABCDEFG")
  TEST_EQUAL((String)(protein_ids[0].getHits()[0].getMetaValue("name")),"ProteinHit")
  //protein hit 2
  TEST_REAL_SIMILAR(protein_ids[0].getHits()[1].getScore(),24.4)
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),"PROT2")
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),"ABCDEFG")

  //peptide id 1
  TEST_EQUAL(peptide_ids[0].getScoreType(),"Mascot:score")
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),false)
  TEST_EQUAL(peptide_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),675.9)
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),1234.5)
  //TEST_EQUAL((UInt)(peptide_ids[0].getMetaValue("spectrum_reference")),17)
  TEST_EQUAL((String)(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
  TEST_EQUAL(peptide_ids[0].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),0.9)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(), AASequence::fromString("PEPTIDER"))
  TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),1)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getAABefore(),'A')
  TEST_EQUAL(peptide_ids[0].getHits()[0].getAAAfter(),'B')
  TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions().size(),2)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[0],"PROT1")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[1],"PROT2")
  TEST_EQUAL((String)(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERR"))
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
  TEST_EQUAL(peptide_ids[0].getHits()[1].getProteinAccessions().size(),0)
  TEST_EQUAL(peptide_ids[0].getHits()[1].getAABefore(),' ')
  TEST_EQUAL(peptide_ids[0].getHits()[1].getAAAfter(),' ')

  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),"Mascot:score")
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(peptide_ids[1].getIdentifier(),"Mascot_2006-01-12T12:13:14")
  TEST_EQUAL(peptide_ids[1].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),44.4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),2)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getProteinAccessions().size(),0)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getAABefore(),' ')
  TEST_EQUAL(peptide_ids[1].getHits()[0].getAAAfter(),' ')

  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),33.3)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERRRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[1].getCharge(),2)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getProteinAccessions().size(),0)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getAABefore(),' ')
  TEST_EQUAL(peptide_ids[1].getHits()[1].getAAAfter(),' ')

  /////////////// protein id 2 //////////////////
  TEST_EQUAL(protein_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(protein_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[1].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[1].getSearchEngineVersion(),"2.1.1")
  TEST_EQUAL(protein_ids[1].getDateTime().getDate(),"2007-01-12")
  TEST_EQUAL(protein_ids[1].getDateTime().getTime(),"12:13:14")
  TEST_EQUAL(protein_ids[1].getIdentifier(),"Mascot_2007-01-12T12:13:14")
  TEST_EQUAL(protein_ids[1].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[1].getSearchParameters().db_version,"1.1")
  TEST_EQUAL(protein_ids[1].getSearchParameters().enzyme,ProteinIdentification::UNKNOWN_ENZYME)
  TEST_EQUAL(protein_ids[1].getSearchParameters().charges,"+1, +2, +3")
  TEST_EQUAL(protein_ids[1].getSearchParameters().mass_type,ProteinIdentification::MONOISOTOPIC)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().peak_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().precursor_tolerance,1.0)
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications.size(),2)
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications[0],"Fixed")
  TEST_EQUAL(protein_ids[1].getSearchParameters().fixed_modifications[1],"Fixed2")
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications.size(),2)
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications[0],"Variable")
  TEST_EQUAL(protein_ids[1].getSearchParameters().variable_modifications[1],"Variable2")
  TEST_EQUAL(protein_ids[1].getHits().size(),1)
  //protein hit 1
  TEST_REAL_SIMILAR(protein_ids[1].getHits()[0].getScore(),100.0)
  TEST_EQUAL(protein_ids[1].getHits()[0].getAccession(),"PROT3")
  TEST_EQUAL(protein_ids[1].getHits()[0].getSequence(),"")
  //peptide id 3
  TEST_EQUAL(peptide_ids[2].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[2].isHigherScoreBetter(),true)
  TEST_EQUAL(peptide_ids[2].getIdentifier(),"Mascot_2007-01-12T12:13:14")
  TEST_EQUAL(peptide_ids[2].getHits().size(),1)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),1.4)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRRRR"))
  TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),1)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions().size(),1)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions()[0],"PROT3")
  TEST_EQUAL(peptide_ids[2].getHits()[0].getAABefore(),' ')
  TEST_EQUAL(peptide_ids[2].getHits()[0].getAAAfter(),' ')
END_SECTION

START_SECTION(void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids) )

  //store and load data
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String document_id, document_id2;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  FuzzyStringComparator fuzzy;
  fuzzy.setWhitelist(ListUtils::create<String>("<?xml-stylesheet"));
  fuzzy.setAcceptableAbsolute(0.0001);
  bool result = fuzzy.compareFiles(input_path, filename);
  TEST_EQUAL(result, true);
END_SECTION

START_SECTION(bool isValid(const String& filename, std::ostream& os = std::cerr))
  MzIdentMLFile file;
    // TODO
END_SECTION

START_SECTION(bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String filename;
  StringList errs,warns;
  MzIdentMLFile f;

  //test if empty file is valid
  NEW_TMP_FILE(filename)
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isSemanticallyValid(filename, errs, warns), true);

  //test if full file is valid
  NEW_TMP_FILE(filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid"), protein_ids2, peptide_ids2);
  protein_ids2[0].setMetaValue("stringvalue",String("bla"));
  protein_ids2[0].setMetaValue("intvalue",4711);
  protein_ids2[0].setMetaValue("floatvalue",5.3);
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isSemanticallyValid(filename, errs, warns), true);

END_SECTION

START_SECTION(([EXTRA] No protein identification bug))
  MzIdentMLFile mzidfile;
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  mzidfile.load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_no_proteinhits.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(), 1)
  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
  TEST_EQUAL(peptide_ids.size(), 10)
  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

  String filename;
  NEW_TMP_FILE(filename)
  mzidfile.store(filename , protein_ids, peptide_ids);

  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;
  mzidfile.load(filename, protein_ids2, peptide_ids2);

  TEST_EQUAL(protein_ids == protein_ids2, true)
  TEST_EQUAL(peptide_ids == peptide_ids2, true)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
