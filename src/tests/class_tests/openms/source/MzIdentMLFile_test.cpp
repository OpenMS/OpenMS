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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>


using namespace OpenMS;
using namespace std;

START_TEST(MzIdentMLFile, "$Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MzIdentMLFile* ptr = nullptr;
MzIdentMLFile* nullPointer = nullptr;
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
  std::vector<String> fm;
  fm.push_back("Carbamidomethyl (C)");
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(protein_ids[0].getHits().size(),2)
  TEST_EQUAL(protein_ids[1].getHits().size(),1)
  TEST_EQUAL(peptide_ids.size(),5)
  TEST_EQUAL(peptide_ids[0].getHits().size(),1)
  TEST_EQUAL(peptide_ids[1].getHits().size(),1)
  TEST_EQUAL(peptide_ids[2].getHits().size(),1)
  TEST_EQUAL(peptide_ids[3].getHits().size(),1)
  TEST_EQUAL(peptide_ids[4].getHits().size(),1)

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getSearchEngine(),"MS-GF+")
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),"Beta (v9979)")
  TEST_NOT_EQUAL(protein_ids[0].getDateTime().getDate(),"0000-00-00")
  TEST_NOT_EQUAL(protein_ids[0].getDateTime().getTime(),"00:00:00")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"database.fasta")
  TEST_EQUAL(protein_ids[0].getSearchParameters().missed_cleavages, 1000)
  TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications.size(), fm.size())
  TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications.back(), fm.back())
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,0)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,20)

  //ProteinGroups not nupported yet, also no ProteinDetection, too few input here
//  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 0);
//  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 0);

  //protein hit 1
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),"sp|P0A9K9|SLYD_ECOLI")
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),"")
  //protein hit 2
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),"sp|P0A786|PYRB_ECOLI")
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),"")

  //peptide id s
  TEST_EQUAL(peptide_ids[0].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),195)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence().toString(),"LATEFSGNVPVLNAGDGSNQHPTQTLLDLFTIQETQGR")
  TEST_EQUAL(peptide_ids[0].getMetaValue("spectrum_reference"),"controllerType=0 controllerNumber=1 scan=32805")
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),182)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[1].getMetaValue("spectrum_reference"),"controllerType=0 controllerNumber=1 scan=26090")
  TEST_EQUAL(peptide_ids[2].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),191)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[2].getMetaValue("spectrum_reference"),"controllerType=0 controllerNumber=1 scan=26157")
  TEST_EQUAL(peptide_ids[3].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[3].getHits()[0].getScore(),211)
  TEST_EQUAL(peptide_ids[3].getHits()[0].getSequence().toString(),"VGAGPFPTELFDETGEFLC(Carbamidomethyl)K")
  TEST_EQUAL(peptide_ids[3].getMetaValue("spectrum_reference"),"controllerType=0 controllerNumber=1 scan=15094")

END_SECTION

START_SECTION(void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids) )

  //store and load data from various sources, starting with idxml, contents already checked above, so checking integrity of the data over repeated r/w
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  MzIdentMLFile().load(filename, protein_ids, peptide_ids);
  TEST_EQUAL(protein_ids.size(),protein_ids2.size())
  TEST_EQUAL(protein_ids[0].getHits().size(),protein_ids2[0].getHits().size())
  TEST_EQUAL(peptide_ids.size(),peptide_ids2.size())
  TEST_EQUAL(peptide_ids[0].getHits().size(),peptide_ids2[0].getHits().size())
  TEST_EQUAL(peptide_ids[1].getHits().size(),peptide_ids2[1].getHits().size())
  TEST_EQUAL(peptide_ids[2].getHits().size(),peptide_ids2[2].getHits().size())

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getSearchEngine(),protein_ids2[0].getSearchEngine())
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),protein_ids2[0].getSearchEngineVersion())
  TEST_EQUAL(protein_ids[0].getDateTime().getDate(),protein_ids2[0].getDateTime().getDate())
  TEST_EQUAL(protein_ids[0].getDateTime().getTime(),protein_ids2[0].getDateTime().getTime())
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,protein_ids2[0].getSearchParameters().db)
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,protein_ids2[0].getSearchParameters().db_version)
  TEST_EQUAL(protein_ids[0].getSearchParameters().digestion_enzyme.getName(),protein_ids2[0].getSearchParameters().digestion_enzyme.getName())
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,protein_ids2[0].getSearchParameters().charges)
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,protein_ids2[0].getSearchParameters().mass_type)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,protein_ids2[0].getSearchParameters().fragment_mass_tolerance)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,protein_ids2[0].getSearchParameters().precursor_mass_tolerance)

  TEST_EQUAL(protein_ids[0].getSearchParameters().variable_modifications.size(),protein_ids2[0].getSearchParameters().variable_modifications.size())
  for (size_t i = 0; i < protein_ids[0].getSearchParameters().variable_modifications.size(); ++i)
  {
    TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().variable_modifications[i],protein_ids2[0].getSearchParameters().variable_modifications[i])
  }
  TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().variable_modifications.back(),"Acetyl (N-term)")
  TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications.size(),protein_ids2[0].getSearchParameters().fixed_modifications.size())
  for (size_t i = 0; i < protein_ids[0].getSearchParameters().fixed_modifications.size(); ++i)
  {
    TEST_STRING_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications[i],protein_ids2[0].getSearchParameters().fixed_modifications[i])
  }

  //ProteinGroups not nupported yet, also no ProteinDetection, too few input here
//  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 0);
//  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 0);

  //protein hit 1
  TEST_EQUAL(protein_ids[0].getHits()[0].getAccession(),protein_ids2[0].getHits()[0].getAccession())
  TEST_EQUAL(protein_ids[0].getHits()[0].getSequence(),protein_ids2[0].getHits()[0].getSequence())
  //protein hit 2
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),protein_ids2[0].getHits()[1].getAccession())
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),protein_ids2[0].getHits()[1].getSequence())

  //peptide id 1
  TEST_EQUAL(peptide_ids[0].getScoreType(),peptide_ids2[0].getScoreType())
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),peptide_ids2[0].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),peptide_ids2[0].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),peptide_ids2[0].getRT())
  TEST_EQUAL(peptide_ids[0].getMetaValue("spectrum_reference"),peptide_ids2[0].getMetaValue("spectrum_reference"))
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),peptide_ids2[0].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(),peptide_ids2[0].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),peptide_ids2[0].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[0].getHits()[0].getPeptideEvidences().size(); ++i){
    //AA before/after tested by peptide evidences vector equality check - not working if the order of proteins is pertubated
    //TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[0].getHits()[0].getPeptideEvidences()[i],true)
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getStart(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getStart())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getEnd(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getEnd())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getAABefore(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getAABefore())
    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getAAAfter(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getAAAfter())
//    TEST_EQUAL(peptide_ids[0].getHits()[0].getPeptideEvidences()[i].getProteinAccession(),peptide_ids2[0].getHits()[0].getPeptideEvidences()[i].getProteinAccession())
  }
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),peptide_ids2[0].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(),peptide_ids2[0].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),peptide_ids2[0].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[0].getHits()[1].getPeptideEvidences().size(); ++i){
//    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[0].getHits()[1].getPeptideEvidences()[i],true)
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getStart(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getStart())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getEnd(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getEnd())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getAABefore(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getAABefore())
    TEST_EQUAL(peptide_ids[0].getHits()[1].getPeptideEvidences()[i].getAAAfter(),peptide_ids2[0].getHits()[1].getPeptideEvidences()[i].getAAAfter())
  }

  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),peptide_ids2[1].getScoreType())
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),peptide_ids2[1].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[1].getMZ(),peptide_ids2[1].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[1].getRT(),peptide_ids2[1].getRT())
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),peptide_ids2[1].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(),peptide_ids2[1].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),peptide_ids2[1].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[1].getHits()[0].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[1].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[1].getHits()[0].getPeptideEvidences()[i],true)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(),peptide_ids2[1].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[1].getHits()[1].getCharge(),peptide_ids2[1].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[1].getHits()[1].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[1].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[1].getHits()[1].getPeptideEvidences()[i],true)
  //peptide id 3
  TEST_EQUAL(peptide_ids[2].getScoreType(),peptide_ids2[2].getScoreType())
  TEST_EQUAL(peptide_ids[2].isHigherScoreBetter(),peptide_ids2[2].isHigherScoreBetter())
  TEST_REAL_SIMILAR(peptide_ids[2].getMZ(),peptide_ids2[2].getMZ())
  TEST_REAL_SIMILAR(peptide_ids[2].getRT(),peptide_ids2[2].getRT())
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),peptide_ids2[2].getHits()[0].getScore())
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(),peptide_ids2[2].getHits()[0].getSequence())
  TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),peptide_ids2[2].getHits()[0].getCharge())
  for (size_t i = 0; i < peptide_ids[2].getHits()[0].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[2].getHits()[0].getPeptideEvidences()[i]==peptide_ids2[2].getHits()[0].getPeptideEvidences()[i],true)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids[2].getHits()[1].getSequence(),peptide_ids2[2].getHits()[1].getSequence())
  TEST_EQUAL(peptide_ids[2].getHits()[1].getCharge(),peptide_ids2[2].getHits()[1].getCharge())
  for (size_t i = 0; i < peptide_ids[2].getHits()[1].getPeptideEvidences().size(); ++i)
    TEST_EQUAL(peptide_ids[2].getHits()[1].getPeptideEvidences()[i]==peptide_ids2[2].getHits()[1].getPeptideEvidences()[i],true)

END_SECTION

START_SECTION(([EXTRA] multiple runs))
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  MzIdentMLFile().load(filename, protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),protein_ids2.size())

  TEST_EQUAL(protein_ids[0].getHits().size(),protein_ids2[0].getHits().size())
  TEST_EQUAL(protein_ids[1].getHits().size(),protein_ids2[1].getHits().size())
  TEST_EQUAL(protein_ids[2].getHits().size(),protein_ids2[2].getHits().size())

  TEST_EQUAL(protein_ids[0].getSearchParameters().precursor_mass_tolerance_ppm, true)

END_SECTION

START_SECTION(([EXTRA] psm ranking))
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    size_t r = 0;
    for (size_t j = 0; j < peptide_ids[i].getHits().size(); ++j)
    {
      TEST_EQUAL(peptide_ids[i].getHits()[j].getRank()>=r,true)
      r = peptide_ids[i].getHits()[j].getRank();
    }
  }
END_SECTION

START_SECTION(([EXTRA] thresholds))
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),1)
  TEST_EQUAL(protein_ids[0].getSignificanceThreshold(),0.5)

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    if (peptide_ids[i].getMetaValue("spectrum_reference") == "17")
    {
      TEST_EQUAL(peptide_ids[i].getHits().size(),2)
      for (size_t j = 0; j < peptide_ids[i].getHits().size(); ++j)
      {
        TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),false)
      }
      PeptideHit x = peptide_ids[i].getHits().back();
      x.removeMetaValue("pass_threshold");
      x.setSequence(AASequence::fromString("TESTER"));
      x.setScore(0.4);
      peptide_ids[i].insertHit(x);
    }
  }

  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  protein_ids.clear();
  peptide_ids.clear();
  MzIdentMLFile().load(filename, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    if (peptide_ids[i].getMetaValue("spectrum_reference") == "17")
    {
      TEST_EQUAL(peptide_ids[i].getHits().size(),3)
      for (size_t j = 0; j < peptide_ids[i].getHits().size(); ++j)
      {
        if (peptide_ids[i].getHits()[j].getScore() > protein_ids[0].getSignificanceThreshold())
        {
          TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),false)
        }
        else
          TEST_EQUAL(peptide_ids[i].getHits()[j].getMetaValue("pass_threshold"),true)
      }
    }
  }


END_SECTION
START_SECTION(([EXTRA] regression test for file loading on example files))
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
//  input_path = OPENMS_GET_TEST_DATA_PATH("Mascot_MSMS_example.mzid");
//  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
  input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);
  input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);

END_SECTION


//START_SECTION(([EXTRA] compability issues))
//  MzIdentMLFile mzidfile;
//  vector<ProteinIdentification> protein_ids;
//  vector<PeptideIdentification> peptide_ids;
//  mzidfile.load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_no_proteinhits.mzid"), protein_ids, peptide_ids);

//  TEST_EQUAL(protein_ids.size(), 1)
//  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
//  TEST_EQUAL(peptide_ids.size(), 10)
//  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

//  String filename;
//  NEW_TMP_FILE(filename)
//  mzidfile.store(filename , protein_ids, peptide_ids);

//  vector<ProteinIdentification> protein_ids2;
//  vector<PeptideIdentification> peptide_ids2;
//  mzidfile.load(filename, protein_ids2, peptide_ids2);

//  TEST_EQUAL(protein_ids == protein_ids2, true)
//  TEST_EQUAL(peptide_ids == peptide_ids2, true)

//  Misplaced Elements ignored in ParamGroup

//  Converting unknown score type to search engine specific score CV. #should not occurr  scoretype is whatever

//  PSM without peptide evidences reigstered in the given search database found. This will cause an invalid MzIdentML file (which OpenMS still can consume). #might occurr when reading idxml. no protein reference accession

//  No RT #might occurr when reading idxml. no rt to peptidehit
//  No MZ #might occurr when reading idxml. no mz to peptidehit

//  PeptideEvidence without reference to the positional in originating sequence found. #will always occurr when reading idxml  no start end positional arguments
//END_SECTION

START_SECTION(([EXTRA] XLMS data labeled cross-linker))
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;

  String input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_labelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids[1].getHits()[1].getMetaValue("xl_pos"), 0)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getMetaValue("xl_term_spec"), "N_TERM")
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence().toString(), "KELLK")

  // Reading and writing
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // parameters from written and reloaded file
  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"), "138.0680796")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass_isoshift"), "12.075321")
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // cross-linking search

  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 4)
  TEST_EQUAL(peptide_ids2[1].getRT(), peptide_ids2[2].getRT())
  TEST_REAL_SIMILAR(peptide_ids2[1].getRT(), 2132.4757)
  TEST_REAL_SIMILAR(peptide_ids2[1].getMZ(), 721.0845)
  TEST_EQUAL(peptide_ids2[1].getMetaValue("spectrum_reference"), peptide_ids2[2].getMetaValue("spectrum_reference"))
  TEST_EQUAL(peptide_ids2[1].getMetaValue("spectrum_reference"), "controllerType=0 controllerNumber=1 scan=3647,controllerType=0 controllerNumber=1 scan=3539")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[3].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 2)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_chain"), "MS:1002509") // XL donor
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getMetaValue("xl_chain"), "MS:1002510") // XL acceptor
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_type"), "cross-link")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("spec_heavy_RT"), 2089.55329999998)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("spec_heavy_MZ"), 725.108947753906)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getScore(), peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getScore(), peptide_ids2[2].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "LM(Oxidation)VEMEKKLEK")
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getSequence().toString(), "KELLK")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_pos"), 6)
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getMetaValue("xl_pos"), 0)
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getMetaValue("xl_term_spec"), "N_TERM")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_mass"), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_mod"), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b8]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[20].annotation, "[alpha|xi$b9]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[25].charge, 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[25].annotation, "[alpha|xi$y8]")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_type"), "loop-link")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_pos"), 7)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_pos2"), 14)
  TEST_EQUAL(peptide_ids2[3].getHits()[0].getMetaValue("xl_type"), "mono-link")

END_SECTION

START_SECTION(([EXTRA] XLMS data unlabeled cross-linker))
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;

  String input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_unlabelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);

  // Reading and writing
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"), "138.0680796")
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // cross-linking search

  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 4)
  TEST_EQUAL(peptide_ids2[0].getRT(), peptide_ids[1].getRT())
  TEST_REAL_SIMILAR(peptide_ids2[0].getRT(), 2132.4757)
  TEST_REAL_SIMILAR(peptide_ids2[0].getMZ(), 721.0845)
  TEST_EQUAL(peptide_ids2[0].getMetaValue("spectrum_reference"), peptide_ids2[1].getMetaValue("spectrum_reference"))
  TEST_EQUAL(peptide_ids2[0].getMetaValue("spectrum_reference"), "controllerType=0 controllerNumber=1 scan=3647")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 2)
  TEST_EQUAL(peptide_ids2[3].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 2)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_chain"), "MS:1002509") // XL donor
  TEST_EQUAL(peptide_ids2[0].getHits()[1].getMetaValue("xl_chain"), "MS:1002510") // XL acceptor
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_type"), "cross-link")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_pos"), 0)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue("xl_term_spec"), "N_TERM")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getScore(), peptide_ids2[0].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getScore(), peptide_ids2[1].getHits()[1].getScore())
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getSequence().toString(), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_pos"), 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[1].getMetaValue("xl_pos"), 4)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_mass"), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue("xl_mod"), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$b3]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|ci$y5]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[20].annotation, "[beta|xi$y6]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[26].charge, 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[26].annotation, "[alpha|xi$y8]")
  TEST_EQUAL(peptide_ids2[3].getHits()[0].getSequence().toString(), "VLVKVHPEGKYVVDISPDIDIK")
  TEST_EQUAL(peptide_ids2[3].getHits()[0].getMetaValue("xl_type"), "loop-link")
  TEST_EQUAL(peptide_ids2[3].getHits()[0].getMetaValue("xl_pos"), 3)
  TEST_EQUAL(peptide_ids2[3].getHits()[0].getMetaValue("xl_pos2"), 9)
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue("xl_type"), "mono-link")



END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
