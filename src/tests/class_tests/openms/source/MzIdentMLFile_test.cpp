// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/CONCEPT/Constants.h>


using namespace OpenMS;
using namespace std;

START_TEST(MzIdentMLFile, "$Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MzIdentMLFile* ptr = nullptr;
MzIdentMLFile* nullPointer = nullptr;

START_SECTION((MzIdentMLFile()))
{
  ptr = new MzIdentMLFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MzIdentMLFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids) )
{
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  std::vector<String> fm{"Carbamidomethyl (C)", "Xlink:DTSSP[88] (Protein N-term)"};
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
  ABORT_IF(protein_ids[0].getSearchParameters().fixed_modifications.size() != fm.size())
  for (size_t i = 0; i < fm.size(); ++i)
  {
    TEST_EQUAL(protein_ids[0].getSearchParameters().fixed_modifications[i], fm[i]);
  }
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
  TEST_EQUAL(peptide_ids[0].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=32805")
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),182)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[1].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=26090")
  TEST_EQUAL(peptide_ids[2].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),191)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence().toString(),"FLAETDQGPVPVEITAVEDDHVVVDGNHMLAGQNLK")
  TEST_EQUAL(peptide_ids[2].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=26157")
  TEST_EQUAL(peptide_ids[3].getScoreType(),"MS-GF:RawScore")
  TEST_REAL_SIMILAR(peptide_ids[3].getHits()[0].getScore(),211)
  TEST_EQUAL(peptide_ids[3].getHits()[0].getSequence().toString(),"VGAGPFPTELFDETGEFLC(Carbamidomethyl)K")
  TEST_EQUAL(peptide_ids[3].getSpectrumReference(),"controllerType=0 controllerNumber=1 scan=15094")
}
END_SECTION

START_SECTION(void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids) )
{
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
  TEST_EQUAL(peptide_ids[0].getSpectrumReference(),peptide_ids2[0].getSpectrumReference())
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
}
END_SECTION

START_SECTION(([EXTRA] multiple runs))
{
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
}
END_SECTION

START_SECTION(([EXTRA] psm ranking))
{
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
}
END_SECTION

START_SECTION(([EXTRA] thresholds))
{
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  String input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),1)
  TEST_EQUAL(protein_ids[0].getSignificanceThreshold(),0.5)

  TEST_EQUAL(peptide_ids.size(),3)
  for (size_t i = 0; i < peptide_ids.size(); ++i)
  {
    if (peptide_ids[i].getSpectrumReference() == "17")
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
    if (peptide_ids[i].getSpectrumReference() == "17")
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

}
END_SECTION

START_SECTION(([EXTRA] regression test for file loading on example files))
{
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
}
END_SECTION


START_SECTION(([EXTRA] compability issues))
{
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

//  TEST_TRUE(protein_ids == protein_ids2)
//  TEST_TRUE(peptide_ids == peptide_ids2)

//  Misplaced Elements ignored in ParamGroup

//  Converting unknown score type to search engine specific score CV. #should not occur  scoretype is whatever

//  PSM without peptide evidences registered in the given search database found. This will cause an invalid MzIdentML file (which OpenMS still can consume). #might occur when reading idxml. no protein reference accession

//  No RT #might occurr when reading idxml. no rt to peptidehit
//  No MZ #might occurr when reading idxml. no mz to peptidehit

//  PeptideEvidence without reference to the positional in originating sequence found. #will always occur when reading idxml  no start end positional arguments
}
END_SECTION

START_SECTION(([EXTRA] XLMS data labeled cross-linker))
{
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;

  String input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_labelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  //
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 3)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")

  // Reading and writing
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // parameters from written and reloaded file
  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_REAL_SIMILAR(String(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass")).toDouble(), 138.0680796)
  TEST_REAL_SIMILAR(String(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass_isoshift")).toDouble(), 12.075321)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("extra_features"), "precursor_mz_error_ppm,\
OpenPepXL:score,isotope_error,OpenPepXL:xquest_score,OpenPepXL:xcorr xlink,\
OpenPepXL:xcorr common,OpenPepXL:match-odds,OpenPepXL:intsum,OpenPepXL:wTIC,OpenPepXL:TIC,OpenPepXL:prescore,OpenPepXL:log_occupancy,\
OpenPepXL:log_occupancy_alpha,OpenPepXL:log_occupancy_beta,matched_xlink_alpha,matched_xlink_beta,matched_linear_alpha,\
matched_linear_beta,ppm_error_abs_sum_linear_alpha,ppm_error_abs_sum_linear_beta,ppm_error_abs_sum_xlinks_alpha,\
ppm_error_abs_sum_xlinks_beta,ppm_error_abs_sum_linear,ppm_error_abs_sum_xlinks,ppm_error_abs_sum_alpha,ppm_error_abs_sum_beta,\
ppm_error_abs_sum,precursor_total_intensity,precursor_target_intensity,precursor_signal_proportion,precursor_target_peak_count,\
precursor_residual_peak_count")
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // crosslinking search


  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 10)
  TEST_EQUAL(peptide_ids2[1].getRT(), peptide_ids2[2].getRT())
  TEST_REAL_SIMILAR(peptide_ids2[1].getRT(), 2132.4757)
  TEST_REAL_SIMILAR(peptide_ids2[1].getMZ(), 721.0845)
  TEST_EQUAL(peptide_ids2[1].getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE), "spectrum=131,spectrum=113")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[3].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT), 2125.5966796875)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ), 725.109252929687841)
  TEST_REAL_SIMILAR(peptide_ids2[1].getHits()[0].getScore(), -0.190406834856118)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 4)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_REAL_SIMILAR(String(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS)).toDouble(), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[beta|ci$y2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b4]")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 5)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
}
END_SECTION

START_SECTION(([EXTRA] XLMS data unlabeled cross-linker))
{
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;

  String input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_unlabelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);

  // Reading and writing
  String filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_EQUAL(String(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass")).toDouble(), 138.0680796)
  TEST_EQUAL(protein_ids[0].getMetaValue("SpectrumIdentificationProtocol"), "MS:1002494") // crosslinking search

  // PeptideIdentification (Indices may change, without making the reading/writing invalid, if e.g. more is added to the test file)
  TEST_EQUAL(peptide_ids2.size(), 3)
  TEST_REAL_SIMILAR(peptide_ids2[0].getRT(), 2175.3003)
  TEST_REAL_SIMILAR(peptide_ids2[0].getMZ(), 787.740356445313)
  TEST_EQUAL(peptide_ids2[0].getMetaValue(Constants::UserParam::SPECTRUM_REFERENCE), "controllerType=0 controllerNumber=1 scan=2395")

  // PeptideHit
  TEST_EQUAL(peptide_ids2[0].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[1].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[2].getHits().size(), 1)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")

  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 5)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")

  TEST_EQUAL(peptide_ids2[1].getHits()[0].getSequence().toString(), "KNVPIEFPVIDR")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "LGCKALHVLFER")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 0)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 3)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations().size(), 5)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$y5]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$y7]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[2].annotation, "[beta|ci$y7]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[3].annotation, "[alpha|ci$y8]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[2].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[4].charge, 2)

  TEST_EQUAL(peptide_ids2[2].getHits()[0].getSequence().toString(), "VEPSWLGPLFPDK(Xlink:DSS[156])TSNLR")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "-")
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 12)
  TEST_EQUAL(peptide_ids2[2].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
