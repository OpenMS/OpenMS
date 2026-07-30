// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/FORMAT/TextFile.h>
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

START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) )
{
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::vector<std::string> fm{"Carbamidomethyl (C)", "Xlink:DTSSP[88] (Protein N-term)"};
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(protein_ids[0].getHits().size(),2)
  TEST_EQUAL(protein_ids[1].getHits().size(),1)
  TEST_EQUAL(peptide_ids.size(),5)

  // loading again into the same containers must replace, not accumulate
  // (regression: the DOM handler only appends; load() now clears first, matching IdXMLFile::load)
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);
  TEST_EQUAL(protein_ids.size(),2)
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

START_SECTION(([EXTRA] read mzIdentML Modification without the optional location attribute (issue #5443)))
{
  // The 'location' attribute of <Modification> is optional in mzIdentML; some search engines omit it for
  // terminal modifications. Reading such a file used to crash with an IndexOverflow exception. The reader
  // must now infer a terminal position (or skip the modification) instead of applying a negative index.
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_missing_mod_location.mzid"), protein_ids, peptide_ids);

  ABORT_IF(peptide_ids.size() != 5)
  for (Size i = 0; i < peptide_ids.size(); ++i)
  {
    ABORT_IF(peptide_ids[i].getHits().empty())
  }

  // 1) N-terminal Acetyl without 'location' -> inferred as N-terminal
  const AASequence& acetyl_seq = peptide_ids[0].getHits()[0].getSequence();
  TEST_EQUAL(acetyl_seq.toUnmodifiedString(), "PEPTIDEK")
  TEST_TRUE(acetyl_seq.hasNTerminalModification())
  TEST_FALSE(acetyl_seq.hasCTerminalModification())
  TEST_EQUAL(acetyl_seq.getNTerminalModificationName(), "Acetyl")

  // 2) Oxidation without 'location' -> position not uniquely terminal -> skipped (but no crash and no misplacement)
  const AASequence& oxidation_seq = peptide_ids[1].getHits()[0].getSequence();
  TEST_EQUAL(oxidation_seq.toUnmodifiedString(), "PEPTIDEM")
  TEST_FALSE(oxidation_seq.isModified())
  TEST_EQUAL(oxidation_seq.toString(), "PEPTIDEM")

  // 3) Amidated without 'location' -> exclusively C-terminal -> inferred as C-terminal
  const AASequence& amidated_seq = peptide_ids[2].getHits()[0].getSequence();
  TEST_EQUAL(amidated_seq.toUnmodifiedString(), "PEPTIDER")
  TEST_TRUE(amidated_seq.hasCTerminalModification())
  TEST_FALSE(amidated_seq.hasNTerminalModification())
  TEST_EQUAL(amidated_seq.getCTerminalModificationName(), "Amidated")

  // 4) Control: a regular internal modification WITH a valid 'location' is still parsed correctly
  const AASequence& carbamidomethyl_seq = peptide_ids[3].getHits()[0].getSequence();
  TEST_EQUAL(carbamidomethyl_seq.toString(), "PEPTIDEC(Carbamidomethyl)K")

  // 5) Acetyl without 'location' but with residues="K" -> residue-specific (internal), position unknown,
  //    so it is skipped rather than forced onto the N-terminus even though Acetyl supports the N-terminus
  const AASequence& acetyl_on_k_seq = peptide_ids[4].getHits()[0].getSequence();
  TEST_EQUAL(acetyl_on_k_seq.toUnmodifiedString(), "PEPTIKDE")
  TEST_FALSE(acetyl_on_k_seq.isModified())
}
END_SECTION

START_SECTION(void store(std::string filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids) )
{
  //store and load data from various sources, starting with idxml, contents already checked above, so checking integrity of the data over repeated r/w
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);

  // Exercise regular (non-XL) fragment annotation persistence. The two known
  // errors and one unknown error share an IonType, requiring separate parallel
  // error-array groups in mzIdentML without inventing a sentinel value.
  vector<PeptideHit::PeakAnnotation> fragment_annotations(3);
  fragment_annotations[0].annotation = "[y3]";
  fragment_annotations[0].charge = 1;
  fragment_annotations[0].mz = 500.125;
  fragment_annotations[0].intensity = 1234.0;
  fragment_annotations[0].theoretical_mz = 500.1;
  fragment_annotations[1].annotation = "[y4]";
  fragment_annotations[1].charge = 1;
  fragment_annotations[1].mz = 600.2;
  fragment_annotations[1].intensity = 2345.0;
  fragment_annotations[2].annotation = "[y5]";
  fragment_annotations[2].charge = 1;
  fragment_annotations[2].mz = 700.3;
  fragment_annotations[2].intensity = 3456.0;
  fragment_annotations[2].theoretical_mz = 700.3;
  peptide_ids2[0].getHits()[0].setPeakAnnotations(fragment_annotations);

  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids2, peptide_ids2);

  // Parallel FragmentArrays must stay aligned with IonType::index. Corrupt the
  // generated two-value error array and verify that loading fails explicitly.
  TextFile malformed_file(filename);
  bool error_array_corrupted = false;
  for (auto& line : malformed_file)
  {
    if (line.find("<FragmentArray measure_ref=\"Measure_error\"") == std::string::npos)
    {
      continue;
    }
    const Size values_attribute = line.find("values=\"");
    if (values_attribute == std::string::npos)
    {
      continue;
    }
    const Size values_begin = values_attribute + 8;
    const Size separator = line.find(' ', values_begin);
    if (separator != std::string::npos)
    {
      const Size values_end = line.find('"', separator);
      if (values_end != std::string::npos)
      {
        line.erase(separator, values_end - separator);
        error_array_corrupted = true;
        break;
      }
    }
  }
  TEST_TRUE(error_array_corrupted)
  std::string malformed_filename;
  NEW_TMP_FILE(malformed_filename)
  malformed_file.store(malformed_filename);
  vector<ProteinIdentification> malformed_protein_ids;
  PeptideIdentificationList malformed_peptide_ids;
  TEST_EXCEPTION(Exception::ParseError, MzIdentMLFile().load(malformed_filename, malformed_protein_ids, malformed_peptide_ids))

  // Measure IDs are arbitrary. Also exercise ppm input and the legacy
  // capitalized Measure_ref spelling used by existing mzIdentML files.
  TextFile ppm_file(filename);
  bool ppm_measure_changed = false;
  bool ppm_values_changed = false;
  for (auto& line : ppm_file)
  {
    if (line.find("<Measure id=\"Measure_error\"") != std::string::npos)
    {
      StringUtils::substitute(line, "Measure_error", "custom_error_ppm");
      ppm_measure_changed = true;
    }
    else if (line.find("accession=\"MS:1001227\"") != std::string::npos)
    {
      StringUtils::substitute(line, "unitAccession=\"MS:1000040\"", "unitAccession=\"UO:0000169\"");
      StringUtils::substitute(line, "unitCvRef=\"PSI-MS\"", "unitCvRef=\"UO\"");
      StringUtils::substitute(line, "unitName=\"m/z\"", "unitName=\"parts per million\"");
    }
    else if (line.find("<FragmentArray measure_ref=\"Measure_error\"") != std::string::npos)
    {
      StringUtils::substitute(line, "measure_ref=\"Measure_error\"", "Measure_ref=\"custom_error_ppm\"");
      const Size values_attribute = line.find("values=\"");
      ABORT_IF(values_attribute == std::string::npos)
      const Size values_begin = values_attribute + 8;
      const Size values_end = line.find('"', values_begin);
      ABORT_IF(values_end == std::string::npos)
      const double y3_error_ppm = (500.125 - 500.1) / 500.1 * 1e6;
      line.replace(values_begin, values_end - values_begin, StringUtils::toStr(y3_error_ppm) + " 0");
      ppm_values_changed = true;
    }
  }
  TEST_TRUE(ppm_measure_changed)
  TEST_TRUE(ppm_values_changed)
  std::string ppm_filename;
  NEW_TMP_FILE(ppm_filename)
  ppm_file.store(ppm_filename);
  vector<ProteinIdentification> ppm_protein_ids;
  PeptideIdentificationList ppm_peptide_ids;
  MzIdentMLFile().load(ppm_filename, ppm_protein_ids, ppm_peptide_ids);
  bool found_ppm_y3 = false;
  for (const auto& annotation : ppm_peptide_ids[0].getHits()[0].getPeakAnnotations())
  {
    if (annotation.annotation == "[y3]")
    {
      found_ppm_y3 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, 500.1)
    }
  }
  TEST_TRUE(found_ppm_y3)

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
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations().size(), 3)
  bool found_y3 = false;
  bool found_y4 = false;
  bool found_y5 = false;
  for (const auto& annotation : peptide_ids[0].getHits()[0].getPeakAnnotations())
  {
    if (annotation.annotation == "[y3]")
    {
      found_y3 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, 500.1)
    }
    else if (annotation.annotation == "[y4]")
    {
      found_y4 = true;
      TEST_FALSE(annotation.theoretical_mz.has_value())
    }
    else if (annotation.annotation == "[y5]")
    {
      found_y5 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, 700.3)
    }
  }
  TEST_TRUE(found_y3)
  TEST_TRUE(found_y4)
  TEST_TRUE(found_y5)
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
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentML_3runs.mzid");
  MzIdentMLFile().load(input_path, protein_ids2, peptide_ids2);
  std::string filename;
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

START_SECTION(([EXTRA] thresholds))
{
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
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

  std::string filename;
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
  PeptideIdentificationList peptide_ids;
  std::string input_path = OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_whole.mzid");
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
//  PeptideIdentificationList peptide_ids;
//  mzidfile.load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_no_proteinhits.mzid"), protein_ids, peptide_ids);

//  TEST_EQUAL(protein_ids.size(), 1)
//  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
//  TEST_EQUAL(peptide_ids.size(), 10)
//  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

//  std::string filename;
//  NEW_TMP_FILE(filename)
//  mzidfile.store(filename , protein_ids, peptide_ids);

//  vector<ProteinIdentification> protein_ids2;
//  PeptideIdentificationList peptide_ids2;
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
  PeptideIdentificationList peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;

  std::string input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_labelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  //
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 3)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "SAVIKTSTR")
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence().toString(), "FIVKASSGPR")

  auto& xl_fragment_annotations = peptide_ids[1].getHits()[0].getPeakAnnotations();
  ABORT_IF(xl_fragment_annotations.size() < 9)
  bool set_alpha_ci_b2 = false;
  bool set_alpha_xi_b4 = false;
  bool set_alpha_xi_y8 = false;
  for (auto& annotation : xl_fragment_annotations)
  {
    if (annotation.annotation == "[alpha|ci$b2]")
    {
      annotation.theoretical_mz = annotation.mz - 0.0125;
      set_alpha_ci_b2 = true;
    }
    else if (annotation.annotation == "[alpha|xi$b4]")
    {
      annotation.theoretical_mz = annotation.mz + 0.025;
      set_alpha_xi_b4 = true;
    }
    else if (annotation.annotation == "[alpha|xi$y8]" && annotation.charge == 3)
    {
      annotation.theoretical_mz = 0.0;
      set_alpha_xi_y8 = true;
    }
  }
  TEST_TRUE(set_alpha_ci_b2)
  TEST_TRUE(set_alpha_xi_b4)
  TEST_TRUE(set_alpha_xi_y8)

  // Reading and writing
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // parameters from written and reloaded file
  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"))), 138.0680796)
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass_isoshift"))), 12.075321)
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
  TEST_REAL_SIMILAR(StringUtils::toDouble(StringUtils::toStr(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS))), 138.0680796)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD), "DSS")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].annotation, "[beta|ci$y2]")
  TEST_FALSE(peptide_ids2[1].getHits()[0].getPeakAnnotations()[1].theoretical_mz.has_value())
  TEST_EQUAL(peptide_ids2[1].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b4]")
  bool found_alpha_ci_b2 = false;
  bool found_alpha_xi_b4 = false;
  bool found_alpha_xi_y8 = false;
  for (const auto& annotation : peptide_ids2[1].getHits()[0].getPeakAnnotations())
  {
    if (annotation.annotation == "[alpha|ci$b2]")
    {
      found_alpha_ci_b2 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, annotation.mz - 0.0125)
    }
    else if (annotation.annotation == "[alpha|xi$b4]")
    {
      found_alpha_xi_b4 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, annotation.mz + 0.025)
    }
    else if (annotation.annotation == "[alpha|xi$y8]" && annotation.charge == 3)
    {
      found_alpha_xi_y8 = true;
      TEST_TRUE(annotation.theoretical_mz.has_value())
      TEST_REAL_SIMILAR(*annotation.theoretical_mz, 0.0)
    }
  }
  TEST_TRUE(found_alpha_ci_b2)
  TEST_TRUE(found_alpha_xi_b4)
  TEST_TRUE(found_alpha_xi_y8)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "mono-link")
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 5)
  TEST_EQUAL(peptide_ids2[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), "-")
}
END_SECTION

START_SECTION(([EXTRA] XLMS data unlabeled cross-linker))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;

  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_unlabelled.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);

  // Reading and writing
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);

  // ProteinIdentification
  TEST_EQUAL(protein_ids2[0].getSearchParameters().fragment_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().precursor_mass_tolerance_ppm, true)
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue1"), "[K, N-term]")
  TEST_EQUAL(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:residue2"), "[K, N-term]")
  TEST_EQUAL(StringUtils::toDouble(StringUtils::toStr(protein_ids2[0].getSearchParameters().getMetaValue("cross_link:mass"))), 138.0680796)
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

START_SECTION(([EXTRA] mzIdentML 1.3 crosslinking scores_and_thresholds))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_crosslinking.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // every parsed identification carries a spectrum reference; collect hits across all ids
  // (a SpectrumIdentificationResult may legitimately yield an id without hits, so guard the hit access)
  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    TEST_EQUAL(pid.getSpectrumReference().empty(), false)
    total_hits += pid.getHits().size();
    if (!pid.getHits().empty())
    {
      TEST_EQUAL(pid.getHits()[0].getSequence().empty(), false)
    }
  }
  TEST_EQUAL(total_hits > 0, true)

  // store/load roundtrip: the writer must default to mzIdentML 1.3.0 and the reader must read it back without loss
  std::string filename;
  NEW_TMP_FILE(filename)
  MzIdentMLFile().store(filename, protein_ids, peptide_ids);

  // the default writer path emits version 1.3.0
  TextFile written(filename);
  bool has_version_130 = false;
  for (TextFile::ConstIterator it = written.begin(); it != written.end(); ++it)
  {
    if (StringUtils::hasSubstring(*it, "version=\"1.3.0\"")) { has_version_130 = true; break; }
  }
  TEST_EQUAL(has_version_130, true)

  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
  MzIdentMLFile().load(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(!peptide_ids2.empty(), true)
  TEST_EQUAL(!protein_ids2.empty(), true)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 noncovalent association))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_noncov_assoc.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // parsed hits must have valid sequences (modification/cvParam fallbacks resolved correctly)
  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    total_hits += pid.getHits().size();
    if (!pid.getHits().empty())
    {
      TEST_EQUAL(pid.getHits()[0].getSequence().empty(), false)
    }
  }
  TEST_EQUAL(total_hits > 0, true)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 EDC crosslinking))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_xlink_edc.mzid");
  // EDC files mix crosslinked and standalone (non-XL) peptides; both must parse without throwing
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  Size total_hits = 0;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    total_hits += pid.getHits().size();
  }
  TEST_EQUAL(total_hits > 0, true)
}
END_SECTION

START_SECTION(([EXTRA] mzIdentML 1.3 multiple spectra per identification))
{
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string input_file = OPENMS_GET_TEST_DATA_PATH("MzIdentML_v1.3_multi_spectra.mzid");
  MzIdentMLFile().load(input_file, protein_ids, peptide_ids);
  TEST_EQUAL(!protein_ids.empty(), true)
  TEST_EQUAL(!peptide_ids.empty(), true)

  // multiple SpectrumIdentificationLists yield several distinct spectrum references
  set<std::string> refs;
  for (const PeptideIdentification& pid : peptide_ids)
  {
    refs.insert(pid.getSpectrumReference());
  }
  TEST_EQUAL(refs.size() > 1, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
