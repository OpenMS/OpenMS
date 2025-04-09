// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

///////////////////////////

START_TEST(IdXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((IdXMLFile()))
  IdXMLFile* ptr = nullptr;
  IdXMLFile* nullPointer = nullptr;
  ptr = new IdXMLFile();
  TEST_NOT_EQUAL(ptr,nullPointer)
  delete ptr;
END_SECTION

START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids) )
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),3)
END_SECTION


START_SECTION(void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, String& document_id) )
  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids, document_id);

  TEST_STRING_EQUAL(document_id,"LSID1234")
  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),3)

  /////////////// protein id 1 //////////////////
  TEST_EQUAL(protein_ids[0].getScoreType(),"MOWSE")
  TEST_EQUAL(protein_ids[0].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[0].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[0].getSearchEngineVersion(),"2.1.0")
  TEST_EQUAL(protein_ids[0].getDateTime().getDate(),"2006-01-12")
  TEST_EQUAL(protein_ids[0].getDateTime().getTime(),"12:13:14")
  TEST_EQUAL(protein_ids[0].getIdentifier().hasPrefix("Mascot_2006-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"1, 2")
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::AVERAGE)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,1.0)
  TEST_EQUAL((String)(protein_ids[0].getMetaValue("name")),"ProteinIdentification")

  TEST_EQUAL(protein_ids[0].getProteinGroups().size(), 1);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].probability, 0.88);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions.size(), 2);
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions[0], "PROT1");
  TEST_EQUAL(protein_ids[0].getProteinGroups()[0].accessions[1], "PROT2");

  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins().size(), 1);
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions.size(),
             2);
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions[0],
             "PROT1");
  TEST_EQUAL(protein_ids[0].getIndistinguishableProteins()[0].accessions[1],
             "PROT2");

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
  TEST_EQUAL(peptide_ids[0].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),false)
  TEST_EQUAL(peptide_ids[0].getIdentifier().hasPrefix("Mascot_2006-01-12T12:13:14"), true)
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),675.9)
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),1234.5)
  TEST_EQUAL((peptide_ids[0].getSpectrumReference()),"17")
  TEST_EQUAL((String)(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
  TEST_EQUAL(peptide_ids[0].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),0.9)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(), AASequence::fromString("PEPTIDER"))
  TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),1)
  vector<PeptideEvidence> pes0 = peptide_ids[0].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes0.size(),2)
  TEST_EQUAL(pes0[0].getProteinAccession(),"PROT1")
  TEST_EQUAL(pes0[1].getProteinAccession(),"PROT2")
  TEST_EQUAL(pes0[0].getAABefore(),'A')
  TEST_EQUAL(pes0[0].getAAAfter(),'B')
  TEST_EQUAL((String)(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
  vector<PeptideEvidence> pes1 = peptide_ids[0].getHits()[1].getPeptideEvidences();
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERR"))
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
  TEST_EQUAL(pes1.size(),0)
  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(peptide_ids[1].getIdentifier().hasPrefix("Mascot_2006-01-12T12:13:14"), true)
  TEST_EQUAL(peptide_ids[1].getHits().size(),2)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),44.4)
  TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),2)
  vector<PeptideEvidence> pes2 = peptide_ids[1].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes2.size(),0)
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),33.3)
  TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERRRR"))
  TEST_EQUAL(peptide_ids[1].getHits()[1].getCharge(),2)
  vector<PeptideEvidence> pes3 = peptide_ids[1].getHits()[1].getPeptideEvidences();
  TEST_EQUAL(pes3.size(),0)

  /////////////// protein id 2 //////////////////
  TEST_EQUAL(protein_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(protein_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(protein_ids[1].getSearchEngine(),"Mascot")
  TEST_EQUAL(protein_ids[1].getSearchEngineVersion(),"2.1.1")
  TEST_EQUAL(protein_ids[1].getDateTime().getDate(),"2007-01-12")
  TEST_EQUAL(protein_ids[1].getDateTime().getTime(),"12:13:14")
  TEST_EQUAL(protein_ids[1].getIdentifier().hasPrefix("Mascot_2007-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[1].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[1].getSearchParameters().db_version,"1.1")
  TEST_EQUAL(protein_ids[1].getSearchParameters().charges,"1, 2, 3")
  TEST_EQUAL(protein_ids[1].getSearchParameters().mass_type,ProteinIdentification::MONOISOTOPIC)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().fragment_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[1].getSearchParameters().precursor_mass_tolerance,1.0)
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
  TEST_EQUAL(peptide_ids[2].getIdentifier().hasPrefix("Mascot_2007-01-12T12:13:14"), true)
  TEST_EQUAL(peptide_ids[2].getHits().size(),1)
  //peptide hit 1
  TEST_REAL_SIMILAR(peptide_ids[2].getHits()[0].getScore(),1.4)
  TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(), AASequence::fromString("PEPTIDERRRRR"))
  TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),1)
  vector<PeptideEvidence> pes4 = peptide_ids[2].getHits()[0].getPeptideEvidences();
  TEST_EQUAL(pes4.size(),1)
  TEST_EQUAL(pes4[0].getProteinAccession(),"PROT3")
  TEST_EQUAL(pes4[0].getAABefore(), PeptideEvidence::UNKNOWN_AA)
  TEST_EQUAL(pes4[0].getAAAfter(), PeptideEvidence::UNKNOWN_AA)
END_SECTION

START_SECTION(void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, const String& document_id="") )

  // load, store, and reload data
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String document_id, document_id2;
  String target_file = OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML");
  IdXMLFile().load(target_file, protein_ids2, peptide_ids2, document_id2);

  String actual_file;
  NEW_TMP_FILE(actual_file)
  IdXMLFile().store(actual_file, protein_ids2, peptide_ids2, document_id2);

  FuzzyStringComparator fuzzy;
  fuzzy.setWhitelist(ListUtils::create<String>("<?xml-stylesheet"));
  fuzzy.setAcceptableAbsolute(0.0001);
  bool result = fuzzy.compareFiles(actual_file, target_file);
  TEST_EQUAL(result, true);
END_SECTION


START_SECTION([EXTRA] static bool isValid(const String& filename))
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
  String filename;
  IdXMLFile f;

  //test if empty file is valid
  NEW_TMP_FILE(filename)
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isValid(filename, std::cerr),true);

  //test if full file is valid
  NEW_TMP_FILE(filename);
  String document_id;
  f.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids2, peptide_ids2, document_id);
  protein_ids2[0].setMetaValue("stringvalue",String("bla"));
  protein_ids2[0].setMetaValue("intvalue",4711);
  protein_ids2[0].setMetaValue("floatvalue",5.3);
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isValid(filename, std::cerr),true);

  //check if meta information can be loaded
  f.load(filename, protein_ids2, peptide_ids2, document_id);
END_SECTION

START_SECTION(([EXTRA] No protein identification bug))
  IdXMLFile id_xmlfile;
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;
  id_xmlfile.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_no_proteinhits.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(), 1)
  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
  TEST_EQUAL(peptide_ids.size(), 10)
  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

  String filename;
  NEW_TMP_FILE(filename)
  id_xmlfile.store(filename , protein_ids, peptide_ids);

  vector<ProteinIdentification> protein_ids2;
  vector<PeptideIdentification> peptide_ids2;
  id_xmlfile.load(filename, protein_ids2, peptide_ids2);

  // identifiers contain a random number when loaded to avoid ambiguities when merging ProtIDs; make them equal for our purposes here
  protein_ids2[0].setIdentifier(protein_ids[0].getIdentifier());
  for (auto& pep : peptide_ids2)
  {
    pep.setIdentifier(peptide_ids[0].getIdentifier());
  }

  TEST_TRUE(protein_ids == protein_ids2)
  TEST_TRUE(peptide_ids == peptide_ids2)

END_SECTION

START_SECTION(([EXTRA] XLMS data labeled cross-linker))
  vector<ProteinIdentification> protein_ids;
  vector<PeptideIdentification> peptide_ids;

  String input_file= OPENMS_GET_TEST_DATA_PATH("IdXML_XLMS_labelled.idXML");
  IdXMLFile().load(input_file, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b8]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[20].annotation, "[alpha|xi$b9]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].charge, 3)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].annotation, "[alpha|xi$y8]")

END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
