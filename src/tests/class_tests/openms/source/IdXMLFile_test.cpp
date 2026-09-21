// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/TestFileValidation.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <fstream>
#include <sstream>

///////////////////////////

using namespace OpenMS;

namespace
{
  // registers a tool-defined modification; ModificationsDB is process-wide, so every section uses its own name
  const ResidueModification* defineMod4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return ModificationsDB::getInstance()->registerDefinition(d);
  }

  // a definition record for a name that is NOT registered in this process
  std::string freshRecord4b(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return d.toDefinitionString();
  }

  std::string slurp4b(const std::string& path)
  {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  bool fileContains4b(const std::string& path, const std::string& needle)
  {
    return slurp4b(path).find(needle) != std::string::npos;
  }

  // first occurrence only; returns false when @p from is absent
  bool replaceInFile4b(const std::string& path, const std::string& from, const std::string& to)
  {
    std::string s = slurp4b(path);
    const std::size_t pos = s.find(from);
    if (pos == std::string::npos) return false;
    s.replace(pos, from.size(), to);
    std::ofstream out(path);
    out << s;
    return true;
  }
}

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

START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids) )
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(),2)
  TEST_EQUAL(peptide_ids.size(),3)
END_SECTION


START_SECTION(void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids, std::string& document_id) )
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  std::string document_id;
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
  TEST_EQUAL(StringUtils::hasPrefix(protein_ids[0].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
  TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"1, 2")
  TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::PeakMassType::AVERAGE)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().fragment_mass_tolerance,0.3)
  TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_mass_tolerance,1.0)
  TEST_EQUAL(std::string(protein_ids[0].getMetaValue("name")),"ProteinIdentification")

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
  TEST_EQUAL(std::string(protein_ids[0].getHits()[0].getMetaValue("name")),"ProteinHit")
  //protein hit 2
  TEST_REAL_SIMILAR(protein_ids[0].getHits()[1].getScore(),24.4)
  TEST_EQUAL(protein_ids[0].getHits()[1].getAccession(),"PROT2")
  TEST_EQUAL(protein_ids[0].getHits()[1].getSequence(),"ABCDEFG")

  //peptide id 1
  TEST_EQUAL(peptide_ids[0].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[0].isHigherScoreBetter(),false)
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[0].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
  TEST_REAL_SIMILAR(peptide_ids[0].getMZ(),675.9)
  TEST_REAL_SIMILAR(peptide_ids[0].getRT(),1234.5)
  TEST_EQUAL((peptide_ids[0].getSpectrumReference()),"17")
  TEST_EQUAL(std::string(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
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
  TEST_EQUAL(std::string(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
  //peptide hit 2
  TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
  vector<PeptideEvidence> pes1 = peptide_ids[0].getHits()[1].getPeptideEvidences();
  TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(), AASequence::fromString("PEPTIDERR"))
  TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
  TEST_EQUAL(pes1.size(),0)
  //peptide id 2
  TEST_EQUAL(peptide_ids[1].getScoreType(),"MOWSE")
  TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[1].getIdentifier(), "Mascot_2006-01-12T12:13:14"), true)
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
  TEST_EQUAL(StringUtils::hasPrefix(protein_ids[1].getIdentifier(), "Mascot_2007-01-12T12:13:14"), true)
  TEST_EQUAL(protein_ids[1].getSearchParameters().db,"MSDB")
  TEST_EQUAL(protein_ids[1].getSearchParameters().db_version,"1.1")
  TEST_EQUAL(protein_ids[1].getSearchParameters().charges,"1, 2, 3")
  TEST_EQUAL(protein_ids[1].getSearchParameters().mass_type,ProteinIdentification::PeakMassType::MONOISOTOPIC)
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
  TEST_EQUAL(StringUtils::hasPrefix(peptide_ids[2].getIdentifier(), "Mascot_2007-01-12T12:13:14"), true)
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

START_SECTION(void store(std::string filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids, const std::string& document_id="") )

  // load, store, and reload data
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string document_id, document_id2;
  std::string target_file = OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML");
  IdXMLFile().load(target_file, protein_ids2, peptide_ids2, document_id2);

  std::string actual_file;
  NEW_TMP_FILE(actual_file)
  IdXMLFile().store(actual_file, protein_ids2, peptide_ids2, document_id2);

  FuzzyStringComparator fuzzy;
  fuzzy.setWhitelist(ListUtils::create<std::string>("<?xml-stylesheet"));
  fuzzy.setAcceptableAbsolute(0.0001);
  bool result = fuzzy.compareFiles(actual_file, target_file);
  TEST_EQUAL(result, true);
END_SECTION


START_SECTION([EXTRA] static bool isValid(const std::string& filename))
  std::vector<ProteinIdentification> protein_ids, protein_ids2;
  PeptideIdentificationList peptide_ids, peptide_ids2;
  std::string filename;
  IdXMLFile f;

  //test if empty file is valid
  NEW_TMP_FILE(filename)
  f.store(filename, protein_ids2, peptide_ids2);
  TEST_EQUAL(f.isValid(filename, std::cerr),true);

  //test if full file is valid
  NEW_TMP_FILE(filename);
  std::string document_id;
  f.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids2, peptide_ids2, document_id);
  protein_ids2[0].setMetaValue("stringvalue",std::string("bla"));
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
  PeptideIdentificationList peptide_ids;
  id_xmlfile.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_no_proteinhits.idXML"), protein_ids, peptide_ids);

  TEST_EQUAL(protein_ids.size(), 1)
  TEST_EQUAL(protein_ids[0].getHits().size(), 0)
  TEST_EQUAL(peptide_ids.size(), 10)
  TEST_EQUAL(peptide_ids[0].getHits().size(), 1)

  std::string filename;
  NEW_TMP_FILE(filename)
  id_xmlfile.store(filename , protein_ids, peptide_ids);

  vector<ProteinIdentification> protein_ids2;
  PeptideIdentificationList peptide_ids2;
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
  PeptideIdentificationList peptide_ids;

  std::string input_file= OPENMS_GET_TEST_DATA_PATH("IdXML_XLMS_labelled.idXML");
  IdXMLFile().load(input_file, protein_ids, peptide_ids);

  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[0].charge, 1)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[1].annotation, "[alpha|ci$b2]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[8].annotation, "[alpha|xi$b8]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[20].annotation, "[alpha|xi$b9]")
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].charge, 3)
  TEST_EQUAL(peptide_ids[0].getHits()[0].getPeakAnnotations()[25].annotation, "[alpha|xi$y8]")

END_SECTION

START_SECTION([EXTRA] Compressed file writing - gzip round-trip)
  // Load reference data
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  // Store as gzip-compressed file
  std::string tmp_gz;
  NEW_TMP_FILE(tmp_gz);
  tmp_gz += ".gz";
  IdXMLFile().store(tmp_gz, protein_ids, peptide_ids);

  // Load back from compressed file
  std::vector<ProteinIdentification> protein_ids_gz;
  PeptideIdentificationList peptide_ids_gz;
  IdXMLFile().load(tmp_gz, protein_ids_gz, peptide_ids_gz);

  // Verify round-trip integrity
  TEST_EQUAL(protein_ids_gz.size(), protein_ids.size())
  TEST_EQUAL(peptide_ids_gz.size(), peptide_ids.size())
  TEST_EQUAL(protein_ids_gz[0].getHits().size(), protein_ids[0].getHits().size())
  TEST_EQUAL(protein_ids_gz[0].getHits()[0].getAccession(), protein_ids[0].getHits()[0].getAccession())
END_SECTION

START_SECTION([EXTRA] Compressed file writing - bzip2 round-trip)
  // Load reference data
  std::vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids);

  // Store as bzip2-compressed file
  std::string tmp_bz2;
  NEW_TMP_FILE(tmp_bz2);
  tmp_bz2 += ".bz2";
  IdXMLFile().store(tmp_bz2, protein_ids, peptide_ids);

  // Load back from compressed file
  std::vector<ProteinIdentification> protein_ids_bz2;
  PeptideIdentificationList peptide_ids_bz2;
  IdXMLFile().load(tmp_bz2, protein_ids_bz2, peptide_ids_bz2);

  // Verify round-trip integrity
  TEST_EQUAL(protein_ids_bz2.size(), protein_ids.size())
  TEST_EQUAL(peptide_ids_bz2.size(), peptide_ids.size())
  TEST_EQUAL(protein_ids_bz2[0].getHits().size(), protein_ids[0].getHits().size())
  TEST_EQUAL(protein_ids_bz2[0].getHits()[0].getAccession(), protein_ids[0].getHits()[0].getAccession())
END_SECTION

START_SECTION([EXTRA] store/load - a tool-defined modification travels with its definition)
{
  TEST_TRUE(defineMod4b("TestIdXML:Adduct", 'K', "C9H11N2O8P") != nullptr)
  ProteinIdentification prot;
  prot.setIdentifier("run4b");
  prot.setDateTime(DateTime::now());
  PeptideIdentification pep;
  pep.setIdentifier("run4b");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AEADNLDDK(TestIdXML:Adduct)K"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_TRUE(fileContains4b(filename, "name=\"modification_definitions\""))
  TEST_TRUE(fileContains4b(filename, "1|TestIdXML:Adduct|TestIdXML:Adduct (K)|"))

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_EQUAL(prots_in.size(), 1)
  TEST_EQUAL(peps_in.size(), 1)
  if (prots_in.size() == 1 && peps_in.size() == 1 && !peps_in[0].getHits().empty())
  {
    TEST_TRUE(prots_in[0].getSearchParameters().metaValueExists(Constants::UserParam::MODIFICATION_DEFINITIONS))
    const AASequence& back = peps_in[0].getHits()[0].getSequence();
    TEST_EQUAL(back.toString(), "AEADNLDDK(TestIdXML:Adduct)K")
    TEST_EQUAL(back.getFormula(Residue::Full, 0).toString(), "C54H86N15O28P1")
  }
}
END_SECTION

START_SECTION([EXTRA] store - runs with equal parameters share one block that carries both definition sets)
{
  TEST_TRUE(defineMod4b("TestIdXML:RunA", 'K', "C2H2O") != nullptr)
  TEST_TRUE(defineMod4b("TestIdXML:RunB", 'R', "CH2") != nullptr)
  std::vector<ProteinIdentification> prots(2);
  prots[0].setIdentifier("runA");
  prots[0].setDateTime(DateTime::now());
  prots[1].setIdentifier("runB");
  prots[1].setDateTime(DateTime::now());
  PeptideIdentificationList peps(2);
  peps[0].setIdentifier("runA");
  peps[1].setIdentifier("runB");
  PeptideHit ha, hb;
  ha.setSequence(AASequence::fromString("PEPK(TestIdXML:RunA)IDE"));
  hb.setSequence(AASequence::fromString("PEPR(TestIdXML:RunB)IDE"));
  peps[0].insertHit(ha);
  peps[1].insertHit(hb);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  const std::string text = slurp4b(filename);
  Size blocks = 0;
  for (std::size_t pos = text.find("<SearchParameters "); pos != std::string::npos; pos = text.find("<SearchParameters ", pos + 1)) ++blocks;
  TEST_EQUAL(blocks, 1)
  TEST_TRUE(text.find("1|TestIdXML:RunA|TestIdXML:RunA (K)|") != std::string::npos)
  TEST_TRUE(text.find("1|TestIdXML:RunB|TestIdXML:RunB (R)|") != std::string::npos)

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_EQUAL(prots_in.size(), 2)
  TEST_EQUAL(peps_in.size(), 2)
  if (peps_in.size() == 2 && !peps_in[0].getHits().empty() && !peps_in[1].getHits().empty())
  {
    TEST_EQUAL(peps_in[0].getHits()[0].getSequence().toString(), "PEPK(TestIdXML:RunA)IDE")
    TEST_EQUAL(peps_in[1].getHits()[0].getSequence().toString(), "PEPR(TestIdXML:RunB)IDE")
  }
}
END_SECTION

START_SECTION([EXTRA] store - a defined variable modification on no hit keeps its definition)
{
  TEST_TRUE(defineMod4b("TestIdXML:VarOnly", 'S', "HPO3") != nullptr)
  ProteinIdentification prot;
  prot.setIdentifier("run4b_var");
  prot.setDateTime(DateTime::now());
  ProteinIdentification::SearchParameters sp;
  sp.variable_modifications.push_back("TestIdXML:VarOnly (S)");
  prot.setSearchParameters(sp);
  PeptideIdentification pep;
  pep.setIdentifier("run4b_var");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_TRUE(fileContains4b(filename, "1|TestIdXML:VarOnly|TestIdXML:VarOnly (S)|"))
}
END_SECTION

START_SECTION([EXTRA] load - definitions are registered before the sequences are parsed)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  TEST_FALSE(db->hasDefinedModification("TestIdXML:Fresh"))
  ProteinIdentification prot;
  prot.setIdentifier("run4b_fresh");
  prot.setDateTime(DateTime::now());
  ProteinIdentification::SearchParameters sp;
  sp.setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, freshRecord4b("TestIdXML:Fresh", 'K', "C2H2O"));
  prot.setSearchParameters(sp);
  PeptideIdentification pep;
  pep.setIdentifier("run4b_fresh");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTKIDE"));
  pep.insertHit(hit);
  std::vector<ProteinIdentification> prots(1, prot);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  std::string filename;
  NEW_TMP_FILE(filename)
  IdXMLFile().store(filename, prots, peps);
  TEST_FALSE(db->hasDefinedModification("TestIdXML:Fresh")) // storing registers nothing
  // a hit using the not-yet-registered name, as a file from another process would carry it
  TEST_TRUE(replaceInFile4b(filename, "sequence=\"PEPTKIDE\"", "sequence=\"PEPTK(TestIdXML:Fresh)IDE\""))

  std::vector<ProteinIdentification> prots_in;
  PeptideIdentificationList peps_in;
  IdXMLFile().load(filename, prots_in, peps_in);
  TEST_TRUE(db->hasDefinedModification("TestIdXML:Fresh"))
  TEST_EQUAL(peps_in.size(), 1)
  if (peps_in.size() == 1 && !peps_in[0].getHits().empty())
  {
    const AASequence& back = peps_in[0].getHits()[0].getSequence();
    TEST_EQUAL(back.toString(), "PEPTK(TestIdXML:Fresh)IDE")
    TEST_REAL_SIMILAR(back.getMonoWeight(), AASequence::fromString("PEPTKIDE").getMonoWeight() + EmpiricalFormula("C2H2O").getMonoWeight())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/// check the temporary files written above against their XML schema (types without a validator are skipped)
VALIDATE_TMP_FILES

END_TEST
