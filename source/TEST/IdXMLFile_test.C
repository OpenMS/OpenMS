// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/IdXMLFile.h>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IdXMLFile* ptr;
START_SECTION((IdXMLFile()))
	ptr = new IdXMLFile();
	TEST_NOT_EQUAL(ptr,0)
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
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids, peptide_ids,document_id);
	
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
	TEST_EQUAL(protein_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
	TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
	TEST_EQUAL(protein_ids[0].getSearchParameters().enzyme,ProteinIdentification::TRYPSIN)
	TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"+1, +2")
	TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::AVERAGE)
	TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().peak_mass_tolerance,0.3)
	TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_tolerance,1.0)
	TEST_EQUAL((String)(protein_ids[0].getMetaValue("name")),"ProteinIdentification")
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
	TEST_EQUAL(peptide_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_REAL_SIMILAR((DoubleReal)(peptide_ids[0].getMetaValue("MZ")),675.9)
	TEST_REAL_SIMILAR((DoubleReal)(peptide_ids[0].getMetaValue("RT")),1234.5)
	TEST_EQUAL((UInt)(peptide_ids[0].getMetaValue("spectrum_reference")),17)
	TEST_EQUAL((String)(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
	TEST_EQUAL(peptide_ids[0].getHits().size(),2)
	//peptide hit 1
	TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),0.9)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(),"PEP1")
	TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),1)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getAABefore(),'A')
	TEST_EQUAL(peptide_ids[0].getHits()[0].getAAAfter(),'B')
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions().size(),2)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[0],"PROT1")
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[1],"PROT2")
	TEST_EQUAL((String)(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
	//peptide hit 2
	TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(),"PEP2")
	TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getProteinAccessions().size(),0)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[0].getHits()[1].getAAAfter(),' ')
	
	//peptide id 2
	TEST_EQUAL(peptide_ids[1].getScoreType(),"MOWSE")
	TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
	TEST_EQUAL(peptide_ids[1].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_EQUAL(peptide_ids[1].getHits().size(),2)
	//peptide hit 1
	TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),44.4)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(),"PEP3")
	TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),2)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getProteinAccessions().size(),0)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[1].getHits()[0].getAAAfter(),' ')
	
	//peptide hit 2
	TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),33.3)
	TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(),"PEP4")
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
	TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(),"PEP5")
	TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),1)
	TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions().size(),1)
	TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions()[0],"PROT3")
	TEST_EQUAL(peptide_ids[2].getHits()[0].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[2].getHits()[0].getAAAfter(),' ')
END_SECTION

START_SECTION(void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids, const String& document_id="") )
	
	//store and load data
	std::vector<ProteinIdentification> protein_ids, protein_ids2;
	std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
	String document_id, document_id2;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids2, peptide_ids2, document_id2);	
	String filename;
	NEW_TMP_FILE(filename)
	IdXMLFile().store(filename, protein_ids2, peptide_ids2, document_id2);	
	IdXMLFile().load(filename, protein_ids, peptide_ids, document_id);
	
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
	TEST_EQUAL(protein_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_EQUAL(protein_ids[0].getSearchParameters().db,"MSDB")
	TEST_EQUAL(protein_ids[0].getSearchParameters().db_version,"1.0")
	TEST_EQUAL(protein_ids[0].getSearchParameters().enzyme,ProteinIdentification::TRYPSIN)
	TEST_EQUAL(protein_ids[0].getSearchParameters().charges,"+1, +2")
	TEST_EQUAL(protein_ids[0].getSearchParameters().mass_type,ProteinIdentification::AVERAGE)
	TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().peak_mass_tolerance,0.3)
	TEST_REAL_SIMILAR(protein_ids[0].getSearchParameters().precursor_tolerance,1.0)
	TEST_EQUAL((String)(protein_ids[0].getMetaValue("name")),"ProteinIdentification")
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
	TEST_EQUAL(peptide_ids[0].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_REAL_SIMILAR((DoubleReal)(peptide_ids[0].getMetaValue("MZ")),675.9)
	TEST_REAL_SIMILAR((DoubleReal)(peptide_ids[0].getMetaValue("RT")),1234.5)
	TEST_EQUAL((UInt)(peptide_ids[0].getMetaValue("spectrum_reference")),17)
	TEST_EQUAL((String)(peptide_ids[0].getMetaValue("name")),"PeptideIdentification")
	TEST_EQUAL(peptide_ids[0].getHits().size(),2)
	//peptide hit 1
	TEST_REAL_SIMILAR(peptide_ids[0].getHits()[0].getScore(),0.9)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getSequence(),"PEP1")
	TEST_EQUAL(peptide_ids[0].getHits()[0].getCharge(),1)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getAABefore(),'A')
	TEST_EQUAL(peptide_ids[0].getHits()[0].getAAAfter(),'B')
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions().size(),2)
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[0],"PROT1")
	TEST_EQUAL(peptide_ids[0].getHits()[0].getProteinAccessions()[1],"PROT2")
	TEST_EQUAL((String)(peptide_ids[0].getHits()[0].getMetaValue("name")),"PeptideHit")
	//peptide hit 2
	TEST_REAL_SIMILAR(peptide_ids[0].getHits()[1].getScore(),1.4)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getSequence(),"PEP2")
	TEST_EQUAL(peptide_ids[0].getHits()[1].getCharge(),1)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getProteinAccessions().size(),0)
	TEST_EQUAL(peptide_ids[0].getHits()[1].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[0].getHits()[1].getAAAfter(),' ')
	//peptide id 2
	TEST_EQUAL(peptide_ids[1].getScoreType(),"MOWSE")
	TEST_EQUAL(peptide_ids[1].isHigherScoreBetter(),true)
	TEST_EQUAL(peptide_ids[1].getIdentifier(),"Mascot_2006-01-12T12:13:14")
	TEST_EQUAL(peptide_ids[1].getHits().size(),2)
	//peptide hit 1
	TEST_REAL_SIMILAR(peptide_ids[1].getHits()[0].getScore(),44.4)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getSequence(),"PEP3")
	TEST_EQUAL(peptide_ids[1].getHits()[0].getCharge(),2)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getProteinAccessions().size(),0)
	TEST_EQUAL(peptide_ids[1].getHits()[0].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[1].getHits()[0].getAAAfter(),' ')
	//peptide hit 2
	TEST_REAL_SIMILAR(peptide_ids[1].getHits()[1].getScore(),33.3)
	TEST_EQUAL(peptide_ids[1].getHits()[1].getSequence(),"PEP4")
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
	TEST_EQUAL(peptide_ids[2].getHits()[0].getSequence(),"PEP5")
	TEST_EQUAL(peptide_ids[2].getHits()[0].getCharge(),1)
	TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions().size(),1)
	TEST_EQUAL(peptide_ids[2].getHits()[0].getProteinAccessions()[0],"PROT3")
	TEST_EQUAL(peptide_ids[2].getHits()[0].getAABefore(),' ')
	TEST_EQUAL(peptide_ids[2].getHits()[0].getAAAfter(),' ')

	TEST_EQUAL(peptide_ids==peptide_ids2,true)
	TEST_EQUAL(protein_ids==protein_ids2,true)
END_SECTION


START_SECTION([EXTRA] static bool isValid(const String& filename))
	std::vector<ProteinIdentification> protein_ids, protein_ids2;
	std::vector<PeptideIdentification> peptide_ids, peptide_ids2;
	String filename;
	IdXMLFile f;

  //test if empty file is valid
	NEW_TMP_FILE(filename)
	f.store(filename, protein_ids2, peptide_ids2);	
  TEST_EQUAL(f.isValid(filename),true);	
	
	//test if full file is valid
	NEW_TMP_FILE(filename);
	String document_id;
	f.load(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML"), protein_ids2, peptide_ids2, document_id);
	protein_ids2[0].setMetaValue("stringvalue",String("bla"));
	protein_ids2[0].setMetaValue("intvalue",4711);
	protein_ids2[0].setMetaValue("floatvalue",5.3);
	f.store(filename, protein_ids2, peptide_ids2);	
  TEST_EQUAL(f.isValid(filename),true);
  
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
	
	TEST_EQUAL(protein_ids == protein_ids2, true)
	TEST_EQUAL(peptide_ids == peptide_ids2, true)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
