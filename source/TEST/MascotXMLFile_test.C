// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>


///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(MascotXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MascotXMLFile xml_file;
MascotXMLFile* ptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
DateTime date;
PeptideHit peptide_hit;
vector<String> references;

date.set("2006-03-09 11:31:52");

MascotXMLFile* nullPointer = 0;
START_SECTION((MascotXMLFile()))
	ptr = new MascotXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((void load(const String &filename, ProteinIdentification &protein_identification, std::vector< PeptideIdentification > &id_data)))
  
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_1.mascotXML"),
							protein_identification, 
				   		peptide_identifications);
  
  {
	ProteinIdentification::SearchParameters search_parameters = protein_identification.getSearchParameters();
	TEST_EQUAL(search_parameters.missed_cleavages, 1);
	TEST_EQUAL(search_parameters.taxonomy, ". . Eukaryota (eucaryotes)");
	TEST_EQUAL(search_parameters.mass_type, ProteinIdentification::AVERAGE);
	TEST_EQUAL(search_parameters.enzyme, ProteinIdentification::TRYPSIN);
	TEST_EQUAL(search_parameters.db, "MSDB_chordata");
	TEST_EQUAL(search_parameters.db_version, "MSDB_chordata_20070910.fasta");
	TEST_EQUAL(search_parameters.peak_mass_tolerance, 0.2);
	TEST_EQUAL(search_parameters.precursor_tolerance, 1.4);
	TEST_EQUAL(search_parameters.charges, "1+, 2+ and 3+");
  TEST_EQUAL(search_parameters.fixed_modifications.size(), 3);
	TEST_EQUAL(search_parameters.fixed_modifications[0], "Carboxymethyl (C)");
	TEST_EQUAL(search_parameters.fixed_modifications[1], "Deamidated (NQ)");
	TEST_EQUAL(search_parameters.fixed_modifications[2], "Guanidinyl (K)");
  TEST_EQUAL(search_parameters.variable_modifications.size(), 3);
	TEST_EQUAL(search_parameters.variable_modifications[0], "Acetyl (Protein N-term)");
	TEST_EQUAL(search_parameters.variable_modifications[1], "Biotin (K)");
	TEST_EQUAL(search_parameters.variable_modifications[2], "Carbamyl (K)");
	TEST_EQUAL(peptide_identifications.size(), 3)
	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(peptide_identifications[0].getMetaValue("MZ"), 789.83)
	TEST_REAL_SIMILAR(peptide_identifications[1].getMetaValue("MZ"), 135.29)
	TEST_REAL_SIMILAR(peptide_identifications[2].getMetaValue("MZ"), 982.58)
	TOLERANCE_ABSOLUTE(0.00001)	
	TEST_EQUAL(protein_identification.getHits().size(), 2)
	TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736")
	TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 619)
	TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 293)
	TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime().get(), "2006-03-09 11:31:52")
	
	TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 31.8621)
	TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
	
	peptide_hit = peptide_identifications[0].getHits()[0];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 2)
	TEST_EQUAL(references[0], "AAN17824")
	TEST_EQUAL(references[1], "GN1736")
	peptide_hit = peptide_identifications[0].getHits()[1];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0], "AAN17824")
	peptide_hit = peptide_identifications[1].getHits()[0];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0], "GN1736")
	
	TEST_EQUAL(peptide_identifications[1].getHits().size(), 1)
	TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_SIMILAR(peptide_identifications[1].getHits()[0].getScore(), 43.9)
	TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime() == date, true)	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATN(MOD:00565)FK(MOD:00445)")
	TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK(MOD:00445)")
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "HSK(MOD:00445)LSAK(MOD:00445)")
  }


  /// for new MascotXML 2.1 as used by Mascot Server 2.3
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_2.mascotXML"),
							protein_identification, 
				   		peptide_identifications);
  {
  ProteinIdentification::SearchParameters search_parameters = protein_identification.getSearchParameters();
  TEST_EQUAL(search_parameters.missed_cleavages, 7);
  TEST_EQUAL(search_parameters.taxonomy, "All entries");
  TEST_EQUAL(search_parameters.mass_type, ProteinIdentification::MONOISOTOPIC);
  TEST_EQUAL(search_parameters.enzyme, ProteinIdentification::TRYPSIN);
  TEST_EQUAL(search_parameters.db, "IPI_human");
  TEST_EQUAL(search_parameters.db_version, "ipi.HUMAN.v3.61.fasta");
  TEST_EQUAL(search_parameters.peak_mass_tolerance, 0.3);
  TEST_EQUAL(search_parameters.precursor_tolerance, 3);
  TEST_EQUAL(search_parameters.charges, "");
  TEST_EQUAL(search_parameters.fixed_modifications.size(), 1);
  TEST_EQUAL(search_parameters.fixed_modifications[0], "Carbamidomethyl (C)");
  TEST_EQUAL(search_parameters.variable_modifications.size(), 3);
  TEST_EQUAL(search_parameters.variable_modifications[0], "Oxidation (M)");
  TEST_EQUAL(search_parameters.variable_modifications[1], "Acetyl (N-term)");
  TEST_EQUAL(search_parameters.variable_modifications[2], "Phospho (Y)");
  // not necessarily equal to numQueries as some hits might not be contained, e.g. peptide's might start with <peptide rank="10"...> so 9 peptides are missing
  // thus empty peptides are removed (see MascotXMLFile.C::load() ) after the handler() call
  TEST_EQUAL(peptide_identifications.size(), 1112)
  TOLERANCE_ABSOLUTE(0.0001)
  TEST_REAL_SIMILAR(peptide_identifications[0].getMetaValue("MZ"), 304.6967)
  TEST_REAL_SIMILAR(peptide_identifications[1].getMetaValue("MZ"), 314.1815)
  TEST_REAL_SIMILAR(peptide_identifications[1111].getMetaValue("MZ"), 583.7948)
  TOLERANCE_ABSOLUTE(0.00001)	
  TEST_EQUAL(protein_identification.getHits().size(), 66)
  TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "IPI00745872")
  TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "IPI00908876")
  TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 122)
  TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 122)
  TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
  TEST_EQUAL(protein_identification.getDateTime().get(), "2011-06-24 19:34:54")

  TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 5)
  TEST_EQUAL(peptide_identifications[0].getHits().size(), 1)

  peptide_hit = peptide_identifications[0].getHits()[0];
  TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 0)
  references = peptide_identifications[34].getHits()[0].getProteinAccessions(); // corresponds to <peptide query="35" ...>
  ABORT_IF(references.size() != 5)
  TEST_EQUAL(references[0], "IPI00745872")
  TEST_EQUAL(references[4], "IPI00878517")

  TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 5.34)
  TEST_REAL_SIMILAR(peptide_identifications[49].getHits()[0].getScore(), 14.83)
  TEST_REAL_SIMILAR(peptide_identifications[49].getHits()[1].getScore(), 17.5)
  TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
  TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
  TEST_EQUAL(protein_identification.getDateTime().get() == "2011-06-24 19:34:54", true)	
  TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "VVFIK")
  TEST_EQUAL(peptide_identifications[49].getHits()[0].getSequence(), "LASYLDK")
  TEST_EQUAL(peptide_identifications[49].getHits()[1].getSequence(), "(Acetyl)AAFESDK")
  //for (int i=520;i<540;++i) std::cerr << "i: " << i << " " << peptide_identifications[i].getHits()[0].getSequence() << "\n";
  TEST_EQUAL(peptide_identifications[522].getHits()[0].getSequence(), "(Acetyl)GALM(Oxidation)NEIQAAK")
  TEST_EQUAL(peptide_identifications[67].getHits()[0].getSequence(), "SHY(Phospho)GGSR")
  }

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_3.mascotXML"),
    protein_identification, 
    peptide_identifications);
  {
    std::vector<ProteinIdentification> pids;
    pids.push_back(protein_identification);
    String filename;
    NEW_TMP_FILE(filename)
    IdXMLFile().store(filename, pids, peptide_identifications);
    FuzzyStringComparator fuzzy;
    fuzzy.setWhitelist(StringList::create("<?xml-stylesheet"));
    fuzzy.setAcceptableAbsolute(0.0001);
    bool result = fuzzy.compareFiles(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_out_3.idXML"), filename);
    TEST_EQUAL(result, true);
  }

END_SECTION

START_SECTION((void load(const String &filename, ProteinIdentification &protein_identification, std::vector< PeptideIdentification > &id_data, std::map< String, std::vector< AASequence > > &peptides)))
	std::map<String, vector<AASequence> > modified_peptides;
	AASequence aa_sequence_1;
	AASequence aa_sequence_2;
	AASequence aa_sequence_3;
	vector<AASequence> temp;
	
	aa_sequence_1.setStringSequence("LHASGITVTEIPVTATNFK");
	aa_sequence_1.setModification(16, "Deamidated");
	aa_sequence_2.setStringSequence("MRSLGYVAVISAVATDTDK");
	aa_sequence_2.setModification(2, "Phospho");
	aa_sequence_3.setStringSequence("HSKLSAK");
	aa_sequence_3.setModification(4, "Phospho");
	temp.push_back(aa_sequence_1);
	temp.push_back(aa_sequence_2);
	modified_peptides.insert(make_pair("789.83", temp));
	temp.clear();
	temp.push_back(aa_sequence_3);
	modified_peptides.insert(make_pair("135.29", temp));
		
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_1.mascotXML"),
								protein_identification, 
					   		peptide_identifications,
					   		modified_peptides);
	
	TEST_EQUAL(peptide_identifications.size(), 3)
	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(peptide_identifications[0].getMetaValue("MZ"), 789.83)
	TEST_REAL_SIMILAR(peptide_identifications[1].getMetaValue("MZ"), 135.29)
	TEST_REAL_SIMILAR(peptide_identifications[2].getMetaValue("MZ"), 982.58)
	TOLERANCE_ABSOLUTE(0.00001)	
	TEST_EQUAL(protein_identification.getHits().size(), 2)
	TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736")
	TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 619)
	TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 293)
	TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime().get(), "2006-03-09 11:31:52")
	
	TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 31.8621)
	TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
	
	peptide_hit = peptide_identifications[0].getHits()[0];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 2)
	TEST_EQUAL(references[0], "AAN17824")
	TEST_EQUAL(references[1], "GN1736")
	peptide_hit = peptide_identifications[0].getHits()[1];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0], "AAN17824")
	peptide_hit = peptide_identifications[1].getHits()[0];
	references = peptide_hit.getProteinAccessions();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0], "GN1736")
	
	TEST_EQUAL(peptide_identifications[1].getHits().size(), 1)
	TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_SIMILAR(peptide_identifications[1].getHits()[0].getScore(), 43.9)
	TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime() == date, true)	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), aa_sequence_1)
	TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), aa_sequence_2)
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), aa_sequence_3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
