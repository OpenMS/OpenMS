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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(MascotXMLFile, "$Id: MascotXMLFile_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

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

START_SECTION((MascotXMLFile()))
	ptr = new MascotXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((void load(const String &filename, ProteinIdentification &protein_identification, std::vector< PeptideIdentification > &id_data)))

	xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_1.mascotXML"),
							protein_identification, 
				   		peptide_identifications);
				   		
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
	TEST_EQUAL(search_parameters.fixed_modifications[0], "Carboxymethyl (C)");
	TEST_EQUAL(search_parameters.fixed_modifications[1], "Deamidated (NQ)");
	TEST_EQUAL(search_parameters.fixed_modifications[2], "Guanidinyl (K)");
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
