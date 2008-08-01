// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;
vector<String> references;

date.set("2006-03-09 11:31:52");

CHECK((MascotXMLFile()))
	ptr = new MascotXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((void load(const String &filename, ProteinIdentification &protein_identification, std::vector<PeptideIdentification> &id_data) const))

	xml_file.load("data/MascotXMLFile_test_1.mascotXML",
							protein_identification, 
				   		peptide_identifications);
	
	TEST_EQUAL(peptide_identifications.size(), 3)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(peptide_identifications[0].getMetaValue("MZ"), 789.83)
	TEST_REAL_EQUAL(peptide_identifications[1].getMetaValue("MZ"), 135.29)
	TEST_REAL_EQUAL(peptide_identifications[2].getMetaValue("MZ"), 982.58)
	PRECISION(0.00001)	
	TEST_EQUAL(protein_identification.getHits().size(), 2)
	TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identification.getHits()[0].getScore(), 619)
	TEST_REAL_EQUAL(protein_identification.getHits()[1].getScore(), 293)
	TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
	
	protein_identification.getDateTime().get(date_string_1);
	TEST_EQUAL(date_string_1, "2006-03-09 11:31:52")
	
	TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 31.8621)
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
	TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(peptide_identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(peptide_identifications[1].getHits()[0].getScore(), 43.9)
	TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime() == date, true)	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "HSKLSAK")
RESULT

CHECK((void load(const String &filename, ProteinIdentification &protein_identification, std::vector<PeptideIdentification> &id_data, std::map<String, std::vector<AASequence> > &peptides) const))	
	std::map<String, vector<AASequence> > modified_peptides;
	AASequence aa_sequence_1;
	AASequence aa_sequence_2;
	AASequence aa_sequence_3;
	vector<AASequence> temp;
	
	aa_sequence_1.setStringSequence("LHASGITVTEIPVTATNFK");
	aa_sequence_1.setModification(6, "Deamidated");
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
		
	xml_file.load("data/MascotXMLFile_test_1.mascotXML",
								protein_identification, 
					   		peptide_identifications,
					   		modified_peptides);
	
	TEST_EQUAL(peptide_identifications.size(), 3)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(peptide_identifications[0].getMetaValue("MZ"), 789.83)
	TEST_REAL_EQUAL(peptide_identifications[1].getMetaValue("MZ"), 135.29)
	TEST_REAL_EQUAL(peptide_identifications[2].getMetaValue("MZ"), 982.58)
	PRECISION(0.00001)	
	TEST_EQUAL(protein_identification.getHits().size(), 2)
	TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identification.getHits()[0].getScore(), 619)
	TEST_REAL_EQUAL(protein_identification.getHits()[1].getScore(), 293)
	TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
	
	protein_identification.getDateTime().get(date_string_1);
	TEST_EQUAL(date_string_1, "2006-03-09 11:31:52")
	
	TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 31.8621)
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
	TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(peptide_identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(peptide_identifications[1].getHits()[0].getScore(), 43.9)
	TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getDateTime() == date, true)	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), aa_sequence_1)
	TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), aa_sequence_2)
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), aa_sequence_3)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
