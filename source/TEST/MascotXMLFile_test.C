// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

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
vector<IdentificationData> identifications; 
vector<IdentificationData> identifications2; 
DateTime date;
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;
vector< pair<String, String> > references;

date.set("2006-03-09 11:31:52");

CHECK((MascotXMLFile()))
	ptr = new MascotXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((void load(const String& filename, ProteinIdentification& protein_identification, std::vector<IdentificationData>& id_data ) const throw(Exception::FileNotFound, Exception::ParseError)))

	xml_file.load("data/MascotXMLFile_test_1.mascotXML",
							protein_identification, 
				   		identifications);
	
	TEST_EQUAL(identifications.size(), 3)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(identifications[0].mz, 789.83)
	TEST_REAL_EQUAL(identifications[1].mz, 135.29)
	TEST_REAL_EQUAL(identifications[2].mz, 982.58)
	PRECISION(0.00001)	
	TEST_EQUAL(protein_identification.getProteinHits().size(), 2)
	TEST_EQUAL(protein_identification.getProteinHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identification.getProteinHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identification.getProteinHits()[0].getScore(), 619)
	TEST_REAL_EQUAL(protein_identification.getProteinHits()[1].getScore(), 293)
	TEST_EQUAL(protein_identification.getProteinHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identification.getProteinHits()[1].getScoreType(), "Mascot")
	
	protein_identification.getDateTime().get(date_string_1);
	TEST_EQUAL(date_string_1, "2006-03-09 11:31:52")
	
	TEST_REAL_EQUAL(identifications[0].id.getPeptideSignificanceThreshold(), 31.8621)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	
	peptide_hit = identifications[0].id.getPeptideHits()[0];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 2)
	TEST_EQUAL(references[0].second, "AAN17824")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	TEST_EQUAL(references[1].second, "GN1736")
	TEST_EQUAL(references[1].first, "2006-03-09 11:31:52")
	peptide_hit = identifications[0].id.getPeptideHits()[1];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0].second, "AAN17824")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	peptide_hit = identifications[1].id.getPeptideHits()[0];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0].second, "GN1736")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 43.9)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getDateTime() == date, true)	
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "HSKLSAK")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
