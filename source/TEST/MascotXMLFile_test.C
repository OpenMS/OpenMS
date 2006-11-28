// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
vector<Identification> identifications; 
vector<float> precursor_retention_times;
vector<float> precursor_mz_values;
vector<Identification> identifications2; 
vector<float> precursor_retention_times2;
vector<float> precursor_mz_values2;
DateTime date;
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;
vector< pair<String, String> > references;

date.set("2006-03-09 11:31:52");

CHECK(MascotXMLFile())
	ptr = new MascotXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((void load(const String& filename, ProteinIdentification* protein_identification, std::vector<Identification>* identifications, std::vector<float>* precursor_retention_times, std::vector<float>* precursor_mz_values) const throw(Exception::FileNotFound, Exception::FileNotReadable, Exception::FileEmpty, Exception::ParseError)))

	xml_file.load("data/MascotXMLFile_test_1.mascotXML",
							&protein_identification, 
				   		&identifications, 
							&precursor_retention_times, 
							&precursor_mz_values);
	
	TEST_EQUAL(identifications.size(), 3)
//	TEST_EQUAL(precursor_retention_times.size(), 2)
	TEST_EQUAL(precursor_mz_values.size(), 3)
//	TEST_EQUAL(precursor_retention_times[0], 120)
//	TEST_EQUAL(precursor_retention_times[1], 120)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(precursor_mz_values[0], 789.83)
	TEST_REAL_EQUAL(precursor_mz_values[1], 135.29)
	TEST_REAL_EQUAL(precursor_mz_values[2], 982.58)
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
	
	TEST_REAL_EQUAL(identifications[0].getPeptideSignificanceThreshold(), 31.8621)
	TEST_EQUAL(identifications[0].getCharge(), 3)
	TEST_EQUAL(identifications[1].getCharge(), 2)
	TEST_EQUAL(identifications[0].getPeptideHits().size(), 2)
	
	peptide_hit = identifications[0].getPeptideHits()[0];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 2)
	TEST_EQUAL(references[0].second, "AAN17824")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	TEST_EQUAL(references[1].second, "GN1736")
	TEST_EQUAL(references[1].first, "2006-03-09 11:31:52")
	peptide_hit = identifications[0].getPeptideHits()[1];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0].second, "AAN17824")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	peptide_hit = identifications[1].getPeptideHits()[0];
	references = peptide_hit.getProteinIndices();
	TEST_EQUAL(references.size(), 1)
	TEST_EQUAL(references[0].second, "GN1736")
	TEST_EQUAL(references[0].first, "2006-03-09 11:31:52")
	
	TEST_EQUAL(identifications[1].getPeptideHits().size(), 1)
	TEST_REAL_EQUAL(identifications[0].getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].getPeptideHits()[0].getScore(), 43.9)
	TEST_EQUAL(identifications[0].getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getDateTime() == date, true)	
	TEST_EQUAL(identifications[0].getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].getPeptideHits()[0].getSequence(), "HSKLSAK")
RESULT

CHECK(~MascotXMLFile())
	delete ptr;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
