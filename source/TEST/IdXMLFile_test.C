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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IdXMLFile* ptr;
IdXMLFile xml_file;

CHECK((IdXMLFile()))
	ptr = new IdXMLFile();
RESULT

CHECK((void load(const String& filename, std::vector<Identification>& protein_identifications, std::vector<PeptideIdentification>& id_data) const throw(Exception::FileNotFound, Exception::ParseError)))

	vector<Identification> protein_identifications; 
	vector<PeptideIdentification> identifications;
	String date_string = "2007-05-02 21:22:05"; 

	xml_file.load("data/IdXMLFile_test.idXML",
							protein_identifications, 
				   		identifications);
	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(DoubleReal(identifications[0].getMetaValue("RT")), 120)
	TEST_EQUAL(DoubleReal(identifications[1].getMetaValue("RT")), 150)
	TEST_EQUAL(DoubleReal(identifications[2].getMetaValue("RT")), 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(DoubleReal(identifications[0].getMetaValue("MZ")), 789.83)
	TEST_REAL_EQUAL(DoubleReal(identifications[1].getMetaValue("MZ")), 135.29)
	TEST_REAL_EQUAL(DoubleReal(identifications[2].getMetaValue("MZ")), 982.58)
	TEST_REAL_EQUAL(identifications[0].getSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].getSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].getSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].getHits().size(), 2)
	TEST_EQUAL(identifications[1].getHits().size(), 1)
	TEST_EQUAL(identifications[2].getHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].getHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].getHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].getHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].getHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].getHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].getHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(DoubleReal(identifications[0].getHits()[0].getMetaValue("predicted_RT")), 122.5)
	TEST_EQUAL(DoubleReal(identifications[0].getHits()[1].getMetaValue("predicted_RT")), 122.5)
	TEST_EQUAL(DoubleReal(identifications[1].getHits()[0].getMetaValue("predicted_RT")), 151.5)
	TEST_EQUAL(DoubleReal(identifications[2].getHits()[0].getMetaValue("predicted_RT")), 159.5)
	TEST_EQUAL(DoubleReal(identifications[2].getHits()[1].getMetaValue("predicted_RT")), 159.5)
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getScoreType(), "Mascot")
	DateTime temp_date = protein_identifications[0].getDateTime();
	String temp_string;
	temp_date.get(temp_string);
	TEST_EQUAL(temp_string, date_string)
	Identification::SearchParameters sp = protein_identifications[0].getSearchParameters();
	TEST_EQUAL(sp.db, "MSDB")
	TEST_EQUAL(sp.db_version, "1.0")
	TEST_EQUAL(sp.taxonomy, "vertebrates")
RESULT

CHECK((void store(String filename, const std::vector<Identification>& protein_identifications, const std::vector<PeptideIdentification>& id_data) const throw(Exception::UnableToCreateFile)))
												
	vector<Identification> protein_identifications; 
	vector<PeptideIdentification> identifications; 

	String temp_filename = "data/IdXMLFile_test_2.idXML";

	NEW_TMP_FILE(temp_filename)
	xml_file.load("data/IdXMLFile_test.idXML", 
							protein_identifications, 
				   		identifications);
	xml_file.store(temp_filename, 
							    protein_identifications, 
				   				identifications);
	xml_file.load(temp_filename, 
							protein_identifications, 
				   		identifications);
	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(DoubleReal(identifications[0].getMetaValue("RT")), 120)
	TEST_EQUAL(DoubleReal(identifications[1].getMetaValue("RT")), 150)
	TEST_EQUAL(DoubleReal(identifications[2].getMetaValue("RT")), 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(DoubleReal(identifications[0].getMetaValue("MZ")), 789.83)
	TEST_REAL_EQUAL(DoubleReal(identifications[1].getMetaValue("MZ")), 135.29)
	TEST_REAL_EQUAL(DoubleReal(identifications[2].getMetaValue("MZ")), 982.58)
	TEST_REAL_EQUAL(identifications[0].getSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].getSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].getSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].getHits().size(), 2)
	TEST_EQUAL(identifications[1].getHits().size(), 1)
	TEST_EQUAL(identifications[2].getHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].getHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].getHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].getHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].getHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].getHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].getHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getScoreType(), "Mascot")
									
RESULT

CHECK((void store(String filename, const std::vector<Identification>& protein_identifications, const std::vector<PeptideIdentification>& id_data, const std::map<String, double>& predicted_retention_times) const throw(Exception::UnableToCreateFile)))
												
	vector<Identification> protein_identifications; 
	vector<PeptideIdentification> identifications; 

	String temp_filename = "data/IdXMLFile_test_2.idXML";
	NEW_TMP_FILE(temp_filename)
	
	xml_file.load("data/IdXMLFile_test.idXML", 
							protein_identifications, 
				   		identifications);
	xml_file.store(temp_filename, 
									protein_identifications, 
				   				identifications);

	xml_file.load(temp_filename,
							protein_identifications, 
				   		identifications);

	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(DoubleReal(identifications[0].getMetaValue("RT")), 120)
	TEST_EQUAL(DoubleReal(identifications[1].getMetaValue("RT")), 150)
	TEST_EQUAL(DoubleReal(identifications[2].getMetaValue("RT")), 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(DoubleReal(identifications[0].getMetaValue("MZ")), 789.83)
	TEST_REAL_EQUAL(DoubleReal(identifications[1].getMetaValue("MZ")), 135.29)
	TEST_REAL_EQUAL(DoubleReal(identifications[2].getMetaValue("MZ")), 982.58)
	TEST_REAL_EQUAL(identifications[0].getSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].getSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].getSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].getHits().size(), 2)
	TEST_EQUAL(identifications[1].getHits().size(), 1)
	TEST_EQUAL(identifications[2].getHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].getHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].getHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].getHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].getHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].getHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].getHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].getHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].getHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getHits().size(), 2)
	TEST_EQUAL(protein_identifications[0].getHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].getHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].getHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].getHits()[1].getSequence(), "MYSTVGPA")
									
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
