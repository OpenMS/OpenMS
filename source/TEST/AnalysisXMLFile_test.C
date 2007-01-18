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
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

AnalysisXMLFile* ptr;
AnalysisXMLFile xml_file;

CHECK((AnalysisXMLFile()))
	ptr = new AnalysisXMLFile();
RESULT

CHECK((void load(const String& filename, std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data) const throw(Exception::FileNotFound, Exception::ParseError)))

	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications;

	xml_file.load("data/AnalysisXMLFile_test.analysisXML",
							protein_identifications, 
				   		identifications);
	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(identifications[0].rt, 120)
	TEST_EQUAL(identifications[1].rt, 150)
	TEST_EQUAL(identifications[2].rt, 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(identifications[0].mz, 789.83)
	TEST_REAL_EQUAL(identifications[1].mz, 135.29)
	TEST_REAL_EQUAL(identifications[2].mz, 982.58)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)
	TEST_EQUAL(identifications[2].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccessionType(), "SwissProt")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccessionType(), "SwissProt")
	
RESULT

CHECK((void load(const String& filename, std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data, std::map<String, double>& predicted_retention_times, DoubleReal& predicted_sigma) const throw(Exception::FileNotFound, Exception::ParseError)))

	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications; 
	map<String, double> predicted_retention_times;
	DoubleReal predicted_sigma = 0.0;

	xml_file.load("data/AnalysisXMLFile_test.analysisXML",
							protein_identifications, 
				   		identifications,
							predicted_retention_times,
							predicted_sigma);
	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(identifications[0].rt, 120)
	TEST_EQUAL(identifications[1].rt, 150)
	TEST_EQUAL(identifications[2].rt, 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(identifications[0].mz, 789.83)
	TEST_REAL_EQUAL(identifications[1].mz, 135.29)
	TEST_REAL_EQUAL(identifications[2].mz, 982.58)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)
	TEST_EQUAL(identifications[2].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccessionType(), "SwissProt")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccessionType(), "SwissProt")
	TEST_REAL_EQUAL(predicted_sigma, 0.0852201)
	TEST_EQUAL(predicted_retention_times.size(), 5)
	TEST_REAL_EQUAL(predicted_retention_times[string("LHASGITVTEIPVTATNFK")], 122.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("MRSLGYVAVISAVATDTDK")], 122.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("HSKLSAK")], 151.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("RASNSPQDPQSATAHSFR")], 159.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("MYSTVGPA")], 159.5)


RESULT

CHECK((void store(String filename, const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data) const throw(Exception::UnableToCreateFile)))
												
	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications; 

	String temp_filename = "data/AnalysisXMLFile_test_2.analysisXML";

	NEW_TMP_FILE(temp_filename)
	xml_file.load("data/AnalysisXMLFile_test.analysisXML", 
							protein_identifications, 
				   		identifications);
	xml_file.store(temp_filename, 
							    protein_identifications, 
				   				identifications);
	xml_file.load(temp_filename, 
							protein_identifications, 
				   		identifications);
	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(identifications[0].rt, 120)
	TEST_EQUAL(identifications[1].rt, 150)
	TEST_EQUAL(identifications[2].rt, 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(identifications[0].mz, 789.83)
	TEST_REAL_EQUAL(identifications[1].mz, 135.29)
	TEST_REAL_EQUAL(identifications[2].mz, 982.58)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)
	TEST_EQUAL(identifications[2].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccessionType(), "SwissProt")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccessionType(), "SwissProt")
									
RESULT

CHECK((void store(String filename, const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data, const std::map<String, double>& predicted_retention_times, DoubleReal predicted_sigma) const throw(Exception::UnableToCreateFile)))
												
	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications; 
	map<String, double> predicted_retention_times;
	DoubleReal predicted_sigma = 0.0;

	String temp_filename = "data/AnalysisXMLFile_test_2.analysisXML";
	NEW_TMP_FILE(temp_filename)
	
	xml_file.load("data/AnalysisXMLFile_test.analysisXML", 
							protein_identifications, 
				   		identifications,
							predicted_retention_times,
							predicted_sigma);
	xml_file.store(temp_filename, 
									protein_identifications, 
				   				identifications,
									predicted_retention_times,
									predicted_sigma);

	xml_file.load(temp_filename,
							protein_identifications, 
				   		identifications,
							predicted_retention_times,
							predicted_sigma);

	TEST_EQUAL(identifications.size(), 3)
	TEST_EQUAL(identifications[0].rt, 120)
	TEST_EQUAL(identifications[1].rt, 150)
	TEST_EQUAL(identifications[2].rt, 160)
	PRECISION(0.0001)
	TEST_REAL_EQUAL(identifications[0].mz, 789.83)
	TEST_REAL_EQUAL(identifications[1].mz, 135.29)
	TEST_REAL_EQUAL(identifications[2].mz, 982.58)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideSignificanceThreshold(), 31.8621)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideSignificanceThreshold(), 12)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideSignificanceThreshold(), 19)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)
	TEST_EQUAL(identifications[2].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 33.85)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 33.12)
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 43.9)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[0].getScore(), 5.41)
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[1].getScore(), 7.87)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "HSKLSAK")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getSequence(), "MYSTVGPA")
	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getProteinHits().size(), 2)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccession(), "AAN17824")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccession(), "GN1736")
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[0].getScore(), 103.55)
	TEST_REAL_EQUAL(protein_identifications[0].getProteinHits()[1].getScore(), 67.85)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getScoreType(), "Mascot")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[0].getAccessionType(), "SwissProt")
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccessionType(), "SwissProt")
	TEST_REAL_EQUAL(predicted_sigma, 0.0852201)
	TEST_EQUAL(predicted_retention_times.size(), 5)
	TEST_REAL_EQUAL(predicted_retention_times[string("LHASGITVTEIPVTATNFK")], 122.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("MRSLGYVAVISAVATDTDK")], 122.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("HSKLSAK")], 151.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("RASNSPQDPQSATAHSFR")], 159.5)
	TEST_REAL_EQUAL(predicted_retention_times[string("MYSTVGPA")], 159.5)
									
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
