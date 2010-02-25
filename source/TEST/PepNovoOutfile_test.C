// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

PepNovoOutfile* ptr = 0;
START_SECTION(PepNovoOutfile())
	ptr = new PepNovoOutfile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PepNovoOutfile())
	delete ptr;
END_SECTION

START_SECTION((PepNovoOutfile& operator=(const PepNovoOutfile &pepnovo_outfile)))
	PepNovoOutfile pepnovo_outfile1;
  PepNovoOutfile pepnovo_outfile2;
	pepnovo_outfile2 = pepnovo_outfile1;
	PepNovoOutfile pepnovo_outfile3;
	pepnovo_outfile1 = PepNovoOutfile();
	TEST_EQUAL(( pepnovo_outfile2 == pepnovo_outfile3 ), true)
END_SECTION

START_SECTION((PepNovoOutfile(const PepNovoOutfile &pepnovo_outfile)))
	PepNovoOutfile pepnovo_outfile1;
	PepNovoOutfile pepnovo_outfile2(pepnovo_outfile1);
	PepNovoOutfile pepnovo_outfile3;
	pepnovo_outfile1 = PepNovoOutfile();
	TEST_EQUAL(( pepnovo_outfile2 == pepnovo_outfile3 ), true)
END_SECTION

START_SECTION((bool operator==(const PepNovoOutfile &pepnovo_outfile) const))
	PepNovoOutfile pepnovo_outfile1;
	PepNovoOutfile pepnovo_outfile2;
	TEST_EQUAL(( pepnovo_outfile1 == pepnovo_outfile2 ), true)
END_SECTION

PepNovoOutfile file;


START_SECTION((void load(const std::string &result_filename, std::vector< PeptideIdentification > &peptide_identifications, ProteinIdentification &protein_identification, const DoubleReal &score_threshold, const std::map< String, std::pair< DoubleReal, DoubleReal > > &dta_filenames_and_precursor_retention_times, const std::map< String, String > &keys_to_id)))	
	std::vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;
	map< String, DoubleReal > filenames_and_precursor_retention_times;

	// test exceptions
	//TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.load("a", peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), "the file 'a' could not be found")
	
	//TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out1"), peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), OPENMS_GET_TEST_DATA_PATH_MESSAGE("", "PepNovoOutfile.out1", " in: Not enough columns in file in line 2 (should be 8)!"))
	
	//TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out2"), peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), OPENMS_GET_TEST_DATA_PATH_MESSAGE("", "PepNovoOutfile.out2", " in: Not enough columns in file in line 7 (should be 8)!" ))
	
	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());
	
	
	// test the actual program
	map<String, String>key_to_mod;
	map< String, pair<DoubleReal, DoubleReal> > rt_and_index;

	key_to_mod["K+42"]="Acetyl (K)";
	key_to_mod["Y+42"]="Acetyl (Y)";

	file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), peptide_identifications, protein_identification, 9.000f, rt_and_index, key_to_mod);

	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		//TEST_EQUAL(peptide_identifications[0].getScoreType(), "PepNovo2")
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 9.0)
		if( peptide_identifications[0].getHits().size() == 2)
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 9.888)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "VKEAMAPK")
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 2)
			
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 9.236)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), "KEAMAPK")
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 2)
		}
	}

	file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), peptide_identifications, protein_identification, 0.000f, rt_and_index, key_to_mod);
	
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 20)
		//TEST_EQUAL(peptide_identifications[0].getScoreType(), "PepNovo2")
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 0.0)
		if( peptide_identifications[0].getHits().size() == 20)
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[11].getScore(),8.045)
			TEST_EQUAL(peptide_identifications[0].getHits()[11].getSequence(), "GK(Acetyl)EAMAPK")
			TEST_EQUAL(peptide_identifications[0].getHits()[11].getRank(), 12)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 2)
		}
	}

	END_SECTION


START_SECTION(void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification))
	ProteinIdentification protein_identification;
	
	// test the actual program
	file.getSearchEngineAndVersion(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), protein_identification);
	TEST_EQUAL(protein_identification.getSearchEngine(), "PepNovo+");
	TEST_EQUAL(protein_identification.getSearchEngineVersion(), "Build 20081230");
	TEST_EQUAL(protein_identification.getSearchParameters().peak_mass_tolerance, 0.5);
	TEST_REAL_SIMILAR(protein_identification.getSearchParameters().peak_mass_tolerance, 0.5);
	TEST_REAL_SIMILAR(protein_identification.getSearchParameters().precursor_tolerance, 2.5);
	TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications.size(), 2);
	if(protein_identification.getSearchParameters().variable_modifications.size()== 2)
	{
	  TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications[0], "K+42");
	  TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications[1], "Y+42");
	}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
