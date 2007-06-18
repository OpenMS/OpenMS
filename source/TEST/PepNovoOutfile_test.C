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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/PepNovoOutfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PepNovoOutfile, "$Id: PepNovoInfile_test.C 2076 2007-05-25 16:20:21Z martinlangwisch $")

/////////////////////////////////////////////////////////////

PepNovoOutfile* ptr = 0;
CHECK(PepNovoOutfile())
	ptr = new PepNovoOutfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PepNovoOutfile())
	delete ptr;
RESULT

PepNovoOutfile file;


CHECK(void load(const std::string& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const DoubleReal& p_value_threshold, const std::map< String, Real >& filenames_and_precursor_retention_times) throw (Exception::FileNotFound, Exception::ParseError))
	
	std::vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;
	map< String, Real > filenames_and_precursor_retention_times;
	
	file.load("data/PepNovoOutfile.out", peptide_identifications, protein_identification, 0.915, filenames_and_precursor_retention_times);
	
	TEST_EQUAL(peptide_identifications.size(), 2)
	if ( peptide_identifications.size() == 2 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().empty(), true)
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "PepNovo2")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.915)
		
		TEST_EQUAL(peptide_identifications[1].getHits().size(), 2)
		TEST_EQUAL(peptide_identifications[1].getScoreType(), "PepNovo2")
		TEST_REAL_EQUAL(peptide_identifications[1].getSignificanceThreshold(), 0.915)
		if( peptide_identifications[1].getHits().size() == 2)
		{
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[0].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "CKPTDN")
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getCharge(), 2)
			
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[1].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getSequence(), "KCPTDN")
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getCharge(), 2)
		}
	}
	
	file.load("data/PepNovoOutfile.out", peptide_identifications, protein_identification, 0.94, filenames_and_precursor_retention_times);
	
	TEST_EQUAL(peptide_identifications.size(), 2)
	if ( peptide_identifications.size() == 2 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().empty(), true)
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "PepNovo2")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.94)
		
		TEST_EQUAL(peptide_identifications[1].getHits().size(), 4)
		TEST_EQUAL(peptide_identifications[1].getScoreType(), "PepNovo2")
		TEST_REAL_EQUAL(peptide_identifications[1].getSignificanceThreshold(), 0.94)
		if( peptide_identifications[1].getHits().size() == 4)
		{
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[0].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "CKPTDN")
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getCharge(), 2)
			
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[1].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getSequence(), "KCPTDN")
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getCharge(), 2)
			
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[2].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getSequence(), "KCPTND")
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getRank(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getCharge(), 2)
			
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[3].getScore(), -12.997)
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getSequence(), "CKPTND")
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getRank(), 4)
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getCharge(), 2)
		}
	}
RESULT

CHECK(void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification) throw (Exception::FileNotFound))
	ProteinIdentification protein_identification;
	file.getSearchEngineAndVersion("data/PepNovoOutfile_version_file.txt", protein_identification);
	TEST_EQUAL(protein_identification.getSearchEngine(), "PepNovo");
	TEST_EQUAL(protein_identification.getSearchEngineVersion(), "v2.00");
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
