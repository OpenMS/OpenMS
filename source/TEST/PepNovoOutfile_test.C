// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
PepNovoOutfile* nullPointer = 0;
START_SECTION(PepNovoOutfile())
	ptr = new PepNovoOutfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
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


START_SECTION((void load(const std::string &result_filename, std::vector< PeptideIdentification > &peptide_identifications, ProteinIdentification &protein_identification, const DoubleReal &score_threshold, const std::map< String, std::pair< DoubleReal, DoubleReal > > &id_rt_mz, const std::map< String, String > &mod_id_map)))	
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
