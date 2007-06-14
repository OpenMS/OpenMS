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

#include <OpenMS/FORMAT/SequestOutfile.h>

using namespace OpenMS;
using namespace std;


START_TEST(String, "$Id$")


SequestOutfile* ptr = 0;

CHECK(SequestOutfile())
	ptr = new SequestOutfile();
	TEST_NOT_EQUAL(ptr, 0)
	delete(ptr);
RESULT

SequestOutfile file;

CHECK(void load(const String& result_filename, vector< PeptideIdentification >&	peptide_identifications, ProteinIdentification&	protein_identification, const Real& p_value_threshold, const vector< Real >& pvalues, const String& database = "", const String& snd_database = "") throw (Exception::FileNotFound, Exception::ParseError))
	vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;
	vector< Real > pvalues;
	
	file.load("data/SequestOutfile2.out", peptide_identifications, protein_identification, 1.0, pvalues);
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 ) TEST_EQUAL(peptide_identifications[0].getHits().empty(), true)
	
	file.load("data/SequestOutfile.out", peptide_identifications, protein_identification, 1.0, pvalues);
	
	TEST_EQUAL(peptide_identifications.size(), 2)
	if ( peptide_identifications.size() == 2 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().empty(), true)
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "SEQUEST prelim.")
		TEST_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 1.0)
		TEST_EQUAL(peptide_identifications[1].getHits().size(), 4)
		TEST_EQUAL(peptide_identifications[1].getScoreType(), "SEQUEST")
		TEST_EQUAL(peptide_identifications[1].getIdentifier(), "TurboSEQUEST_2004-03-16")
		TEST_REAL_EQUAL(peptide_identifications[1].getSignificanceThreshold(), 1.0)
		if ( peptide_identifications[1].getHits().size() == 4 )
		{
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[0].getScore(), 0.05)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "ETQAWSIATILETLYDL")
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getAABefore(), 'C')
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getAAAfter(), '-')
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[0].getProteinAccessions().size(), 3)
			if ( peptide_identifications[1].getHits()[0].getProteinAccessions().size() == 3 )
			{
				TEST_EQUAL(peptide_identifications[1].getHits()[0].getProteinAccessions()[0], "P35574")
				TEST_EQUAL(peptide_identifications[1].getHits()[0].getProteinAccessions()[1], "2136928")
				TEST_EQUAL(peptide_identifications[1].getHits()[0].getProteinAccessions()[2], "L10605")
			}
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[1].getScore(), 0.04)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getSequence(), "QVLNPLLVLIFIPLFDL")
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getAABefore(), 'M')
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getAAAfter(), 'V')
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[1].getProteinAccessions().size(), 2)
			if ( peptide_identifications[1].getHits()[1].getProteinAccessions().size() == 2 )
			{
				TEST_EQUAL(peptide_identifications[1].getHits()[1].getProteinAccessions()[0], "P46029")
				TEST_EQUAL(peptide_identifications[1].getHits()[1].getProteinAccessions()[1], "U32507")
			}
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[2].getScore(), 0.02)
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getSequence(), "WVELGPSVLAGVGVMVLLI")
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getAABefore(), 'L')
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getAAAfter(), 'P')
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getRank(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[2].getProteinAccessions().size(), 1)
			if ( peptide_identifications[1].getHits()[2].getProteinAccessions().size() == 1 )
			{
				TEST_EQUAL(peptide_identifications[1].getHits()[2].getProteinAccessions()[0], "e148876")
			}
			TEST_REAL_EQUAL(peptide_identifications[1].getHits()[3].getScore(), 0.14)
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getSequence(), "FDEITAMTGDGVNDAPALK")
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getAABefore(), 'S')
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getRank(), 4)
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions().size(), 8)
			if ( peptide_identifications[1].getHits()[3].getProteinAccessions().size() == 8 )
			{
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[0], "P20647")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[1], "P04192")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[2], "67962")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[3], "67961")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[4], "109166")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[5], "224621")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[6], "X02814")
				TEST_EQUAL(peptide_identifications[1].getHits()[3].getProteinAccessions()[7], "J04703")
			}
		}
	}

	peptide_identifications.clear();
	pvalues.push_back(0.001);
	pvalues.push_back(0.01);
	pvalues.push_back(0.05);
	pvalues.push_back(0.5);
	file.load("data/SequestOutfile.out", peptide_identifications, protein_identification, 0.01, pvalues);
	
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "SEQUEST")
		TEST_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.01)
		if ( peptide_identifications[0].getHits().size() == 2 )
		{
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), 0.05)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "ETQAWSIATILETLYDL")
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAABefore(), 'C')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAAAfter(), '-')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 3)
			if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 3 )
			{
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P35574")
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[1], "2136928")
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[2], "L10605")
			}
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[1].getScore(), 0.04)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), "QVLNPLLVLIFIPLFDL")
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getAABefore(), 'M')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getAAAfter(), 'V')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 3)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getProteinAccessions().size(), 2)
			if ( peptide_identifications[0].getHits()[1].getProteinAccessions().size() == 2 )
			{
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getProteinAccessions()[0], "P46029")
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getProteinAccessions()[1], "U32507")
			}
		}
		TEST_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
	}
	
	TEST_EQUAL(protein_identification.getSearchEngine(), "TurboSEQUEST")
	TEST_EQUAL(protein_identification.getSearchEngineVersion(), "v.27 (rev. 12)")
	TEST_EQUAL(protein_identification.getIdentifier(), "TurboSEQUEST_2004-03-16")
RESULT

CHECK(bool getColumns(const String& line, vector< String >& substrings, UInt number_of_columns, UInt reference_column))
	String line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT   +2   C.ETQAWSIATILETLYDL.-";
	vector< String > substrings, columns;
	columns.push_back("1.");
	columns.push_back("1/80");
	columns.push_back("0");
	columns.push_back("1967.0013");
	columns.push_back("0.0000");
	columns.push_back("1.5789");
	columns.push_back("310.3");
	columns.push_back("0.05");
	columns.push_back("0");
	columns.push_back("18/64");
	columns.push_back("gi|544379|sp|P35574|GDE RABIT+2");
	columns.push_back("C.ETQAWSIATILETLYDL.-");
	TEST_EQUAL(file.getColumns("", substrings, 12, 10), false)
	TEST_EQUAL(file.getColumns(line, substrings, 12, 10), true)
	TEST_EQUAL((columns == substrings), true)
	
	line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT+2   C.ETQAWSIATILETLYDL.-";
	TEST_EQUAL(file.getColumns(line, substrings, 12, 10), true)
	TEST_EQUAL((columns == substrings), true)
	
	line = "  1.   1/80          0 1967.0013  0.0000  1.5789   310.3 0.05    0 18/64 gi|544379|sp|P35574|GDE RABIT   +X   C.ETQAWSIATILETLYDL.-";
	TEST_EQUAL(file.getColumns(line, substrings, 12, 10), true)
	columns[10] = "gi|544379|sp|P35574|GDE RABIT +X";
	TEST_EQUAL((columns == substrings), true)
RESULT

CHECK(void getSequences(const String& database_filename, const map< String, UInt >& ac_position_map, vector< String >& sequences, vector< pair< String, UInt > >& found, map< String, UInt >& not_found) throw (Exception::FileNotFound))
	map< String, UInt > ac_position_map, not_found;
	ac_position_map["P02666"] = 0;
	ac_position_map["Q9CQV8"] = 1;
	ac_position_map["Q5EEQ7"] = 2;
	ac_position_map["P68509"] = 3;
	vector< String > sequences, found_sequences;
	sequences.push_back("MKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV");
	sequences.push_back("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
	vector< pair< String, UInt > > found;
	
	file.getSequences("data/Sequest_test.fasta", ac_position_map, found_sequences, found, not_found);
	TEST_EQUAL(found.size(), 2)
	TEST_EQUAL(not_found.size(), 2)
	if ( found.size() == 2 && not_found.size() == 2 )
	{
		map< String, UInt >::const_iterator i = ++ac_position_map.begin();
		TEST_EQUAL(i->first, found[0].first)
		TEST_EQUAL(i->second, found[0].second)

		i = --ac_position_map.end();
		TEST_EQUAL(i->first, found[1].first)
		TEST_EQUAL(i->second, found[1].second)
		
		TEST_EQUAL(sequences[3], found_sequences[0])
		TEST_EQUAL(sequences[1], found_sequences[1])

		file.getSequences("data/Sequest_test2.fasta", not_found, found_sequences, found, not_found);
		TEST_EQUAL(found.size(), 4)
		TEST_EQUAL(not_found.size(), 0)

		i = ac_position_map.begin();
		TEST_EQUAL(i->first, found[2].first)
		TEST_EQUAL(i->second, found[2].second)

		++i; ++i;
		TEST_EQUAL(i->first, found[3].first)
		TEST_EQUAL(i->second, found[3].second)
		
		TEST_EQUAL(sequences[0], found_sequences[2])
		TEST_EQUAL(sequences[2], found_sequences[3])
	}
	ac_position_map.clear();
	sequences.clear();
	found_sequences.clear();
	found.clear();
	not_found.clear();
RESULT

CHECK(void getACAndACType(String line, String& accession, String& accession_type))
	String accession, accession_type;
	file.getACAndACType(">sp|P02666|CASB_BOVIN Beta-casein precursor - Bos taurus (Bovine).", accession, accession_type);
	TEST_EQUAL(accession, "P02666")
	TEST_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType(">tr|Q5EEQ7|Q5EEQ7_BOVIN Beta-casein (Fragment) - Bos taurus (Bovine).", accession, accession_type);
	TEST_EQUAL(accession, "Q5EEQ7")
	TEST_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType("gi|110174602|gb|DQ660451.1|", accession, accession_type);
	TEST_EQUAL("DQ660451.1", accession)
	TEST_EQUAL("GenBank", accession_type)
	
	file.getACAndACType("gi|1655698|emb|Y07752|VCPHEROPH", accession, accession_type);
	TEST_EQUAL(accession, "Y07752")
	TEST_EQUAL(accession_type, "EMBL")
	
	file.getACAndACType("gi|10038695|dbj|BAB12730|", accession, accession_type);
	TEST_EQUAL(accession, "BAB12730")
	TEST_EQUAL(accession_type, "DDBJ")
	
	file.getACAndACType("gi|9628804|ref|NP_043835|", accession, accession_type);
	TEST_EQUAL(accession, "NP_043835")
	TEST_EQUAL(accession_type, "NCBI")
	
	file.getACAndACType("gi|21362794|sp|P58858.0|", accession, accession_type);
	TEST_EQUAL(accession, "P58858.0")
	TEST_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType("gi|21362794|tr|P58858.0|", accession, accession_type);
	TEST_EQUAL(accession, "P58858.0")
	TEST_EQUAL(accession_type, "SwissProt")

	file.getACAndACType("gi|1619818|gnl|PID|d1013471|", accession, accession_type);
	TEST_EQUAL(accession, "d1013471")
	TEST_EQUAL(accession_type, "PID")

	file.getACAndACType("Q30DX2 Gamma-gliadin/LMW-glutenin chimera Ch7 (Fragment).", accession, accession_type);
	TEST_EQUAL(accession, "Q30DX2")
	TEST_EQUAL(accession_type, "SwissProt")

	file.getACAndACType(">P68509|1433F_BOVIN", accession, accession_type);
	TEST_EQUAL(accession, "P68509")
	TEST_EQUAL(accession_type, "SwissProt")

	file.getACAndACType(">ACBLA (P68509) F_BOVIN", accession, accession_type);
	TEST_EQUAL(accession, "P68509")
	TEST_EQUAL(accession_type, "SwissProt")
RESULT

END_TEST
