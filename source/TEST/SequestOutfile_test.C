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

CHECK(void load(const string& result_filename, vector< IdentificationData >&	identifications, ProteinIdentification&	protein_identification, const DoubleReal& p_value_threshold, const vector< Real >& pvalues, const string& database = "", const string& snd_database = "") throw (Exception::FileNotFound, Exception::ParseError))
	vector< IdentificationData > identifications;
	ProteinIdentification protein_identification;
	vector< Real > pvalues;
	
	file.load("data/SequestOutfile2.out", identifications, protein_identification, 1.0, pvalues);
	TEST_EQUAL(identifications.empty(), 1)
	
	file.load("data/SequestOutfile.out", identifications, protein_identification, 1.0, pvalues);
	
	TEST_EQUAL(identifications.size(), 1)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 4)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 0.05)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "ETQAWSIATILETLYDL")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 0.04)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "QVLNPLLVLIFIPLFDL")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getRank(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[2].getScore(), 0.02)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[2].getSequence(), "WVELGPSVLAGVGVM*VLLI")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[2].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[2].getRank(), 3)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[3].getScore(), 0.14)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[3].getSequence(), "FDEITAMTGDGVNDAPALK")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[3].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[3].getRank(), 4)

	identifications.clear();
	pvalues.push_back(0.001);
	pvalues.push_back(0.01);
	pvalues.push_back(0.05);
	pvalues.push_back(0.5);
	file.load("data/SequestOutfile.out", identifications, protein_identification, 0.01, pvalues);
	
	TEST_EQUAL(identifications.size(), 1)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 0.05)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "ETQAWSIATILETLYDL")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), 0.04)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "QVLNPLLVLIFIPLFDL")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "SEQUEST")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getRank(), 2)
RESULT

CHECK(bool getColumns(const String& line, vector< String >& substrings, UnsignedInt number_of_columns, UnsignedInt reference_column))
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

CHECK(void getSequences(const String& database_filename, const map< String, UnsignedInt >& ac_position_map, vector< String >& sequences, vector< pair< String, UnsignedInt > >& found, map< String, UnsignedInt >& not_found) throw (Exception::FileNotFound))
	map< String, UnsignedInt > ac_position_map, not_found;
	ac_position_map["P02666"] = 0;
	ac_position_map["Q9CQV8"] = 1;
	ac_position_map["Q5EEQ7"] = 2;
	ac_position_map["P68509"] = 3;
	vector< String > sequences, found_sequences;
	sequences.push_back("MKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV");
	sequences.push_back("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
	vector< pair< String, UnsignedInt > > found;
	
	file.getSequences("data/Sequest_test.fasta", ac_position_map, found_sequences, found, not_found);
	TEST_EQUAL(found.size(), 2)
	TEST_EQUAL(not_found.size(), 2)
	if ( found.size() == 2 && not_found.size() == 2 )
	{
		map< String, UnsignedInt >::const_iterator i = ++ac_position_map.begin();
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

CHECK(void getACAndACType(String line, string& accession, string& accession_type))
	string accession, accession_type;
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

CHECK(bool updatePeptideHits(PeptideHit& peptide_hit, vector< PeptideHit >& peptide_hits))
	DateTime date;
	date.now();
	vector< PeptideHit > peptide_hits;
	PeptideHit peptide_hit1, peptide_hit2, peptide_hit3;
	peptide_hit1.setScore(12.1);
	peptide_hit1.setScoreType("SEQUEST prelim.");
	peptide_hit1.setSequence("SEQVEST");
	peptide_hit1.setRank(2);
	peptide_hit1.addProteinIndex(date, "P02666");
	TEST_EQUAL(file.updatePeptideHits(peptide_hit1, peptide_hits), true)
	TEST_EQUAL(peptide_hits.size(), 1)
	if ( 1 == peptide_hits.size() )
	{
		TEST_EQUAL(peptide_hit1.getScore(), peptide_hits.front().getScore())
		TEST_EQUAL(peptide_hit1.getScoreType(), peptide_hits.front().getScoreType())
		TEST_EQUAL(peptide_hit1.getSequence(), peptide_hits.front().getSequence())
		TEST_EQUAL(peptide_hit1.getRank(), peptide_hits.front().getRank())

		// the peptides within a vector have to have the same scoreType
		peptide_hit2 = peptide_hit1;
		peptide_hit2.setScoreType("SEQUEST");
		TEST_EQUAL(file.updatePeptideHits(peptide_hit2, peptide_hits), false)
		TEST_EQUAL(peptide_hits.size(), 1)
		
		TEST_EQUAL(peptide_hit1.getScore(), peptide_hits.front().getScore())
		TEST_EQUAL(peptide_hit1.getScoreType(), peptide_hits.front().getScoreType())
		TEST_EQUAL(peptide_hit1.getSequence(), peptide_hits.front().getSequence())
		TEST_EQUAL(peptide_hit1.getRank(), peptide_hits.front().getRank())

		// if the identical peptide hit has already been inserted, nothing changes
		peptide_hit2.setScoreType("SEQUEST prelim.");
		TEST_EQUAL(file.updatePeptideHits(peptide_hit2, peptide_hits), true)
		TEST_EQUAL(peptide_hits.size(), 1)
			
		TEST_EQUAL(peptide_hit1.getScore(), peptide_hits.front().getScore())
		TEST_EQUAL(peptide_hit1.getScoreType(), peptide_hits.front().getScoreType())
		TEST_EQUAL(peptide_hit1.getSequence(), peptide_hits.front().getSequence())
		TEST_EQUAL(peptide_hit1.getRank(), peptide_hits.front().getRank())

		// two peptide hits are considered not equal if either their sequence or their score differs
		peptide_hit2.setSequence("SEQVESTT");
		TEST_EQUAL(file.updatePeptideHits(peptide_hit2, peptide_hits), true)
		TEST_EQUAL(peptide_hits.size(), 2)
		
		TEST_EQUAL(peptide_hit2.getScore(), peptide_hits.back().getScore())
		TEST_EQUAL(peptide_hit2.getScoreType(), peptide_hits.back().getScoreType())
		TEST_EQUAL(peptide_hit2.getSequence(), peptide_hits.back().getSequence())
		TEST_EQUAL(peptide_hit2.getRank(), peptide_hits.back().getRank())

		peptide_hit2.setScore(23.42);
		TEST_EQUAL(file.updatePeptideHits(peptide_hit2, peptide_hits), true)
		TEST_EQUAL(peptide_hits.size(), 3)
		
		TEST_EQUAL(peptide_hit2.getScore(), peptide_hits.back().getScore())
		TEST_EQUAL(peptide_hit2.getScoreType(), peptide_hits.back().getScoreType())
		TEST_EQUAL(peptide_hit2.getSequence(), peptide_hits.back().getSequence())
		TEST_EQUAL(peptide_hit2.getRank(), peptide_hits.back().getRank())

		// if the peptide hit has already been inserted, add additional protein hits
		peptide_hit2.addProteinIndex(date, "Q5EEQ7");
		TEST_EQUAL(file.updatePeptideHits(peptide_hit2, peptide_hits), true)
		TEST_EQUAL(peptide_hits.size(), 3)
		
		TEST_EQUAL(peptide_hit2.getScore(), peptide_hits.back().getScore())
		TEST_EQUAL(peptide_hit2.getScoreType(), peptide_hits.back().getScoreType())
		TEST_EQUAL(peptide_hit2.getSequence(), peptide_hits.back().getSequence())
		TEST_EQUAL(peptide_hit2.getRank(), peptide_hits.back().getRank())

		// but the peptide hit (outside the vector) remains unchanged
		peptide_hit3 = peptide_hit2;
		peptide_hit3.getProteinIndices().clear();
		peptide_hit3.addProteinIndex(date, "Q6UN63");
		TEST_EQUAL(file.updatePeptideHits(peptide_hit3, peptide_hits), true)
		TEST_EQUAL(peptide_hits.size(), 3)
		peptide_hit2.addProteinIndex(date, "Q6UN63");
		
		TEST_EQUAL(peptide_hit2.getScore(), peptide_hits.back().getScore())
		TEST_EQUAL(peptide_hit2.getScoreType(), peptide_hits.back().getScoreType())
		TEST_EQUAL(peptide_hit2.getSequence(), peptide_hits.back().getSequence())
		TEST_EQUAL(peptide_hit2.getRank(), peptide_hits.back().getRank())
		
		TEST_EQUAL(peptide_hit3.getProteinIndices().size(), 1)
	}
RESULT


END_TEST
