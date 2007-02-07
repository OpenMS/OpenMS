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

#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////


InspectOutfile* ptr = 0;
CHECK(InspectOutfile())
	ptr = new InspectOutfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~InspectOutfile())
	delete ptr;
RESULT

InspectOutfile file;

CHECK(vector< UnsignedInt > load(const string& result_filename, vector< IdentificationData >&	identifications, ProteinIdentification&	protein_identification, const DoubleReal& p_value_threshold, const string& database_filename) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument))
	vector< IdentificationData > identifications;
	ProteinIdentification protein_identification;
	file.load("data/InspectOutfile.out", identifications, protein_identification, 0.001);
	
	TEST_EQUAL(identifications.size(), 1)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 1)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), -257)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "KKLE")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Inspect")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)

	identifications.clear();
	file.load("data/InspectOutfile.out", identifications, protein_identification, 0.01);
	
	TEST_EQUAL(identifications.size(), 1)
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 2)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), -1456)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "EKIE")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Inspect")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[1].getScore(), -257)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getSequence(), "KKLE")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getScoreType(), "Inspect")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getRank(), 2)
RESULT

CHECK(void generateTrieDB(const std::string& source_database_filename, const std::string& database_filename, const std::string& index_filename, bool append = false, const std::string species = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile))
	file.generateTrieDB("data/Inspect_test.fasta", "InspectOutfile_test.trie", "InspectOutfile_test.index");
	TEST_FILE("InspectOutfile_test.trie", "data/Inspect_test.trie")
	TEST_FILE("InspectOutfile_test.index", "data/Inspect_test.index")
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
	
	file.generateTrieDB("data/Inspect_test2.fasta", "InspectOutfile_test.trie", "InspectOutfile_test.index");
	TEST_FILE("InspectOutfile_test.trie", "data/Inspect_test2.trie")
	TEST_FILE("InspectOutfile_test.index", "data/Inspect_test2.index")
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
RESULT

CHECK(void compressTrieDB(const std::string& database_filename, const std::string& index_filename, const std::vector< UnsignedInt >& wanted_records, const std::string& snd_database_filename, const std::string& snd_index_filename, bool append = false) throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile))
	vector< UnsignedInt > wanted_records(1, 0);
	file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", true);

	bool append = true;
	wanted_records.clear();
	file.compressTrieDB("data/Inspect_test2.trie", "data/Inspect_test2.index", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);

	wanted_records.push_back(1);
	file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);

	wanted_records.clear();
	wanted_records.push_back(0);
	wanted_records.push_back(3);
	file.compressTrieDB("InspectOutfile_test.trie", "InspectOutfile_test.index", wanted_records, "InspectOutfile_test2.trie", "InspectOutfile_test2.index");

	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");

	TEST_FILE("InspectOutfile_test2.trie", "data/Inspect_test.trie")
	TEST_FILE("InspectOutfile_test2.index", "data/Inspect_test.index")


	remove("InspectOutfile_test2.trie");
	remove("InspectOutfile_test2.index");
RESULT


CHECK(void getSequences(const String& database_filename, const map< Unsigned, UnsignedInt >& rt_position_map, vector< String >& sequences) throw (Exception::FileNotFound, Exception::ParseError))
	map< UnsignedInt, UnsignedInt > rn_position_map;
	rn_position_map[0] = 0;
	rn_position_map[1] = 1;
	vector< String > sequences, found_sequences;
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");

	file.getSequences("data/Inspect_test.trie", rn_position_map, found_sequences);
	if ( found_sequences.size() == 1 )
	{
		TEST_EQUAL(found_sequences.front(), sequences.front())

		file.getSequences("data/Inspect_test2.trie", rn_position_map, found_sequences);

		TEST_EQUAL((sequences == found_sequences), true)
	}
	rn_position_map.clear();
	sequences.clear();
	found_sequences.clear();
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
RESULT

CHECK(void getPrecursorRTandMZ(const vector< pair< String, vector< UnsignedInt > > >& files_and_scan_numbers, std::vector< IdentificationData >& ids))
	vector< pair< String, vector< UnsignedInt > > > files_and_scan_numbers;
	files_and_scan_numbers.push_back(make_pair("data/InspectOutfile_test_1.mzXML", vector< UnsignedInt >(1, 4)));
	files_and_scan_numbers.push_back(make_pair("data/InspectOutfile_test_2.mzXML", vector< UnsignedInt >(1, 4)));
	
	vector< IdentificationData > ids, ids_found;
	ids_found.push_back(IdentificationData());
	ids_found.push_back(IdentificationData());
	
	ids.push_back(IdentificationData());
	ids.back().rt = -1;
	ids.back().mz = 123.456;
	ids.push_back(IdentificationData());
	ids.back().rt = 180;
	ids.back().mz = 123.456;
	
	file.getPrecursorRTandMZ(files_and_scan_numbers, ids_found);
	
	TEST_REAL_EQUAL(ids_found.front().rt, ids.front().rt);
	TEST_REAL_EQUAL(ids_found.front().mz, ids.front().mz);
	TEST_REAL_EQUAL(ids_found.back().rt, ids.back().rt);
	TEST_REAL_EQUAL(ids_found.back().mz, ids.back().mz);
RESULT

CHECK(void getLabels(const string& source_database_filename, string& ac_label, string& sequence_start_label, string& sequence_end_label, string& comment_label, string& species_label) throw (Exception::FileNotFound, Exception::ParseError))
	string ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
	file.getLabels("data/Inspect_test.fasta", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
	TEST_EQUAL(ac_label, ">")
	TEST_EQUAL(sequence_start_label, ">")
	TEST_EQUAL(sequence_end_label, ">")
	TEST_EQUAL(comment_label, ";")
	TEST_EQUAL(species_label, ">")
RESULT

CHECK(vector< UnsignedInt > getWantedRecords(const string& result_filename, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument))
	vector< UnsignedInt > wanted_records = file.getWantedRecords("data/InspectOutfile.out", 0.01);
	TEST_EQUAL (wanted_records.size(), 1)
	if ( !wanted_records.empty() ) TEST_EQUAL (wanted_records.front(), 0)
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
