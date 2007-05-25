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

CHECK(vector< UInt > load(const String& result_filename, vector<PeptideIdentification>&	peptide_identifications, ProteinIdentification&	protein_identification, const DoubleReal& p_value_threshold, const String& database_filename) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument))
	vector<PeptideIdentification> peptide_identifications;
	ProteinIdentification protein_identification;
	
	file.load("data/InspectOutfile.out", peptide_identifications, protein_identification, 0.001);
	
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 1)
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.001)
		if( peptide_identifications[0].getHits().size() == 1)
		{
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), -257)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "KKLE")
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAABefore(), 'E')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 1)
			if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
			{
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
			}
		}
	}

	peptide_identifications.clear();
	file.load("data/InspectOutfile.out", peptide_identifications, protein_identification, 0.01);
	
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		TEST_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.01)
		if( peptide_identifications[0].getHits().size() == 2 )
		{
	// 		if( peptide_identifications[0].getHits().size() == 1)
	// 		{
				TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), -257)
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "KKLE")
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getAABefore(), 'E')
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getAAAfter(), 'K')
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
				TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 1)
				if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
				{
					TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
				}
				TEST_REAL_EQUAL(peptide_identifications[0].getHits()[1].getScore(), -1456)
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), "EKIE")
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getAABefore(), 'R')
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getAAAfter(), 'K')
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 1)
				TEST_EQUAL(peptide_identifications[0].getHits()[1].getProteinAccessions().size(), 1)
				if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
				{
					TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
				}
	// 		}
		}
	}
RESULT

CHECK(void generateTrieDB(const std::String& source_database_filename, const std::String& database_filename, const std::String& index_filename, bool append = false, const std::String species = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile))
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

CHECK(void compressTrieDB(const std::String& database_filename, const std::String& index_filename, const std::vector< UInt >& wanted_records, const std::String& snd_database_filename, const std::String& snd_index_filename, bool append = false) throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile))
	vector< UInt > wanted_records(1, 0);
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


CHECK(void getSequences(const String& database_filename, const map< Unsigned, UInt >& getMetaValue("RT")_position_map, vector< String >& sequences) throw (Exception::FileNotFound, Exception::ParseError))
	map< UInt, UInt > rn_position_map;
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

CHECK(void getPrecursorRTandMZ(const vector< pair< String, vector< UInt > > >& files_and_scan_numbers, std::vector< PeptideIdentification >& ids))
	vector< pair< String, vector< UInt > > > files_and_scan_numbers;
	files_and_scan_numbers.push_back(make_pair("data/InspectOutfile_test_1.mzXML", vector< UInt >(1, 4)));
	files_and_scan_numbers.push_back(make_pair("data/InspectOutfile_test_2.mzXML", vector< UInt >(1, 4)));
	
	vector< PeptideIdentification > ids, ids_found;
	ids_found.push_back(PeptideIdentification());
	ids_found.push_back(PeptideIdentification());
	
	ids.push_back(PeptideIdentification());
	ids.back().setMetaValue("RT", -1);
	ids.back().setMetaValue("MZ", 123.456);
	ids.push_back(PeptideIdentification());
	ids.back().setMetaValue("RT", 180);
	ids.back().setMetaValue("MZ", 123.456);
	
	file.getPrecursorRTandMZ(files_and_scan_numbers, ids_found);
	
	TEST_REAL_EQUAL(ids_found.front().getMetaValue("RT"), ids.front().getMetaValue("RT"));
	TEST_REAL_EQUAL(ids_found.front().getMetaValue("MZ"), ids.front().getMetaValue("MZ"));
	TEST_REAL_EQUAL(ids_found.back().getMetaValue("RT"), ids.back().getMetaValue("RT"));
	TEST_REAL_EQUAL(ids_found.back().getMetaValue("MZ"), ids.back().getMetaValue("MZ"));
RESULT

CHECK(void getLabels(const String& source_database_filename, String& ac_label, String& sequence_start_label, String& sequence_end_label, String& comment_label, String& species_label) throw (Exception::FileNotFound, Exception::ParseError))
	String ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
	file.getLabels("data/Inspect_test.fasta", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
	TEST_EQUAL(ac_label, ">")
	TEST_EQUAL(sequence_start_label, ">")
	TEST_EQUAL(sequence_end_label, ">")
	TEST_EQUAL(comment_label, ";")
	TEST_EQUAL(species_label, ">")
RESULT

CHECK(vector< UInt > getWantedRecords(const String& result_filename, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument))
	vector< UInt > wanted_records = file.getWantedRecords("data/InspectOutfile.out", 0.01);
	TEST_EQUAL (wanted_records.size(), 1)
	if ( !wanted_records.empty() ) TEST_EQUAL (wanted_records.front(), 0)
RESULT

CHECK(template<typename PeakT> void getExperiment(MSExperiment<PeakT>& exp, String& type, const String& in_filename) throw(Exception::ParseError))
	MSExperiment< > exp;
	String type;
	file.getExperiment(exp, type, "TOPP/Inspect_Sequest.mzXML");
	TEST_EQUAL(type, "mzXML")
	file.getExperiment(exp, type, "TOPP/Inspect_Sequest.mzData");
	TEST_EQUAL(type, "mzData")
RESULT

CHECK(void getSearchEngineAndVersion(const String& inspect_output_without_parameters_filename, ProteinIdentification& protein_identification) throw (Exception::FileNotFound))
	ProteinIdentification protein_identification;
	file.getSearchEngineAndVersion("data/InspectOutfile_version_file.txt", protein_identification);
	TEST_EQUAL(protein_identification.getSearchEngine(), "InsPecT");
	TEST_EQUAL(protein_identification.getSearchEngineVersion(), "20060907");
RESULT

CHECK(void readOutHeader(const String& filename, const String& header_line, Int& spectrum_file_column, Int& scan_column, Int& peptide_column, Int& protein_column, Int& charge_column, Int& MQ_score_column, Int& p_value_column, Int& record_number_column, Int& DB_file_pos_column, Int& spec_file_pos_column) throw (Exception::ParseError))
		
	String header_line = "#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos";
	
	Int spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column;
	UInt number_of_columns;
	
	file.readOutHeader("dummy_testfile", header_line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns);
	TEST_EQUAL(spectrum_file_column, 0)
	TEST_EQUAL(scan_column, 1)
	TEST_EQUAL(peptide_column, 2)
	TEST_EQUAL(protein_column, 3)
	TEST_EQUAL(charge_column, 4)
	TEST_EQUAL(MQ_score_column, 5)
	TEST_EQUAL(p_value_column, 13)
	TEST_EQUAL(record_number_column, 17)
	TEST_EQUAL(DB_file_pos_column, 18)
	TEST_EQUAL(spec_file_pos_column, 19)
	TEST_EQUAL(number_of_columns, 20)
	
	header_line = "#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	";
	
	TEST_EXCEPTION(Exception::ParseError, file.readOutHeader("dummy_testfile", header_line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns));
RESULT

END_TEST
