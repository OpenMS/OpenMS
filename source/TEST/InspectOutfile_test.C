// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
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

CHECK((InspectOutfile& operator=(const InspectOutfile &inspect_outfile)))
	InspectOutfile inspect_outfile1;
  InspectOutfile inspect_outfile2;
	inspect_outfile2 = inspect_outfile1;
	InspectOutfile inspect_outfile3;
	inspect_outfile1 = InspectOutfile();
	TEST_EQUAL( (inspect_outfile2 == inspect_outfile3 ), true)
RESULT

CHECK((InspectOutfile(const InspectOutfile &inspect_outfile)))
	InspectOutfile inspect_outfile1;
	InspectOutfile inspect_outfile2(inspect_outfile1);
	InspectOutfile inspect_outfile3;
	inspect_outfile1 = InspectOutfile();
	TEST_EQUAL( (inspect_outfile2 == inspect_outfile3 ), true)
RESULT

CHECK((bool operator==(const InspectOutfile &inspect_outfile) const))
	InspectOutfile inspect_outfile1;
	InspectOutfile inspect_outfile2;
	TEST_EQUAL(( inspect_outfile1 == inspect_outfile2 ), true)
RESULT

InspectOutfile file;

CHECK(std::vector< UInt > load(const String& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const Real p_value_threshold, const String& database_filename = ""))
	vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;
	
	// test exceptions
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification, 2.0))
	
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification, -1.0))
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.load("a", peptide_identifications, protein_identification, 0.01), "the file 'a' could not be found")
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileEmpty, file.load("data/Inspect_empty_file.txt", peptide_identifications, protein_identification, 0.01), "the file 'data/Inspect_empty_file.txt' is empty")
	
	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());
	
	
	file.load("data/InspectOutfile.out", peptide_identifications, protein_identification, 0.001);
	
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 1)
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.001)
		if( peptide_identifications[0].getHits().size() == 1)
		{
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), -257)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "KKLE")
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAABefore(), 'E')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 1)
			if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
			{
				TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
			}
		}
	}
	peptide_identifications.clear();
std::cout << "MARTIN" << std::endl;
	file.load("data/InspectOutfile.out", peptide_identifications, protein_identification, 0.01);
std::cout << "MARTIN" << std::endl;
	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_EQUAL(peptide_identifications[0].getSignificanceThreshold(), 0.01)
		if( peptide_identifications[0].getHits().size() == 2 )
		{
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[0].getScore(), -257)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "KKLE")
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAABefore(), 'E')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions().size(), 1)
			if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
			{
				TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
			}
			TEST_REAL_EQUAL(peptide_identifications[0].getHits()[1].getScore(), -1456)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getSequence().toString(), "EKIE")
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getAABefore(), 'R')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getProteinAccessions().size(), 1)
			if ( peptide_identifications[0].getHits()[0].getProteinAccessions().size() == 1 )
			{
				TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getProteinAccessions()[0], "P68509")
			}
		}
	}
	
	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());
	TEST_EQUAL(file.load("data/InspectOutfile.out1", peptide_identifications, protein_identification, 0.01).size(), 1)
	if ( file.load("data/InspectOutfile.out1", peptide_identifications, protein_identification, 0.01).size() == 1 )
	{
		TEST_EQUAL(file.load("data/InspectOutfile.out1", peptide_identifications, protein_identification, 0.01)[0], 2)
	}
RESULT

CHECK(void generateTrieDB(const std::String& source_database_filename, const std::String& database_filename, const std::String& index_filename, bool append = false, const std::String species = ""))
	// test exceptions
	// test file not found for input file
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.generateTrieDB("a", "", ""), "the file 'a' could not be found")
	
	// test unable to create file
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.generateTrieDB("data/Inspect_test.fasta", "data/Inspect_unreadable_unwriteable.txt", ""), "the file 'data/Inspect_unreadable_unwriteable.txt' could not be created")
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.generateTrieDB("data/Inspect_test.fasta", "InspectOutfile_test.trie", "data/Inspect_unreadable_unwriteable.txt"), "the file 'data/Inspect_unreadable_unwriteable.txt' could not be created")

	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");


	// test the actual program
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

CHECK(void compressTrieDB(const String& database_filename, const String& index_filename, std::vector< UInt >& wanted_records, const String& snd_database_filename, const String& snd_index_filename, bool append = false))
	vector< UInt > wanted_records(1, 0);

	// test exceptions
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
	remove("InspectOutfile_test2.trie");
	remove("InspectOutfile_test2.index");
	// test for equal filenames
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "data/Inspect_test.trie", ""), "data/Inspect_test.trie in: Same filename can not be used for original and second database!")
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "", "data/Inspect_test.index"), "data/Inspect_test.index in: Same filename can not be used for original and second database!")

	// test file not found for input files (using empty filenames)
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.compressTrieDB("a", "", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index"), "the file 'a' could not be found")
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.compressTrieDB("data/Inspect_test.trie", "b", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index"), "the file 'b' could not be found")

	// test for unable to create file
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "data/Inspect_unreadable_unwriteable.txt", "", true), "the file 'data/Inspect_unreadable_unwriteable.txt' could not be created")
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "InspectOutfile_test.trie", "data/Inspect_unreadable_unwriteable.txt", true), "the file 'data/Inspect_unreadable_unwriteable.txt' could not be created")

	// test for parse error
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_empty_file.txt", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", true), "data/Inspect_empty_file.txt in: index file is too short!")
	
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
	wanted_records.clear();


	// test the actual program
	bool append = true;
	wanted_records.push_back(0);
	file.compressTrieDB("data/Inspect_test.trie", "data/Inspect_test.index", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);
	wanted_records.clear();

	file.compressTrieDB("data/Inspect_test2.trie", "data/Inspect_test2.index", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);

	wanted_records.clear();
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


CHECK(std::vector< UInt > getSequences(const String& database_filename, const std::map< UInt, UInt >& wanted_records, std::vector< String >& sequences))
	map< UInt, UInt > rn_position_map;
	rn_position_map[0] = 0;
	rn_position_map[1] = 1;
	vector< String > sequences, found_sequences;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getSequences("a", rn_position_map, found_sequences), "the file 'a' could not be found")
	
	found_sequences.clear();

	// test actual program
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");
	
	file.getSequences("data/Inspect_test.trie", rn_position_map, found_sequences);
	if ( found_sequences.size() == 1 )
	{
		TEST_STRING_EQUAL(found_sequences.front(), sequences.front())

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
	TEST_STRING_EQUAL(accession, "P02666")
	TEST_STRING_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType(">tr|Q5EEQ7|Q5EEQ7_BOVIN Beta-casein (Fragment) - Bos taurus (Bovine).", accession, accession_type);
	TEST_STRING_EQUAL(accession, "Q5EEQ7")
	TEST_STRING_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType("gi|110174602|gb|DQ660451.1|", accession, accession_type);
	TEST_STRING_EQUAL("DQ660451.1", accession)
	TEST_STRING_EQUAL("GenBank", accession_type)
	
	file.getACAndACType("gi|1655698|emb|Y07752|VCPHEROPH", accession, accession_type);
	TEST_STRING_EQUAL(accession, "Y07752")
	TEST_STRING_EQUAL(accession_type, "EMBL")
	
	file.getACAndACType("gi|10038695|dbj|BAB12730|", accession, accession_type);
	TEST_STRING_EQUAL(accession, "BAB12730")
	TEST_STRING_EQUAL(accession_type, "DDBJ")
	
	file.getACAndACType("gi|9628804|ref|NP_043835|", accession, accession_type);
	TEST_STRING_EQUAL(accession, "NP_043835")
	TEST_STRING_EQUAL(accession_type, "NCBI")
	
	file.getACAndACType("gi|21362794|sp|P58858.0|", accession, accession_type);
	TEST_STRING_EQUAL(accession, "P58858.0")
	TEST_STRING_EQUAL(accession_type, "SwissProt")
	
	file.getACAndACType("gi|21362794|tr|P58858.0|", accession, accession_type);
	TEST_STRING_EQUAL(accession, "P58858.0")
	TEST_STRING_EQUAL(accession_type, "SwissProt")

	file.getACAndACType("gi|1619818|gnl|PID|d1013471|", accession, accession_type);
	TEST_STRING_EQUAL(accession, "d1013471")
	TEST_STRING_EQUAL(accession_type, "PID")

	file.getACAndACType("Q30DX2 Gamma-gliadin/LMW-glutenin chimera Ch7 (Fragment).", accession, accession_type);
	TEST_STRING_EQUAL(accession, "Q30DX2")
	TEST_STRING_EQUAL(accession_type, "SwissProt")

	file.getACAndACType(">P68509|1433F_BOVIN", accession, accession_type);
	TEST_STRING_EQUAL(accession, "P68509")
	TEST_STRING_EQUAL(accession_type, "SwissProt")

	file.getACAndACType(">ACBLA (P68509) F_BOVIN", accession, accession_type);
	TEST_STRING_EQUAL(accession, "P68509")
	TEST_STRING_EQUAL(accession_type, "SwissProt")
RESULT

CHECK(void getPrecursorRTandMZ(const vector< pair< String, vector< pair < UInt, UInt > > > >& files_and_peptide_identification_with_scan_number, std::vector< PeptideIdentification >& ids))
	vector< pair< String, vector< pair< UInt, UInt > > > > files_and_peptide_identification_with_scan_number;
	vector< PeptideIdentification > ids, ids_found;

	// test exceptions
	files_and_peptide_identification_with_scan_number.push_back(make_pair("data/InspectOutfile_test_1.mzXML", vector< pair< UInt, UInt > >(1, make_pair(0, 10))));
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.getPrecursorRTandMZ(files_and_peptide_identification_with_scan_number, ids_found), "data/InspectOutfile_test_1.mzXML in: Not enought scans in file! (4 available, should be at least 10)")
	
	files_and_peptide_identification_with_scan_number.clear();
	ids.clear();
	ids_found.clear();
	
	files_and_peptide_identification_with_scan_number.push_back(make_pair("data/InspectOutfile_test_1.mzXML", vector< pair < UInt, UInt > >(1, make_pair(0, 4))));
	files_and_peptide_identification_with_scan_number.push_back(make_pair("data/InspectOutfile_test_2.mzXML", vector< pair < UInt, UInt > >(1, make_pair(1, 4))));
	ids_found.push_back(PeptideIdentification());
	ids_found.push_back(PeptideIdentification());
	
	ids.push_back(PeptideIdentification());
	ids.back().setMetaValue("RT", -1);
	ids.back().setMetaValue("MZ", 123.456);
	ids.push_back(PeptideIdentification());
	ids.back().setMetaValue("RT", 180);
	ids.back().setMetaValue("MZ", 123.456);
	
	file.getPrecursorRTandMZ(files_and_peptide_identification_with_scan_number, ids_found);
	
	TEST_REAL_EQUAL(ids_found.front().getMetaValue("RT"), ids.front().getMetaValue("RT"));
	TEST_REAL_EQUAL(ids_found.front().getMetaValue("MZ"), ids.front().getMetaValue("MZ"));
	TEST_REAL_EQUAL(ids_found.back().getMetaValue("RT"), ids.back().getMetaValue("RT"));
	TEST_REAL_EQUAL(ids_found.back().getMetaValue("MZ"), ids.back().getMetaValue("MZ"));
RESULT

CHECK(void getLabels(const String& source_database_filename, String& ac_label, String& sequence_start_label, String& sequence_end_label, String& comment_label, String& species_label))
	String ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getLabels("a", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label), "the file 'a' could not be found")
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.getLabels("data/Inspect_test1.fasta", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label), "data/Inspect_test1.fasta in: database has unknown file format (neither trie nor FASTA nor swissprot)")
	
	// test the actual program
	file.getLabels("data/Inspect_test.fasta", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
	TEST_STRING_EQUAL(ac_label, ">")
	TEST_STRING_EQUAL(sequence_start_label, ">")
	TEST_STRING_EQUAL(sequence_end_label, ">")
	TEST_STRING_EQUAL(comment_label, ";")
	TEST_STRING_EQUAL(species_label, ">")
RESULT

CHECK(vector< UInt > getWantedRecords(const String& result_filename, Real p_value_threshold))

	// test exceptions
	TEST_EXCEPTION(Exception::IllegalArgument, file.getWantedRecords("", 2.0))

	TEST_EXCEPTION(Exception::IllegalArgument, file.getWantedRecords("", -1.0))
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getWantedRecords("a", 0.01), "the file 'a' could not be found")
	
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileEmpty, file.getWantedRecords("data/Inspect_empty_file.txt", 0.01), "the file 'data/Inspect_empty_file.txt' is empty")
	

	// test the actual program
	vector< UInt > wanted_records = file.getWantedRecords("data/InspectOutfile.out", 0.01);
	TEST_EQUAL (wanted_records.size(), 1)
	if ( !wanted_records.empty() ) TEST_EQUAL (wanted_records.front(), 0)
RESULT

CHECK(template< typename PeakT > void getExperiment(MSExperiment< PeakT >& exp, String& type, const String& in_filename))
	MSExperiment< > exp;
	String type;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.getExperiment(exp, type, "data/InspectOutfile_version_file.txt"), "data/InspectOutfile_version_file.txt in: Could not determine type of the file. Aborting!")

	
	// test the actual program
	file.getExperiment(exp, type, "TOPP/Inspect.mzXML");
	TEST_STRING_EQUAL(type, "mzXML")
	file.getExperiment(exp, type, "TOPP/Inspect.mzData");
	TEST_STRING_EQUAL(type, "mzData")
RESULT

CHECK(void getSearchEngineAndVersion(const String& inspect_output_without_parameters_filename, ProteinIdentification& protein_identification) )
	ProteinIdentification protein_identification;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getSearchEngineAndVersion("a", protein_identification), "the file 'a' could not be found")

	protein_identification.setHits(vector< ProteinHit >());
	
	
	// test the actual program
	file.getSearchEngineAndVersion("data/InspectOutfile_version_file.txt", protein_identification);
	TEST_STRING_EQUAL(protein_identification.getSearchEngine(), "InsPecT");
	TEST_STRING_EQUAL(protein_identification.getSearchEngineVersion(), "20060907");
RESULT

CHECK(void readOutHeader(const String& filename, const String& header_line, Int& spectrum_file_column, Int& scan_column, Int& peptide_column, Int& protein_column, Int& charge_column, Int& MQ_score_column, Int& p_value_column, Int& record_number_column, Int& DB_file_pos_column, Int& spec_file_pos_column, UInt &number_of_columns))
		
	String header_line = "#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	";
	
	Int spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column;
	UInt number_of_columns;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.readOutHeader("dummy_testfile", header_line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns), "dummy_testfile in: at least one of the columns '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge', 'MQScore', 'p-value', 'RecordNumber', 'DBFilePos' or 'SpecFilePos' is missing!")

	
	// test the actual program
	header_line = "#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos";
	
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
RESULT

END_TEST
