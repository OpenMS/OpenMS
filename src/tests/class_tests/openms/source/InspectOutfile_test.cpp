// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

String spectrum_file1 = OPENMS_GET_TEST_DATA_PATH("InspectOutfile_test_1.mzXML");
String spectrum_file2 = OPENMS_GET_TEST_DATA_PATH("InspectOutfile_test_2.mzXML");

//create input file (they contain absolute paths...)
String input_file_name;
NEW_TMP_FILE(input_file_name);
TextFile outfile_content;
outfile_content.addLine("#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	CutScore	IntenseBY	BYPresent	Unused	p-value	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos");
outfile_content.addLine(spectrum_file1 + "	4	N.EER.N	Q9CQV8|1433B_MOUSE	1	-6889	1150	200	500	0	0.95000	0	0	1	287	270451");
outfile_content.addLine(spectrum_file1 + "	4	R.EKIE.K	P68509|1433F_BOVIN	1	-1456	3388	300	667	0	0.01000	-1199	-1199	0	87	276369");
outfile_content.addLine(spectrum_file1 + "	4	E.KKLE.K	P68509|1433F_BOVIN	1	-257	3820	450	333	0	0.00001	1199	1199	0	78	276369");
outfile_content.addLine(spectrum_file1 + "	25	D.DAIAE.L	P68509|1433F_BOVIN	1	-6526	354	250	375	0	0.97500	-6269	-6269	0	205	276369");
outfile_content.addLine(spectrum_file1 + "	37	D.EAIAEL.D	Q9CQV8|1433B_MOUSE	1	-2807	2464	350	300	0	0.98536	0	0	1	446	387302");
outfile_content.addLine(spectrum_file1 + "	49	D.KFLIK.N	P68509|1433F_BOVIN	1	-5308	1813	250	375	0	0.98758	0	0	0	106	473207");
outfile_content.addLine(spectrum_file1 + "	120	N.EKKLEKVKA.Y	P68509|1433F_BOVIN	2	-9115	-1187	185	125	0	0.96296	0	0	0	77	1062539");
outfile_content.addLine(spectrum_file1 + "	217	N.EDRNLL.S	P68509|1433F_BOVIN	1	-11835	140	100	0	0	0.97143	0	0	0	41	1720635");
outfile_content.addLine(spectrum_file1 + "	249	A.LLDKF.L	P68509|1433F_BOVIN	1	-4503	2554	300	250	0	0.99043	0	0	0	103	1943362");
outfile_content.addLine(spectrum_file1 + "	501	F.DEAIAELDTLNEE.S	Q9CQV8|1433B_MOUSE	2	-7874	-1829	231	0	0	0.95652	0	0	1	452	3860094");
outfile_content.addLine(spectrum_file1 + "	667	K.LAEQAERYDDMAA.A	Q9CQV8|1433B_MOUSE	2	-9787	-1362	77	125	0	0.96552	0	0	1	266	5279013");
outfile_content.addLine(spectrum_file1 + "	685	N.LTLWTSENQGDEGDAG.E	Q9CQV8|1433B_MOUSE	2	-9056	-2174	125	0	0	0.96296	0	0	1	479	5448607");
outfile_content.addLine(spectrum_file1 + "	736	Y.QEAFEIS.K	Q9CQV8|1433B_MOUSE	1	-11507	478	0	250	0	0.98789	17	17	1	399	6018155");
outfile_content.addLine(spectrum_file1 + "	736	Y.KEAFEIS.K	P68509|1433F_BOVIN	1	-11524	457	0	250	0	0.97059	-17	0	0	156	6018155");
outfile_content.addLine(spectrum_file1 + "	758	S.NEDRNLLSVAYKN.V	P68509|1433F_BOVIN	2	-4841	-481	256	208	0	0.98948	0	0	0	43	6167475");
outfile_content.addLine(spectrum_file1 + "	764	L.LAKQAFDDAIAELDTLNED.S	P68509|1433F_BOVIN	2	-4474	-560	140	250	0	0.99043	0	0	0	201	6208454");
outfile_content.addLine(spectrum_file1 + "	786	T.MDKSELV.Q	Q9CQV8|1433B_MOUSE	1	-12849	-1629	48	167	0	0.97368	0	0	1	253	6399323");
outfile_content.addLine(spectrum_file1 + "	1045	S.VFYYEI.Q	P68509|1433F_BOVIN	1	-15475	-2579	0	100	0	0.97500	0	0	0	184	9842931");
outfile_content.addLine(spectrum_file1 + "	1962	E.AFEIS.K	P68509|1433F_BOVIN	1	-10496	-742	100	375	0	0.96774	0	0	0	159	19098337");
outfile_content.store(input_file_name);

String input_file_name2;
NEW_TMP_FILE(input_file_name2);
TextFile outfile_content2;
outfile_content2.addLine("#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	CutScore	IntenseBY	BYPresent	Unused	p-value	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos");
outfile_content2.addLine(spectrum_file1 + "	N.EER.N	Q9CQV8|1433B_MOUSE	1	-6889	1150	200	500	0	0.95000	0	0	1	287	270451");
outfile_content2.store(input_file_name2);


InspectOutfile* ptr = nullptr;
InspectOutfile* nullPointer = nullptr;
START_SECTION(InspectOutfile())
	ptr = new InspectOutfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~InspectOutfile())
	delete ptr;
END_SECTION

START_SECTION((InspectOutfile& operator=(const InspectOutfile &inspect_outfile)))
	InspectOutfile inspect_outfile1;
  InspectOutfile inspect_outfile2;
	inspect_outfile2 = inspect_outfile1;
	InspectOutfile inspect_outfile3;
	inspect_outfile1 = InspectOutfile();
	TEST_EQUAL( (inspect_outfile2 == inspect_outfile3 ), true)
END_SECTION

START_SECTION((InspectOutfile(const InspectOutfile &inspect_outfile)))
	InspectOutfile inspect_outfile1;
	InspectOutfile inspect_outfile2(inspect_outfile1);
	InspectOutfile inspect_outfile3;
	inspect_outfile1 = InspectOutfile();
	TEST_EQUAL( (inspect_outfile2 == inspect_outfile3 ), true)
END_SECTION

START_SECTION((bool operator==(const InspectOutfile &inspect_outfile) const))
	InspectOutfile inspect_outfile1;
	InspectOutfile inspect_outfile2;
	TEST_EQUAL(( inspect_outfile1 == inspect_outfile2 ), true)
END_SECTION

InspectOutfile file;

START_SECTION(std::vector< Size > load(const String& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const double p_value_threshold, const String& database_filename = ""))
	vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;

	// test exceptions
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification, 2.0))
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification, -1.0))
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.load("a", peptide_identifications, protein_identification, 0.01), "the file 'a' could not be found")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileEmpty, file.load(OPENMS_GET_TEST_DATA_PATH("Inspect_empty_file.txt"), peptide_identifications, protein_identification, 0.01), OPENMS_GET_TEST_DATA_PATH_MESSAGE("the file '","Inspect_empty_file.txt","' is empty"))

	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());

	file.load(input_file_name, peptide_identifications, protein_identification, 0.001);

	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 1)
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 0.001)
		if( peptide_identifications[0].getHits().size() == 1)
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), -257)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "KKLE")
            vector<PeptideEvidence> pes = peptide_identifications[0].getHits()[0].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'E')
            TEST_EQUAL(pes[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
            TEST_EQUAL(pes.size(), 1)
            TEST_STRING_EQUAL(pes[0].getProteinAccession(), "P68509")
		}
	}
	peptide_identifications.clear();

	file.load(input_file_name, peptide_identifications, protein_identification, 0.01);

    TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "Inspect")
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 0.01)
		if( peptide_identifications[0].getHits().size() == 2 )
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), -257)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "KKLE")
            vector<PeptideEvidence> pes = peptide_identifications[0].getHits()[0].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'E')
            TEST_EQUAL(pes[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 1)
            TEST_EQUAL(pes.size(), 1)
            TEST_STRING_EQUAL(pes[0].getProteinAccession(), "P68509")
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), -1456)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getSequence().toString(), "EKIE")
            vector<PeptideEvidence> pes1 = peptide_identifications[0].getHits()[1].getPeptideEvidences();
            TEST_EQUAL(pes1[0].getAABefore(), 'R')
            TEST_EQUAL(pes1[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 1)
            TEST_EQUAL(pes1.size(), 1)
            TEST_STRING_EQUAL(pes1[0].getProteinAccession(), "P68509")
		}
	}

	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());
	TEST_EQUAL(file.load(input_file_name2, peptide_identifications, protein_identification, 0.01).size(), 1)
	if ( file.load(input_file_name2, peptide_identifications, protein_identification, 0.01).size() == 1 )
	{
		TEST_EQUAL(file.load(input_file_name2, peptide_identifications, protein_identification, 0.01)[0], 2)
	}
END_SECTION

START_SECTION(void generateTrieDB(const std::String& source_database_filename, const std::String& database_filename, const std::String& index_filename, bool append = false, const std::String species = ""))
	// test exceptions
	// test file not found for input file
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.generateTrieDB("a", "", ""), "the file 'a' could not be found")

	// test unable to create file
	TEST_EXCEPTION(Exception::UnableToCreateFile, file.generateTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.fasta"), OPENMS_GET_TEST_DATA_PATH("Inspect_unreadable_unwriteable.txt"), ""))

	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");


	// test the actual program
	file.generateTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.fasta"), "InspectOutfile_test.trie", "InspectOutfile_test.index");
	TEST_FILE_EQUAL("InspectOutfile_test.trie", OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"))
	TEST_FILE_EQUAL("InspectOutfile_test.index", OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"))
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");

	file.generateTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test2.fasta"), "InspectOutfile_test.trie", "InspectOutfile_test.index");
	TEST_FILE_EQUAL("InspectOutfile_test.trie", OPENMS_GET_TEST_DATA_PATH("Inspect_test2.trie"))
	TEST_FILE_EQUAL("InspectOutfile_test.index", OPENMS_GET_TEST_DATA_PATH("Inspect_test2.index"))
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
END_SECTION

START_SECTION(void compressTrieDB(const String& database_filename, const String& index_filename, std::vector< Size >& wanted_records, const String& snd_database_filename, const String& snd_index_filename, bool append = false))
	vector< Size > wanted_records(1, 0);

	// test exceptions
	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
	remove("InspectOutfile_test2.trie");
	remove("InspectOutfile_test2.index");
	// test for equal filenames
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"), wanted_records, OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), ""), OPENMS_GET_TEST_DATA_PATH_MESSAGE("","Inspect_test.trie"," in: Same filename can not be used for original and second database!"))
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"), wanted_records, "", OPENMS_GET_TEST_DATA_PATH("Inspect_test.index")), OPENMS_GET_TEST_DATA_PATH_MESSAGE("","Inspect_test.index"," in: Same filename can not be used for original and second database!"))

	// test file not found for input files (using empty filenames)
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.compressTrieDB("a", "", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index"), "the file 'a' could not be found")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), "b", wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index"), "the file 'b' could not be found")

	// test for unable to create file
	TEST_EXCEPTION(Exception::UnableToCreateFile, file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"), wanted_records, OPENMS_GET_TEST_DATA_PATH("Inspect_unreadable_unwriteable.txt"), "", true))

	// test for parse error
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_empty_file.txt"), wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", true),OPENMS_GET_TEST_DATA_PATH_MESSAGE("", "Inspect_empty_file.txt"," in: index file is too short!"))

	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");
	wanted_records.clear();


	// test the actual program
	bool append = true;
	wanted_records.push_back(0);
	file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"), wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);
	wanted_records.clear();

	file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test2.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test2.index"), wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);

	wanted_records.clear();
	wanted_records.push_back(1);
	file.compressTrieDB(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"), wanted_records, "InspectOutfile_test.trie", "InspectOutfile_test.index", append);

	wanted_records.clear();
	wanted_records.push_back(0);
	wanted_records.push_back(3);
	file.compressTrieDB("InspectOutfile_test.trie", "InspectOutfile_test.index", wanted_records, "InspectOutfile_test2.trie", "InspectOutfile_test2.index");

	remove("InspectOutfile_test.trie");
	remove("InspectOutfile_test.index");

	TEST_FILE_EQUAL("InspectOutfile_test2.trie", OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"))
	TEST_FILE_EQUAL("InspectOutfile_test2.index", OPENMS_GET_TEST_DATA_PATH("Inspect_test.index"))

	remove("InspectOutfile_test2.trie");
	remove("InspectOutfile_test2.index");
END_SECTION


START_SECTION(std::vector< Size > getSequences(const String& database_filename, const std::map< Size, Size >& wanted_records, std::vector< String >& sequences))
	map< Size, Size > rn_position_map;
	rn_position_map[0] = 0;
	rn_position_map[1] = 1;
	vector< String > sequences, found_sequences;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getSequences("a", rn_position_map, found_sequences), "the file 'a' could not be found")

	found_sequences.clear();

	// test actual program
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");

	file.getSequences(OPENMS_GET_TEST_DATA_PATH("Inspect_test.trie"), rn_position_map, found_sequences);
	if ( found_sequences.size() == 1 )
	{
		TEST_STRING_EQUAL(found_sequences.front(), sequences.front())

		file.getSequences(OPENMS_GET_TEST_DATA_PATH("Inspect_test2.trie"), rn_position_map, found_sequences);

		TEST_EQUAL((sequences == found_sequences), true)
	}
	rn_position_map.clear();
	sequences.clear();
	found_sequences.clear();
END_SECTION

START_SECTION(void getACAndACType(String line, String& accession, String& accession_type))
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
END_SECTION

START_SECTION(void getPrecursorRTandMZ(const vector< pair< String, vector< pair < Size, Size > > > >& files_and_peptide_identification_with_scan_number, std::vector< PeptideIdentification >& ids))
	vector< pair< String, vector< pair< Size, Size > > > > files_and_peptide_identification_with_scan_number;
	vector< PeptideIdentification > ids, ids_found;

	// test exceptions
	files_and_peptide_identification_with_scan_number.push_back(make_pair(spectrum_file1, vector< pair< Size, Size > >(1, make_pair(0, 10))));
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.getPrecursorRTandMZ(files_and_peptide_identification_with_scan_number, ids_found), OPENMS_GET_TEST_DATA_PATH_MESSAGE("","InspectOutfile_test_1.mzXML"," in: Not enough scans in file! (4 available, should be at least 10)"))

	files_and_peptide_identification_with_scan_number.clear();
	ids.clear();
	ids_found.clear();

	files_and_peptide_identification_with_scan_number.push_back(make_pair(spectrum_file1, vector< pair < Size, Size > >(1, make_pair(0, 4))));
	files_and_peptide_identification_with_scan_number.push_back(make_pair(spectrum_file2, vector< pair < Size, Size > >(1, make_pair(1, 4))));
	ids_found.push_back(PeptideIdentification());
	ids_found.push_back(PeptideIdentification());

	ids.push_back(PeptideIdentification());
	ids.back().setRT(-1);
	ids.back().setMZ(123.456);
	ids.push_back(PeptideIdentification());
	ids.back().setRT(180);
	ids.back().setMZ(123.456);

	file.getPrecursorRTandMZ(files_and_peptide_identification_with_scan_number, ids_found);

	TEST_REAL_SIMILAR(ids_found.front().getRT(), ids.front().getRT());
	TEST_REAL_SIMILAR(ids_found.front().getMZ(), ids.front().getMZ());
	TEST_REAL_SIMILAR(ids_found.back().getRT(), ids.back().getRT());
	TEST_REAL_SIMILAR(ids_found.back().getMZ(), ids.back().getMZ());
END_SECTION

START_SECTION(void getLabels(const String& source_database_filename, String& ac_label, String& sequence_start_label, String& sequence_end_label, String& comment_label, String& species_label))
	String ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getLabels("a", ac_label, sequence_start_label, sequence_end_label, comment_label, species_label), "the file 'a' could not be found")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.getLabels(OPENMS_GET_TEST_DATA_PATH("Inspect_test1.fasta"), ac_label, sequence_start_label, sequence_end_label, comment_label, species_label), OPENMS_GET_TEST_DATA_PATH_MESSAGE("","Inspect_test1.fasta"," in: database has unknown file format (neither trie nor FASTA nor swissprot)"))

	// test the actual program
	file.getLabels(OPENMS_GET_TEST_DATA_PATH("Inspect_test.fasta"), ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
	TEST_STRING_EQUAL(ac_label, ">")
	TEST_STRING_EQUAL(sequence_start_label, ">")
	TEST_STRING_EQUAL(sequence_end_label, ">")
	TEST_STRING_EQUAL(comment_label, ";")
	TEST_STRING_EQUAL(species_label, ">")
END_SECTION

START_SECTION(vector< Size > getWantedRecords(const String& result_filename, double p_value_threshold))

	// test exceptions
	TEST_EXCEPTION(Exception::IllegalArgument, file.getWantedRecords("", 2.0))
	TEST_EXCEPTION(Exception::IllegalArgument, file.getWantedRecords("", -1.0))
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getWantedRecords("a", 0.01), "the file 'a' could not be found")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileEmpty, file.getWantedRecords(OPENMS_GET_TEST_DATA_PATH("Inspect_empty_file.txt"), 0.01), OPENMS_GET_TEST_DATA_PATH_MESSAGE("the file '","Inspect_empty_file.txt","' is empty"))


	// test the actual program
	vector< Size > wanted_records = file.getWantedRecords(input_file_name, 0.01);
	TEST_EQUAL (wanted_records.size(), 1)
	if ( !wanted_records.empty() ) TEST_EQUAL (wanted_records.front(), 0)
END_SECTION

START_SECTION(template< typename PeakT > void getExperiment(MSExperiment< PeakT >& exp, String& type, const String& in_filename))
	PeakMap exp;
	String type;

	// test the actual program
	file.getExperiment(exp, type, OPENMS_GET_TEST_DATA_PATH("Inspect.mzXML"));
	TEST_STRING_EQUAL(type, "mzXML")
	file.getExperiment(exp, type, OPENMS_GET_TEST_DATA_PATH("Inspect.mzData"));
	TEST_STRING_EQUAL(type, "mzData")
END_SECTION

START_SECTION(bool getSearchEngineAndVersion(const String& cmd_output, ProteinIdentification& protein_identification) )
	ProteinIdentification protein_identification;

	protein_identification.setHits(vector< ProteinHit >());

  String output = "\
InsPecT vesrion 20060907\
  Interpretation of Peptides with Post-translational Modifications.\
  Copyright 2006, The Regents of the University of California\
  [See Docs directory for usage manual and copyright information]\
\
\
Sample command-line:\
Inspect.exe -i Foo.in -o Foo.txt -e ErrorsFoo.txt\
Command-line arguments:\
 -i InputFileName: Path to a config file specifying search parameters.\
 -o OutputFileName: Output file for match results.  If not\
          specified, output goes to stdout.\
 -e ErrorFileName: Output file for errors and warnings, if any.  If not\
          specified, any errors go to Inspect.err; if there are no errors.\
          or warnings reported, this file will be erased at end of run.\
 -r ResourceDir: Directory for resource files (such \
     as AminoAcidMasses.txt).  Defaults to current directory. \
  Consult the documentation (Inspect.html) for further details.\
";

	// test the actual program
	TEST_EQUAL(file.getSearchEngineAndVersion(output, protein_identification), true);
	TEST_STRING_EQUAL(protein_identification.getSearchEngine(), "InsPecT");
	TEST_STRING_EQUAL(protein_identification.getSearchEngineVersion(), "20060907");

  output = "\
InsPecT version 20100331\
  Interpretation of Peptides with Post-translational Modifications.\
  Copyright 2007,2008,2009 The Regents of the University of California\
  [See Docs directory for usage manual and copyright information]\
\
\
Sample command-line:\
Inspect.exe -i Foo.in -o Foo.txt -e ErrorsFoo.txt\
Command-line arguments:\
 -i InputFileName: Path to a config file specifying search parameters.\
 -o OutputFileName: Output file for match results.  If not\
          specified, output goes to stdout.\
 -e ErrorFileName: Output file for errors and warnings, if any.  If not\
          specified, any errors go to Inspect.err; if there are no errors.\
          or warnings reported, this file will be erased at end of run.\
 -r ResourceDir: Directory for resource files (such\
     as AminoAcidMasses.txt).  Defaults to current directory.\
 -a AminoAcidMassesFile: Specify a file containing non-standard amino acid masses.\
  Consult the documentation (Inspect.html) for further details.";

	// test the actual program
	TEST_EQUAL(file.getSearchEngineAndVersion(output, protein_identification), true);
	TEST_STRING_EQUAL(protein_identification.getSearchEngine(), "InsPecT");
	TEST_STRING_EQUAL(protein_identification.getSearchEngineVersion(), "20100331");

END_SECTION

START_SECTION(void readOutHeader(const String& filename, const String& header_line, Int& spectrum_file_column, Int& scan_column, Int& peptide_column, Int& protein_column, Int& charge_column, Int& MQ_score_column, Int& p_value_column, Int& record_number_column, Int& DB_file_pos_column, Int& spec_file_pos_column, Size &number_of_columns))

	String header_line = "#SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	";

	Int spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column;
	Size number_of_columns;

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
END_SECTION

END_TEST
