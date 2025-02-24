// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;


START_TEST(String, "$Id$")


SequestOutfile* ptr = nullptr;
SequestOutfile* nullPointer = nullptr;
START_SECTION(SequestOutfile())
	ptr = new SequestOutfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~SequestOutfile())
	delete ptr;
END_SECTION

START_SECTION((SequestOutfile& operator=(const SequestOutfile &sequest_outfile)))
	SequestOutfile sequest_outfile1;
  SequestOutfile sequest_outfile2;
	sequest_outfile2 = sequest_outfile1;
	SequestOutfile sequest_outfile3;
	sequest_outfile1 = SequestOutfile();
	TEST_EQUAL(( sequest_outfile2 == sequest_outfile3 ), true)
END_SECTION

START_SECTION((SequestOutfile(const SequestOutfile &sequest_outfile)))
	SequestOutfile sequest_outfile1;
	SequestOutfile sequest_outfile2(sequest_outfile1);
	SequestOutfile sequest_outfile3;
	sequest_outfile1 = SequestOutfile();
	TEST_EQUAL(( sequest_outfile2 == sequest_outfile3 ), true)
END_SECTION

START_SECTION((bool operator==(const SequestOutfile &sequest_outfile) const))
	SequestOutfile sequest_outfile1;
	SequestOutfile sequest_outfile2;
	TEST_EQUAL(( sequest_outfile1 == sequest_outfile2 ), true)
END_SECTION


SequestOutfile file;

START_SECTION(void load(const String& result_filename, std::vector< PeptideIdentification >& peptide_identifications, ProteinIdentification& protein_identification, const double p_value_threshold, std::vector< double >& pvalues, const String& database="", const bool ignore_proteins_per_peptide=false))
	vector< PeptideIdentification > peptide_identifications;
	ProteinIdentification protein_identification;
	vector< double > pvalues;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.load("a", peptide_identifications, protein_identification, 0.01, pvalues), "the file 'a' could not be found")
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("SequestOutfile.out1"), peptide_identifications, protein_identification, 0.01, pvalues),  OPENMS_GET_TEST_DATA_PATH_MESSAGE("","SequestOutfile.out1", " in: Wrong number of columns in line 16! (11 present, should be 12)"))
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification, 2.0, pvalues))
	TEST_EXCEPTION(Exception::IllegalArgument, file.load("", peptide_identifications, protein_identification,-1.0, pvalues))

	peptide_identifications.clear();
	protein_identification.setHits(vector< ProteinHit >());
	pvalues.clear();


	// test the actual program
	file.load(OPENMS_GET_TEST_DATA_PATH("SequestOutfile2.out"), peptide_identifications, protein_identification, 1.0, pvalues);
	TEST_EQUAL(peptide_identifications.size(), 0)

	file.load(OPENMS_GET_TEST_DATA_PATH("SequestOutfile.out"), peptide_identifications, protein_identification, 1.0, pvalues);

	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 4)
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "SEQUEST")
		TEST_STRING_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 1.0)
		if ( peptide_identifications[0].getHits().size() == 4 )
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 0.05)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "ETQAWSIATILETLYDL")
            vector<PeptideEvidence> pes = peptide_identifications[0].getHits()[0].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'C')
            TEST_EQUAL(pes[0].getAAAfter(), '-')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[0].getMetaValue("RankSp")), "1/80")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[0].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("MH")), 1967.0013)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("DeltCn")), 0.0000)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("XCorr")), 1.5789)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("Sp")), 310.3)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("Sf")), 0.05)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[0].getMetaValue("Ions")), "18/64")
            set<String> protein_accessions = peptide_identifications[0].getHits()[0].extractProteinAccessionsSet();
            TEST_EQUAL(protein_accessions.size(), 3)
            if (protein_accessions.size() == 3 )
			{
                set<String>::const_iterator s_it = protein_accessions.begin();
                TEST_STRING_EQUAL(*s_it++, "2136928")
                TEST_STRING_EQUAL(*s_it++, "L10605")
                TEST_STRING_EQUAL(*s_it++, "P35574")
			}
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 0.04)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getSequence().toString(), "QVLNPLLVLIFIPLFDL")
            PeptideEvidence pe = peptide_identifications[0].getHits()[1].getPeptideEvidences()[0];
            TEST_EQUAL(pe.getAABefore(), 'M')
            TEST_EQUAL(pe.getAAAfter(), 'V')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[1].getMetaValue("RankSp")), "2/85")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[1].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("MH")), 1967.1985)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("DeltCn")), 0.0390)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("XCorr")), 1.5173)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("Sp")), 308.3)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("Sf")), 0.04)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[1].getMetaValue("Ions")), "19/64")
            TEST_EQUAL(peptide_identifications[0].getHits()[1].getPeptideEvidences().size(), 2)
            if ( peptide_identifications[0].getHits()[1].getPeptideEvidences().size() == 2 )
			{
                TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getPeptideEvidences()[0].getProteinAccession(), "P46029")
                TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getPeptideEvidences()[1].getProteinAccession(), "U32507")
			}
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[2].getScore(), 0.02)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[2].getSequence().toString(), "WVELGPSVLAGVGVMVLLI")
            pes = peptide_identifications[0].getHits()[2].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'L')
            TEST_EQUAL(pes[0].getAAAfter(), 'P')
			TEST_EQUAL(peptide_identifications[0].getHits()[2].getRank(), 3)
			TEST_EQUAL(peptide_identifications[0].getHits()[2].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[2].getMetaValue("RankSp")), "3/117")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[2].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[2].getMetaValue("MH")), 1968.1244)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[2].getMetaValue("DeltCn")), 0.0501)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[2].getMetaValue("XCorr")), 1.4998)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[2].getMetaValue("Sp")), 292.4)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[2].getMetaValue("Sf")), 0.02)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[2].getMetaValue("Ions")), "17/72")
            TEST_EQUAL(pes.size(), 1)
            if ( pes.size() == 1 )
			{
                TEST_STRING_EQUAL(pes[0].getProteinAccession(), "e148876")
			}
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[3].getScore(), 0.14)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[3].getSequence().toString(), "FDEITAMTGDGVNDAPALK")
            pes = peptide_identifications[0].getHits()[3].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'S')
            TEST_EQUAL(pes[0].getAAAfter(), 'K')
			TEST_EQUAL(peptide_identifications[0].getHits()[3].getRank(), 4)
			TEST_EQUAL(peptide_identifications[0].getHits()[3].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[3].getMetaValue("RankSp")), "4/1")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[3].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[3].getMetaValue("MH")), 1964.9275)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[3].getMetaValue("DeltCn")), 0.0627)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[3].getMetaValue("XCorr")), 1.4799)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[3].getMetaValue("Sp")), 530.9)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[3].getMetaValue("Sf")), 0.14)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[3].getMetaValue("Ions")), "24/72")
            TEST_EQUAL(pes.size(), 8)
            if ( pes.size() == 8 )
			{
                TEST_STRING_EQUAL(pes[0].getProteinAccession(), "P20647")
                TEST_STRING_EQUAL(pes[1].getProteinAccession(), "P04192")
                TEST_STRING_EQUAL(pes[2].getProteinAccession(), "67962")
                TEST_STRING_EQUAL(pes[3].getProteinAccession(), "67961")
                TEST_STRING_EQUAL(pes[4].getProteinAccession(), "109166")
                TEST_STRING_EQUAL(pes[5].getProteinAccession(), "224621")
                TEST_STRING_EQUAL(pes[6].getProteinAccession(), "X02814")
                TEST_STRING_EQUAL(pes[7].getProteinAccession(), "J04703")
			}
		}
	}

	peptide_identifications.clear();
	pvalues.push_back(0.001);
	pvalues.push_back(0.01);
	pvalues.push_back(0.05);
	pvalues.push_back(0.5);
	file.load(OPENMS_GET_TEST_DATA_PATH("SequestOutfile.out"), peptide_identifications, protein_identification, 0.01, pvalues);

	TEST_EQUAL(peptide_identifications.size(), 1)
	if ( peptide_identifications.size() == 1 )
	{
		TEST_STRING_EQUAL(peptide_identifications[0].getScoreType(), "SEQUEST")
		TEST_STRING_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
		TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)
		TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 0.01)
		if ( peptide_identifications[0].getHits().size() == 2 )
		{
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 0.05)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[0].getSequence().toString(), "ETQAWSIATILETLYDL")
            vector<PeptideEvidence> pes = peptide_identifications[0].getHits()[0].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'C')
            TEST_EQUAL(pes[0].getAAAfter(), '-')
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
			TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[0].getMetaValue("RankSp")), "1/80")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[0].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("MH")), 1967.0013)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("DeltCn")), 0.0000)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("XCorr")), 1.5789)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("Sp")), 310.3)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[0].getMetaValue("Sf")), 0.05)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[0].getMetaValue("Ions")), "18/64")
            TEST_EQUAL(pes.size(), 3)
            if ( pes.size() == 3 )
			{
                TEST_STRING_EQUAL(pes[0].getProteinAccession(), "P35574")
                TEST_STRING_EQUAL(pes[1].getProteinAccession(), "2136928")
                TEST_STRING_EQUAL(pes[2].getProteinAccession(), "L10605")
			}
			TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 0.04)
			TEST_STRING_EQUAL(peptide_identifications[0].getHits()[1].getSequence().toString(), "QVLNPLLVLIFIPLFDL")
            pes = peptide_identifications[0].getHits()[1].getPeptideEvidences();
            TEST_EQUAL(pes[0].getAABefore(), 'M')
            TEST_EQUAL(pes[0].getAAAfter(), 'V')
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
			TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 3)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[1].getMetaValue("RankSp")), "2/85")
			TEST_EQUAL(static_cast<Int>(peptide_identifications[0].getHits()[1].getMetaValue("SequestId")), 0)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("MH")), 1967.1985)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("DeltCn")), 0.0390)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("XCorr")), 1.5173)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("Sp")), 308.3)
			TEST_REAL_SIMILAR(static_cast<double>(peptide_identifications[0].getHits()[1].getMetaValue("Sf")), 0.04)
			TEST_STRING_EQUAL(static_cast<String>(peptide_identifications[0].getHits()[1].getMetaValue("Ions")), "19/64")
            TEST_EQUAL(pes.size(), 2)
            if ( pes.size() == 2 )
			{
                TEST_STRING_EQUAL(pes[0].getProteinAccession(), "P46029")
                TEST_STRING_EQUAL(pes[1].getProteinAccession(), "U32507")
			}
		}
		TEST_STRING_EQUAL(peptide_identifications[0].getIdentifier(), "TurboSEQUEST_2004-03-16")
	}

	TEST_STRING_EQUAL(protein_identification.getSearchEngine(), "TurboSEQUEST")
	TEST_STRING_EQUAL(protein_identification.getSearchEngineVersion(), "v.27 (rev. 12)")
	TEST_STRING_EQUAL(protein_identification.getIdentifier(), "TurboSEQUEST_2004-03-16")
END_SECTION

START_SECTION(bool getColumns(const String& line, vector< String >& substrings, Size number_of_columns, Size reference_column))
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
END_SECTION

START_SECTION(void getSequences(const String& database_filename, const map< String, Size >& ac_position_map, vector< String >& sequences, vector< pair< String, Size > >& found, map< String, Size >& not_found))
	map< String, Size > ac_position_map, not_found;
	vector< String > sequences, found_sequences;
	vector< pair< String, Size > > found;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.getSequences("a", not_found, found_sequences, found, not_found), "the file 'a' could not be found")

	// test the actual program
	ac_position_map["P02666"] = 0;
	ac_position_map["Q9CQV8"] = 1;
	ac_position_map["Q5EEQ7"] = 2;
	ac_position_map["P68509"] = 3;

	sequences.push_back("MKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV");
	sequences.push_back("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN");
	sequences.push_back("SAPPSLLVLYFGKKELRAMKVLILACLVALALARELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPP");
	sequences.push_back("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");

	ABORT_IF(ac_position_map.size() != 4)
	file.getSequences(OPENMS_GET_TEST_DATA_PATH("Sequest_test.fasta"), ac_position_map, found_sequences, found, not_found);
	ABORT_IF(ac_position_map.size() != 4)
	TEST_EQUAL(found.size(), 2)
	TEST_EQUAL(not_found.size(), 2)
	ABORT_IF( found.size() != 2 || not_found.size() != 2 )

	TEST_EQUAL(String("P68509"), found[0].first)
	TEST_EQUAL(ac_position_map["P68509"], found[0].second)
	TEST_EQUAL(sequences[ac_position_map["P68509"]], found_sequences[0])

	TEST_EQUAL(String("Q9CQV8"), found[1].first)
	TEST_EQUAL(ac_position_map["Q9CQV8"], found[1].second)
	TEST_EQUAL(sequences[ac_position_map["Q9CQV8"]], found_sequences[1])

	// create a copy as getSequences() does some weird things with the actual map
	map< String, Size > ac_position_map_subset = not_found;
	file.getSequences(OPENMS_GET_TEST_DATA_PATH("Sequest_test2.fasta"), ac_position_map_subset, found_sequences, found, not_found);
	TEST_EQUAL(found.size(), 4)
	TEST_EQUAL(not_found.size(), 0)
	ABORT_IF(found.size() != 4 || !not_found.empty())

	TEST_EQUAL(String("P02666"), found[2].first)
	TEST_EQUAL(ac_position_map["P02666"], found[2].second)
	TEST_EQUAL(sequences[ac_position_map["P02666"]], found_sequences[2])

	TEST_EQUAL(String("Q5EEQ7"), found[3].first)
	TEST_EQUAL(ac_position_map["Q5EEQ7"], found[3].second)
	TEST_EQUAL(sequences[ac_position_map["Q5EEQ7"]], found_sequences[3])

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

START_SECTION(void readOutHeader(const String& result_filename, DateTime& datetime, double& precursor_mz_value, Int& charge, Size& precursor_mass_type, Size& ion_mass_type, Size& displayed_peptides, String& sequest, String& sequest_version, String& database_type, Int& number_column, Int& rank_sp_column, Int& id_column, Int& mh_column, Int& delta_cn_column, Int& xcorr_column, Int& sp_column, Int& sf_column, Int& ions_column, Int& reference_column, Int& peptide_column, Int& score_column, Size& number_of_columns))

	String result_filename = OPENMS_GET_TEST_DATA_PATH("Sequest.mzXML.13.1.d.out");
	DateTime datetime;

	double precursor_mz_value(0.0);
	Int
	 charge(-1),
	 number_column(-1),
	 rank_sp_column(-1),
	 id_column(-1),
	 mh_column(-1),
	 delta_cn_column(-1),
	 xcorr_column(-1),
	 sp_column(-1),
	 sf_column(-1),
	 ions_column(-1),
	 reference_column(-1),
	 peptide_column(-1),
	 score_column(-1);

	Size
	 precursor_mass_type(0),
	 ion_mass_type(0),
	 displayed_peptides(0),
	 number_of_columns(0);

	String
	 sequest,
	 sequest_version,
	 database_type;

	// test exceptions
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.readOutHeader("a", datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, displayed_peptides, sequest, sequest_version, database_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns), "the file 'a' could not be found")

	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.readOutHeader(OPENMS_GET_TEST_DATA_PATH("SequestOutfile_headerfile.txt"), datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, displayed_peptides, sequest, sequest_version, database_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns),  OPENMS_GET_TEST_DATA_PATH_MESSAGE("","SequestOutfile_headerfile.txt", " in: No Sequest version found!"))


	// test the actual program
	file.readOutHeader(result_filename, datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, displayed_peptides, sequest, sequest_version, database_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns);

	TOLERANCE_ABSOLUTE(0.0001)
	TEST_REAL_SIMILAR(precursor_mz_value, 866.606)

	TEST_STRING_EQUAL(sequest, "TurboSEQUEST")
	TEST_STRING_EQUAL(sequest_version, "v.27 (rev. 12)")
	TEST_STRING_EQUAL(database_type, "amino acids")
	TEST_STRING_EQUAL(datetime.get(), "2007-01-17 17:29:00")

	TEST_EQUAL(charge, 2)
	TEST_EQUAL(number_column, 0)
	TEST_EQUAL(rank_sp_column, 1)
	TEST_EQUAL(id_column, 2)
	TEST_EQUAL(mh_column, 3)
	TEST_EQUAL(delta_cn_column, 4)
	TEST_EQUAL(xcorr_column, 5)
	TEST_EQUAL(sp_column, 6)
	TEST_EQUAL(sf_column, 7)
	TEST_EQUAL(ions_column, 9)
	TEST_EQUAL(reference_column, 10)
	TEST_EQUAL(peptide_column, 11)
	TEST_EQUAL(score_column, 7)
	TEST_EQUAL(number_of_columns, 12)
	TEST_EQUAL(precursor_mass_type, 0)
	TEST_EQUAL(ion_mass_type, 0)
	TEST_EQUAL(displayed_peptides, 2)
END_SECTION

END_TEST
