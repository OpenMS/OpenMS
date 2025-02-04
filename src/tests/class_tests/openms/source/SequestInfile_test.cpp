// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/SequestInfile.h>

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SequestInfile, "$Id$")

/////////////////////////////////////////////////////////////

SequestInfile* ptr = nullptr;
SequestInfile* nullPointer = nullptr;
START_SECTION(SequestInfile())
	ptr = new SequestInfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~SequestInfile())
	delete ptr;
END_SECTION

START_SECTION((SequestInfile& operator=(const SequestInfile &sequest_infile)))
	SequestInfile sequest_infile1;
	sequest_infile1.setDatabase("dummy");
	SequestInfile sequest_infile2 = sequest_infile1;
  TEST_EQUAL(( sequest_infile1 == sequest_infile2 ), true)

  SequestInfile sequest_infile3;
  TEST_EQUAL(( sequest_infile2 == sequest_infile3 ), false)
END_SECTION

START_SECTION((SequestInfile(const SequestInfile &sequest_infile)))
	SequestInfile sequest_infile1;
	sequest_infile1.setDatabase("dummy");
	SequestInfile sequest_infile2(sequest_infile1);
  TEST_EQUAL(( sequest_infile1 == sequest_infile2 ), true)

	SequestInfile sequest_infile3;
  TEST_EQUAL(( sequest_infile2 == sequest_infile3 ), false)
END_SECTION

START_SECTION((bool operator==(const SequestInfile &sequest_infile) const))
	SequestInfile sequest_infile1;
	sequest_infile1.setDatabase("dummy");
	SequestInfile sequest_infile2;
  sequest_infile2.setDatabase("dummy");
	TEST_EQUAL(( sequest_infile1 == sequest_infile2 ), true)

  SequestInfile sequest_infile3;
  sequest_infile3.setDatabase("another dummy");
  TEST_EQUAL(( sequest_infile1 == sequest_infile3 ), false)
END_SECTION

SequestInfile file;

stringstream ss;
ss << "[SEQUEST_ENZYME_INFO]" << endl;
ss << "0.  AspN                    0     D             -" << endl;
ss << "1.  AspN_DE                 0     DE            -" << endl;
ss << "2.  Chymotrypsin            1     FWYL          -" << endl;
ss << "3.  Chymotrypsin_WYF        1     FWY           -" << endl;
ss << "4.  Clostripain             1     R             -" << endl;
ss << "5.  Cyanogen_Bromide        1     M             -" << endl;
ss << "6.  Elastase                1     ALIV          P" << endl;
ss << "7.  Elastase/Tryp/Chymo     1     ALIVKRWFY     P" << endl;
ss << "8.  GluC                    1     E             -" << endl;
ss << "9.  GluC_ED                 1     ED            -" << endl;
ss << "10.  IodosoBenzoate          1     W             -" << endl;
ss << "11.  LysC                    1     K             -" << endl;
ss << "12.  No_Enzyme               0     -             -" << endl;
ss << "13.  Proline_Endopept        1     P             -" << endl;
ss << "14.  Trypsin                 1     KRLNH         -" << endl;
ss << "15.  Trypsin/Chymo           1     KRLFWYN       -" << endl;
ss << "16.  Trypsin_Strict          1     KR            -" << endl;

START_SECTION((const String getEnzymeInfoAsString() const))
	TEST_STRING_EQUAL(file.getEnzymeInfoAsString(), ss.str())
END_SECTION

START_SECTION(void addEnzymeInfo(std::vector< String >& enzyme_info))
	std::vector< String > e_info;
	e_info.push_back("Z_TestEnzyme");
	e_info.push_back("1");
	e_info.push_back("RMW");
	e_info.push_back("-");
	file.addEnzymeInfo(e_info);
	e_info.clear();
	ss << "17.  Z_TestEnzyme            1     RMW           -" << endl;
	TEST_STRING_EQUAL(file.getEnzymeInfoAsString(), ss.str())
END_SECTION

START_SECTION(void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic))

	// test exceptions
	String modification_line = "Phosphorylation";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.handlePTMs(modification_line, "a", true), "the file 'a' could not be found")
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotReadable, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Sequest_unreadable_unwriteable.txt"), true), "the file `data/Sequest_unreadable_unwriteable.txt' is not readable for the current user")
	
	modification_line = "2H20,KRLNH,fix";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), true), "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix")

	modification_line = "10.3+";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), true), "No residues for modification given. Aborting! in: 10.3+")

	modification_line = "10.3+,KRLNH,stat,PTM_0";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), true), "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0")

	modification_line = "Phosphorylation:Phosphorylation";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), true), "There's already a modification with this name. Aborting! in: Phosphorylation")
	

	// test the actual program
	modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";

	// average masses
  file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), false);

	map< String, vector< String > > modifications;
	modifications["PTM_0"] = vector< String >(3);
	modifications["PTM_0"][0] = "KRLNH";
	modifications["PTM_0"][1] = "10.3";
	modifications["PTM_0"][2] = "FIX";
//	modifications["Phosphorylation"] = vector< String >(3);
//	modifications["Phosphorylation"][0] = "STYDHCR";
//	modifications["Phosphorylation"][1] = "79.97990";
//	modifications["Phosphorylation"][2] = "OPT";
	modifications["PTM_1"] = vector< String >(3);
	modifications["PTM_1"][0] = "C";
	modifications["PTM_1"][1] = "16";
	modifications["PTM_1"][2] = "OPT";
// 	modifications["Carbamylation"] = vector< String >(3);
// 	modifications["Carbamylation"][0] = "NTERM";
// 	modifications["Carbamylation"][1] = "43.02474";
// 	modifications["Carbamylation"][2] = "OPT";
//	modifications["Methylation"] = vector< String >(3);
//	modifications["Methylation"][0] = "CHKNQRILDEST";
//	modifications["Methylation"][1] = "14.02658";
//	modifications["Methylation"][2] = "OPT";
// 	modifications["PTM_5"] = vector< String >(3);
// 	modifications["PTM_5"][0] = "CTERM";
// 	modifications["PTM_5"][1] = "-16";
// 	modifications["PTM_5"][2] = "OPT";
// 	modifications["PTM_6"] = vector< String >(3);
// 	modifications["PTM_6"][0] = "NTERM";
// 	modifications["PTM_6"][1] = "-16";
// 	modifications["PTM_6"][2] = "FIX";
// 	modifications["PTM_7"] = vector< String >(3);
// 	modifications["PTM_7"][0] = "CTERM_PROT";
// 	modifications["PTM_7"][1] = "-17";
// 	modifications["PTM_7"][2] = "OPT";
// 	modifications["PTM_8"] = vector< String >(3);
// 	modifications["PTM_8"][0] = "NTERM_PROT";
// 	modifications["PTM_8"][1] = "-17";
// 	modifications["PTM_8"][2] = "FIX";
	modifications["PTM_2"] = vector< String >(3);
	modifications["PTM_2"][0] = "CTERM";
	modifications["PTM_2"][1] = "-16";
	modifications["PTM_2"][2] = "OPT";
	modifications["PTM_3"] = vector< String >(3);
	modifications["PTM_3"][0] = "NTERM";
	modifications["PTM_3"][1] = "-16";
	modifications["PTM_3"][2] = "FIX";
	modifications["PTM_4"] = vector< String >(3);
	modifications["PTM_4"][0] = "CTERM_PROT";
	modifications["PTM_4"][1] = "-17";
	modifications["PTM_4"][2] = "OPT";
	modifications["PTM_5"] = vector< String >(3);
	modifications["PTM_5"][0] = "NTERM_PROT";
	modifications["PTM_5"][1] = "-17";
	modifications["PTM_5"][2] = "FIX";

	map< String, vector< String > >::const_iterator result_mod_i = file.getModifications().begin();
	TEST_EQUAL(file.getModifications().size(), modifications.size())
	if ( file.getModifications().size() == modifications.size() )
	{
		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
		{
			TEST_STRING_EQUAL(result_mod_i->first, mod_i->first)
			TEST_EQUAL(result_mod_i->second.size(), 3)
			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
			if ( result_mod_i->second.size() == mod_i->second.size() )
			{
				TEST_STRING_EQUAL(result_mod_i->second[0], mod_i->second[0])
				TEST_STRING_EQUAL(result_mod_i->second[1], mod_i->second[1])
				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}

	// monoisotopic masses
	file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), true);

//	modifications["Phosphorylation"][1] = "79.96634";
// 	modifications["Carbamylation"][1] = "43.00582";
//	modifications["Methylation"][1] = "14.01565";

	result_mod_i = file.getModifications().begin();
	TEST_EQUAL(file.getModifications().size(), modifications.size())
	if ( file.getModifications().size() == modifications.size() )
	{
		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
		{
			TEST_STRING_EQUAL(result_mod_i->first, mod_i->first)
			TEST_EQUAL(result_mod_i->second.size(), 3)
			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
			if ( result_mod_i->second.size() == mod_i->second.size() )
			{
				TEST_STRING_EQUAL(result_mod_i->second[0], mod_i->second[0])
				TEST_STRING_EQUAL(result_mod_i->second[1], mod_i->second[1])
				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}
END_SECTION

START_SECTION((const std::map< String, std::vector< String > >& getModifications() const))
	String modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";

	// average masses
	file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Sequest_PTMs.xml"), false);

	map< String, vector< String > > modifications;
	modifications["PTM_0"] = vector< String >(3);
	modifications["PTM_0"][0] = "KRLNH";
	modifications["PTM_0"][1] = "10.3";
	modifications["PTM_0"][2] = "FIX";
//	modifications["Phosphorylation"] = vector< String >(3);
//	modifications["Phosphorylation"][0] = "STYDHCR";
//	modifications["Phosphorylation"][1] = "79.97990";
//	modifications["Phosphorylation"][2] = "OPT";
	modifications["PTM_1"] = vector< String >(3);
	modifications["PTM_1"][0] = "C";
	modifications["PTM_1"][1] = "16";
	modifications["PTM_1"][2] = "OPT";
// 	modifications["Carbamylation"] = vector< String >(3);
// 	modifications["Carbamylation"][0] = "NTERM";
// 	modifications["Carbamylation"][1] = "43.02474";
// 	modifications["Carbamylation"][2] = "OPT";
//	modifications["Methylation"] = vector< String >(3);
//	modifications["Methylation"][0] = "CHKNQRILDEST";
//	modifications["Methylation"][1] = "14.02658";
//	modifications["Methylation"][2] = "OPT";
// 	modifications["PTM_5"] = vector< String >(3);
// 	modifications["PTM_5"][0] = "CTERM";
// 	modifications["PTM_5"][1] = "-16";
// 	modifications["PTM_5"][2] = "OPT";
// 	modifications["PTM_6"] = vector< String >(3);
// 	modifications["PTM_6"][0] = "NTERM";
// 	modifications["PTM_6"][1] = "-16";
// 	modifications["PTM_6"][2] = "FIX";
// 	modifications["PTM_7"] = vector< String >(3);
// 	modifications["PTM_7"][0] = "CTERM_PROT";
// 	modifications["PTM_7"][1] = "-17";
// 	modifications["PTM_7"][2] = "OPT";
// 	modifications["PTM_8"] = vector< String >(3);
// 	modifications["PTM_8"][0] = "NTERM_PROT";
// 	modifications["PTM_8"][1] = "-17";
// 	modifications["PTM_8"][2] = "FIX";
	modifications["PTM_2"] = vector< String >(3);
	modifications["PTM_2"][0] = "CTERM";
	modifications["PTM_2"][1] = "-16";
	modifications["PTM_2"][2] = "OPT";
	modifications["PTM_3"] = vector< String >(3);
	modifications["PTM_3"][0] = "NTERM";
	modifications["PTM_3"][1] = "-16";
	modifications["PTM_3"][2] = "FIX";
	modifications["PTM_4"] = vector< String >(3);
	modifications["PTM_4"][0] = "CTERM_PROT";
	modifications["PTM_4"][1] = "-17";
	modifications["PTM_4"][2] = "OPT";
	modifications["PTM_5"] = vector< String >(3);
	modifications["PTM_5"][0] = "NTERM_PROT";
	modifications["PTM_5"][1] = "-17";
	modifications["PTM_5"][2] = "FIX";

	map< String, vector< String > >::const_iterator result_mod_i = file.getModifications().begin();
	TEST_EQUAL(file.getModifications().size(), modifications.size())
	if ( file.getModifications().size() == modifications.size() )
	{
		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
		{
			TEST_STRING_EQUAL(result_mod_i->first, mod_i->first)
			TEST_EQUAL(result_mod_i->second.size(), 3)
			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
			if ( result_mod_i->second.size() == mod_i->second.size() )
			{
				TEST_STRING_EQUAL(result_mod_i->second[0], mod_i->second[0])
				TEST_STRING_EQUAL(result_mod_i->second[1], mod_i->second[1])
				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}
END_SECTION


START_SECTION(void setDatabase(const String& database))
	file.setDatabase(R"(\\bude\langwisc\sequest_test\Analysis.mzXML.fasta)");
	TEST_STRING_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
END_SECTION

START_SECTION((const String& getDatabase() const))
	TEST_STRING_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
END_SECTION

START_SECTION(void setNeutralLossesForIons(const String& neutral_losses_for_ions))
	file.setNeutralLossesForIons("0 1 1");
	TEST_STRING_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
END_SECTION

START_SECTION((const String& getNeutralLossesForIons() const))
	TEST_STRING_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
END_SECTION

START_SECTION(void setIonSeriesWeights(const String& ion_series_weights))
	file.setIonSeriesWeights("0 1.0 0 0 0 0 0 1.0 0");
	TEST_STRING_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
END_SECTION

START_SECTION((const String& getIonSeriesWeights() const))
	TEST_STRING_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
END_SECTION

START_SECTION(void setPartialSequence(const String& partial_sequence))
	file.setPartialSequence("SEQVEST TEST");
	TEST_STRING_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
END_SECTION

START_SECTION((const String& getPartialSequence() const))
	TEST_STRING_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
END_SECTION

START_SECTION(void setSequenceHeaderFilter(const String& sequence_header_filter))
	file.setSequenceHeaderFilter("homo~sapiens !mus musculus");
	TEST_STRING_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
END_SECTION

START_SECTION((const String& getSequenceHeaderFilter() const))
	TEST_STRING_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
END_SECTION


START_SECTION(void setPrecursorMassTolerance(float precursor_mass_tolerance))
	file.setPrecursorMassTolerance(1.3f);
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance() , 1.3)
END_SECTION

START_SECTION((float getPrecursorMassTolerance() const))
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance() , 1.3)
END_SECTION

START_SECTION(void setPeakMassTolerance(float peak_mass_tolerance))
	file.setPeakMassTolerance(0.3f);
	TEST_REAL_SIMILAR(file.getPeakMassTolerance() , 0.3)
END_SECTION

START_SECTION((float getPeakMassTolerance() const))
	TEST_REAL_SIMILAR(file.getPeakMassTolerance() , 0.3)
END_SECTION

START_SECTION(void setMatchPeakTolerance(float match_peak_tolerance))
	file.setMatchPeakTolerance(1.2f);
	TEST_REAL_SIMILAR(file.getMatchPeakTolerance() , 1.2)
END_SECTION

START_SECTION((float getMatchPeakTolerance() const))
	TEST_REAL_SIMILAR(file.getMatchPeakTolerance() , 1.2)
END_SECTION

START_SECTION(void setIonCutoffPercentage(float ion_cutoff_percentage))
	file.setIonCutoffPercentage(0.3f);
	TEST_REAL_SIMILAR(file.getIonCutoffPercentage() , 0.3)
END_SECTION

START_SECTION((float getIonCutoffPercentage() const))
	TEST_REAL_SIMILAR(file.getIonCutoffPercentage() , 0.3)
END_SECTION

START_SECTION(void setProteinMassFilter(const String& protein_mass_filter))
	file.setProteinMassFilter("30.2 0");
	TEST_STRING_EQUAL(file.getProteinMassFilter() , "30.2 0")
END_SECTION

START_SECTION((const String& getProteinMassFilter() const))
	TEST_STRING_EQUAL(file.getProteinMassFilter() , "30.2 0")
END_SECTION

START_SECTION(void setPeptideMassUnit(Size peptide_mass_unit))
	file.setPeptideMassUnit(0);
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
END_SECTION

START_SECTION((Size getPeptideMassUnit() const))
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
END_SECTION

START_SECTION(void setOutputLines(Size output_lines))
	file.setOutputLines(10);
	TEST_EQUAL(file.getOutputLines() , 10)
END_SECTION

START_SECTION((Size getOutputLines() const))
	TEST_EQUAL(file.getOutputLines() , 10)
END_SECTION

START_SECTION(Size setEnzyme(String enzyme_name))
	TEST_EQUAL(file.setEnzyme("i_dont_exist_enzyme"), 18)
	TEST_EQUAL(file.setEnzyme("Trypsin"), 0)
	TEST_EQUAL(file.getEnzymeNumber() , 14)
END_SECTION

START_SECTION((String getEnzymeName() const))
	TEST_STRING_EQUAL(file.getEnzymeName(), "Trypsin")
END_SECTION

START_SECTION((Size getEnzymeNumber() const))
	TEST_EQUAL(file.getEnzymeNumber() , 14)
END_SECTION

START_SECTION(void setMaxAAPerModPerPeptide(Size max_aa_per_mod_per_peptide))
	file.setMaxAAPerModPerPeptide(4);
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
END_SECTION

START_SECTION((Size getMaxAAPerModPerPeptide() const))
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
END_SECTION

START_SECTION(void setMaxModsPerPeptide(Size max_mods_per_peptide))
	file.setMaxModsPerPeptide(3);
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
END_SECTION

START_SECTION((Size getMaxModsPerPeptide() const))
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
END_SECTION

START_SECTION(void setNucleotideReadingFrame(Size nucleotide_reading_frame))
	file.setNucleotideReadingFrame(0);
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
END_SECTION

START_SECTION((Size getNucleotideReadingFrame() const))
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
END_SECTION

START_SECTION(void setMaxInternalCleavageSites(Size max_internal_cleavage_sites))
	file.setMaxInternalCleavageSites(2);
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
END_SECTION

START_SECTION((Size getMaxInternalCleavageSites() const))
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
END_SECTION

START_SECTION(void setMatchPeakCount(Size match_peak_count))
	file.setMatchPeakCount(5);
	TEST_EQUAL(file.getMatchPeakCount() , 5)
END_SECTION

START_SECTION((Size getMatchPeakCount() const))
	TEST_EQUAL(file.getMatchPeakCount() , 5)
END_SECTION

START_SECTION(void setMatchPeakAllowedError(Size match_peak_allowed_error))
	file.setMatchPeakAllowedError(4);
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
END_SECTION

START_SECTION((Size getMatchPeakAllowedError() const))
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
END_SECTION

START_SECTION(void setShowFragmentIons(bool show_fragments))
	file.setShowFragmentIons(true);
	TEST_EQUAL(file.getShowFragmentIons() , true)
END_SECTION

START_SECTION((bool getShowFragmentIons() const))
	TEST_EQUAL(file.getShowFragmentIons() , true)
END_SECTION

START_SECTION(void setPrintDuplicateReferences(bool print_duplicate_references))
	file.setPrintDuplicateReferences(true);
	TEST_EQUAL(file.getPrintDuplicateReferences() , true)
END_SECTION

START_SECTION((bool getPrintDuplicateReferences() const))
	TEST_EQUAL(file.getPrintDuplicateReferences() , true)
END_SECTION

START_SECTION(void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks))
	file.setRemovePrecursorNearPeaks(true);
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
END_SECTION

START_SECTION((bool getRemovePrecursorNearPeaks() const))
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
END_SECTION

START_SECTION(void setMassTypeParent(bool mass_type_parent))
	file.setMassTypeParent(true);
	TEST_EQUAL(file.getMassTypeParent() , true)
END_SECTION

START_SECTION((bool getMassTypeParent() const))
	TEST_EQUAL(file.getMassTypeParent() , true)
END_SECTION

START_SECTION(void setMassTypeFragment(bool mass_type_fragment))
	file.setMassTypeFragment(true);
	TEST_EQUAL(file.getMassTypeFragment() , true)
END_SECTION

START_SECTION((bool getMassTypeFragment() const))
	TEST_EQUAL(file.getMassTypeFragment() , true)
END_SECTION

START_SECTION(void setNormalizeXcorr(bool normalize_xcorr))
	file.setNormalizeXcorr(true);
	TEST_EQUAL(file.getNormalizeXcorr() , true)
END_SECTION

START_SECTION((bool getNormalizeXcorr() const))
	TEST_EQUAL(file.getNormalizeXcorr() , true)
END_SECTION

START_SECTION(void setResiduesInUpperCase(bool residues_in_upper_case))
	file.setResiduesInUpperCase(true);
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
END_SECTION

START_SECTION((bool getResiduesInUpperCase() const))
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
END_SECTION

START_SECTION(void store(const String& filename))
	String filename;
	NEW_TMP_FILE(filename)

	// test exceptions
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.store(OPENMS_GET_TEST_DATA_PATH("Sequest_unreadable_unwriteable.txt")), "the file `data/Sequest_unreadable_unwriteable.txt' could not be created")

	// test actual program
	file.store(filename);
	TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("SequestInfile_test_template1.txt"));
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
