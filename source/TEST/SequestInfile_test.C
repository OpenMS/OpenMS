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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/SequestInfile.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SequestInfile, "$Id$")

/////////////////////////////////////////////////////////////

SequestInfile* ptr = 0;
CHECK(SequestInfile())
	ptr = new SequestInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SequestInfile())
	delete ptr;
RESULT

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

CHECK(const String getEnzymeInfoAsString())
	TEST_EQUAL(file.getEnzymeInfoAsString(), ss.str())
RESULT

CHECK(void addEnzymeInfo(String& value))
	std::vector< String > e_info;
	e_info.push_back("Z_TestEnzyme");
	e_info.push_back("1");
	e_info.push_back("RMW");
	e_info.push_back("-");
	file.addEnzymeInfo(e_info);
	e_info.clear();
	ss << "17.  Z_TestEnzyme            1     RMW           -" << endl;
	TEST_EQUAL(file.getEnzymeInfoAsString(), ss.str())
RESULT

CHECK(String handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic) throw (Exception::FileNotReadable, Exception::FileNotFound, Exception::ParseError))
	String modification_line = "10.3+,KRLNH,fix:Phosphorylation:+16,C:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
	TEST_EXCEPTION(Exception::FileNotFound, file.handlePTMs(modification_line, "", true))
	
	modification_line = "2H20,KRLNH,fix";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_EQUAL(String(p_e.getMessage()), "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix")
	}

	modification_line = "10.3+";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_EQUAL(String(p_e.getMessage()), "No residues for modification given. Aborting! in: 10.3+")
	}

	modification_line = "10.3+,KRLNH,stat,PTM_0";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_EQUAL(String(p_e.getMessage()), "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0")
	}

	modification_line = "Phosphorylation:Phosphorylation";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_EQUAL(String(p_e.getMessage()), "There's already a modification with this name. Aborting! in: Phosphorylation")
	}

	modification_line = "10.3+,KRLNH,fix:Phosphorylation:+16,C:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm,opt:-16,nterm,fix:17-,cterm_prot:-17,nterm_prot,fix";

	// average masses
	file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", false);

	map< String, vector< String > > modifications;
	modifications["PTM_0"] = vector< String >(3);
	modifications["PTM_0"][0] = "KRLNH";
	modifications["PTM_0"][1] = "10.3";
	modifications["PTM_0"][2] = "FIX";
	modifications["Phosphorylation"] = vector< String >(3);
	modifications["Phosphorylation"][0] = "STYDHCR";
	modifications["Phosphorylation"][1] = "79.9799";
	modifications["Phosphorylation"][2] = "OPT";
	modifications["PTM_2"] = vector< String >(3);
	modifications["PTM_2"][0] = "C";
	modifications["PTM_2"][1] = "16";
	modifications["PTM_2"][2] = "OPT";
// 	modifications["Carbamylation"] = vector< String >(3);
// 	modifications["Carbamylation"][0] = "NTERM";
// 	modifications["Carbamylation"][1] = "43.02474";
// 	modifications["Carbamylation"][2] = "OPT";
	modifications["Methylation"] = vector< String >(3);
	modifications["Methylation"][0] = "CHKNQRILDEST";
	modifications["Methylation"][1] = "14.02658";
	modifications["Methylation"][2] = "OPT";
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
	modifications["PTM_4"] = vector< String >(3);
	modifications["PTM_4"][0] = "CTERM";
	modifications["PTM_4"][1] = "-16";
	modifications["PTM_4"][2] = "OPT";
	modifications["PTM_5"] = vector< String >(3);
	modifications["PTM_5"][0] = "NTERM";
	modifications["PTM_5"][1] = "-16";
	modifications["PTM_5"][2] = "FIX";
	modifications["PTM_6"] = vector< String >(3);
	modifications["PTM_6"][0] = "CTERM_PROT";
	modifications["PTM_6"][1] = "-17";
	modifications["PTM_6"][2] = "OPT";
	modifications["PTM_7"] = vector< String >(3);
	modifications["PTM_7"][0] = "NTERM_PROT";
	modifications["PTM_7"][1] = "-17";
	modifications["PTM_7"][2] = "FIX";

	map< String, vector< String > >::const_iterator result_mod_i = file.getModifications().begin();
	TEST_EQUAL(file.getModifications().size(), modifications.size())
	if ( file.getModifications().size() == modifications.size() )
	{
		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
		{
			TEST_EQUAL(result_mod_i->first, mod_i->first)
			TEST_EQUAL(result_mod_i->second.size(), 3)
			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
			if ( result_mod_i->second.size() == mod_i->second.size() )
			{
				TEST_EQUAL(result_mod_i->second[0], mod_i->second[0])
				TEST_EQUAL(result_mod_i->second[1], mod_i->second[1])
				TEST_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}

	// monoisotopic masses
	file.handlePTMs(modification_line, "TOPP/Sequest_PTMs.xml", true);

	modifications["Phosphorylation"][1] = "79.96635";
// 	modifications["Carbamylation"][1] = "43.00582";
	modifications["Methylation"][1] = "14.01565";

	result_mod_i = file.getModifications().begin();
	TEST_EQUAL(file.getModifications().size(), modifications.size())
	if ( file.getModifications().size() == modifications.size() )
	{
		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
		{
			TEST_EQUAL(result_mod_i->first, mod_i->first)
			TEST_EQUAL(result_mod_i->second.size(), 3)
			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
			if ( result_mod_i->second.size() == mod_i->second.size() )
			{
				TEST_EQUAL(result_mod_i->second[0], mod_i->second[0])
				TEST_EQUAL(result_mod_i->second[1], mod_i->second[1])
				TEST_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}
RESULT

CHECK(void setDatabase(const String& value))
	file.setDatabase("\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta");
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

CHECK(const String& getDatabase())
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

CHECK(void setNeutralLossesForIons(const String& value))
	file.setNeutralLossesForIons("0 1 1");
	TEST_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
RESULT

CHECK(const String& getNeutralLossesForIons())
	TEST_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
RESULT

CHECK(void setIonSeriesWeights(const String& value))
	file.setIonSeriesWeights("0 1.0 0 0 0 0 0 1.0 0");
	TEST_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
RESULT

CHECK(const String& getIonSeriesWeights())
	TEST_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
RESULT

CHECK(void setPartialSequence(const String& value))
	file.setPartialSequence("SEQVEST TEST");
	TEST_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
RESULT

CHECK(const String& getPartialSequence())
	TEST_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
RESULT

CHECK(void setSequenceHeaderFilter(const String& value))
	file.setSequenceHeaderFilter("homo~sapiens !mus musculus");
	TEST_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
RESULT

CHECK(const String& getSequenceHeaderFilter())
	TEST_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
RESULT


CHECK(void setPrecursorMassTolerance(Real value))
	file.setPrecursorMassTolerance(1.3);
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 1.3)
RESULT

CHECK(Real getPrecursorMassTolerance())
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 1.3)
RESULT

CHECK(void setPeakMassTolerance(Real value))
	file.setPeakMassTolerance(0.3);
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 0.3)
RESULT

CHECK(Real getPeakMassTolerance())
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 0.3)
RESULT

CHECK(void setMatchPeakTolerance(Real value))
	file.setMatchPeakTolerance(1.2);
	TEST_REAL_EQUAL(file.getMatchPeakTolerance() , 1.2)
RESULT

CHECK(Real getMatchPeakTolerance())
	TEST_REAL_EQUAL(file.getMatchPeakTolerance() , 1.2)
RESULT

CHECK(void setIonCutoffPercentage(Real value))
	file.setIonCutoffPercentage(0.3);
	TEST_REAL_EQUAL(file.getIonCutoffPercentage() , 0.3)
RESULT

CHECK(Real getIonCutoffPercentage())
	TEST_REAL_EQUAL(file.getIonCutoffPercentage() , 0.3)
RESULT

CHECK(void setProteinMassFilter(const String& protein_mass_filter))
	file.setProteinMassFilter("30.2 0");
	TEST_EQUAL(file.getProteinMassFilter() , "30.2 0")
RESULT

CHECK(Real getProteinMassFilter())
	TEST_EQUAL(file.getProteinMassFilter() , "30.2 0")
RESULT

CHECK(void setPeptideMassUnit(Int value))
	file.setPeptideMassUnit(0);
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK(Int getPeptideMassUnit())
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK(void setOutputLines(Int value))
	file.setOutputLines(10);
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK(Int getOutputLines())
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK(Int setEnzymeNumber(Int value))
	TEST_EQUAL(file.setEnzyme("i_dont_exist_enzyme"), 18)
	TEST_EQUAL(file.setEnzyme("Trypsin"), 0)
	TEST_EQUAL(file.getEnzymeNumber() , 14)
RESULT

CHECK(Int getEnzymeNumber())
	TEST_EQUAL(file.getEnzymeNumber() , 14)
RESULT

CHECK(void setMaxAAPerModPerPeptide(Int value))
	file.setMaxAAPerModPerPeptide(4);
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK(Int getMaxAAPerModPerPeptide())
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK(void setMaxModsPerPeptide(Int value))
	file.setMaxModsPerPeptide(3);
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK(Int getMaxModsPerPeptide())
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK(void setNucleotideReadingFrame(Int value))
	file.setNucleotideReadingFrame(0);
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK(Int getNucleotideReadingFrame())
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK(void setMaxInternalCleavageSites(Int value))
	file.setMaxInternalCleavageSites(2);
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK(Int getMaxInternalCleavageSites())
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK(void setMatchPeakCount(Int value))
	file.setMatchPeakCount(5);
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK(Int getMatchPeakCount())
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK(void setMatchPeakAllowedError(Int value))
	file.setMatchPeakAllowedError(4);
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
RESULT

CHECK(Int getMatchPeakAllowedError())
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
RESULT


CHECK(void setShowFragmentIons(bool value))
	file.setShowFragmentIons(true);
	TEST_EQUAL(file.getShowFragmentIons() , true)
RESULT

CHECK(bool getShowFragmentIons())
	TEST_EQUAL(file.getShowFragmentIons() , true)
RESULT

CHECK(void setPrintDuplicateReferences(bool value))
	file.setPrintDuplicateReferences(true);
	TEST_EQUAL(file.getPrintDuplicateReferences() , true)
RESULT

CHECK(bool getPrintDuplicateReferences())
	TEST_EQUAL(file.getPrintDuplicateReferences() , true)
RESULT

CHECK(void setRemovePrecursorNearPeaks(bool value))
	file.setRemovePrecursorNearPeaks(true);
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
RESULT

CHECK(bool getRemovePrecursorNearPeaks())
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
RESULT

CHECK(void setMassTypeParent(bool value))
	file.setMassTypeParent(true);
	TEST_EQUAL(file.getMassTypeParent() , true)
RESULT

CHECK(bool getMassTypeParent())
	TEST_EQUAL(file.getMassTypeParent() , true)
RESULT

CHECK(void setMassTypeFragment(bool value))
	file.setMassTypeFragment(true);
	TEST_EQUAL(file.getMassTypeFragment() , true)
RESULT

CHECK(bool getMassTypeFragment())
	TEST_EQUAL(file.getMassTypeFragment() , true)
RESULT

CHECK(void setNormalizeXcorr(bool value))
	file.setNormalizeXcorr(true);
	TEST_EQUAL(file.getNormalizeXcorr() , true)
RESULT

CHECK(bool getNormalizeXcorr())
	TEST_EQUAL(file.getNormalizeXcorr() , true)
RESULT

CHECK(void setResiduesInUpperCase(bool value))
	file.setResiduesInUpperCase(true);
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
RESULT

CHECK(bool getResiduesInUpperCase())
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
RESULT

CHECK(void store(const String& filename))
	String filename;
	NEW_TMP_FILE(filename)
	file.store(filename);
	TEST_FILE(filename.c_str(), "data/SequestInfile_test_template1.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
