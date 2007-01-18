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
ss << "0.  No_Enzyme              0      -                        -" << endl;
ss << "1.  Trypsin_Strict         1      KR                       -" << endl;
ss << "2.  Trypsin                1      KRLNH                    -" << endl;
ss << "3.  Chymotrypsin           1      FWYL                     -" << endl;
ss << "4.  Chymotrypsin_WYF       1      FWY                      -" << endl;
ss << "5.  Clostripain            1      R                        -" << endl;
ss << "6.  Cyanogen_Bromide       1      M                        -" << endl;
ss << "7.  IodosoBenzoate         1      W                        -" << endl;
ss << "8.  Proline_Endopept       1      P                        -" << endl;
ss << "9.  GluC                   1      E                        -" << endl;
ss << "10. GluC_ED                1      ED                       -" << endl;
ss << "11. LysC                   1      K                        -" << endl;
ss << "12. AspN                   0      D                        -" << endl;
ss << "13. AspN_DE                0      DE                       -" << endl;
ss << "14. Elastase               1      ALIV                     P" << endl;
ss << "15. Elastase/Tryp/Chymo    1      ALIVKRWFY                P" << endl;
ss << "16. Trypsin/Chymo          1      KRLFWYN                  -" << endl;

CHECK(const String& getEnzymeInfo())
	TEST_EQUAL(file.getEnzymeInfo(), ss.str())
RESULT

CHECK(void addEnzymeInfo(String& value))
	std::vector< String > e_info;
	e_info.push_back("TestEnzyme");
	e_info.push_back("1");
	e_info.push_back("RMW");
	e_info.push_back("-");
	file.addEnzymeInfo(e_info);
	e_info.clear();
	ss << "17. TestEnzyme             1      RMW                      -" << endl;
	TEST_EQUAL(file.getEnzymeInfo(), ss.str())
RESULT

CHECK(void setDatabase(const String& value))
	file.setDatabase("\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta");
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

CHECK(const String& getDatabase())
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

CHECK(void setSndDatabase(const String& value))
	file.setSndDatabase("\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1");
	TEST_EQUAL(file.getSndDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1")
RESULT

CHECK(const String& getSndDatabase())
	TEST_EQUAL(file.getSndDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1")
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

CHECK(void setDynMods(const String& value))
	file.setDynMods("57.8 T 78 YW 0 X");
	TEST_EQUAL(file.getDynMods() , "57.8 T 78 YW 0 X")
RESULT

CHECK(const String& getDynMods())
	TEST_EQUAL(file.getDynMods() , "57.8 T 78 YW 0 X")
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


CHECK(void setPeptideMassTolerance(Real value))
	file.setPeptideMassTolerance(1.3);
	TEST_REAL_EQUAL(file.getPeptideMassTolerance() , 1.3)
RESULT

CHECK(Real getPeptideMassTolerance())
	TEST_REAL_EQUAL(file.getPeptideMassTolerance() , 1.3)
RESULT

CHECK(void setFragmentIonTolerance(Real value))
	file.setFragmentIonTolerance(0.3);
	TEST_REAL_EQUAL(file.getFragmentIonTolerance() , 0.3)
RESULT

CHECK(Real getFragmentIonTolerance())
	TEST_REAL_EQUAL(file.getFragmentIonTolerance() , 0.3)
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

CHECK(void setDynNTermMod(Real value))
	file.setDynNTermMod(32.4);
	TEST_REAL_EQUAL(file.getDynNTermMod() , 32.4)
RESULT

CHECK(Real getDynNTermMod())
	TEST_REAL_EQUAL(file.getDynNTermMod() , 32.4)
RESULT

CHECK(void setDynCTermMod(Real value))
	file.setDynCTermMod(32.0);
	TEST_REAL_EQUAL(file.getDynCTermMod() , 32.0)
RESULT

CHECK(Real getDynCTermMod())
	TEST_REAL_EQUAL(file.getDynCTermMod() , 32.0)
RESULT

CHECK(void setStatNTermMod(Real value))
	file.setStatNTermMod(32.1);
	TEST_REAL_EQUAL(file.getStatNTermMod() , 32.1)
RESULT

CHECK(Real getStatNTermMod())
	TEST_REAL_EQUAL(file.getStatNTermMod() , 32.1)
RESULT

CHECK(void setStatCTermMod(Real value))
	file.setStatCTermMod(32.2);
	TEST_REAL_EQUAL(file.getStatCTermMod() , 32.2)
RESULT

CHECK(Real getStatCTermMod())
	TEST_REAL_EQUAL(file.getStatCTermMod() , 32.2)
RESULT

CHECK(void setStatNTermProtMod(Real value))
	file.setStatNTermProtMod(32.3);
	TEST_REAL_EQUAL(file.getStatNTermProtMod() , 32.3)
RESULT

CHECK(Real getStatNTermProtMod())
	TEST_REAL_EQUAL(file.getStatNTermProtMod() , 32.3)
RESULT

CHECK(void setStatCTermProtMod(Real value))
	file.setStatCTermProtMod(32.5);
	TEST_REAL_EQUAL(file.getStatCTermProtMod() , 32.5)
RESULT

CHECK(Real getStatCTermProtMod())
	TEST_REAL_EQUAL(file.getStatCTermProtMod() , 32.5)
RESULT


CHECK(void setPeptideMassUnit(SignedInt value))
	file.setPeptideMassUnit(0);
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK(SignedInt getPeptideMassUnit())
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK(void setOutputLines(SignedInt value))
	file.setOutputLines(10);
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK(SignedInt getOutputLines())
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK(SignedInt setEnzymeNumber(SignedInt value))
	TEST_EQUAL(file.setEnzymeNumber(-1), 17)
	TEST_EQUAL(file.setEnzymeNumber(20), 17)
	TEST_EQUAL(file.setEnzymeNumber(1), 0)
	file.setEnzymeNumber(1);
	TEST_EQUAL(file.getEnzymeNumber() , 1)
RESULT

CHECK(SignedInt getEnzymeNumber())
	TEST_EQUAL(file.getEnzymeNumber() , 1)
RESULT

CHECK(void setMaxAAPerModPerPeptide(SignedInt value))
	file.setMaxAAPerModPerPeptide(4);
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK(SignedInt getMaxAAPerModPerPeptide())
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK(void setMaxModsPerPeptide(SignedInt value))
	file.setMaxModsPerPeptide(3);
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK(SignedInt getMaxModsPerPeptide())
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK(void setNucleotideReadingFrame(SignedInt value))
	file.setNucleotideReadingFrame(0);
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK(SignedInt getNucleotideReadingFrame())
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK(void setMaxInternalCleavageSites(SignedInt value))
	file.setMaxInternalCleavageSites(2);
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK(SignedInt getMaxInternalCleavageSites())
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK(void setMatchPeakCount(SignedInt value))
	file.setMatchPeakCount(5);
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK(SignedInt getMatchPeakCount())
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK(void setMatchPeakAllowedError(SignedInt value))
	file.setMatchPeakAllowedError(4);
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
RESULT

CHECK(SignedInt getMatchPeakAllowedError())
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

// CHECK(void setUsePhosphoFragmentation(bool value))
// 	file.setUsePhosphoFragmentation(true);
// 	TEST_EQUAL(file.getUsePhosphoFragmentation() , true)
// RESULT
// 
// CHECK(bool getUsePhosphoFragmentation())
// 	TEST_EQUAL(file.getUsePhosphoFragmentation() , true)
// RESULT

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

vector< Real > masses(4, 0.1);
String aas = "GASPVTCLIXNOBDQKZEMHFRYW";
CHECK(char setStatMod(String& amino_acids, Real mass));
	Real mass = 0.1;
	String string_buffer = aas.substr(0,4);
	TEST_EQUAL(file.setStatMod(string_buffer, mass), 0)
	for ( String::const_iterator s_i = aas.begin()+4; s_i != aas.end(); ++s_i )
	{
		masses.push_back(++mass);
		string_buffer = *s_i;
		TEST_EQUAL(file.setStatMod(string_buffer, mass), 0)
	}
	
	map< char, Real > mod_map = file.getStatMods();
	size_t pos;
	for ( map< char, Real >::const_iterator m_i = mod_map.begin(); m_i != mod_map.end(); ++m_i )
	{
		pos = aas.find(m_i->first);
		TEST_NOT_EQUAL(pos, string::npos)
		if ( pos != string::npos ) TEST_REAL_EQUAL(masses[pos], m_i->second)
	}
RESULT

CHECK((const map< char, Real >& getStatMods()))
	map< char, Real > mod_map = file.getStatMods();
	size_t pos;
	for ( map< char, Real >::const_iterator m_i = mod_map.begin(); m_i != mod_map.end(); ++m_i )
	{
		pos = aas.find(m_i->first);
		TEST_NOT_EQUAL(pos, string::npos)
		if ( pos != string::npos ) TEST_REAL_EQUAL(masses[pos], m_i->second)
	}
RESULT
masses.clear();

CHECK(void store(const String& filename))
	String filename;
	NEW_TMP_FILE(filename)
	file.store(filename);
	TEST_FILE(filename.c_str(), "data/SequestInfile_test_template1.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
