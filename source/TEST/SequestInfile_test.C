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
CHECK((SequestInfile()))
	ptr = new SequestInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~SequestInfile()))
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

CHECK((const String getEnzymeInfoAsString() const))
	TEST_EQUAL(file.getEnzymeInfoAsString(), ss.str())
RESULT

CHECK((void addEnzymeInfo(std::vector< String > &enzyme_info)))
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

CHECK((void setDatabase(const String &database)))
	file.setDatabase("\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta");
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

CHECK((const String& getDatabase() const))
	TEST_EQUAL(file.getDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta")
RESULT

// CHECK(void setSndDatabase(const String& value))
// 	file.setSndDatabase("\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1");
// 	TEST_EQUAL(file.getSndDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1")
// RESULT
// 
// CHECK(const String& getSndDatabase())
// 	TEST_EQUAL(file.getSndDatabase() , "\\\\bude\\langwisc\\sequest_test\\Analysis.mzXML.fasta1")
// RESULT

CHECK((void setNeutralLossesForIons(const String &neutral_losses_for_ions)))
	file.setNeutralLossesForIons("0 1 1");
	TEST_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
RESULT

CHECK((const String& getNeutralLossesForIons() const))
	TEST_EQUAL(file.getNeutralLossesForIons() , "0 1 1")
RESULT

CHECK((void setIonSeriesWeights(const String &ion_series_weights)))
	file.setIonSeriesWeights("0 1.0 0 0 0 0 0 1.0 0");
	TEST_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
RESULT

CHECK((const String& getIonSeriesWeights() const))
	TEST_EQUAL(file.getIonSeriesWeights() , "0 1.0 0 0 0 0 0 1.0 0")
RESULT

CHECK((void setDynMods(const String &dyn_mods)))
	file.setDynMods("57.8 T 78 YW 0 X");
	TEST_EQUAL(file.getDynMods() , "57.8 T 78 YW 0 X")
RESULT

CHECK((const String& getDynMods() const))
	TEST_EQUAL(file.getDynMods() , "57.8 T 78 YW 0 X")
RESULT

CHECK((void setPartialSequence(const String &partial_sequence)))
	file.setPartialSequence("SEQVEST TEST");
	TEST_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
RESULT

CHECK((const String& getPartialSequence() const))
	TEST_EQUAL(file.getPartialSequence() , "SEQVEST TEST")
RESULT

CHECK((void setSequenceHeaderFilter(const String &sequence_header_filter)))
	file.setSequenceHeaderFilter("homo~sapiens !mus musculus");
	TEST_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
RESULT

CHECK((const String& getSequenceHeaderFilter() const))
	TEST_EQUAL(file.getSequenceHeaderFilter() , "homo~sapiens !mus musculus")
RESULT


CHECK((void setPrecursorMassTolerance(Real precursor_mass_tolerance)))
	file.setPrecursorMassTolerance(1.3);
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 1.3)
RESULT

CHECK((Real getPrecursorMassTolerance() const))
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 1.3)
RESULT

CHECK((void setPeakMassTolerance(Real peak_mass_tolerance)))
	file.setPeakMassTolerance(0.3);
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 0.3)
RESULT

CHECK((Real getPeakMassTolerance() const))
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 0.3)
RESULT

CHECK((void setMatchPeakTolerance(Real match_peak_tolerance)))
	file.setMatchPeakTolerance(1.2);
	TEST_REAL_EQUAL(file.getMatchPeakTolerance() , 1.2)
RESULT

CHECK((Real getMatchPeakTolerance() const))
	TEST_REAL_EQUAL(file.getMatchPeakTolerance() , 1.2)
RESULT

CHECK((void setIonCutoffPercentage(Real cutoff_percentage)))
	file.setIonCutoffPercentage(0.3);
	TEST_REAL_EQUAL(file.getIonCutoffPercentage() , 0.3)
RESULT

CHECK((Real getIonCutoffPercentage() const))
	TEST_REAL_EQUAL(file.getIonCutoffPercentage() , 0.3)
RESULT

CHECK((void setProteinMassFilter(const String& protein_mass_filter)))
	file.setProteinMassFilter("30.2 0");
	TEST_EQUAL(file.getProteinMassFilter() , "30.2 0")
RESULT

CHECK((const String& getProteinMassFilter() const))
	TEST_EQUAL(file.getProteinMassFilter() , "30.2 0")
RESULT

CHECK((void setDynNTermMod(Real dyn_n_term_mod)))
	file.setDynNTermMod(32.4);
	TEST_REAL_EQUAL(file.getDynNTermMod() , 32.4)
RESULT

CHECK((Real getDynNTermMod() const))
	TEST_REAL_EQUAL(file.getDynNTermMod() , 32.4)
RESULT

CHECK((void setDynCTermMod(Real dyn_c_term_mod)))
	file.setDynCTermMod(32.0);
	TEST_REAL_EQUAL(file.getDynCTermMod() , 32.0)
RESULT

CHECK((Real getDynCTermMod() const))
	TEST_REAL_EQUAL(file.getDynCTermMod() , 32.0)
RESULT

CHECK((void setStatNTermMod(Real stat_n_term_mod)))
	file.setStatNTermMod(32.1);
	TEST_REAL_EQUAL(file.getStatNTermMod() , 32.1)
RESULT

CHECK((Real getStatNTermMod() const))
	TEST_REAL_EQUAL(file.getStatNTermMod() , 32.1)
RESULT

CHECK((void setStatCTermMod(Real stat_c_term_mod)))
	file.setStatCTermMod(32.2);
	TEST_REAL_EQUAL(file.getStatCTermMod() , 32.2)
RESULT

CHECK((Real getStatCTermMod() const))
	TEST_REAL_EQUAL(file.getStatCTermMod() , 32.2)
RESULT

CHECK((void setStatNTermProtMod(Real stat_n_term_prot_mod)))
	file.setStatNTermProtMod(32.3);
	TEST_REAL_EQUAL(file.getStatNTermProtMod() , 32.3)
RESULT

CHECK((Real getStatNTermProtMod() const))
	TEST_REAL_EQUAL(file.getStatNTermProtMod() , 32.3)
RESULT

CHECK((void setStatCTermProtMod(Real stat_c_term_prot_mod)))
	file.setStatCTermProtMod(32.5);
	TEST_REAL_EQUAL(file.getStatCTermProtMod() , 32.5)
RESULT

CHECK((Real getStatCTermProtMod() const))
	TEST_REAL_EQUAL(file.getStatCTermProtMod() , 32.5)
RESULT


CHECK((void setPeptideMassUnit(Int peptide_mass_unit)))
	file.setPeptideMassUnit(0);
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK((Int getPeptideMassUnit() const))
	TEST_EQUAL(file.getPeptideMassUnit() , 0)
RESULT

CHECK((void setOutputLines(Int output_lines)))
	file.setOutputLines(10);
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK((Int getOutputLines() const))
	TEST_EQUAL(file.getOutputLines() , 10)
RESULT

CHECK((Int getEnzymeNumber() const))
	TEST_EQUAL(file.setEnzyme("i_dont_exist_enzyme"), 18)
	TEST_EQUAL(file.setEnzyme("Trypsin"), 0)
	TEST_EQUAL(file.getEnzymeNumber() , 14)
RESULT

CHECK((Int getEnzymeNumber() const))
	TEST_EQUAL(file.getEnzymeNumber() , 14)
RESULT

CHECK((void setMaxAAPerModPerPeptide(Int max_aa_per_mod_per_peptide)))
	file.setMaxAAPerModPerPeptide(4);
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK((Int getMaxAAPerModPerPeptide() const))
	TEST_EQUAL(file.getMaxAAPerModPerPeptide() , 4)
RESULT

CHECK((void setMaxModsPerPeptide(Int max_mods_per_peptide)))
	file.setMaxModsPerPeptide(3);
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK((Int getMaxModsPerPeptide() const))
	TEST_EQUAL(file.getMaxModsPerPeptide() , 3)
RESULT

CHECK((void setNucleotideReadingFrame(Int nucleotide_reading_frame)))
	file.setNucleotideReadingFrame(0);
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK((Int getNucleotideReadingFrame() const))
	TEST_EQUAL(file.getNucleotideReadingFrame() , 0)
RESULT

CHECK((void setMaxInternalCleavageSites(Int max_internal_cleavage_sites)))
	file.setMaxInternalCleavageSites(2);
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK((Int getMaxInternalCleavageSites() const))
	TEST_EQUAL(file.getMaxInternalCleavageSites() , 2)
RESULT

CHECK((void setMatchPeakCount(Int match_peak_count)))
	file.setMatchPeakCount(5);
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK((Int getMatchPeakCount() const))
	TEST_EQUAL(file.getMatchPeakCount() , 5)
RESULT

CHECK((void setMatchPeakAllowedError(Int match_peak_allowed_error)))
	file.setMatchPeakAllowedError(4);
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
RESULT

CHECK((Int getMatchPeakAllowedError() const))
	TEST_EQUAL(file.getMatchPeakAllowedError() , 4)
RESULT


CHECK((void setShowFragmentIons(bool show_fragments)))
	file.setShowFragmentIons(true);
	TEST_EQUAL(file.getShowFragmentIons() , true)
RESULT

CHECK((bool getShowFragmentIons() const))
	TEST_EQUAL(file.getShowFragmentIons() , true)
RESULT

CHECK((void setPrintDuplicateReferences(bool print_duplicate_references)))
	file.setPrintDuplicateReferences(true);
	TEST_EQUAL(file.getPrintDuplicateReferences() , true)
RESULT

CHECK((bool getPrintDuplicateReferences() const))
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

CHECK((void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks)))
	file.setRemovePrecursorNearPeaks(true);
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
RESULT

CHECK((bool getRemovePrecursorNearPeaks() const))
	TEST_EQUAL(file.getRemovePrecursorNearPeaks() , true)
RESULT

CHECK((void setMassTypeParent(bool mass_type_parent)))
	file.setMassTypeParent(true);
	TEST_EQUAL(file.getMassTypeParent() , true)
RESULT

CHECK((bool getMassTypeParent() const))
	TEST_EQUAL(file.getMassTypeParent() , true)
RESULT

CHECK((void setMassTypeFragment(bool mass_type_fragment)))
	file.setMassTypeFragment(true);
	TEST_EQUAL(file.getMassTypeFragment() , true)
RESULT

CHECK((bool getMassTypeFragment() const))
	TEST_EQUAL(file.getMassTypeFragment() , true)
RESULT

CHECK((void setNormalizeXcorr(bool normalize_xcorr)))
	file.setNormalizeXcorr(true);
	TEST_EQUAL(file.getNormalizeXcorr() , true)
RESULT

CHECK((bool getNormalizeXcorr() const))
	TEST_EQUAL(file.getNormalizeXcorr() , true)
RESULT

CHECK((void setResiduesInUpperCase(bool residues_in_upper_case)))
	file.setResiduesInUpperCase(true);
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
RESULT

CHECK((bool getResiduesInUpperCase() const))
	TEST_EQUAL(file.getResiduesInUpperCase() , true)
RESULT

vector< Real > masses(4, 0.1);
String aas = "GASPVTCLIXNOBDQKZEMHFRYW";
CHECK((char setStatMod(String amino_acid, Real mass)))
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

CHECK((const std::map< char, Real >& getStatMods() const))
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

CHECK((void store(const String &filename) throw (Exception::UnableToCreateFile)))
	String filename;
	NEW_TMP_FILE(filename)
	file.store(filename);
	TEST_FILE(filename.c_str(), "data/SequestInfile_test_template1.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
