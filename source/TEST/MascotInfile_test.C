// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MascotInfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(MascotInfile, "$Id$")

/////////////////////////////////////////////////////////////

PeakSpectrum spec;
Peak1D tmp;
vector<Int> charges;
charges.push_back(2);
for (Size i=1;i<10;i+=1)
{
	tmp.setPosition(DPosition<1>(i));
	tmp.setIntensity(i * i);
	spec.push_back(tmp);	
}

MascotInfile* ptr = 0;
START_SECTION((MascotInfile()))
	ptr = new MascotInfile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MascotInfile()))
	delete ptr;
END_SECTION

MascotInfile file;
file.setCharges(charges);

START_SECTION((void setBoundary(const String &boundary)))
	file.setBoundary("ABCDEFGHIJKMNOPQRSTUVWXYZ");
	TEST_EQUAL(file.getBoundary() , "ABCDEFGHIJKMNOPQRSTUVWXYZ")
END_SECTION

START_SECTION((const String& getBoundary()))
	TEST_EQUAL(file.getBoundary() , "ABCDEFGHIJKMNOPQRSTUVWXYZ")
END_SECTION

START_SECTION((void store(const String &filename, const PeakSpectrum& spec, DoubleReal mz, DoubleReal retention_time, String search_title)))

	// here a fixed name has to be used as it has to be in the template
	file.store("MascotInfile_test.txt", spec, 1998.0f, 25.379, "TestTitle");
	TEST_FILE_EQUAL("MascotInfile_test.txt", OPENMS_GET_TEST_DATA_PATH("MascotInfile_test_template1.txt"));
	remove("MascotInfile_test.txt");
END_SECTION

START_SECTION((void setDB(const String &db)))
	file.setDB("DB_TEST");
	TEST_EQUAL(file.getDB() , "DB_TEST")
END_SECTION

START_SECTION((const String& getDB()))
	TEST_EQUAL(file.getDB() , "DB_TEST")
END_SECTION

START_SECTION((void setSearchType(const String &search_type)))
	file.setSearchType("SearchType_TEST");
	TEST_EQUAL(file.getSearchType() , "SearchType_TEST")
END_SECTION

START_SECTION((const String& getSearchType()))
	TEST_EQUAL(file.getSearchType() , "SearchType_TEST")
END_SECTION

START_SECTION((void setHits(const String &hits)))
	file.setHits("Hits_TEST");
	TEST_EQUAL(file.getHits() , "Hits_TEST")
END_SECTION

START_SECTION((const String& getHits()))
	TEST_EQUAL(file.getHits() , "Hits_TEST")
END_SECTION

START_SECTION((void setCleavage(const String &cleavage)))
	file.setCleavage("Cleavage_TEST");
	TEST_EQUAL(file.getCleavage() , "Cleavage_TEST")
END_SECTION

START_SECTION((const String& getCleavage()))
	TEST_EQUAL(file.getCleavage() , "Cleavage_TEST")
END_SECTION

START_SECTION((void setMassType(const String &mass_type)))
	file.setMassType("MassType_TEST");
	TEST_EQUAL(file.getMassType() , "MassType_TEST")
END_SECTION

START_SECTION((const String& getMassType()))
	TEST_EQUAL(file.getMassType() , "MassType_TEST")
END_SECTION

START_SECTION((void setInstrument(const String &instrument)))
	file.setInstrument("Instrument_TEST");
	TEST_EQUAL(file.getInstrument() , "Instrument_TEST")
END_SECTION

START_SECTION((const String& getInstrument()))
	TEST_EQUAL(file.getInstrument() , "Instrument_TEST")
END_SECTION

START_SECTION((void setMissedCleavages(UInt missed_cleavages)))
	file.setMissedCleavages(4711);
	TEST_EQUAL(file.getMissedCleavages() , 4711)
END_SECTION

START_SECTION((UInt getMissedCleavages()))
	TEST_EQUAL(file.getMissedCleavages() , 4711)
END_SECTION

START_SECTION((void setPrecursorMassTolerance(Real precursor_mass_tolerance)))
	file.setPrecursorMassTolerance(4711.1f);
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance() , 4711.1f)
END_SECTION

START_SECTION((Real getPrecursorMassTolerance()))
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance() , 4711.1f)
END_SECTION

START_SECTION((void setPeakMassTolerance(Real ion_mass_tolerance)))
	file.setPeakMassTolerance(4711.2f);
	TEST_REAL_SIMILAR(file.getPeakMassTolerance() , 4711.2f)
END_SECTION

START_SECTION((Real getPeakMassTolerance()))
	TEST_REAL_SIMILAR(file.getPeakMassTolerance() , 4711.2f)
END_SECTION

START_SECTION((void setTaxonomy(const String &taxonomy)))
	file.setTaxonomy("Taxonomy_TEST");
	TEST_EQUAL(file.getTaxonomy() , "Taxonomy_TEST")
END_SECTION

START_SECTION((const String& getTaxonomy()))
	TEST_EQUAL(file.getTaxonomy() , "Taxonomy_TEST")
END_SECTION

START_SECTION((void setFormVersion(const String &form_version)))
	file.setFormVersion("FormVersion_TEST");
	TEST_EQUAL(file.getFormVersion() , "FormVersion_TEST")
END_SECTION

START_SECTION((const String& getFormVersion()))
	TEST_EQUAL(file.getFormVersion() , "FormVersion_TEST")
END_SECTION

vector<String> mods;
mods.push_back("Modifiactions_TEST_1");
mods.push_back("Modifiactions_TEST_2");
vector<String> vmods;
vmods.push_back("Variable_Modifiactions_TEST_1");
vmods.push_back("Variable_Modifiactions_TEST_2");

START_SECTION((void setModifications(const std::vector<String>& mods)))
	file.setModifications(mods);
	TEST_EQUAL(file.getModifications() == mods, true)
END_SECTION

START_SECTION((const std::vector<String>& getModifications()))
	TEST_EQUAL(file.getModifications() == mods, true)
END_SECTION

START_SECTION((void setVariableModifications(const std::vector<String>& mods)))
	file.setVariableModifications(vmods);
	TEST_EQUAL(file.getVariableModifications() == vmods, true)
END_SECTION

START_SECTION((const std::vector<String>& getVariableModifications()))
	TEST_EQUAL(file.getVariableModifications() == vmods, true)
END_SECTION

START_SECTION([EXTRA] void store(const std::string& filename, const PeakSpectrum& spec, double mz, double retention_time, std::string search_title))
	// here a fixed name has to be used as it has to be in the tamplate
	file.store("MascotInfile_test.txt", spec, 1998.0f, 25.379, "TestTitle");
	TEST_FILE_EQUAL("MascotInfile_test.txt", OPENMS_GET_TEST_DATA_PATH("MascotInfile_test_template2.txt"));
	remove("MascotInfile_test.txt");
END_SECTION

START_SECTION((void setCharges(std::vector<Int>& charges)))
	charges.push_back(3);
	charges.push_back(1);
	file.setCharges(charges);
	TEST_EQUAL(file.getCharges(), "1+, 2+ and 3+")
END_SECTION

START_SECTION((const String& getCharges()))
	TEST_EQUAL(file.getCharges(), "1+, 2+ and 3+")
END_SECTION

START_SECTION((void store(const String &filename, const MSExperiment<> &experiment, String search_title)))
	MSExperiment<> exp;
	MSExperiment<>::SpectrumType spec;
	MSExperiment<>::PeakType peak;

	// first spectrum (MS)
	spec.setRT(11.1);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 5;
	peak.setIntensity(47.11f);
	spec.push_back(peak);
	peak.getPosition()[0] = 10;
	peak.setIntensity(48.11f);
	spec.push_back(peak);
	peak.getPosition()[0] = 15;
	spec.push_back(peak);
	exp.push_back(spec);

	// second spectrum (MS/MS)
	spec.clear(true);
	spec.setRT(11.5);
	spec.getPrecursors().resize(1);
	spec.getPrecursors()[0].setMZ(11.4);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 6;
	spec.push_back(peak);
	peak.getPosition()[0] = 11;
	spec.push_back(peak);
	exp.push_back(spec);	

	// third spectrum (MS)
	spec.clear(true);
	spec.setRT(12.2);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 20;
	spec.push_back(peak);
	peak.getPosition()[0] = 25;
	spec.push_back(peak);
	exp.push_back(spec);	

	// forth spectrum (MS/MS)
	spec.clear(true);
	spec.setRT(12.5);
	spec.getPrecursors().resize(1);
	spec.getPrecursors()[0].setMZ(21.4);
	spec.setMSLevel(2);
	peak.getPosition()[0] = 21;
	spec.push_back(peak);
	peak.getPosition()[0] = 26;
	spec.push_back(peak);
	peak.getPosition()[0] = 31;
	spec.push_back(peak);
	exp.push_back(spec);	

	file.store("MascotInfile_test.txt", exp, "Experiment");
	TEST_FILE_EQUAL("MascotInfile_test.txt", OPENMS_GET_TEST_DATA_PATH("MascotInfile_test_template3.txt"));
	remove("MascotInfile_test.txt");
END_SECTION

START_SECTION(template <typename MapType> void load(const String &filename, MapType &exp))
	MascotInfile infile;
	MSExperiment<> experiment;
	MSExperiment<>::SpectrumType spectrum;
	
	infile.load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), experiment);
	spectrum = experiment[0];
	TEST_REAL_SIMILAR(spectrum.getRT(), 25.379)
	TEST_EQUAL(spectrum.getPrecursors().size(),1)
	TEST_REAL_SIMILAR(spectrum.getPrecursors()[0].getMZ(), 1998) 
	TEST_EQUAL(spectrum.getMetaValue("TITLE"), "Testtitle");

	TEST_REAL_SIMILAR(spectrum[0].getIntensity(), 1)
	TEST_REAL_SIMILAR(spectrum[0].getMZ(), 1)
	TEST_REAL_SIMILAR(spectrum[1].getIntensity(), 4)
	TEST_REAL_SIMILAR(spectrum[1].getMZ(), 2)
	TEST_REAL_SIMILAR(spectrum[2].getIntensity(), 9)
	TEST_REAL_SIMILAR(spectrum[2].getMZ(), 3)
	TEST_REAL_SIMILAR(spectrum[3].getIntensity(), 16)
	TEST_REAL_SIMILAR(spectrum[3].getMZ(), 4)
	TEST_REAL_SIMILAR(spectrum[4].getIntensity(), 25)
	TEST_REAL_SIMILAR(spectrum[4].getMZ(), 5)
	TEST_REAL_SIMILAR(spectrum[5].getIntensity(), 36)
	TEST_REAL_SIMILAR(spectrum[5].getMZ(), 6)
	TEST_REAL_SIMILAR(spectrum[6].getIntensity(), 49)
	TEST_REAL_SIMILAR(spectrum[6].getMZ(), 7)
	TEST_REAL_SIMILAR(spectrum[7].getIntensity(), 64)
	TEST_REAL_SIMILAR(spectrum[7].getMZ(), 8)
	TEST_REAL_SIMILAR(spectrum[8].getIntensity(), 81)
	TEST_REAL_SIMILAR(spectrum[8].getMZ(), 9)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
