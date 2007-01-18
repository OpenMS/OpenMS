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
// $Maintainer: Nico Pfeifer $
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

//DPeakArray (dummy for spectrum)
DPeakArray<1> spec;
DPeak<1> tmp;
vector<SignedInt> charges;
charges.push_back(2);
for (UnsignedInt i=1;i<10;i+=1)
{
	tmp.setPosition(DPosition<1>(i));
	tmp.setIntensity(i*i);
	spec.push_back(tmp);	
}

MascotInfile* ptr = 0;
CHECK((MascotInfile()))
	ptr = new MascotInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MascotInfile()))
	delete ptr;
RESULT

MascotInfile file;
file.setCharges(charges);

CHECK((void setBoundary(const std::string& boundary)))
	file.setBoundary("ABCDEFGHIJKMNOPQRSTUVWXYZ");
	TEST_EQUAL(file.getBoundary() , "ABCDEFGHIJKMNOPQRSTUVWXYZ")
RESULT

CHECK((const std::string& getBoundary()))
	TEST_EQUAL(file.getBoundary() , "ABCDEFGHIJKMNOPQRSTUVWXYZ")
RESULT

CHECK((void store(const std::string& filename, const DPeakArray<1>& spec, double mz, double retention_time, std::string search_title)))

	// here a fixed name has to be used as it has to be in the template
	file.store("MascotInfile_test.txt", spec, 1998.0f, 25.379, "TestTitle");
	TEST_FILE("MascotInfile_test.txt", "data/MascotInfile_test_template1.txt");
	remove("MascotInfile_test.txt");
RESULT

CHECK((void setDB(const std::string& db)))
	file.setDB("DB_TEST");
	TEST_EQUAL(file.getDB() , "DB_TEST")
RESULT

CHECK((const std::string& getDB()))
	TEST_EQUAL(file.getDB() , "DB_TEST")
RESULT

CHECK((void setSearchType(const std::string& search_type)))
	file.setSearchType("SearchType_TEST");
	TEST_EQUAL(file.getSearchType() , "SearchType_TEST")
RESULT

CHECK((const std::string& getSearchType()))
	TEST_EQUAL(file.getSearchType() , "SearchType_TEST")
RESULT

CHECK((void setHits(const std::string& hits)))
	file.setHits("Hits_TEST");
	TEST_EQUAL(file.getHits() , "Hits_TEST")
RESULT

CHECK((const std::string& getHits()))
	TEST_EQUAL(file.getHits() , "Hits_TEST")
RESULT

CHECK((void setCleavage(const std::string& cleavage)))
	file.setCleavage("Cleavage_TEST");
	TEST_EQUAL(file.getCleavage() , "Cleavage_TEST")
RESULT

CHECK((const std::string& getCleavage()))
	TEST_EQUAL(file.getCleavage() , "Cleavage_TEST")
RESULT

CHECK((void setMassType(const std::string& mass_type)))
	file.setMassType("MassType_TEST");
	TEST_EQUAL(file.getMassType() , "MassType_TEST")
RESULT

CHECK((const std::string& getMassType()))
	TEST_EQUAL(file.getMassType() , "MassType_TEST")
RESULT

CHECK((void setInstrument(const std::string& instrument)))
	file.setInstrument("Instrument_TEST");
	TEST_EQUAL(file.getInstrument() , "Instrument_TEST")
RESULT

CHECK((const std::string& getInstrument()))
	TEST_EQUAL(file.getInstrument() , "Instrument_TEST")
RESULT

CHECK((void setMissedCleavages(UnsignedInt missed_cleavages)))
	file.setMissedCleavages(4711);
	TEST_EQUAL(file.getMissedCleavages() , 4711)
RESULT

CHECK((UnsignedInt getMissedCleavages()))
	TEST_EQUAL(file.getMissedCleavages() , 4711)
RESULT

CHECK((void setPrecursorMassTolerance(float precursor_mass_tolerance)))
	file.setPrecursorMassTolerance(4711.1f);
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 4711.1f)
RESULT

CHECK((float getPrecursorMassTolerance()))
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance() , 4711.1f)
RESULT

CHECK((void setPeakMassTolerance(float ion_mass_tolerance)))
	file.setPeakMassTolerance(4711.2f);
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 4711.2f)
RESULT

CHECK((float getPeakMassTolerance()))
	TEST_REAL_EQUAL(file.getPeakMassTolerance() , 4711.2f)
RESULT

CHECK((void setTaxonomy(const std::string& taxonomy)))
	file.setTaxonomy("Taxonomy_TEST");
	TEST_EQUAL(file.getTaxonomy() , "Taxonomy_TEST")
RESULT

CHECK((const std::string& getTaxonomy()))
	TEST_EQUAL(file.getTaxonomy() , "Taxonomy_TEST")
RESULT

CHECK((void setFormVersion(const std::string& form_version)))
	file.setFormVersion("FormVersion_TEST");
	TEST_EQUAL(file.getFormVersion() , "FormVersion_TEST")
RESULT

CHECK((const std::string& getFormVersion()))
	TEST_EQUAL(file.getFormVersion() , "FormVersion_TEST")
RESULT

vector<String> mods;
mods.push_back("Modifiactions_TEST_1");
mods.push_back("Modifiactions_TEST_2");
vector<String> vmods;
vmods.push_back("Variable_Modifiactions_TEST_1");
vmods.push_back("Variable_Modifiactions_TEST_2");

CHECK((void setModifications(const std::vector<String>& mods)))
	file.setModifications(mods);
	TEST_EQUAL(file.getModifications() == mods, true)
RESULT

CHECK((const std::vector<String>& getModifications()))
	TEST_EQUAL(file.getModifications() == mods, true)
RESULT

CHECK((void setVariableModifications(const std::vector<String>& mods)))
	file.setVariableModifications(vmods);
	TEST_EQUAL(file.getVariableModifications() == vmods, true)
RESULT

CHECK((const std::vector<String>& getVariableModifications()))
	TEST_EQUAL(file.getVariableModifications() == vmods, true)
RESULT

CHECK([EXTRA] void store(const std::string& filename, const DPeakArray<1>& spec, double mz, double retention_time, std::string search_title))
	// here a fixed name has to be used as it has to be in the tamplate
	file.store("MascotInfile_test.txt", spec, 1998.0f, 25.379, "TestTitle");
	TEST_FILE("MascotInfile_test.txt", "data/MascotInfile_test_template2.txt");
	remove("MascotInfile_test.txt");
RESULT

CHECK((void setCharges(std::vector<SignedInt>& charges)))
	charges.push_back(3);
	charges.push_back(1);
	file.setCharges(charges);
	TEST_EQUAL(file.getCharges(), "1+, 2+ and 3+")
RESULT

CHECK((const std::string& getCharges()))
	TEST_EQUAL(file.getCharges(), "1+, 2+ and 3+")
RESULT

CHECK((void store(const std::string& filename, const MSExperiment< DPeak<1> >& experiment, std::string search_title)))
	MSExperiment<> exp;
	MSExperiment<>::SpectrumType spec;
	MSExperiment<>::PeakType peak;

	// first spectrum (MS)
	spec.setRetentionTime(11.1);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 5;
	peak.setIntensity(47.11);
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 10;
	peak.setIntensity(48.11);
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 15;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);

	// second spectrum (MS/MS)
	spec.getContainer().clear();
	spec.setRetentionTime(11.5);
	spec.getPrecursorPeak().getPosition()[0] = 11.4;
	spec.setMSLevel(2);
	peak.getPosition()[0] = 6;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 11;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// third spectrum (MS)
	spec.getContainer().clear();
	spec.setRetentionTime(12.2);
	spec.setMSLevel(1);
	peak.getPosition()[0] = 20;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 25;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	// forth spectrum (MS/MS)
	spec.getContainer().clear();
	spec.setRetentionTime(12.5);
	spec.getPrecursorPeak().getPosition()[0] = 21.4;
	spec.setMSLevel(2);
	peak.getPosition()[0] = 21;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 26;
	spec.getContainer().push_back(peak);
	peak.getPosition()[0] = 31;
	spec.getContainer().push_back(peak);
	exp.push_back(spec);	

	file.store("MascotInfile_test.txt", exp, "Experiment");
	TEST_FILE("MascotInfile_test.txt", "data/MascotInfile_test_template3.txt");
	remove("MascotInfile_test.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
