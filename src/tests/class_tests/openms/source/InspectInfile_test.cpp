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

#include <OpenMS/FORMAT/InspectInfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(InspectInfile, "$Id$")

/////////////////////////////////////////////////////////////

InspectInfile* ptr = nullptr;
InspectInfile* nullPointer = nullptr;
START_SECTION(InspectInfile())
	ptr = new InspectInfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~InspectInfile())
	delete ptr;
END_SECTION

START_SECTION((InspectInfile& operator=(const InspectInfile &inspect_infile)))
	InspectInfile inspect_infile1;
	inspect_infile1.setSpectra("dummy");
	InspectInfile inspect_infile2 = inspect_infile1;
	InspectInfile inspect_infile3;
	inspect_infile3.setSpectra("dummy");
	inspect_infile1 = InspectInfile();
	TEST_EQUAL(( inspect_infile2 == inspect_infile3 ), true)
	InspectInfile inspect_infile4;
	TEST_EQUAL(( inspect_infile1 == inspect_infile4 ), true)
END_SECTION

START_SECTION((InspectInfile(const InspectInfile &inspect_infile)))
	InspectInfile inspect_infile1;
	inspect_infile1.setSpectra("dummy");
	InspectInfile inspect_infile2(inspect_infile1);
	InspectInfile inspect_infile3;
	inspect_infile3.setSpectra("dummy");
	inspect_infile1 = InspectInfile();
	TEST_EQUAL(( inspect_infile2 == inspect_infile3 ), true)
	InspectInfile inspect_infile4;
	TEST_EQUAL(( inspect_infile1 == inspect_infile4 ), true)
END_SECTION

START_SECTION((bool operator==(const InspectInfile &inspect_infile) const))
	InspectInfile inspect_infile1;
	inspect_infile1.setSpectra("dummy");
	InspectInfile inspect_infile2;
	inspect_infile2.setSpectra("dummy");
	TEST_EQUAL(( inspect_infile1 == inspect_infile2 ), true)
END_SECTION

InspectInfile file;

START_SECTION(void setSpectra(const String& spectra))
	file.setSpectra("dummy4712");
	TEST_STRING_EQUAL(file.getSpectra(), "dummy4712")
END_SECTION

START_SECTION((const String& getSpectra() const))
	TEST_STRING_EQUAL(file.getSpectra(), "dummy4712")
END_SECTION


START_SECTION(void setDb(const String& db))
	file.setDb("dummy4711");
	TEST_STRING_EQUAL(file.getDb(), "dummy4711");
END_SECTION

START_SECTION((const String& getDb() const))
	TEST_STRING_EQUAL(file.getDb(), "dummy4711");
END_SECTION


START_SECTION(void setEnzyme(const String& enzyme))
	file.setEnzyme("Trypsin");
	TEST_STRING_EQUAL(file.getEnzyme(), "Trypsin")
END_SECTION

START_SECTION((const String& getEnzyme() const))
	TEST_STRING_EQUAL(file.getEnzyme(), "Trypsin")
END_SECTION

START_SECTION(void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic))

	// test exceptions
	String modification_line = "Phosphorylation";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.handlePTMs(modification_line, "a", true), "the file 'a' could not be found")
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotReadable, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Inspect_unreadable_unwriteable.txt"), true), "the file `data/Inspect_unreadable_unwriteable.txt' is not readable for the current user")
	
	modification_line = "2H20,KRLNH,fix";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), true), "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix")
	
	modification_line = "10.3+";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), true), "No residues for modification given. Aborting! in: 10.3+")
	
	modification_line = "10.3+,KRLNH,stat,PTM_0";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), true), "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0")
	
	modification_line = "Phosphorylation:Phosphorylation";
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), true), "There's already a modification with this name. Aborting! in: Phosphorylation")
	
	
	// test the actual program
	modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";

	// average masses
  file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), false);

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
// 	modifications["PTM_6"][2] = "OPT";
	modifications["PTM_2"] = vector< String >(3);
	modifications["PTM_2"][0] = "CTERM";
	modifications["PTM_2"][1] = "-16";
	modifications["PTM_2"][2] = "OPT";
	modifications["PTM_3"] = vector< String >(3);
	modifications["PTM_3"][0] = "NTERM";
	modifications["PTM_3"][1] = "-16";
	modifications["PTM_3"][2] = "OPT";

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
  file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), true);

//	modifications["Phosphorylation"][1] = "79.96635";
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
				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}
END_SECTION

START_SECTION((const Map< String, std::vector< String > >& getModifications() const))
	String modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";

	// average masses
	file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("Inspect_PTMs.xml"), false);

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
// 	modifications["PTM_6"][2] = "OPT";
	modifications["PTM_2"] = vector< String >(3);
	modifications["PTM_2"][0] = "CTERM";
	modifications["PTM_2"][1] = "-16";
	modifications["PTM_2"][2] = "OPT";
	modifications["PTM_3"] = vector< String >(3);
	modifications["PTM_3"][0] = "NTERM";
	modifications["PTM_3"][1] = "-16";
	modifications["PTM_3"][2] = "OPT";

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


START_SECTION(void setModificationsPerPeptide(Int modifications_per_peptide))
	file.setModificationsPerPeptide(2);
	TEST_EQUAL(file.getModificationsPerPeptide(), 2)
END_SECTION

START_SECTION((Int getModificationsPerPeptide() const))
	TEST_EQUAL(file.getModificationsPerPeptide(), 2)
END_SECTION


START_SECTION(void setBlind(UInt blind))
	file.setBlind(1);
	TEST_EQUAL(file.getBlind(), 1)
END_SECTION

START_SECTION((UInt getBlind() const))
	TEST_EQUAL(file.getBlind(), 1)
END_SECTION


START_SECTION(void setMaxPTMsize(float maxptmsize))
	file.setMaxPTMsize(250);
	TEST_EQUAL(file.getMaxPTMsize(), 250)
END_SECTION

START_SECTION((float getMaxPTMsize() const))
	TEST_EQUAL(file.getMaxPTMsize(), 250)
END_SECTION


START_SECTION(void setPrecursorMassTolerance(float precursor_mass_tolerance))
	file.setPrecursorMassTolerance(1.3f);
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance(), 1.3f)
END_SECTION

START_SECTION((float getPrecursorMassTolerance() const))
	TEST_REAL_SIMILAR(file.getPrecursorMassTolerance(), 1.3f)
END_SECTION


START_SECTION(void setPeakMassTolerance(float peak_mass_tolerance))
	file.setPeakMassTolerance(0.3f);
	TEST_REAL_SIMILAR(file.getPeakMassTolerance(), 0.3f)
END_SECTION

START_SECTION((float getPeakMassTolerance() const))
	TEST_REAL_SIMILAR(file.getPeakMassTolerance(), 0.3f)
END_SECTION


START_SECTION(void setMulticharge(UInt multicharge))
	file.setMulticharge(1);
	TEST_EQUAL(file.getMulticharge(), 1)
END_SECTION

START_SECTION((UInt getMulticharge() const))
	TEST_EQUAL(file.getMulticharge(), 1)
END_SECTION


START_SECTION(void setInstrument(const String& instrument))
	file.setInstrument("ESI-ION-TRAP");
	TEST_STRING_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
END_SECTION

START_SECTION((const String& getInstrument() const))
	TEST_STRING_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
END_SECTION


START_SECTION(void setTagCount(Int TagCount))
	file.setTagCount(1);
	TEST_EQUAL(file.getTagCount(), 1)
END_SECTION

START_SECTION((Int getTagCount() const))
	TEST_EQUAL(file.getTagCount(), 1)
END_SECTION


START_SECTION(void store(const String& filename))
	String filename;
	NEW_TMP_FILE(filename)
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.store(OPENMS_GET_TEST_DATA_PATH("Inspect_unreadable_unwriteable.txt")), "the file `data/Inspect_unreadable_unwriteable.txt' could not be created")
	
	file.store(filename);
	TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("InspectInfile_test_template1.txt"))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
