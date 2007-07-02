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

#include <OpenMS/FORMAT/InspectInfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(InspectInfile, "$Id$")

/////////////////////////////////////////////////////////////

InspectInfile* ptr = 0;
CHECK(InspectInfile())
	ptr = new InspectInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~InspectInfile())
	delete ptr;
RESULT

InspectInfile file;

CHECK(void setSpectra(const string& spectra))
	file.setSpectra("dummy4712");
	TEST_STRING_EQUAL(file.getSpectra(), "dummy4712")
RESULT

CHECK(const string& getSpectra() const)
	TEST_STRING_EQUAL(file.getSpectra(), "dummy4712")
RESULT


CHECK(void setDb(const String& db))
	file.setDb("dummy4711");
	TEST_STRING_EQUAL(file.getDb(), "dummy4711");
RESULT

CHECK(const String& getDb() const)
	TEST_STRING_EQUAL(file.getDb(), "dummy4711");
RESULT


CHECK(void setEnzyme(const String& protease))
	file.setEnzyme("Trypsin");
	TEST_STRING_EQUAL(file.getEnzyme(), "Trypsin")
RESULT

CHECK(const String& getEnzyme() const)
	TEST_STRING_EQUAL(file.getEnzyme(), "Trypsin")
RESULT

CHECK(String handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic) throw (Exception::FileNotReadable, Exception::FileNotFound, Exception::ParseError))
	String modification_line = "10.3+,KRLNH,fix:Phosphorylation:+16,C:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
	
	TEST_EXCEPTION(Exception::FileNotFound, file.handlePTMs(modification_line, "", true))
	
	modification_line = "2H20,KRLNH,fix";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_STRING_EQUAL(String(p_e.getMessage()), "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix")
	}

	modification_line = "10.3+";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_STRING_EQUAL(String(p_e.getMessage()), "No residues for modification given. Aborting! in: 10.3+")
	}

	modification_line = "10.3+,KRLNH,stat,PTM_0";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_STRING_EQUAL(String(p_e.getMessage()), "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0")
	}

//<<<<<<< .working
//	modification_line = "Phosphorylation:Phosphorylation";
//	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true))
//	try
//	{
//		file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
//	}
//	catch ( Exception::ParseError p_e )
//	{
//		TEST_STRING_EQUAL(String(p_e.getMessage()), "There's already a modification with this name. Aborting! in: Phosphorylation")
//	}
//
//	modification_line = "10.3+,KRLNH,fix:Phosphorylation:+16,C:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
//// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
//
//	// average masses
//	file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", false);
//
//	map< String, vector< String > > modifications;
//	modifications["PTM_0"] = vector< String >(3);
//	modifications["PTM_0"][0] = "KRLNH";
//	modifications["PTM_0"][1] = "10.3";
//	modifications["PTM_0"][2] = "FIX";
//	modifications["Phosphorylation"] = vector< String >(3);
//	modifications["Phosphorylation"][0] = "STYDHCR";
//	modifications["Phosphorylation"][1] = "79.9799";
//	modifications["Phosphorylation"][2] = "OPT";
//	modifications["PTM_2"] = vector< String >(3);
//	modifications["PTM_2"][0] = "C";
//	modifications["PTM_2"][1] = "16";
//	modifications["PTM_2"][2] = "OPT";
//// 	modifications["Carbamylation"] = vector< String >(3);
//// 	modifications["Carbamylation"][0] = "NTERM";
//// 	modifications["Carbamylation"][1] = "43.02474";
//// 	modifications["Carbamylation"][2] = "OPT";
//	modifications["Methylation"] = vector< String >(3);
//	modifications["Methylation"][0] = "CHKNQRILDEST";
//	modifications["Methylation"][1] = "14.02658";
//	modifications["Methylation"][2] = "OPT";
//// 	modifications["PTM_5"] = vector< String >(3);
//// 	modifications["PTM_5"][0] = "CTERM";
//// 	modifications["PTM_5"][1] = "-16";
//// 	modifications["PTM_5"][2] = "OPT";
//// 	modifications["PTM_6"] = vector< String >(3);
//// 	modifications["PTM_6"][0] = "NTERM";
//// 	modifications["PTM_6"][1] = "-16";
//// 	modifications["PTM_6"][2] = "OPT";
//	modifications["PTM_4"] = vector< String >(3);
//	modifications["PTM_4"][0] = "CTERM";
//	modifications["PTM_4"][1] = "-16";
//	modifications["PTM_4"][2] = "OPT";
//	modifications["PTM_5"] = vector< String >(3);
//	modifications["PTM_5"][0] = "NTERM";
//	modifications["PTM_5"][1] = "-16";
//	modifications["PTM_5"][2] = "OPT";
//
//	map< String, vector< String > >::const_iterator result_mod_i = file.getModifications().begin();
//	TEST_EQUAL(file.getModifications().size(), modifications.size())
//	if ( file.getModifications().size() == modifications.size() )
//	{
//		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
//		{
//			TEST_STRING_EQUAL(result_mod_i->first, mod_i->first)
//			TEST_EQUAL(result_mod_i->second.size(), 3)
//			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
//			if ( result_mod_i->second.size() == mod_i->second.size() )
//			{
//				TEST_STRING_EQUAL(result_mod_i->second[0], mod_i->second[0])
//				TEST_STRING_EQUAL(result_mod_i->second[1], mod_i->second[1])
//				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
//			}
//		}
//	}
//
//	// monoisotopic masses
//	file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
//
//	modifications["Phosphorylation"][1] = "79.96635";
//// 	modifications["Carbamylation"][1] = "43.00582";
//	modifications["Methylation"][1] = "14.01565";
//
//	result_mod_i = file.getModifications().begin();
//	TEST_EQUAL(file.getModifications().size(), modifications.size())
//	if ( file.getModifications().size() == modifications.size() )
//	{
//		for ( map< String, vector< String > >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i, ++result_mod_i )
//		{
//			TEST_STRING_EQUAL(result_mod_i->first, mod_i->first)
//			TEST_EQUAL(result_mod_i->second.size(), 3)
//			TEST_EQUAL(result_mod_i->second.size(), mod_i->second.size())
//			if ( result_mod_i->second.size() == mod_i->second.size() )
//			{
//				TEST_STRING_EQUAL(result_mod_i->second[0], mod_i->second[0])
//				TEST_STRING_EQUAL(result_mod_i->second[1].substr(0, 7), mod_i->second[1].substr(0, 7))
//				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
//			}
//		}
//	}
//=======
	modification_line = "Phosphorylation:Phosphorylation";
	TEST_EXCEPTION(Exception::ParseError, file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true))
	try
	{
		file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);
	}
	catch ( Exception::ParseError p_e )
	{
		TEST_STRING_EQUAL(String(p_e.getMessage()), "There's already a modification with this name. Aborting! in: Phosphorylation")
	}

	modification_line = "10.3+,KRLNH,fix:Phosphorylation:+16,C:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";

	// average masses
	file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", false);

	map< String, vector< String > > modifications;
	modifications["PTM_0"] = vector< String >(3);
	modifications["PTM_0"][0] = "KRLNH";
	modifications["PTM_0"][1] = "10.3";
	modifications["PTM_0"][2] = "FIX";
	modifications["Phosphorylation"] = vector< String >(3);
	modifications["Phosphorylation"][0] = "STYDHCR";
	modifications["Phosphorylation"][1] = "79.97990108";
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
	modifications["Methylation"][1] = "14.02658033";
	modifications["Methylation"][2] = "OPT";
// 	modifications["PTM_5"] = vector< String >(3);
// 	modifications["PTM_5"][0] = "CTERM";
// 	modifications["PTM_5"][1] = "-16";
// 	modifications["PTM_5"][2] = "OPT";
// 	modifications["PTM_6"] = vector< String >(3);
// 	modifications["PTM_6"][0] = "NTERM";
// 	modifications["PTM_6"][1] = "-16";
// 	modifications["PTM_6"][2] = "OPT";
	modifications["PTM_4"] = vector< String >(3);
	modifications["PTM_4"][0] = "CTERM";
	modifications["PTM_4"][1] = "-16";
	modifications["PTM_4"][2] = "OPT";
	modifications["PTM_5"] = vector< String >(3);
	modifications["PTM_5"][0] = "NTERM";
	modifications["PTM_5"][1] = "-16";
	modifications["PTM_5"][2] = "OPT";

	PRECISION(0.0001)

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
	file.handlePTMs(modification_line, "TOPP/Inspect_PTMs.xml", true);

	modifications["Phosphorylation"][1] = "79.96634495";
// 	modifications["Carbamylation"][1] = "43.00582";
	modifications["Methylation"][1] = "14.01565";

	PRECISION(0.0001)

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
				TEST_REAL_EQUAL(result_mod_i->second[1].toDouble(), mod_i->second[1].toDouble())
				TEST_STRING_EQUAL(result_mod_i->second[2], mod_i->second[2])
			}
		}
	}
//>>>>>>> .merge-right.r2320
RESULT


CHECK(void setModificationsPerPeptide(int modifications_per_peptide))
	file.setModificationsPerPeptide(2);
	TEST_EQUAL(file.getModificationsPerPeptide(), 2)
RESULT

CHECK(const int getModificationsPerPeptide() const)
	TEST_EQUAL(file.getModificationsPerPeptide(), 2)
RESULT


CHECK(void setBlind(unsigned int blind))
	file.setBlind(1);
	TEST_EQUAL(file.getBlind(), 1)
RESULT

CHECK(const unsigned int getBlind() const)
	TEST_EQUAL(file.getBlind(), 1)
RESULT


CHECK(void setMaxPTMsize(DoubleReal maxptmsize))
	file.setMaxPTMsize(250);
	TEST_EQUAL(file.getMaxPTMsize(), 250)
RESULT

CHECK(const DoubleReal getMaxPTMsize() const)
	TEST_EQUAL(file.getMaxPTMsize(), 250)
RESULT


CHECK(void setPrecursorMassTolerance(DoubleReal precursor_mass_tolerance))
	file.setPrecursorMassTolerance(1.3);
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance(), 1.3)
RESULT

CHECK(const DoubleReal getPrecursorMassTolerance() const)
	TEST_REAL_EQUAL(file.getPrecursorMassTolerance(), 1.3)
RESULT


CHECK(void setPeakMassTolerance(DoubleReal peak_mass_tolerance))
	file.setPeakMassTolerance(0.3);
	TEST_REAL_EQUAL(file.getPeakMassTolerance(), 0.3)
RESULT

CHECK(const DoubleReal getPeakMassTolerance() const)
	TEST_REAL_EQUAL(file.getPeakMassTolerance(), 0.3)
RESULT


CHECK(void setMulticharge(unsigned int multicharge))
	file.setMulticharge(1);
	TEST_EQUAL(file.getMulticharge(), 1)
RESULT

CHECK(const unsigned int getMulticharge() const)
	TEST_EQUAL(file.getMulticharge(), 1)
RESULT


CHECK(void setInstrument(const String& instrument))
	file.setInstrument("ESI-ION-TRAP");
	TEST_STRING_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
RESULT

CHECK(const String& getInstrument() const)
	TEST_STRING_EQUAL(file.getInstrument(), "ESI-ION-TRAP")
RESULT


CHECK(void setTagCount(int TagCount))
	file.setTagCount(1);
	TEST_EQUAL(file.getTagCount(), 1)
RESULT

CHECK(const int getTagCount() const)
	TEST_EQUAL(file.getTagCount(), 1)
RESULT


CHECK(void store(const String& filename) throw (Exception::UnableToCreateFile))
	String filename;
	NEW_TMP_FILE(filename)
	file.store(filename);
	TEST_FILE(filename.c_str(), "data/InspectInfile_test_template1.txt");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
