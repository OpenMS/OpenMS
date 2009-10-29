// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/PepNovoInfile.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PepNovoInfile, "$Id$")

/////////////////////////////////////////////////////////////

PepNovoInfile* ptr = 0;
START_SECTION(PepNovoInfile())
	ptr = new PepNovoInfile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PepNovoInfile())
	delete ptr;
END_SECTION

START_SECTION((PepNovoInfile& operator=(const PepNovoInfile& pepnovo_infile)))
	PepNovoInfile pepnovo_infile1;
  PepNovoInfile pepnovo_infile2;
	pepnovo_infile2 = pepnovo_infile1;
	PepNovoInfile pepnovo_infile3;
	pepnovo_infile1 = PepNovoInfile();
	TEST_EQUAL(( pepnovo_infile2 == pepnovo_infile3 ), true)
END_SECTION

START_SECTION((PepNovoInfile(const PepNovoInfile &pepnovo_infile)))
	PepNovoInfile pepnovo_infile1;
	PepNovoInfile pepnovo_infile2(pepnovo_infile1);
	PepNovoInfile pepnovo_infile3;
	pepnovo_infile1 = PepNovoInfile();
	TEST_EQUAL(( pepnovo_infile2 == pepnovo_infile3 ), true)
END_SECTION

START_SECTION((bool operator==(const PepNovoInfile &pepnovo_infile) const))
	PepNovoInfile pepnovo_infile1;
	PepNovoInfile pepnovo_infile2;
	TEST_EQUAL(( pepnovo_infile1 == pepnovo_infile2 ), true)
END_SECTION

PepNovoInfile file;

START_SECTION(void handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic))
	//TODO add real test
/*
	// test exceptions
	String modification_line = "Phosphorylation";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.handlePTMs(modification_line, "a", true), "the file 'a' could not be found")
	
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotReadable, file.handlePTMs(modification_line, OPENMS_GET_TEST_DATA_PATH("PepNovo_unreadable_unwriteable.txt"), true), "the file `data/PepNovo_unreadable_unwriteable.txt' is not readable for the current user")
	
	modification_line = "2H20,KRLNH,fix";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), true), "There's something wrong with this modification. Aborting! in: 2H20,KRLNH,fix")
	
	modification_line = "10.3+";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), true), "No residues for modification given. Aborting! in: 10.3+")
	
	modification_line = "10.3+,KRLNH,stat,PTM_0";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), true), "There's something wrong with the type of this modification. Aborting! in: 10.3+,KRLNH,stat,PTM_0")
	
	modification_line = "Phosphorylation:Phosphorylation";
	TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), true), "There's already a modification with this name. Aborting! in: Phosphorylation")
	
	
	// test the actual program
	modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";

	// average masses
	file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), false);
*/
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
/*
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
*/
	// monoisotopic masses
	//file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), true);

//	modifications["Phosphorylation"][1] = "79.96635";
// 	modifications["Carbamylation"][1] = "43.00582";
//	modifications["Methylation"][1] = "14.01565";
/*
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
	*/
END_SECTION

START_SECTION((const std::map< String, std::vector< String > >& getModifications() const))
//TODO add real test
	String modification_line = "10.3+,KRLNH,fix:+16,C:16-,cterm:-16,nterm";
// 	"10.3+,KRLNH,fix:Phosphorylation:+16,C:HCNO,nterm,Carbamylation:H2C,CHKNQRILDEST,opt,Methylation:16-,cterm:-16,nterm";

	// average masses
	//file.handlePTMs(modification_line,  OPENMS_GET_TEST_DATA_PATH("../TOPP/PepNovo_PTMs.xml"), false);

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
/*
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
	*/
END_SECTION

START_SECTION(String store(const String& filename))
//TODO add real test
	String filename;
	NEW_TMP_FILE(filename)

	// test exceptions
// 	TEST_EXCEPTION_WITH_MESSAGE(Exception::UnableToCreateFile, file.store(OPENMS_GET_TEST_DATA_PATH("PepNovo_unreadable_unwriteable.txt")), "the file `data/PepNovo_unreadable_unwriteable.txt' could not be created")
	
	
	// test actual program
	file.store(filename);
	//TEST_FILE_EQUAL(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("PepNovoInfile_test_template.txt"));
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
