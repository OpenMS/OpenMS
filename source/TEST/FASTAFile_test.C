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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FASTAFile* ptr;
START_SECTION((FASTAFile()))
	ptr = new FASTAFile();
	TEST_EQUAL(ptr == 0, false)
END_SECTION

START_SECTION((~FASTAFile()))
  delete(ptr);
END_SECTION

FASTAFile file;
vector< FASTAFile::FASTAEntry > sequences;
vector< FASTAFile::FASTAEntry >::const_iterator sequences_iterator;
FASTAFile::FASTAEntry temp_entry;

START_SECTION([FASTAFile::FASTAEntry] FASTAEntry())
	FASTAFile::FASTAEntry * ptr_e;
	ptr_e = new FASTAFile::FASTAEntry();
	TEST_EQUAL(ptr_e == 0, false)
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] FASTAEntry(String id, String desc, String seq))
FASTAFile::FASTAEntry entry("ID", "DESC", "DAVLDELNER");
TEST_EQUAL(entry.identifier, "ID")
TEST_EQUAL(entry.description, "DESC")
TEST_EQUAL(entry.sequence, "DAVLDELNER")
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] bool operator==(const FASTAEntry &rhs) const)
FASTAFile::FASTAEntry entry1("ID", "DESC", "DAVLDELNER");
FASTAFile::FASTAEntry entry2("ID", "DESC", "DAVLDELNER");
FASTAFile::FASTAEntry entry3("ID2", "DESC", "DAVLDELNER");
TEST_EQUAL(entry1==entry2, true)
TEST_EQUAL(entry1==entry3, false)
END_SECTION
	
temp_entry.identifier = String("P68509|1433F_BOVIN");
temp_entry.description = String("This is the description of the first protein");
temp_entry.sequence = 		String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") + 
		String("NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
		String("NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
		String("KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
		String("LRDNLTLWTSDQQDEEAGEGN");
		
sequences.push_back(temp_entry);	


START_SECTION((void load(const String& filename, std::vector< FASTAEntry > &data)))
	vector<FASTAFile::FASTAEntry> data;
	FASTAFile file;
	
	TEST_EXCEPTION(Exception::FileNotFound, file.load("FASTAFile_test_this_file_does_not_exist",data))
	
	file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"),data);
	sequences_iterator = data.begin();
	TEST_EQUAL(data.size(), 2)
	TEST_EQUAL(sequences_iterator->identifier, String("P68509|1433F_BOVIN"))
	TEST_EQUAL(sequences_iterator->description, String("This is the description of the first protein"))
	TEST_EQUAL(sequences_iterator->sequence, String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") + 
		String("NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
		String("NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
		String("KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
		String("LRDNLTLWTSDQQDEEAGEGN"))
	sequences_iterator++;
	TEST_EQUAL(sequences_iterator->identifier, "Q9CQV8|1433B_MOUSE")
	TEST_EQUAL(sequences_iterator->sequence, String("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTE") + 
		String("QGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICND") + 
		String("VLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQ") + 
		String("PTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLT") + 
		String("LWTSENQGDEGDAGEGEN"))
END_SECTION

START_SECTION((void store(const String& filename, const std::vector< FASTAEntry > &data) const))
	vector<FASTAFile::FASTAEntry> data, data2;
	String tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	FASTAFile file;
	
	file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"),data);
	TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.txt",data))
	
	file.store(tmp_filename,data);
	file.load(tmp_filename,data2);
	TEST_EQUAL(data==data2,true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
