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
CHECK((FASTAFile()))
	ptr = new FASTAFile();
	TEST_EQUAL(ptr == 0, false)
RESULT

CHECK((~FASTAFile()))
  delete(ptr);
RESULT

FASTAFile file();
vector< pair < String, String > > sequences;
vector< pair < String, String > >::const_iterator sequences_iterator;
	
sequences.push_back(make_pair(String("P68509|1433F_BOVIN"), 
		String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") + 
		String("NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
		String("NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
		String("KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
		String("LRDNLTLWTSDQQDEEAGEGN")));	


CHECK((void load(const String& filename, FASTAType& data) throw(Exception::FileNotFound, Exception::ParseError)))
	FASTAFile::FASTAType data;
	FASTAFile file;
	
	TEST_EXCEPTION(Exception::FileNotFound, file.load("bla",data))
	
	file.load("data/FASTAFile_test.fasta",data);
	sequences_iterator = data.begin();
	TEST_EQUAL(data.size(), 2)
	TEST_EQUAL(sequences_iterator->first, String("P68509|1433F_BOVIN"))
	TEST_EQUAL(sequences_iterator->second, String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") + 
		String("NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
		String("NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
		String("KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
		String("LRDNLTLWTSDQQDEEAGEGN"))
	sequences_iterator++;
	TEST_EQUAL(sequences_iterator->first, "Q9CQV8|1433B_MOUSE")
	TEST_EQUAL(sequences_iterator->second, String("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTE") + 
		String("QGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICND") + 
		String("VLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQ") + 
		String("PTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLT") + 
		String("LWTSENQGDEGDAGEGEN"))
RESULT

CHECK((void store(const String& filename, const FASTAType& data) const throw(Exception::UnableToCreateFile)))
	FASTAFile::FASTAType data, data2;
	String tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	FASTAFile file;
	
	file.load("data/FASTAFile_test.fasta",data);
	TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.txt",data))
	
	file.store(tmp_filename,data);
	file.load(tmp_filename,data2);
	TEST_EQUAL(data==data2,true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
