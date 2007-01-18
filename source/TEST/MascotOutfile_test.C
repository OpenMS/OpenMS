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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

///////////////////////////

#include <OpenMS/FORMAT/MascotOutfile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

DateTime date;
date.set("27.01.2005 17:47:41");
MascotOutfile* ptr = 0;

CHECK(MascotOutfile())
	ptr = new MascotOutfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((void load(String filename, std::vector<IdentificationData>& identifications, Real p = 0.05) throw(Exception::ParseError)))
	ptr = new MascotOutfile();
	vector<IdentificationData> identifications;

	ptr->load("data/MascotOutfile2.txt",
						identifications);

	TEST_EQUAL(identifications.size(), 4)	
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 1)	
	TEST_EQUAL(identifications[1].id.getPeptideHits().size(), 1)	
	TEST_EQUAL(identifications[2].id.getPeptideHits().size(), 10)	
	TEST_EQUAL(identifications[3].id.getPeptideHits().size(), 10)	
	TEST_REAL_EQUAL(identifications[0].id.getPeptideHits()[0].getScore(), 19.1f)	
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "NSSEA")	
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(identifications[1].id.getPeptideHits()[0].getScore(), 0.93f)	
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getSequence(), "FGASK")	
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(identifications[1].id.getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[0].getScore(), 9.72f)	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getSequence(), "AGGNAK")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[1].getScore(), 8.77f)	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getSequence(), "KGANK")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getScoreType(), "Mascot")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[1].getRank(), 2)	
	TEST_REAL_EQUAL(identifications[2].id.getPeptideHits()[2].getScore(), 8.77f)	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[2].getSequence(), "KXANK")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[2].getScoreType(), "Mascot")	
	TEST_EQUAL(identifications[2].id.getPeptideHits()[2].getRank(), 3)	
		
  TEST_REAL_EQUAL(identifications[0].rt, 88.3466f)
  TEST_REAL_EQUAL(identifications[1].rt, 96.9993f)
  TEST_REAL_EQUAL(identifications[2].rt, 105.615f)
  TEST_REAL_EQUAL(identifications[3].rt, 105.615f)

	TEST_REAL_EQUAL(identifications[0].mz, 508.119f)	
	TEST_REAL_EQUAL(identifications[1].mz, 508.458f)	
	TEST_REAL_EQUAL(identifications[2].mz, 517.267f)	
	TEST_REAL_EQUAL(identifications[3].mz, 517.324f)	

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
