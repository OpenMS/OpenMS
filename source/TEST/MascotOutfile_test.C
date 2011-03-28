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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
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
MascotOutfile* nullPointer = 0;

START_SECTION((MascotOutfile()))
	ptr = new MascotOutfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((void load(String filename, ProteinIdentification &protein_identification, std::vector< PeptideIdentification > &peptide_identifications, Real p=0.05)))
	ptr = new MascotOutfile();
	vector<PeptideIdentification> peptide_identifications;
	ProteinIdentification protein_identification;

	ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotOutfile2.txt"), protein_identification, peptide_identifications);

	TEST_EQUAL(peptide_identifications.size(), 4)	
	TEST_EQUAL(peptide_identifications[0].getHits().size(), 1)	
	TEST_EQUAL(peptide_identifications[1].getHits().size(), 1)	
	TEST_EQUAL(peptide_identifications[2].getHits().size(), 10)	
	TEST_EQUAL(peptide_identifications[3].getHits().size(), 10)	
	TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 19.1f)	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), "NSSEA")	
	TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)	
	TEST_REAL_SIMILAR(peptide_identifications[1].getHits()[0].getScore(), 0.93f)	
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), "FGASK")	
	TEST_EQUAL(peptide_identifications[1].getHits()[0].getRank(), 1)	
	TEST_REAL_SIMILAR(peptide_identifications[2].getHits()[0].getScore(), 9.72f)	
	TEST_EQUAL(peptide_identifications[2].getHits()[0].getSequence(), "AGGNAK")	
	TEST_EQUAL(peptide_identifications[2].getHits()[0].getRank(), 1)	
	TEST_REAL_SIMILAR(peptide_identifications[2].getHits()[1].getScore(), 8.77f)	
	TEST_EQUAL(peptide_identifications[2].getHits()[1].getSequence(), "KGANK")	
	TEST_EQUAL(peptide_identifications[2].getHits()[1].getRank(), 2)	
	TEST_REAL_SIMILAR(peptide_identifications[2].getHits()[2].getScore(), 8.77f)	
	TEST_EQUAL(peptide_identifications[2].getHits()[2].getSequence(), "KXANK")	
	TEST_EQUAL(peptide_identifications[2].getScoreType(), "Mascot")	
	TEST_EQUAL(peptide_identifications[2].getHits()[2].getRank(), 3)	
		
  TEST_REAL_SIMILAR((float)peptide_identifications[0].getMetaValue("RT"), 88.3466f)
  TEST_REAL_SIMILAR((float)peptide_identifications[1].getMetaValue("RT"), 96.9993f)
  TEST_REAL_SIMILAR((float)peptide_identifications[2].getMetaValue("RT"), 105.615f)
  TEST_REAL_SIMILAR((float)peptide_identifications[3].getMetaValue("RT"), 105.615f)

	TEST_REAL_SIMILAR((float)peptide_identifications[0].getMetaValue("MZ"), 508.119f)	
	TEST_REAL_SIMILAR((float)peptide_identifications[1].getMetaValue("MZ"), 508.458f)	
	TEST_REAL_SIMILAR((float)peptide_identifications[2].getMetaValue("MZ"), 517.267f)	
	TEST_REAL_SIMILAR((float)peptide_identifications[3].getMetaValue("MZ"), 517.324f)	

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
