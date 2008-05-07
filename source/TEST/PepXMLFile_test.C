// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/PepXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PepXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepXMLFile* ptr = 0;
PepXMLFile file;
CHECK(PepXMLFile())
{
        ptr = new PepXMLFile();
        TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~PepXMLFile())
{
        delete ptr;
}
RESULT

CHECK(void load(const String& filename,  std::map<String, std::vector<AASequence> >& peptides))
{
	std::map<String, std::vector<AASequence> > peptides;
	std::map<String, std::vector<AASequence> >::iterator it;
	std::vector<AASequence> temp_sequences;
	String filename = "data/PepXMLFile_test.pepXML";
	
	file.load(filename, peptides);
	it = peptides.begin();
	
	TEST_EQUAL(peptides.size(), 3)

	TEST_EQUAL(it->first, "402.42_14_3688.53")
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 4)
	TEST_EQUAL(temp_sequences[0], "LAPSAAEDGAFR")			
	TEST_EQUAL(temp_sequences[1], "GQPLGQLEAHR")			
	TEST_EQUAL(temp_sequences[2], "SPAPLVPVRLR")			
	TEST_EQUAL(temp_sequences[3], "QPYLHPSFSK")			
	++it;
	TEST_EQUAL(it->first, "404.875_14_333.442")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)			
	TEST_EQUAL(temp_sequences[0], "AQAVAAEIRER")			
	TEST_EQUAL(temp_sequences[1], "AALNAADIVTVR")			
	TEST_EQUAL(temp_sequences[2], "AEPAAELALEAK")			
	TEST_EQUAL(temp_sequences[3], "LANAASPAITQR")			
	TEST_EQUAL(temp_sequences[4], "SGHSAVLLQDGK")			
	TEST_EQUAL(temp_sequences[5], "AGAGIANVQAAIR")			
	TEST_EQUAL(temp_sequences[6], "VTAPAARSAALGK")			
	++it;
	TEST_EQUAL(it->first, "411.766_14_3724.98")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)
	TEST_EQUAL(temp_sequences[0], "LLAWMGRTER")			
	TEST_EQUAL(temp_sequences[1], "VLALYRAAQAR")			
	TEST_EQUAL(temp_sequences[2], "RTLLMSLTGLK")			
	TEST_EQUAL(temp_sequences[3], "LLGLSRFGLQK")			
	TEST_EQUAL(temp_sequences[4], "MGGIALLDEIGK")			
	TEST_EQUAL(temp_sequences[5], "DQMDNALRIR")			
	TEST_EQUAL(temp_sequences[6], "QTLAGRMVVQK")				
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

