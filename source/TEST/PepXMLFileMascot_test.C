// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/PepXMLFileMascot.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PepXMLFileMascot, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepXMLFileMascot* ptr = 0;
PepXMLFileMascot file;
START_SECTION(PepXMLFileMascot())
	ptr = new PepXMLFileMascot();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PepXMLFileMascot())
	delete ptr;
END_SECTION

START_SECTION(void load(const String& filename,  std::map<String, std::vector<AASequence> >& peptides))
	std::map<String, std::vector<AASequence> > peptides;
	std::map<String, std::vector<AASequence> >::iterator it;
	std::vector<AASequence> temp_sequences;
	String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFileMascot_test.pepXML");
	
	file.load(filename, peptides);
	it = peptides.begin();
	
	TEST_EQUAL(peptides.size(), 3)

	TEST_EQUAL(it->first, "402.42_14_3688.53")
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 4)
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "LAPSAAEDGAFR")
	TEST_EQUAL(temp_sequences[0].toString(), "LAPSAAEDGAFR")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "GQPLGQLEAHR")			
	TEST_EQUAL(temp_sequences[1].toString(), "GQPLGQLEAHR")
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "SPAPLVPVRLR")			
	TEST_EQUAL(temp_sequences[2].toString(), "SPAPLVPVRLR")
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "QPYLHPSFSK")
	TEST_EQUAL(temp_sequences[3].toString(), "Q(Deamidated)PYLHPSFSK")
	++it;
	TEST_EQUAL(it->first, "404.875_14_333.442")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)			
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "AQAVAAEIRER")			
	TEST_EQUAL(temp_sequences[0].toString(), "AQAVAAEIRER")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "AALNAADIVTVR")			
	TEST_EQUAL(temp_sequences[1].toString(), "AALNAADIVTVR")
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "AEPAAELALEAK")			
	TEST_EQUAL(temp_sequences[2].toString(), "AEPAAELALEAK")
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "LANAASPAITQR")			
	TEST_EQUAL(temp_sequences[3].toString(), "LANAASPAITQ(Deamidated)R")
	TEST_EQUAL(temp_sequences[4].toUnmodifiedString(), "SGHSAVLLQDGK")			
	TEST_EQUAL(temp_sequences[4].toString(), "SGHSAVLLQ(Deamidated)DGK")
	TEST_EQUAL(temp_sequences[5].toUnmodifiedString(), "AGAGIANVQAAIR")			
	TEST_EQUAL(temp_sequences[5].toString(), "AGAGIAN(Deamidated)VQ(Deamidated)AAIR")
	TEST_EQUAL(temp_sequences[6].toUnmodifiedString(), "VTAPAARSAALGK")
	TEST_EQUAL(temp_sequences[6].toString(), "VTAPAARSAALGK")
	++it;
	TEST_EQUAL(it->first, "411.766_14_3724.98")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "LLAWMGRTER")
	TEST_EQUAL(temp_sequences[0].toString(), "LLAWMGRTER")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "VLALYRAAQAR")
	TEST_EQUAL(temp_sequences[1].toString(), "VLALYRAAQ(Deamidated)AR")			
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "RTLLMSLTGLK")
	TEST_EQUAL(temp_sequences[2].toString(), "RTLLMSLTGLK")			
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "LLGLSRFGLQK")
	TEST_EQUAL(temp_sequences[3].toString(), "LLGLSRFGLQ(Deamidated)K")			
	TEST_EQUAL(temp_sequences[4].toUnmodifiedString(), "MGGIALLDEIGK")
	TEST_EQUAL(temp_sequences[4].toString(), "M(Oxidation)GGIALLDEIGK")			
	TEST_EQUAL(temp_sequences[5].toUnmodifiedString(), "DQMDNALRIR")
	TEST_EQUAL(temp_sequences[5].toString(), "DQMDN(Deamidated)ALRIR")			
	TEST_EQUAL(temp_sequences[6].toUnmodifiedString(), "QTLAGRMVVQK")
	TEST_EQUAL(temp_sequences[6].toString(), "Q(Deamidated)TLAGRMVVQ(Deamidated)K")				
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

